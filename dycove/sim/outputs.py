###############################################################
#  outputs.py
###############################################################

from pathlib import Path
import numpy as np
from dataclasses import asdict
import xarray as xr
import json

from dycove.utils.simulation_reporting import Reporter

r = Reporter()


class OutputManager:
    COHORT_INDEX_FNAME = "_cohort_files_ets_index.json"
    METADATA_FNAME = "_eco_time_vars.json"

    """ For saving :class:`~dycove.sim.vegetation_data.VegCohort` instances to output files. """
    def __init__(self, engine, save_freq=1, save_mort=True):
        self.engine = engine
        self.veg = engine.veg
        self.veg_dir   = Path(engine.model_dir) / "veg_output"
        self.save_freq = save_freq
        self.save_mort = save_mort

        if self.veg:
            self.veg_dir.mkdir(parents=True, exist_ok=True)
        
        self.n_cohort_steps = []  # for numbering output files
        self.cohort_index = {}  # nested dict index of output files for each year-ets combo


    def save_vegetation_step(self, simstate, vts):
        """ Save vegetation cohort state for a given ecological timestep """
        for c in self.veg.cohorts:

            fname_base = f"cohort{c.cohort_id}_{c.n_ets:03d}"
            fname = (fname_base + f"_proc{self.engine.get_rank()}" 
                     if self.engine.is_parallel() 
                     else fname_base
                     )
            
            if vts % self.save_freq == 0:            
                self.save_netcdf(self.veg_dir, 
                                 fname, 
                                 asdict(c), 
                                 eco_year = simstate.eco_year, 
                                 ets = simstate.ets, 
                                 cohort_id = c.cohort_id,
                                 )
                self.cohort_indexing(fname_base, simstate.eco_year, simstate.ets)
                self.save_cohort_index()  # saves every step in case of incomplete simulation


    def cohort_indexing(self, fname_base, year, ets):
        """ Log the year/ets combo and create (or append to) list for cohort file names """
        if f"{year}" not in self.cohort_index:
            self.cohort_index[f"{year}"] = {}
        if f"{ets}" not in self.cohort_index[f"{year}"]:
            self.cohort_index[f"{year}"][f"{ets}"] = [fname_base]
        else:
            self.cohort_index[f"{year}"][f"{ets}"].append(fname_base)


    def save_simulation_metadata(self, simstate):
        """ Save record of ecological time inputs for easier post-processing """
        sim_dict = {"n_ets": simstate.n_ets,
                    "veg_interval": simstate.veg_interval,
                    "ecofac": simstate.ecofac,
                    "save_frequency": self.save_freq,
                    "save_mortality": self.save_mort,
                    }
        with open(self.veg_dir / self.METADATA_FNAME, "w") as f:
            json.dump(sim_dict, f)


    def save_cohort_index(self):
        c_index_as_str = {str(key): value for key, value in self.cohort_index.items()}
        with open(self.veg_dir / self.COHORT_INDEX_FNAME, "w") as f:
            json.dump(c_index_as_str, f)


    def save_netcdf(self,
                    directory: Path,
                    filename: str,
                    data: dict,
                    eco_year: int | None = None,
                    ets: int | None = None,
                    cohort_id: int | None = None,
                    saved_attrs: dict | None = None
                    ):
        """ Save dict-like data as a NetCDF file using xarray """

        data_vars = {}
        attrs = {}

        for key, value in data.items():
            # Pass arrays and scalars separately
            if isinstance(value, np.ndarray):
                data_vars[key] = xr.DataArray(value)
                if not self.save_mort and "mort" in key:
                    data_vars.pop(key)
            else:
                attrs[key] = value

        ds = xr.Dataset(data_vars=data_vars)

        if saved_attrs is None:  # do this when call comes from save_vegetation_step()
            ds.attrs.update(attrs)
            ds.attrs.update(eco_year=eco_year, ets=ets, cohort=cohort_id)
        else:  # do this when call comes from merge_parallel_veg() at end of simulation
            ds.attrs.update(saved_attrs)

        # Need engine='scipy' because the DFM BMI ctypes wrapper has a path conflict with netCDF
        ds.to_netcdf(directory / (filename + ".nc"), engine="scipy")


    def reconcile_vegetation_output(self, simstate):
        """ 
        Merge vegetation outputs across MPI subdomains into single files.

        We want one output file per cohort, per ecological timestep.

        This executes on main processor only and can take a while for giant domains.
        
        Method is part of the engine class because it uses domain-specific info to 
        re-map files, but we should keep this access point in OutputManager for now 
        because we may add other similar tasks here later (and b/c it involves I/O 
        and directory access).
        """
        
        if self.veg and self.engine.is_parallel() and self.engine.get_rank() == 0:
            self.file_index = self.load_cohort_file_index()
            self.engine.merge_parallel_veg(self)  # requires self object as argument
                             

    def load_cohort_file_index(self):
        # Load as separate method for testing purposes
        with open(self.veg_dir / self.COHORT_INDEX_FNAME, "r") as f:
            return json.load(f)