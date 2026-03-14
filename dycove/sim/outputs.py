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
    """ For saving :class:`~dycove.sim.vegetation_data.VegCohort` instances to output files. """
    def __init__(self, engine, save_freq=1, save_mort=True):
        self.engine = engine
        self.veg = engine.veg
        self.veg_dir   = Path(engine.model_dir) / "veg_output"
        self.save_freq = save_freq
        self.save_mort = save_mort

        self.veg_dir.mkdir(parents=True, exist_ok=True)
        
        self.n_cohort_steps = []  # for numbering output files
        self.cohort_index = {}  # nested dict index of output files for each year-ets combo


    def save_vegetation_step(self, simstate, vts):
        """ Save vegetation cohort state for a given ecological timestep """
        self.update_file_counts()
        for i, cohort in enumerate(self.veg.cohorts):

            self.fname_base = f"cohort{i}_{self.n_cohort_steps[i]:02d}"
            fname = (self.fname_base + f"_proc{self.engine.get_rank()}" 
                     if self.engine.is_parallel() 
                     else self.fname_base
                     )
            
            if vts % self.save_freq == 0:            
                self.save_netcdf(self.veg_dir, 
                                 fname, 
                                 asdict(cohort), 
                                 eco_year = simstate.eco_year, 
                                 ets = simstate.ets, 
                                 cohort_id = i,
                                 )
                self.cohort_indexing(simstate.eco_year, simstate.ets)

            self.n_cohort_steps[i] += 1


    def update_file_counts(self):
        """ Update file count for each cohort based on current number of cohorts. """
        if len(self.veg.cohorts) > len(self.n_cohort_steps):
            self.n_cohort_steps.extend([0]*(len(self.veg.cohorts) - len(self.n_cohort_steps)))


    def cohort_indexing(self, year, ets):
        """ Log the year/ets combo and create (or append to) list for cohort file names """
        if f"{year}" not in self.cohort_index:
            self.cohort_index[f"{year}"] = {}
        if f"{ets}" not in self.cohort_index[f"{year}"]:
            self.cohort_index[f"{year}"][f"{ets}"] = [self.fname_base]
        else:
            self.cohort_index[f"{year}"][f"{ets}"].append(self.fname_base)


    def save_simulation_indices(self, simstate):
        self.save_simulation_metadata(simstate)
        self.save_cohort_index()


    def save_simulation_metadata(self, simstate):
        """ Save record of ecological time inputs for easier post-processing """
        sim_dict = {"n_ets": simstate.n_ets,
                    "veg_interval": simstate.veg_interval,
                    "ecofac": simstate.ecofac,
                    "save_frequency": self.save_freq,
                    "save_mortality": self.save_mort,
                    }
        with open(self.veg_dir / "_eco_time_vars.json", "w") as f:
            json.dump(sim_dict, f)


    def save_cohort_index(self):
        c_index_as_str = {str(key): value for key, value in self.cohort_index.items()}
        with open(self.veg_dir / "_cohort_files_ets_index.json", "w") as f:
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
        
        if self.veg:
            if not self.engine.is_parallel() or self.engine.get_rank() == 0:
                self.save_simulation_indices(simstate)
                if self.engine.is_parallel():
                    self.engine.merge_parallel_veg(self)  # requires self object as argument