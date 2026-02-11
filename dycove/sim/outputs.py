###############################################################
#  outputs.py
###############################################################

from pathlib import Path
import numpy as np
from dataclasses import asdict
import xarray as xr

from dycove.utils.simulation_reporting import Reporter

r = Reporter()


class OutputManager:
    """ For saving :class:`~dycove.sim.vegetation_data.VegCohort` instances to output files. """
    def __init__(self, engine):
        self.engine = engine
        self.veg = engine.veg
        self.veg_dir   = Path(engine.model_dir) / "veg_output"
        self.veg_dir.mkdir(parents=True, exist_ok=True)
        self.n_cohort_steps = []  # for numbering output files
        self.fname_base = "cohort0_01"

    def save_vegetation_step(self, sim):
        """ Save vegetation cohort state for a given ecological timestep """
        self.update_file_counts()
        for i, cohort in enumerate(self.veg.cohorts):

            self.fname_base = f"cohort{i}_{self.n_cohort_steps[i]:02d}"
            fname = (self.fname_base + f"_proc{self.engine.get_rank()}" 
                     if self.engine.is_parallel() 
                     else self.fname_base
                     )
            
            self.save_netcdf(self.veg_dir, 
                             fname, 
                             asdict(cohort), 
                             eco_year = sim.eco_year, 
                             ets = sim.ets, 
                             cohort_id = i,
                             )
            self.n_cohort_steps[i] += 1

    def update_file_counts(self):
        """ Update file count for each cohort based on current number of cohorts. """
        if len(self.veg.cohorts) > len(self.n_cohort_steps):
            self.n_cohort_steps.extend([0]*(len(self.veg.cohorts) - len(self.n_cohort_steps)))

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
            else:
                attrs[key] = value

        ds = xr.Dataset(data_vars=data_vars)

        if saved_attrs is None:  # do this when call comes from save_vegetation_step()
            ds.attrs.update(attrs)
            ds.attrs.update(eco_year=eco_year, ets=ets, cohort=cohort_id)
        else:                    # do this when call comes from merge_parallel_veg()
            ds.attrs.update(saved_attrs)

        # Need engine='scipy' because the DFM BMI ctypes wrapper has a path conflict with netCDF
        ds.to_netcdf(directory / (filename + ".nc"), engine="scipy")


    def reconcile_vegetation_output(self):
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
            self.engine.merge_parallel_veg(self)