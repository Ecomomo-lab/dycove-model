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
        self.cohort_count = []  # for numbering output files

    def save_vegetation_step(self, sim):
        """ Save vegetation cohort state for a given ecological timestep """
        self.update_file_counts()
        for i, cohort in enumerate(self.veg.cohorts):

            fname = (f"cohort{i}_{self.cohort_count[i]:02d}_proc{self.engine.get_rank()}" 
                     if self.engine.is_parallel() 
                     else f"cohort_data_{self.cohort_count[i]:02d}"
                     )
            
            self.save_netcdf(self.veg_dir, 
                             fname, 
                             asdict(cohort), 
                             eco_year = sim.eco_year, 
                             ets = sim.ets, 
                             cohort_id = i,
                             )
            self.cohort_count[i] += 1

    def update_file_counts(self):
        """ Update file count for each cohort based on current number of cohorts. """
        if len(self.veg.cohorts) > len(self.cohort_count):
            self.cohort_count.extend([0]*(len(self.veg.cohorts) - len(self.cohort_count)))

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

        #for key, value in data_flat.items():
        for key, value in data.items():
            if isinstance(value, np.ndarray):
                # Let xarray infer dimensions
                data_vars[key] = xr.DataArray(value)
            else:
                # Scalars / metadata go into attributes
                attrs[key] = value

        ds = xr.Dataset(data_vars=data_vars)

        if saved_attrs is None:  # do this when call comes from save_vegetation_step()
            ds.attrs.update(attrs)
            ds.attrs.update(eco_year=eco_year, ets=ets, cohort=cohort_id)
        else:                    # do this when call comes from merge_parallel_veg()
            ds.attrs.update(saved_attrs)

        ds.to_netcdf(directory / (filename + ".nc"))


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