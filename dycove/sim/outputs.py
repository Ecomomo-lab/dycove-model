###############################################################
#  outputs.py
###############################################################

from pathlib import Path
import numpy as np
from dataclasses import asdict


class OutputManager:
    """ For saving :class:`~dycove.sim.vegetation_data.VegCohort` instances to output files. """
    def __init__(self, engine):
        self.engine = engine
        self.veg = engine.veg
        self.veg_dir   = Path(engine.model_dir) / "veg_output"
        self.veg_dir.mkdir(parents=True, exist_ok=True)


    def save_vegetation_step(self, year: int, ets: int):
        """ Save vegetation cohort state for a given ecological timestep """
        for n, cohort in enumerate(self.veg.cohorts):
            fname = f"cohort{n+1}_proc{self.engine.get_rank()}" if self.engine.is_parallel() else f"cohort{n+1}"
            self.save_binary(self.veg_dir, fname, asdict(cohort), year, ets)

    def save_binary(self, directory: Path, filename: str, data: dict, year: int, ets: int):
        """ Save dict-like data as compressed numpy .npz file """
        data_flat = self.make_numpy_friendly(data)
        fname = f"{filename}_year{year}_ets{ets}.npz"
        np.savez_compressed(directory / fname, **data_flat)

    def make_numpy_friendly(self, obj):
        """ Recursively convert arrays/lists to numpy arrays where possible """
        if isinstance(obj, np.ndarray):
            return obj
        elif isinstance(obj, dict):
            return {k: self.make_numpy_friendly(v) for k, v in obj.items()}
        elif isinstance(obj, (list, tuple)):
            try:
                return np.asarray(obj)
            except Exception:
                return [self.make_numpy_friendly(v) for v in obj]
        else:
            return obj

    def reconcile_vegetation_output(self):
        """ 
        Merge vegetation outputs across MPI subdomains into single files.

        We want one output file per cohort, per ecological timestep.

        This executes on main processor only and can take a while for giant domains.
        
        Method is part of the engine class because it uses domain-specific info to 
        re-map files, but we should keep this in OutputManager for now because we may 
        add other similar tasks here later (and b/c it involves I/O and directory 
        access).
        """
        
        if self.veg and self.engine.is_parallel() and self.engine.get_rank() == 0:
            self.engine.merge_parallel_veg(self)