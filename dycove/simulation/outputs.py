from pathlib import Path
import numpy as np
from dataclasses import asdict


""" For saving variables to output files """
class OutputManager:
    # vegetation quantities after being spatially averaged and applied to hydrodynamic model
    VEG_MODEL_VARS = ["stemdensity", "stemdiameter", "stemheight"]
    # # vegetation quantities as they are stored in the vegetation module within VegCohort objects
    # VEG_COHORT_VARS = ["fraction", "height", "diameter", "density"]
    # vegetation mortaility arrays as they are stored in the vegetation module (cumulatively)
    VEG_MORTALITY_VARS = ["from_flood", "from_desic", "from_uproot", "from_burial", "from_scour", "total"]
    
    def __init__(self, engine):
        self.engine = engine
        self.veg = engine.veg
        #self.hydro_dir = Path(engine.model_dir) / "hydro_output"
        self.veg_dir   = Path(engine.model_dir) / "veg_output"
        #self.hydro_dir.mkdir(parents=True, exist_ok=True)
        self.veg_dir.mkdir(parents=True, exist_ok=True)


    def save_vegetation_step(self, year: int, ets: int):
        for n, cohort in enumerate(self.veg.cohorts):
            fname = f"cohort{n+1}_proc{self.engine.get_rank()}" if self.engine.is_parallel() else f"cohort{n+1}"
            self.save_binary(self.veg_dir, fname, asdict(cohort), year, ets)

    # def save_final_vegetation(self):
    #     # mortality variables which are cumulative can be saved at the end
    #     for mort_var in self.VEG_MORTALITY_VARS:
    #         self.save_json(self.veg_dir, f"mort_{mort_var}", getattr(self.veg, f"mort_{mort_var}"))

    # @staticmethod
    # def save_array(directory: Path, name: str, arr, i=None):
    #     fname = f"{name}_{i}.npy" if i is not None else f"{name}.npy"
    #     np.save(directory / fname, arr)

    def save_binary(self, directory: Path, name: str, data, year, ets):
        """Save dict-like data as compressed NumPy .npz file"""
        fname = f"{name}_year{year}_ets{ets}.npz"  #if i is not None else f"{name}.npz"

        # Flatten dataclass → dict of NumPy arrays/scalars
        flat = self.make_numpy_friendly(data)
        np.savez_compressed(directory / fname, **flat)

    def make_numpy_friendly(self, obj):
        """Recursively convert arrays/lists to NumPy arrays where possible"""
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
        This executes on the main processor only and can take a while for giant parallelized domains.
        Method is part of engine because it uses domain-specific info to re-map files.
        Keeping this in OutputManager for now because we may add other similar tasks here later 
          (and b/c it involves I/O and directory access).
        """
        if self.veg and self.engine.is_parallel() and self.engine.get_rank() == 0:
            self.engine.merge_parallel_veg(self)