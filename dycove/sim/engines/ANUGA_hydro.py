###############################################################
#  ANUGA_hydro.py
###############################################################

from datetime import datetime
import numpy as np
import xarray as xr
from pathlib import Path


from dycove.sim.base import HydroSimulationBase, HydroEngineBase
from dycove.sim.engines.ANUGA_baptist import Baptist_operator
from dycove.utils.simulation_reporting import Reporter
from dycove.constants import H_LIM_VELOCITY


r = Reporter()

def _import_anuga():
    """ Lazy loading of ANUGA parallel methods to avoid import errors when ANUGA will not be tested/used. """
    try:
        from anuga import myid, numprocs, finalize, barrier
        return myid, numprocs, finalize, barrier
    except ImportError:
        msg = ("The `anuga` package is not installed. "
               "Refer to the documentation for installation instructions.")
        r.report(msg, level="ERROR")
        raise ImportError(msg)


class ANUGA(HydroSimulationBase):
    """
    Hydrodynamic simulation wrapper for the ANUGA model.

    This class connects the generic :class:`~dycove.sim.base.HydroSimulationBase` 
    interface to the ANUGA hydrodynamic engine 
    :class:`~dycove.sim.base.engines.ANUGA_hydro.AnugaEngine`.

    Notes
    -----
    - The ANUGA ``domain`` object is expected to be pre-constructed by the user.
      This preserves the typical ANUGA workflow where domains are created
      directly in Python scripts.
    - All higher-level logic that can be abstracted from the engine classes is
      handled in :class:`~dycove.sim.base.HydroSimulationBase`; all low-level 
      model interactions are delegated to 
      :class:`~dycove.sim.base.engines.ANUGA_hydro.AnugaEngine`.

    """

    def __init__(self, anuga_domain, vegetation=None):
        # build ANUGA engine
        engine = AnugaEngine(anuga_domain, vegetation)
        # pass ANUGA engine to the base class
        super().__init__(engine)
    

class AnugaEngine(HydroEngineBase):
    """ 
    Engine interface for ANUGA hydrodynamic model.

    This engine:
    
    - Holds the ANUGA domain object.
    - Reads/writes flow and vegetation state directly through Python objects
      (unlike BMI-based engine used for Delft3D FM).
    - Contains all methods that are specific to the ANUGA model.

    Parameters
    ----------
    anuga_domain : anuga.domain
        The pre-constructed computational domain.
    vegetation : VegetationSpecies or MultipleVegetationSpecies, optional
        Vegetation object passed down from the base simulation.

    Notes
    -----
    - ANUGA runs in many short `domain.evolve()` loops (see 
      :meth:`~dycove.sim.engines.ANUGA_hydro.AnugaEngine.step`) 
      rather than one long simulation call; the ``skip_step`` mechanism 
      prevents duplicate first timesteps across loops.
    - Parallel execution (if enabled) requires merging local vegetation states
      after simulation; see 
      :meth:`~dycove.sim.engines.ANUGA_hydro.AnugaEngine.merge_parallel_veg`.
    """

    def __init__(self, anuga_domain, vegetation=None):
        # Lazy anuga loading
        self.myid, self.numprocs, self.finalize, self.barrier = _import_anuga()

        self.domain = anuga_domain
        self.model_dir = self.domain.get_datadir()

        # Passing vegetation as attribute of the engine
        self.veg = vegetation

        # With DYCOVE-ANUGA, we run many consecutive "domain.evolve" loops, rather than just one big loop.
        # We don't want to skip the first step for the first of these loops, but for every other loop, we do.
        # Otherwise, we get repeated steps.
        self.skip_step = False

        # Interval (seconds) for saving ANUGA output, this is just a placeholder.
        # Actual value is set via run_simulation() and used as argument to domain.evolve()
        self.save_interval = 3600


    def initialize(self):
        # ANUGA doesn't have an "initialize" method like DFM.
        # But we can include some required steps here rather than just having an empty method.
        if self.veg is not None:
            self.Baptist = Baptist_operator(self.domain, 
                                            veg_diameter=0, 
                                            veg_density=0, 
                                            veg_height=0, 
                                            drag=self.veg.get_drag())
        self.morphology = False
        

    def step(self, seconds):
        # Normally, all processes in ANUGA would be performed within the domain.evolve() loop at each
        #   yieldstep, but for consistency across all potential models, we wrap it up here under the 
        #   "step" method.
        # If performing a "big" step, reduce yieldstep so it equals save_interval
        yieldstep = min(seconds, self.save_interval)
        self.barrier()
        for t in self.domain.evolve(yieldstep=yieldstep, 
                                    outputstep=self.save_interval,
                                    duration=seconds,
                                    skip_initial_step=self.skip_step):
            if self.myid == 0:
                r.report(self.domain.timestepping_statistics())    
            # Skip initial yieldstep for all future loops, to avoid rerunning yieldsteps at "restart"
            self.skip_step = True

        # Enforces wait time for all cores so they catch up to each other when this is called (ignored if not parallel)
        self.barrier()


    def cleanup(self):
        if self.is_parallel():
            self.domain.sww_merge(delete_old=True)
            self.finalize()
        else:
            pass


    def get_refdate(self):
        # Hardcoded for now, because DFM requires this in input file. TODO: fix or remove
        return datetime(2001, 1, 1)   
    

    def get_cell_count(self):
        return len(self.get_elevation())  


    def get_elevation(self):
        return self.domain.quantities["elevation"].centroid_values


    def get_velocity_and_depth(self):
        stage = self.domain.quantities["stage"].centroid_values
        depth = stage - self.get_elevation()
        xmom = self.domain.quantities["xmomentum"].centroid_values
        ymom = self.domain.quantities["ymomentum"].centroid_values
        # Convert depth-averaged momentum to velocity
        with np.errstate(divide="ignore", invalid="ignore"):
            # Ignore velocities where depth is insufficient
            xvel = np.where(depth < H_LIM_VELOCITY, 0., xmom/depth)
            yvel = np.where(depth < H_LIM_VELOCITY, 0., ymom/depth)
        velocity = np.sqrt(xvel**2 + yvel**2)
        return velocity, depth


    def get_vegetation(self):
        # Pull directly from Baptist operator
        stemdensity = self.Baptist.veg_density.centroid_values
        stemdiameter = self.Baptist.veg_diameter.centroid_values
        stemheight = self.Baptist.veg_height.centroid_values
        return stemdensity, stemdiameter, stemheight
    

    def set_vegetation(self, stemdensity, stemdiameter, stemheight):
        # Update Baptist operator with new quantities
        self.Baptist.set_vegetation(veg_diameter=stemdiameter,
                                    veg_density=stemdensity,
                                    veg_height=stemheight,
                                    )
            

    def check_simulation_inputs(self, simstate):
        # Nothing implemented yet
        pass


    # --------------------------------------------------------
    # Parallel methods
    # --------------------------------------------------------
    
    def get_rank(self):
        return self.myid
    
    
    def is_parallel(self):
        return True if self.numprocs > 1 else False


    def merge_parallel_veg(self, OutputManager):

        outputdir = OutputManager.veg_dir
        n_cohort_steps = OutputManager.n_cohort_steps

        sww_name = self.domain.get_name()
        base_name = sww_name[:-2]
        sww_dir = Path(self.domain.get_datadir())

        sww_file_list = [f"{sww_dir}/{base_name}_{p}.sww" for p in range(self.numprocs)]

        tri_l2g, tri_full_flag = [], []
        for sww_file in sww_file_list:
            with xr.open_dataset(sww_file) as sww:
                tri_l2g.append(sww["tri_l2g"].values)
                tri_full_flag.append(sww["tri_full_flag"].values)

        n_global = int(max(tri.max() for tri in tri_l2g) + 1)

        f_ids = [np.where(tri_full_flag[p] == 1)[0] for p in range(self.numprocs)]
        f_gids = [tri_l2g[p][f_ids[p]] for p in range(self.numprocs)]

        for cohort_id in range(len(self.veg.cohorts)):
            for file_num in range(n_cohort_steps[cohort_id]):

                local_data = []
                c_merged_attrs = None

                for p in range(self.numprocs):
                    c_file_sub = outputdir / f"cohort{cohort_id}_{file_num:02d}_proc{p}.nc"

                    if not c_file_sub.exists():
                        msg = ("No individual processor (proc) output files found for this "
                            "cohort, year, and ETS combination")
                        r.report(msg, level="ERROR")
                        raise FileNotFoundError(msg)

                    with xr.open_dataset(c_file_sub) as c_sub:
                        if p == 0:
                            c_merged_attrs = dict(c_sub.attrs)

                        local_data.append(
                            {k: v.values for k, v in c_sub.data_vars.items()}
                        )

                    c_file_sub.unlink()

                c_merged = self.merge_local_to_global(local_data=local_data,
                                                      f_gids=f_gids,
                                                      f_ids=f_ids,
                                                      n_global=n_global,
                                                      )

                OutputManager.save_netcdf(outputdir,
                                          f"cohort{cohort_id}_{file_num:02d}",
                                          c_merged,
                                          saved_attrs=c_merged_attrs,
                                          )


    @staticmethod
    def merge_local_to_global(local_data: list[dict[str, np.ndarray]],
                              f_gids: list[np.ndarray],
                              f_ids: list[np.ndarray],
                              n_global: int,
                              ) -> dict[str, np.ndarray]:
        """
        Merge per-processor local arrays into global arrays using
        local-to-global index mappings.
        """
        merged: dict[str, np.ndarray] = {}

        for p, local_vars in enumerate(local_data):
            for key, values in local_vars.items():
                if key not in merged:
                    merged[key] = np.zeros(n_global, dtype=float)

                merged[key][f_gids[p]] = values[f_ids[p]]

        return merged
