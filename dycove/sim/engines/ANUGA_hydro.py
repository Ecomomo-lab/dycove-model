###############################################################
#  ANUGA_hydro.py
###############################################################

from datetime import datetime
import numpy as np
from netCDF4 import Dataset  # type: ignore
from pathlib import Path
import re

from dycove.sim.base import HydroSimulationBase, HydroEngineBase
from dycove.sim.engines.ANUGA_baptist import Baptist_operator
from dycove.utils.simulation_reporting import Reporter

from anuga import myid, numprocs, finalize, barrier

r = Reporter()


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
      :meth:`~dycove.sim.engines.ANUGA_hydro.AnugaEngine.merge_parallel_veg`) 
      rather than one long simulation call; the ``skip_step`` mechanism 
      prevents duplicate first timesteps across loops.
    - Parallel execution (if enabled) requires merging local vegetation states
      after simulation; see 
      :meth:`~dycove.sim.engines.ANUGA_hydro.AnugaEngine.merge_parallel_veg`.
    """

    def __init__(self, anuga_domain, vegetation=None):
        self.domain = anuga_domain
        self.model_dir = self.domain.get_datadir()

        # passing vegetation as attribute of the engine
        self.veg = vegetation

        # With DYCOVE-ANUGA, we run many consecutive "domain.evolve" loops, rather than just one big loop.
        # We don't want to skip the first step for the first of these loops, but for every other loop, we do.
        # Otherwise, we get repeated steps.
        self.skip_step = False

        # interval (seconds) for saving ANUGA output, this is just a placeholder
        # actual value is set via run_simulation() and used as argument to domain.evolve()
        self.save_interval = 3600


    def initialize(self):
        # ANUGA doesn't have an "initialize" method like DFM, 
        # but we can include some required steps here rather than just having an empty method
        self.Baptist = Baptist_operator(self.domain, veg_diameter=0, veg_density=0, veg_height=0)
        self.morphology = False
        
    def step(self, seconds):
        # normally, all processes would be performed within this domain.evolve() loop, 
        # but for consistency across all potential models, we wrap it up here under the "step" method.
        yieldstep = min(seconds, self.save_interval)  # if performing a "big" step, reduce yieldstep so it equals save_interval
        barrier()
        for t in self.domain.evolve(yieldstep=yieldstep, 
                                    outputstep=self.save_interval,
                                    duration=seconds,
                                    skip_initial_step=self.skip_step):
            if myid == 0:
                r.report(self.domain.timestepping_statistics())
                
            # skip initial yieldstep for all future loops, to avoid rerunning yieldsteps at restart
            self.skip_step = True
        # barrier() enforces wait time for all cores so they catch up to each other when this is called (ignored if not parallel)
        barrier()

    def cleanup(self):
        # cleanup only required for parallel models
        if self.is_parallel():
            self.domain.sww_merge(delete_old=True)
            finalize()
        else:
            pass

    def get_rank(self):
        return myid

    def get_cell_count(self):
        # get number of cells (not nodes) in numerical model grid
        return len(self.get_elevation())

    def get_refdate(self):
        # hardcoded for now, because DFM requires this in input file
        # TODO: fix or remove
        return datetime(2001, 1, 1)     

    def get_elevation(self):
        # return array of elevation values at cell centers
        return self.domain.quantities["elevation"].centroid_values

    def get_velocity_and_depth(self):
        # get time-varying quantities available in anuga and convert to velocity and depth
        stage = self.domain.quantities["stage"].centroid_values
        depth = stage - self.get_elevation()
        xmom = self.domain.quantities["xmomentum"].centroid_values
        ymom = self.domain.quantities["ymomentum"].centroid_values
        with np.errstate(divide="ignore", invalid="ignore"):
            h_lim = self.domain.get_minimum_storable_height()
            xvel = np.where(depth < h_lim, 0., xmom/depth)
            yvel = np.where(depth < h_lim, 0., ymom/depth)
        velocity = np.sqrt(xvel**2 + yvel**2)
        return velocity, depth

    def get_vegetation(self):
        # pull directly from Baptist operator
        stemdensity = self.Baptist.veg_density.centroid_values
        stemdiameter = self.Baptist.veg_diameter.centroid_values
        stemheight = self.Baptist.veg_height.centroid_values
        return stemdensity, stemdiameter, stemheight
    
    def set_vegetation(self, stemdensity, stemdiameter, stemheight):
        self.Baptist.set_vegetation(veg_diameter=stemdiameter,
                                    veg_density=stemdensity,
                                    veg_height=stemheight)

    def check_simulation_inputs(self, simstate):
        # Nothing implemented yet
        pass

    def is_parallel(self):
        return True if numprocs > 1 else False
    
    def merge_parallel_veg(self, OutputManager):
        """
        Merge vegetation output files from parallel ANUGA runs into single files per cohort.
        
        This involves reading local-to-global index mapping from each ANUGA subdomain's 
        ``sww`` file, then using that mapping to merge local vegetation arrays into global 
        arrays.
        
        Finally, save merged arrays using :class:`~dycove.sim.outputs.OutputManager` and 
        delete local files.
        """
        outputdir = OutputManager.veg_dir

        sww_name = self.domain.get_name()  # this is the specific name for myid 0, so it has suffix "_0"
        base_name = sww_name[:-2]  # remove the suffix for myid 0
        sww_dir = Path(self.domain.get_datadir())

        # put sww variables from each processor into list
        sww_file_list = [f"{sww_dir}/{base_name}_{p}.sww" for p in range(numprocs)]

        # get subdomain index mapping arrays
        tri_l2g, tri_full_flag = [], []
        for sww_file in sww_file_list:
            with Dataset(sww_file) as sww:
                tri_l2g.append(sww["tri_l2g"][:])
                tri_full_flag.append(sww["tri_full_flag"][:])

        # get no. global grid cells -- this must happen prior to loops b/c we need to know n_global beforehand
        # get_cell_count() method would return only local processor cell count
        n_global = max([tri.max() for tri in tri_l2g]) + 1

        # get indices for full (non-ghost) cells
        f_ids = [np.where(tri_full_flag[p]==1)[0] for p in range(numprocs)]
        f_gids = [tri_l2g[p][f_ids[p]] for p in range(numprocs)]

        # loop over vegetation cohorts (we want one file per cohort
        for cohort_id, cohort in enumerate(self.veg.cohorts, start=1):
            
            # all files for this cohort
            c_files = sorted(outputdir.glob(f"cohort{cohort_id}_proc*.npz"),
                             key=lambda f: int(re.search(r"_(\d+)\.npz", f.name).group(1)))
            
            if len(c_files) == 0:
                msg = (f"Warning: No output files found for cohort {cohort_id}, skipping merge. "
                        "This is likely due to the simulation ending in the middle of an ETS, "
                        "where the last cohort has not been written to an output file yet.")
                r.report(msg, level="WARNING")
                continue

            # first ETS for this cohort, from first filename
            ets_0 = int(Path(c_files[0]).stem.split("_")[-1])

            # number of ecological timesteps experienced by this cohort
            n_ets = len(c_files)//numprocs

            for ets in range(ets_0, ets_0 + n_ets):
                # for storing merged arrays
                merged = {}
                
                for p in range(numprocs):
                    # get local cohort output file (don't use previous file list because proc ID may be unsorted)
                    fname = outputdir / f"cohort{cohort_id}_proc{p}_{ets}.npz"

                    # load numpy and convert to dictionary
                    cohort_local = dict(np.load(fname, allow_pickle=True))

                    # loop over cohort quantities that are array-like
                    for key, values in cohort_local.items():
                        if np.ndim(values) > 0:
                            if key not in merged:  # allocate global array only once
                                merged[key] = np.zeros(n_global, dtype=float)

                            # merge into global array using local-to-global index mapping
                            merged[key][f_gids[p]] = values[f_ids[p]]

                        # for scalar quantities, just pass the output (only need to do it once since it doesn not vary spatially)
                        elif p == 0:
                            merged[key] = values

                    # remove local file
                    fname.unlink()

                # save merged output (basically just remove the "proc" part of the filename
                OutputManager.save_binary(outputdir, f"cohort{cohort_id}", merged, i=ets)

