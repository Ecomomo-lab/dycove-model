###############################################################
#  base.py
###############################################################

from abc import ABC, abstractmethod
import numpy as np

from dycove.sim.coupler import VegetationCoupler
from dycove.sim.simulation_data import SimulationTimeState, HydrodynamicStats
from dycove.sim.outputs import OutputManager
from dycove.utils.simulation_reporting import print_model_time_info, print_runtime_updates
from dycove.utils.log import Reporter

r = Reporter()

"""
A base class for creating and running hydrodynamic simulations. The class is model-agnostic; all model-specific
methods are called via the respective model "engine" (for which the base class is also included in this file)
"""
class HydroSimulationBase(ABC):

    def __init__(self, engine):
        self.engine = engine
        self.veg_coupler = VegetationCoupler(engine) if engine.veg else None
        self.outputs = OutputManager(engine)

    def run_simulation(self, 
                       sim_time,
                       sim_time_unit='eco-morphodynamic years',
                       n_ets=14, 
                       veg_interval=43200, 
                       hydro_interval=900, 
                       save_interval=3600,
                       vegfac=None,
                       fl_dr=0.15):
        """
        sim_time      : how long to run the simulation, in days if sim_time_units='hydrodynamic days', or years if 
                        sim_time_units='eco-morphodynamic years'
        sim_time_unit : 'hydrodynamic days' or 'eco-morphodynamic years', where eco time is equal to hydro time times morfac/vegfac.
                        Note that if morphology is turned on, vegfac must equal morfac (there is a check for this).
        n_ets         : no. eco timesteps per veg year
        veg_interval  : seconds between eco timesteps. Default is one half day, ~ 1 tidal cycle
        hydro_interval: seconds between internal hydrodynamic timesteps, needed for capturing representative water levels 
                        within eco timesteps so the min water level per eco timestep can be calculated. Also for 
                        calculating average velocities at every hydro_interval, the maximum of which is used for uprooting
                        calculations. So the interval needs to be small enough so that high velocities can be detected.
        save_interval : interval in seconds for writing output to files. This value is not used in DFM because the map
                        output parameter in the MDU file is used by the model internally
        vegfac        : vegetation time factor, i.e. how much faster vegetation processes occur relative to hydrodynamics.
                        Equivalent to morfac when morphology is turned on. Default is None, but will be set to morfac if morfac
                        exists, and if not, will be computed from veg_interval and n_ets inputs.
        fl_dr         : wet/dry threshold (meters) to avoid issue where cells maintain thin film and never dry
        """

        # initialize model
        self.engine.initialize()

        # If engine has save_interval attribute (ANUGA), pass it here
        if hasattr(self.engine, "save_interval"):
            self.engine.save_interval = save_interval

        # define object to hold all static and dynamic variables related to the simulation time frame
        self.simstate = SimulationTimeState(eco_year=1, ets=1, sim_time=sim_time, sim_time_unit=sim_time_unit,
                                            n_ets=n_ets, veg_interval=veg_interval, hydro_interval=hydro_interval,
                                            morfac=int(self.engine.morph_vars["MorFac"]) if self.engine.morphology else None,
                                            vegfac=vegfac,
                                            refdate=self.engine.get_refdate(),
                                            )

        # define object to hold all hydrodynamic statistics relevant for vegetation processes
        self.hydrostats = HydrodynamicStats(fl_dr=fl_dr)

        # perform some checks on simulation inputs
        self.engine.check_simulation_inputs(self.simstate)
        print_model_time_info(self.simstate)

        # loop over vegetation time steps
        for vts in range(self.simstate.n_veg_steps):
            # print runtime stats to the screen
            print_runtime_updates(self.simstate, vts)
            if self.engine.veg is None:
                self.step_hydrodynamics(self.simstate.veg_interval)
            else:
                # get flooding, drying, velocity, bed level changes to be used for colonization and mortality calcs
                self.loop_hydrodynamics()
                # update vegetation and inject those updates back into hydrodynamic model
                self.veg_coupler.update(self.simstate, self.hydrostats)
                # write vegetation variables to output
                self.outputs.save_vegetation_step(self.simstate.eco_year, self.simstate.ets)
        
        r.report("Merging outputs, cleaning up, and finalizing simulation...")
        self.outputs.reconcile_vegetation_output()
        self.engine.cleanup()
        r.report("Simulation complete!")

    def loop_hydrodynamics(self):
        """Advance hydrodynamics by a number of sub-iterations to get finer-scale hydrodynamic statistics"""
        # add empty placeholders for hydro stats like hmin, vmax, etc
        self.hydrostats.reset(self.engine.get_cell_count())
        # get bed level before hydro loop (TODO: morphodynamic simulations only, right now is irrelevent if mor=0)
        self.hydrostats.bedlevel_0 = self.engine.get_elevation()

        # do inner loop of smaller hydrodynamic intervals
        for hts in range(self.simstate.n_hydro_steps):
            # run hydrodynamic model for specified time interval
            self.step_hydrodynamics(self.simstate.hydro_interval)
            # get mean velocity and depth using model-specific methods
            velocity, depth = self.engine.get_velocity_and_depth()
            # update hydrostats counter
            self.hydrostats.update(velocity, depth)

        # get bed level changes (TODO: morphodynamic simulations only, right now is irrelevent if mor=0)
        self.hydrostats.bedlevel_f = self.engine.get_elevation()

    def step_hydrodynamics(self, seconds):
        """Advance hydrodynamics by a specified interval"""
        # By wrapping ANUGA's domain.evolve loop inside the engine.step method, step_hydrodynamics is model-agnostic
        self.engine.step(seconds)
        self.simstate.advance_time(seconds)


"""Abstract interface for hydrodynamic engines."""
class HydroEngineBase(ABC):
    #DAYS_PER_YEAR = 350

    @abstractmethod
    def initialize(self):
        """Prepare engine to run (allocate memory, print start, etc.)."""
        pass

    @abstractmethod
    def step(self, seconds: int):
        """
        Advance the hydrodynamic model forward by a given number of seconds.
        ANUGA might implement this by looping internally with evolve(),
        while DFM can call dimr.update().
        """
        pass

    @abstractmethod
    def cleanup(self):
        """Clean up resources."""
        pass

    @abstractmethod
    def get_rank(self):
        """Get current parallel processor/rank."""
        pass

    @abstractmethod
    def get_cell_count(self) -> int:
        """Return number of active grid cells in the model."""
        pass

    @abstractmethod
    def get_elevation(self) -> np.ndarray:
        """Return static bed elevation array (cell-centered)."""
        pass

    @abstractmethod
    def get_velocity_and_depth(self) -> tuple[np.ndarray, np.ndarray]:
        """
        Return velocity magnitude and depth arrays for each cell.
        If applicable, velocity and depth arrays must be trimmed to active cells.
        """
        pass

    @abstractmethod
    def get_vegetation(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return stem density, diameter, and height arrays."""
        pass

    @abstractmethod
    def set_vegetation(self, stemdens, stemdiam, stemheight):
        """Push vegetation arrays back into the hydro model."""
        pass

    @abstractmethod
    def check_simulation_inputs(self, simstate):
        """Perform checks on simulation time inputs"""
        pass

    @abstractmethod
    def is_parallel(self):
        """Return True if model is running in parallel based on number of active processors"""
        pass

    @abstractmethod
    def merge_parallel_veg(self, OutputManager):
        """Merge vegetation output files across MPI subdomains into single files."""
        pass