###############################################################
#  base.py
###############################################################

from abc import ABC, abstractmethod
import numpy as np

from dycove.sim.coupler import VegetationCoupler
from dycove.sim.simulation_data import SimulationTimeState, HydrodynamicStats
from dycove.sim.outputs import OutputManager
from dycove.utils.log import Reporter
from dycove.utils.simulation_reporting import print_model_time_info, print_runtime_updates


r = Reporter()

class HydroSimulationBase(ABC):
    """
    Base class for running hydrodynamic simulations.

    This class is model-agnostic: it is a high-level wrapper that manages the time 
    stepping logic, both for the hydrodynamic and vegetation components.
    
    Model-specific computations are delegated to the appropriate 
    :class:`~dycove.sim.base.HydroEngineBase`.
    Tracking of simulation time state variables is delegated to the 
    :class:`~dycove.sim.simulation_data.SimulationTimeState` dataclass.
    Calculation of hydrodynamic statistics relevant for vegetation processes is 
    handled by the :class:`~dycove.sim.simulation_data.HydrodynamicStats` dataclass.
    Coupling logic between vegetation and hydrodynamic models is managed by 
    :class:`~dycove.sim.coupler.VegetationCoupler`.
    Writing of output files is performed by the 
    :class:`~dycove.sim.outputs.OutputManager` class.
    """

    def __init__(self, engine):
        # Initialize the simulation with a hydrodynamic engine
        self.engine = engine
        # Initialize coupler and output manager
        self.veg_coupler = VegetationCoupler(engine) if engine.veg else None
        self.outputs = OutputManager(engine)

    def run_simulation(self, 
                       sim_time,
                       sim_time_unit='eco-morphodynamic years',
                       n_ets=14, 
                       veg_interval=43200, 
                       hydro_interval=900, 
                       save_interval=3600,
                       ecofac=None,
                       fl_dr=0.15):
        """
        Run a full simulation over the specified time period.

        This method handles the main simulation loop, including:
        
        - Initializing the engine
        - Creating :class:`~dycove.sim.simulation_data.SimulationTimeState` and
          :class:`~dycove.sim.simulation_data.HydrodynamicStats` objects
        - Running vegetation and hydrodynamic loops
        - Saving outputs

        Parameters
        ----------
        sim_time : float
            Simulation duration in units of ``sim_time_unit``.
        sim_time_unit : str, optional
            Either ``'hydrodynamic days'`` or ``'eco-morphodynamic years'``.
            Determines interpretation of ``sim_time``.
        n_ets : int, optional
            Number of ecological time steps per vegetation year. Default is 14.
        veg_interval : int, optional
            Time interval in seconds between vegetation updates.
            Default is 43200 (12 hours).
        hydro_interval : int, optional
            Time interval in seconds between hydrodynamic substeps.
            Default is 900 (15 minutes).
        save_interval : int, optional
            Interval in seconds for writing outputs. Default is 3600 (1 hour).
        ecofac : int, optional
            Ecological (vegetation) acceleration factor relative to hydrodynamics.
            Must equal ``MORFAC`` if Delft3D morphology is enabled. If ``None``,
            it is computed automatically.
        fl_dr : float, optional
            Wet/dry threshold in meters. Cells below this depth are considered dry
            to avoid thin-film instability issues.
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
                                            ecofac=ecofac,
                                            refdate=self.engine.get_refdate()
                                            )

        # define object to hold all hydrodynamic statistics relevant for vegetation processes
        self.hydrostats = HydrodynamicStats(fl_dr=fl_dr)

        # perform some checks on simulation inputs
        self.engine.check_simulation_inputs(self.simstate)
        print_model_time_info(self.simstate)

        # loop over all vegetation time steps
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
        """
        Advance hydrodynamics over all substeps within the current vegetation step.

        Computes min/max water depths, maximum velocities, and flooding
        statistics (HydrodynamicStats) required for vegetation updates.
        """
                
        # add empty placeholders for hydro stats like hmin, vmax, etc
        self.hydrostats.reset(self.engine.get_cell_count())
        # get bed level before hydro loop 
        # (TODO: morphodynamic simulations only, right now is irrelevent if mor=0)
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


class HydroEngineBase(ABC):
    """
    Abstract interface for hydrodynamic engines.
    
    All hydrodynamic models must implement this interface to be compatible with the base 
    hydrodynamic simulation class.
    """

    @abstractmethod
    def initialize(self):
        """ Prepare engine to run (allocate memory, print start, etc.). """
        pass

    @abstractmethod
    def step(self, seconds: int):
        """
        Advance the hydrodynamic model forward by a given number of seconds.
        
        ANUGA implements this by looping with domain.evolve(), while DFM calls dimr.update().
        """
        pass

    @abstractmethod
    def cleanup(self):
        """ Clean up resources. """
        pass

    @abstractmethod
    def get_rank(self):
        """ Get current parallel processor/rank. """
        pass

    @abstractmethod
    def get_cell_count(self) -> int:
        """ Return number of active grid cells in the model. """
        pass

    @abstractmethod
    def get_elevation(self) -> np.ndarray:
        """ Return static bed elevation array (cell-centered). """
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
        """ Return stem density, diameter, and height arrays. """
        pass

    @abstractmethod
    def set_vegetation(self, stemdens, stemdiam, stemheight):
        """ Push vegetation arrays back into the hydro model. """
        pass

    @abstractmethod
    def check_simulation_inputs(self, simstate):
        """ Perform checks on simulation time inputs. """
        pass

    @abstractmethod
    def is_parallel(self):
        """ Return True if model is running in parallel based on number of active processors. """
        pass

    @abstractmethod
    def merge_parallel_veg(self, OutputManager):
        """ Merge vegetation output files across MPI subdomains into single files. """
        pass