###############################################################
#  base.py
###############################################################

from abc import ABC, abstractmethod
import numpy as np

from dycove.sim.coupler import VegetationCoupler
from dycove.sim.simulation_data import SimulationTimeState, HydrodynamicStats
from dycove.sim.outputs import OutputManager
from dycove.utils.simulation_reporting import Reporter

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
        self.engine = engine
        self.veg_active = True if self.engine.veg else False
        self.veg_coupler = VegetationCoupler(engine) if self.veg_active else None
        self.outputs = OutputManager(engine)

    def run_simulation(self, 
                       sim_time,
                       sim_time_unit,
                       n_ets=14, 
                       veg_interval=43200, 
                       hydro_interval=900, 
                       save_interval=3600,
                       ecofac=None
                       ):
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
        sim_time_unit : str
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
        """

        eco_args = {"n_ets": (n_ets, 14),  # actual vs. default arguments
                    "veg_interval": (veg_interval, 43200),
                    "ecofac": (ecofac, None)
                    }
        self.verify_eco_args(eco_args)       

        self.engine.initialize()

        # If engine has save_interval attribute (ANUGA), pass it here
        if hasattr(self.engine, "save_interval"):
            self.engine.save_interval = save_interval

        # Define object to hold all static and dynamic variables related to the simulation time frame
        self.simstate = SimulationTimeState(eco_year=1, ets=1, sim_time=sim_time, sim_time_unit=sim_time_unit,
                                            n_ets=n_ets, veg_interval=veg_interval, hydro_interval=hydro_interval,
                                            veg_active=self.veg_active,
                                            morfac=int(self.engine.morph_vars["MorFac"]) if self.engine.morphology else None,
                                            ecofac=ecofac,
                                            refdate=self.engine.get_refdate()
                                            )

        # Define object to hold all hydrodynamic statistics relevant for vegetation processes
        self.hydrostats = HydrodynamicStats(n_hydro_substeps=self.simstate.n_hydro_substeps,
                                            n_cells=self.engine.get_cell_count()
                                            )

        self.engine.check_simulation_inputs(self.simstate)

        r.print_model_time_info(self.simstate, self.veg_active)

        # Loop over all Ecological Time Steps
        if not self.veg_active:
            self.run_simulation_without_vegetation()
        else:
            self.run_simulation_with_vegetation()


    def verify_eco_args(self, inputs: dict):
        default_inputs = [values[0] != values[1] for name, values in inputs.items()]
        if any(default_inputs) and not self.veg_active:
            msg = ("No vegetation object was specified for this simulation but one or "
                   "vegetation-related inputs to run_simulation() were given non-default "
                   "values. If the intent was to run a non-vegetation simulation, please "
                   "leave vegetation-related inputs as their default values.")
            r.report(msg, level="ERROR")
            raise ValueError(msg)   
        

    def run_simulation_without_vegetation(self):
        for hts in range(self.simstate.n_hydro_steps):
            r.print_runtime_updates(self.simstate, hts, self.veg_active)
            self.hydro_step(self.simstate.veg_interval)  # run a big step
        self.finalize_simulation()


    def run_simulation_with_vegetation(self):
        for vts in range(self.simstate.n_veg_steps):
            r.print_runtime_updates(self.simstate, vts, self.veg_active)
            self.eco_hydro_loop(self.simstate.n_hydro_substeps, self.simstate.hydro_interval)
            self.veg_coupler.update(self.simstate, self.hydrostats)
            self.outputs.save_vegetation_step(self.simstate)
        self.finalize_simulation()


    def eco_hydro_loop(self, n_substeps, interval):
        """
        Advance hydrodynamics over all substeps within the current vegetation step.

        Computes min/max water depths, maximum velocities, and flooding
        statistics (HydrodynamicStats) required for vegetation updates.
        """
        
        # Add empty placeholders for hydro stats like hmin, vmax, etc
        self.hydrostats.reset()

        # Get bed level before hydro loop (right now is irrelevent if mor=0)
        self.hydrostats.bedlevel_0 = self.engine.get_elevation()

        # Run inner loop of smaller hydrodynamic intervals
        for hts in range(n_substeps):
            self.hydro_step(interval)
            velocity, depth = self.engine.get_velocity_and_depth()
            self.hydrostats.update(hts, velocity, depth)

        # Get bed level changes (right now is irrelevent if mor=0)
        self.hydrostats.bedlevel_f = self.engine.get_elevation()


    def hydro_step(self, seconds):
        """ Advance hydrodynamics by a specified interval """
        self.engine.step(seconds)
        self.simstate.advance_time(seconds)


    def finalize_simulation(self):
        r.report("Merging outputs, cleaning up, and finalizing simulation...")
        self.outputs.reconcile_vegetation_output()
        self.engine.cleanup()
        r.report("Simulation complete!")


class HydroEngineBase(ABC):
    """
    Abstract interface for hydrodynamic engines.
    
    All hydrodynamic models must implement this interface to be compatible with the base 
    hydrodynamic simulation class.
    """

    @abstractmethod
    def initialize(self):
        """ Prepare engine to run (allocate memory, read/check inputs, print start, etc.) """
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
    def get_refdate(self):
        """ 
        Get model reference date as a datetime object. 
        Required for now, only really used by DFM. 
        """
        pass
    
    @abstractmethod
    def get_cell_count(self) -> int:
        """ Return number of active grid cells in the model. """
        pass

    @abstractmethod
    def get_elevation(self) -> np.ndarray:
        """ Return array of bed elevations at cell-centers. """
        pass

    @abstractmethod
    def get_velocity_and_depth(self) -> tuple[np.ndarray, np.ndarray]:
        """
        Return velocity magnitude and depth arrays at cell-centers.
        """
        pass

    @abstractmethod
    def get_vegetation(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """ Return stem dens., diam., and height arrays from numerical model. """
        pass

    @abstractmethod
    def set_vegetation(self, stemdens, stemdiam, stemheight):
        """ Push vegetation arrays back into the numerical model. """
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