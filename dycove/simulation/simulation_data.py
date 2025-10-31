import datetime as dt
import numpy as np
import time
from typing import Optional
from dataclasses import dataclass, field

from dycove.utils.log import Reporter

r = Reporter()


@dataclass
class SimulationTimeState:
    """
    Tracks and manages simulation time state variables.

    This class handles both the hydrodynamic and eco-morphodynamic 
    time scales, including conversions between them.

    Attributes:
        eco_year (int): Current ecological year in the simulation.
        ets (int): Current ecological time step.
        sim_time (float): Simulation time in units given by `sim_time_unit`.
        sim_time_unit (str): Unit of `sim_time`. 
            Must be either 'hydrodynamic days' or 'eco-morphodynamic years'.
        n_ets (int): Number of ecological time steps per year.
        veg_interval (int): Time interval (in seconds) between vegetation updates.
        hydro_interval (int): Time interval (in seconds) between hydrodynamic updates 
            (updates to HydrodynamicStats class).
        morfac (int): Morphological acceleration factor from morphodynamic model.
        vegfac (Optional[int]): Vegetation acceleration factor. If None, it is computed.
        refdate (datetime): Reference start date of the simulation.
        hydrotime_seconds (int): Accumulated hydrodynamic simulation time in seconds.

    Computed Attributes (set in __post_init__):
        days_per_year (float): Implied number of days per ecological year.
        hydro_sim_days (float): Total hydrodynamic days elapsed in simulation.
        veg_sim_years (float): Total eco-morphodynamic years elapsed.
        time_0 (float): Wall-clock timestamp of simulation start.
        times_elapsed (list): List of elapsed wall-clock times at updates.
        n_veg_steps (int): Total number of eco steps in simulation.
        n_hydro_steps (int): Total number of hydrodynamic steps per eco step.
        hydrotime_date (datetime): Current simulation date based on hydrodynamic time.
    """
    
    eco_year: int
    ets: int
    sim_time: float
    sim_time_unit: str  # 'hydrodynamic days' or 'eco-morphodynamic years'
    n_ets: int
    veg_interval: int
    hydro_interval: int
    morfac: int
    vegfac: Optional[int] = None
    refdate: dt.datetime = field(default_factory=lambda: dt.datetime.now())
    hydrotime_seconds: int = 0
    
    # Computed attributes
    days_per_year: float = field(init=False)
    hydro_sim_days: float = field(init=False)
    veg_sim_years: float = field(init=False)
    time_0: float = field(init=False)
    times_elapsed: list = field(init=False)
    n_veg_steps: int = field(init=False)
    n_hydro_steps: int = field(init=False)
    hydrotime_date: dt.datetime = field(init=False)

    def __post_init__(self):
        """Initialize computed attributes and validate inputs"""
        self.time_0 = time.time()
        self.times_elapsed = []
        
        # Compute/validate vegfac first
        self.compute_vegfac()
        
        # Compute days_per_year based on vegfac
        self.compute_days_per_year()
        
        # Use days_per_year to compute time conversions
        self.compute_time_conversions()
        
        # Compute step counts
        self.n_veg_steps = int(self.veg_sim_years * self.n_ets)
        self.n_hydro_steps = int(self.veg_interval / self.hydro_interval)

    @property
    def vegtime_date(self) -> dt.datetime:
        return self.refdate + dt.timedelta(seconds=self.hydrotime_seconds * self.vegfac)

    def advance_time(self, seconds: int):
        self.hydrotime_seconds += seconds

    def update_ets(self):
        if self.ets < self.n_ets:
            self.ets += 1
        else:
            self.ets = 1
            self.eco_year += 1

    def compute_vegfac(self):
        """
        Determine vegetation acceleration factor, 
          prioritizing: model morfac > input value > derived value.
        """
        # Check that input vegfac matches morfac if morphology is on
        if self.vegfac is not None:
            if self.morfac and self.morfac != self.vegfac:
                msg = ("MORFAC mismatch: when morphology is on (.mor file exists), input vegfac must "
                       "either equal the DFM morfac value or be left as None.")
                r.report(msg, level="ERROR")
                raise ValueError(msg)
        else:
            # Use morfac from model if morphology is on
            if self.morfac:
                msg = f"Using morfac = {self.morfac:d} from morphology file as vegfac."
                r.report(msg, level="INFO")
                self.vegfac = self.morfac
            # Otherwise, compute vegfac implied from vegetation coupling times
            else:
                # Compute ideal (continuous) vegfac
                vegfac_est = (86400 * 365.) / (self.veg_interval * self.n_ets)
                # Round to nearest whole number
                self.vegfac = int(round(vegfac_est))
                
                msg = (f"Computed vegfac = {self.vegfac:d}, derived from veg_interval and n_ets "
                       f"(rounded from {vegfac_est:.3f})")
                r.report(msg, level="INFO")
    
    def compute_days_per_year(self):
        """
        Compute the implied 'days per year' from vegfac, veg_interval, and n_ets.
        Validates that the result is reasonable (350-380 days).
        """
        self.days_per_year = (self.vegfac * self.veg_interval * self.n_ets) / 86400.
        
        if 350 <= self.days_per_year <= 380:
            msg = f"vegfac = {self.vegfac:d} corresponds to {self.days_per_year:.2f} days per year."
            r.report(msg, level="INFO")
        else:
            msg = (f"Input parameters (vegfac = {self.vegfac:d}, veg_interval = {self.veg_interval}, "
                   f"and n_ets = {self.n_ets}) correspond to {self.days_per_year:.2f} days per year. "
                   "Please check these three inputs for consistency. "
                   "Note that the vegfac input can be left out (None) to have it computed automatically "
                   "based on veg_interval and n_ets.")
            r.report(msg, level="ERROR")
            raise ValueError(msg)
    
    def compute_time_conversions(self):
        """Convert between hydrodynamic days and eco-morphodynamic years"""
        if self.sim_time_unit == 'hydrodynamic days':
            self.hydro_sim_days = self.sim_time
            self.veg_sim_years = self.sim_time * self.vegfac / self.days_per_year
        elif self.sim_time_unit == 'eco-morphodynamic years':
            self.veg_sim_years = self.sim_time
            self.hydro_sim_days = self.sim_time * self.days_per_year / self.vegfac
        else:
            msg = f"Invalid sim_time_unit: '{self.sim_time_unit}'. Must be 'hydrodynamic days' or 'eco-morphodynamic years'."
            r.report(msg, level="ERROR")
            raise ValueError(msg)


@dataclass
class HydrodynamicStats:
    """
    Stores hydrodynamic model outputs used to compute vegetation responses.

    Tracks water depths, velocities, and flood statistics across the computational grid.

    Attributes:
        fl_dr (float): Flooding threshold (depth above which a cell is considered flooded).
        h_min (np.ndarray): Minimum water depth observed in each cell during hydro_interval.
        h_max (np.ndarray): Maximum water depth observed in each cell during hydro_interval.
        v_max (np.ndarray): Maximum velocity observed in each cell during hydro_interval.
        flood_counts (np.ndarray): Number of timesteps each cell was flooded during hydro_interval.
        bedlevel_0 (np.ndarray): Bed elevations per cell before hydro_interval.
        bedlevel_f (np.ndarray): Bed elevations per cell after hydro_interval.

    Properties:
        bedlevel_diff (np.ndarray): Change in bed elevation (bedlevel_f - bedlevel_0).

    Methods:
        flood_frac(n_substeps): Fraction of time each cell was flooded over `n_substeps`.
        dry_frac(n_substeps): Fraction of time each cell was dry over `n_substeps`.
        reset(n_cells): Initialize/reset all arrays for `n_cells` grid cells.
        update(vel, depth): Update min/max values and flood counts from new hydrodynamic step.
    """
    
    fl_dr: float
    h_min: Optional[np.ndarray] = None
    h_max: Optional[np.ndarray] = None
    v_max: Optional[np.ndarray] = None
    flood_counts: Optional[np.ndarray] = None
    bedlevel_0: Optional[np.ndarray] = None
    bedlevel_f: Optional[np.ndarray] = None

    @property
    def bedlevel_diff(self) -> np.ndarray:
        return self.bedlevel_f - self.bedlevel_0
    
    def flood_frac(self, n_substeps: int) -> np.ndarray:
        # fraction of time cells were flooded over ets
        return self.flood_counts/n_substeps
    
    def dry_frac(self, n_substeps: int) -> np.ndarray:
        # fraction of time cells were dry over ets
        return (n_substeps - self.flood_counts)/n_substeps
    
    def reset(self, n_cells: int):
        self.h_min = np.full(n_cells, np.inf)
        self.h_max = np.full(n_cells, -np.inf)
        self.v_max = np.full(n_cells, -np.inf)
        self.flood_counts = np.zeros(n_cells)     

    def update(self, vel, depth):
        # update arrays of min/max hydro variables
        np.maximum(self.v_max, vel, out=self.v_max)
        np.minimum(self.h_min, depth, out=self.h_min)
        np.maximum(self.h_max, depth, out=self.h_max)
        # add to count of flooded/dry cells
        self.flood_counts[depth >= self.fl_dr] += 1
