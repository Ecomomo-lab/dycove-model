import numpy as np
from typing import Optional
from dataclasses import dataclass, field


@dataclass
class VegetationAttributes:
    """
    Container for vegetation trait parameters and life-cycle settings.

    This class stores all vegetation input attributes and computes
    linear growth rates for height, stem diameter, and root length
    across life stages. Values are typically read from the vegetation
    configuration text file used by the eco-morphodynamic model.
    """

    age_max: int  # max age
    seed_dispersal: int  # amount of ets seed dispersal
    col_method: int  # colonisation method (1 = on bare substrate between max and min water levels, 2 = on bare substrate with mud content
    growth_method: int  # growth method (2 = linear seasonal growth within year; e.g. saltmarsh)
    veg_formula: int  # vegetation formula (is fixed with Baptiste 2007: 154, for different formula other parameters are required)
    nls: int  # number of life stages 
    rootlength_0: float  # initial root length in m
    shootlength_0: float  # initial shoot length in m
    stemdiam_0: float  # initial stem diameter in m
    veg_type: int  # vegetation type (1 = saltmarsh)
    ### TODO: implement growth strategy 1
    start_growth_ets: int  # ecological timestep start growth shoot (in case of growth method 2)
    end_growth_ets: int  # ecological timestep end growth shoot (in case of growth method 2) 
    winter_ets: int  # ecological timestep  step start winter period
    start_col_ets: int  # ecological timestep at which colonisation starts 
    end_col_ets: int  # ecological timestep of last colonisation. Right now, is non-inclusive. TODO: decide on this
    # each attribute below is a list with length equal to number of life stages (nls)
    ht_max: list[float]  # max plant height growth [m]
    diam_max: list[float]  # max stem diameter at ets "end growth shoot"
    rootlength_max: list[float]  # max root length [m] at ets "end growth root"
    years_max: list[int]  # max number of years in life-stage
    stemdens: list[float]  # number of stems per m2 
    fraction_0: float  # initial colonization fraction (0-1)... values are redundant, TODO: move this to gen_veg_attr
    drag: list[float]  # drag coefficient
    mud_percent: list[float]  # mud percentage for colonization
    desic_no_mort: list[float]  # dry fraction below which there is no mortality (zero turns dessication mortality off)
    desic_all_mort: list[float]  # dry fraction above which there total mortality
    flood_no_mort: list[float]  # flood fraction below which there is no mortality (zero turns flooding mortality off)
    flood_all_mort: list[float]  # flood fraction above which there total no mortality
    uproot_no_mort: list[float]  # flow velocity below which there is no mortality (zero turns uprooting mortality off)
    uproot_all_mort: list[float]  # flow velocity above which there is total mortality
    ht_winter_max: list[float]  # max height during winter time
    ### TODO: implement salinity tolerance
    # computed below, after reading in above attributes
    ht_growth_rates: list[float] = field(init=False)
    diam_growth_rates: list[float] = field(init=False)
    root_growth_rates: list[float] = field(init=False)

    def __post_init__(self):
        """Initialize computed growth rate attributes"""
        self.compute_growth_rates()

    def compute_growth_rates(self):
        """
        Compute seasonal linear growth rates for each life stage.

        Stage 0 grows from initial colonization values.
        Later stages grow from previous winter values.
        """
        self.ht_growth_rates = []
        self.diam_growth_rates = []
        self.root_growth_rates = []
        
        for n in range(self.nls):
            if n == 0:
                rates = self._compute_stage_0_rates()
            else:
                rates = self._compute_stage_n_rates(n)
            
            self.ht_growth_rates.append(rates['height'])
            self.diam_growth_rates.append(rates['diameter'])
            self.root_growth_rates.append(rates['root'])

    def _compute_stage_0_rates(self) -> dict:
        """
        Compute growth rates for first life stage starting from
        colonization geometry (height, diameter, roots).
        """
        # TODO: verify we want all 3 growth rates based on start_growth_ets ,
        #       which was technically described initially as being for shoot growth.
        #       I think it makes sense: shoot growth -> diameter/root growth.
        return {
            'height': (self.ht_max[0] - self.shootlength_0) / (self.end_growth_ets - self.start_growth_ets),
            'diameter': (self.diam_max[0] - self.stemdiam_0) / (self.winter_ets - self.start_growth_ets) / self.years_max[0],
            'root': (self.rootlength_max[0] - self.rootlength_0) / (self.winter_ets - self.start_growth_ets) / self.years_max[0]
        }

    def _compute_stage_n_rates(self, n: int) -> dict:
        """
        Compute growth rates for life stage `n` using previous
        winter geometry as starting point.
        """
        return {
            'height': (self.ht_max[n] - self.ht_winter_max[n-1]) / (self.end_growth_ets - self.start_growth_ets),
            'diameter': (self.diam_max[n] - self.diam_max[n-1]) / (self.winter_ets - self.start_growth_ets) / self.years_max[n],
            'root': (self.rootlength_max[n] - self.rootlength_max[n-1]) / (self.winter_ets - self.start_growth_ets) / self.years_max[n]
        }
    

@dataclass
class VegCohort:
    """
    Represents a vegetation cohort (single colonization event) tracked through time.

    A cohort stores plant geometry, density, life-stage progression, and
    mortality contributions across the model domain.

    Attributes
    ----------
    fraction : ndarray
        Vegetation fractional cover per cell (0–1).
    density : float
        Stem density (stems/m²).
    diameter, height, rootlength : float
        Plant geometry (m).
    lifestage : int
        Current life-stage index (0 .. nls-1).
    lifestage_year : int
        Years elapsed in the current life stage.

    potential_mort_* : ndarray
        Mortality potential based solely on environmental stress.
        Represents vulnerability prior to application to vegetated fraction.

    applied_mort_* : ndarray
        Mortality applied to vegetation fraction.
        Represents actual loss per cell.
    """

    fraction: np.ndarray  # array, length = n_cells
    density: float        # stems/m²
    diameter: float       # m
    height: float         # m
    rootlength: float     # m
    lifestage: int        # current life stage index (1..nls)
    lifestage_year: int   # years spent in this stage

    # tracking mortality causes; potential = between 0 and 1, based on stressor not veg fraction
    potential_mort_flood: Optional[np.ndarray] = None
    potential_mort_desic: Optional[np.ndarray] = None
    potential_mort_uproot: Optional[np.ndarray] = None
    potential_mort_burial: Optional[np.ndarray] = None
    potential_mort_scour: Optional[np.ndarray] = None

    # tracking mortality causes; applied = (veg fraction) * (potential mortality)
    applied_mort_flood: Optional[np.ndarray] = None
    applied_mort_desic: Optional[np.ndarray] = None
    applied_mort_uproot: Optional[np.ndarray] = None
    applied_mort_burial: Optional[np.ndarray] = None
    applied_mort_scour: Optional[np.ndarray] = None
    applied_mort_total: Optional[np.ndarray] = None
