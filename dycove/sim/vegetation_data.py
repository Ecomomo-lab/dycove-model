###############################################################
#  vegetation_data.py
###############################################################

import numpy as np
from dataclasses import dataclass, field, fields
from typing import List, get_type_hints


@dataclass
class VegetationAttributes:
    """
    Container for vegetation trait parameters and life-cycle settings.

    This class stores all vegetation input attributes and computes
    linear growth rates for height, stem diameter, and root length
    across life stages. Values are read from a vegetation configuration 
    JSON file (e.g., `veg1.json`). The JSON file must have a number of 
    dictionaries of life-stage attributes equal to the number of life 
    stages indicated by ``nls`` in the general parameters. Note that 
    the type hints for life-stage attributes are lists, because 
    :meth:`~dycove.sim.vegetation.VegetationSpecies.load_vegetation_attributes`
    converts them from a series of dictionaries to one list for each
    attribute.

    Parameters (General)
    --------------------
    age_max : int
        Maximum plant age (in years).
    nls : int
        Number of life stages.
    fraction_0 : float
        Initial colonization fraction (0–1).
    stemht_0 : float
        Initial stem height (m).
    rootlength_0 : float
        Initial root length (m).
    stemdiam_0 : float
        Initial stem diameter (m).
    drag : float
        Drag coefficient C_d.
    start_growth_ets : int
        Ecological timestep at which shoot growth starts.
    end_growth_ets : int
        Ecological timestep at which shoot growth ends.
    winter_ets : int
        Ecological timestep marking the start of the winter period.
    start_col_ets : int
        Ecological timestep at which colonization starts.
    end_col_ets : int
        Ecological timestep at which colonization ends (for the year; currently non-inclusive).
    col_method : int
        Colonization method (1, 2, 3, or 4, see Documentation for details)

    Parameters (Life Stage)
    -----------------------
    stemht_max : list[float]
        Maximum stem height per life stage (m).
    rootlength_max : list[float]
        Maximum root length per life stage (m).
    stemdiam_max : list[float]
        Maximum stem diameter per life stage (m).
    years_max : list[int]
        Maximum number of years spent in each life stage.
    stemdens : list[float]
        Stem density per life stage (number of stems per m²).
    desic_no_mort : list[float]
        Dry fraction below which there is no desiccation mortality.
        (Zero disables desiccation mortality.)
    desic_all_mort : list[float]
        Dry fraction above which there is total desiccation mortality.
    flood_no_mort : list[float]
        Flood fraction below which there is no flooding mortality.
        (Zero disables flooding mortality.)
    flood_all_mort : list[float]
        Flood fraction above which there is total flooding mortality.
    uproot_no_mort : list[float]
        Flow velocity below which there is no uprooting mortality.
        (Zero disables uprooting mortality.)
    uproot_all_mort : list[float]
        Flow velocity above which there is total uprooting mortality.
    stemht_winter_max : list[float]
        Maximum stem height during winter for each life stage.

    Attributes (Computed)
    ---------------------
    ht_growth_rates : list[float]
        Linear growth rate of plant height per life stage (computed in :meth:`__post_init__`).
    diam_growth_rates : list[float]
        Linear growth rate of stem diameter per life stage (computed in :meth:`__post_init__`).
    root_growth_rates : list[float]
        Linear growth rate of root length per life stage (computed in :meth:`__post_init__`).
    """

    age_max: int
    nls: int
    fraction_0: float
    stemht_0: float
    rootlength_0: float
    stemdiam_0: float
    drag: float
    start_growth_ets: int
    end_growth_ets: int
    winter_ets: int
    start_col_ets: int
    end_col_ets: int
    col_method: int
    
    stemht_max: list[float]
    rootlength_max: list[float]
    stemdiam_max: list[float]
    years_max: list[int]
    stemdens: list[float]
    desic_no_mort: list[float]
    desic_all_mort: list[float]
    flood_no_mort: list[float]
    flood_all_mort: list[float]
    uproot_no_mort: list[float]
    uproot_all_mort: list[float]
    stemht_winter_max: list[float]
    
    ht_growth_rates: list[float] = field(init=False)
    diam_growth_rates: list[float] = field(init=False)
    root_growth_rates: list[float] = field(init=False)


    def __post_init__(self):
        # Initialize computed growth rate attributes
        self.compute_growth_rates()

        # Check validity of inputs
        self.validate_input_list_lengths()
        self.validate_mortality_limits()


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
        # TODO: verify we want all 3 growth rates based on start_growth_ets,
        #       which was technically described initially as being for shoot growth.
        #       I think it makes sense: shoot growth -> diameter/root growth.
        return {
            'height': (self.stemht_max[0] - self.stemht_0) / (self.end_growth_ets - self.start_growth_ets),
            'diameter': (self.stemdiam_max[0] - self.stemdiam_0) / (self.winter_ets - self.start_growth_ets) / self.years_max[0],
            'root': (self.rootlength_max[0] - self.rootlength_0) / (self.winter_ets - self.start_growth_ets) / self.years_max[0]
        }


    def _compute_stage_n_rates(self, n: int) -> dict:
        return {
            'height': (self.stemht_max[n] - self.stemht_winter_max[n-1]) / (self.end_growth_ets - self.start_growth_ets),
            'diameter': (self.stemdiam_max[n] - self.stemdiam_max[n-1]) / (self.winter_ets - self.start_growth_ets) / self.years_max[n],
            'root': (self.rootlength_max[n] - self.rootlength_max[n-1]) / (self.winter_ets - self.start_growth_ets) / self.years_max[n]
        }
    

    def validate_input_list_lengths(self):
        """ Check all list fields have length == nls """
        for f in fields(self):
            val = getattr(self, f.name)
            if isinstance(val, list):
                if len(val) != self.nls:
                    msg = (f"One of more life stage attributes have length {len(val)}, but nls == {self.nls}. "
                           "Check that the number of life stage dictionaries in the input .json file equals "
                           "the number of life stages 'nls'.")
                    raise ValueError(msg)
                
    def validate_mortality_limits(self):
        """ Check <stressor>_no_mort < <stressor>_all_mort for each stressor """
        zero_mort_fields = {f.name for f in fields(self) if f.name.endswith('_no_mort')}
        for zero_field in zero_mort_fields:
            all_field = zero_field.replace('_no_mort', '_all_mort')
            zero_vals = getattr(self, zero_field)
            all_vals = getattr(self, all_field)
            for i, (zero, all_) in enumerate(zip(zero_vals, all_vals)):
                if zero != 0 and zero >= all_:
                    msg = (f"'{zero_field}[{i}]' (= {zero}) must be < '{all_field}[{i}]' (= {all_}) "
                        f"or equal to 0 (disabled)")
                    raise ValueError(msg)


@dataclass
class VegCohort:
    """
    Represents a vegetation cohort (single colonization event) tracked through time.

    A cohort stores plant geometry, density, life-stage progression, and
    mortality contributions across the model domain.

    Attributes
    ----------
    name : str
        Name of species.
    fraction : numpy.ndarray
        Vegetation fractional cover per cell (0–1).
    density : float
        Stem density (stems/m²).
    diameter, height, rootlength : float
        Plant geometry (m).
    lifestage : int
        Current life-stage index (0 .. nls-1).
    lifestage_year : int
        Years elapsed in the current life stage.

    potential_mort_* : numpy.ndarray
        Mortality potential based solely on environmental stress.
        Represents vulnerability prior to application to vegetated fraction.

    applied_mort_* : numpy.ndarray
        Mortality applied to vegetation fraction.
        Represents actual loss per cell.
    """

    name: str
    fraction: np.ndarray
    density: float
    diameter: float
    height: float
    rootlength: float
    lifestage: int
    lifestage_year: int

    # Tracking mortality causes; potential = between 0 and 1, based on stressor not veg fraction
    potential_mort_flood: np.ndarray | None = None
    potential_mort_desic: np.ndarray | None = None
    potential_mort_uproot: np.ndarray | None = None
    potential_mort_burial: np.ndarray | None = None
    potential_mort_scour: np.ndarray | None = None

    # Tracking mortality causes; applied = (veg fraction) * (potential mortality)
    applied_mort_flood: np.ndarray | None = None
    applied_mort_desic: np.ndarray | None = None
    applied_mort_uproot: np.ndarray | None = None
    applied_mort_burial: np.ndarray | None = None
    applied_mort_scour: np.ndarray | None = None
    applied_mort_total: np.ndarray | None = None
