from pytest import approx, mark
import numpy as np
from dycove.sim.vegetation import VegetationSpecies
from dycove.sim.vegetation_data import VegetationAttributes
from dycove.sim.vegetation_data import VegCohort


def test_colonization(constants, attributes):
    c = constants
    a = attributes

    def mock_attrs():
        class MockAttr:
            fraction_0    = a["fraction_0"]
            start_col_ets = a["start_col_ets"]
            end_col_ets   = a["end_col_ets"]
            stemdens      = [c["m"]]  # list because life stage attribute
            stemdiam_0    = c["D"]
            stemht_0      = c["hv"]
            rootlength_0  = a["rootlength_0"]
        return MockAttr()

    # Create species object without running __init__, and mock relevant attributes
    veg = VegetationSpecies.__new__(VegetationSpecies)
    veg.attrs = mock_attrs()
    veg.name = "veg"
    veg.seed_frac = 1.0
    veg.seed_method = "deterministic"
    veg.cohorts = []
    
    # Define test inundation limits
    min_depths = np.array([1.0, 0.2, 0.1, 0.0])
    max_depths = np.array([1.0, 1.0, 1.0, 0.1])
    fl_dr = 0.15

    # Test one ETS each inside and outside the range of colonization based on input file
    for ets in [1, 2]:
        veg.colonization(ets, min_depths, max_depths, fl_dr)

    # Check only one cohort colonized
    assert len(veg.cohorts) == 1
    # Check colonized fractions
    fractions = veg.cohorts[0].fraction
    expected = np.array([0., 0., a["fraction_0"], 0.])
    assert (fractions == expected).all()


def test_hydrodynamic_mortality(attributes):
    a = attributes

    def mock_attrs():
        class MockAttr:
            nls             = a["nls"]
            flood_no_mort   = [a["flood_no_mort"]]
            flood_all_mort  = [a["flood_all_mort"]]
            desic_no_mort   = [a["desic_no_mort"]]
            desic_all_mort  = [a["desic_all_mort"]]
            uproot_no_mort  = [a["uproot_no_mort"]]
            uproot_all_mort = [a["uproot_all_mort"]]
        return MockAttr()
    
    # Create species object without running __init__, and mock relevant attributes
    veg = VegetationSpecies.__new__(VegetationSpecies)
    veg.attrs = mock_attrs()
    
    # Define wet/dry fractions and max velocity
    vel_max  = np.array([0.8, 1.1, 1.3, 1.5, 1.8])
    fld_frac = np.array([0.1, 0.4, 0.5, 0.6, 0.9])
    #dry_frac = np.array([1 - f for f in fld_frac])
    dry_frac = np.array([0.9, 0.6, 0.5, 0.4, 0.1])

    # Test method
    veg.mortality_hydrodynamic(fld_frac, dry_frac, vel_max)

    # Check against expected values
    mort_fld = veg.potential_mort_flood
    mort_dry = veg.potential_mort_desic
    mort_upr = veg.potential_mort_uproot
    expected_fld = np.array([0., 0., 0.5, 1., 1.])
    expected_dry = np.array([0.5, 0., 0., 0., 0.])
    expected_upr = np.array([0., 0.2, 0.6, 1., 1.])

    assert mort_fld[0] == approx(expected_fld, rel=1e-6)
    assert mort_dry[0] == approx(expected_dry, rel=1e-6)
    assert mort_upr[0] == approx(expected_upr, rel=1e-6)


def test_morphodynamic_mortality(constants, attributes):
    a = attributes
    c = constants

    def mock_cohort():
        class MockCohort:
            height     = c["hv"]
            rootlength = a["rootlength_0"]
            potential_mort_burial = None,
            potential_mort_scour  = None,
        return MockCohort()
    
    # Create species object without running __init__, and mock relevant attributes
    veg = VegetationSpecies.__new__(VegetationSpecies)
    veg.cohorts = [mock_cohort()]
    veg.mor = 1
    
    # Define wet/dry fractions and max velocity
    bl_diff  = np.array([0.7, 0.5, 0.1, -0.1, -0.5])  # positive = deposition, negative = erosion

    # Test method
    veg.mortality_morphodynamic(bl_diff)

    # Check against expected values
    mort_bur = veg.cohorts[0].potential_mort_burial
    mort_sco = veg.cohorts[0].potential_mort_scour
    expected_bur = np.array([1., 1., 0., 0., 0.])
    expected_sco = np.array([0., 0., 0., 1., 1.])

    assert mort_bur == approx(expected_bur, rel=1e-6)
    assert mort_sco == approx(expected_sco, rel=1e-6)


def test_stemheight_growth(constants, attributes):
    a = attributes
    c = constants

    def mock_attrs():
        class MockAttr:
            start_growth_ets = a["start_growth_ets"]
            end_growth_ets = a["end_growth_ets"]
            winter_ets = a["winter_ets"]
            stemht_max = [a["stemht_max"]]
            ht_growth_rates = [a["ht_growth_rates"]]
            stemht_winter_max = [a["stemht_winter_max"]]
        return MockAttr()
    
    def mock_cohort():
        class MockCohort:
            lifestage = 1
            height = c["hv"]
        return MockCohort()

    # Create species object without running __init__, and mock relevant attributes
    veg = VegetationSpecies.__new__(VegetationSpecies)
    veg.attrs = mock_attrs()
    veg.cohorts = [mock_cohort()]

    # Test growth over several ETS
    for ets in range(2, 13):
        veg.stemheight_growth(ets)
        if ets < a["winter_ets"]:
            height = c["hv"] + max(a["ht_growth_rates"]*(ets - a["start_growth_ets"]), 0)
        else:
            height = a["stemht_winter_max"]
        assert veg.cohorts[0].height == min(a["stemht_max"], height)