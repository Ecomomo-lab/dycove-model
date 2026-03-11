import pytest
from types import SimpleNamespace
import numpy as np
from dycove.sim.vegetation import VegetationSpecies
from dycove.sim.vegetation_data import VegetationAttributes
from dycove.sim.vegetation_data import VegCohort


# ------------------------------------------------------------ #
# tests for colonization
# ------------------------------------------------------------ #

class TestColonization:

    @staticmethod
    def mock_attrs(c, a):
        class MockAttr:
            fraction_0    = a["fraction_0"]
            start_col_ets = a["start_col_ets"]
            end_col_ets   = a["end_col_ets"]
            stemdens      = [c["m"]]  # list because life stage attribute
            stemdiam_0    = c["D"]
            stemht_0      = c["hv"]
            rootlength_0  = a["rootlength_0"]
        return MockAttr()
    

    @pytest.mark.unit
    def test_colonization_outside_ets_window(self, constants, attributes):
        veg = VegetationSpecies.__new__(VegetationSpecies)
        veg.attrs = self.mock_attrs(constants, attributes)

        veg.cohorts = []
        veg.name = "test"

        veg.colonization(ets=1,  # just before window
                         min_depths=np.zeros(5),
                         max_depths=np.ones(5),
                         fl_dr=0.15,
                         )
        veg.colonization(ets=4,  # just after window
                         min_depths=np.zeros(5),
                         max_depths=np.ones(5),
                         fl_dr=0.15,
                         )
        
        assert len(veg.cohorts) == 0


    @pytest.mark.unit
    def test_colonization_inside_ets_window(self, constants, attributes):
        veg = VegetationSpecies.__new__(VegetationSpecies)
        veg.attrs = self.mock_attrs(constants, attributes)

        veg.cohorts = []
        veg.name = "test"
        veg.seed_frac = 1.0
        veg.seed_method = "deterministic"

        veg.colonization(ets=2,  # inside window
                         min_depths=np.zeros(5),
                         max_depths=np.ones(5),
                         fl_dr=0.15,
                         )
        veg.colonization(ets=3,  # inside window
                         min_depths=np.zeros(5),
                         max_depths=np.ones(5),
                         fl_dr=0.15,
                         )
        
        assert len(veg.cohorts) == 2


    @pytest.mark.unit
    def test_colonization_mask_logic(self, monkeypatch, constants, attributes):
        veg = VegetationSpecies.__new__(VegetationSpecies)
        veg.attrs = self.mock_attrs(constants, attributes)

        veg.cohorts = []
        veg.name = "test"

        min_depths = np.array([0.0, 0.2, 0.0])
        max_depths = np.array([0.2, 0.2, 0.1])

        # Deterministic seed mask
        monkeypatch.setattr(veg,
                            "create_seed_fraction_mask",
                            lambda n: np.array([True, False, True])
                            )

        captured = {}
        def fake_compute(existing, inds):
            captured["inds"] = inds.copy()
            return np.zeros(len(min_depths))

        monkeypatch.setattr(veg, "compute_new_fraction", fake_compute)

        veg.colonization(ets=2, min_depths=min_depths, max_depths=max_depths, fl_dr=0.15)

        # Only index 0 satisfies dry & flood & seed
        assert np.array_equal(captured["inds"], np.array([True, False, False]))


    @pytest.mark.unit
    def test_colonization_creates_new_cohort(self, monkeypatch, constants, attributes):
        veg = VegetationSpecies.__new__(VegetationSpecies)
        veg.attrs = self.mock_attrs(constants, attributes)

        veg.cohorts = []
        veg.name = "speciesA"

        min_depths = np.array([0.0, 0.2, 0.0])
        max_depths = np.array([0.2, 0.2, 0.1])

        monkeypatch.setattr(veg,
                            "create_seed_fraction_mask",
                            lambda n: np.ones(n, dtype=bool)
                            )

        monkeypatch.setattr(veg,
                            "compute_new_fraction",
                            lambda existing, inds: np.array([0.5, 0.0, 0.0])
                            )

        veg.colonization(ets=2,
                         min_depths=min_depths,
                         max_depths=max_depths,
                         fl_dr=0.15,
                         )

        assert len(veg.cohorts) == 1
        c = veg.cohorts[0]

        assert c.name == "speciesA"
        assert np.allclose(c.fraction, [0.5, 0.0, 0.0])
        assert c.density == constants["m"]
        assert c.lifestage == 1


    @pytest.mark.unit
    def test_colonization_uses_combined_cohorts(self, monkeypatch, constants, attributes):
        veg = VegetationSpecies.__new__(VegetationSpecies)
        veg.attrs = self.mock_attrs(constants, attributes)

        veg.cohorts = ["speciesA_obj"]
        veg.name = "speciesA"

        combined_cohorts = ["speciesA_obj", "speciesB_obj"]

        def fake_compute(existing, inds):
            assert existing is combined_cohorts
            return np.zeros(3)

        monkeypatch.setattr(veg, "compute_new_fraction", fake_compute)
        monkeypatch.setattr(veg, "create_seed_fraction_mask", lambda n: np.ones(n))

        veg.colonization(ets=1,
                         min_depths=np.zeros(3),
                         max_depths=np.ones(3),
                         fl_dr=0.15,
                         combined_cohorts=combined_cohorts,
                         )


# ------------------------------------------------------------ #
# tests for create_seed_fraction_mask()
# ------------------------------------------------------------ #

class TestCreateSeedFractionMask:

    @pytest.mark.unit
    def test_seed_mask_deterministic(self, array_len):
        veg = VegetationSpecies.__new__(VegetationSpecies)
        veg.seed_method = "deterministic"
        veg.seed_frac = 0.3
        
        mask1 = veg.create_seed_fraction_mask(array_len)
        mask2 = veg.create_seed_fraction_mask(array_len)

        assert mask1.dtype == bool
        assert np.sum(mask1) == int(np.floor(0.3 * array_len))
        assert np.array_equal(mask1, mask2)


    @pytest.mark.unit
    def test_seed_mask_random(self, array_len):
        veg = VegetationSpecies.__new__(VegetationSpecies)
        veg.seed_method = "random"
        veg.seed_frac = 0.3
        
        mask = veg.create_seed_fraction_mask(array_len)

        assert mask.shape == (array_len,)
        assert mask.dtype == bool
        assert mask.sum() == int(np.floor(0.3 * array_len))


    @pytest.mark.unit
    def test_seed_mask_zero_fraction(self, array_len):
        veg = VegetationSpecies.__new__(VegetationSpecies)
        veg.seed_method = "random"
        veg.seed_frac = 0.0

        mask = veg.create_seed_fraction_mask(array_len)

        assert not mask.any()
        assert mask.sum() == 0


    @pytest.mark.unit
    def test_seed_mask_full_fraction(self, array_len):
        veg = VegetationSpecies.__new__(VegetationSpecies)
        veg.seed_method = "random"
        veg.seed_frac = 1.0

        mask = veg.create_seed_fraction_mask(array_len)

        assert mask.all()
        assert mask.sum() == array_len


# ------------------------------------------------------------ #
# tests for compute_new_fraction()
# ------------------------------------------------------------ #

class TestComputeNewFraction:

    @staticmethod
    def mock_attrs(a):
        class MockAttr:  # only need this one attribute
            fraction_0 = a["fraction_0"]
        return MockAttr()
    
    @staticmethod
    def mock_cohort(frac):
        class MockCohort:
            fraction = frac
        return MockCohort()
    
    
    @pytest.mark.unit
    def test_compute_new_fraction_no_existing_cohorts(self, attributes):
        veg = VegetationSpecies.__new__(VegetationSpecies)
        veg.attrs = self.mock_attrs(attributes)

        inds2colonize = np.array([True, False, True, False, False])
        existing_cohorts = []

        new_frac = veg.compute_new_fraction(
            existing_cohorts=existing_cohorts,
            inds2colonize=inds2colonize,
        )

        expected = np.zeros(len(inds2colonize))
        expected[inds2colonize] = veg.attrs.fraction_0

        assert np.allclose(new_frac, expected)
        

    @pytest.mark.unit
    def test_compute_new_fraction_with_existing_cohorts(self, attributes):
        veg = VegetationSpecies.__new__(VegetationSpecies)
        veg.attrs = self.mock_attrs(attributes)

        inds2colonize = np.array([True, True, False, True])

        # Existing fractions sum to less than 1 everywhere
        existing_cohorts = [
            self.mock_cohort(np.array([0.2, 0.1, 0.3, 0.4])),
            self.mock_cohort(np.array([0.3, 0.2, 0.1, 0.1])),
        ]

        new_frac = veg.compute_new_fraction(
            existing_cohorts=existing_cohorts,
            inds2colonize=inds2colonize,
        )

        fraction_uncovered = 1 - (
            existing_cohorts[0].fraction + existing_cohorts[1].fraction
        )
        candidate = np.minimum(fraction_uncovered, veg.attrs.fraction_0)

        expected = np.zeros(len(inds2colonize))
        expected[inds2colonize] = candidate[inds2colonize]

        assert np.allclose(new_frac, expected)


    @pytest.mark.unit
    def test_compute_new_fraction_no_available_space(self, attributes):    
        veg = VegetationSpecies.__new__(VegetationSpecies)
        veg.attrs = self.mock_attrs(attributes)

        inds2colonize = np.array([True, True, True])

        # Fractions sum to 1 everywhere
        existing_cohorts = [
            self.mock_cohort(np.array([0.6, 0.5, 0.7])),
            self.mock_cohort(np.array([0.4, 0.5, 0.3])),
        ]

        new_frac = veg.compute_new_fraction(
            existing_cohorts=existing_cohorts,
            inds2colonize=inds2colonize,
        )

        assert np.allclose(new_frac, 0.0)


# ------------------------------------------------------------ #
# tests for linear_mortality_func()
# ------------------------------------------------------------ #

class TestLinearMortalityFunc:

    @pytest.mark.unit
    def test_linear_mortality_mixed_values(self):
        stressor = np.array([0.0, 0.5, 0.75, 1.25, 1.5, 2.0])
        th_min = 0.5
        th_max = 1.5

        result = VegetationSpecies.linear_mortality_func(stressor, th_min, th_max)

        expected = np.array([0.0,   # below th_min
                            0.0,   # exactly th_min
                            0.25,  # linear region (1/4)
                            0.75,  # linear region (3/4)
                            1.0,   # exactly th_max
                            1.0,   # above th_max
                            ])

        assert np.allclose(result, expected)


    @pytest.mark.unit
    def test_linear_mortality_clipping(self):
        stressor = np.array([-10.0, 0.0, 10.0])
        th_min = 0.0
        th_max = 1.0

        result = VegetationSpecies.linear_mortality_func(stressor, th_min, th_max)

        assert np.all(result >= 0.0)
        assert np.all(result <= 1.0)


# ------------------------------------------------------------ #
# tests for mortality_hydrodynamic()
# ------------------------------------------------------------ #

class TestMortalityHydrodynamic:

    @pytest.mark.unit
    def test_mortality_hydro_all_disabled(self):
        veg = VegetationSpecies.__new__(VegetationSpecies)

        veg.attrs = SimpleNamespace(nls=2,
                                    flood_no_mort=[0, 0],
                                    desic_no_mort=[0, 0],
                                    uproot_no_mort=[0, 0],
                                    )

        fld = np.array([0.1, 0.2])  # fraction of time flooded
        dry = np.array([0.9, 0.8])  # inverse of fld
        vel = np.array([0.5, 0.6])  # in m/s

        veg.mortality_hydrodynamic(fld, dry, vel)

        for n in range(veg.attrs.nls):
            assert np.all(veg.potential_mort_flood[n] == 0)
            assert np.all(veg.potential_mort_desic[n] == 0)
            assert np.all(veg.potential_mort_uproot[n] == 0)


    @pytest.mark.unit
    def test_mortality_hydro_calls_linear(self, monkeypatch):
        veg = VegetationSpecies.__new__(VegetationSpecies)

        veg.attrs = SimpleNamespace(nls=1,
                                    flood_no_mort=[0.3],
                                    flood_all_mort=[0.5],
                                    desic_no_mort=[0.8],
                                    desic_all_mort=[0.9],
                                    uproot_no_mort=[1.0],
                                    uproot_all_mort=[1.5],
                                    )

        fld = np.array([0.5])
        dry = np.array([0.5])
        vel = np.array([0.5])

        calls = []

        def fake_linear_func(stressor, th_min, th_max):
            calls.append((th_min, th_max))
            return stressor

        monkeypatch.setattr(veg, "linear_mortality_func", fake_linear_func)

        veg.mortality_hydrodynamic(fld, dry, vel)

        assert len(calls) == 3
        assigned_calls = [(0.3, 0.5), (0.8, 0.9), (1., 1.5)]  # (th_min, th_max) pairs
        assert all(call in assigned_calls for call in calls)


    @pytest.mark.unit
    def test_mortality_hydro_mixed_enable(self, monkeypatch):
        veg = VegetationSpecies.__new__(VegetationSpecies)

        veg.attrs = SimpleNamespace(nls=2,
                                    flood_no_mort=[0, 0.3],
                                    flood_all_mort=[1, 0.5],
                                    desic_no_mort=[0.8, 0],
                                    desic_all_mort=[0.9, 1],
                                    uproot_no_mort=[0, 1.0],
                                    uproot_all_mort=[1, 1.5],
                                    )

        arr = np.array([0.5])

        monkeypatch.setattr(veg, 
                            "linear_mortality_func",
                            lambda *args: np.array([0.7])
                            )

        veg.mortality_hydrodynamic(arr, arr, arr)

        # First life stage
        assert veg.potential_mort_flood[0][0] == 0
        assert veg.potential_mort_desic[0][0] == 0.7
        assert veg.potential_mort_uproot[0][0] == 0
        # Second life stage
        assert veg.potential_mort_flood[1][0] == 0.7
        assert veg.potential_mort_desic[1][0] == 0
        assert veg.potential_mort_uproot[1][0] == 0.7


    @pytest.mark.unit
    def test_mortality_hydro_array_length(self):
        veg = VegetationSpecies.__new__(VegetationSpecies)

        veg.attrs = SimpleNamespace(nls=1,
                                    flood_no_mort=[0.3],
                                    flood_all_mort=[0.5],
                                    desic_no_mort=[0.8],
                                    desic_all_mort=[0.9],
                                    uproot_no_mort=[1.0],
                                    uproot_all_mort=[1.5],
                                    )

        arr = np.array([0.5, 0.5, 0.5, 0.5])

        veg.mortality_hydrodynamic(arr, arr, arr)

        outputs = [veg.potential_mort_flood,
                   veg.potential_mort_desic,
                   veg.potential_mort_uproot,
                   ]
        assert all(len(a) == 1 for a in outputs)  # b/c nls == 1
        assert all(len(a[0]) == len(arr) for a in outputs)  # array shapes must be 1D and match


# ------------------------------------------------------------ #
# tests for mortality_morphodynamic()
# ------------------------------------------------------------ #

class TestMortalityMorphodynamic:
    
    @staticmethod
    def mock_cohort(c, a):
        class MockCohort:
            height     = c["hv"]
            rootlength = a["rootlength_0"]
            fraction   = np.array([0., 0.4, 0.6, 1.0])
        return MockCohort()
    
    @pytest.mark.unit
    def test_mortality_morpho_disabled(self, constants, attributes):
        veg = VegetationSpecies.__new__(VegetationSpecies)

        veg.mor = 0  # disabled

        c = self.mock_cohort(constants, attributes)
        veg.cohorts = [c]

        bedlevel_change = np.array([10, 10, 10, 10])

        veg.mortality_morphodynamic(bedlevel_change)

        assert np.all(c.potential_mort_burial == 0)
        assert np.all(c.potential_mort_scour == 0)


    @pytest.mark.unit
    def test_mortality_morpho_enabled(self, constants, attributes):
        
        veg = VegetationSpecies.__new__(VegetationSpecies)

        veg.mor = 1  # enabled

        c = self.mock_cohort(constants, attributes)
        veg.cohorts = [c]

        bedlevel_change = np.array([0.2, -0.3, 2, -2])  # positive is deposition --> burial

        veg.mortality_morphodynamic(bedlevel_change, burial_frac=1.0, scour_frac=1.0)

        expected_burial = np.where(bedlevel_change >= constants["hv"], 1, 0)
        assert np.allclose(c.potential_mort_burial, expected_burial)
        expected_scour = np.where(bedlevel_change <= -attributes["rootlength_0"], 1, 0)
        assert np.allclose(c.potential_mort_scour, expected_scour)


class TestApplyMortality:

    @staticmethod
    def mock_cohort(ls):
        class MockCohort:
            fraction = np.array([0.4, 0.6])
            lifestage = ls  # life stage numbering starts at 1, not 0
            potential_mort_burial = np.array([0.1, 0.1])
            potential_mort_scour  = np.array([0.1, 0.1])
        return MockCohort()


    @staticmethod
    def veg_obj():
        veg = VegetationSpecies.__new__(VegetationSpecies)
        # two life stages, one array per life stage
        veg.potential_mort_flood = [np.array([0.1, 0.1]),
                                    np.array([0.2, 0.2]),
                                    np.array([1.0, 1.0]),
                                    ]
        veg.potential_mort_desic = [np.array([0.1, 0.1]),
                                    np.array([0.2, 0.2]),
                                    np.array([1.0, 1.0]),
                                    ]
        veg.potential_mort_uproot = [np.array([0.0, 0.0]),
                                     np.array([0.5, 0.5]),
                                     np.array([1.0, 1.0]),
                                     ]
        return veg
    

    @pytest.mark.unit
    def test_apply_mortality_lifestage_indexing(self):
        veg = self.veg_obj()
        c = self.mock_cohort(2)
        veg.cohorts = [c]

        veg.apply_mortality()

        assert np.all(c.potential_mort_flood == veg.potential_mort_flood[1])
        assert np.all(c.potential_mort_desic == veg.potential_mort_desic[1])
        assert np.all(c.potential_mort_uproot == veg.potential_mort_uproot[1])


    @pytest.mark.unit
    def test_apply_mortality_multiplication(self):
        veg = self.veg_obj()
        c = self.mock_cohort(1)
        veg.cohorts = [c]

        expected_fld = c.fraction * veg.potential_mort_flood[0]
        expected_des = c.fraction * veg.potential_mort_desic[0]
        expected_upr = c.fraction * veg.potential_mort_uproot[0]
        expected_bur = c.fraction * c.potential_mort_burial
        expected_sco = c.fraction * c.potential_mort_scour

        veg.apply_mortality()

        assert np.allclose(c.applied_mort_flood, expected_fld)
        assert np.allclose(c.applied_mort_desic, expected_des)
        assert np.allclose(c.applied_mort_uproot, expected_upr)
        assert np.allclose(c.applied_mort_burial, expected_bur)
        assert np.allclose(c.applied_mort_scour, expected_sco)


    @pytest.mark.unit
    def test_apply_mortality_total_additive(self):
        veg = self.veg_obj()
        c = self.mock_cohort(1)
        veg.cohorts = [c]

        expected_fld = c.fraction * veg.potential_mort_flood[0]
        expected_des = c.fraction * veg.potential_mort_desic[0]
        expected_upr = c.fraction * veg.potential_mort_uproot[0]
        expected_bur = c.fraction * c.potential_mort_burial
        expected_sco = c.fraction * c.potential_mort_scour

        veg.apply_mortality()

        expected_total = expected_fld + expected_des + expected_upr + expected_bur + expected_sco
        assert np.allclose(c.applied_mort_total, expected_total)


    @pytest.mark.unit
    def test_apply_mortality_fully_removes(self):
        veg = self.veg_obj()
        c = self.mock_cohort(3)  # 100% mortality for all hydrodynamic modes
        veg.cohorts = [c]
        c.potential_mort_burial = np.ones(2)
        c.potential_mort_scour = np.ones(2)

        veg.apply_mortality()
        assert np.all(c.fraction == 0)


    @pytest.mark.unit
    def test_apply_mortality_identity(self):
        veg = self.veg_obj()
        veg.potential_mort_flood = [np.zeros(2)]
        veg.potential_mort_desic = [np.zeros(2)]
        veg.potential_mort_uproot = [np.zeros(2)]

        c = self.mock_cohort(1)
        veg.cohorts = [c]
        c.potential_mort_burial = np.zeros(2)
        c.potential_mort_scour = np.zeros(2)

        before = c.fraction.copy()

        veg.apply_mortality()
        assert np.allclose(c.fraction, before)


    @pytest.mark.unit
    def test_apply_mortality_multiple_cohorts(self):
        veg = self.veg_obj()
        c1 = self.mock_cohort(1)
        c2 = self.mock_cohort(3)
        veg.cohorts = [c1, c2]

        orig_fraction = c2.fraction.copy()
        veg.apply_mortality()

        # No applied mortality
        assert np.all(c1.applied_mort_uproot == 0.0)
        # 100% applied mortality, max is actual fraction present
        assert np.allclose(c2.applied_mort_uproot, orig_fraction)


# ------------------------------------------------------------ #
# tests for growth methods
# ------------------------------------------------------------ #

class TestGrowthMethods:

    @staticmethod
    def mock_attrs():
        class MockAttr:
            start_growth_ets = 4
            end_growth_ets = 8
            winter_ets = 10
            stemht_max = [1.0, 2.5]
            stemht_winter_max = [0.1, 0.2]
            stemdiam_max = [0.01, 0.02]
            rootlength_max = [1.0, 1.5]
            ht_growth_rates = [0.1, 0.2]
            diam_growth_rates = [0.001, 0.002]
            root_growth_rates = [0.1, 0.2]
        return MockAttr()


    @staticmethod
    def mock_cohort(ls):
        class MockCohort:
            lifestage  = ls  # life stage numbering starts at 1, not 0
            height     = 0.5
            diameter   = 0.005
            rootlength = 0.5
        return MockCohort()
    

    @pytest.mark.unit
    def test_ets_at_start_growth_no_change(self):
        """ ets == start_growth_ets: loop does not run, height unchanged """
        veg = VegetationSpecies.__new__(VegetationSpecies)
        veg.attrs = self.mock_attrs()

        c = self.mock_cohort(1)
        veg.cohorts = [c]

        veg.stemheight_growth(ets=4)
        veg.stemdiam_growth(ets=4)
        veg.root_growth(ets=4)

        assert c.height == 0.5
        assert c.diameter == 0.005
        assert c.rootlength == 0.5


    @pytest.mark.unit
    def test_ets_one_above_start_growth_triggers_growth(self):
        """ ets == start_growth_ets + 1: first ets where growth runs """
        veg = VegetationSpecies.__new__(VegetationSpecies)
        veg.attrs = self.mock_attrs()

        c = self.mock_cohort(1)
        veg.cohorts = [c]

        veg.stemheight_growth(ets=5)
        veg.stemdiam_growth(ets=5)
        veg.root_growth(ets=5)

        assert c.height == pytest.approx(0.6)
        assert c.diameter == pytest.approx(0.006)
        assert c.rootlength == pytest.approx(0.6)


    @pytest.mark.unit
    def test_ets_at_end_growth_still_grows(self):
        """ ets == end_growth_ets: growth branch is inclusive, height should increase """
        veg = VegetationSpecies.__new__(VegetationSpecies)
        veg.attrs = self.mock_attrs()

        c = self.mock_cohort(1)
        veg.cohorts = [c]

        veg.stemheight_growth(ets=8)
        veg.stemdiam_growth(ets=8)
        veg.root_growth(ets=8)

        assert c.height == pytest.approx(0.6)
        assert c.diameter == pytest.approx(0.006)
        assert c.rootlength == pytest.approx(0.6)


    @pytest.mark.unit
    def test_ets_after_end_growth_no_change(self):
        """ ets between end_growth_ets and winter_ets: neither branch fires, height unchanged """
        veg = VegetationSpecies.__new__(VegetationSpecies)
        veg.attrs = self.mock_attrs()

        c = self.mock_cohort(1)
        veg.cohorts = [c]

        veg.stemheight_growth(ets=9)
        veg.stemdiam_growth(ets=9)
        veg.root_growth(ets=9)

        assert c.height == 0.5
        assert c.diameter == 0.005
        assert c.rootlength == 0.5


    @pytest.mark.unit
    def test_growth_at_max_no_change(self):
        """Height == stemht_max: growth is blocked by the < check."""
        veg = VegetationSpecies.__new__(VegetationSpecies)
        veg.attrs = self.mock_attrs()

        c = self.mock_cohort(1)
        c.height = 1.0
        c.diameter = 0.01
        c.rootlength = 1.0
        veg.cohorts = [c]

        veg.stemheight_growth(ets=6)
        veg.stemdiam_growth(ets=6)
        veg.root_growth(ets=6)

        assert c.height == 1.0
        assert c.diameter == 0.01
        assert c.rootlength == 1.0


    @pytest.mark.unit
    def test_growth_above_max_no_change(self):
        """
        Height > stemht_max: growth is blocked, height stays the same.
        Theoretically should not occur, but this is handled by VegetationAttributes.
        """
        veg = VegetationSpecies.__new__(VegetationSpecies)
        veg.attrs = self.mock_attrs()

        c = self.mock_cohort(1)
        c.height = 2.0
        c.diameter = 0.1
        c.rootlength = 2.0
        veg.cohorts = [c]

        veg.stemheight_growth(ets=6)
        veg.stemdiam_growth(ets=6)
        veg.root_growth(ets=6)

        assert c.height == pytest.approx(2.0)
        assert c.diameter == pytest.approx(0.1)
        assert c.rootlength == pytest.approx(2.0)


    @pytest.mark.unit
    def test_winter_ets_above_max_reduced(self):
        """ Height > stemht_winter_max: height drops to winter max """
        veg = VegetationSpecies.__new__(VegetationSpecies)
        veg.attrs = self.mock_attrs()

        c = self.mock_cohort(1)
        veg.cohorts = [c]

        veg.stemheight_growth(ets=10)

        assert c.height == pytest.approx(0.1)


    @pytest.mark.unit
    def test_winter_ets_remains_constant(self):
        """ Diameter/Rootlength: remains constant """
        veg = VegetationSpecies.__new__(VegetationSpecies)
        veg.attrs = self.mock_attrs()

        c = self.mock_cohort(1)
        veg.cohorts = [c]

        veg.stemdiam_growth(ets=10)
        veg.root_growth(ets=10)

        assert c.diameter == pytest.approx(0.005)
        assert c.rootlength == pytest.approx(0.5)


    @pytest.mark.unit
    def test_multiple_cohorts_each_use_own_lifestage_index(self):
        """ Two cohorts with different lifestages each grow by their own rate """
        veg = VegetationSpecies.__new__(VegetationSpecies)
        veg.attrs = self.mock_attrs()

        c1 = self.mock_cohort(1)
        c2 = self.mock_cohort(2)
        veg.cohorts = [c1, c2]

        veg.stemheight_growth(ets=6)
        veg.stemdiam_growth(ets=6)
        veg.root_growth(ets=6)

        assert c1.height == pytest.approx(0.6)  # grew by 0.1
        assert c2.height == pytest.approx(0.7)  # grew by 0.2
        assert c1.diameter == pytest.approx(0.006)
        assert c2.diameter == pytest.approx(0.007)
        assert c1.rootlength == pytest.approx(0.6)
        assert c2.rootlength == pytest.approx(0.7)


    @pytest.mark.unit
    def test_empty_cohorts_no_error(self):
        """ Empty cohorts list: method completes without error """
        veg = VegetationSpecies.__new__(VegetationSpecies)
        veg.attrs = self.mock_attrs()

        veg.cohorts = []

        veg.stemheight_growth(ets=6)  # should not raise



class TestUpdateLifestageAndStemdensity:

    @staticmethod
    def mock_attrs(num_lifestage=3):
        class MockAttr:
            nls = num_lifestage
            years_max = [2, 3, 4]  # max years per lifestage
            stemdens = [500, 300, 100]  # density per lifestage
        return MockAttr()

    @staticmethod
    def mock_cohort(ls, ls_year, dens=500):
        class MockCohort:
            lifestage = ls
            lifestage_year = ls_year
            density = dens
        return MockCohort()


    @pytest.mark.unit
    def test_mid_lifestage_year_increments(self):
        """ cohort not at years_max: lifestage_year increments, nothing else changes """
        veg = VegetationSpecies.__new__(VegetationSpecies)
        veg.attrs = self.mock_attrs()

        c = self.mock_cohort(ls=1, ls_year=1)
        veg.cohorts = [c]
        veg.remove_old_cohorts = lambda inds: None

        veg.update_lifestage_and_stemdensity()

        assert c.lifestage == 1
        assert c.lifestage_year == 2
        assert c.density == 500


    @pytest.mark.unit
    def test_advances_to_next_lifestage_at_years_max(self):
        """ cohort at years_max but not final lifestage: lifestage advances, year resets to 1 """
        veg = VegetationSpecies.__new__(VegetationSpecies)
        veg.attrs = self.mock_attrs()

        c = self.mock_cohort(ls=1, ls_year=2)  # years_max[0] == 2
        veg.cohorts = [c]
        veg.remove_old_cohorts = lambda inds: None

        veg.update_lifestage_and_stemdensity()

        assert c.lifestage == 2
        assert c.lifestage_year == 1


    @pytest.mark.unit
    def test_density_updated_to_next_lifestage_on_advance(self):
        """ On lifestage advance, density is set to stemdens of the NEW lifestage """
        veg = VegetationSpecies.__new__(VegetationSpecies)
        veg.attrs = self.mock_attrs()

        c = self.mock_cohort(ls=1, ls_year=2)
        veg.cohorts = [c]
        veg.remove_old_cohorts = lambda inds: None

        veg.update_lifestage_and_stemdensity()

        assert c.density == 300  # stemdens[1], i.e. lifestage 2


    @pytest.mark.unit
    def test_cohort_at_final_lifestage_and_years_max_is_removed(self):
        """ cohort at final lifestage and years_max: remove_old_cohorts is called with its index """
        veg = VegetationSpecies.__new__(VegetationSpecies)
        veg.attrs = self.mock_attrs()
        c = self.mock_cohort(ls=3, ls_year=4)  # nls=3, years_max[2]==4
        veg.cohorts = [c]

        removed = []
        veg.remove_old_cohorts = lambda inds: removed.extend(inds)

        veg.update_lifestage_and_stemdensity()

        assert removed == [0]


    @pytest.mark.unit
    def test_remove_old_cohorts_not_called_when_no_cohorts_expire(self):
        """ No cohorts at end of final lifestage: remove_old_cohorts is never called """
        veg = VegetationSpecies.__new__(VegetationSpecies)
        veg.attrs = self.mock_attrs()
        c = self.mock_cohort(ls=1, ls_year=1)
        veg.cohorts = [c]

        called = []
        veg.remove_old_cohorts = lambda inds: called.append(inds)

        veg.update_lifestage_and_stemdensity()

        assert called == []


    @pytest.mark.unit
    def test_mixed_cohorts_each_handled_independently(self):
        """ One cohort mid-lifestage, one advancing, one expiring: all handled correctly """
        veg = VegetationSpecies.__new__(VegetationSpecies)
        veg.attrs = self.mock_attrs()

        c_mid      = self.mock_cohort(ls=1, ls_year=1)  # will increment year
        c_advance  = self.mock_cohort(ls=1, ls_year=2)  # will advance lifestage
        c_expire   = self.mock_cohort(ls=3, ls_year=4)  # will be removed
        veg.cohorts = [c_mid, c_advance, c_expire]

        removed = []
        veg.remove_old_cohorts = lambda inds: removed.extend(inds)

        veg.update_lifestage_and_stemdensity()

        assert c_mid.lifestage_year == 2         # incremented
        assert c_advance.lifestage == 2          # advanced
        assert c_advance.lifestage_year == 1     # reset
        assert removed == [2]                    # only the expiring cohort's index


    @pytest.mark.unit
    def test_remove_old_cohorts_removes_correct_cohorts(self):
        """ Only cohorts at the specified indices are removed; others are retained """
        veg = VegetationSpecies.__new__(VegetationSpecies)
        c0 = self.mock_cohort(ls=1, ls_year=1)
        c1 = self.mock_cohort(ls=1, ls_year=1)
        c2 = self.mock_cohort(ls=1, ls_year=1)
        veg.cohorts = [c0, c1, c2]
        veg.remove_old_cohorts([0, 2])
        assert veg.cohorts == [c1]