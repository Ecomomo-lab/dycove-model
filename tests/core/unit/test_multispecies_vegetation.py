import pytest
from unittest.mock import MagicMock
from dycove.sim.vegetation import MultipleVegetationSpecies


# ------------------------------------------------------------ #
# tests for multiple species delegation
# ------------------------------------------------------------ #

class TestMultipleVegetationSpecies:

    @staticmethod
    def mock_species(mor=1):
        sp = MagicMock()
        sp.mor = mor
        sp.cohorts = []
        return sp


    @pytest.mark.unit
    def test_init_raises_if_mor_values_differ(self):
        """ Mixed mor values across species should raise ValueError on construction """
        sp1 = self.mock_species(mor=1)
        sp2 = self.mock_species(mor=0)
        with pytest.raises(ValueError):
            MultipleVegetationSpecies([sp1, sp2])


    @pytest.mark.unit
    def test_init_sets_mor_from_first_species(self):
        """ self.mor should be inherited from the first species when all agree """
        sp1 = self.mock_species(mor=1)
        sp2 = self.mock_species(mor=1)
        mvs = MultipleVegetationSpecies([sp1, sp2])
        assert mvs.mor == 1


    @pytest.mark.unit
    def test_cohorts_flattens_across_species(self):
        """ cohorts property should return a flat list of all cohorts from all species """
        sp1 = self.mock_species()
        sp2 = self.mock_species()
        sp1.cohorts = ["c1", "c2"]
        sp2.cohorts = ["c3"]
        mvs = MultipleVegetationSpecies([sp1, sp2])
        assert mvs.cohorts == ["c1", "c2", "c3"]


    @pytest.mark.unit
    def test_cohorts_empty_when_no_cohorts(self):
        """ cohorts property returns empty list when all species have no cohorts """
        sp1 = self.mock_species()
        sp2 = self.mock_species()
        mvs = MultipleVegetationSpecies([sp1, sp2])
        assert mvs.cohorts == []


    @pytest.mark.unit
    def test_stemheight_growth_delegates_to_all_species(self):
        sp1, sp2 = self.mock_species(), self.mock_species()
        mvs = MultipleVegetationSpecies([sp1, sp2])
        mvs.stemheight_growth(ets=5)
        sp1.stemheight_growth.assert_called_once_with(5)
        sp2.stemheight_growth.assert_called_once_with(5)


    @pytest.mark.unit
    def test_stemdiam_growth_delegates_to_all_species(self):
        sp1, sp2 = self.mock_species(), self.mock_species()
        mvs = MultipleVegetationSpecies([sp1, sp2])
        mvs.stemdiam_growth(ets=5)
        sp1.stemdiam_growth.assert_called_once_with(5)
        sp2.stemdiam_growth.assert_called_once_with(5)


    @pytest.mark.unit
    def test_root_growth_delegates_to_all_species(self):
        sp1, sp2 = self.mock_species(), self.mock_species()
        mvs = MultipleVegetationSpecies([sp1, sp2])
        mvs.root_growth(ets=5)
        sp1.root_growth.assert_called_once_with(5)
        sp2.root_growth.assert_called_once_with(5)


    @pytest.mark.unit
    def test_mortality_delegates_to_all_species(self):
        sp1, sp2 = self.mock_species(), self.mock_species()
        mvs = MultipleVegetationSpecies([sp1, sp2])
        mvs.mortality(hydro_vars="h", morpho_vars="m")
        sp1.mortality.assert_called_once_with("h", "m")
        sp2.mortality.assert_called_once_with("h", "m")


    @pytest.mark.unit
    def test_update_lifestage_delegates_to_all_species(self):
        sp1, sp2 = self.mock_species(), self.mock_species()
        mvs = MultipleVegetationSpecies([sp1, sp2])
        mvs.update_lifestage_and_stemdensity()
        sp1.update_lifestage_and_stemdensity.assert_called_once()
        sp2.update_lifestage_and_stemdensity.assert_called_once()


    @pytest.mark.unit
    def test_colonization_delegates_with_combined_cohorts(self):
        """colonization passes the flattened cohort list as combined_cohorts to each species."""
        sp1, sp2 = self.mock_species(), self.mock_species()
        sp1.cohorts = ["c1"]
        sp2.cohorts = ["c2"]
        mvs = MultipleVegetationSpecies([sp1, sp2])
        mvs.colonization(ets=5, min_depths="mn", max_depths="mx", fl_dr="fl")
        sp1.colonization.assert_called_once_with(5, "mn", "mx", "fl", combined_cohorts=["c1", "c2"])
        sp2.colonization.assert_called_once_with(5, "mn", "mx", "fl", combined_cohorts=["c1", "c2"])


    @pytest.mark.unit
    def test_get_drag_returns_mean_of_species_drags(self):
        """get_drag should return the mean of each species' get_drag value."""
        sp1, sp2 = self.mock_species(), self.mock_species()
        sp1.get_drag.return_value = 0.4
        sp2.get_drag.return_value = 0.6
        mvs = MultipleVegetationSpecies([sp1, sp2])
        assert mvs.get_drag() == pytest.approx(0.5)