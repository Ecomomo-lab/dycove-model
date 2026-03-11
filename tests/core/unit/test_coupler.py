import pytest
from unittest.mock import MagicMock, patch
import numpy as np

from dycove.sim.coupler import VegetationCoupler


class TestVegetationCoupler:

    @staticmethod
    def mock_engine():
        engine = MagicMock()
        engine.veg = MagicMock()
        return engine


    @staticmethod
    def mock_simstate(ets=3):
        simstate = MagicMock()
        simstate.ets = ets
        def set_ets_to_1():
            simstate.ets = 1
        simstate.update_ets = MagicMock()
        return simstate


    @staticmethod
    def mock_hydro():
        hydro = MagicMock()
        hydro.h_min = np.array([0.0, 0.2])
        hydro.h_max = np.array([0.4, 0.5])
        hydro.fl_dr = 0.1
        hydro.flood_frac.return_value = np.array([0.4, 0.6])
        hydro.dry_frac.return_value = np.array([0.6, 0.4])
        hydro.v_max_95th.return_value = np.array([0.8, 1.2])
        hydro.bedlevel_diff = np.array([-0.2, 0.2])
        return hydro


    @pytest.mark.unit
    def test_update_calls_all_steps_in_order(self):
        """ All six update steps are called, in the correct order """
        coupler = VegetationCoupler(self.mock_engine())

        # Replace all methods with mocks and track call order
        call_order = []
        methods_in_order = ["lifestage_update", 
                            "apply_growth", 
                            "do_colonization",
                            "compute_mortality", 
                            "compute_veg_model_quantities", 
                            "push_veg_to_hydro"]
        for method in methods_in_order:
            def mock_method(*args, name=method, **kwargs):
                call_order.append(name)
            setattr(coupler, method, MagicMock(side_effect=mock_method))

        coupler.update(MagicMock(), MagicMock())  # mock both dataclasses, contents don't matter
        
        assert call_order == methods_in_order


    @pytest.mark.unit
    def test_lifestage_update_calls_update_ets(self):
        """ simstate.update_ets() is always called """
        coupler = VegetationCoupler(self.mock_engine())
        simstate = self.mock_simstate(ets=3)

        coupler.lifestage_update(simstate)

        simstate.update_ets.assert_called_once()


    @pytest.mark.unit
    def test_lifestage_update_triggers_lifestage_when_ets_is_1(self):
        """ update_lifestage_and_stemdensity is called when ets == 1 (new eco year) """
        coupler = VegetationCoupler(self.mock_engine())
        simstate = self.mock_simstate(ets=1)

        # update_ets is mocked, so no progression of ets
        coupler.lifestage_update(simstate)

        coupler.veg.update_lifestage_and_stemdensity.assert_called_once()

    
    @pytest.mark.unit
    def test_lifestage_update_skips_lifestage_when_ets_is_not_1(self):
        """ update_lifestage_and_stemdensity is not called mid-year """
        coupler = VegetationCoupler(self.mock_engine())
        simstate = self.mock_simstate(ets=3)

        coupler.lifestage_update(simstate)

        coupler.veg.update_lifestage_and_stemdensity.assert_not_called()


    @pytest.mark.unit
    def test_apply_growth_delegates_correctly(self):
        coupler = VegetationCoupler(self.mock_engine())
        simstate = self.mock_simstate(ets=5)

        coupler.apply_growth(simstate)

        coupler.veg.stemheight_growth.assert_called_once_with(5)
        coupler.veg.stemdiam_growth.assert_called_once_with(5)
        coupler.veg.root_growth.assert_called_once_with(5)


    @pytest.mark.unit
    def test_do_colonization_delegates_correctly(self):
        coupler = VegetationCoupler(self.mock_engine())
        simstate = self.mock_simstate(ets=5)
        hydro = self.mock_hydro()

        coupler.do_colonization(simstate, hydro)

        coupler.veg.colonization.assert_called_once_with(5, hydro.h_min, hydro.h_max, hydro.fl_dr)


    @pytest.mark.unit
    def test_compute_mortality_passes_correct_values(self):
        """ mortality() receives values from the correct hydrostats methods """
        coupler = VegetationCoupler(self.mock_engine())
        hydro = self.mock_hydro()

        coupler.compute_mortality(hydro)

        call_kwargs = coupler.veg.mortality.call_args
        hydro_vars  = call_kwargs[0][0]
        morpho_vars = call_kwargs[0][1]

        assert np.allclose(hydro_vars["fld_frac"], hydro.flood_frac())
        assert np.allclose(hydro_vars["dry_frac"], hydro.dry_frac())
        assert np.allclose(hydro_vars["vel_max"], hydro.v_max_95th())
        assert np.allclose(morpho_vars["bl_diff"], hydro.bedlevel_diff)


    @pytest.mark.unit
    def test_compute_veg_model_quantities_delegates_correctly(self):
        coupler = VegetationCoupler(self.mock_engine())
        coupler.compute_veg_model_quantities()
        coupler.veg.compute_veg_model_quantities.assert_called_once()


    @pytest.mark.unit
    def test_push_veg_to_hydro_skips_when_stemdensity_is_none(self):
        """ If veg.stemdensity is None, engine.get_vegetation is never called """
        engine = self.mock_engine()
        coupler = VegetationCoupler(engine)
        coupler.veg.stemdensity = None

        coupler.push_veg_to_hydro()

        engine.get_vegetation.assert_not_called()
        engine.set_vegetation.assert_not_called()

    
    @pytest.mark.unit
    def test_push_veg_to_hydro_sets_vegetation(self):
        """ When stemdensity is set, engine.set_vegetation is called with updated arrays """
        engine = self.mock_engine()
        coupler = VegetationCoupler(engine)

        coupler.veg.stemdensity  = np.array([100, 200])
        coupler.veg.stemdiameter = np.array([0.01, 0.02])
        coupler.veg.stemheight   = np.array([0.5, 0.6])

        # Would be equal to the above in practice, but for the test, just use zeros
        engine.get_vegetation.return_value = (np.zeros(2),
                                              np.zeros(2),
                                              np.zeros(2),
                                              )
        coupler.push_veg_to_hydro()

        engine.set_vegetation.assert_called_once()


    @pytest.mark.unit
    def test_push_veg_to_hydro_only_overwrites_active_portion(self):
        """ Only the first N elements of the hydro arrays are replaced, rest stay zero """
        engine = self.mock_engine()
        coupler = VegetationCoupler(engine)

        coupler.veg.stemdensity  = np.array([100, 200])
        coupler.veg.stemdiameter = np.array([0.01, 0.02])
        coupler.veg.stemheight   = np.array([0.5, 0.6])

        # Larger arrays, in case of DFM which includes boundary cells when pulled
        engine.get_vegetation.return_value = (np.zeros(5),
                                              np.zeros(5),
                                              np.zeros(5),
                                              )
        coupler.push_veg_to_hydro()

        passed_density = engine.set_vegetation.call_args[0][0]
        passed_diameter = engine.set_vegetation.call_args[0][1]
        passed_height = engine.set_vegetation.call_args[0][2]
        assert list(passed_density[:2]) == [100, 200]    # overwritten
        assert list(passed_diameter[:2]) == [0.01, 0.02]  # overwritten
        assert list(passed_height[:2]) == [0.5, 0.6]    # overwritten
        assert list(passed_density[2:]) == [0.0, 0.0, 0.0]  # untouched