import pytest
from unittest.mock import MagicMock, patch
import numpy as np
import sys, os

sys.path.append(os.path.realpath(os.path.dirname(__file__)+"/../.."))
#from conftest import make_anuga_domain_and_engine as make_anuga
from dycove.sim.engines.ANUGA_hydro import ANUGA, AnugaEngine


@pytest.mark.anuga
@pytest.mark.unit
def test_anuga_hydro_components():
    assert callable(ANUGA)
    assert callable(AnugaEngine)


@pytest.mark.anuga
@pytest.mark.unit
def test_engine_abstractmethods(required_engine_methods):
    # Check for methods defined by abstract class base.HydroEngineBase
    for method in required_engine_methods:
        assert hasattr(AnugaEngine, method)


class TestAnugaEngine:

    @staticmethod
    def mock_domain():
        domain = MagicMock()
        domain.get_datadir.return_value = "/tmp/test_model"
        domain.evolve.return_value = iter([0])  # simulate one yieldstep
        return domain
    

    @staticmethod
    def mock_veg():
        veg = MagicMock()
        veg.get_drag.return_value = 1.0
        return veg


    @staticmethod
    def make_engine(domain, vegetation=None, save_interval=3600, 
                    myid=1, nprocs=2,  # use as default to avoid r.report call under step()
                    barrier=MagicMock()  # one test will implement something here
                    ):
        with patch("dycove.sim.engines.ANUGA_hydro._import_anuga",
                return_value=(myid, nprocs, MagicMock(), barrier)):
            engine = AnugaEngine(domain, vegetation=vegetation)
        engine.save_interval = save_interval
        return engine


    @pytest.mark.anuga
    @pytest.mark.unit
    def test_anuga_obj_has_anuga_engine_attrs_after_init(self):
        domain = self.mock_domain()

        # patch separately here because testing ANUGA class not AnugaEngine directly
        with patch("dycove.sim.engines.ANUGA_hydro._import_anuga", 
                return_value=(1, 2, MagicMock(), MagicMock())):
            model = ANUGA(domain)

        assert isinstance(model.engine, AnugaEngine)
        assert model.engine.domain is domain


    @pytest.mark.anuga
    @pytest.mark.unit
    def test_initialize_sets_morphology_false(self):
        """ morphology is always set to False regardless of vegetation """
        domain = self.mock_domain()
        engine = self.make_engine(domain)
        engine.initialize()
        assert engine.morphology == False


    @pytest.mark.anuga
    @pytest.mark.unit
    def test_initialize_creates_baptist_operator_when_veg_present(self):
        """ Baptist_operator is constructed with drag from veg.get_drag() when veg is set """
        domain = self.mock_domain()
        veg = self.mock_veg()
        engine = self.make_engine(domain, vegetation=veg)

        with patch("dycove.sim.engines.ANUGA_hydro.Baptist_operator") as mock_baptist:
            engine.initialize()

        mock_baptist.assert_called_once_with(domain, 
                                             veg_diameter=0, veg_density=0, veg_height=0, 
                                             drag=1.0)


    @pytest.mark.anuga
    @pytest.mark.unit
    def test_initialize_skips_baptist_operator_when_no_veg(self):
        """ Baptist_operator is not created when vegetation is None """
        domain = self.mock_domain()
        engine = self.make_engine(domain)

        with patch("dycove.sim.engines.ANUGA_hydro.Baptist_operator") as mock_baptist:
            engine.initialize()

        mock_baptist.assert_not_called()


    @pytest.mark.anuga
    @pytest.mark.unit
    def test_step_yieldstep_is_min_of_seconds_and_save_interval(self):
        """ yieldstep passed to domain.evolve is capped at save_interval """
        domain = self.mock_domain()
        engine = self.make_engine(domain)

        engine.step(seconds=43200)

        call_kwargs = domain.evolve.call_args[1]
        assert call_kwargs["yieldstep"] == 3600  # capped at save_interval, not 43200


    @pytest.mark.anuga
    @pytest.mark.unit
    def test_step_yieldstep_smaller_than_save_interval(self):
        """ yieldstep equals `seconds` when seconds < save_interval """
        domain = self.mock_domain()
        engine = self.make_engine(domain)

        engine.step(seconds=900)

        call_kwargs = domain.evolve.call_args[1]
        assert call_kwargs["yieldstep"] == 900


    @pytest.mark.anuga
    @pytest.mark.unit
    def test_step_skip_step_is_true_after_first_call(self):
        """ skip_step is False on first call but True after, preventing duplicate yieldsteps """
        domain = self.mock_domain()
        engine = self.make_engine(domain)

        assert engine.skip_step == False
        engine.step(seconds=3600)
        assert engine.skip_step == True


    @pytest.mark.anuga
    @pytest.mark.unit
    def test_step_passes_skip_step_to_evolve(self):
        """ domain.evolve receives the current value of skip_step """
        domain = self.mock_domain()
        engine = self.make_engine(domain)

        engine.step(seconds=3600)

        call_kwargs = domain.evolve.call_args[1]
        assert call_kwargs["skip_initial_step"] == False


    @pytest.mark.anuga
    @pytest.mark.unit
    def test_step_calls_barrier_twice(self):
        """ barrier() is called at start and end to synchronize parallel processes """
        domain = self.mock_domain()

        domain.evolve.side_effect = lambda **kwargs: call_order.append("evolve") or iter([0])
        barrier = MagicMock(side_effect=lambda: call_order.append("barrier"))
                
        engine = self.make_engine(domain, barrier=barrier)

        call_order = []
        engine.step(seconds=3600)

        assert call_order == ["barrier", "evolve", "barrier"]






    # @pytest.mark.anuga
    # @pytest.mark.unit
    # def test_cell_count():
    #     domain, engine = make_anuga()

    #     # get_cell_count() uses elevation array to return number of cells, so do a cross-check
    #     n_cells = engine.get_cell_count()
    #     elevation_c = engine.get_elevation()  # centroids
    #     assert n_cells == len(elevation_c)


    # @pytest.mark.anuga
    # @pytest.mark.unit
    # def test_get_velocity_and_depth():
    #     domain, engine = make_anuga(momentum=True)

    #     # Get velocity and depth arrays from domain quantities
    #     stage = domain.quantities["stage"].centroid_values
    #     elev = domain.quantities["elevation"].centroid_values
    #     xmom = domain.quantities["xmomentum"].centroid_values
    #     ymom = domain.quantities["ymomentum"].centroid_values
    #     dep = stage - elev
    #     vel = np.sqrt((xmom/dep)**2 + (ymom/dep)**2)

    #     # Get from engine method
    #     velocity, depth = engine.get_velocity_and_depth()

    #     # Cross-check values
    #     assert velocity == pytest.approx(vel, rel=1e-6)
    #     assert depth == pytest.approx(dep, rel=1e-6)


    # @pytest.mark.anuga
    # @pytest.mark.unit
    # def test_get_set_vegetation(constants):
    #     from dycove.sim.engines.ANUGA_baptist import Baptist_operator

    #     c = constants
    #     domain, engine = make_anuga()

    #     engine.Baptist = Baptist_operator(domain)
    #     engine.set_vegetation(c["m"], c["D"], c["hv"])

    #     # Check that vegetation was set correctly in domain/Baptist
    #     dens, diam, ht = engine.get_vegetation()

    #     assert np.all(dens == pytest.approx(c["m"], rel=1e-6))
    #     assert np.all(diam == pytest.approx(c["D"], rel=1e-6))
    #     assert np.all(ht == pytest.approx(c["hv"], rel=1e-6))