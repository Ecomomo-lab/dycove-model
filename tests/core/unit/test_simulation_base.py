import pytest
from unittest.mock import MagicMock, patch
from dycove.sim.base import HydroSimulationBase


@pytest.mark.unit
def test_sim_base_has_required_methods():
    assert hasattr(HydroSimulationBase, "run_simulation")
    assert hasattr(HydroSimulationBase, "run_simulation_without_vegetation")
    assert hasattr(HydroSimulationBase, "run_simulation_with_vegetation")
    assert hasattr(HydroSimulationBase, "eco_hydro_loop")
    assert hasattr(HydroSimulationBase, "hydro_step")


class ConcreteHydroSimulation(HydroSimulationBase):
    pass  # ABC requires no additional methods to be implemented

class TestHydroSimulationBase:

    @staticmethod
    def mock_engine(veg=True, morphology=False):
        engine = MagicMock()
        engine.veg = MagicMock() if veg else None
        engine.model_dir = "/tmp/test_model"
        engine.morphology = morphology
        engine.morph_vars = {"MorFac": 1}
        engine.get_refdate.return_value = "2000-01-01"
        engine.get_cell_count.return_value = 10
        return engine


    @staticmethod
    def make_sim(engine):
        """ Instantiate HydroSimulationBase with all heavy dependencies mocked out """
        with patch("dycove.sim.coupler.VegetationCoupler"), \
             patch("dycove.sim.outputs.OutputManager"), \
             patch("dycove.sim.simulation_data.SimulationTimeState"), \
             patch("dycove.sim.simulation_data.HydrodynamicStats"):
            sim = ConcreteHydroSimulation(engine)  # a minimal concrete subclass of HydroSimulationBase
        return sim


    @pytest.mark.unit
    def test_run_simulation_routes_correctly_when_veg_active(self):
        """ When veg is active, run_simulation_with_vegetation is called, not the other """
        engine = self.mock_engine(veg=True)
        sim = self.make_sim(engine)
        sim.run_simulation_with_vegetation = MagicMock()
        sim.run_simulation_without_vegetation = MagicMock()

        with patch("dycove.sim.simulation_data.SimulationTimeState"), \
             patch("dycove.sim.simulation_data.HydrodynamicStats"), \
             patch("dycove.sim.outputs.OutputManager"):
            sim.run_simulation(sim_time=1, sim_time_unit="hydrodynamic days")

        sim.run_simulation_with_vegetation.assert_called_once()
        sim.run_simulation_without_vegetation.assert_not_called()


    @pytest.mark.unit
    def test_run_simulation_routes_correctly_when_veg_inactive(self):
        engine = self.mock_engine(veg=False)
        sim = self.make_sim(engine)
        sim.run_simulation_with_vegetation = MagicMock()
        sim.run_simulation_without_vegetation = MagicMock()

        with patch("dycove.sim.simulation_data.SimulationTimeState"), \
             patch("dycove.sim.simulation_data.HydrodynamicStats"), \
             patch("dycove.sim.outputs.OutputManager"):
            sim.run_simulation(sim_time=1, sim_time_unit="hydrodynamic days")

        sim.run_simulation_without_vegetation.assert_called_once()
        sim.run_simulation_with_vegetation.assert_not_called()


    @pytest.mark.unit
    def test_verify_eco_args_raises_if_veg_inactive_and_nondefault_args(self):
        engine = self.mock_engine(veg=False)
        sim = self.make_sim(engine)
        with pytest.raises(ValueError):
            sim.run_simulation(sim_time=1, sim_time_unit="hydrodynamic days", n_ets=7)


    @pytest.mark.unit
    def test_verify_eco_args_does_not_raise_if_veg_active(self):
        engine = self.mock_engine(veg=True)
        sim = self.make_sim(engine)
        sim.run_simulation_with_vegetation = MagicMock()

        with patch("dycove.sim.simulation_data.SimulationTimeState"), \
             patch("dycove.sim.simulation_data.HydrodynamicStats"), \
             patch("dycove.sim.outputs.OutputManager"):
            sim.run_simulation(sim_time=1, sim_time_unit="hydrodynamic days", n_ets=7)  # should not raise b/c veg=True


    @pytest.mark.unit
    def test_verify_eco_args_does_not_raise_if_default_args(self):
        engine = self.mock_engine(veg=False)
        sim = self.make_sim(engine)
        sim.run_simulation_with_vegetation = MagicMock()

        with patch("dycove.sim.simulation_data.SimulationTimeState"), \
             patch("dycove.sim.simulation_data.HydrodynamicStats"), \
             patch("dycove.sim.outputs.OutputManager"):
            sim.run_simulation(sim_time=1, sim_time_unit="hydrodynamic days")  # should not raise b/c default


    @pytest.mark.unit
    def test_engine_initialize_always_called(self):
        engine = self.mock_engine(veg=False)
        sim = self.make_sim(engine)
        sim.run_simulation_without_vegetation = MagicMock()

        with patch("dycove.sim.simulation_data.SimulationTimeState"), \
             patch("dycove.sim.simulation_data.HydrodynamicStats"), \
             patch("dycove.sim.outputs.OutputManager"):
            sim.run_simulation(sim_time=1, sim_time_unit="hydrodynamic days")

        engine.initialize.assert_called_once()


    @pytest.mark.unit
    def test_veg_loop_calls_steps_n_times(self):
        """eco_hydro_loop, veg_coupler.update, and save_vegetation_step each called n_veg_steps times."""
        engine = self.mock_engine(veg=True)
        sim = self.make_sim(engine)
        sim.eco_hydro_loop = MagicMock()
        sim.veg_coupler = MagicMock()
        sim.outputs = MagicMock()
        sim.finalize_simulation = MagicMock()
        sim.simstate = MagicMock()
        sim.hydrostats = MagicMock()
        sim.simstate.n_veg_steps = 3

        with patch("dycove.sim.base.r.print_runtime_updates"):
            sim.run_simulation_with_vegetation()

        assert sim.eco_hydro_loop.call_count == 3
        assert sim.veg_coupler.update.call_count == 3
        assert sim.outputs.save_vegetation_step.call_count == 3
        sim.finalize_simulation.assert_called_once()


    @pytest.mark.unit
    def test_eco_hydro_loop_calls_correct_methods(self):
        """ hydro_step and hydrostats.update are called once per substep """
        engine = self.mock_engine(veg=True)
        sim = self.make_sim(engine)
        sim.hydro_step = MagicMock()
        sim.hydrostats = MagicMock()
        engine.get_velocity_and_depth.return_value = ("vel", "depth")
        engine.get_elevation.return_value = 0.0

        sim.eco_hydro_loop(n_substeps=3, interval=900)

        assert sim.hydro_step.call_count == 3
        assert sim.hydrostats.update.call_count == 3
        sim.hydrostats.reset.assert_called_once()
        assert sim.engine.get_elevation.call_count == 2