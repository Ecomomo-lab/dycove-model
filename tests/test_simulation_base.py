from dycove.sim.base import HydroSimulationBase

def test_sim_base_has_required_methods():
    assert hasattr(HydroSimulationBase, "run_simulation")
    assert hasattr(HydroSimulationBase, "loop_hydrodynamics")
    assert hasattr(HydroSimulationBase, "step_hydrodynamics")