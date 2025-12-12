from pytest import approx, mark
import shutil
import numpy as np
import sys, os

sys.path.append(os.path.realpath(os.path.dirname(__file__)+"/.."))
from conftest import anuga_domain_and_engine


@mark.anuga
def test_anuga_hydro_components():
    from dycove import ANUGA_hydro
    assert hasattr(ANUGA_hydro, "ANUGA")
    assert hasattr(ANUGA_hydro, "AnugaEngine")


@mark.anuga
def test_engine_abstractmethods(required_engine_methods):
    from dycove import ANUGA_hydro
    # Check for methods defined by abstract class base.HydroEngineBase
    for method in required_engine_methods:
        assert hasattr(ANUGA_hydro.AnugaEngine, method)
    

@mark.anuga
def test_step():
    domain, engine = anuga_domain_and_engine(friction=True, dirichlet=True)

    # Timestamp and step model by small dt
    t0 = domain.get_time()
    dt = 1
    engine.step(dt)

    # Check domain time state
    assert domain.get_time() - dt == approx(t0, rel=1e-6)


@mark.anuga
def test_cell_count():
    domain, engine = anuga_domain_and_engine()

    # get_cell_count() uses elevation array to return number of cells, so do a cross-check
    n_cells = engine.get_cell_count()
    elevation_c = engine.get_elevation()  # centroids
    assert n_cells == len(elevation_c)


@mark.anuga
def test_get_velocity_and_depth():
    domain, engine = anuga_domain_and_engine(momentum=True)

    # Get velocity and depth arrays from domain quantities
    stage = domain.quantities["stage"].centroid_values
    elev = domain.quantities["elevation"].centroid_values
    xmom = domain.quantities["xmomentum"].centroid_values
    ymom = domain.quantities["ymomentum"].centroid_values
    dep = stage - elev
    vel = np.sqrt((xmom/dep)**2 + (ymom/dep)**2)

    # Get from engine method
    velocity, depth = engine.get_velocity_and_depth()

    # Cross-check values
    assert velocity == approx(vel, rel=1e-6)
    assert depth == approx(dep, rel=1e-6)


@mark.anuga
def test_get_set_vegetation(constants):
    from dycove.sim.engines.ANUGA_baptist import Baptist_operator

    c = constants
    domain, engine = anuga_domain_and_engine()

    engine.Baptist = Baptist_operator(domain)
    engine.set_vegetation(c["m"], c["D"], c["hv"])

    # Check that vegetation was set correctly in domain/Baptist
    dens, diam, ht = engine.get_vegetation()

    assert np.all(dens == approx(c["m"], rel=1e-6))
    assert np.all(diam == approx(c["D"], rel=1e-6))
    assert np.all(ht == approx(c["hv"], rel=1e-6))


@mark.anuga
def test_merge_parallel_veg(sample_data_path, tmp_data_path, monkeypatch):
    # Define fake VegetationSpecies class
    def mock_veg():
        class MockVeg:
            cohorts = [1, 2, 3, 4]
        return MockVeg()
    
    from dycove.sim.outputs import OutputManager

    domain, engine = anuga_domain_and_engine()

    # Copy existing sample output files so we can place them back cleanly afterwards
    tmp_data_path.mkdir()
    shutil.copytree(sample_data_path / "veg_output", tmp_data_path / "veg_output")
    shutil.copy(sample_data_path / "domain_P2_0.sww", tmp_data_path)
    shutil.copy(sample_data_path / "domain_P2_1.sww", tmp_data_path)

    # So that merge_parallel_veg() can find the files
    domain.set_name("domain_P2_0")  # this is the specific name for myid 0, so it has suffix "_0"
    domain.set_datadir(tmp_data_path)

    # Mock some values to the engine for the sake of the test
    engine.model_dir = domain.get_datadir()
    engine.numprocs = 2
    engine.veg = mock_veg()

    # Normally created upon creation of HydroSimulationBase object
    om = OutputManager(engine)

    # Shortcutting to this method allows bypass of creating VegetationSpecies object
    engine.merge_parallel_veg(om)

    # Remove temporary directory
    shutil.rmtree(tmp_data_path)
