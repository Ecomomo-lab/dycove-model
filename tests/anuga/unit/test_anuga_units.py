import pytest
import numpy as np
import sys, os

sys.path.append(os.path.realpath(os.path.dirname(__file__)+"/../.."))
from conftest import make_anuga_domain_and_engine as make_anuga


@pytest.mark.anuga
@pytest.mark.unit
def test_anuga_hydro_components():
    from dycove import ANUGA_hydro
    assert callable(ANUGA_hydro.ANUGA)
    assert callable(ANUGA_hydro.AnugaEngine)


@pytest.mark.anuga
@pytest.mark.unit
def test_engine_abstractmethods(required_engine_methods):
    from dycove import ANUGA_hydro
    # Check for methods defined by abstract class base.HydroEngineBase
    for method in required_engine_methods:
        assert hasattr(ANUGA_hydro.AnugaEngine, method)


@pytest.mark.anuga
@pytest.mark.unit
def test_cell_count():
    domain, engine = make_anuga()

    # get_cell_count() uses elevation array to return number of cells, so do a cross-check
    n_cells = engine.get_cell_count()
    elevation_c = engine.get_elevation()  # centroids
    assert n_cells == len(elevation_c)


@pytest.mark.anuga
@pytest.mark.unit
def test_get_velocity_and_depth():
    domain, engine = make_anuga(momentum=True)

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
    assert velocity == pytest.approx(vel, rel=1e-6)
    assert depth == pytest.approx(dep, rel=1e-6)


@pytest.mark.anuga
@pytest.mark.unit
def test_get_set_vegetation(constants):
    from dycove.sim.engines.ANUGA_baptist import Baptist_operator

    c = constants
    domain, engine = make_anuga()

    engine.Baptist = Baptist_operator(domain)
    engine.set_vegetation(c["m"], c["D"], c["hv"])

    # Check that vegetation was set correctly in domain/Baptist
    dens, diam, ht = engine.get_vegetation()

    assert np.all(dens == pytest.approx(c["m"], rel=1e-6))
    assert np.all(diam == pytest.approx(c["D"], rel=1e-6))
    assert np.all(ht == pytest.approx(c["hv"], rel=1e-6))