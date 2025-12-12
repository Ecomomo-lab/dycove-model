
import numpy as np
from pytest import approx, mark

from conftest import anuga_domain_and_engine


@mark.anuga
def test_baptist_quantities():
    from dycove.sim.engines.ANUGA_baptist import Baptist_operator

    domain, engine = anuga_domain_and_engine()

    Baptist = Baptist_operator(domain)

    # Test for quantities added to domain object
    assert "veg_diameter" in domain.quantities
    assert "veg_density" in domain.quantities
    assert "veg_height" in domain.quantities

    # Test for Baptist attributes in __init__
    assert hasattr(Baptist, "g")
    assert hasattr(Baptist, "K")
    assert hasattr(Baptist, "Cd")
    assert hasattr(Baptist, "depth")
    assert hasattr(Baptist, "bed_friction")


@mark.anuga
def test_baptist_Cv_calculation(constants):
    from dycove.sim.engines.ANUGA_baptist import Baptist_operator

    c = constants
    domain, engine = anuga_domain_and_engine()

    Baptist = Baptist_operator(domain,
        veg_diameter=c["D"],
        veg_density=c["m"],
        veg_height=c["hv"],
        bed_friction_const=c["C"],
    )

    # Test internal precomputed coefficients
    a1_expected = c["C"]**-2
    a2_expected = c["m"]*c["D"]/(2*c["g"])
    a3_expected = np.sqrt(c["g"])/c["K"]
    assert Baptist.a1 == approx(a1_expected, rel=1e-6)
    assert Baptist.a2 == approx(a2_expected, rel=1e-2)  # lower threshold in case g differs in sig figs
    assert Baptist.a3 == approx(a3_expected, rel=1e-1)  # lower threshold in case both g and K differ in sig figs


@mark.anuga
def test_baptist_set_vegetation(constants):
    from dycove.sim.engines.ANUGA_baptist import Baptist_operator

    c = constants
    domain, engine = anuga_domain_and_engine()
    
    Baptist = Baptist_operator(domain)

    # Check that external setting of vegetation works
    Baptist.set_vegetation(veg_diameter=c["D"], veg_density=c["m"], veg_height=c["hv"])
    
    assert np.all(domain.quantities['veg_diameter'].get_values() == approx(c["D"], rel=1e-6))
    assert np.all(domain.quantities['veg_density'].get_values() == approx(c["m"], rel=1e-6))
    assert np.all(domain.quantities['veg_height'].get_values() == approx(c["hv"], rel=1e-6))


@mark.anuga
def test_baptist_update_quantities_reduces_momentum(constants):
    from dycove.sim.engines.ANUGA_baptist import Baptist_operator

    c = constants
    domain, engine = anuga_domain_and_engine(momentum=True)

    Baptist = Baptist_operator(domain,
        veg_diameter=c["D"],
        veg_density=c["m"],
        veg_height=c["hv"],
        bed_friction_const=c["C"],
    )

    # Force deterministic timestep
    Baptist.dt = 0.1

    # Snapshot pre-update momenta (centroid arrays)
    xmom_before = domain.quantities["xmomentum"].centroid_values.copy()
    ymom_before = domain.quantities["ymomentum"].centroid_values.copy()

    # Masks where vegetation exists and depth > 0
    inds = (Baptist.veg_density.centroid_values > 0) & (Baptist.depth > 0.01)

    # Perform update
    Baptist.update_quantities(inds)

    xmom_after = domain.quantities["xmomentum"].centroid_values
    ymom_after = domain.quantities["ymomentum"].centroid_values

    # Momentum should decrease for vegetated cells
    assert np.all(xmom_after[inds] < xmom_before[inds])
    assert np.all(ymom_after[inds] < ymom_before[inds])

    # Non-vegetated cells should remain unchanged
    assert np.all(xmom_after[~inds] == approx(xmom_before[~inds]))
    assert np.all(ymom_after[~inds] == approx(ymom_before[~inds]))

    # Momentum should not drop below zero
    assert np.all(xmom_after >= 0)
    assert np.all(ymom_after >= 0)
