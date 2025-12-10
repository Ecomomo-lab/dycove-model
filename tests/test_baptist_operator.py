
import numpy as np
from pytest import approx
from dycove.sim.engines.ANUGA_baptist import Baptist_operator
import anuga

def test_baptist_cv_computation():
    # minimal fake ANUGA domain object
    class FakeDomain: pass
    fake = FakeDomain()

    ht, diam, dens, C = 0.5, 0.005, 100, 65
    g, K = 9.81, 0.41

    op = Baptist_operator(
        domain=fake,
        veg_diameter=diam,
        veg_density=dens,
        veg_height=ht,
        bed_friction_const=C,
    )

    # test internal precomputed coefficients
    assert op.a1 == approx(C**-2)
    assert op.a2 == approx(dens*diam/(2*g), rel=1e-2)
    assert op.a3 == approx(np.sqrt(g)/K, rel=1e-1)
