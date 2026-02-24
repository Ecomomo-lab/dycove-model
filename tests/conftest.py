from pytest import fixture
from pathlib import Path
import numpy as np


@fixture
def dfm_dll_path():
    # Update this if necessary/desired, to run Delft3D FM tests
    return Path(
        "C:/Program Files (x86)/"
        "Deltares/"
        "Delft3D Flexible Mesh Suite HM (2021.03)/"
        "plugins/"
        "DeltaShell.Dimr/"
        "kernels/"
        "x64"
        )


@fixture
def array_len():
    return 100


@fixture
def constants():
    # Vegetation test parameters (Baptist)
    # Don't update these; many tests depend on them
    return {
        "D" : 0.01,
        "m" : 100,
        "hv": 0.5,
        "C" : 65,
        "g" : 9.81,
        "K" : 0.41,
    }


@fixture
def attributes():
    # Vegetation test parameters (attributes)
    # Don't update these; many tests depend on them
    return {
        "fraction_0"     : 0.5,
        "start_col_ets"  : 2,
        "end_col_ets"    : 4,
        "rootlength_0"   : 0.5,
        "nls"            : 1,
        "flood_no_mort"  : 0.4,
        "flood_all_mort" : 0.6,
        "desic_no_mort"  : 0.85,
        "desic_all_mort" : 0.95,
        "uproot_no_mort" : 1.0,
        "uproot_all_mort": 1.5,
        "start_growth_ets": 4,
        "end_growth_ets": 8,
        "winter_ets": 10,
        "stemht_max": 1.0,
        "ht_growth_rates": 0.125,
        "stemht_winter_max": 0.1,
    }


@fixture
def required_engine_methods():
    return (
        "initialize",
        "step",
        "cleanup",
        "get_rank",
        "get_cell_count",
        "get_elevation",
        "get_velocity_and_depth",
        "get_vegetation",
        "set_vegetation",
        "check_simulation_inputs",
        "is_parallel",
        "merge_parallel_veg",
    )


@fixture
def sample_data_path():
    return Path(__file__).parent / "sample_data"


def make_anuga_domain_and_engine(*, friction=False, momentum=False, tide=False, datadir=None):
    import anuga
    from dycove import ANUGA_hydro

    # Minimal ANUGA domain object (very large elements to speed up time steps)
    domain = anuga.rectangular_cross_domain(5, 5, 100, 100)
    domain.set_datadir(str(datadir))

    # Set uniform depth and momentum for tests that check quantities
    domain.set_quantity("stage", 0.5)
    domain.set_quantity("elevation", 0.0)
    if friction:
        domain.set_quantity("friction", 0.3)
    if momentum:
        domain.set_quantity("xmomentum", 0.2)
        domain.set_quantity("ymomentum", 0.1)

    # Set all reflective boundaries for tests that evolve domain
    Br = anuga.Reflective_boundary(domain)
    if tide:
        def tide_func(t, A=0.5, T=43200):
            return A*np.cos(2*np.pi * t/T)
        Bt = anuga.Time_boundary(domain, function=lambda t: [tide_func(t), 0.0, 0.0])
        domain.set_boundary({'left': Bt, 'right': Bt, 'top': Br, 'bottom': Br}) 
    else:   
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

    # Create ANUGA engine object
    engine = ANUGA_hydro.AnugaEngine(domain)

    return domain, engine
