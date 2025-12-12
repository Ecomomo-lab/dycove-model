from pytest import fixture
from pathlib import Path

@fixture
def constants():
    # Vegetation test parameters (Baptist)
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
    return {
        "fraction_0"     : 0.5,
        "start_col_ets"  : 2,
        "end_col_ets"    : 3,
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

@fixture
def tmp_data_path():
    return Path(__file__).parent / "tmp_data"


def anuga_domain_and_engine(*, friction=False, momentum=False, dirichlet=False):
    import anuga
    from dycove import ANUGA_hydro

    # Minimal ANUGA domain object
    domain = anuga.rectangular_cross_domain(5, 5, 1000, 1000)

    # Set uniform depth and momentum for tests that check quantities
    domain.set_quantity("stage", 1.0)
    domain.set_quantity("elevation", 0.0)
    if friction:
        domain.set_quantity("friction", 0.3)
    if momentum:
        domain.set_quantity("xmomentum", 1.0)
        domain.set_quantity("ymomentum", 0.5)

    # Set all reflective boundaries for tests that evolve domain
    Br = anuga.Reflective_boundary(domain)
    if dirichlet:
        Bd_1 = anuga.Dirichlet_boundary([1.0, 0, 0])
        Bd_2 = anuga.Dirichlet_boundary([0.5, 0, 0])
        domain.set_boundary({'left': Bd_1, 'right': Bd_2, 'top': Br, 'bottom': Br})     
    else:   
        domain.set_boundary({'left': Br, 'right': Br, 'top': Br, 'bottom': Br})

    # Create ANUGA engine object
    engine = ANUGA_hydro.AnugaEngine(domain)

    return domain, engine
