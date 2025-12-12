from pytest import mark
import sys, os

sys.path.append(os.path.realpath(os.path.dirname(__file__)+"/.."))
from conftest import anuga_domain_and_engine


@mark.anuga
def test_anuga_dycove_serial(sample_data_path):
    from dycove import ANUGA_hydro, VegetationSpecies, MultipleVegetationSpecies

    domain, engine = anuga_domain_and_engine(friction=True, dirichlet=True)

    veg_1 = VegetationSpecies(sample_data_path / "NL_ex.json", "veg1")
    veg_2 = VegetationSpecies(sample_data_path / "CE_ex.json", "veg2")
    multi = MultipleVegetationSpecies([veg_1, veg_2])

    model = ANUGA_hydro.ANUGA(domain, vegetation=multi)

    model.run_simulation(1.5, sim_time_unit="hydrodynamic days")

