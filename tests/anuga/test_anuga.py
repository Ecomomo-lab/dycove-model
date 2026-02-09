from pytest import mark
import sys, os
import shutil

sys.path.append(os.path.realpath(os.path.dirname(__file__)+"/.."))
from conftest import make_anuga_domain_and_engine as make_anuga


@mark.anuga
def test_anuga_dycove_serial(sample_data_path, tmp_data_path):
    from dycove import ANUGA_hydro, VegetationSpecies, MultipleVegetationSpecies

    domain, engine = make_anuga(friction=True, dirichlet=True)

    veg_1 = VegetationSpecies(sample_data_path / "NL_ex.json", "veg1")
    veg_2 = VegetationSpecies(sample_data_path / "CE_ex.json", "veg2")
    multi = MultipleVegetationSpecies([veg_1, veg_2])

    tmp_data_path.mkdir(exist_ok=True)
    domain.set_datadir(str(tmp_data_path))

    model = ANUGA_hydro.ANUGA(domain, vegetation=multi)

    model.run_simulation(1.5, sim_time_unit="hydrodynamic days")

    shutil.rmtree(tmp_data_path)

