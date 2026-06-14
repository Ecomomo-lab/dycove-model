from pytest import approx, mark
import shutil
import sys, os

sys.path.append(os.path.realpath(os.path.dirname(__file__)+"/../.."))
from conftest import make_anuga_domain_and_engine as make_anuga
    
    
@mark.anuga
@mark.integration
def test_anuga_dycove_serial(sample_data_path, tmp_path):
    from dycove import ANUGA_hydro, VegetationSpecies, MultipleVegetationSpecies

    domain, engine = make_anuga(friction=True, tide=True, datadir=tmp_path)

    veg_1 = VegetationSpecies(sample_data_path / "NL_ex.json", "veg1")
    veg_2 = VegetationSpecies(sample_data_path / "CE_ex.json", "veg2")
    multi = MultipleVegetationSpecies([veg_1, veg_2])

    model = ANUGA_hydro.ANUGA(domain, vegetation=multi)

    model.run_simulation(1, sim_time_unit="hydrodynamic days")


@mark.anuga
@mark.integration
def test_step(tmp_path):
    domain, engine = make_anuga(friction=True, tide=True, datadir=tmp_path)

    # Timestamp and step model by small dt
    t0 = domain.get_time()
    dt = 1
    engine.step(dt)

    # Check domain time state
    assert domain.get_time() - dt == approx(t0, rel=1e-6)


@mark.anuga
@mark.integration
def test_merge_parallel_veg(sample_data_path, tmp_path):
    
    from dycove.sim.outputs import OutputManager

    domain, engine = make_anuga(datadir=tmp_path)

    # Copy existing sample output files so we can place them back cleanly afterwards
    shutil.copytree(sample_data_path / "veg_output", tmp_path / "veg_output")
    shutil.copy(sample_data_path / "domain_P2_0.sww", tmp_path)
    shutil.copy(sample_data_path / "domain_P2_1.sww", tmp_path)

    # So that merge_parallel_veg() can find the files
    domain.set_name("domain_P2_0")  # this is the specific name for myid 0, so it has suffix "_0"
    engine.model_dir = domain.get_datadir()
    engine.numprocs = 2

    # Normally created upon creation of HydroSimulationBase object
    om = OutputManager(engine)

    # Load example cohort file index
    om.file_index = om.load_cohort_file_index()

    # Shortcutting to this method allows bypass of creating VegetationSpecies object
    engine.merge_parallel_veg(om)
