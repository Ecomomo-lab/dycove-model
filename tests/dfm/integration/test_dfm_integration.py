from pytest import mark
import shutil



@mark.dfm
@mark.integration
def test_dfm_dycove_serial(sample_data_path, dfm_dll_path, tmp_path):
    from dycove import DFM_hydro, VegetationSpecies, MultipleVegetationSpecies

    config_file = "dimr_config.xml"
    vegfile_1 = "NL_ex.json"
    vegfile_2 = "CE_ex.json"

    run_dir = tmp_path / "run"
    run_dir.mkdir()

    shutil.copytree(sample_data_path / "dflowfm", run_dir / "dflowfm")
    shutil.copy(sample_data_path / config_file, run_dir)
    shutil.copy(sample_data_path / vegfile_1, run_dir)
    shutil.copy(sample_data_path / vegfile_2, run_dir)

    veg_1 = VegetationSpecies(run_dir / vegfile_1)
    veg_2 = VegetationSpecies(run_dir / vegfile_2)

    model = DFM_hydro.DFM(dfm_dll_path, 
                          run_dir / config_file, 
                          run_dir / "dflowfm/FlowFM.mdu", 
                          vegetation=MultipleVegetationSpecies([veg_1, veg_2]))

    model.run_simulation(1, sim_time_unit="hydrodynamic days")
