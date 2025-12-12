from pytest import mark
import sys, os



# @mark.dfm
# def test_dfm_dycove_serial(sample_data_path):
#     from dycove import DFM_hydro, VegetationSpecies, MultipleVegetationSpecies

#     veg_1 = VegetationSpecies(sample_data_path / "NL_ex.json", "veg1")
#     veg_2 = VegetationSpecies(sample_data_path / "CE_ex.json", "veg2")
#     multi = MultipleVegetationSpecies([veg_1, veg_2])

#     model = DFM_hydro.DFM(DFM_DLLs, config_file, mdu_file, vegetation=multi)

#     model.run_simulation(1.5, sim_time_unit="hydrodynamic days")