import pytest
from unittest.mock import MagicMock, patch
from dataclasses import asdict
import numpy as np
import xarray as xr

from dycove.sim.outputs import OutputManager
from dycove.sim.vegetation_data import VegCohort


class TestOutputManager:

    @staticmethod
    def mock_engine(is_parallel=False, rank=0):
        engine = MagicMock()
        engine.veg = MagicMock()
        engine.model_dir = "/tmp/test_model"
        engine.is_parallel.return_value = is_parallel
        engine.get_rank.return_value = rank
        return engine


    @staticmethod
    def mock_cohort():
        """ Create full VegCohort object to test output writing """
        return VegCohort(
            name = "test_species",
            fraction=np.array([0.1, 0.2]),
            density = 100.0,
            diameter = 0.01,
            height = 0.5,
            rootlength = 0.3,
            lifestage = 1,
            lifestage_year = 2,
            potential_mort_flood = np.array([0.5, 0.1]),
            potential_mort_desic = np.array([0.5, 0.1]),
            potential_mort_uproot = np.array([0.5, 0.1]),
            potential_mort_burial = np.array([0.5, 0.1]),
            potential_mort_scour = np.array([0.5, 0.1]),
            applied_mort_flood = np.array([0.1, 0.1]),
            applied_mort_desic = np.array([0.1, 0.1]),
            applied_mort_uproot = np.array([0.1, 0.1]),
            applied_mort_burial = np.array([0.1, 0.1]),
            applied_mort_scour = np.array([0.1, 0.1]),
            applied_mort_total = np.array([0.1, 0.1]),
        )


    @staticmethod
    def mock_simstate(eco_year=1, ets=3):
        simstate = MagicMock()
        simstate.eco_year = eco_year
        simstate.ets = ets
        simstate.n_ets = 14
        simstate.veg_interval = 43200
        simstate.ecofac = 10
        return simstate


    def make_output_manager(self, engine, save_freq=1, save_mort=True):
        with patch("pathlib.Path.mkdir"):
            return OutputManager(engine, save_freq, save_mort)


    @pytest.mark.unit
    def test_update_file_counts_extends_with_new_cohorts(self):
        """ If cohorts outnumber tracked steps, n_cohort_steps is extended """
        engine = self.mock_engine()
        engine.veg.cohorts = [self.mock_cohort(), self.mock_cohort()]
        om = self.make_output_manager(engine)
        om.n_cohort_steps = [2]  # only tracking one cohort so far
        om.update_file_counts()
        assert len(om.n_cohort_steps) == 2


    @pytest.mark.unit
    def test_update_file_counts_no_change_when_counts_match(self):
        """ n_cohort_steps is not changed if already the right length """
        engine = self.mock_engine()
        engine.veg.cohorts = [self.mock_cohort()]
        om = self.make_output_manager(engine)
        om.n_cohort_steps = [3]
        om.update_file_counts()
        assert om.n_cohort_steps == [3]


    @pytest.mark.unit
    def test_cohort_indexing_creates_new_entry(self):
        """ First call for a year/ets combo creates nested dict and list """
        engine = self.mock_engine()
        om = self.make_output_manager(engine)
        om.fname_base = "cohort0_01"
        om.cohort_indexing(year=1, ets=3)
        assert om.cohort_index == {"1": {"3": ["cohort0_01"]}}


    @pytest.mark.unit
    def test_cohort_indexing_appends_to_existing_entry(self):
        """ Second cohort at same year/ets appends to existing list """
        engine = self.mock_engine()
        om = self.make_output_manager(engine)
        om.fname_base = "cohort0_01"
        om.cohort_indexing(year=1, ets=3)
        om.fname_base = "cohort1_01"
        om.cohort_indexing(year=1, ets=3)
        assert om.cohort_index["1"]["3"] == ["cohort0_01", "cohort1_01"]


    @pytest.mark.unit
    def test_save_vegetation_step_passes_args_to_save_netcdf(self):
        """ save_vegetation_step passes eco_year, ets, cohort_id to save_netcdf """
        engine = self.mock_engine()
        engine.veg.cohorts = [self.mock_cohort()]
        om = self.make_output_manager(engine)
        om.n_cohort_steps = [0]
        om.save_netcdf = MagicMock()

        simstate = self.mock_simstate(eco_year=2, ets=5)
        om.save_vegetation_step(simstate, vts=0)

        call_kwargs = om.save_netcdf.call_args.kwargs
        assert call_kwargs["eco_year"] == 2
        assert call_kwargs["ets"] == 5
        assert call_kwargs["cohort_id"] == 0
        assert "saved_attrs" not in call_kwargs


    @pytest.mark.unit
    def test_save_netcdf_excludes_mort_when_save_mort_false(self, tmp_path):
        engine = self.mock_engine()
        om = self.make_output_manager(engine, save_mort=False)
        data = {"height": 0.5, "mort_frac": np.array([0.1, 0.2])}

        om.save_netcdf(tmp_path, "test_file", data, eco_year=1, ets=1, cohort_id=0)

        ds = xr.open_dataset(tmp_path / "test_file.nc", engine="scipy")
        assert "mort_frac" not in ds.data_vars
        assert "height" in ds.attrs
        ds.close()


    @pytest.mark.unit
    def test_save_vegetation_step_skips_save_when_not_save_freq(self):
        """ When vts % save_freq != 0, save_netcdf is not called """
        engine = self.mock_engine()
        engine.veg.cohorts = [self.mock_cohort()]
        om = self.make_output_manager(engine, save_freq=2)
        om.n_cohort_steps = [0]
        om.save_netcdf = MagicMock()

        om.save_vegetation_step(self.mock_simstate(), vts=1)  # 1 % 2 != 0

        om.save_netcdf.assert_not_called()


    @pytest.mark.unit
    def test_save_vegetation_step_saves_when_save_freq(self):
        """ When vts % save_freq == 0, save_netcdf is called once per cohort """
        engine = self.mock_engine()
        engine.veg.cohorts = [self.mock_cohort(), self.mock_cohort()]
        om = self.make_output_manager(engine, save_freq=2)
        om.n_cohort_steps = [0, 0]
        om.save_netcdf = MagicMock()

        om.save_vegetation_step(self.mock_simstate(), vts=2)  # 2 % 2 == 0

        assert om.save_netcdf.call_count == 2


    @pytest.mark.unit
    def test_save_vegetation_step_uses_parallel_filename(self):
        """ In parallel mode, filenames must include processor rank """
        engine = self.mock_engine(is_parallel=True, rank=2)
        engine.veg.cohorts = [self.mock_cohort()]
        om = self.make_output_manager(engine)
        om.n_cohort_steps = [0]
        om.save_netcdf = MagicMock()

        om.save_vegetation_step(self.mock_simstate(), vts=0)

        saved_filename = om.save_netcdf.call_args[0][1]  # second positional arg
        assert "proc2" in saved_filename


    @pytest.mark.unit
    def test_n_cohort_steps_increments_regardless_of_save_freq(self):
        """ Step counter increments even when the save is skipped """
        engine = self.mock_engine()
        engine.veg.cohorts = [self.mock_cohort()]
        om = self.make_output_manager(engine, save_freq=2)
        om.n_cohort_steps = [0]
        om.save_netcdf = MagicMock()

        om.save_vegetation_step(self.mock_simstate(), vts=1)  # skipped

        assert om.n_cohort_steps[0] == 1  # still incremented


    @pytest.mark.unit
    def test_reconcile_calls_merge_when_parallel_and_rank0_and_veg_active(self):
        engine = self.mock_engine(is_parallel=True, rank=0)
        om = self.make_output_manager(engine)
        om.save_simulation_indices = MagicMock()

        om.reconcile_vegetation_output(self.mock_simstate())

        engine.merge_parallel_veg.assert_called_once()


    @pytest.mark.unit
    def test_reconcile_does_not_merge_when_not_parallel(self):
        engine = self.mock_engine(is_parallel=False, rank=0)
        om = self.make_output_manager(engine)
        om.save_simulation_indices = MagicMock()

        om.reconcile_vegetation_output(self.mock_simstate())

        engine.merge_parallel_veg.assert_not_called()


    @pytest.mark.unit
    def test_reconcile_does_not_merge_when_parallel_but_not_rank_0(self):
        """ Non-root processors should not trigger a merge """
        engine = self.mock_engine(is_parallel=True, rank=1)
        om = self.make_output_manager(engine)
        om.save_simulation_indices = MagicMock()

        om.reconcile_vegetation_output(self.mock_simstate())

        engine.merge_parallel_veg.assert_not_called()


    @pytest.mark.unit
    def test_reconcile_does_not_merge_when_not_veg_active(self):
        """ Non-vegetation simulation should not trigger anything (b/c no files) """
        engine = self.mock_engine(is_parallel=True, rank=0)
        om = self.make_output_manager(engine)
        om.save_simulation_indices = MagicMock()

        om.veg = None
        om.reconcile_vegetation_output(self.mock_simstate())

        engine.merge_parallel_veg.assert_not_called()


    @pytest.mark.unit
    def test_reconcile_passes_self_to_merge_parallel(self):
        """ merge_parallel_veg requires the OutputManager instance as argument """
        engine = self.mock_engine(is_parallel=True, rank=0)
        om = self.make_output_manager(engine)
        om.save_simulation_indices = MagicMock()

        om.reconcile_vegetation_output(self.mock_simstate())


        assert engine.merge_parallel_veg.call_args[0][0] is om


    @pytest.mark.unit
    def test_save_netcdf_field_names_match_vegcohort(self, tmp_path):
        """ Field names in saved NetCDF match VegCohort, catching any renames that break post-processing """
        engine = self.mock_engine()
        om = self.make_output_manager(engine)

        cohort = self.mock_cohort()
        om.save_netcdf(tmp_path, "test_cohort", asdict(cohort), eco_year=1, ets=1, cohort_id=0)

        ds = xr.open_dataset(tmp_path / "test_cohort.nc", engine="scipy")

        # Arrays should be data variables
        expected_data_vars = {"fraction", 
                              "potential_mort_flood", 
                              "potential_mort_desic",
                              "potential_mort_uproot",
                              "potential_mort_burial",
                              "potential_mort_scour",
                              "applied_mort_flood",
                              "applied_mort_desic",
                              "applied_mort_uproot",
                              "applied_mort_scour",
                              "applied_mort_total"}
        assert expected_data_vars.issubset(set(ds.data_vars))

        # Scalars should be attributes
        expected_attrs = {"name", "density", "diameter", "height", "rootlength", "lifestage", "lifestage_year"}
        assert expected_attrs.issubset(set(ds.attrs))

        ds.close()