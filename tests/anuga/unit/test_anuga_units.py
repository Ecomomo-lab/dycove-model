import pytest
from unittest.mock import MagicMock, patch
import numpy as np
import xarray as xr
import sys, os

sys.path.append(os.path.realpath(os.path.dirname(__file__)+"/../.."))
#from conftest import make_anuga_domain_and_engine as make_anuga
from dycove.sim.engines.ANUGA_hydro import ANUGA, AnugaEngine


@pytest.mark.anuga
@pytest.mark.unit
def test_anuga_hydro_components():
    assert callable(ANUGA)
    assert callable(AnugaEngine)


@pytest.mark.anuga
@pytest.mark.unit
def test_engine_abstractmethods(required_engine_methods):
    # Check for methods defined by abstract class base.HydroEngineBase
    for method in required_engine_methods:
        assert hasattr(AnugaEngine, method)


class TestAnugaEngine:

    @staticmethod
    def mock_domain():
        domain = MagicMock()
        domain.get_datadir.return_value = "/tmp/test_model"
        domain.evolve.return_value = iter([0])  # simulate one yieldstep
        return domain
    

    @staticmethod
    def mock_veg():
        veg = MagicMock()
        veg.get_drag.return_value = 1.0
        return veg


    @staticmethod
    def make_engine(domain, vegetation=None, save_interval=3600, 
                    myid=1, nprocs=2,  # use as default to avoid r.report call under step()
                    barrier=MagicMock(),  # one test will implement something here
                    elevation=np.array([0.0, 0.0]),
                    ):
        with patch("dycove.sim.engines.ANUGA_hydro._import_anuga",
                return_value=(myid, nprocs, MagicMock(), barrier)):
            engine = AnugaEngine(domain, vegetation=vegetation)
        engine.save_interval = save_interval
        engine.get_elevation = MagicMock(return_value=elevation)
        return engine


    @pytest.mark.anuga
    @pytest.mark.unit
    def test_anuga_obj_has_anuga_engine_attrs_after_init(self):
        domain = self.mock_domain()

        # patch separately here because testing ANUGA class not AnugaEngine directly
        with patch("dycove.sim.engines.ANUGA_hydro._import_anuga", 
                return_value=(1, 2, MagicMock(), MagicMock())):
            model = ANUGA(domain)

        assert isinstance(model.engine, AnugaEngine)
        assert model.engine.domain is domain


    @pytest.mark.anuga
    @pytest.mark.unit
    def test_initialize_sets_morphology_false(self):
        """ morphology is always set to False regardless of vegetation """
        domain = self.mock_domain()
        engine = self.make_engine(domain)
        engine.initialize()
        assert engine.morphology == False


    @pytest.mark.anuga
    @pytest.mark.unit
    def test_initialize_creates_baptist_operator_when_veg_present(self):
        """ Baptist_operator is constructed with drag from veg.get_drag() when veg is set """
        domain = self.mock_domain()
        veg = self.mock_veg()
        engine = self.make_engine(domain, vegetation=veg)

        with patch("dycove.sim.engines.ANUGA_hydro.Baptist_operator") as mock_baptist:
            engine.initialize()

        mock_baptist.assert_called_once_with(domain, 
                                             veg_diameter=0, veg_density=0, veg_height=0, 
                                             drag=1.0)


    @pytest.mark.anuga
    @pytest.mark.unit
    def test_initialize_skips_baptist_operator_when_no_veg(self):
        """ Baptist_operator is not created when vegetation is None """
        domain = self.mock_domain()
        engine = self.make_engine(domain)

        with patch("dycove.sim.engines.ANUGA_hydro.Baptist_operator") as mock_baptist:
            engine.initialize()

        mock_baptist.assert_not_called()


    @pytest.mark.anuga
    @pytest.mark.unit
    def test_step_yieldstep_is_min_of_seconds_and_save_interval(self):
        """ yieldstep passed to domain.evolve is capped at save_interval """
        domain = self.mock_domain()
        engine = self.make_engine(domain)

        engine.step(seconds=43200)

        call_kwargs = domain.evolve.call_args[1]
        assert call_kwargs["yieldstep"] == 3600  # capped at save_interval, not 43200


    @pytest.mark.anuga
    @pytest.mark.unit
    def test_step_yieldstep_smaller_than_save_interval(self):
        """ yieldstep equals `seconds` when seconds < save_interval """
        domain = self.mock_domain()
        engine = self.make_engine(domain)

        engine.step(seconds=900)

        call_kwargs = domain.evolve.call_args[1]
        assert call_kwargs["yieldstep"] == 900


    @pytest.mark.anuga
    @pytest.mark.unit
    def test_step_skip_step_is_true_after_first_call(self):
        """ skip_step is False on first call but True after, preventing duplicate yieldsteps """
        domain = self.mock_domain()
        engine = self.make_engine(domain)

        assert engine.skip_step == False
        engine.step(seconds=3600)
        assert engine.skip_step == True


    @pytest.mark.anuga
    @pytest.mark.unit
    def test_step_passes_skip_step_to_evolve(self):
        """ domain.evolve receives the current value of skip_step """
        domain = self.mock_domain()
        engine = self.make_engine(domain)

        engine.step(seconds=3600)

        call_kwargs = domain.evolve.call_args[1]
        assert call_kwargs["skip_initial_step"] == False


    @pytest.mark.anuga
    @pytest.mark.unit
    def test_step_calls_barrier_twice(self):
        """ barrier() is called at start and end to synchronize parallel processes """
        domain = self.mock_domain()

        domain.evolve.side_effect = lambda **kwargs: call_order.append("evolve") or iter([0])
        barrier = MagicMock(side_effect=lambda: call_order.append("barrier"))
                
        engine = self.make_engine(domain, barrier=barrier)

        call_order = []
        engine.step(seconds=3600)

        assert call_order == ["barrier", "evolve", "barrier"]


    @pytest.mark.anuga
    @pytest.mark.unit
    def test_cleanup_does_parallel_ops_when_parallel(self):
        """ cleanup() calls finalization methods when parallel """
        domain = self.mock_domain()
        engine = self.make_engine(domain)  # parallel arguments by default
        domain.sww_merge = MagicMock()

        engine.cleanup()
            
        domain.sww_merge.assert_called_once()
        engine.finalize.assert_called_once()


    @pytest.mark.anuga
    @pytest.mark.unit
    def test_get_cell_count_equals_len_elevation_array(self):
        """ get_cell_count should call get_elevation corrcetly """
        domain = self.mock_domain()
        engine = self.make_engine(domain)  # parallel arguments by default

        n_cells = engine.get_cell_count()

        assert n_cells == 2  # based on length of MagicMock array in make_engine


    @pytest.mark.anuga
    @pytest.mark.unit
    def test_get_velocity_and_depth_computes_depth_correctly(self):
        """ depth = stage - elevation """
        domain = self.mock_domain()
        engine = self.make_engine(domain, elevation=np.array([1.0, 2.0]))
        domain.quantities = {
            "stage"     : MagicMock(centroid_values=np.array([3.0, 4.0])),
            "xmomentum" : MagicMock(centroid_values=np.array([0.0, 0.0])),
            "ymomentum" : MagicMock(centroid_values=np.array([0.0, 0.0])),
        }
        _, depth = engine.get_velocity_and_depth()
        assert np.allclose(depth, np.array([2.0, 2.0]))


    @pytest.mark.anuga
    @pytest.mark.unit
    def test_get_velocity_and_depth_computes_velocity_magnitude(self):
        """ velocity = sqrt(xvel^2 + yvel^2) from momentum/depth """
        domain = self.mock_domain()
        engine = self.make_engine(domain, elevation=np.array([0.0]))
        domain.quantities = {
            "stage"     : MagicMock(centroid_values=np.array([2.0])),
            "xmomentum" : MagicMock(centroid_values=np.array([3.0])),  # xvel = 3/2 = 1.5
            "ymomentum" : MagicMock(centroid_values=np.array([4.0])),  # yvel = 4/2 = 2.0
        }
        velocity, _ = engine.get_velocity_and_depth()
        assert np.allclose(velocity, np.array([2.5]))  # sqrt(1.5^2 + 2.0^2)


    @pytest.mark.anuga
    @pytest.mark.unit
    def test_get_velocity_and_depth_zeros_velocity_below_depth_limit(self):
        """ Velocity is set to zero where depth < H_LIM_VELOCITY """
        domain = self.mock_domain()
        engine = self.make_engine(domain, elevation=np.array([0.0, 0.0]))
        domain.quantities = {
            "stage"     : MagicMock(centroid_values=np.array([0.0001, 2.0])),  # first cell below H_LIM
            "xmomentum" : MagicMock(centroid_values=np.array([1.0,    1.0])),
            "ymomentum" : MagicMock(centroid_values=np.array([1.0,    1.0])),
        }
        velocity, _ = engine.get_velocity_and_depth()
        assert velocity[0] == 0.0       # zeroed out
        assert velocity[1] > 0.0       # normal cell unaffected


    @pytest.mark.anuga
    @pytest.mark.unit
    def test_get_vegetation_returns_correct_baptist_arrays(self):
        """ stemdensity, stemdiameter, stemheight are pulled from the correct Baptist attributes """
        domain = self.mock_domain()
        engine = self.make_engine(domain)
        engine.Baptist = MagicMock()
        engine.Baptist.veg_density.centroid_values  = np.array([100.0])
        engine.Baptist.veg_diameter.centroid_values = np.array([0.01])
        engine.Baptist.veg_height.centroid_values   = np.array([0.5])

        stemdensity, stemdiameter, stemheight = engine.get_vegetation()

        assert np.allclose(stemdensity,  np.array([100.0]))
        assert np.allclose(stemdiameter, np.array([0.01]))
        assert np.allclose(stemheight,   np.array([0.5]))


    @pytest.mark.anuga
    @pytest.mark.unit
    def test_set_vegetation_passes_correct_keyword_arguments_to_baptist(self):
        """Arguments are passed to Baptist.set_vegetation with the correct keyword names."""
        domain = self.mock_domain()
        engine = self.make_engine(domain)
        engine.Baptist = MagicMock()

        engine.set_vegetation(stemdensity=np.array([100.0]),
                              stemdiameter=np.array([0.01]),
                              stemheight=np.array([0.5]))

        engine.Baptist.set_vegetation.assert_called_once_with(
            veg_diameter=np.array([0.01]),
            veg_density=np.array([100.0]),
            veg_height=np.array([0.5]),
        )



class TestMergeLocalToGlobal:
    """ @staticmethod, no mocking needed """

    @pytest.mark.anuga
    @pytest.mark.unit
    def test_basic_merge_maps_local_to_global_indices(self):
        """ Values are placed at correct global indices via f_gids mapping """
        local_data = [
            {"fraction": np.array([0.1, 0.2, 0.3])},  # proc 0
            {"fraction": np.array([0.4, 0.5, 0.6])},  # proc 1
        ]
        f_ids  = [np.array([0, 1, 2]), np.array([0, 1, 2])]  # all local cells are "full"
        f_gids = [np.array([0, 1, 2]), np.array([3, 4, 5])]  # proc 0 owns globals 0-2, proc 1 owns 3-5
        n_global = 6

        merged = AnugaEngine.merge_local_to_global(local_data, f_gids, f_ids, n_global)

        assert np.allclose(merged["fraction"], [0.1, 0.2, 0.3, 0.4, 0.5, 0.6])


    @pytest.mark.anuga
    @pytest.mark.unit
    def test_merge_only_uses_full_non_ghost_cells(self):
        """ Only cells where tri_full_flag == 1 (f_ids) contribute to the global array """
        local_data = [
            {"fraction": np.array([0.1, 0.2, 0.3])},  # proc 0 has 3 cells but only 2 are "full"
        ]
        f_ids  = [np.array([0, 2])]  # cells 0 and 2 are full, cell 1 is ghost
        f_gids = [np.array([0, 1])]  # full cells map to global 0 and 1
        n_global = 2

        merged = AnugaEngine.merge_local_to_global(local_data, f_gids, f_ids, n_global)

        assert np.allclose(merged["fraction"], [0.1, 0.3])  # cell 1 (0.2) excluded


    @pytest.mark.anuga
    @pytest.mark.unit
    def test_merge_handles_multiple_variables(self):
        """ All variables in local_data are merged into the output dict """
        local_data = [{"fraction": np.array([0.5]), "height": np.array([1.0])}]
        f_ids  = [np.array([0])]
        f_gids = [np.array([0])]
        n_global = 1

        merged = AnugaEngine.merge_local_to_global(local_data, f_gids, f_ids, n_global)

        assert set(merged.keys()) == {"fraction", "height"}


    @pytest.mark.anuga
    @pytest.mark.unit
    def test_merge_output_array_length_equals_n_global(self):
        """ Output arrays have length n_global regardless of local array sizes """
        local_data = [{"fraction": np.array([0.5])}]
        f_ids  = [np.array([0])]
        f_gids = [np.array([2])]  # maps to global index 2
        n_global = 5

        merged = AnugaEngine.merge_local_to_global(local_data, f_gids, f_ids, n_global)

        assert len(merged["fraction"]) == 5


class TestMergeParallelVeg:    

    @staticmethod
    def make_engine(tmp_path, numprocs=2, n_cohorts=1):
        domain = MagicMock()
        domain.get_name.return_value = f"test_domain_P{numprocs}_0"  # 0 for myid == 0
        domain.get_datadir.return_value = str(tmp_path)

        with patch("dycove.sim.engines.ANUGA_hydro._import_anuga",
                   return_value=(0, numprocs, MagicMock(), MagicMock())), \
             patch("dycove.sim.engines.ANUGA_hydro.r"):
            engine = AnugaEngine(domain)

        engine.veg = MagicMock()
        engine.veg.cohorts = [MagicMock()] * n_cohorts
        engine.numprocs = numprocs
        return engine


    @staticmethod
    def write_sww_file(path, tri_l2g, tri_full_flag):
        """ Write a minimal fake .sww file with just the fields merge_parallel_veg needs """
        ds = xr.Dataset({
            "tri_l2g"       : xr.DataArray(tri_l2g),
            "tri_full_flag" : xr.DataArray(tri_full_flag),
        })
        ds.to_netcdf(path)


    @staticmethod
    def write_cohort_file(path, fraction, attrs=None):
        """ Write a minimal vegetation cohort NetCDF file """
        ds = xr.Dataset({"fraction": xr.DataArray(fraction)})
        if attrs:
            ds.attrs.update(attrs)
        ds.to_netcdf(path, engine="scipy")


    @staticmethod
    def mock_output_manager(tmp_path, n_cohort_steps=[1]):
        om = MagicMock()
        om.veg_dir = tmp_path
        om.n_cohort_steps = n_cohort_steps
        return om


    @pytest.mark.anuga
    @pytest.mark.unit
    def test_raises_filenotfounderror_when_cohort_file_missing(self, tmp_path):
        """ FileNotFoundError is raised when a processor cohort file does not exist """
        engine = self.make_engine(tmp_path)

        # Write SWW files
        for p in range(engine.numprocs):
            self.write_sww_file(
                tmp_path / f"test_domain_P{engine.numprocs}_{p}.sww",
                tri_l2g=np.array([0, 1]),
                tri_full_flag=np.array([1, 1])
            )

        # Write only proc 0 cohort file, not proc 1
        self.write_cohort_file(tmp_path / "cohort0_00_proc0.nc", fraction=np.array([0.1, 0.2]))

        om = self.mock_output_manager(tmp_path)

        with patch("dycove.sim.engines.ANUGA_hydro.r"):
            with pytest.raises(FileNotFoundError):
                engine.merge_parallel_veg(om)


    @pytest.mark.anuga
    @pytest.mark.unit
    def test_merged_file_created(self, tmp_path):
        """ A merged output NetCDF file is created after merging processor files """
        engine = self.make_engine(tmp_path)

        for p in range(engine.numprocs):
            self.write_sww_file(
                tmp_path / f"test_domain_P{engine.numprocs}_{p}.sww",
                tri_l2g=np.array([p]),       # proc 0 owns global 0, proc 1 owns global 1
                tri_full_flag=np.array([1])
            )
            self.write_cohort_file(tmp_path / f"cohort0_00_proc{p}.nc",
                                   fraction=np.array([0.1 * (p + 1)]))

        om = self.mock_output_manager(tmp_path)
        om.save_netcdf = MagicMock()

        with patch("dycove.sim.engines.ANUGA_hydro.r"):
            engine.merge_parallel_veg(om)

        om.save_netcdf.assert_called_once()


    @pytest.mark.anuga
    @pytest.mark.unit
    def test_processor_files_deleted_after_merge(self, tmp_path):
        """ Individual processor cohort files are removed after merging """
        engine = self.make_engine(tmp_path)

        for p in range(engine.numprocs):
            self.write_sww_file(
                tmp_path / f"test_domain_P{engine.numprocs}_{p}.sww",
                tri_l2g=np.array([p]),
                tri_full_flag=np.array([1])
            )
            self.write_cohort_file(tmp_path / f"cohort0_00_proc{p}.nc",
                                   fraction=np.array([0.5]))

        om = self.mock_output_manager(tmp_path)
        om.save_netcdf = MagicMock()

        with patch("dycove.sim.engines.ANUGA_hydro.r"):
            engine.merge_parallel_veg(om)

        assert not (tmp_path / "cohort0_00_proc0.nc").exists()
        assert not (tmp_path / "cohort0_00_proc1.nc").exists()


    @pytest.mark.anuga
    @pytest.mark.unit
    def test_attrs_taken_from_proc0_file(self, tmp_path):
        """ Merged file attrs come from processor 0's cohort file """
        engine = self.make_engine(tmp_path)

        for p in range(engine.numprocs):
            self.write_sww_file(
                tmp_path / f"test_domain_P{engine.numprocs}_{p}.sww",
                tri_l2g=np.array([p]),
                tri_full_flag=np.array([1])
            )
            self.write_cohort_file(tmp_path / f"cohort0_00_proc{p}.nc",
                                   fraction=np.array([0.5]),
                                   attrs={"eco_year": 1, "ets": 3, "cohort": 0})

        om = self.mock_output_manager(tmp_path)
        om.save_netcdf = MagicMock()

        with patch("dycove.sim.engines.ANUGA_hydro.r"):
            engine.merge_parallel_veg(om)

        saved_attrs = om.save_netcdf.call_args[1]["saved_attrs"]
        assert saved_attrs["eco_year"] == 1
        assert saved_attrs["ets"] == 3


    @pytest.mark.anuga
    @pytest.mark.unit
    def test_merged_values_correctly_assembled_from_processors(self, tmp_path):
        """ Global array contains correct values from each processor at the right indices """
        engine = self.make_engine(tmp_path)

        # proc 0 owns global cells 0,1 — proc 1 owns global cells 2,3
        for p in range(engine.numprocs):
            self.write_sww_file(
                tmp_path / f"test_domain_P{engine.numprocs}_{p}.sww",
                tri_l2g=np.array([p*2, p*2+1]),
                tri_full_flag=np.array([1, 1])
            )
            self.write_cohort_file(tmp_path / f"cohort0_00_proc{p}.nc",
                                   fraction=np.array([0.1*(p*2+1), 0.1*(p*2+2)]))

        om = self.mock_output_manager(tmp_path)
        om.save_netcdf = MagicMock()

        with patch("dycove.sim.engines.ANUGA_hydro.r"):
            engine.merge_parallel_veg(om)

        merged_data = om.save_netcdf.call_args[0][2]  # third positional arg is the data dict
        assert np.allclose(merged_data["fraction"], [0.1, 0.2, 0.3, 0.4])


    @pytest.mark.anuga
    @pytest.mark.unit
    def test_save_netcdf_has_correct_filename(self, tmp_path):
        """ Global array contains correct values from each processor at the right indices """
        engine = self.make_engine(tmp_path)

        # proc 0 owns global cells 0,1 — proc 1 owns global cells 2,3
        for p in range(engine.numprocs):
            self.write_sww_file(
                tmp_path / f"test_domain_P{engine.numprocs}_{p}.sww",
                tri_l2g=np.array([p]),
                tri_full_flag=np.array([1])
            )
            self.write_cohort_file(tmp_path / f"cohort0_00_proc{p}.nc",
                                   fraction=np.array([0.5]))

        om = self.mock_output_manager(tmp_path)
        om.save_netcdf = MagicMock()

        with patch("dycove.sim.engines.ANUGA_hydro.r"):
            engine.merge_parallel_veg(om)

        fname_stem = om.save_netcdf.call_args[0][1]  # third positional arg is the data dict
        assert fname_stem == "cohort0_00"