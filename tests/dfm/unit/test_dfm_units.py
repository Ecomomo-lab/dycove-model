import numpy as np
import pytest
import xarray as xr

from unittest.mock import MagicMock
from pytest import mark


@mark.dfm
@mark.unit
def test_dfm_engine_imports():
    from dycove import DFM_hydro
    assert callable(DFM_hydro.DFM)
    assert callable(DFM_hydro.DFMEngine)


@mark.dfm
@mark.unit
def test_engine_abstractmethods(required_engine_methods):
    from dycove import DFM_hydro
    for method in required_engine_methods:
        assert hasattr(DFM_hydro.DFMEngine, method)


class TestDFMParallelVegetationMerge:
    """
    Unit tests for reconstruction of rank-local DYCOVE vegetation
    output onto the native merged DFM map ordering.
    """

    @staticmethod
    def make_engine(tmp_path):
        from dycove import DFM_hydro

        # Avoid DFM/BMI initialization; construct only the state required
        # by merge_parallel_veg().
        engine = object.__new__(DFM_hydro.DFMEngine)

        engine.model_dir = tmp_path
        engine.mdu_path = tmp_path / "FlowFM.mdu"

        # Local BMI order on each rank:
        #
        # rank 0:
        #   gid 3 = owned
        #   gid 2 = ghost
        #
        # rank 1:
        #   gid 1 = owned
        #   gid 2 = owned
        #   gid 3 = ghost
        #
        # Complete physical domain = gids 1, 2, 3.
        engine.parallel_partition_maps = [
            (
                np.array([3, 2], dtype=np.int64),
                np.array([True, False]),
            ),
            (
                np.array([1, 2, 3], dtype=np.int64),
                np.array([True, True, False]),
            ),
        ]

        # merge_parallel_veg() is a rank-0 operation.
        engine.get_rank = MagicMock(return_value=0)

        return engine


    @staticmethod
    def write_merged_map(tmp_path):
        output_dir = tmp_path / "output"
        output_dir.mkdir()

        # Deliberately NOT ascending global-ID order.
        # This verifies vegetation follows the native merged DFM map,
        # rather than assuming output_index == global_id - 1.
        ds = xr.Dataset({
            "mesh2d_flowelem_globalnr": xr.DataArray(
                np.array([3, 1, 2], dtype=np.int64)
            )
        })

        ds.to_netcdf(
            output_dir / "FlowFM_map.nc",
            engine="scipy",
        )


    @staticmethod
    def write_cohort_file(path, fraction, attrs=None):
        ds = xr.Dataset({
            "fraction": xr.DataArray(
                np.asarray(fraction, dtype=float)
            )
        })

        if attrs:
            ds.attrs.update(attrs)

        ds.to_netcdf(path, engine="scipy")


    @staticmethod
    def make_output_manager(tmp_path):
        om = MagicMock()

        om.veg_dir = tmp_path / "veg_output"
        om.veg_dir.mkdir()

        om.file_index = {
            "1": {
                "2": ["cohort0_000"]
            }
        }

        return om


    @pytest.mark.dfm
    @pytest.mark.unit
    def test_merge_uses_owned_cells_and_native_dfm_map_order(
        self,
        tmp_path,
    ):
        """
        Vegetation is reconstructed from owned cells only and written
        in the exact global-cell order of the native merged DFM map.
        """
        engine = self.make_engine(tmp_path)
        self.write_merged_map(tmp_path)

        om = self.make_output_manager(tmp_path)

        # rank 0 local BMI order = [gid3 owned, gid2 ghost]
        self.write_cohort_file(
            om.veg_dir / "cohort0_000_proc0.nc",
            fraction=[0.30, 9.99],
            attrs={"eco_year": 1, "ets": 2, "cohort": 0},
        )

        # rank 1 local BMI order = [gid1 owned, gid2 owned, gid3 ghost]
        self.write_cohort_file(
            om.veg_dir / "cohort0_000_proc1.nc",
            fraction=[0.10, 0.20, 8.88],
            attrs={"eco_year": 1, "ets": 2, "cohort": 0},
        )

        engine.merge_parallel_veg(om)

        # Native merged DFM map order = [gid3, gid1, gid2]
        # Therefore expected vegetation = [0.30, 0.10, 0.20].
        merged_data = om.save_netcdf.call_args[0][2]

        assert np.allclose(
            merged_data["fraction"],
            [0.30, 0.10, 0.20],
        )

        # Ghost values 9.99 and 8.88 must not appear.
        assert 9.99 not in merged_data["fraction"]
        assert 8.88 not in merged_data["fraction"]
