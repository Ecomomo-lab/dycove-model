from types import SimpleNamespace
from unittest.mock import MagicMock

import numpy as np
import pytest
import xarray as xr

from dycove.sim.outputs import OutputManager


def make_output_manager(tmp_path, storage):
    engine = MagicMock()
    engine.veg = MagicMock()
    engine.organic = SimpleNamespace(
        sediment_storage=storage,
        target_layer=0,
        target_sediment=2,
    )
    engine.model_dir = tmp_path
    engine.is_parallel.return_value = False
    engine.get_rank.return_value = 0
    return OutputManager(engine)


def cell_records():
    return [{
        "species": "species_1",
        "fveg": np.array([0.2, 0.4]),
        "agb_estimate": np.array([0.01, 0.02]),
        "delta_msed": np.array([0.001, 0.002]),
    }]


@pytest.mark.unit
def test_organic_netcdf_uses_neutral_name_and_records_units(tmp_path):
    outputs = make_output_manager(tmp_path, "msed")
    outputs.save_organic_step(
        SimpleNamespace(eco_year=1, ets=2), records=[],
        cell_records=cell_records(),
    )

    path = tmp_path / "om_output" / "om_year0001_ets002.nc"
    with xr.open_dataset(path, engine="scipy") as ds:
        assert "organic_increment" in ds
        assert "delta_msed" not in ds
        np.testing.assert_allclose(ds["organic_increment"], [[0.001, 0.002]])
        assert ds["fveg"].attrs["units"] == "1"
        assert ds["agb_estimate"].attrs["units"] == "kg m-2"
        assert ds["organic_increment"].attrs["units"] == "kg m-2"
        assert ds["cell"].attrs["units"] == "1"
        assert ds.attrs["storage"] == "msed"
        assert ds.attrs["target_layer"] == 0
        assert ds.attrs["target_sediment"] == 2


@pytest.mark.unit
def test_organic_netcdf_marks_layer_not_applicable_for_bodsed(tmp_path):
    outputs = make_output_manager(tmp_path, "bodsed")
    outputs.save_organic_step(
        SimpleNamespace(eco_year=1, ets=2), records=[],
        cell_records=cell_records(),
    )

    path = tmp_path / "om_output" / "om_year0001_ets002.nc"
    with xr.open_dataset(path, engine="scipy") as ds:
        assert ds.attrs["storage"] == "bodsed"
        assert ds.attrs["target_layer"] == "not_applicable"
