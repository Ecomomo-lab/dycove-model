import pytest
from unittest.mock import MagicMock, patch

from dycove.utils.simulation_reporting import Reporter


@pytest.mark.unit
def test_reporter_is_singleton():
    r1 = Reporter()
    r2 = Reporter()
    assert r1 is r2


@pytest.mark.unit
def test_print_model_time_info_veg_active_reports_all_three_messages():
    r = Reporter()
    simstate = MagicMock()
    simstate.n_ets = 14
    simstate.veg_interval = 43200
    simstate.hydro_sim_days = 365
    simstate.veg_sim_years = 10

    with patch.object(r, "report") as mock_report:
        r.print_model_time_info(simstate, veg_active=True)
        assert mock_report.call_count == 3


@pytest.mark.unit
def test_print_model_time_info_veg_inactive_reports_only_hydro():
    r = Reporter()
    simstate = MagicMock()
    simstate.hydro_sim_days = 365

    with patch.object(r, "report") as mock_report:
        r.print_model_time_info(simstate, veg_active=False)
        assert mock_report.call_count == 1