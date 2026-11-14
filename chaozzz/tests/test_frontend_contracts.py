"""Session-state contracts consumed by the Streamlit integration pages."""
from __future__ import annotations

from modules import state


def test_state_has_contracts_for_optional_science_integrations():
    defaults = state._defaults()
    expected_keys = {
        "selected_epoch_ids",
        "multi_epoch_result",
        "temporal_pair_comparison",
        "sdt_comparison",
        "injection_recovery_report",
        "jpl_candidate_result",
        "archive_query_history",
    }
    assert expected_keys.issubset(defaults)
    assert defaults["selected_epoch_ids"] == []
    assert defaults["archive_query_history"] == []


def test_result_containers_start_empty_without_fake_measurements():
    defaults = state._defaults()
    for key in ("multi_epoch_result", "temporal_pair_comparison", "sdt_comparison", "injection_recovery_report", "jpl_candidate_result"):
        assert defaults[key] is None
