"""Offline tests for the optional Source Discovery Tool adapter."""
from modules import spherex_sdt_adapter


def test_independent_catalogue_comparison_reports_agreement_and_disagreement():
    result = spherex_sdt_adapter.compare_with_photutils(
        [
            {"x_pixel": 10.0, "y_pixel": 20.0},
            {"x_pixel": 40.0, "y_pixel": 40.0},
        ],
        [
            {"xcentroid": 10.5, "ycentroid": 20.1},
            {"xcentroid": 70.0, "ycentroid": 70.0},
        ],
        match_radius_px=2.0,
    )
    assert result["matched_pairs"] == 1
    assert result["photutils_only_count"] == 1
    assert result["sdt_only_count"] == 1
    assert "not a calibrated" in result["interpretation"]

