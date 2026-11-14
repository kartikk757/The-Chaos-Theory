"""Offline tests for uncertainty-aware 2MASS and AllWISE matching."""
from types import SimpleNamespace

from modules import catalogue_matching


def _candidate(ra=10.0, dec=20.0):
    return SimpleNamespace(candidate_id="CAND-001", ra_deg=ra, dec_deg=dec)


def test_uncertainty_radius_and_multiple_counterparts_are_reported():
    rows = [
        {"designation": "near", "ra": "10.0001", "dec": "20.0", "err_maj": "1.0", "err_min": "1.0", "cc_flg": "000"},
        {"designation": "alt", "ra": "10.0005", "dec": "20.0", "err_maj": "1.0", "err_min": "1.0", "cc_flg": "000"},
    ]
    result = catalogue_matching.cross_match_candidates(
        [_candidate()], rows, radius_arcsec=1.0, candidate_sigma_arcsec=1.0, maximum_radius_arcsec=10.0
    )["CAND-001"]

    assert result["status"] == "matched"
    assert result["multiple_counterparts"] is True
    assert result["all_match_count"] == 2
    assert result["matches"][0]["separation_over_sigma"] < result["matches"][1]["separation_over_sigma"]
    assert result["matches"][0]["match_radius_arcsec"] > 1.0


def test_allwise_apparent_motion_is_diagnostic_not_solar_system_identity():
    row = {
        "designation": "WISE source",
        "ra": "10.0",
        "dec": "20.0",
        "sigra": "0.1",
        "sigdec": "0.1",
        "ra_pm": "10.0",
        "dec_pm": "20.0",
        "pmra": "500",
        "pmdec": "0",
        "sigpmra": "50",
        "sigpmdec": "50",
    }
    result = catalogue_matching.cross_match_candidates([_candidate()], [row], radius_arcsec=1.0)["CAND-001"]
    match = result["matches"][0]
    assert "AllWISE apparent-motion fit" in match["classification_hint"]
    assert "includes parallax" in match["allwise_motion_diagnostic"]["interpretation"]


def test_quality_flags_are_carried_with_2mass_match():
    row = {
        "designation": "confused",
        "ra": "10.0",
        "dec": "20.0",
        "err_maj": "1.0",
        "err_min": "0.5",
        "cc_flg": "c00",
        "prox": "4.0",
        "use_src": "0",
        "dup_src": "2",
    }
    match = catalogue_matching.cross_match_candidates([_candidate()], [row], radius_arcsec=1.0)["CAND-001"]["matches"][0]
    assert len(match["catalogue_quality"]) >= 3
