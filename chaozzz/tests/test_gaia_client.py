"""Offline tests for Gaia DR3 epoch propagation and uncertainty-aware matching."""
from __future__ import annotations

import math

from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.time import Time
import pytest

from modules import gaia_client


def _row(**overrides):
    row = {
        "source_id": "1234567890123456789",
        "ra": "359.999",
        "dec": "20.0",
        "ref_epoch": "2016.0",
        "ra_error": "10.0",
        "dec_error": "20.0",
        "ra_dec_corr": "0.0",
        "pmra": "1000.0",
        "pmra_error": "10.0",
        "pmdec": "0.0",
        "pmdec_error": "20.0",
        "ra_pmra_corr": "0.0",
        "ra_pmdec_corr": "0.0",
        "dec_pmra_corr": "0.0",
        "dec_pmdec_corr": "0.0",
        "pmra_pmdec_corr": "0.0",
        "parallax": "2.0",
        "parallax_error": "0.1",
        "astrometric_params_solved": "31",
        "ruwe": "1.0",
        "phot_g_mean_mag": "15.0",
        "phot_bp_mean_mag": "15.5",
        "phot_rp_mean_mag": "14.4",
    }
    row.update(overrides)
    return row


def test_gaia_source_position_is_propagated_with_pmra_star():
    target_time = Time(2020.0, format="jyear", scale="tcb")
    position = gaia_client.propagate_gaia_source(_row(), target_time.utc.mjd)
    expected = SkyCoord(359.999 * u.deg, 20.0 * u.deg).directional_offset_by(90 * u.deg, 4 * u.arcsec)

    assert position["status"] == "ok"
    assert SkyCoord(position["ra_deg"] * u.deg, position["dec_deg"] * u.deg).separation(expected).arcsec < 0.01
    assert position["delta_years"] == pytest.approx(4.0, abs=1e-6)
    assert position["proper_motion_available"] is True
    assert position["sigma_arcsec"] == pytest.approx(math.sqrt(8500) / 1000, rel=1e-6)


def test_gaia_source_without_proper_motion_is_skipped_for_long_baseline():
    target_time = Time(2026.0, format="jyear", scale="tcb")
    result = gaia_client.propagate_gaia_source(_row(pmra="", pmdec=""), target_time.utc.mjd)
    assert result["status"] == "unavailable"
    assert "proper-motion" in result["reason"]


def test_uncertainty_aware_crossmatch_reports_epoch_and_radius():
    target_time = Time(2020.0, format="jyear", scale="tcb")
    propagated = gaia_client.propagate_gaia_source(_row(), target_time.utc.mjd)
    candidate = {
        "candidate_id": "CAND-001",
        "ra_deg": propagated["ra_deg"],
        "dec_deg": propagated["dec_deg"],
        "mjd": target_time.utc.mjd,
        "position_sigma_arcsec": 1.0,
    }
    matches = gaia_client.cross_match_candidates([candidate], [_row()], maximum_radius_arcsec=10.0)

    match = matches["CAND-001"]
    assert match["status"] == "matched"
    assert match["matches"][0]["source_id"] == "1234567890123456789"
    assert match["matches"][0]["separation_arcsec"] < 0.01
    assert match["matches"][0]["match_radius_arcsec"] >= 3.0
    assert match["matches"][0]["propagated_to_jyear_tcb"] == pytest.approx(2020.0, abs=1e-6)


def test_gaia_tap_query_is_read_only_and_records_provenance(monkeypatch):
    class FakeResponse:
        text = "source_id,ra,dec\n1,10.0,20.0\n"

        def raise_for_status(self):
            return None

    captured = {}

    def fake_post(url, *, data, timeout, headers):
        captured.update(url=url, data=data, timeout=timeout)
        return FakeResponse()

    monkeypatch.setattr(gaia_client.requests, "post", fake_post)
    rows, diagnostics = gaia_client.fetch_gaia_dr3(10.0, 20.0, 0.05, 20)

    assert rows[0]["source_id"] == "1"
    assert "SELECT TOP 20" in diagnostics["query"]
    assert diagnostics["endpoint"] == gaia_client.GAIA_TAP_SYNC
    assert captured["data"]["REQUEST"] == "doQuery"
    assert captured["data"]["LANG"] == "ADQL"

