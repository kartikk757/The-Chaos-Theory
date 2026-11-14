"""Offline tests for known-object ephemeris comparison."""
from astropy.coordinates import SkyCoord
import astropy.units as u

from modules import jpl_crossmatch


def test_horizons_csv_parser_reads_degree_columns():
    payload = {"result": "header\nDate__(UT), , R.A.__(ICRF), DEC__(ICRF)\n$$SOE\n2026-Jan-01, , 10.0, 20.0\n$$EOE\n"}
    rows = jpl_crossmatch._parse_horizons_csv(payload, 1)
    assert rows[0]["ra_deg"] == 10.0
    assert rows[0]["dec_deg"] == 20.0


def test_known_object_investigation_compares_exact_epochs_without_network(monkeypatch):
    candidate = {
        "candidate_id": "TRACK-1",
        "detections": [
            {"observation_id": "E1", "mjd": 60000.0, "ra_deg": 10.0, "dec_deg": 20.0},
            {"observation_id": "E2", "mjd": 60001.0, "ra_deg": 10.001, "dec_deg": 20.0},
        ],
    }
    sbdb = {"status": "ok", "object": {"name": "Known object"}, "orbit": {"elements": []}}
    monkeypatch.setattr(jpl_crossmatch, "fetch_sbdb", lambda target: sbdb)

    def fake_horizons(target, jd_values, center):
        return {
            "positions": [
                {"requested_jd_utc": jd, "ra_deg": ra, "dec_deg": 20.0}
                for jd, ra in zip(jd_values, (10.0, 10.001))
            ],
            "endpoint": "offline fake",
            "center": center,
            "timescale": "UT",
            "reference_frame": "ICRF",
        }

    monkeypatch.setattr(jpl_crossmatch, "_fetch_horizons_positions", fake_horizons)
    result = jpl_crossmatch.investigate_known_object("known designation", [candidate], match_radius_arcsec=1.0)

    assert result["status"] == "ok"
    assert result["summary"]["epochs_within_radius"] == 2
    assert result["summary"]["hypothesis_consistent_with_positions"] is True
    assert all(row["angular_separation_arcsec"] < 1e-6 for row in result["matches"])
    assert result["sbdb"] == sbdb


def test_ephemeris_mismatch_is_not_treated_as_unknown_object(monkeypatch):
    candidate = {"candidate_id": "C1", "mjd": 60000.0, "ra_deg": 10.0, "dec_deg": 20.0}
    monkeypatch.setattr(jpl_crossmatch, "fetch_sbdb", lambda target: {"status": "ok", "object": {}, "orbit": {}})
    monkeypatch.setattr(jpl_crossmatch, "_fetch_horizons_positions", lambda target, jds, center: {
        "positions": [{"requested_jd_utc": jds[0], "ra_deg": 11.0, "dec_deg": 20.0}],
        "endpoint": "offline fake", "center": center, "timescale": "UT", "reference_frame": "ICRF",
    })
    result = jpl_crossmatch.investigate_known_object("known", [candidate], match_radius_arcsec=1.0)
    assert result["summary"]["hypothesis_consistent_with_positions"] is False
    assert "does not imply" in result["limitations"][-1]

