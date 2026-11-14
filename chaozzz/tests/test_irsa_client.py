"""Offline tests for coordinate parsing, query validation and status plumbing.

No network access is required: these tests cover the deterministic parts of the
archive client. Live retrieval is exercised separately by the smoke script.
"""
from __future__ import annotations

import pytest

import config
from modules import irsa_client


def test_parse_decimal_and_sexagesimal_ra():
    assert irsa_client.parse_coordinate("210.80227", "ra") == pytest.approx(210.80227)
    assert irsa_client.parse_coordinate("14:03:12.65", "ra") == pytest.approx(210.8027, abs=1e-3)
    assert irsa_client.parse_coordinate("359.999", "ra") == pytest.approx(359.999)


def test_parse_declination_and_wrap():
    assert irsa_client.parse_coordinate("+54:20:56.2", "dec") == pytest.approx(54.34895, abs=1e-4)
    assert irsa_client.parse_coordinate("-12.5", "dec") == pytest.approx(-12.5)
    with pytest.raises(ValueError):
        irsa_client.parse_coordinate("120.0", "dec")


def test_parse_rejects_empty_and_garbage():
    with pytest.raises(ValueError):
        irsa_client.parse_coordinate("", "ra")
    with pytest.raises(ValueError):
        irsa_client.parse_coordinate("not-a-coordinate", "dec")


def test_mjd_round_trip():
    mjd = irsa_client.mjd_from_date("2026-05-17")
    assert irsa_client.mjd_to_iso(mjd).startswith("2026-05-17")
    assert irsa_client.mjd_to_iso("") == ""


def test_validate_query_accepts_documented_default():
    validated = irsa_client.validate_query(**config.DEFAULT_QUERY)
    assert validated["ra_deg"] == pytest.approx(config.DEFAULT_QUERY["ra_deg"])
    assert validated["end_mjd"] > validated["start_mjd"]
    assert validated["band"] == "SPHEREx-D2"


@pytest.mark.parametrize(
    "overrides, message",
    [
        ({"radius_deg": 5.0}, "Search radius"),
        ({"radius_deg": 0.0001}, "Search radius"),
        ({"size_deg": 2.0}, "Cutout width"),
        ({"band": "SPHEREx-D9"}, "detector band"),
        ({"release": "qr2"}, "release"),
        ({"survey": "medium"}, "survey"),
        ({"max_rows": 5000}, "Maximum observation rows"),
        ({"end_date": "2026-01-01", "start_date": "2026-06-01"}, "End date"),
        ({"ra_deg": "abc"}, "RA"),
        ({"dec_deg": "abc"}, "DEC"),
    ],
)
def test_validate_query_rejects_invalid_input(overrides, message):
    params = dict(config.DEFAULT_QUERY)
    params.update(overrides)
    with pytest.raises(ValueError) as error:
        irsa_client.validate_query(**params)
    assert message.lower() in str(error.value).lower()


def test_build_observation_record_uses_documented_fields():
    row = {
        "uri": "https://irsa.ipac.caltech.edu/ibe/data/spherex/qr3/level2/2026W20_2A/l2b-v27-2026-267/2/xD2_spx.fits?center=210.8,54.3d&size=0.1",
        "raw_uri": "ibe/data/spherex/qr3/level2/2026W20_2A/l2b-v27-2026-267/2/level2_2026W20_2A_0464_1D2_spx_l2b-v27-2026-267.fits",
        "time_bounds_lower": "61177.453993",
        "energy_bandpassname": "SPHEREx-D2",
        "energy_bounds_lower": "1.10e-6",
        "energy_bounds_upper": "1.62e-6",
        "calibrationlevel": "2",
        "quality_flag": "",
        "provenance_version": "7.0.5",
        "productid": "D2",
        "obsid": "1234",
        "planeid": "5678",
        "contentlength": "4000320",
        "producttype": "science",
        "contenttype": "image/fits",
        "obs_publisher_did": "ivo://nasa.heasarc/spherex/2026W20_2A_0464_1",
    }
    record = irsa_client.build_observation_record(row)
    assert record["band"] == "SPHEREx-D2"
    assert record["release"] == "qr3"
    assert record["is_real_spherex"] is True
    assert record["cutout_url"].endswith("size=0.1")
    assert record["full_fits_url"].startswith("https://irsa.ipac.caltech.edu/ibe/data/spherex/qr3/")
    assert "?center=" not in record["full_fits_url"]
    assert record["date_utc"].startswith("2026-05-17")
    assert record["wavelength_range_um"] == (1.1, 1.62)
    assert record["content_length"] == 4000320
    assert any("quality_flag" in warning for warning in record["warnings"])


def test_sia2_enrichment_marks_unmatched_rows():
    observation = irsa_client.build_observation_record({"raw_uri": "ibe/data/spherex/qr3/a/b/fileD2.fits", "time_bounds_lower": "61000.0", "energy_bandpassname": "SPHEREx-D2"})
    before = list(observation["warnings"])
    enriched = irsa_client.enrich_observations_with_sia2([observation], [], 210.8, 54.3)
    assert enriched[0]["ra_deg"] == 210.8
    assert "requested field centre" in enriched[0]["coordinates_source"]
    assert enriched[0]["coordinates_role"] == "requested_field_center"
    assert enriched[0]["metadata"]["requested_field_center"] == {
        "ra_deg": 210.8,
        "dec_deg": 54.3,
        "coordinate_frame": "ICRS",
        "source": "IRSA archive query / cutout request",
    }
    assert len(enriched[0]["warnings"]) == len(before) + 1


def test_sia2_collection_alias_mapping():
    assert irsa_client._sia2_collections("all", "wide") == ["spherex_dr1", "spherex_qr3"]
    assert irsa_client._sia2_collections("qr3", "deep") == ["spherex_qr3_deep"]
