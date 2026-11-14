"""Known-object investigation through JPL SBDB and Horizons.

This module only investigates a user-supplied known-object identifier. It does
not perform blind solar-system object discovery or infer that an unmatched
candidate is unknown.
"""
from __future__ import annotations

import csv
import re
import time
from typing import Any

import astropy.units as u
import numpy as np
import requests
import streamlit as st
from astropy.coordinates import SkyCoord

import config
from modules.network_guard import wait_for_service


class JPLQueryError(RuntimeError):
    """Raised when JPL's public small-body services fail or return no data."""


def _get_json(url: str, params: dict[str, Any], service: str) -> dict:
    last_error: Exception | None = None
    for attempt in range(2):
        wait_for_service(f"jpl_{service}", 1.0)
        try:
            response = requests.get(url, params=params, timeout=(10, 45), headers={"User-Agent": "TheChaosTheory/1.1 (public astronomy data client)"})
            response.raise_for_status()
            payload = response.json()
            if payload.get("error"):
                raise JPLQueryError(str(payload["error"]))
            return payload
        except (requests.RequestException, ValueError, JPLQueryError) as exc:
            last_error = exc
            if isinstance(exc, JPLQueryError) or attempt == 1:
                break
            time.sleep(0.4)
    raise JPLQueryError(f"JPL {service} request failed: {last_error}") from last_error


@st.cache_data(ttl=3600, max_entries=32, show_spinner=False)
def fetch_sbdb(target: str) -> dict:
    """Resolve a known small-body identifier and return SBDB object/orbit data."""
    query = str(target).strip()
    if not query or len(query) > 100:
        raise JPLQueryError("Enter a known small-body name, designation, or SPK identifier (1–100 characters).")
    payload = _get_json(
        config.JPL_SBDB_API,
        {"sstr": query, "alt-des": "true", "alt-spk": "true", "cd-epoch": "true", "full-prec": "true"},
        "sbdb",
    )
    obj = payload.get("object")
    orbit = payload.get("orbit")
    if not obj:
        raise JPLQueryError("SBDB did not resolve this identifier to a known object.")
    return {
        "status": "ok",
        "query": query,
        "object": obj,
        "orbit": orbit,
        "signature": payload.get("signature"),
        "endpoint": config.JPL_SBDB_API,
        "provenance": "JPL SBDB API response; object and orbit identifiers are retained verbatim.",
    }


def _candidate_positions(candidate: Any) -> list[dict]:
    if isinstance(candidate, dict):
        detections = candidate.get("detections")
        if isinstance(detections, list):
            return [dict(row) for row in detections if isinstance(row, dict)]
        return [candidate]
    metadata = getattr(candidate, "metadata", {}) or {}
    output = []
    for name in ("epoch1", "epoch2"):
        row = metadata.get(name, {}) if isinstance(metadata, dict) else {}
        if row.get("mjd") is not None:
            output.append({"ra_deg": row.get("ra_deg", getattr(candidate, "ra_deg", None)), "dec_deg": row.get("dec_deg", getattr(candidate, "dec_deg", None)), "mjd": row.get("mjd"), "observation_id": row.get("id")})
    if not output:
        output.append({"ra_deg": getattr(candidate, "ra_deg", None), "dec_deg": getattr(candidate, "dec_deg", None), "mjd": getattr(candidate, "mjd", None)})
    return output


def _parse_horizons_csv(payload: dict, requested_count: int) -> list[dict]:
    text = str(payload.get("result", ""))
    if "$$SOE" not in text or "$$EOE" not in text:
        message = re.search(r"(?:INPUT ERROR|Cannot|No ephemeris).*", text, flags=re.IGNORECASE)
        raise JPLQueryError(f"Horizons did not return an ephemeris table: {message.group(0) if message else text[:500]}")
    before, remainder = text.split("$$SOE", 1)
    body, _ = remainder.split("$$EOE", 1)
    header = None
    for line in reversed(before.splitlines()):
        cells = next(csv.reader([line]))
        lowered = [cell.strip().lower() for cell in cells]
        if any("r.a." in cell or cell == "ra" for cell in lowered) and any("dec" in cell for cell in lowered):
            header = [cell.strip() for cell in cells]
            break
    if not header:
        raise JPLQueryError("Horizons returned a table without recognizable RA/Dec column names.")
    ra_index = next(i for i, cell in enumerate(header) if "r.a." in cell.lower() or cell.lower() == "ra")
    dec_index = next(i for i, cell in enumerate(header) if "dec" in cell.lower())
    rows = []
    for line in body.splitlines():
        if not line.strip():
            continue
        cells = next(csv.reader([line]))
        if max(ra_index, dec_index) >= len(cells):
            continue
        try:
            ra, dec = float(cells[ra_index].strip()), float(cells[dec_index].strip())
        except ValueError:
            continue
        rows.append({"ra_deg": ra, "dec_deg": dec, "raw_csv_row": cells})
    if len(rows) != requested_count:
        raise JPLQueryError(f"Horizons returned {len(rows)} of {requested_count} requested ephemeris epoch(s).")
    return rows


@st.cache_data(ttl=1800, max_entries=32, show_spinner=False)
def _fetch_horizons_positions(target: str, jd_values: tuple[float, ...], center: str = "500@399") -> dict:
    if not jd_values or len(jd_values) > 50:
        raise JPLQueryError("Horizons lookup requires 1–50 observation times.")
    # Quotes are required by Horizons' TLIST syntax; requests performs URL encoding.
    tlist = " ".join(f"'{float(jd):.9f}'" for jd in jd_values)
    payload = _get_json(
        config.JPL_HORIZONS_API,
        {
            "format": "json", "COMMAND": f"'{str(target).strip()}'", "OBJ_DATA": "YES",
            "MAKE_EPHEM": "YES", "EPHEM_TYPE": "OBSERVER", "CENTER": center,
            "TLIST": tlist, "TLIST_TYPE": "JD", "TIME_TYPE": "UT", "TIME_DIGITS": "SECONDS",
            "QUANTITIES": "1", "ANG_FORMAT": "DEG", "CSV_FORMAT": "YES", "EXTRA_PREC": "YES",
        },
        "horizons",
    )
    positions = _parse_horizons_csv(payload, len(jd_values))
    for position, jd in zip(positions, jd_values):
        position["requested_jd_utc"] = float(jd)
        position["requested_mjd_utc"] = float(jd - 2400000.5)
    return {
        "positions": positions, "endpoint": config.JPL_HORIZONS_API, "signature": payload.get("signature"),
        "center": center, "timescale": "UT observer ephemeris; input times are UTC-derived JD",
        "reference_frame": "ICRF, as returned by Horizons",
        "raw_result_header": str(payload.get("result", "")).split("$$SOE", 1)[0][-1600:],
    }


def investigate_known_object(target: str, candidates: list[Any], match_radius_arcsec: float = 60.0, center: str = "500@399") -> dict:
    """Resolve a user-entered known object and compare its ephemeris to positions."""
    if not 0.1 <= float(match_radius_arcsec) <= 3600.0:
        raise JPLQueryError("Match radius must be between 0.1 and 3600 arcsec.")
    sbdb = fetch_sbdb(target)
    observations = []
    for candidate_index, candidate in enumerate(candidates):
        for position in _candidate_positions(candidate):
            try:
                ra, dec, mjd = float(position["ra_deg"]), float(position["dec_deg"]), float(position["mjd"])
                if not (0.0 <= ra < 360.0 and -90.0 <= dec <= 90.0 and np.isfinite(mjd)):
                    continue
            except (KeyError, TypeError, ValueError):
                continue
            observations.append({**position, "candidate_index": candidate_index, "ra_deg": ra, "dec_deg": dec, "mjd": mjd})
    if not observations:
        return {"status": "insufficient_candidate_epochs", "sbdb": sbdb, "matches": [], "detail": "The known object was resolved, but candidates did not include both sky coordinates and UTC MJD observation times."}
    # De-duplicate requested epochs while preserving the exact response mapping.
    unique_times = sorted({round(row["mjd"] + 2400000.5, 9) for row in observations})
    ephemeris = _fetch_horizons_positions(target, tuple(unique_times), center)
    predicted_by_jd = {row["requested_jd_utc"]: row for row in ephemeris["positions"]}
    matches = []
    for row in observations:
        jd = round(row["mjd"] + 2400000.5, 9)
        predicted = predicted_by_jd[jd]
        candidate_coord = SkyCoord(row["ra_deg"] * u.deg, row["dec_deg"] * u.deg, frame="icrs")
        predicted_coord = SkyCoord(predicted["ra_deg"] * u.deg, predicted["dec_deg"] * u.deg, frame="icrs")
        separation = float(candidate_coord.separation(predicted_coord).arcsec)
        matches.append({
            "candidate_index": row["candidate_index"], "observation_id": row.get("observation_id"), "mjd_utc": row["mjd"],
            "candidate_ra_deg": row["ra_deg"], "candidate_dec_deg": row["dec_deg"],
            "predicted_ra_deg": predicted["ra_deg"], "predicted_dec_deg": predicted["dec_deg"],
            "angular_separation_arcsec": separation, "match_radius_arcsec": float(match_radius_arcsec),
            "within_match_radius": separation <= float(match_radius_arcsec),
        })
    near = [row for row in matches if row["within_match_radius"]]
    return {
        "status": "ok", "query_type": "user-specified known-object hypothesis only", "sbdb": sbdb,
        "ephemeris": ephemeris, "matches": matches,
        "summary": {"epochs_compared": len(matches), "epochs_within_radius": len(near),
                    "minimum_separation_arcsec": min(row["angular_separation_arcsec"] for row in matches),
                    "hypothesis_consistent_with_positions": bool(near),
                    "criteria": f"At least one candidate epoch lies within {float(match_radius_arcsec):g} arcsec of the Horizons ephemeris."},
        "limitations": [
            "A positional match investigates a known-object hypothesis; it is not a new-object discovery or orbital determination.",
            "The ephemeris uses the specified observer centre and Horizons time/reference settings; it does not model SPHEREx-specific observing geometry unless supplied.",
            "Candidate coordinate uncertainties are not folded into the user-selected angular match radius.",
            "No match is inconclusive and does not imply that a candidate is unknown.",
        ],
    }
