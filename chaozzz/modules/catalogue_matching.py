"""Catalogue cross-matching against IRSA-hosted survey tables."""
from __future__ import annotations

import csv
import io
import time
from typing import Any

import astropy.units as u
import numpy as np
import streamlit as st
from astropy.coordinates import SkyCoord

import config
from modules.irsa_client import ArchiveError, _http_get_text
from modules.network_guard import wait_for_service

CATALOGUES = {
    "2MASS PSC": {
        "table": "fp_psc",
        "label": "2MASS Point Source Catalog (IRSA TAP fp_psc)",
        "columns": "designation, ra, dec, err_maj, err_min, err_ang, date, jdate, j_m, h_m, k_m, cc_flg, prox, use_src, dup_src, mp_flg",
        "id_column": "designation",
        "note": "J/H/K point sources with 1-sigma position ellipse, observation date and quality/confusion fields. 2MASS positions are not proper-motion propagated.",
    },
    "AllWISE": {
        "table": "allwise_p3as_psd",
        "label": "AllWISE Source Catalog (IRSA TAP allwise_p3as_psd)",
        "columns": "designation, ra, dec, sigra, sigdec, sigradec, ra_pm, dec_pm, sigra_pm, sigdec_pm, sigradec_pm, pmra, sigpmra, pmdec, sigpmdec, w1mpro, w2mpro",
        "id_column": "designation",
        "note": "Mid-infrared sources with stationary and motion-fit positions. AllWISE apparent motion includes parallax and is not equivalent to proper motion.",
    },
}


@st.cache_data(ttl=config.CACHE_TTL_SECONDS, max_entries=16, show_spinner=False)
def fetch_catalogue(table_key: str, ra_deg: float, dec_deg: float, radius_deg: float, max_rows: int = 400) -> tuple[list[dict], dict]:
    spec = CATALOGUES.get(table_key)
    if spec is None:
        raise ArchiveError(f"Unknown catalogue selection: {table_key}")
    query = (
        f"SELECT TOP {int(max_rows)} {spec['columns']} FROM {spec['table']} "
        f"WHERE 1 = CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra_deg:.8f}, {dec_deg:.8f}, {radius_deg:.6f}))"
    )
    started = time.time()
    wait_for_service("irsa_tap", 0.5)
    try:
        text, transport = _http_get_text(
            config.IRSA_TAP_SYNC,
            {"QUERY": query, "FORMAT": "csv", "MAXREC": int(max_rows)},
            timeout=config.ARCHIVE_TIMEOUT_SECONDS,
        )
    except ArchiveError as exc:
        # Keep the existing catalogue selector functional on older IRSA table
        # versions while showing that astrometric enrichment was unavailable.
        basic_columns = "designation, ra, dec, j_m, h_m, k_m" if table_key == "2MASS PSC" else "designation, ra, dec, w1mpro, w2mpro"
        fallback_query = (
            f"SELECT TOP {int(max_rows)} {basic_columns} FROM {spec['table']} "
            f"WHERE 1 = CONTAINS(POINT('ICRS', ra, dec), CIRCLE('ICRS', {ra_deg:.8f}, {dec_deg:.8f}, {radius_deg:.6f}))"
        )
        wait_for_service("irsa_tap", 0.5)
        text, transport = _http_get_text(
            config.IRSA_TAP_SYNC,
            {"QUERY": fallback_query, "FORMAT": "csv", "MAXREC": int(max_rows)},
            timeout=config.ARCHIVE_TIMEOUT_SECONDS,
        )
        query = fallback_query
        enrichment_warning = f"Optional astrometric/quality columns were unavailable; used basic catalogue columns. Original query error: {exc}"
    else:
        enrichment_warning = None
    rows = list(csv.DictReader(io.StringIO(text)))
    diagnostics = {
        "catalogue": spec["label"],
        "table": spec["table"],
        "query": query,
        "rows": len(rows),
        "radius_deg": radius_deg,
        "elapsed_seconds": round(time.time() - started, 2),
        "endpoint": config.IRSA_TAP_SYNC,
        "note": spec["note"],
        "transport": transport,
        "astrometric_enrichment": "full" if enrichment_warning is None else "basic columns only",
        "warning": enrichment_warning,
        "epoch_note": "2MASS is matched at its catalogued position without propagation; AllWISE apparent-motion fields are diagnostic only because the measured motion includes parallax and is not guaranteed to extrapolate linearly.",
    }
    return rows, diagnostics


def cross_match_candidates(
    candidates: list[Any],
    rows: list[dict],
    ra_column: str = "ra",
    dec_column: str = "dec",
    radius_arcsec: float = 3.0,
    candidate_sigma_arcsec: float = 3.1,
    sigma_multiplier: float = 3.0,
    maximum_radius_arcsec: float = 15.0,
) -> dict[str, dict]:
    results: dict[str, dict] = {}
    usable = []
    for row in rows:
        try:
            usable.append((row, SkyCoord(float(row[ra_column]) * u.deg, float(row[dec_column]) * u.deg, frame="icrs")))
        except (KeyError, TypeError, ValueError):
            continue
    for candidate in candidates:
        if candidate.ra_deg is None or candidate.dec_deg is None:
            results[candidate.candidate_id] = {"status": "unavailable", "matches": [], "detail": "Candidate sky position is unavailable."}
            continue
        try:
            target = SkyCoord(float(candidate.ra_deg) * u.deg, float(candidate.dec_deg) * u.deg, frame="icrs")
        except Exception:
            results[candidate.candidate_id] = {"status": "unavailable", "matches": [], "detail": "Candidate coordinate could not be parsed."}
            continue
        candidate_sigma = max(float(candidate_sigma_arcsec), 0.1)
        matches = []
        for row, coordinate in usable:
            separation = float(coordinate.separation(target).arcsec)
            cat_sigma = _catalogue_position_sigma(row)
            match_radius = min(
                float(maximum_radius_arcsec),
                max(float(radius_arcsec), float(sigma_multiplier) * np.hypot(candidate_sigma, cat_sigma)),
            )
            if separation <= match_radius:
                motion = _allwise_motion_diagnostic(row) if "ra_pm" in row else None
                quality = _catalogue_quality(row)
                matches.append({
                    **row,
                    "separation_arcsec": round(separation, 3),
                    "candidate_sigma_arcsec": candidate_sigma,
                    "catalogue_sigma_arcsec": cat_sigma,
                    "combined_sigma_arcsec": float(np.hypot(candidate_sigma, cat_sigma)),
                    "match_radius_arcsec": match_radius,
                    "separation_over_sigma": separation / max(float(np.hypot(candidate_sigma, cat_sigma)), 1e-6),
                    "classification_hint": _classification_hint(row, motion),
                    "allwise_motion_diagnostic": motion,
                    "catalogue_quality": quality,
                    "epoch_handling": "No catalogue-position propagation was applied; see epoch note in query provenance.",
                })
        matches.sort(key=lambda item: (item["separation_over_sigma"], item["separation_arcsec"]))
        if matches:
            results[candidate.candidate_id] = {
                "status": "matched",
                "matches": matches[:5],
                "all_match_count": len(matches),
                "detail": f"{len(matches)} catalogue source(s) within a combined-uncertainty radius (configured floor {radius_arcsec:.1f} arcsec; cap {maximum_radius_arcsec:.1f} arcsec).",
                "multiple_counterparts": len(matches) > 1,
                "candidate_sigma_arcsec": candidate_sigma,
                "quality_explanation": "Match ordering uses separation divided by the combined positional uncertainty; review all alternatives and source-quality flags.",
            }
        else:
            results[candidate.candidate_id] = {
                "status": "no match",
                "matches": [],
                "detail": f"No catalogue source within the uncertainty-scaled search radius (floor {radius_arcsec:.1f} arcsec); this does not establish a new object.",
                "candidate_sigma_arcsec": candidate_sigma,
                "multiple_counterparts": False,
            }
    return results


def _number(row: dict, name: str) -> float | None:
    try:
        value = float(row.get(name))
        return value if np.isfinite(value) else None
    except (TypeError, ValueError):
        return None


def _catalogue_position_sigma(row: dict) -> float:
    """Return a conservative radial 1-sigma from catalogue position columns."""
    major = _number(row, "err_maj")
    minor = _number(row, "err_min")
    if major is not None and minor is not None and major > 0 and minor > 0:
        return float(np.hypot(major, minor) / np.sqrt(2.0))
    sigra = _number(row, "sigra")
    sigdec = _number(row, "sigdec")
    if sigra is not None and sigdec is not None and sigra > 0 and sigdec > 0:
        # AllWISE RA uncertainty is in angular units for the coordinate; use
        # the documented coordinate ellipse as supplied, without dec scaling.
        return float(np.hypot(sigra, sigdec) / np.sqrt(2.0))
    return 1.0


def _allwise_motion_diagnostic(row: dict) -> dict | None:
    ra, dec = _number(row, "ra_pm"), _number(row, "dec_pm")
    pmra, pmdec = _number(row, "pmra"), _number(row, "pmdec")
    err_ra, err_dec = _number(row, "sigpmra"), _number(row, "sigpmdec")
    if None in (ra, dec, pmra, pmdec, err_ra, err_dec) or err_ra < 0 or err_dec < 0:
        return None
    motion_snr = float(np.hypot(pmra, pmdec) / max(np.hypot(err_ra, err_dec), 1e-9))
    return {
        "position_ra_deg_at_mjd_55400": ra,
        "position_dec_deg_at_mjd_55400": dec,
        "reference_mjd_utc": 55400.0,
        "apparent_motion_ra_mas_per_year": pmra,
        "apparent_motion_dec_mas_per_year": pmdec,
        "apparent_motion_uncertainty_ra_mas_per_year": err_ra,
        "apparent_motion_uncertainty_dec_mas_per_year": err_dec,
        "motion_signal_to_uncertainty": motion_snr,
        "interpretation": "AllWISE apparent motion includes parallax and can be affected by source confusion, artifacts and noise; this is not a proper-motion or solar-system-object identification.",
    }


def _catalogue_quality(row: dict) -> list[str]:
    flags = []
    cc_flg = str(row.get("cc_flg", "")).strip()
    if cc_flg and set(cc_flg) - {"0", "", "n", "N"}:
        flags.append(f"2MASS contamination/confusion flag: {cc_flg}")
    prox = _number(row, "prox")
    if prox is not None and prox < 6.0:
        flags.append(f"2MASS nearest-neighbour distance {prox:.1f} arcsec indicates possible confusion")
    if str(row.get("use_src", "1")).strip() == "0":
        flags.append("2MASS source is excluded by use_src quality selection")
    dup = _number(row, "dup_src")
    if dup is not None and dup > 0:
        flags.append(f"2MASS duplicate/overlap flag {int(dup)}")
    return flags


def _classification_hint(row: dict, motion: dict | None) -> str:
    if motion and motion["motion_signal_to_uncertainty"] >= 3.0:
        return "AllWISE apparent-motion fit reported (investigate; includes parallax)"
    if str(row.get("mp_flg", "0")).strip() not in ("", "0", "--", "null"):
        return "2MASS solar-system/asteroid association flag present; verify the linked association"
    return "Catalogue counterpart at catalogued position; stationary interpretation is provisional"
