"""Catalogue cross-matching against IRSA-hosted survey tables."""
from __future__ import annotations

import csv
import io
import time
from typing import Any

import astropy.units as u
import streamlit as st
from astropy.coordinates import SkyCoord

import config
from modules.irsa_client import ArchiveError, _http_get_text

CATALOGUES = {
    "2MASS PSC": {
        "table": "fp_psc",
        "label": "2MASS Point Source Catalog (IRSA TAP fp_psc)",
        "columns": "designation, ra, dec, j_m, h_m, k_m",
        "id_column": "designation",
        "note": "J/H/K near-infrared point sources; useful for separating known stars from unmatched residuals.",
    },
    "AllWISE": {
        "table": "allwise_p3as_psd",
        "label": "AllWISE Source Catalog (IRSA TAP allwise_p3as_psd)",
        "columns": "designation, ra, dec, w1mpro, w2mpro",
        "id_column": "designation",
        "note": "Mid-infrared sources (W1/W2); a stronger overlap with SPHEREx wavelengths.",
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
    text, transport = _http_get_text(
        config.IRSA_TAP_SYNC,
        {"QUERY": query, "FORMAT": "csv", "MAXREC": int(max_rows)},
        timeout=config.ARCHIVE_TIMEOUT_SECONDS,
    )
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
    }
    return rows, diagnostics


def cross_match_candidates(
    candidates: list[Any],
    rows: list[dict],
    ra_column: str = "ra",
    dec_column: str = "dec",
    radius_arcsec: float = 3.0,
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
        matches = []
        for row, coordinate in usable:
            separation = float(coordinate.separation(target).arcsec)
            if separation <= radius_arcsec:
                matches.append({**row, "separation_arcsec": round(separation, 3)})
        matches.sort(key=lambda item: item["separation_arcsec"])
        if matches:
            results[candidate.candidate_id] = {
                "status": "matched",
                "matches": matches[:3],
                "detail": f"{len(matches)} catalogue source(s) within {radius_arcsec:.1f} arcsec.",
            }
        else:
            results[candidate.candidate_id] = {
                "status": "no match",
                "matches": [],
                "detail": f"No catalogue source within {radius_arcsec:.1f} arcsec; this does not establish a new object.",
            }
    return results
