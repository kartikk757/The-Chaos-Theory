"""Gaia DR3 cone search and uncertainty-aware epoch cross-match."""
from __future__ import annotations

import csv
import io
import math
import time
from typing import Any

import astropy.units as u
import numpy as np
import requests
import streamlit as st
from astropy.coordinates import SkyCoord
from astropy.time import Time

import config
from modules.network_guard import wait_for_service

GAIA_TAP_SYNC = config.GAIA_TAP_SYNC
GAIA_TABLE = "gaiadr3.gaia_source"
GAIA_COLUMNS = (
    "source_id, ra, dec, ref_epoch, ra_error, dec_error, ra_dec_corr, "
    "pmra, pmra_error, pmdec, pmdec_error, ra_pmra_corr, ra_pmdec_corr, "
    "dec_pmra_corr, dec_pmdec_corr, pmra_pmdec_corr, parallax, parallax_error, "
    "astrometric_params_solved, ruwe, phot_g_mean_mag, phot_bp_mean_mag, phot_rp_mean_mag"
)


class GaiaQueryError(RuntimeError):
    """Raised when the public Gaia TAP query cannot return a usable response."""


@st.cache_data(ttl=1800, max_entries=12, show_spinner=False)
def fetch_gaia_dr3(ra_deg: float, dec_deg: float, radius_deg: float, max_rows: int = 10000) -> tuple[list[dict], dict]:
    """Query the ESA Gaia DR3 source table within an ICRS cone."""
    if not (0.0 <= float(ra_deg) < 360.0 and -90.0 <= float(dec_deg) <= 90.0):
        raise GaiaQueryError("The Gaia cone-search centre is outside valid ICRS coordinates.")
    if not (0.0 < float(radius_deg) <= 2.0):
        raise GaiaQueryError("The Gaia cone-search radius must be greater than 0 and no more than 2 degrees.")
    row_limit = max(1, min(int(max_rows), 10000))
    query = (
        f"SELECT TOP {row_limit} {GAIA_COLUMNS} FROM {GAIA_TABLE} "
        f"WHERE 1 = CONTAINS(POINT('ICRS', ra, dec), "
        f"CIRCLE('ICRS', {float(ra_deg):.9f}, {float(dec_deg):.9f}, {float(radius_deg):.7f}))"
    )
    started = time.monotonic()
    response = None
    last_error = None
    for attempt in range(2):
        wait_for_service("gaia_tap", 1.0)
        try:
            response = requests.post(
                GAIA_TAP_SYNC,
                data={"REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "csv", "MAXREC": row_limit, "QUERY": query},
                timeout=(10, 60),
                headers={"User-Agent": "TheChaosTheory/1.0 (Gaia DR3 cross-match)"},
            )
            response.raise_for_status()
            break
        except requests.RequestException as exc:
            last_error = exc
            if attempt == 1 or (getattr(exc.response, "status_code", 0) or 0) < 500:
                raise GaiaQueryError(f"Gaia DR3 TAP request failed: {exc}") from exc
            time.sleep(0.5)
    if response is None or last_error is not None and not response.ok:
        raise GaiaQueryError(f"Gaia DR3 TAP request failed: {last_error}") from last_error
    body = response.text.strip()
    if not body or body.lower().startswith(("error", "<!doctype html", "<html")):
        raise GaiaQueryError(f"Gaia DR3 TAP returned an error response: {body[:400] or 'empty response'}")
    try:
        rows = list(csv.DictReader(io.StringIO(body)))
        if rows and not any(key in rows[0] for key in ("source_id", "SOURCE_ID")):
            raise ValueError("source_id column is missing")
    except Exception as exc:
        raise GaiaQueryError(f"Could not parse Gaia DR3 CSV response: {exc}") from exc
    diagnostics = {
        "status": "ok",
        "catalogue": "Gaia DR3 gaia_source",
        "release": "DR3",
        "endpoint": GAIA_TAP_SYNC,
        "table": GAIA_TABLE,
        "query": query,
        "center_icrs_deg": [float(ra_deg), float(dec_deg)],
        "radius_deg": float(radius_deg),
        "rows": len(rows),
        "row_limit": row_limit,
        "possibly_truncated": len(rows) >= row_limit,
        "elapsed_seconds": round(time.monotonic() - started, 2),
        "propagation_note": "Gaia positions refer to ref_epoch in Julian years TCB; source positions are propagated to the observation epoch before matching.",
    }
    return rows, diagnostics


def _number(row: dict, name: str) -> float | None:
    try:
        value = float(row.get(name))
        return value if np.isfinite(value) else None
    except (TypeError, ValueError):
        return None


def _corr(row: dict, name: str) -> float:
    value = _number(row, name)
    return float(np.clip(value, -1.0, 1.0)) if value is not None else 0.0


def propagate_gaia_source(row: dict, target_mjd: float, max_no_pm_delta_years: float = 0.5) -> dict:
    """Propagate Gaia DR3 ICRS position using pmra=mu_alpha*cos(dec).

    Positional covariance is propagated with the catalogue standard errors and
    available correlations in the local tangent plane. Radial velocity and
    perspective acceleration are not modelled and are called out in metadata.
    """
    ra, dec = _number(row, "ra"), _number(row, "dec")
    ref_epoch = _number(row, "ref_epoch")
    if ra is None or dec is None or ref_epoch is None:
        return {"status": "unavailable", "reason": "Gaia source lacks a valid position or ref_epoch."}
    target_year = Time(float(target_mjd), format="mjd", scale="utc").tcb.jyear
    delta_years = float(target_year - ref_epoch)
    pmra, pmdec = _number(row, "pmra"), _number(row, "pmdec")
    if pmra is None or pmdec is None:
        if abs(delta_years) > float(max_no_pm_delta_years):
            return {
                "status": "unavailable",
                "reason": "Gaia source has no complete proper-motion solution for this epoch baseline.",
                "ref_epoch": ref_epoch,
                "target_jyear_tcb": target_year,
                "delta_years": delta_years,
            }
        coord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")
        pm_available = False
    else:
        coord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame="icrs")
        pm_total = math.hypot(pmra, pmdec)
        if pm_total > 0 and delta_years != 0:
            position_angle = math.atan2(pmra, pmdec) * u.rad
            if delta_years < 0:
                position_angle += 180 * u.deg
            coord = coord.directional_offset_by(position_angle, pm_total * abs(delta_years) * u.mas)
        pm_available = True

    ra_error = _number(row, "ra_error") or 0.0
    dec_error = _number(row, "dec_error") or 0.0
    pmra_error = _number(row, "pmra_error") or 0.0
    pmdec_error = _number(row, "pmdec_error") or 0.0
    var_x = ra_error**2 + delta_years**2 * pmra_error**2 + 2 * delta_years * _corr(row, "ra_pmra_corr") * ra_error * pmra_error
    var_y = dec_error**2 + delta_years**2 * pmdec_error**2 + 2 * delta_years * _corr(row, "dec_pmdec_corr") * dec_error * pmdec_error
    cov_xy = (
        _corr(row, "ra_dec_corr") * ra_error * dec_error
        + delta_years * _corr(row, "ra_pmdec_corr") * ra_error * pmdec_error
        + delta_years * _corr(row, "dec_pmra_corr") * dec_error * pmra_error
        + delta_years**2 * _corr(row, "pmra_pmdec_corr") * pmra_error * pmdec_error
    )
    sigma_arcsec = math.sqrt(max(var_x + var_y, 0.0)) / 1000.0
    return {
        "status": "ok",
        "ra_deg": float(coord.icrs.ra.deg),
        "dec_deg": float(coord.icrs.dec.deg),
        "sigma_arcsec": sigma_arcsec,
        "sigma_ra_star_mas": math.sqrt(max(var_x, 0.0)),
        "sigma_dec_mas": math.sqrt(max(var_y, 0.0)),
        "position_pm_covariance_mas2": cov_xy,
        "ref_epoch": ref_epoch,
        "target_jyear_tcb": float(target_year),
        "delta_years": delta_years,
        "pmra_mas_per_year": pmra,
        "pmdec_mas_per_year": pmdec,
        "proper_motion_available": pm_available,
        "parallax_mas": _number(row, "parallax"),
        "propagation_model": "uniform linear proper motion on the celestial sphere; parallax and radial-velocity perspective acceleration are not applied",
    }


def _candidate_value(candidate: Any, name: str, default=None):
    if isinstance(candidate, dict):
        return candidate.get(name, default)
    return getattr(candidate, name, default)


def _candidate_epoch(candidate: Any, fallback_mjd: float | None) -> float | None:
    direct = _candidate_value(candidate, "mjd")
    if direct is not None:
        return _number({"v": direct}, "v")
    metadata = _candidate_value(candidate, "metadata", {}) or {}
    for epoch_name in ("epoch1", "epoch2"):
        epoch = metadata.get(epoch_name, {}) if isinstance(metadata, dict) else {}
        if epoch.get("mjd") is not None:
            return _number({"v": epoch.get("mjd")}, "v")
    return fallback_mjd


def _candidate_epoch_positions(candidate: Any, fallback_mjd: float | None) -> list[dict]:
    """Collect per-epoch positions where those measurements are available."""
    if isinstance(candidate, dict):
        detections = candidate.get("detections")
        if isinstance(detections, list):
            rows = [dict(row) for row in detections if isinstance(row, dict)]
            return [row for row in rows if row.get("mjd") is not None and row.get("ra_deg") is not None and row.get("dec_deg") is not None]
        if candidate.get("mjd") is not None:
            return [candidate]
    metadata = _candidate_value(candidate, "metadata", {}) or {}
    rows = []
    if isinstance(metadata, dict):
        for epoch_name in ("epoch1", "epoch2"):
            epoch = metadata.get(epoch_name, {})
            if not isinstance(epoch, dict) or epoch.get("mjd") is None:
                continue
            ra = epoch.get("ra_deg", _candidate_value(candidate, "ra_deg"))
            dec = epoch.get("dec_deg", _candidate_value(candidate, "dec_deg"))
            if ra is not None and dec is not None:
                rows.append({**epoch, "ra_deg": ra, "dec_deg": dec, "epoch_name": epoch_name})
    if rows:
        return rows
    ra = _candidate_value(candidate, "ra_deg", _candidate_value(candidate, "ra"))
    dec = _candidate_value(candidate, "dec_deg", _candidate_value(candidate, "dec"))
    epoch_mjd = _candidate_epoch(candidate, fallback_mjd)
    return [{"ra_deg": ra, "dec_deg": dec, "mjd": epoch_mjd}] if ra is not None and dec is not None and epoch_mjd is not None else []


def _candidate_sigma(candidate: Any, fallback_arcsec: float) -> tuple[float, str]:
    direct = _candidate_value(candidate, "position_sigma_arcsec")
    value = _number({"v": direct}, "v") if direct is not None else None
    if value is not None and value > 0:
        return value, "candidate-provided astrometric uncertainty"
    quality = _candidate_value(candidate, "quality", {}) or {}
    if isinstance(quality, dict):
        direct = quality.get("position_sigma_arcsec")
        value = _number({"v": direct}, "v") if direct is not None else None
        if value is not None and value > 0:
            return value, "candidate quality metadata"
    return max(float(fallback_arcsec), 0.1), "conservative half-pixel/IRSA position-resolution floor"


def cross_match_candidates(
    candidates: list[Any],
    rows: list[dict],
    *,
    fallback_mjd: float | None = None,
    fallback_position_sigma_arcsec: float = 3.1,
    sigma_radius: float = 3.0,
    minimum_radius_arcsec: float = 1.0,
    maximum_radius_arcsec: float = 15.0,
    max_no_pm_delta_years: float = 0.5,
) -> dict[str, dict]:
    """Match candidate detections after propagating Gaia to each epoch."""
    results: dict[str, dict] = {}
    for index, candidate in enumerate(candidates, start=1):
        candidate_id = str(_candidate_value(candidate, "candidate_id", f"CAND-{index:03d}"))
        epoch_positions = _candidate_epoch_positions(candidate, fallback_mjd)
        if not epoch_positions:
            results[candidate_id] = {"status": "unavailable", "matches": [], "detail": "Candidate position or observation epoch is unavailable."}
            continue
        grouped: dict[str, dict] = {}
        epoch_match_rows = []
        no_pm_skipped = 0
        for epoch_position in epoch_positions:
            try:
                ra = float(epoch_position["ra_deg"])
                dec = float(epoch_position["dec_deg"])
                epoch_mjd = float(epoch_position["mjd"])
                if not (0.0 <= ra < 360.0 and -90.0 <= dec <= 90.0 and np.isfinite(epoch_mjd)):
                    continue
            except (KeyError, TypeError, ValueError):
                continue
            target = SkyCoord(ra * u.deg, dec * u.deg, frame="icrs")
            epoch_sigma = _number({"v": epoch_position.get("position_sigma_arcsec")}, "v")
            if epoch_sigma is not None and epoch_sigma > 0:
                candidate_sigma, sigma_source = epoch_sigma, "per-epoch candidate astrometric uncertainty"
            else:
                candidate_sigma, sigma_source = _candidate_sigma(candidate, fallback_position_sigma_arcsec)
            for source in rows:
                position = propagate_gaia_source(source, epoch_mjd, max_no_pm_delta_years=max_no_pm_delta_years)
                if position.get("status") != "ok":
                    if "proper-motion" in str(position.get("reason", "")):
                        no_pm_skipped += 1
                    continue
                coordinate = SkyCoord(position["ra_deg"] * u.deg, position["dec_deg"] * u.deg, frame="icrs")
                separation = float(target.separation(coordinate).arcsec)
                combined_sigma = math.hypot(candidate_sigma, float(position["sigma_arcsec"]))
                radius = min(float(maximum_radius_arcsec), max(float(minimum_radius_arcsec), float(sigma_radius) * combined_sigma))
                if separation > radius:
                    continue
                source_id = str(source.get("source_id", source.get("SOURCE_ID", "")))
                match = {
                    **source,
                    "propagated_ra_deg": position["ra_deg"],
                    "propagated_dec_deg": position["dec_deg"],
                    "propagated_sigma_arcsec": position["sigma_arcsec"],
                    "propagated_from_ref_epoch": position["ref_epoch"],
                    "propagated_to_jyear_tcb": position["target_jyear_tcb"],
                    "delta_years": position["delta_years"],
                    "proper_motion_available": position["proper_motion_available"],
                    "propagation_model": position["propagation_model"],
                    "separation_arcsec": separation,
                    "combined_sigma_arcsec": combined_sigma,
                    "match_radius_arcsec": radius,
                    "separation_over_sigma": separation / combined_sigma if combined_sigma > 0 else None,
                    "matched_observation_id": epoch_position.get("observation_id", epoch_position.get("id")),
                    "matched_epoch_name": epoch_position.get("epoch_name"),
                    "matched_epoch_mjd": epoch_mjd,
                    "candidate_position_icrs_deg": [ra, dec],
                }
                epoch_match_rows.append(match)
                previous = grouped.get(source_id)
                if previous is None:
                    grouped[source_id] = {**match, "matched_epochs": [match]}
                else:
                    previous["matched_epochs"].append(match)
                    if match["separation_over_sigma"] < previous["separation_over_sigma"]:
                        retained_epochs = previous["matched_epochs"]
                        grouped[source_id] = {**match, "matched_epochs": retained_epochs}
        propagated = sorted(grouped.values(), key=lambda row: (row["separation_over_sigma"], row["separation_arcsec"]))
        candidate_sigma, sigma_source = _candidate_sigma(candidate, fallback_position_sigma_arcsec)
        if propagated:
            results[candidate_id] = {
                "status": "matched",
                "matches": propagated[:3],
                "epoch_matches": epoch_match_rows,
                "detail": f"{len(propagated)} Gaia DR3 source(s) match one or more measured epochs after proper-motion propagation.",
                "candidate_epochs_mjd": [float(row["mjd"]) for row in epoch_positions if row.get("mjd") is not None],
                "candidate_position_sigma_arcsec": candidate_sigma,
                "candidate_position_sigma_source": sigma_source,
                "match_radius_sigma": float(sigma_radius),
                "no_proper_motion_rows_skipped": no_pm_skipped,
            }
        else:
            results[candidate_id] = {
                "status": "no match",
                "matches": [],
                "epoch_matches": [],
                "detail": "No propagated Gaia DR3 source fell within the uncertainty-scaled radius at the available candidate epochs; this is not evidence of a new object.",
                "candidate_epochs_mjd": [float(row["mjd"]) for row in epoch_positions if row.get("mjd") is not None],
                "candidate_position_sigma_arcsec": candidate_sigma,
                "candidate_position_sigma_source": sigma_source,
                "no_proper_motion_rows_skipped": no_pm_skipped,
            }
    return results


def cross_match_with_gaia(
    candidates: list[Any],
    *,
    fallback_mjd: float | None,
    center_ra_deg: float,
    center_dec_deg: float,
    fallback_position_sigma_arcsec: float,
    cfg: dict | None = None,
) -> tuple[dict[str, dict], dict]:
    settings = {
        "max_rows": 10000,
        "query_min_radius_deg": 0.05,
        "query_motion_buffer_mas_per_year": 10000.0,
        "query_max_radius_deg": 2.0,
        "match_sigma_radius": 3.0,
        "match_min_radius_arcsec": 1.0,
        "match_max_radius_arcsec": 15.0,
        "max_no_pm_delta_years": 0.5,
    }
    settings.update(cfg or {})
    epochs = [row.get("mjd") for candidate in candidates for row in _candidate_epoch_positions(candidate, fallback_mjd)]
    target_mjds = [float(value) for value in epochs if value is not None]
    if not target_mjds:
        raise GaiaQueryError("No observation epoch is available for Gaia position propagation.")
    max_delta_years = max(abs(Time(value, format="mjd", scale="utc").tcb.jyear - 2016.0) for value in target_mjds)
    motion_margin_deg = float(settings["query_motion_buffer_mas_per_year"]) * max_delta_years / 3_600_000.0
    radius = min(float(settings["query_max_radius_deg"]), max(float(settings["query_min_radius_deg"]), motion_margin_deg + 0.02))
    rows, diagnostics = fetch_gaia_dr3(
        float(center_ra_deg), float(center_dec_deg), radius, int(settings["max_rows"])
    )
    matches = cross_match_candidates(
        candidates,
        rows,
        fallback_mjd=fallback_mjd,
        fallback_position_sigma_arcsec=fallback_position_sigma_arcsec,
        sigma_radius=float(settings["match_sigma_radius"]),
        minimum_radius_arcsec=float(settings["match_min_radius_arcsec"]),
        maximum_radius_arcsec=float(settings["match_max_radius_arcsec"]),
        max_no_pm_delta_years=float(settings["max_no_pm_delta_years"]),
    )
    diagnostics["proper_motion_search_margin_deg"] = motion_margin_deg
    diagnostics["query_radius_capped"] = radius >= float(settings["query_max_radius_deg"])
    diagnostics["position_epoch_basis"] = "per-epoch measured candidate coordinates where available; otherwise selected pair epoch"
    diagnostics["maximum_epoch_delta_from_gaia_ref_epoch_years"] = float(max_delta_years)
    diagnostics["match_method"] = "propagate Gaia DR3 position with pmra=mu_alpha*cos(dec), then compare against a radius based on combined positional uncertainty"
    return matches, diagnostics

