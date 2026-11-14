"""NASA/IPAC IRSA SPHEREx archive client.

Verified against the official SPHEREx Data Access guide, IRSA TAP schema and
IRSA SPHEREx cutout tutorial on 2026-11-13. This module never invents metadata:
requests either return real archive rows or raise a visible ArchiveError.
"""
from __future__ import annotations

import csv
import hashlib
import io
import json
import re
import time
from datetime import date
from html import escape
from pathlib import Path
from typing import Any
from urllib.parse import urlparse

import astropy.units as u
import requests
import streamlit as st
from astropy.coordinates import Angle, SkyCoord
from astropy.time import Time

import config


class ArchiveError(RuntimeError):
    """Raised when a documented archive service cannot provide valid data."""


def _cache_path(key_material: dict) -> Path:
    serialized = json.dumps(key_material, sort_keys=True, default=str)
    digest = hashlib.sha256(serialized.encode("utf-8")).hexdigest()[:24]
    directory = config.CACHE_DIR / "http"
    directory.mkdir(parents=True, exist_ok=True)
    return directory / f"{digest}.json"


def _http_get_text(url: str, params: dict, timeout: int, cache_ttl: int = config.CACHE_TTL_SECONDS) -> tuple[str, dict]:
    key = {"url": url, "params": params}
    path = _cache_path(key)
    now = time.time()
    if cache_ttl > 0 and path.exists():
        try:
            payload = json.loads(path.read_text(encoding="utf-8"))
            if now - float(payload.get("timestamp", 0)) < cache_ttl:
                return payload["text"], {"cache": "file", "age_seconds": round(now - payload["timestamp"], 1)}
        except (OSError, ValueError, KeyError, TypeError):
            pass
    started = time.time()
    response = requests.get(url, params=params, timeout=timeout)
    elapsed = round(time.time() - started, 3)
    try:
        response.raise_for_status()
    except requests.HTTPError as exc:
        raise ArchiveError(
            f"{urlparse(url).netloc} returned HTTP {response.status_code}: {response.text[:300]}"
        ) from exc
    text = response.text
    if not text.strip():
        raise ArchiveError(f"{url} returned an empty response")
    if text.lstrip().startswith(("<html", "<HTML", "<VOTABLE")) and "csv" not in url:
        raise ArchiveError(f"{url} returned markup instead of the requested CSV table")
    if cache_ttl > 0:
        path.write_text(
            json.dumps({"timestamp": now, "url": url, "params": params, "text": text}),
            encoding="utf-8",
        )
    return text, {"cache": "network", "elapsed_seconds": elapsed, "status_code": response.status_code, "bytes": len(text.encode("utf-8"))}


def parse_coordinate(value: Any, kind: str = "ra") -> float:
    """Parse degrees or sexagesimal ICRS text into degrees."""
    text = str(value).strip()
    if not text:
        raise ValueError(f"{kind.upper()} cannot be empty")
    try:
        degrees = float(text)
    except ValueError:
        unit = u.hourangle if kind.lower() == "ra" else u.deg
        try:
            degrees = float(Angle(text, unit=unit).to_value(u.deg))
        except Exception as exc:
            raise ValueError(f"Cannot parse {kind.upper()} coordinate '{text}'") from exc
    if kind.lower() == "ra":
        degrees = degrees % 360.0
    elif not -90.0 <= degrees <= 90.0:
        raise ValueError(f"Declination must be between -90° and +90°, got {degrees}°")
    return float(degrees)


def mjd_from_date(value: Any) -> float:
    try:
        return float(Time(str(value), format="iso").mjd)
    except Exception as exc:
        raise ValueError(f"Cannot parse date '{value}'") from exc


def mjd_to_iso(value: Any) -> str:
    if value in (None, ""):
        return ""
    try:
        return Time(float(value), format="mjd").isot[:19]
    except Exception:
        return ""


def _float_or_none(value: Any) -> float | None:
    if value is None or str(value).strip() == "":
        return None
    try:
        result = float(value)
        return result if result == result and abs(result) != float("inf") else None
    except (TypeError, ValueError):
        return None


def _int_or_none(value: Any) -> int | None:
    number = _float_or_none(value)
    return int(number) if number is not None else None


def validate_query(
    ra_deg: Any,
    dec_deg: Any,
    radius_deg: Any,
    size_deg: Any,
    start_date: Any,
    end_date: Any,
    band: str,
    release: str,
    survey: str,
    max_rows: Any,
) -> dict:
    ra_value = parse_coordinate(ra_deg, "ra")
    dec_value = parse_coordinate(dec_deg, "dec")
    radius = float(radius_deg)
    size = float(size_deg)
    if not 0.005 <= radius <= 0.25:
        raise ValueError("Search radius must be between 0.005° and 0.25° to keep cutouts small")
    if not 0.01 <= size <= 0.5:
        raise ValueError("Cutout width must be between 0.01° and 0.5°")
    if abs(size - 2 * radius) > 0.05:
        # Not an error: documented because users may deliberately choose a different cutout width.
        pass
    start_mjd = mjd_from_date(start_date)
    end_mjd = mjd_from_date(end_date) + 1.0  # include the whole end calendar day
    if end_mjd <= start_mjd:
        raise ValueError("End date must be on or after start date")
    if band != "ALL" and band not in config.DETECTOR_BANDS:
        raise ValueError(f"Unknown SPHEREx detector band: {band}")
    if release not in {"all", "dr1", "qr3"}:
        raise ValueError(f"Unknown release selection: {release}")
    if survey not in {"wide", "deep", "all"}:
        raise ValueError(f"Unknown survey selection: {survey}")
    rows = int(max_rows)
    if not 1 <= rows <= 300:
        raise ValueError("Maximum observation rows must be between 1 and 300")
    return {
        "ra_deg": ra_value,
        "dec_deg": dec_value,
        "radius_deg": radius,
        "size_deg": size,
        "start_date": str(start_date),
        "end_date": str(end_date),
        "start_mjd": start_mjd,
        "end_mjd": end_mjd,
        "band": band,
        "release": release,
        "survey": survey,
        "max_rows": rows,
    }


def _derive_release(raw_uri: str) -> str:
    match = re.search(r"/spherex/([^/]+)/", "/" + raw_uri.lstrip("/"))
    return match.group(1) if match else "unknown"


def _derive_obs_id(raw_uri: str, publisher_did: str = "") -> str:
    match = re.search(r"\?([^/]+)/", publisher_did or "")
    if match:
        return match.group(1)
    filename = Path(raw_uri).name
    match = re.search(r"level2_(.+?)D\d+_spx", filename)
    return match.group(1) if match else filename.replace(".fits", "")


@st.cache_data(ttl=config.CACHE_TTL_SECONDS, max_entries=32, show_spinner=False)
def query_tap_observations(
    ra_deg: float,
    dec_deg: float,
    size_deg: float,
    start_mjd: float,
    end_mjd: float,
    band: str,
    release: str,
    max_rows: int,
) -> tuple[list[dict], dict]:
    """Query spherex.plane + spherex.artifact and construct documented cutout URLs."""
    filters = [
        f"1 = CONTAINS(POINT('ICRS', {ra_deg:.8f}, {dec_deg:.8f}), p.poly)",
        f"p.time_bounds_lower >= {start_mjd:.6f}",
        f"p.time_bounds_lower < {end_mjd:.6f}",
        "a.producttype = 'science'",
        "a.contenttype = 'image/fits'",
    ]
    if band != "ALL":
        filters.append(f"p.energy_bandpassname = '{band}'")
    if release in {"dr1", "qr3"}:
        filters.append(f"a.uri LIKE 'ibe/data/spherex/{release}/%'")
    where = "\n        AND ".join(filters)
    cutout_suffix = f"?center={ra_deg:.8f},{dec_deg:.8f}d&size={size_deg:.5f}"
    query = f"""
SELECT TOP {max_rows}
    'https://irsa.ipac.caltech.edu/' || a.uri || '{cutout_suffix}' AS uri,
    a.uri AS raw_uri,
    p.time_bounds_lower, p.time_bounds_upper, p.time_exposure,
    p.energy_bandpassname, p.energy_bounds_lower, p.energy_bounds_upper,
    p.calibrationlevel, p.quality_flag, p.provenance_version, p.productid,
    p.obsid, p.planeid, p.obs_publisher_did, p.position_resolution,
    p.position_samplesize, p.metrics_background, p.metrics_backgroundstddev,
    a.contentlength, a.producttype, a.contenttype, a.contentchecksum,
    a.ipac_pub_date, a.lastmodified
FROM spherex.artifact a
JOIN spherex.plane p ON a.planeid = p.planeid
WHERE {where}
ORDER BY p.time_bounds_lower
""".strip()
    text, transport = _http_get_text(
        config.IRSA_TAP_SYNC,
        {"QUERY": query, "FORMAT": "csv", "MAXREC": max_rows},
        timeout=config.ARCHIVE_TIMEOUT_SECONDS,
    )
    reader = csv.DictReader(io.StringIO(text))
    raw_rows = list(reader)
    if reader.fieldnames and any("error" in name.lower() for name in reader.fieldnames):
        raise ArchiveError(f"IRSA TAP returned an error table: {text[:500]}")
    records = [build_observation_record(row) for row in raw_rows]
    diagnostics = {
        "endpoint": config.IRSA_TAP_SYNC,
        "query": query,
        "rows": len(records),
        "transport": transport,
        "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    return records, diagnostics


def build_observation_record(row: dict) -> dict:
    raw_uri = row.get("raw_uri", "") or ""
    full_fits_url = row.get("uri", "") or ""
    # The TAP query aliases a.uri to include the cutout suffix; reconstruct the exact direct URL.
    full_fits_url = full_fits_url.split("?center=")[0] if "?center=" in full_fits_url else full_fits_url
    if not full_fits_url and raw_uri:
        full_fits_url = config.IRSA_CUTOUT_ROOT + raw_uri.lstrip("/")
    cutout_url = row.get("uri", "") or ""
    if not cutout_url and full_fits_url:
        cutout_url = full_fits_url
    mjd = _float_or_none(row.get("time_bounds_lower"))
    publisher = str(row.get("obs_publisher_did", "") or "")
    band = str(row.get("energy_bandpassname", "") or "unknown")
    release = _derive_release(raw_uri)
    detector = band.split("-D")[-1] if "-D" in band else ""
    obs_id = _derive_obs_id(raw_uri, publisher)
    date_utc = mjd_to_iso(mjd)
    stable = hashlib.sha256((raw_uri + str(mjd)).encode("utf-8")).hexdigest()[:12].upper()
    lower_m = _float_or_none(row.get("energy_bounds_lower"))
    upper_m = _float_or_none(row.get("energy_bounds_upper"))
    quality = str(row.get("quality_flag", "") or "").strip()
    warnings = []
    if not quality:
        warnings.append("Archive quality_flag was not supplied for this plane.")
    if not publisher:
        warnings.append("Observation ID derived from the documented filename convention.")
    record = {
        "id": f"SPX-{stable}",
        "label": f"{date_utc[:10] or 'Unknown date'} · {band} · {obs_id}",
        "origin": "NASA/IPAC IRSA · SPHEREx TAP",
        "data_mode": config.DATA_MODES["LIVE"],
        "ra_deg": None,  # enriched by SIA2 when a footprint is available
        "dec_deg": None,
        "date_utc": date_utc or None,
        "mjd": mjd,
        "band": band,
        "detector": detector,
        "release": release,
        "obs_id": obs_id,
        "cutout_url": cutout_url,
        "full_fits_url": full_fits_url,
        "content_length": _int_or_none(row.get("contentlength")),
        "calibration_level": _int_or_none(row.get("calibrationlevel")),
        "quality_flag": quality,
        "provenance_version": str(row.get("provenance_version", "") or ""),
        "exposure_s": _float_or_none(row.get("time_exposure")),
        "energy_bounds_lower_m": lower_m,
        "energy_bounds_upper_m": upper_m,
        "wavelength_range_um": (
            round(lower_m * 1e6, 5) if lower_m is not None else None,
            round(upper_m * 1e6, 5) if upper_m is not None else None,
        ),
        "product_type": str(row.get("producttype", "") or ""),
        "content_type": str(row.get("contenttype", "") or ""),
        "plane_id": str(row.get("planeid", "") or ""),
        "obs_publisher_did": publisher,
        "position_resolution_arcsec": _float_or_none(row.get("position_resolution")),
        "position_samplesize_arcsec": _float_or_none(row.get("position_samplesize")),
        "metrics_background": _float_or_none(row.get("metrics_background")),
        "metrics_backgroundstddev": _float_or_none(row.get("metrics_backgroundstddev")),
        "contentchecksum": str(row.get("contentchecksum", "") or ""),
        "ipac_pub_date": str(row.get("ipac_pub_date", "") or ""),
        "is_real_spherex": True,
        "is_demo": False,
        "central_wavelength_um": None,
        "central_bandwidth_um": None,
        "coordinates_source": "requested field centre (no SIA2 footprint match)",
        "metadata": {"tap_row": dict(row), "sia2_match": False},
        "warnings": warnings,
    }
    return record


def _parse_s_region(value: str, center_ra: float) -> list[tuple[float, float]]:
    if not value or not value.upper().startswith("POLYGON"):
        return []
    tokens = value.replace("POLYGON", "").replace("ICRS", "").strip().split()
    try:
        numbers = [float(token) for token in tokens]
    except ValueError:
        return []
    points: list[tuple[float, float]] = []
    for ra, dec in zip(numbers[0::2], numbers[1::2]):
        while ra - center_ra > 180:
            ra -= 360
        while ra - center_ra < -180:
            ra += 360
        points.append((ra, dec))
    return points


def _sia2_collections(release: str, survey: str) -> list[str]:
    releases = ["DR1", "QR3"] if release == "all" else [release.upper()]
    surveys = ["wide", "deep"] if survey == "all" else [survey]
    return [config.SIA2_COLLECTIONS[f"{r} {s}"] for r in releases for s in surveys]


@st.cache_data(ttl=900, max_entries=32, show_spinner=False)
def query_sia2_coverage(
    ra_deg: float,
    dec_deg: float,
    radius_deg: float,
    release: str,
    survey: str,
    band: str,
    max_records: int = 120,
) -> tuple[list[dict], dict]:
    collections = _sia2_collections(release, survey)
    rows: list[dict] = []
    collection_status: dict[str, dict] = {}
    for collection in collections:
        params = {
            "COLLECTION": collection,
            "POS": f"circle {ra_deg:.8f} {dec_deg:.8f} {radius_deg:.6f}",
            "RESPONSEFORMAT": "CSV",
            "MAXREC": int(min(max_records, 250)),
        }
        try:
            text, transport = _http_get_text(config.IRSA_SIA2, params, timeout=config.ARCHIVE_TIMEOUT_SECONDS, cache_ttl=900)
            parsed = list(csv.DictReader(io.StringIO(text)))
            for row in parsed:
                if band != "ALL" and row.get("energy_bandpassname") != band:
                    continue
                s_ra = _float_or_none(row.get("s_ra"))
                s_dec = _float_or_none(row.get("s_dec"))
                row["collection"] = collection
                row["polygon"] = _parse_s_region(str(row.get("s_region", "")), ra_deg)
                row["s_ra_float"] = s_ra
                row["s_dec_float"] = s_dec
                row["t_min_float"] = _float_or_none(row.get("t_min"))
                rows.append(row)
            collection_status[collection] = {"status": "ok", "rows": len(parsed), "transport": transport}
        except Exception as exc:
            collection_status[collection] = {"status": "failed", "error": repr(exc)}
    ok = sum(1 for value in collection_status.values() if value["status"] == "ok")
    status = "ok" if ok == len(collections) else "partial" if ok else "failed"
    detail = f"{ok}/{len(collections)} SIA2 collections responded; {len(rows)} band-filtered footprints."
    coverage_status = {
        "status": status,
        "detail": detail,
        "collections": collection_status,
        "endpoint": config.IRSA_SIA2,
        "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }
    return rows, coverage_status


def enrich_observations_with_sia2(observations: list[dict], coverage_rows: list[dict], ra_deg: float, dec_deg: float) -> list[dict]:
    by_key: dict[tuple[str, str], list[dict]] = {}
    for row in coverage_rows:
        key = (str(row.get("access_url", "")), str(row.get("energy_bandpassname", "")))
        by_key.setdefault(key, []).append(row)
    for obs in observations:
        obs["ra_deg"] = ra_deg
        obs["dec_deg"] = dec_deg
        candidates = by_key.get((obs.get("full_fits_url", ""), obs.get("band", "")), [])
        match = None
        if candidates:
            if obs.get("mjd") is not None:
                match = min(candidates, key=lambda r: abs((r.get("t_min_float") or 0.0) - float(obs["mjd"])))
            else:
                match = candidates[0]
        if match:
            s_ra = match.get("s_ra_float")
            s_dec = match.get("s_dec_float")
            if s_ra is not None and s_dec is not None:
                obs["ra_deg"], obs["dec_deg"] = s_ra, s_dec
                obs["coordinates_source"] = "SIA2 s_ra/s_dec (observation footprint centre)"
            obs["obs_id"] = str(match.get("obs_id") or obs["obs_id"])
            obs["obs_collection"] = str(match.get("obs_collection") or "")
            obs["s_fov_deg"] = _float_or_none(match.get("s_fov"))
            obs["s_pixel_scale_arcsec"] = _float_or_none(match.get("s_pixel_scale"))
            obs["s_region"] = str(match.get("s_region") or "")
            obs["cloud_access"] = str(match.get("cloud_access") or "")
            if match.get("t_exptime"):
                obs["exposure_s"] = _float_or_none(match.get("t_exptime")) or obs.get("exposure_s")
            obs["metadata"]["sia2_match"] = True
            obs["metadata"]["sia2"] = {k: v for k, v in match.items() if k not in {"polygon", "cloud_access"}}
        else:
            obs["warnings"].append("No SIA2 footprint match; coordinates are the requested field centre, not the observation centre.")
    return observations


def query_archive(query: dict) -> dict:
    validated = validate_query(**query)
    observations, tap_diag = query_tap_observations(
        validated["ra_deg"], validated["dec_deg"], validated["size_deg"],
        validated["start_mjd"], validated["end_mjd"], validated["band"],
        validated["release"], validated["max_rows"],
    )
    coverage_rows, coverage_status = query_sia2_coverage(
        validated["ra_deg"], validated["dec_deg"], validated["radius_deg"],
        validated["release"], validated["survey"], validated["band"],
        max_records=validated["max_rows"],
    )
    observations = enrich_observations_with_sia2(observations, coverage_rows, validated["ra_deg"], validated["dec_deg"])
    archive_status = {
        "status": "online",
        "detail": f"TAP returned {len(observations)} rows; {coverage_status['detail']}",
        "utc": tap_diag["utc"],
    }
    return {
        "validated_query": validated,
        "observations": observations,
        "coverage_rows": coverage_rows,
        "coverage_status": coverage_status,
        "tap_diagnostics": tap_diag,
        "archive_status": archive_status,
    }


@st.cache_data(ttl=300, max_entries=4, show_spinner=False)
def check_archive_status() -> dict:
    started = time.time()
    results: dict[str, dict] = {}
    try:
        text, transport = _http_get_text(
            config.IRSA_TAP_SYNC,
            {"QUERY": "SELECT TOP 1 planeid FROM spherex.plane", "FORMAT": "csv"},
            timeout=30,
            cache_ttl=0,
        )
        rows = list(csv.DictReader(io.StringIO(text)))
        results["tap"] = {"status": "ok", "rows": len(rows), "transport": transport}
    except Exception as exc:
        results["tap"] = {"status": "failed", "error": repr(exc)}
    try:
        text, transport = _http_get_text(
            config.IRSA_SIA2,
            {"COLLECTION": "spherex_qr3", "POS": f"circle {config.DEFAULT_QUERY['ra_deg']} {config.DEFAULT_QUERY['dec_deg']} 0.01", "RESPONSEFORMAT": "CSV", "MAXREC": 1},
            timeout=30,
            cache_ttl=0,
        )
        rows = list(csv.DictReader(io.StringIO(text)))
        results["sia2"] = {"status": "ok", "rows": len(rows), "transport": transport}
    except Exception as exc:
        results["sia2"] = {"status": "failed", "error": repr(exc)}
    ok = sum(1 for value in results.values() if value["status"] == "ok")
    status = "online" if ok == 2 else "degraded" if ok == 1 else "unavailable"
    return {
        "status": status,
        "detail": f"{ok}/2 official IRSA endpoints verified",
        "services": results,
        "latency_seconds": round(time.time() - started, 3),
        "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
    }


def verify_vo_clients() -> dict:
    """Optional deep verification using the officially recommended VO Python clients."""
    result: dict[str, dict] = {}
    query = "SELECT TOP 1 planeid FROM spherex.plane"
    try:
        import pyvo

        started = time.time()
        service = pyvo.dal.TAPService("https://irsa.ipac.caltech.edu/TAP")
        rows = service.search(query)
        result["pyvo"] = {"status": "ok", "rows": len(rows), "seconds": round(time.time() - started, 2)}
    except Exception as exc:
        result["pyvo"] = {"status": "failed", "error": repr(exc)}
    try:
        from astroquery.ipac.irsa import Irsa

        Irsa.TIMEOUT = 60
        started = time.time()
        table = Irsa.query_tap(query)
        result["astroquery"] = {"status": "ok", "rows": len(table), "seconds": round(time.time() - started, 2)}
    except Exception as exc:
        result["astroquery"] = {"status": "failed", "error": repr(exc)}
    return result
