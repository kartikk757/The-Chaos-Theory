"""Coordinate roles for archive records and delivered FITS cutouts."""
from __future__ import annotations

import math
import re
from typing import Any
from urllib.parse import parse_qs, urlparse


def _finite_degrees(ra: Any, dec: Any) -> tuple[float, float] | None:
    try:
        ra_value, dec_value = float(ra), float(dec)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(ra_value) or not math.isfinite(dec_value):
        return None
    if not 0.0 <= ra_value < 360.0 or not -90.0 <= dec_value <= 90.0:
        return None
    return ra_value, dec_value


def _url_center(url: str) -> tuple[float, float] | None:
    """Read the requested cutout centre only from an IRSA URL."""
    try:
        parsed = urlparse(str(url))
        if parsed.scheme != "https" or parsed.hostname != "irsa.ipac.caltech.edu":
            return None
        value = parse_qs(parsed.query).get("center", [""])[0]
        match = re.fullmatch(r"\s*([+-]?[\d.]+)\s*,\s*([+-]?[\d.]+)\s*d?\s*", value, re.IGNORECASE)
        if not match:
            return None
        return _finite_degrees(match.group(1), match.group(2))
    except (TypeError, ValueError):
        return None


def coordinate_context(record: dict[str, Any]) -> dict[str, Any]:
    """Normalize coordinate meaning without conflating footprints and cutouts.

    ``ra_deg``/``dec_deg`` remain the catalogue or observation-record position.
    The cutout centre is always measured independently from the FITS WCS.
    """
    source = str(record.get("coordinates_source") or "unspecified")
    role = str(record.get("coordinates_role") or "").strip()
    source_lower = source.lower()
    if not role:
        if "sia2" in source_lower or "footprint" in source_lower:
            role = "parent_observation_footprint_center"
        elif "requested" in source_lower:
            role = "requested_field_center"
        elif "synthetic" in source_lower:
            role = "synthetic_image_reference_center"
        elif "fits wcs" in source_lower:
            role = "fits_wcs_image_center"
        else:
            role = "recorded_observation_position"

    record_position = _finite_degrees(record.get("ra_deg"), record.get("dec_deg"))
    metadata = record.get("metadata") if isinstance(record.get("metadata"), dict) else {}
    requested = metadata.get("requested_field_center")
    requested_position = None
    if isinstance(requested, dict):
        requested_position = _finite_degrees(requested.get("ra_deg"), requested.get("dec_deg"))
    if requested_position is None:
        requested_position = _finite_degrees(record.get("requested_ra_deg"), record.get("requested_dec_deg"))
    if requested_position is None:
        requested_position = _url_center(record.get("cutout_url") or record.get("download_url") or "")
    if requested_position is None and role == "requested_field_center":
        requested_position = record_position

    return {
        "record_position_icrs_deg": list(record_position) if record_position else None,
        "record_position_role": role,
        "record_position_source": source,
        "requested_field_center_icrs_deg": list(requested_position) if requested_position else None,
        "requested_field_center_source": (
            "record metadata" if isinstance(requested, dict) and requested_position
            else "IRSA cutout URL" if requested_position and _url_center(record.get("cutout_url") or record.get("download_url") or "") == requested_position
            else "observation record" if requested_position and role == "requested_field_center"
            else "not available"
        ),
        "coordinate_frame": "ICRS",
    }


def attach_requested_field_center(record: dict[str, Any], ra_deg: float, dec_deg: float) -> dict[str, Any]:
    """Store the query/cutout centre separately from the record position."""
    requested = _finite_degrees(ra_deg, dec_deg)
    if requested is None:
        return record
    metadata = record.setdefault("metadata", {})
    if not isinstance(metadata, dict):
        metadata = {}
        record["metadata"] = metadata
    metadata["requested_field_center"] = {
        "ra_deg": requested[0],
        "dec_deg": requested[1],
        "coordinate_frame": "ICRS",
        "source": "IRSA archive query / cutout request",
    }
    return record

