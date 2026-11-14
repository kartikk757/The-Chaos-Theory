"""Optional SPHEREx HiPS tile discovery and Aladin Lite embed helpers."""
from __future__ import annotations

from typing import Any

import numpy as np
import requests
import streamlit as st

import config
from modules.network_guard import wait_for_service


class HipsUnavailable(RuntimeError):
    """Raised when HiPS metadata, optional HEALPix code, or tiles are unavailable."""


def parse_properties(text: str) -> dict[str, str]:
    properties = {}
    for line in str(text).splitlines():
        line = line.strip()
        if not line or line.startswith("#") or "=" not in line:
            continue
        key, value = line.split("=", 1)
        properties[key.strip()] = value.strip()
    return properties


@st.cache_data(ttl=3600, max_entries=2, show_spinner=False)
def fetch_properties() -> dict[str, str]:
    wait_for_service("irsa_hips", 0.5)
    try:
        response = requests.get(f"{config.SPHEREX_HIPS_ROOT}/properties", timeout=(10, 30))
        response.raise_for_status()
    except requests.RequestException as exc:
        raise HipsUnavailable(f"SPHEREx HiPS properties are unavailable: {exc}") from exc
    return parse_properties(response.text)


def healpix_nested_index(ra_deg: float, dec_deg: float, order: int) -> int:
    """Return a NESTED HEALPix index using an installed scientific library."""
    order = int(order)
    if not 0 <= order <= 12:
        raise ValueError("HEALPix order must be between 0 and 12.")
    if not 0 <= float(ra_deg) < 360 or not -90 <= float(dec_deg) <= 90:
        raise ValueError("Position is outside valid ICRS coordinates.")
    nside = 2**order
    try:
        import hpgeom
        return int(hpgeom.angle_to_pixel(nside, float(ra_deg), float(dec_deg), nest=True, lonlat=True))
    except ImportError:
        try:
            from astropy.coordinates import SkyCoord
            from astropy_healpix import HEALPix
            import astropy.units as u
            coordinate = SkyCoord(float(ra_deg) * u.deg, float(dec_deg) * u.deg, frame="icrs")
            return int(HEALPix(nside=nside, order="nested", frame="icrs").lonlat_to_healpix(coordinate.ra, coordinate.dec))
        except ImportError as exc:
            raise HipsUnavailable("HEALPix indexing requires optional hpgeom or astropy-healpix; FITS cube upload remains available.") from exc


def hips_cube_url(ra_deg: float, dec_deg: float, order: int | None = None, properties: dict[str, str] | None = None) -> dict[str, Any]:
    properties = properties or fetch_properties()
    order = int(order if order is not None else properties.get("hips_order", 5))
    minimum = int(properties.get("hips_order_min", 0))
    maximum = int(properties.get("hips_order", 5))
    if order < minimum or order > maximum:
        raise ValueError(f"This HiPS supports orders {minimum}–{maximum}.")
    ipix = healpix_nested_index(ra_deg, dec_deg, order)
    data_root = properties.get("hips_data_url", config.SPHEREX_HIPS_DATA_ROOT).rstrip("/")
    directory = (ipix // 10000) * 10000
    url = f"{data_root}/Norder{order}/Dir{directory}/Npix{ipix}_cube.fits"
    return {
        "url": url,
        "order": order,
        "ipix_nested": ipix,
        "ra_deg": float(ra_deg),
        "dec_deg": float(dec_deg),
        "data_root": data_root,
        "properties": properties,
        "note": "NESTED HEALPix tile lookup; order and data-root version are read from HiPS properties when available. HiPS cube is visualization context, not time-domain evidence.",
    }


@st.cache_data(ttl=3600, max_entries=1, show_spinner=False)
def fetch_cube_bytes(url: str, max_bytes: int = 200 * 1024 * 1024) -> tuple[bytes, dict]:
    wait_for_service("irsa_hips", 0.5)
    try:
        with requests.get(url, timeout=(10, 90), stream=True) as response:
            response.raise_for_status()
            length = response.headers.get("Content-Length")
            if length and int(length) > max_bytes:
                raise HipsUnavailable(f"HiPS cube exceeds the {max_bytes // (1024 * 1024)} MiB download limit.")
            chunks, total = [], 0
            for chunk in response.iter_content(1024 * 1024):
                if not chunk:
                    continue
                total += len(chunk)
                if total > max_bytes:
                    raise HipsUnavailable(f"HiPS cube exceeded the {max_bytes // (1024 * 1024)} MiB download limit.")
                chunks.append(chunk)
    except (requests.RequestException, ValueError) as exc:
        raise HipsUnavailable(f"HiPS cube download failed: {exc}") from exc
    payload = b"".join(chunks)
    import hashlib
    return payload, {"url": url, "bytes": total, "sha256": hashlib.sha256(payload).hexdigest(), "status_code": response.status_code}


def aladin_lite_html(ra_deg: float, dec_deg: float, fov_deg: float, hips_url: str = config.SPHEREX_HIPS_ROOT) -> str:
    """Return a constrained Aladin Lite v3 embed for the official SPHEREx HiPS."""
    ra = float(np.clip(float(ra_deg), 0.0, 359.999999))
    dec = float(np.clip(float(dec_deg), -90.0, 90.0))
    fov = float(np.clip(float(fov_deg), 0.02, 180.0))
    safe_url = str(hips_url).replace("&", "&amp;").replace("'", "")
    return f"""<!doctype html><html><head><meta name="viewport" content="width=device-width, initial-scale=1.0" />
<style>html,body,#aladin-lite-div{{margin:0;width:100%;height:500px;background:#07111d;}}</style></head><body>
<div id="aladin-lite-div"></div>
<script src="https://aladin.cds.unistra.fr/AladinLite/api/v3/latest/aladin.js" charset="utf-8"></script>
<script>A.init.then(()=>{{const aladin=A.aladin('#aladin-lite-div',{{target:'{ra:.7f} {dec:.7f}',fov:{fov:.5f},cooFrame:'equatorial',showReticle:true,showCooGridControl:true,showSimbadPointerControl:false}});const survey=aladin.newImageSurvey('{safe_url}');aladin.setBaseImageLayer(survey);}});</script>
</body></html>"""
