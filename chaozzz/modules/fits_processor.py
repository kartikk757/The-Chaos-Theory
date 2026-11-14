"""FITS header inspection, WCS extraction, FLAGS-bit parsing and validation.

This layer contains no Streamlit state and no archive/network code.
"""
from __future__ import annotations

import logging
import warnings
from typing import Any

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.nddata import Cutout2D
from astropy.wcs import WCS

import config


def header_cards(header: fits.Header) -> list[dict[str, str]]:
    cards = []
    for card in header.cards:
        cards.append({
            "keyword": str(card.keyword),
            "value": str(card.value),
            "comment": str(card.comment),
        })
    return cards


def header_dict(header: fits.Header) -> dict[str, Any]:
    result: dict[str, Any] = {}
    for key, value in header.items():
        if isinstance(value, (np.generic,)):
            result[str(key)] = value.item()
        elif isinstance(value, (str, int, float, bool, type(None))):
            result[str(key)] = value
        else:
            result[str(key)] = str(value)
    return result


def find_extension(hdul: fits.HDUList, name: str) -> fits.ImageHDU | None:
    try:
        return hdul[name]
    except KeyError:
        pass
    for hdu in hdul:
        if str(getattr(hdu, "name", "")).strip().upper() == name.upper():
            return hdu
    return None


def find_image_hdu(hdul: fits.HDUList) -> fits.ImageHDU:
    named = find_extension(hdul, "IMAGE")
    if named is not None and getattr(named, "data", None) is not None and np.asarray(named.data).ndim == 2:
        return named
    for hdu in hdul:
        data = getattr(hdu, "data", None)
        if data is not None and np.asarray(data).ndim == 2:
            return hdu
    raise ValueError("No 2-D IMAGE extension was found in this FITS file")


def parse_flag_definitions(header: fits.Header | None) -> dict[str, dict[str, Any]]:
    """Read MSKN/MSKM/MSKD (DR1/QR3) or MP_* (older QR) definitions from the header."""
    if header is None:
        return {}
    definitions: dict[str, dict[str, Any]] = {}
    for key in list(header.keys()):
        if not key.startswith("MSKN"):
            continue
        index = key[4:]
        name = str(header.get(key, "")).strip()
        if not name:
            continue
        mask = header.get("MSKM" + index, 1 << int(index)) if index.isdigit() else 0
        try:
            mask = int(mask)
        except (TypeError, ValueError):
            continue
        try:
            description = str(header.comments["MSKD" + index])
        except (KeyError, AttributeError):
            description = ""
        definitions[name] = {"bit_mask": mask, "description": description, "source": key}
    if not definitions:
        for key, value in header.items():
            if not key.startswith("MP_"):
                continue
            name = key[3:].strip().upper()
            try:
                mask = int(value)
            except (TypeError, ValueError):
                continue
            definitions[name] = {
                "bit_mask": mask,
                "description": str(header.comments.get(key, "")) if hasattr(header, "comments") else "",
                "source": key,
            }
    return definitions


def build_flag_masks(flags: np.ndarray | None, definitions: dict[str, dict[str, Any]]) -> tuple[np.ndarray, np.ndarray, dict]:
    if flags is None:
        return np.zeros((1, 1), dtype=bool), np.zeros((1, 1), dtype=bool), {"available": False}
    values = np.asarray(flags).astype(np.int64, copy=False)
    bad = np.zeros(values.shape, dtype=bool)
    informational = np.zeros(values.shape, dtype=bool)
    used: list[str] = []
    unknown: list[str] = []
    for name, definition in definitions.items():
        bit = int(definition["bit_mask"])
        hit = (values & bit) != 0
        if name in config.BAD_FLAG_NAMES:
            bad |= hit
            used.append(name)
        elif name in config.INFO_FLAG_NAMES:
            informational |= hit
        else:
            # Conservative: an unrecognised flag is excluded from scientific pixels.
            bad |= hit
            used.append(name)
            unknown.append(name)
    source_mask = np.zeros(values.shape, dtype=bool)
    if "SOURCE" in definitions:
        source_mask = (values & int(definitions["SOURCE"]["bit_mask"])) != 0
    diagnostics = {
        "available": True,
        "bad_flag_names_used": used,
        "unknown_flag_names_masked": unknown,
        "informational_flag_names": sorted(config.INFO_FLAG_NAMES.intersection(definitions)),
        "bad_pixel_fraction": float(np.mean(bad)),
        "source_flag_fraction": float(np.mean(source_mask)),
        "nonzero_flag_fraction": float(np.mean(values != 0)),
    }
    return bad, source_mask, diagnostics


def build_spatial_wcs(header: fits.Header) -> tuple[WCS | None, list[str]]:
    warnings_out: list[str] = []
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            previous = logging.getLogger("astropy").level
            logging.getLogger("astropy").setLevel(logging.ERROR)
            wcs = WCS(header, naxis=2)
            logging.getLogger("astropy").setLevel(previous)
        if not wcs.has_celestial:
            warnings_out.append("FITS WCS has no celestial axes; alignment is unavailable.")
            return None, warnings_out
        return wcs, warnings_out
    except Exception as exc:
        warnings_out.append(f"Spatial WCS could not be constructed: {exc}")
        return None, warnings_out


def pixel_scale_arcsec(wcs: WCS | None) -> float | None:
    if wcs is None:
        return None
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            # Astropy's WCS method returns a Quantity; the utils helper returns a bare
            # ndarray in recent versions, so the WCS method is used here.
            # Compute from the pixel scale matrix (degrees per pixel) because
            # proj_plane_pixel_scales returns bare ndarrays in recent Astropy versions.
            matrix = np.asarray(wcs.celestial.pixel_scale_matrix, dtype=float)
            values = np.sqrt(np.sum(matrix ** 2, axis=0)) * 3600.0
        value = float(np.mean(np.abs(values)))
        return value if np.isfinite(value) and value > 0 else None
    except Exception:
        return None


def center_sky(wcs: WCS | None, shape: tuple[int, int], center_ra: float | None = None, center_dec: float | None = None) -> tuple[float, float] | None:
    """Return the world coordinate at the geometric centre of the FITS array.

    ``center_ra`` and ``center_dec`` are retained for call compatibility, but
    archive/query coordinates must not override the delivered image WCS centre.
    Those coordinates are stored separately as observation metadata.
    """
    if wcs is None:
        return None
    try:
        x = (shape[1] - 1) / 2.0
        y = (shape[0] - 1) / 2.0
        coord = wcs.pixel_to_world(x, y)
        return float(coord.ra.deg), float(coord.dec.deg)
    except Exception:
        return None


def center_pixel(wcs: WCS | None, shape: tuple[int, int], center_ra: float | None = None, center_dec: float | None = None) -> tuple[float, float]:
    if wcs is not None and center_ra is not None and center_dec is not None:
        try:
            coordinate = SkyCoord(ra=center_ra * u.deg, dec=center_dec * u.deg, frame="icrs")
            x, y = wcs.world_to_pixel(coordinate)
            return float(x), float(y)
        except Exception:
            pass
    return (shape[1] - 1) / 2.0, (shape[0] - 1) / 2.0


def wavelength_quantity_to_um(value: Any) -> float | None:
    try:
        if hasattr(value, "to_value"):
            return float(value.to_value(u.micron))
        return float(value)
    except Exception:
        return None


def spectral_measurement(
    hdul: fits.HDUList,
    image_header: fits.Header,
    spatial_wcs: WCS | None,
    x: float,
    y: float,
) -> tuple[float | None, float | None, list[str]]:
    """Measure WCS-WAVE central wavelength and bandwidth at the requested sky pixel."""
    warnings_out: list[str] = []
    if find_extension(hdul, "WCS-WAVE") is None:
        return None, None, ["No WCS-WAVE lookup table; central wavelength is unavailable."]
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            previous = logging.getLogger("astropy").level
            logging.getLogger("astropy").setLevel(logging.ERROR)
            spectral_wcs = WCS(image_header, fobj=hdul, key="W")
            spectral_wcs.sip = None
            result = spectral_wcs.pixel_to_world(x, y)
            logging.getLogger("astropy").setLevel(previous)
        if result is None:
            raise ValueError("spectral WCS returned no coordinate")
        values = list(result) if isinstance(result, (list, tuple)) else [result]
        if len(values) < 2:
            raise ValueError("spectral WCS did not return wavelength and bandwidth")
        wavelength = wavelength_quantity_to_um(values[0])
        bandwidth = wavelength_quantity_to_um(values[1])
        if wavelength is None or not np.isfinite(wavelength):
            raise ValueError("non-finite wavelength")
        return float(wavelength), (float(bandwidth) if bandwidth is not None and np.isfinite(bandwidth) else None), warnings_out
    except Exception as exc:
        warnings_out.append(f"WCS-WAVE could not be evaluated: {exc}")
        return None, None, warnings_out


def crop_to_limit(
    arrays: dict[str, np.ndarray | None],
    wcs: WCS | None,
    max_pixels: int = config.MAX_FITS_PIXELS,
    crop_side: int = config.DEFAULT_CROP_SIDE,
) -> tuple[dict[str, np.ndarray | None], WCS | None, tuple[float, float], list[str]]:
    """Crop very large local FITS files around the field centre for limited hardware."""
    image = arrays.get("image")
    if image is None:
        return arrays, wcs, (0.0, 0.0), []
    if image.size <= max_pixels:
        return arrays, wcs, ((image.shape[1] - 1) / 2.0, (image.shape[0] - 1) / 2.0), []
    side = int(min(crop_side, min(image.shape)))
    position = ((image.shape[1] - 1) / 2.0, (image.shape[0] - 1) / 2.0)
    cutout = Cutout2D(image, position, side, wcs=wcs, mode="partial", fill_value=np.nan)
    center_in_cutout = cutout.position_cutout if hasattr(cutout, "position_cutout") else (side / 2.0, side / 2.0)
    cropped: dict[str, np.ndarray | None] = {"image": np.asarray(cutout.data, dtype=float)}
    for key, array in arrays.items():
        if key == "image" or array is None:
            continue
        if np.asarray(array).shape == image.shape:
            cropped[key] = np.asarray(Cutout2D(np.asarray(array), position, side, mode="partial", fill_value=0).data)
        else:
            cropped[key] = None
    return cropped, cutout.wcs, (float(center_in_cutout[0]), float(center_in_cutout[1])), [
        f"Input exceeded {max_pixels:,} pixels; a {side}×{side} centre crop was used to protect memory. Full-resolution analysis is not represented here."
    ]


def validate_arrays(
    image: np.ndarray,
    flags: np.ndarray | None,
    variance: np.ndarray | None,
    zodi: np.ndarray | None,
    wcs: WCS | None,
    bunit: str,
    bad_mask: np.ndarray,
    flag_definitions: dict[str, dict[str, Any]],
    center_wavelength_um: float | None,
    pixel_scale: float | None,
) -> list[dict[str, str]]:
    checks: list[dict[str, str]] = []
    finite = np.isfinite(image)
    checks.append({
        "name": "2-D IMAGE array",
        "status": "pass" if image.ndim == 2 and image.size > 0 else "fail",
        "detail": f"shape {image.shape}, dtype {image.dtype}, finite fraction {np.mean(finite):.3f}",
    })
    checks.append({
        "name": "IMAGE units",
        "status": "pass" if bunit else "warn",
        "detail": f"BUNIT = {bunit or 'not supplied'}; all differences inherit this unit.",
    })
    checks.append({
        "name": "Celestial WCS",
        "status": "pass" if wcs is not None and wcs.has_celestial else "fail",
        "detail": "valid celestial WCS" if wcs is not None else "WCS unavailable; alignment will be rejected",
    })
    checks.append({
        "name": "Pixel scale",
        "status": "pass" if pixel_scale is not None else "warn",
        "detail": f"{pixel_scale:.3f} arcsec/pixel" if pixel_scale else "could not derive pixel scale",
    })
    for name, array in [("FLAGS", flags), ("VARIANCE", variance), ("ZODI", zodi)]:
        if array is None:
            checks.append({"name": f"{name} extension", "status": "warn", "detail": "not present; related quality information is unavailable"})
        elif np.asarray(array).shape != image.shape:
            checks.append({"name": f"{name} shape", "status": "fail", "detail": f"{np.asarray(array).shape} does not match {image.shape}"})
        else:
            checks.append({"name": f"{name} extension", "status": "pass", "detail": f"shape {np.asarray(array).shape}"})
    if variance is not None:
        positive = np.mean(np.asarray(variance) > 0)
        checks.append({"name": "Variance values", "status": "pass" if positive > 0.5 else "warn", "detail": f"positive fraction {positive:.3f}"})
    checks.append({
        "name": "FLAGS definitions",
        "status": "pass" if flag_definitions else "warn",
        "detail": f"{len(flag_definitions)} named bit planes read from the FITS header; severe-flag fraction {np.mean(bad_mask):.3f}",
    })
    checks.append({
        "name": "WCS-WAVE wavelength",
        "status": "pass" if center_wavelength_um is not None else "warn",
        "detail": f"{center_wavelength_um:.5f} µm at field centre" if center_wavelength_um is not None else "not available; spectral compatibility cannot be fully verified",
    })
    return checks
