"""SPHEREx ePSF zone selection and descriptive diagnostics.

This module inspects the official FITS ePSF table. It does not perform PSF
photometry or claim a fitted detection; oversampled models require correct
sub-pixel shifting and integration before native-pixel fitting.
"""
from __future__ import annotations

import math
from pathlib import Path
from typing import Any

import numpy as np
from astropy.io import fits


class EPSFError(RuntimeError):
    """Raised when a usable ePSF zone cannot be read from a FITS product."""


def _column_name(data, name: str) -> str | None:
    names = list(getattr(data, "names", []) or [])
    if not names:
        names = list(getattr(getattr(data, "array", None), "names", []) or [])
    if not names:
        names = list(getattr(getattr(data, "dtype", None), "names", []) or [])
    return next((item for item in names if item.upper() == name.upper()), None)


def _scalar(row, names: list[str], default=None):
    for name in names:
        key = _column_name(row, name)
        if key is None:
            continue
        try:
            value = np.asarray(row[key]).item()
            return value.item() if isinstance(value, np.generic) else value
        except Exception:
            continue
    return default


def epsf_for_position(product, x_pixel: float, y_pixel: float) -> dict[str, Any]:
    """Select the nearest official ePSF zone for an image pixel.

    Parameters are zero-based coordinates in the loaded product array. Any
    centre crop applied by the loader is mapped back to the original cutout,
    then CRPIX1A/CRPIX2A are used to recover the parent-detector coordinates as
    documented by the IRSA SPHEREx ePSF tutorial.
    """
    path = Path(product.path)
    if not path.exists():
        raise EPSFError(f"FITS product is no longer available at {path}.")
    crop_origin = product.metadata.get("crop_origin_pixel") or [0.0, 0.0]
    original_cutout_x = float(x_pixel) + float(crop_origin[0])
    original_cutout_y = float(y_pixel) + float(crop_origin[1])

    try:
        with fits.open(path, memmap=False) as hdul:
            hdu = hdul["EPSF"] if "EPSF" in hdul else None
            if hdu is None:
                if "PSF" in hdul:
                    return {
                        "status": "unsupported_legacy_psf",
                        "detail": "This product has a legacy PSF extension but no EPSF table. The app does not substitute or reinterpret legacy PSFs.",
                        "product_version": str(hdul[0].header.get("VERSION", "unknown")),
                    }
                return {"status": "unavailable", "detail": "No EPSF extension is present in this FITS product."}
            image_header = hdul["IMAGE"].header if "IMAGE" in hdul else hdul[1].header
            parent_x = 1.0 + original_cutout_x - float(image_header["CRPIX1A"])
            parent_y = 1.0 + original_cutout_y - float(image_header["CRPIX2A"])
            data = hdu.data
            if data is None or len(data) == 0:
                return {"status": "unavailable", "detail": "The EPSF extension contains no zone rows."}
            names = {name.upper(): name for name in data.names or []}
            required = {"BINX", "BINY", "XCENTER", "YCENTER", "EPSF"}
            if not required.issubset(names):
                return {
                    "status": "unsupported_table", 
                    "detail": f"EPSF table is missing required columns: {', '.join(sorted(required - set(names)))}.",
                    "columns": list(data.names or []),
                }
            xcentres = np.asarray(data[names["XCENTER"]], dtype=float)
            ycentres = np.asarray(data[names["YCENTER"]], dtype=float)
            finite = np.isfinite(xcentres) & np.isfinite(ycentres)
            if not np.any(finite):
                return {"status": "unavailable", "detail": "EPSF zone centres are not finite."}
            distances = np.hypot(xcentres - parent_x, ycentres - parent_y)
            distances[~finite] = np.inf
            row_index = int(np.argmin(distances))
            row = data[row_index]
            model = np.asarray(row[names["EPSF"]], dtype=float)
            if model.ndim == 1:
                dim = int(round(math.sqrt(model.size)))
                if dim * dim != model.size:
                    return {"status": "unsupported_table", "detail": f"EPSF cell has {model.size} values, not a square grid."}
                model = model.reshape((dim, dim))
            model = np.asarray(model, dtype=float)
            finite_model = model[np.isfinite(model)]
            if finite_model.size == 0:
                return {"status": "unavailable", "detail": "Selected EPSF model contains no finite values."}
            model_sum = float(np.sum(finite_model))
            if model_sum != 0:
                normalized = np.where(np.isfinite(model), model / model_sum, 0.0)
            else:
                normalized = model.copy()
            header = hdu.header
            image_wavelength = product.center_wavelength_um
            wavelength_sample_pixel = [original_cutout_x, original_cutout_y]
            try:
                from modules import fits_processor

                spatial_wcs, _ = fits_processor.build_spatial_wcs(image_header)
                measured, _, _ = fits_processor.spectral_measurement(
                    hdul, image_header, spatial_wcs, original_cutout_x, original_cutout_y
                )
                if measured is not None:
                    image_wavelength = measured
            except Exception:
                pass
            epsf_wavelength = _scalar(row, ["CWAVE"])
            epsf_bandwidth = _scalar(row, ["CBAND"])
            wavelength_delta = (
                abs(float(image_wavelength) - float(epsf_wavelength))
                if image_wavelength is not None and epsf_wavelength is not None else None
            )
            warnings_out = []
            if wavelength_delta is None:
                warnings_out.append("Wavelength consistency could not be checked for this ePSF zone.")
            elif epsf_bandwidth and wavelength_delta > float(epsf_bandwidth) / 2.0:
                warnings_out.append("The selected ePSF zone's CWAVE differs from the image centre wavelength by more than half its CBAND.")
            if int(header.get("OVSMPX", 0) or 0) <= 0 or int(header.get("OVSMPY", 0) or 0) <= 0:
                warnings_out.append("EPSF oversampling metadata is missing or invalid.")
            return {
                "status": "ok",
                "product_version": str(hdul[0].header.get("VERSION", "unknown")),
                "detector": str(header.get("DETECTOR", image_header.get("DETECTOR", "unknown"))),
                "source_pixel_loaded_cutout": [float(x_pixel), float(y_pixel)],
                "crop_origin_pixel": [float(crop_origin[0]), float(crop_origin[1])],
                "source_pixel_original_cutout": [original_cutout_x, original_cutout_y],
                "parent_detector_pixel": [parent_x, parent_y],
                "zone_row": row_index,
                "zone_id": [int(_scalar(row, ["BINX"], -1)), int(_scalar(row, ["BINY"], -1))],
                "zone_center_parent_pixel": [float(xcentres[row_index]), float(ycentres[row_index])],
                "zone_distance_parent_pixel": float(distances[row_index]),
                "epsf_shape_oversampled": list(model.shape),
                "oversampling_x": int(header.get("OVSMPX", 0) or 0),
                "oversampling_y": int(header.get("OVSMPY", 0) or 0),
                "epsf_sum_raw": model_sum,
                "epsf_peak_raw": float(np.nanmax(model)),
                "epsf_min_raw": float(np.nanmin(model)),
                "epsf_nonfinite_fraction": float(np.mean(~np.isfinite(model))),
                "epsf_normalized": normalized,
                "epsf_effective_area": _scalar(row, ["NEFF_MEAN"]),
                "epsf_star_count": _scalar(row, ["NSTAR"]),
                "epsf_center_wavelength_um": float(epsf_wavelength) if epsf_wavelength is not None else None,
                "epsf_bandwidth_um": float(epsf_bandwidth) if epsf_bandwidth is not None else None,
                "image_wavelength_at_source_um": image_wavelength,
                "image_wavelength_sample_pixel": wavelength_sample_pixel,
                "wavelength_delta_um": wavelength_delta,
                "zone_grid_shape": [int(np.nanmax(data[names["BINX"]])) + 1, int(np.nanmax(data[names["BINY"]])) + 1],
                "warnings": warnings_out,
                "interpretation": "Zone-selection and model diagnostics only. No PSF fit, matched-filter detection, or validated photometry was performed.",
                "source": "EPSF binary-table extension from this exact FITS product",
            }
    except EPSFError:
        raise
    except Exception as exc:
        raise EPSFError(f"Could not inspect SPHEREx EPSF extension: {exc}") from exc

