"""SPHEREx ePSF zone selection, diagnostics, and supplemental source fits."""
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
                "epsf_model_raw": model.copy(),
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
                "interpretation": "Zone-selection and model diagnostics. Fits use the oversampled model at native-pixel centres with its FITS normalization; they are supplemental measurements, not validated sensitivity gains.",
                "source": "EPSF binary-table extension from this exact FITS product",
            }
    except EPSFError:
        raise
    except Exception as exc:
        raise EPSFError(f"Could not inspect SPHEREx EPSF extension: {exc}") from exc


def _detector_number(value: Any) -> int | None:
    digits = "".join(character for character in str(value or "") if character.isdigit())
    try:
        number = int(digits[-1])
    except (TypeError, ValueError, IndexError):
        return None
    return number if 1 <= number <= 6 else None


def _flux_conversion_to_jy(product) -> tuple[float | None, str]:
    """Return a pixel-integrated surface-brightness conversion when defined."""
    unit = "".join(str(product.bunit or "").lower().split())
    if product.wcs is None or not product.wcs.has_celestial:
        return None, f"{product.bunit or 'image unit'} pixel-sum"
    try:
        from astropy import units as u
        from astropy.wcs.utils import proj_plane_pixel_area

        area_sr = float((proj_plane_pixel_area(product.wcs) * u.deg**2).to_value(u.sr))
    except Exception:
        return None, f"{product.bunit or 'image unit'} pixel-sum"
    if "mjy/sr" in unit or "mjy/steradian" in unit:
        return area_sr * 1e6, "Jy"
    if "jy/sr" in unit or "jy/steradian" in unit:
        return area_sr, "Jy"
    if unit == "mjy":
        return 1e6, "Jy"
    if unit == "jy":
        return 1.0, "Jy"
    return None, f"{product.bunit or 'image unit'} pixel-sum"


def fit_source(
    product,
    x_pixel: float,
    y_pixel: float,
    *,
    fit_shape: int = 9,
    background_inner_radius_px: float = 4.0,
    background_outer_radius_px: float = 7.0,
    max_position_shift_px: float = 1.5,
) -> dict[str, Any]:
    """Fit one selected source with the product's spatially appropriate ePSF.

    The FITS EPSF array is passed to Photutils without renormalizing it. The
    returned flux is labelled in integrated image units, or converted to Jy
    only when the FITS BUNIT and celestial pixel area make that conversion
    explicit. This is a supplemental isolated-source fit: it does not deblend
    neighbors, recalibrate the PSF, classify candidates, or claim a sensitivity
    improvement.
    """
    try:
        x0, y0 = float(x_pixel), float(y_pixel)
    except (TypeError, ValueError):
        return {"status": "invalid_position", "detail": "Finite image-pixel coordinates are required."}
    if not np.isfinite(x0) or not np.isfinite(y0):
        return {"status": "invalid_position", "detail": "Finite image-pixel coordinates are required."}
    try:
        fit_shape = int(fit_shape)
        background_inner_radius_px = float(background_inner_radius_px)
        background_outer_radius_px = float(background_outer_radius_px)
        max_position_shift_px = float(max_position_shift_px)
    except (TypeError, ValueError, OverflowError):
        return {"status": "invalid_configuration", "detail": "The fit window and background/position limits must be numeric."}
    if (
        fit_shape < 5
        or fit_shape % 2 != 1
        or not np.isfinite(background_inner_radius_px)
        or not np.isfinite(background_outer_radius_px)
        or background_inner_radius_px < 0
        or background_outer_radius_px <= background_inner_radius_px
        or not np.isfinite(max_position_shift_px)
        or max_position_shift_px <= 0
    ):
        return {"status": "invalid_configuration", "detail": "The PSF fit window must be an odd size of at least 5 pixels."}
    try:
        info = epsf_for_position(product, x0, y0)
    except EPSFError as exc:
        return {"status": "unavailable", "detail": str(exc)}
    if info.get("status") != "ok":
        return {"status": info.get("status", "unavailable"), "detail": info.get("detail", "No compatible ePSF model is available."), "diagnostics": {key: value for key, value in info.items() if key not in {"epsf_model_raw", "epsf_normalized"}}}
    if info.get("warnings"):
        severe_warnings = [item for item in info["warnings"] if "differs from the image" in item or "oversampling metadata" in item]
        if severe_warnings:
            return {"status": "incompatible", "detail": "The selected ePSF is not demonstrably compatible with this image position.", "warnings": list(info["warnings"])}

    image_detector = _detector_number(product.header.get("DETECTOR", ""))
    epsf_detector = _detector_number(info.get("detector", ""))
    if image_detector is not None and epsf_detector is not None and image_detector != epsf_detector:
        return {"status": "incompatible_detector", "detail": f"IMAGE detector {image_detector} does not match ePSF detector {epsf_detector}.", "warnings": list(info.get("warnings", []))}

    image = np.asarray(product.science, dtype=float)
    if image.ndim != 2 or min(image.shape) < fit_shape:
        return {"status": "invalid_input", "detail": "A two-dimensional science image larger than the fit window is required."}
    margin = max(fit_shape // 2 + 1, int(np.ceil(background_outer_radius_px)))
    if x0 < margin or y0 < margin or x0 > image.shape[1] - 1 - margin or y0 > image.shape[0] - 1 - margin:
        return {"status": "edge_source", "detail": "The selected source is too close to the image edge for the fit window and background annulus."}

    try:
        from astropy.stats import sigma_clipped_stats
        from astropy.table import QTable
        from photutils.psf import ImagePSF, PSFPhotometry
    except Exception as exc:
        return {"status": "unavailable", "detail": f"Photutils ePSF fitting is unavailable: {exc}"}

    finite = np.isfinite(image)
    bad_flags = np.asarray(getattr(product, "bad_flag_mask", np.zeros(image.shape, dtype=bool)), dtype=bool)
    if bad_flags.shape != image.shape:
        bad_flags = np.zeros(image.shape, dtype=bool)
    mask = ~finite | bad_flags
    variance = getattr(product, "variance", None)
    if variance is not None and np.asarray(variance).shape == image.shape:
        variance = np.asarray(variance, dtype=float)
        mask |= ~np.isfinite(variance) | (variance <= 0)
        error = np.sqrt(np.where(mask, 1.0, variance))
        error_source = "SPHEREx VARIANCE extension; flagged/non-positive pixels masked"
    else:
        variance = None
        error_source = "robust sigma-clipped image noise; no usable variance extension"
    if np.count_nonzero(~mask) < fit_shape * fit_shape:
        return {"status": "invalid_input", "detail": "Too few unmasked image pixels remain for an ePSF fit."}

    yy, xx = np.indices(image.shape, dtype=float)
    radius = np.hypot(xx - x0, yy - y0)
    annulus = (radius >= background_inner_radius_px) & (radius <= background_outer_radius_px) & ~mask
    annulus_values = image[annulus]
    if annulus_values.size < 8:
        return {"status": "invalid_background", "detail": "Fewer than 8 valid pixels are available in the local background annulus."}
    background = float(np.median(annulus_values))
    fit_image = np.where(mask, 0.0, image - background)
    if variance is None:
        _, _, noise = sigma_clipped_stats(image[~mask], sigma=3.0, maxiters=5)
        if not np.isfinite(noise) or noise <= 0:
            return {"status": "invalid_noise", "detail": "A positive robust image-noise estimate could not be measured."}
        error = np.full(image.shape, float(noise), dtype=float)

    raw_model = np.asarray(info.get("epsf_model_raw"), dtype=float)
    oversampling_x = int(info.get("oversampling_x") or 0)
    oversampling_y = int(info.get("oversampling_y") or 0)
    if raw_model.ndim != 2 or min(raw_model.shape) < 4 or not np.all(np.isfinite(raw_model)) or np.any(raw_model < 0):
        return {"status": "invalid_epsf", "detail": "The selected ePSF model must be a finite, non-negative two-dimensional array."}
    if oversampling_x <= 0 or oversampling_y <= 0 or float(np.sum(raw_model)) <= 0:
        return {"status": "invalid_epsf", "detail": "The selected ePSF model has invalid oversampling or normalization metadata."}
    expected_oversampled_sum = float(oversampling_x * oversampling_y)
    raw_sum = float(np.sum(raw_model))
    normalization_ratio = raw_sum / expected_oversampled_sum
    try:
        psf_model = ImagePSF(
            data=raw_model,
            oversampling=(oversampling_y, oversampling_x),
            x_0=x0,
            y_0=y0,
        )
        initial_flux = float(np.nansum(fit_image[(radius <= 2.5) & ~mask]))
        if not np.isfinite(initial_flux) or initial_flux <= 0:
            initial_flux = max(float(np.nanmax(fit_image[~mask])) * float(np.sum(raw_model)) / max(float(np.nanmax(raw_model)), 1e-30), 1e-12)
        init = QTable({"x": [x0], "y": [y0], "flux": [initial_flux]})
        photometry = PSFPhotometry(
            psf_model,
            fit_shape=fit_shape,
            xy_bounds=float(max_position_shift_px),
            fitter_maxiters=150,
        )
        results = photometry(fit_image, error=error, mask=mask, init_params=init)
        row = results[0]
        fitted_flux = float(row["flux_fit"])
        fitted_x = float(row["x_fit"])
        fitted_y = float(row["y_fit"])
        flux_error_value = row["flux_err"] if "flux_err" in results.colnames else None
        fitted_error = float(flux_error_value) if flux_error_value is not None and np.isfinite(float(flux_error_value)) else None
        fitted_model = photometry.make_model_image(image.shape, psf_shape=fit_shape)
        xfit_i, yfit_i = int(round(fitted_x)), int(round(fitted_y))
        half = fit_shape // 2
        y_slice = slice(max(0, yfit_i - half), min(image.shape[0], yfit_i + half + 1))
        x_slice = slice(max(0, xfit_i - half), min(image.shape[1], xfit_i + half + 1))
        valid_fit = ~mask[y_slice, x_slice]
        residual = fit_image[y_slice, x_slice] - fitted_model[y_slice, x_slice]
        standardized = residual[valid_fit] / error[y_slice, x_slice][valid_fit]
        dof = max(1, standardized.size - 3)
        reduced_chi2 = float(np.sum(standardized**2) / dof)
        photometry_flags = int(row["flags"]) if "flags" in results.colnames else 0
        conversion, flux_unit = _flux_conversion_to_jy(product)
        warnings_out = list(info.get("warnings", []))
        if not 0.8 <= normalization_ratio <= 1.2:
            warnings_out.append(
                f"Finite ePSF sum / oversampling area is {normalization_ratio:.3f}; the fitted amplitude retains the FITS normalization and its total-flux scale is not validated."
            )
            conversion = None
            flux_unit = f"{product.bunit or 'image unit'} ePSF amplitude (scale unvalidated)"
        flux_out = fitted_flux * conversion if conversion is not None else fitted_flux
        error_out = fitted_error * conversion if conversion is not None and fitted_error is not None else fitted_error
        if variance is not None:
            warnings_out.append("Flux errors use per-pixel variance and omit background-estimation and inter-pixel covariance terms.")
        else:
            warnings_out.append("Flux errors use an image-wide robust-noise approximation; background structure and correlated noise are not modelled.")
        warnings_out.append("This isolated-source ePSF fit is supplemental; no candidate classification or sensitivity-improvement claim is derived from it.")
        return {
            "status": "ok",
            "x_fit": fitted_x,
            "y_fit": fitted_y,
            "x_initial": x0,
            "y_initial": y0,
            "flux_fit": float(flux_out),
            "flux_error": float(error_out) if error_out is not None else None,
            "flux_unit": flux_unit,
            "flux_fit_native": fitted_flux,
            "flux_error_native": fitted_error,
            "native_flux_unit": f"{product.bunit or 'image unit'} pixel-sum",
            "background_level": background,
            "background_annulus_px": [float(background_inner_radius_px), float(background_outer_radius_px)],
            "error_source": error_source,
            "fit_shape": [fit_shape, fit_shape],
            "fit_pixels": int(standardized.size),
            "reduced_chi_square": reduced_chi2,
            "photometry_flags": photometry_flags,
            "epsf_zone_id": list(info["zone_id"]),
            "epsf_oversampling_xy": [oversampling_x, oversampling_y],
            "epsf_sum_finite": raw_sum,
            "epsf_expected_oversampled_sum": expected_oversampled_sum,
            "epsf_normalization_ratio": normalization_ratio,
            "epsf_zone_wavelength_um": info.get("epsf_center_wavelength_um"),
            "image_wavelength_at_source_um": info.get("image_wavelength_at_source_um"),
            "wavelength_delta_um": info.get("wavelength_delta_um"),
            "detector": info.get("detector"),
            "product_version": info.get("product_version"),
            "product_sha256": str(product.sha256),
            "pixel_solid_angle_conversion_to_jy": conversion,
            "warnings": warnings_out,
            "interpretation": "Supplemental single-source fit to the exact product ePSF using Photutils ImagePSF and PSFPhotometry; not a calibrated completeness or sensitivity estimate.",
        }
    except Exception as exc:
        return {
            "status": "fit_failed",
            "detail": f"Photutils ePSF fit failed: {exc}",
            "epsf_zone_id": list(info["zone_id"]),
            "product_sha256": str(product.sha256),
            "warnings": list(info.get("warnings", [])),
        }

