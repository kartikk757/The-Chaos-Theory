"""Difference imaging and propagated uncertainty.

D(x, y) = I2(x, y) - I1(x, y) on the common reprojected grid.
Noise comes from the official VARIANCE extensions when available, otherwise
from a robust estimate of the difference image itself.
"""
from __future__ import annotations

import numpy as np
from astropy.stats import mad_std, sigma_clipped_stats

from modules.models import DifferenceResult, NormalizedPair


def global_difference_noise(difference: np.ndarray, valid_mask: np.ndarray, cfg: dict) -> tuple[float, dict]:
    values = difference[valid_mask & np.isfinite(difference)]
    if values.size < 100:
        return float("nan"), {"method": "insufficient pixels", "count": int(values.size)}
    _, median, std = sigma_clipped_stats(values, sigma=cfg["background_sigma"], maxiters=cfg["background_iterations"])
    mad = float(mad_std(values))
    chosen = float(std) if np.isfinite(std) and std > 0 else mad
    return chosen, {
        "method": "sigma-clipped standard deviation of the difference image",
        "count": int(values.size),
        "median": float(median),
        "sigma_clipped_std": float(std),
        "mad_std": mad,
        "chosen": chosen,
    }


def compute_difference(normalized: NormalizedPair, cfg: dict, units: str = "MJy/sr") -> DifferenceResult:
    warnings_out: list[str] = []
    valid = normalized.valid_mask & np.isfinite(normalized.image1) & np.isfinite(normalized.image2)
    difference = np.subtract(normalized.image2, normalized.image1, dtype=float)
    difference[~valid] = np.nan

    if normalized.variance1 is not None and normalized.variance2 is not None:
        sigma = np.sqrt(np.clip(normalized.variance1, 0, None) + np.clip(normalized.variance2, 0, None))
        sigma[~valid] = np.nan
        sigma_values = sigma[valid]
        global_noise = float(np.nanmedian(sigma_values[np.isfinite(sigma_values)])) if sigma_values.size else float("nan")
        noise_method = "VARIANCE extension propagation (sigma = sqrt(var1 + var2))"
        variance_fraction = float(np.mean(np.isfinite(sigma_values) & (sigma_values > 0))) if sigma_values.size else 0.0
        if not np.isfinite(global_noise) or global_noise <= 0:
            noise_method = "VARIANCE propagation returned no usable uncertainty; fell back to robust difference noise"
            global_noise, noise_detail = global_difference_noise(difference, valid, cfg)
            sigma = np.full(difference.shape, global_noise, dtype=float)
            sigma[~valid] = np.nan
        else:
            noise_detail = {"method": noise_method, "median_sigma": global_noise, "valid_sigma_fraction": variance_fraction}
            if variance_fraction < 0.5:
                warnings_out.append("Fewer than half of valid pixels carry a positive VARIANCE value; SNR is partly extrapolated.")
    else:
        global_noise, noise_detail = global_difference_noise(difference, valid, cfg)
        sigma = np.full(difference.shape, global_noise, dtype=float)
        sigma[~valid] = np.nan
        noise_method = "robust global difference noise (no VARIANCE extension)"
        warnings_out.append("Official VARIANCE extensions were unavailable for at least one epoch; a global noise estimate replaces propagated per-pixel uncertainty.")

    with np.errstate(invalid="ignore", divide="ignore"):
        snr = np.where((sigma > 0) & valid, difference / sigma, np.nan)
    if not np.isfinite(global_noise) or global_noise <= 0:
        warnings_out.append("Difference noise could not be estimated; SNR is undefined for this pair.")

    diagnostics = {
        "equation": "D(x, y) = I2(x, y) - I1(x, y)",
        "sign_convention": "positive residuals indicate the second epoch is brighter than the first at that pixel",
        "noise_model": noise_method,
        "noise_detail": noise_detail,
        "global_noise": global_noise,
        "units": units,
        "difference_pixels_valid": int(np.sum(valid)),
        "difference_pixels_nan": int(np.sum(~valid)),
        "difference_min": float(np.nanmin(difference)) if np.any(valid) else None,
        "difference_max": float(np.nanmax(difference)) if np.any(valid) else None,
        "difference_median": float(np.nanmedian(difference[valid])) if np.any(valid) else None,
        "robust_noise_of_difference": global_difference_noise(difference, valid, cfg)[0],
    }
    warnings_out.append("A difference image can contain residuals from WCS resampling, flux scaling and background mismatch; a residual alone does not establish a change in the sky.")
    return DifferenceResult(
        difference=difference,
        snr=snr,
        sigma_map=sigma,
        valid_mask=valid,
        global_noise=float(global_noise) if np.isfinite(global_noise) else float("nan"),
        diagnostics=diagnostics,
        warnings=warnings_out,
        units=units,
    )
