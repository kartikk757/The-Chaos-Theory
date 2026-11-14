"""Background estimation, robust normalization and noise-model construction."""
from __future__ import annotations

import numpy as np
from astropy.stats import sigma_clipped_stats

from modules.models import AlignedPair, NormalizedPair


def robust_stats(values: np.ndarray, sigma: float = 3.0, max_iterations: int = 5) -> dict:
    array = np.asarray(values, dtype=float)
    array = array[np.isfinite(array)]
    if array.size == 0:
        return {"median": float("nan"), "mean": float("nan"), "std": float("nan"), "count": 0}
    mean, median, std = sigma_clipped_stats(array, sigma=sigma, maxiters=max_iterations)
    return {"median": float(median), "mean": float(mean), "std": float(std), "count": int(array.size)}


def fit_background_plane(
    image: np.ndarray,
    valid_mask: np.ndarray,
    excluded_mask: np.ndarray | None = None,
    sigma: float = 3.0,
    iterations: int = 5,
) -> dict:
    """Fit a robust tilted plane to unmasked sky pixels (instrumental background drift)."""
    rows, columns = np.indices(image.shape)
    usable = valid_mask & np.isfinite(image)
    if excluded_mask is not None:
        usable &= ~excluded_mask
    values = image[usable]
    if values.size < 50:
        return {"plane": None, "offset": 0.0, "gradient_x": 0.0, "gradient_y": 0.0, "pixels_used": int(values.size), "note": "too few usable pixels"}
    keep = np.ones(values.shape, dtype=bool)
    coefficients = None
    design_full = np.column_stack([np.ones(values.size), columns[usable], rows[usable]])
    for _ in range(max(1, iterations)):
        if keep.sum() < 50:
            break
        coefficients, *_ = np.linalg.lstsq(design_full[keep], values[keep], rcond=None)
        residual = values - design_full @ coefficients
        _, median, std = sigma_clipped_stats(residual[keep], sigma=sigma, maxiters=3)
        if not np.isfinite(std) or std <= 0:
            break
        keep = np.abs(residual - median) < sigma * std
    if coefficients is None:
        return {"plane": None, "offset": 0.0, "gradient_x": 0.0, "gradient_y": 0.0, "pixels_used": 0, "note": "fit failed"}
    plane = coefficients[0] + coefficients[1] * columns + coefficients[2] * rows
    return {
        "plane": plane,
        "offset": float(coefficients[0]),
        "gradient_x": float(coefficients[1]),
        "gradient_y": float(coefficients[2]),
        "pixels_used": int(keep.sum()),
        "note": "robust least-squares tilted plane",
    }


def normalize_pair(aligned: AlignedPair, cfg: dict) -> NormalizedPair:
    warnings_out: list[str] = []
    valid = aligned.valid_mask.copy()
    image1 = np.asarray(aligned.science1, dtype=float).copy()
    image2 = np.asarray(aligned.science2, dtype=float).copy()
    source_free = valid & ~aligned.source_mask1 & ~aligned.source_mask2
    # The SPHEREx SOURCE flag is informational and marks most pixels in a real cutout,
    # so it is only used as an exclusion mask when it is actually informative.
    if source_free.sum() >= max(200, int(0.3 * valid.sum())):
        sky_pixels = source_free
    else:
        sky_pixels = valid
        warnings_out.append(
            "The informational SOURCE flag marks most pixels in this cutout, so sigma-clipping (not flag exclusion) "
            "is used to reject sources from the background fit."
        )

    stats1 = robust_stats(image1[sky_pixels], cfg["background_sigma"], cfg["background_iterations"])
    stats2 = robust_stats(image2[sky_pixels], cfg["background_sigma"], cfg["background_iterations"])
    background1 = {"method": "sigma-clipped constant", **stats1}
    background2 = {"method": "sigma-clipped constant", **stats2}

    if cfg.get("fit_background_plane"):
        exclusion1 = aligned.source_mask1 if sky_pixels is not valid else None
        exclusion2 = aligned.source_mask2 if sky_pixels is not valid else None
        plane1 = fit_background_plane(image1, valid, exclusion1)
        plane2 = fit_background_plane(image2, valid, exclusion2)
        if plane1["plane"] is not None and plane2["plane"] is not None:
            image1 = image1 - plane1["plane"]
            image2 = image2 - plane2["plane"]
            background1 = {
                "method": "robust tilted plane",
                "offset": plane1["offset"],
                "gradient_x_per_px": plane1["gradient_x"],
                "gradient_y_per_px": plane1["gradient_y"],
                "pixels_used": plane1["pixels_used"],
            }
            background2 = {
                "method": "robust tilted plane",
                "offset": plane2["offset"],
                "gradient_x_per_px": plane2["gradient_x"],
                "gradient_y_per_px": plane2["gradient_y"],
                "pixels_used": plane2["pixels_used"],
            }
        else:
            warnings_out.append("Background plane fit was not possible; fell back to sigma-clipped constants.")

    post1 = robust_stats(image1[sky_pixels], cfg["background_sigma"], cfg["background_iterations"])
    post2 = robust_stats(image2[sky_pixels], cfg["background_sigma"], cfg["background_iterations"])
    if not np.isfinite(post1["median"]) or not np.isfinite(post2["median"]):
        warnings_out.append("Background statistics are not finite; normalization scaling was not applied.")
        scale = 1.0
    else:
        image1 = image1 - post1["median"]
        image2 = image2 - post2["median"]
        scale = 1.0

    residual1 = robust_stats(image1[sky_pixels], cfg["background_sigma"], cfg["background_iterations"])
    residual2 = robust_stats(image2[sky_pixels], cfg["background_sigma"], cfg["background_iterations"])
    noise_reference = max(residual1["std"], residual2["std"], 0.0)
    linear_coefficients = None
    gain = 1.0
    if noise_reference > 0 and np.isfinite(noise_reference):
        samples = int(sky_pixels.sum())
        rng = np.random.default_rng(20260517)
        indices = rng.choice(samples, size=min(cfg["robust_fit_samples"], samples), replace=False)
        x = image1[sky_pixels][indices]
        y = image2[sky_pixels][indices]
        good = np.isfinite(x) & np.isfinite(y)
        if good.sum() > 200:
            design = np.column_stack([x[good], np.ones(good.sum())])
            coefficients, *_ = np.linalg.lstsq(design, y[good], rcond=None)
            gain, offset = float(coefficients[0]), float(coefficients[1])
            if 0.2 <= gain <= 5.0 and np.isfinite(gain):
                image2 = (image2 - offset) / gain
                linear_coefficients = {"gain": gain, "offset": offset, "samples": int(good.sum())}
            else:
                linear_coefficients = {"gain": gain, "offset": offset, "samples": int(good.sum()), "applied": False}
                warnings_out.append(f"Robust linear gain {gain:.3f} is implausible for identical-band data and was not applied.")
        else:
            warnings_out.append("Too few sky samples for a robust linear scaling coefficient.")
    final1 = robust_stats(image1[sky_pixels], cfg["background_sigma"], cfg["background_iterations"])
    final2 = robust_stats(image2[sky_pixels], cfg["background_sigma"], cfg["background_iterations"])

    diagnostics = {
        "normalization_mode": "background plane + sigma-clipped constant + robust linear sky scaling" if linear_coefficients else "background-only",
        "sky_pixels_used": int(sky_pixels.sum()),
        "background_epoch1": background1,
        "background_epoch2": background2,
        "residual_noise_epoch1": residual1["std"],
        "residual_noise_epoch2": residual2["std"],
        "post_normalization_noise_epoch1": final1["std"],
        "post_normalization_noise_epoch2": final2["std"],
        "linear_scaling": linear_coefficients,
        "global_scale_factor": scale,
        "units_after_normalization": "MJy/sr relative to each epoch's own sky background",
        "note": "Identical-band epochs are expected to share the same photometric scale; only sky background and a robust linear term are removed.",
    }
    variance1 = aligned.variance1
    variance2 = aligned.variance2
    if linear_coefficients and linear_coefficients.get("gain") and variance2 is not None:
        variance2 = variance2 / (linear_coefficients["gain"] ** 2)
    return NormalizedPair(
        image1=image1,
        image2=image2,
        variance1=variance1,
        variance2=variance2,
        valid_mask=valid,
        diagnostics=diagnostics,
        warnings=warnings_out,
    )
