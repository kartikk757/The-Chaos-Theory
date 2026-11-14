"""Configurable, controlled source-injection/recovery evaluation utilities.

This evaluates the multi-epoch peak detector on a caller-supplied control
image and kernel. It does not calibrate the science pipeline or define
probabilities. A blank-control assumption is recorded explicitly.
"""
from __future__ import annotations

import hashlib
import io
from typing import Any

import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from scipy.ndimage import maximum_filter
from scipy.optimize import linear_sum_assignment


def load_2d_fits_image(payload: bytes, max_pixels: int = 4_000_000) -> tuple[np.ndarray, dict]:
    raw = bytes(payload)
    digest = hashlib.sha256(raw).hexdigest()
    with fits.open(io.BytesIO(raw), memmap=False) as hdul:
        hdu = next((item for item in hdul if str(item.name).upper() == "IMAGE" and item.data is not None), None)
        if hdu is None:
            hdu = next((item for item in hdul if item.data is not None and getattr(item.data, "ndim", 0) == 2), None)
        if hdu is None or hdu.data is None:
            raise ValueError("FITS file does not contain a two-dimensional image.")
        image = np.asarray(hdu.data, dtype=float)
        if image.ndim != 2 or int(image.size) > max_pixels:
            raise ValueError(f"Expected a 2-D image with no more than {max_pixels:,} pixels.")
        return image.copy(), {"filename": "user-upload.fits", "extension": hdu.name, "shape": list(image.shape), "bunit": hdu.header.get("BUNIT", "not supplied"), "sha256": digest}


def _robust_noise(image: np.ndarray) -> tuple[float, float]:
    finite = np.asarray(image, dtype=float)[np.isfinite(image)]
    if finite.size < 100:
        raise ValueError("Control image has fewer than 100 finite pixels.")
    mean, median, std = sigma_clipped_stats(finite, sigma=3.0, maxiters=5)
    if not np.isfinite(std) or std <= 0:
        raise ValueError("Control image noise could not be estimated; supply a non-zero noise scale.")
    return float(median), float(std)


def _detect_peaks(image: np.ndarray, background: float, noise: float, threshold_sigma: float, border: int, min_distance: float) -> list[dict]:
    array = np.asarray(image, dtype=float)
    finite = np.isfinite(array)
    significance = np.where(finite, (array - background) / noise, -np.inf)
    maxima = finite & (significance >= threshold_sigma) & (significance == maximum_filter(significance, size=3, mode="constant", cval=-np.inf))
    if border > 0:
        maxima[:border] = False
        maxima[-border:] = False
        maxima[:, :border] = False
        maxima[:, -border:] = False
    ys, xs = np.where(maxima)
    order = np.argsort(significance[ys, xs])[::-1]
    selected = []
    for idx in order:
        x, y = int(xs[idx]), int(ys[idx])
        if any(np.hypot(x - other["x_pixel"], y - other["y_pixel"]) < min_distance for other in selected):
            continue
        selected.append({"x_pixel": x, "y_pixel": y, "snr": float(significance[y, x])})
    return selected


def _match_detections(detections: list[dict], expected: list[dict], radius_px: float) -> tuple[list[dict], set[int], set[int]]:
    if not detections or not expected:
        return [], set(), set()
    cost = np.full((len(expected), len(detections)), 1e9, dtype=float)
    for row, source in enumerate(expected):
        for column, detection in enumerate(detections):
            distance = float(np.hypot(source["x_pixel"] - detection["x_pixel"], source["y_pixel"] - detection["y_pixel"]))
            if distance <= radius_px:
                cost[row, column] = distance
    rows, cols = linear_sum_assignment(cost)
    matches, expected_used, detections_used = [], set(), set()
    for row, column in zip(rows, cols):
        if cost[row, column] >= 1e9:
            continue
        expected_used.add(int(row))
        detections_used.add(int(column))
        matches.append({"injection_index": row, "detection_index": column, "separation_px": float(cost[row, column]), "recovered_snr": detections[column]["snr"]})
    return matches, expected_used, detections_used


def run_injection_recovery(
    control_image: np.ndarray,
    empirical_kernel: np.ndarray,
    *,
    amplitudes_sigma: tuple[float, ...] | list[float] = (4.0, 5.0, 7.0),
    injections_per_level: int = 12,
    threshold_sigma: float = 5.0,
    match_radius_px: float = 2.0,
    min_peak_distance_px: float = 3.0,
    border_px: int = 5,
    random_seed: int = 20261114,
    noise_sigma: float | None = None,
    control_assumed_source_free: bool = False,
) -> dict[str, Any]:
    """Inject a user-supplied empirical kernel and measure detector recovery."""
    image = np.asarray(control_image, dtype=float)
    kernel = np.asarray(empirical_kernel, dtype=float)
    if image.ndim != 2 or min(image.shape) < 16 or not np.any(np.isfinite(image)):
        raise ValueError("Control image must be a finite two-dimensional array at least 16 pixels on each side.")
    if kernel.ndim != 2 or min(kernel.shape) < 3 or any(size % 2 == 0 for size in kernel.shape):
        raise ValueError("Empirical source kernel must be a two-dimensional odd-sized image at least 3×3 pixels.")
    if not np.all(np.isfinite(kernel)) or np.any(kernel < 0) or float(np.sum(kernel)) <= 0:
        raise ValueError("Kernel must contain finite, non-negative values with positive sum.")
    kernel = kernel / float(np.sum(kernel))
    if not amplitudes_sigma or any(not np.isfinite(value) or value <= 0 for value in amplitudes_sigma):
        raise ValueError("Injection amplitudes must be positive finite values expressed in estimated-noise sigma units.")
    if not 1 <= int(injections_per_level) <= 1000:
        raise ValueError("Injection count per amplitude must be between 1 and 1000.")
    if not np.isfinite(threshold_sigma) or threshold_sigma <= 0 or match_radius_px <= 0:
        raise ValueError("Detection threshold and match radius must be positive.")
    background, estimated_noise = _robust_noise(image)
    noise = float(noise_sigma) if noise_sigma is not None else estimated_noise
    if not np.isfinite(noise) or noise <= 0:
        raise ValueError("Noise scale must be positive and finite.")
    py, px = kernel.shape[0] // 2, kernel.shape[1] // 2
    required_border = max(int(border_px), py + 2, px + 2)
    if image.shape[0] <= 2 * required_border or image.shape[1] <= 2 * required_border:
        raise ValueError("Control image is too small for the selected kernel and border exclusion.")
    rng = np.random.default_rng(int(random_seed))
    candidate_locations = []
    attempts = 0
    while len(candidate_locations) < int(injections_per_level) and attempts < 5000:
        attempts += 1
        y = int(rng.integers(required_border, image.shape[0] - required_border))
        x = int(rng.integers(required_border, image.shape[1] - required_border))
        if all(np.hypot(x - other["x_pixel"], y - other["y_pixel"]) >= max(kernel.shape) * 2 for other in candidate_locations):
            candidate_locations.append({"x_pixel": x, "y_pixel": y})
    if len(candidate_locations) < int(injections_per_level):
        raise ValueError("Could not place the requested non-overlapping injections; reduce injection count or use a larger control image.")
    image_hash = hashlib.sha256(np.ascontiguousarray(image).tobytes()).hexdigest()
    kernel_hash = hashlib.sha256(np.ascontiguousarray(kernel).tobytes()).hexdigest()
    baseline_detections = _detect_peaks(image, background, noise, threshold_sigma, required_border, min_peak_distance_px)
    level_reports = []
    for amplitude in amplitudes_sigma:
        injected = image.copy()
        expected = []
        for location in candidate_locations:
            x, y = location["x_pixel"], location["y_pixel"]
            amp = float(amplitude) * noise / float(np.max(kernel))
            injected[y - py:y + py + 1, x - px:x + px + 1] += kernel * amp
            expected.append({**location, "requested_peak_snr": float(amplitude)})
        detections = _detect_peaks(injected, background, noise, threshold_sigma, required_border, min_peak_distance_px)
        matches, expected_used, detections_used = _match_detections(detections, expected, float(match_radius_px))
        extras = [row for index, row in enumerate(detections) if index not in detections_used]
        level_reports.append({
            "injected_peak_snr": float(amplitude),
            "injected_sources": len(expected),
            "recovered_sources": len(matches),
            "missed_sources": len(expected) - len(matches),
            "empirical_recovery_fraction": len(matches) / len(expected),
            "matches": matches,
            "unassociated_peaks_in_injected_image": extras,
            "detections_total": len(detections),
        })
    return {
        "status": "ok",
        "evaluation": "controlled empirical-kernel injection/recovery using a local-maximum detector",
        "configuration": {
            "amplitudes_sigma": [float(value) for value in amplitudes_sigma],
            "injections_per_level": int(injections_per_level),
            "detection_threshold_sigma": float(threshold_sigma),
            "match_radius_px": float(match_radius_px),
            "minimum_peak_distance_px": float(min_peak_distance_px),
            "border_pixels_excluded": required_border,
            "random_seed": int(random_seed),
            "noise_sigma_used": noise,
            "noise_sigma_robust_estimate": estimated_noise,
            "background_median": background,
        },
        "control_image": {"shape": list(image.shape), "sha256": image_hash, "assumed_source_free_by_user": bool(control_assumed_source_free)},
        "empirical_kernel": {"shape": list(kernel.shape), "sha256": kernel_hash, "normalization": "kernel sum normalized to 1; requested peak amplitude scaled in noise-sigma units"},
        "no_injection_control": {
            "detected_peak_count": len(baseline_detections),
            "peaks": baseline_detections,
            "interpretation": "Count can be called false positives only if the caller's source-free control assumption is valid; software cannot verify that assumption.",
        },
        "levels": level_reports,
        "scientific_limitations": [
            "Recovery fractions apply only to this control image, empirical kernel, noise estimate, injection placement and detector configuration.",
            "They are not a calibrated detection completeness curve, false-positive probability, or sensitivity improvement claim.",
            "A supplied kernel is not assumed to be a validated SPHEREx ePSF. Oversampled ePSFs must not be injected without validated detector sampling and subpixel shifting.",
            "Unassociated peaks in a real image may be astrophysical sources rather than false positives.",
        ],
    }
