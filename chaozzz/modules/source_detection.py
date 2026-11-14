"""Statistically significant residual-source detection.

Method note: the SPHEREx PSF is strongly undersampled in full-frame spectral
images (FWHM is comparable to the ~6.1 arcsec pixel), so a DAO-style point-spread
model that requires several pixels is not appropriate. Detection is therefore
performed with a box-car peak finder (Photutils ``find_peaks``) applied to the
SNR map, with bright neighbours above the threshold suppressed by the box size.
"""
from __future__ import annotations

import warnings
from typing import Any

import numpy as np
from astropy.stats import sigma_clipped_stats
from photutils.detection import find_peaks

from modules.models import AlignedPair, DifferenceResult


def _column(row, *names, default=float("nan")):
    """Read a Photutils output column across versions (x_peak vs xcentroid)."""
    for name in names:
        try:
            return row[name]
        except (KeyError, ValueError):
            continue
    return default


def estimate_noise_map(difference: DifferenceResult, aligned: AlignedPair, block: int = 16) -> tuple[np.ndarray, dict]:
    """Block-wise robust noise, used to describe local quality (not to rescale SNR)."""
    image = difference.difference
    valid = difference.valid_mask
    rows, columns = image.shape
    noise = np.full(image.shape, np.nan, dtype=float)
    for row in range(0, rows, block):
        for column in range(0, columns, block):
            tile = image[row:row + block, column:column + block]
            mask = valid[row:row + block, column:column + block] & np.isfinite(tile)
            values = tile[mask]
            if values.size >= 40:
                _, _, std = sigma_clipped_stats(values, sigma=3.0, maxiters=3)
                if np.isfinite(std) and std > 0:
                    noise[row:row + block, column:column + block] = std
    finite = np.isfinite(noise)
    fill = float(np.nanmedian(noise[finite])) if finite.any() else float(difference.global_noise)
    noise[~finite] = fill
    return noise, {
        "method": f"sigma-clipped standard deviation in {block}x{block} tiles",
        "tiles": int(np.ceil(rows / block) * np.ceil(columns / block)),
        "median_local_noise": fill,
    }


def _peaks(array: np.ndarray, threshold: float, border: int, sign: int) -> list[dict]:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table = find_peaks(array, threshold=threshold, box_size=3, border_width=border)
    if table is None:
        return []
    output = []
    for row in table:
        output.append({
            "xcentroid": float(_column(row, "x_peak", "xcentroid")),
            "ycentroid": float(_column(row, "y_peak", "ycentroid")),
            "peak_snr": float(_column(row, "peak_value", "peak")),
            "sign": sign,
        })
    return output


def detect_residual_sources(
    aligned: AlignedPair,
    normalized,
    difference: DifferenceResult,
    cfg: dict,
) -> tuple[list[dict[str, Any]], dict]:
    """Detect significant positive and negative residuals and merge ± dipoles.

    A real brightness change at a fixed position produces a single-signed residual.
    A positional change (or a mis-registration) produces a positive/negative dipole.
    Both are detected; a neighbouring pair of opposite sign is reported as one
    dipole candidate so a single object is not counted twice.
    """
    valid = difference.valid_mask & np.isfinite(difference.difference) & np.isfinite(difference.snr)
    if not np.isfinite(difference.global_noise) or difference.global_noise <= 0:
        return [], {"status": "skipped", "reason": "difference noise is undefined", "sources": 0}
    snr = np.where(valid, difference.snr, np.nan)
    work = np.nan_to_num(snr, nan=0.0, posinf=0.0, neginf=0.0)
    local_noise, noise_diag = estimate_noise_map(difference, aligned)
    border = int(cfg["border_px"])
    threshold = float(cfg["detection_sigma"])
    positives = sorted(_peaks(work, threshold, border, +1), key=lambda item: item["peak_snr"], reverse=True)
    negatives = sorted(_peaks(-work, threshold, border, -1), key=lambda item: item["peak_snr"], reverse=True)
    dipole_distance_px = float(cfg.get("dipole_distance_px", 4.0))
    used_negative: set[int] = set()
    merged: list[dict] = []

    def window_count(x: float, y: float) -> int:
        xi = int(np.clip(round(x), 0, snr.shape[1] - 1))
        yi = int(np.clip(round(y), 0, snr.shape[0] - 1))
        window = np.abs(snr[max(yi - 1, 0):yi + 2, max(xi - 1, 0):xi + 2])
        return int(np.sum(np.isfinite(window) & (window >= threshold)))

    for peak in positives:
        partner = None
        for index, negative in enumerate(negatives):
            if index in used_negative:
                continue
            distance = float(np.hypot(negative["xcentroid"] - peak["xcentroid"], negative["ycentroid"] - peak["ycentroid"]))
            if distance <= dipole_distance_px:
                partner = (index, negative, distance)
                break
        if partner is not None:
            used_negative.add(partner[0])
        merged.append({**peak, "dipole_partner": partner[1] if partner else None, "dipole_separation_px": partner[2] if partner else None})
    for index, negative in enumerate(negatives):
        if index in used_negative:
            continue
        merged.append({**negative, "dipole_partner": None, "dipole_separation_px": None})
    merged.sort(key=lambda item: item["peak_snr"], reverse=True)

    sources: list[dict[str, Any]] = []
    for candidate in merged:
        x, y = candidate["xcentroid"], candidate["ycentroid"]
        xi = int(np.clip(round(x), 0, snr.shape[1] - 1))
        yi = int(np.clip(round(y), 0, snr.shape[0] - 1))
        if not valid[yi, xi]:
            continue
        sources.append({
            "xcentroid": x,
            "ycentroid": y,
            "peak_snr": float(candidate["peak_snr"]),
            "sign": int(candidate["sign"]),
            "local_noise": float(local_noise[yi, xi]),
            "pixels_above_threshold_in_window": window_count(x, y),
            "dipole_partner": candidate["dipole_partner"],
            "dipole_separation_px": candidate["dipole_separation_px"],
        })
        if len(sources) >= int(cfg["max_candidates"]):
            break
    return sources, {
        "algorithm": "Photutils find_peaks on the SNR map and on its inverse; ± pairs merged as dipoles",
        "method_note": "SPHEREx full-frame pixels (~6.1 arcsec) undersample the PSF, so peak detection is used instead of a multi-pixel DAO model.",
        "threshold_sigma": threshold,
        "dipole_merge_radius_px": dipole_distance_px,
        "positive_peaks": len(positives),
        "negative_peaks": len(negatives),
        "merged_candidates": len(sources),
        "noise": noise_diag,
        "max_candidates": cfg["max_candidates"],
        "status": "ok",
        "sources": len(sources),
        "note": "A detection is a statistically significant residual pixel group, not an astrophysical identification.",
    }
