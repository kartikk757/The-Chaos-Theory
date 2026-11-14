"""Candidate measurement, epoch association and conservative classification.

Classifications are limited to observational categories. The code never labels a
residual as a planet, asteroid or discovery.
"""
from __future__ import annotations

import warnings
from typing import Any

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from photutils.aperture import CircularAnnulus, CircularAperture, aperture_photometry

from modules.models import AlignedPair, Candidate, DifferenceResult, FitsProduct, NormalizedPair


def _column(row, *names, default=float("nan")):
    """Read a photutils output column across versions (x_centroid vs xcentroid)."""
    for name in names:
        try:
            return row[name]
        except (KeyError, ValueError):
            continue
    return default

CLASS_MOVING = "Potential moving-source candidate"
CLASS_VARIABLE = "Potential variable-source candidate"
CLASS_UNMATCHED = "Unmatched residual"
CLASS_ARTIFACT = "Likely image artifact (user review required)"
CLASS_UNCLASSIFIED = "Unclassified detection"

CATEGORY_DEFINITIONS = {
    CLASS_MOVING: "A significant residual with a compact counterpart in both epochs whose measured positions differ by more than the alignment tolerance over the known time baseline. It is a screening category only: asteroids, satellites, cosmic rays and optical ghosts can all mimic it.",
    CLASS_VARIABLE: "A significant residual at a position where both epochs show a detectable source, consistent with a brightness change at fixed sky position. Variability can also be produced by residual background or PSF differences.",
    CLASS_UNMATCHED: "A significant residual with no confident counterpart in one or both epochs. Common causes include cosmic-ray hits, detector artifacts and processing differences.",
    CLASS_ARTIFACT: "The residual overlaps severe FLAGS bits, a masked pixel, or the frame border, so it is more plausibly instrumental than astrophysical.",
    CLASS_UNCLASSIFIED: "The residual passed the SNR threshold but the available measurements do not support any of the defined categories.",
}


def _pixel_solid_angle_sr(pixel_scale_arcsec: float | None) -> float | None:
    if pixel_scale_arcsec is None or not np.isfinite(pixel_scale_arcsec):
        return None
    return float((pixel_scale_arcsec * u.arcsec).to_value(u.rad) ** 2)


def _detect_epoch_sources(image: np.ndarray, valid: np.ndarray, cfg: dict) -> list[dict]:
    values = image[valid & np.isfinite(image)]
    if values.size < 200:
        return []
    from astropy.stats import sigma_clipped_stats
    from photutils.detection import find_peaks

    _, median, std = sigma_clipped_stats(values, sigma=3.0, maxiters=5)
    if not np.isfinite(std) or std <= 0:
        return []
    significance = np.where(valid & np.isfinite(image), (image - median) / std, 0.0)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        table = find_peaks(significance, threshold=5.0, box_size=3, border_width=2)
    if table is None:
        return []
    return [
        {
            "x": float(_column(row, "x_peak", "xcentroid")),
            "y": float(_column(row, "y_peak", "ycentroid")),
            "peak": float(_column(row, "peak_value", "peak")),
            "flux": float(_column(row, "peak_value", "flux")),
            "detection_threshold": float(5.0 * std),
        }
        for row in table
    ]


def _nearest_source(sources: list[dict], x: float, y: float, max_pixels: float) -> dict | None:
    best = None
    best_distance = float("inf")
    for source in sources:
        distance = float(np.hypot(source["x"] - x, source["y"] - y))
        if distance < best_distance:
            best, best_distance = source, distance
    if best is not None and best_distance <= max_pixels:
        return {**best, "distance_px": best_distance}
    return None


def _aperture_flux(image: np.ndarray, x: float, y: float, cfg: dict) -> tuple[float, float]:
    """Aperture sum and annulus-background-subtracted aperture sum in IMAGE units.

    Non-finite pixels are excluded from both the aperture and its effective area.
    """
    aperture = CircularAperture([(x, y)], r=float(cfg["aperture_radius_px"]))
    annulus = CircularAnnulus([(x, y)], r_in=float(cfg["annulus_inner_px"]), r_out=float(cfg["annulus_outer_px"]))
    data = np.nan_to_num(np.asarray(image, dtype=float), nan=0.0)
    mask = ~np.isfinite(np.asarray(image, dtype=float))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        phot = aperture_photometry(data, [aperture, annulus], mask=mask, method="exact")
        annulus_area = float(annulus.area_overlap(data, mask=mask, method="exact"))
        aperture_area = float(aperture.area_overlap(data, mask=mask, method="exact"))
    raw_sum = float(phot["aperture_sum_0"][0])
    annulus_sum = float(phot["aperture_sum_1"][0])
    local_background = annulus_sum / annulus_area if annulus_area > 0 else 0.0
    net = raw_sum - local_background * aperture_area
    return raw_sum, float(net)


def _sky_offset_arcsec(aligned: AlignedPair, x1: float, y1: float, x2: float, y2: float) -> tuple[float | None, float | None]:
    try:
        coord1 = aligned.target_wcs.pixel_to_world(x1, y1)
        coord2 = aligned.target_wcs.pixel_to_world(x2, y2)
        separation = float(coord1.separation(coord2).arcsec)
        position_angle = float(coord1.position_angle(coord2).to_value(u.deg))
        return separation, position_angle
    except Exception:
        pixel_scale = aligned.pixel_scale_arcsec
        if pixel_scale and np.isfinite(pixel_scale):
            return float(np.hypot(x2 - x1, y2 - y1) * pixel_scale), None
        return None, None


def analyse_candidates(
    aligned: AlignedPair,
    normalized: NormalizedPair,
    difference: DifferenceResult,
    detection_sources: list[dict],
    product1: FitsProduct,
    product2: FitsProduct,
    record1: dict,
    record2: dict,
    cfg: dict,
    provenance_run_id: str,
    data_mode: str,
) -> tuple[list[Candidate], dict]:
    warnings_out: list[str] = []
    sources1 = _detect_epoch_sources(aligned.science1, difference.valid_mask, cfg)
    sources2 = _detect_epoch_sources(aligned.science2, difference.valid_mask, cfg)
    solid_angle = _pixel_solid_angle_sr(aligned.pixel_scale_arcsec)
    delta_days = None
    if record1.get("mjd") is not None and record2.get("mjd") is not None:
        delta_days = abs(float(record2["mjd"]) - float(record1["mjd"]))
    search_radius_px = float(cfg["annulus_outer_px"]) + 2.0
    candidates: list[Candidate] = []

    for index, source in enumerate(detection_sources, start=1):
        x, y = float(source["xcentroid"]), float(source["ycentroid"])
        xi, yi = int(round(x)), int(round(y))
        window_pixels = int(source.get("pixels_above_threshold_in_window", 0))
        sign = int(source.get("sign", 1))
        dipole = source.get("dipole_partner")
        dipole_separation_px = source.get("dipole_separation_px")
        local_noise = float(source.get("local_noise") or difference.global_noise)
        reasons: list[str] = []
        candidate_warnings: list[str] = []
        quality: dict[str, Any] = {}

        ra = dec = None
        try:
            world = aligned.target_wcs.pixel_to_world(x, y)
            ra, dec = float(world.ra.deg), float(world.dec.deg)
        except Exception:
            candidate_warnings.append("Pixel-to-sky conversion failed for this position.")

        flag_window1 = product1.bad_flag_mask[max(yi - 1, 0):yi + 2, max(xi - 1, 0):xi + 2]
        flag_window2 = product2.bad_flag_mask[max(yi - 1, 0):yi + 2, max(xi - 1, 0):xi + 2]
        flagged1 = float(np.mean(flag_window1)) if flag_window1.size else 0.0
        flagged2 = float(np.mean(flag_window2)) if flag_window2.size else 0.0

        raw1, net1 = _aperture_flux(aligned.science1, x, y, cfg)
        raw2, net2 = _aperture_flux(aligned.science2, x, y, cfg)
        # 1 MJy/sr = 1e6 Jy/sr, so F[Jy] = SB[MJy/sr] × Ω[pixel, sr] × 1e6.
        to_jy = (solid_angle * 1e6) if solid_angle else None
        flux1_jy = net1 * to_jy if to_jy is not None else float("nan")
        flux2_jy = net2 * to_jy if to_jy is not None else float("nan")
        delta_flux = flux2_jy - flux1_jy
        if to_jy is not None:
            quality["pixel_solid_angle_sr"] = solid_angle
            quality["flux_units"] = "Jy (aperture sum of MJy/sr surface brightness × pixel solid angle)"
        else:
            candidate_warnings.append("Pixel scale unavailable; aperture brightness is reported in surface-brightness units only.")

        match1 = _nearest_source(sources1, x, y, search_radius_px)
        match2 = _nearest_source(sources2, x, y, search_radius_px)
        displacement_arcsec = None
        position_angle_deg = None
        if match1 and match2:
            displacement_arcsec, position_angle_deg = _sky_offset_arcsec(aligned, match1["x"], match1["y"], match2["x"], match2["y"])
        alignment_tolerance_arcsec = 1.5 * (aligned.pixel_scale_arcsec or 0.0)

        classification = CLASS_UNCLASSIFIED
        if flagged1 >= 0.34 or flagged2 >= 0.34:
            classification = CLASS_ARTIFACT
            reasons.append(f"{flagged1:.0%} / {flagged2:.0%} of the 3×3 aperture neighbourhood carries severe FLAGS bits in epoch 1 / 2.")
        elif window_pixels <= 1 and not (match1 and match2) and source["peak_snr"] >= 20:
            classification = CLASS_ARTIFACT
            reasons.append(
                "The residual is a single-pixel spike an order of magnitude above the sampled point-spread function and has no "
                "compact counterpart in either epoch; cosmic-ray-like hits and resampling residues produce this signature. "
                "This is a morphological heuristic, not a definitive diagnosis."
            )
        elif not np.isfinite(source["peak_snr"]):
            classification = CLASS_UNCLASSIFIED
            reasons.append("Peak SNR at this position is not finite.")
        elif displacement_arcsec is not None and displacement_arcsec > alignment_tolerance_arcsec and (match1 or match2):
            classification = CLASS_MOVING
            reasons.append(
                f"A compact source is detected near this position in both epochs, separated by {displacement_arcsec:.2f} arcsec "
                f"({displacement_arcsec / (aligned.pixel_scale_arcsec or 1):.1f} px) over Δt = {delta_days:.2f} d."
            )
        elif match1 and match2:
            ratio = (abs(delta_flux) / max(abs(flux1_jy), 1e-30)) if np.isfinite(flux1_jy) else float("nan")
            if np.isfinite(ratio) and ratio > 0.25:
                classification = CLASS_VARIABLE
                reasons.append(
                    f"The same sky position has detectable flux in both epochs, with a measured change of {delta_flux * 1e6:.3f} µJy "
                    f"({ratio * 100:.1f}% of the epoch-1 aperture flux)."
                )
            else:
                classification = CLASS_UNCLASSIFIED
                reasons.append("Both epochs show detections here but the aperture flux difference is small relative to the epoch-1 flux.")
        elif match1 or match2:
            classification = CLASS_UNMATCHED
            reasons.append("A counterpart is present in only one epoch within the search radius.")
        else:
            classification = CLASS_UNMATCHED
            reasons.append("No compact counterpart is detected in either epoch within the association radius.")

        reasons.insert(
            0,
            f"Residual peak reaches {source['peak_snr']:.1f}σ against a local difference noise of {local_noise:.4g} "
            f"(detection threshold {cfg['detection_sigma']:.1f}σ) and spans {window_pixels} pixel(s) above threshold in a 3×3 window.",
        )
        reasons.insert(
            1,
            "Sign of the residual: the second epoch is brighter at this pixel."
            if sign > 0 else
            "Sign of the residual: the first epoch is brighter at this pixel.",
        )
        if dipole is not None and dipole_separation_px is not None:
            dipole_arcsec = dipole_separation_px * aligned.pixel_scale_arcsec if aligned.pixel_scale_arcsec else None
            reasons.append(
                f"The residual is one lobe of an opposite-sign dipole separated by {dipole_separation_px:.1f} px"
                + (f" ({dipole_arcsec:.1f} arcsec)" if dipole_arcsec else "")
                + ". A dipole is the signature of a positional shift or a registration difference, not of a pure brightness change."
            )
        if window_pixels <= 1:
            candidate_warnings.append(
                "The residual is confined to a single pixel. With an undersampled SPHEREx PSF this can be a real source, "
                "a cosmic-ray-like defect, or a resampling residue, so it cannot be classified as astrophysical on its own."
            )
        elif window_pixels >= 8:
            candidate_warnings.append("The residual covers most of a 3×3 window, suggesting an extended or poorly subtracted feature rather than a point source.")
        if displacement_arcsec is not None and alignment_tolerance_arcsec and displacement_arcsec <= alignment_tolerance_arcsec:
            candidate_warnings.append("Any apparent offset is within the WCS/alignment tolerance; motion is not established.")
        if displacement_arcsec is None:
            reasons.append("Angular displacement was not measurable: a matched compact counterpart was not available in both epochs.")

        quality.update({
            "peak_snr": float(source["peak_snr"]),
            "sign": sign,
            "dipole_partner_pixel": [dipole["xcentroid"], dipole["ycentroid"]] if dipole else None,
            "dipole_separation_px": dipole_separation_px,
            "pixels_above_threshold_in_window": window_pixels,
            "local_noise": local_noise,
            "severe_flag_fraction_epoch1": flagged1,
            "severe_flag_fraction_epoch2": flagged2,
            "alignment_tolerance_arcsec": alignment_tolerance_arcsec,
            "matched_epoch1": match1 is not None,
            "matched_epoch2": match2 is not None,
            "offset_epoch1_px": match1["distance_px"] if match1 else None,
            "offset_epoch2_px": match2["distance_px"] if match2 else None,
            "epoch_sources_found": [len(sources1), len(sources2)],
            "delta_days": delta_days,
            "aperture_radius_px": cfg["aperture_radius_px"],
        })
        candidate = Candidate(
            candidate_id=f"CAND-{index:03d}",
            pixel_x=x,
            pixel_y=y,
            ra_deg=ra if ra is not None else float("nan"),
            dec_deg=dec if dec is not None else float("nan"),
            classification=classification,
            snr=float(source["peak_snr"]),
            flux1_jy=float(flux1_jy),
            flux2_jy=float(flux2_jy),
            delta_flux_jy=float(delta_flux),
            surface_brightness_diff=float(np.nan_to_num(difference.difference[yi, xi], nan=float("nan"))),
            local_noise=local_noise,
            quality=quality,
            reasons=reasons,
            warnings=candidate_warnings,
            data_mode=data_mode,
            provenance_run_id=provenance_run_id,
            metadata={
                "epoch1": {
                    "id": record1.get("id"),
                    "date_utc": record1.get("date_utc"),
                    "mjd": record1.get("mjd"),
                    "band": record1.get("band"),
                    "origin": record1.get("origin"),
                    "release": record1.get("release"),
                },
                "epoch2": {
                    "id": record2.get("id"),
                    "date_utc": record2.get("date_utc"),
                    "mjd": record2.get("mjd"),
                    "band": record2.get("band"),
                    "origin": record2.get("origin"),
                    "release": record2.get("release"),
                },
                "residual_sign": sign,
                "dipole_separation_px": dipole_separation_px,
                "angular_displacement_arcsec": displacement_arcsec,
                "position_angle_deg": position_angle_deg,
                "delta_time_days": delta_days,
                "units": {
                    "flux": "Jy",
                    "surface_brightness": "MJy/sr",
                    "angular_displacement": "arcsec",
                    "position_angle": "degrees east of north",
                },
            },
        )
        candidates.append(candidate)

    candidates.sort(key=lambda item: item.snr, reverse=True)
    for rank, candidate in enumerate(candidates, start=1):
        candidate.candidate_id = f"CAND-{rank:03d}"
    summary = {
        "counts": {name: sum(1 for c in candidates if c.classification == name) for name in CATEGORY_DEFINITIONS},
        "total": len(candidates),
        "epoch_sources": [len(sources1), len(sources2)],
        "association_radius_px": search_radius_px,
        "classification_definitions": CATEGORY_DEFINITIONS,
        "note": "Classifications are screening labels derived from measurements on two real observations; they are not discovery claims.",
    }
    if not candidates:
        warnings_out.append("No residual exceeded the configured detection threshold; this is a valid null result for this pair.")
    return candidates, summary
