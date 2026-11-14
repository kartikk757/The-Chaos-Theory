"""Multi-epoch source association and motion screening.

Tracks are observational leads only. Motion uncertainties use a stated
conservative per-epoch astrometric floor and are not a calibrated SPHEREx error
model. A two-position displacement is always marked ``2-EPOCH-ONLY``.
"""
from __future__ import annotations

import itertools
import math
from typing import Any

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord, SkyOffsetFrame
from astropy.stats import sigma_clipped_stats
from astropy.time import Time
from scipy.ndimage import maximum_filter
from scipy.optimize import linear_sum_assignment

import config
from modules import coordinates, provenance, wcs_alignment

DEFAULT_MULTI_EPOCH_CONFIG = {
    "detection_sigma": 5.0,
    "border_px": 2,
    "max_candidates": 80,
    "multi_epoch_association_radius_arcsec": 12.0,
    "multi_epoch_max_motion_arcsec_per_day": 30.0,
    "multi_epoch_min_epochs_for_motion": 3,
    "multi_epoch_min_motion_snr": 3.0,
    "multi_epoch_max_reduced_chi2": 25.0,
    "min_overlap_fraction": 0.55,
}

MAX_MULTI_EPOCHS = 20


def _time_mjd(record: dict) -> float | None:
    value = record.get("mjd")
    try:
        if value is not None and np.isfinite(float(value)):
            return float(value)
    except (TypeError, ValueError):
        pass
    date_value = record.get("date_utc")
    if not date_value:
        return None
    try:
        return float(Time(str(date_value), format="isot", scale="utc").mjd)
    except Exception:
        try:
            return float(Time(str(date_value), scale="utc").mjd)
        except Exception:
            return None


def _signed_area(points: list[tuple[float, float]]) -> float:
    return 0.5 * sum(
        points[i][0] * points[(i + 1) % len(points)][1]
        - points[(i + 1) % len(points)][0] * points[i][1]
        for i in range(len(points))
    ) if len(points) >= 3 else 0.0


def _footprint_polygon(product, origin: SkyCoord) -> list[tuple[float, float]]:
    if product.wcs is None or not product.wcs.has_celestial:
        return []
    height, width = product.shape
    x = np.array([-0.5, width - 0.5, width - 0.5, -0.5], dtype=float)
    y = np.array([-0.5, -0.5, height - 0.5, height - 0.5], dtype=float)
    sky = product.wcs.pixel_to_world(x, y)
    offset = sky.transform_to(SkyOffsetFrame(origin=origin))
    points = list(zip(offset.lon.to_value(u.arcsec).tolist(), offset.lat.to_value(u.arcsec).tolist()))
    if _signed_area(points) < 0:
        points.reverse()
    return points


def _line_intersection(
    start: tuple[float, float], end: tuple[float, float],
    edge_start: tuple[float, float], edge_end: tuple[float, float],
) -> tuple[float, float]:
    sx, sy = start
    dx, dy = end[0] - sx, end[1] - sy
    ex, ey = edge_end[0] - edge_start[0], edge_end[1] - edge_start[1]
    denominator = dx * ey - dy * ex
    if abs(denominator) < 1e-14:
        return end
    t = ((edge_start[0] - sx) * ey - (edge_start[1] - sy) * ex) / denominator
    return sx + t * dx, sy + t * dy


def _inside(point: tuple[float, float], edge_start: tuple[float, float], edge_end: tuple[float, float]) -> bool:
    return (edge_end[0] - edge_start[0]) * (point[1] - edge_start[1]) - (edge_end[1] - edge_start[1]) * (point[0] - edge_start[0]) >= -1e-8


def _polygon_intersection(subject: list[tuple[float, float]], clip: list[tuple[float, float]]) -> list[tuple[float, float]]:
    output = subject
    for index, edge_start in enumerate(clip):
        edge_end = clip[(index + 1) % len(clip)]
        incoming = output
        output = []
        if not incoming:
            break
        previous = incoming[-1]
        for current in incoming:
            current_inside = _inside(current, edge_start, edge_end)
            previous_inside = _inside(previous, edge_start, edge_end)
            if current_inside:
                if not previous_inside:
                    output.append(_line_intersection(previous, current, edge_start, edge_end))
                output.append(current)
            elif previous_inside:
                output.append(_line_intersection(previous, current, edge_start, edge_end))
            previous = current
    return output


def _spatial_overlap(product1, product2) -> dict[str, Any]:
    if product1.wcs is None or product2.wcs is None:
        return {"fraction_of_smaller_footprint": 0.0, "area_arcsec2": 0.0, "method": "unavailable: celestial WCS required"}
    try:
        centre = product1.wcs.pixel_to_world((product1.shape[1] - 1) / 2, (product1.shape[0] - 1) / 2)
        p1 = _footprint_polygon(product1, centre)
        p2 = _footprint_polygon(product2, centre)
        area1, area2 = abs(_signed_area(p1)), abs(_signed_area(p2))
        intersection = _polygon_intersection(p1, p2)
        area = abs(_signed_area(intersection))
        fraction = min(1.0, max(0.0, area / min(area1, area2))) if min(area1, area2) > 0 else 0.0
        return {
            "fraction_of_smaller_footprint": fraction,
            "area_arcsec2": area,
            "footprint1_area_arcsec2": area1,
            "footprint2_area_arcsec2": area2,
            "method": "FITS WCS pixel-edge corners projected to a shared tangent plane; convex polygon intersection",
        }
    except Exception as exc:
        return {"fraction_of_smaller_footprint": 0.0, "area_arcsec2": 0.0, "method": f"failed: {exc}"}


def compatibility_matrix(records: list[dict], products: dict[str, Any], cfg: dict) -> list[dict]:
    """Report photometric and astrometric compatibility for every epoch pair."""
    ordered = sorted(records, key=lambda row: (_time_mjd(row) is None, _time_mjd(row) or 0.0, str(row.get("id", ""))))
    rows = []
    minimum_overlap = float(cfg.get("min_overlap_fraction", 0.55))
    for index1, index2 in itertools.combinations(range(len(ordered)), 2):
        record1, record2 = ordered[index1], ordered[index2]
        product1, product2 = products[record1["id"]], products[record2["id"]]
        check = wcs_alignment.check_observation_compatibility(record1, record2, product1, product2, cfg)
        overlap = _spatial_overlap(product1, product2)
        mode1 = record1.get("data_mode") or product1.data_mode
        mode2 = record2.get("data_mode") or product2.data_mode
        mjd1, mjd2 = _time_mjd(record1), _time_mjd(record2)
        time_baseline = abs(mjd2 - mjd1) if mjd1 is not None and mjd2 is not None else None
        same_mode = mode1 == mode2 == product1.data_mode == product2.data_mode
        has_wcs = bool(product1.wcs and product2.wcs and product1.wcs.has_celestial and product2.wcs.has_celestial)
        overlap_fraction = float(overlap["fraction_of_smaller_footprint"])
        time_ok = time_baseline is not None and time_baseline > 0
        astro_reasons = []
        if not same_mode:
            astro_reasons.append("Data modes differ; mixing real, local, and synthetic inputs is refused.")
        if not has_wcs:
            astro_reasons.append("Both products need valid celestial WCS.")
        if overlap_fraction < minimum_overlap:
            astro_reasons.append(f"Spatial overlap {overlap_fraction:.3f} is below the {minimum_overlap:.3f} requirement.")
        if not time_ok:
            astro_reasons.append("Distinct observation times are required to fit motion.")
        photometric_reasons = list(check["reasons"])
        if overlap_fraction < minimum_overlap:
            photometric_reasons.append(f"Spatial overlap {overlap_fraction:.3f} is below the {minimum_overlap:.3f} requirement.")
        rows.append({
            "observation_id_1": record1.get("id"),
            "observation_id_2": record2.get("id"),
            "mjd_1": mjd1,
            "mjd_2": mjd2,
            "time_separation_days": time_baseline,
            "astrometric_compatible": not astro_reasons,
            "astrometric_reasons": astro_reasons,
            "photometric_compatible": bool(check["compatible"] and overlap_fraction >= minimum_overlap),
            "photometric_reasons": photometric_reasons,
            "photometric_warnings": list(check["warnings"]),
            "wcs_wave_delta_um": check.get("wavelength_delta_um"),
            "wcs_wave_tolerance_um": check.get("wavelength_tolerance_um"),
            "spatial_overlap_fraction": overlap_fraction,
            "spatial_overlap_area_arcsec2": overlap["area_arcsec2"],
            "spatial_coverage_method": overlap["method"],
        })
    return rows


def _position_uncertainty_arcsec(product, record: dict) -> float:
    scale = float(product.pixel_scale_arcsec) if product.pixel_scale_arcsec and product.pixel_scale_arcsec > 0 else 6.15
    try:
        archive_resolution = float(record.get("position_resolution_arcsec"))
        if not np.isfinite(archive_resolution) or archive_resolution <= 0:
            archive_resolution = 0.0
    except (TypeError, ValueError):
        archive_resolution = 0.0
    # Half a pixel is a conservative floor for this undersampled peak-based
    # centroid; archive position_resolution is also retained when available.
    return max(scale / 2.0, archive_resolution)


def _detect_sources(product, record: dict, cfg: dict, reference: SkyCoord) -> tuple[list[dict], dict]:
    image = np.asarray(product.science, dtype=float)
    finite = np.isfinite(image)
    valid_values = image[finite]
    if valid_values.size < 16:
        return [], {"status": "skipped", "detail": "Too few finite image pixels."}
    _, background, robust_std = sigma_clipped_stats(valid_values, sigma=3.0, maxiters=5)
    if not np.isfinite(robust_std) or robust_std <= 0:
        return [], {"status": "skipped", "detail": "A robust per-image noise estimate could not be measured."}

    noise_source = "robust sigma-clipped image noise"
    noise = np.full(image.shape, float(robust_std), dtype=float)
    if product.variance is not None and np.asarray(product.variance).shape == image.shape:
        variance = np.asarray(product.variance, dtype=float)
        usable = np.isfinite(variance) & (variance > 0)
        if np.mean(usable & finite) >= 0.5:
            noise = np.where(usable, np.sqrt(np.maximum(variance, 0.0)), np.nan)
            noise_source = "per-pixel SPHEREx VARIANCE extension"
    with np.errstate(divide="ignore", invalid="ignore"):
        significance = (image - float(background)) / noise
    valid = finite & np.isfinite(significance) & np.isfinite(noise) & (noise > 0)
    threshold = float(cfg.get("detection_sigma", 5.0))
    border = int(cfg.get("border_px", 2))
    max_filter = maximum_filter(np.where(valid, significance, -np.inf), size=3, mode="constant", cval=-np.inf)
    peaks = valid & (significance >= threshold) & (significance == max_filter)
    if border > 0:
        peaks[:border] = False
        peaks[-border:] = False
        peaks[:, :border] = False
        peaks[:, -border:] = False

    peak_y, peak_x = np.where(peaks)
    order = np.argsort(significance[peak_y, peak_x])[::-1]
    min_distance = max(1.0, float(cfg.get("multi_epoch_min_peak_distance_px", 2.0)))
    selected: list[tuple[int, int]] = []
    for idx in order:
        x, y = int(peak_x[idx]), int(peak_y[idx])
        if any(math.hypot(x - sx, y - sy) < min_distance for sx, sy in selected):
            continue
        selected.append((x, y))
        if len(selected) >= int(cfg.get("max_candidates", 80)):
            break

    bad_mask = np.asarray(product.bad_flag_mask, dtype=bool) if product.bad_flag_mask is not None else np.zeros(image.shape, bool)
    pixel_scale = float(product.pixel_scale_arcsec) if product.pixel_scale_arcsec else 6.15
    position_sigma = _position_uncertainty_arcsec(product, record)
    detections = []
    for number, (x, y) in enumerate(selected, start=1):
        y0, y1 = max(0, y - 1), min(image.shape[0], y + 2)
        x0, x1 = max(0, x - 1), min(image.shape[1], x + 2)
        local = image[y0:y1, x0:x1]
        local_valid = valid[y0:y1, x0:x1]
        weights = np.where(local_valid, np.maximum(local - float(background), 0.0), 0.0)
        if weights.sum() > 0:
            yy, xx = np.indices(local.shape)
            cx = float(x0 + np.sum((xx + 0.0) * weights) / weights.sum())
            cy = float(y0 + np.sum((yy + 0.0) * weights) / weights.sum())
        else:
            cx, cy = float(x), float(y)
        neighborhood_bad = bad_mask[y0:y1, x0:x1]
        severe_fraction = float(np.mean(neighborhood_bad)) if neighborhood_bad.size else 0.0
        support = int(np.sum((significance[y0:y1, x0:x1] >= threshold * 0.25) & local_valid))
        artifact_flags = []
        if severe_fraction > 0:
            artifact_flags.append("NEAR_SEVERE_FLAG")
        if support <= 1:
            artifact_flags.append("SINGLE_PIXEL_LIKE_PEAK")
        if min(x, y, image.shape[1] - 1 - x, image.shape[0] - 1 - y) <= border:
            artifact_flags.append("NEAR_IMAGE_EDGE")
        try:
            sky = product.wcs.pixel_to_world(cx, cy)
            offset = sky.transform_to(SkyOffsetFrame(origin=reference))
            xi = float(offset.lon.to_value(u.arcsec))
            eta = float(offset.lat.to_value(u.arcsec))
        except Exception:
            continue
        detections.append({
            "observation_id": record.get("id"),
            "observation_date": record.get("date_utc"),
            "mjd": _time_mjd(record),
            "x_pixel": cx,
            "y_pixel": cy,
            "ra_deg": float(sky.icrs.ra.deg),
            "dec_deg": float(sky.icrs.dec.deg),
            "xi_arcsec": xi,
            "eta_arcsec": eta,
            "peak_snr": float(significance[y, x]),
            "peak_value": float(image[y, x]),
            "position_sigma_arcsec": position_sigma,
            "severe_flag_fraction": severe_fraction,
            "support_pixels": support,
            "artifact_flags": artifact_flags,
            "wavelength_um": product.center_wavelength_um,
            "data_mode": record.get("data_mode") or product.data_mode,
            "release": record.get("release"),
            "sha256": product.sha256,
        })
    return detections, {
        "status": "ok",
        "detection_threshold_sigma": threshold,
        "noise_method": noise_source,
        "background": float(background),
        "robust_image_sigma": float(robust_std),
        "position_uncertainty_floor_arcsec": position_sigma,
        "position_uncertainty_assumption": "max(half the WCS pixel scale, archive position_resolution_arcsec when available)",
        "flagged_pixel_fraction": float(np.mean(bad_mask)),
        "sources_detected": len(detections),
        "note": "Peak centroids are exploratory; the uncertainty floor is not a calibrated SPHEREx astrometric covariance.",
    }


def _predict_track(track: dict, mjd: float | None) -> tuple[float, float] | None:
    detections = track["detections"]
    if not detections:
        return None
    last = detections[-1]
    if len(detections) < 2 or mjd is None or last.get("mjd") is None:
        return last["xi_arcsec"], last["eta_arcsec"]
    previous = detections[-2]
    dt = last["mjd"] - previous.get("mjd", last["mjd"])
    if not np.isfinite(dt) or abs(dt) < 1e-9:
        return last["xi_arcsec"], last["eta_arcsec"]
    future = mjd - last["mjd"]
    vx = (last["xi_arcsec"] - previous["xi_arcsec"]) / dt
    vy = (last["eta_arcsec"] - previous["eta_arcsec"]) / dt
    return last["xi_arcsec"] + vx * future, last["eta_arcsec"] + vy * future


def _associate_epochs(epoch_detections: list[list[dict]], cfg: dict, compatible_pairs: set[tuple[str, str]]) -> list[dict]:
    tracks: list[dict] = []
    base_radius = float(cfg.get("multi_epoch_association_radius_arcsec", 12.0))
    max_speed = float(cfg.get("multi_epoch_max_motion_arcsec_per_day", 30.0))
    for epoch_index, detections in enumerate(epoch_detections):
        if not tracks:
            tracks = [{"detections": [detection], "quality_flags": [], "association_ambiguous": False} for detection in detections]
            continue
        active = [track for track in tracks if track["detections"] and track["detections"][-1].get("mjd") is not None]
        if not active or not detections:
            tracks.extend({"detections": [detection], "quality_flags": [], "association_ambiguous": False} for detection in detections)
            continue
        cost = np.full((len(active), len(detections)), 1e12, dtype=float)
        gates = np.zeros(len(active), dtype=float)
        ambiguities: dict[int, bool] = {}
        for row, track in enumerate(active):
            previous_id = str(track["detections"][-1].get("observation_id"))
            current_id = str(detections[0].get("observation_id"))
            if (previous_id, current_id) not in compatible_pairs:
                continue
            last_mjd = track["detections"][-1].get("mjd")
            current_mjd = detections[0].get("mjd")
            elapsed = max(0.0, float(current_mjd - last_mjd)) if current_mjd is not None and last_mjd is not None else 0.0
            prediction = _predict_track(track, current_mjd)
            if prediction is None:
                continue
            uncertainty = max(d["position_sigma_arcsec"] for d in track["detections"][-2:] + [detections[0]])
            gate = base_radius + max_speed * elapsed + 3.0 * uncertainty
            gates[row] = gate
            distances = np.array([
                math.hypot(detection["xi_arcsec"] - prediction[0], detection["eta_arcsec"] - prediction[1])
                for detection in detections
            ])
            eligible = np.where(distances <= gate)[0]
            if eligible.size > 1:
                nearest = np.sort(distances[eligible])
                ambiguities[row] = (nearest[1] - nearest[0]) <= max(2.0 * uncertainty, 1.0)
            for column in eligible:
                cost[row, column] = distances[column] / max(uncertainty, 0.1)

        assigned_detections: set[int] = set()
        if cost.size and np.any(cost < 1e12):
            rows, columns = linear_sum_assignment(cost)
            for row, column in zip(rows, columns):
                if cost[row, column] >= 1e12:
                    continue
                active[row]["detections"].append(detections[column])
                active[row]["association_gates_arcsec"] = active[row].get("association_gates_arcsec", []) + [float(gates[row])]
                if ambiguities.get(row):
                    active[row]["association_ambiguous"] = True
                    if "ASSOCIATION_AMBIGUOUS" not in active[row]["quality_flags"]:
                        active[row]["quality_flags"].append("ASSOCIATION_AMBIGUOUS")
                assigned_detections.add(int(column))
        tracks.extend(
            {"detections": [detection], "quality_flags": [], "association_ambiguous": False}
            for column, detection in enumerate(detections) if column not in assigned_detections
        )
    return tracks


def _fit_motion(detections: list[dict]) -> dict:
    usable = [d for d in detections if d.get("mjd") is not None]
    times = np.array([float(d["mjd"]) for d in usable], dtype=float)
    if len(usable) < 2 or len(np.unique(times)) < 2:
        return {"status": "unavailable", "reason": "At least two detections at distinct observation times are required for a displacement fit."}
    sigmas = np.array([max(float(d["position_sigma_arcsec"]), 1e-6) for d in usable])
    weights = 1.0 / sigmas**2
    t0 = float(np.average(times, weights=weights))
    design = np.column_stack([np.ones(len(times)), times - t0])
    normal = design.T @ (weights[:, None] * design)
    if np.linalg.matrix_rank(normal) < 2:
        return {"status": "unavailable", "reason": "Motion fit is singular; observation times are not sufficiently distinct."}
    covariance = np.linalg.inv(normal)
    outputs = {}
    residuals = {}
    chi2 = 0.0
    for axis in ("xi_arcsec", "eta_arcsec"):
        values = np.array([float(d[axis]) for d in usable])
        parameters = covariance @ (design.T @ (weights * values))
        residual = values - design @ parameters
        outputs[axis] = {"intercept_arcsec": float(parameters[0]), "slope_arcsec_per_day": float(parameters[1]), "slope_sigma_arcsec_per_day": float(np.sqrt(max(covariance[1, 1], 0.0)))}
        residuals[axis] = residual
        chi2 += float(np.sum(weights * residual**2))
    dof = max(0, 2 * len(usable) - 4)
    vx = outputs["xi_arcsec"]["slope_arcsec_per_day"]
    vy = outputs["eta_arcsec"]["slope_arcsec_per_day"]
    speed = math.hypot(vx, vy)
    if speed > 0:
        speed_sigma = math.sqrt(
            (vx / speed * outputs["xi_arcsec"]["slope_sigma_arcsec_per_day"]) ** 2
            + (vy / speed * outputs["eta_arcsec"]["slope_sigma_arcsec_per_day"]) ** 2
        )
    else:
        speed_sigma = math.hypot(outputs["xi_arcsec"]["slope_sigma_arcsec_per_day"], outputs["eta_arcsec"]["slope_sigma_arcsec_per_day"])
    return {
        "status": "ok",
        "reference_mjd": t0,
        "xi_arcsec_per_day": vx,
        "eta_arcsec_per_day": vy,
        "speed_arcsec_per_day": speed,
        "speed_sigma_arcsec_per_day": speed_sigma,
        "motion_signal_to_uncertainty": speed / speed_sigma if speed_sigma > 0 else None,
        "position_angle_deg_east_of_north": float(np.degrees(np.arctan2(vx, vy)) % 360.0),
        "reduced_chi2": chi2 / dof if dof > 0 else None,
        "degrees_of_freedom": dof,
        "fit_axes": outputs,
        "residuals_arcsec": {axis: values.tolist() for axis, values in residuals.items()},
        "assumed_per_detection_sigma_arcsec": sigmas.tolist(),
        "error_model": "weighted linear fit using a conservative position floor; formal covariance is conditional on this uncalibrated floor",
    }


def _track_report(track_number: int, track: dict, searchable_epochs: int, cfg: dict) -> dict:
    detections = track["detections"]
    detections = sorted(detections, key=lambda d: (d.get("mjd") is None, d.get("mjd") or 0.0, d["observation_id"]))
    motion = _fit_motion(detections)
    flags = list(track.get("quality_flags", []))
    evidence = [f"Associated {len(detections)} peak detection(s) across {len({d['observation_id'] for d in detections})} observation(s)."]
    limitations = [
        "Per-epoch positions come from peak-based centroids on undersampled images; they are not PSF/ePSF astrometry.",
        "The positional uncertainty floor is conservative but not calibrated against SPHEREx astrometric residuals.",
        "Sparse observations can alias unrelated sources into an apparently smooth trajectory; review detections and association gates.",
    ]
    flagged = [d for d in detections if d["severe_flag_fraction"] > 0.25 or "NEAR_IMAGE_EDGE" in d["artifact_flags"]]
    if flagged:
        flags.append("SEVERE_FLAG_OR_EDGE_PROXIMITY")
        classification = "instrumental_artifact"
        evidence.append(f"{len(flagged)} detection(s) overlap severe FLAGS or the image edge.")
    elif len(detections) == 1:
        flags.append("SINGLE_EPOCH_ONLY")
        classification = "insufficient_evidence"
        evidence.append("Only one epoch has a detection; a transient cannot be distinguished from a threshold or coverage effect.")
    elif len(detections) < int(cfg.get("multi_epoch_min_epochs_for_motion", 3)):
        if len(detections) == 2:
            flags.append("2-EPOCH-ONLY")
            evidence.append("2-EPOCH-ONLY: two positions define only a displacement; they cannot establish a fitted motion track.")
        classification = "insufficient_evidence"
    elif motion.get("status") != "ok":
        classification = "insufficient_evidence"
        evidence.append(motion.get("reason", "Motion fit unavailable."))
    elif track.get("association_ambiguous"):
        classification = "insufficient_evidence"
        evidence.append("At least one association had a nearby alternative within the astrometric ambiguity tolerance.")
    elif motion.get("reduced_chi2") is not None and motion["reduced_chi2"] > float(cfg.get("multi_epoch_max_reduced_chi2", 25.0)):
        flags.append("POOR_LINEAR_MOTION_FIT")
        classification = "insufficient_evidence"
        evidence.append(f"The weighted linear motion fit has reduced chi-square {motion['reduced_chi2']:.2f}.")
    elif float(motion.get("motion_signal_to_uncertainty") or 0.0) >= float(cfg.get("multi_epoch_min_motion_snr", 3.0)):
        classification = "potential_moving_source"
        evidence.append(
            f"The fitted tangent-plane speed is {motion['speed_arcsec_per_day']:.4g} ± "
            f"{motion['speed_sigma_arcsec_per_day']:.3g} arcsec/day (uncertainty conditional on the stated floor)."
        )
    else:
        classification = "stationary_consistent_with_zero"
        evidence.append("The fitted motion is not significant relative to the stated positional uncertainty floor.")
    if any("SINGLE_PIXEL_LIKE_PEAK" in d["artifact_flags"] for d in detections):
        flags.append("SINGLE_PIXEL_LIKE_PEAK")
        limitations.append("At least one peak is single-pixel-like; undersampling prevents morphology alone from confirming a point source.")
    if len(detections) < searchable_epochs:
        limitations.append("Non-detections are not treated as transients; depth, crowding, wavelength response and coverage vary by epoch.")
    return {
        "track_id": f"TRACK-{track_number:04d}",
        "classification": classification,
        "status_label": {
            "potential_moving_source": "Potential moving-source candidate (investigative lead)",
            "stationary_consistent_with_zero": "Motion consistent with zero at this uncertainty floor",
            "instrumental_artifact": "Instrumental/artifact-suspect detections",
            "insufficient_evidence": "Insufficient evidence",
        }.get(classification, classification),
        "detection_count": len(detections),
        "required_epochs_for_motion": int(cfg.get("multi_epoch_min_epochs_for_motion", 3)),
        "motion": motion,
        "quality_flags": sorted(set(flags)),
        "association_ambiguous": bool(track.get("association_ambiguous")),
        "association_gates_arcsec": track.get("association_gates_arcsec", []),
        "detections": detections,
        "evidence": evidence,
        "limitations": limitations,
        "disclaimer": "This is an observational screening record, not an object identification or discovery.",
    }


def analyse_multi_epoch(records: list[dict], products: dict[str, Any] | list[Any], cfg: dict | None = None, use_cache: bool = False) -> dict:
    """Associate detections across 2+ observations and fit ≥3-epoch motion.

    Two-epoch input is allowed so its limitation can be reported explicitly.
    A moving-source classification is only possible with at least three epochs.
    """
    if not records:
        raise ValueError("At least one observation record is required.")
    if len(records) > MAX_MULTI_EPOCHS:
        raise ValueError(f"Select no more than {MAX_MULTI_EPOCHS} epochs per analysis to keep coverage checks bounded.")
    cfg_merged = dict(config.PROCESSING_DEFAULTS)
    cfg_merged.update(DEFAULT_MULTI_EPOCH_CONFIG)
    if cfg:
        cfg_merged.update(cfg)
    if isinstance(products, dict):
        product_map = products
    else:
        product_map = {product.product_id: product for product in products}
    ordered_records = sorted(records, key=lambda row: (_time_mjd(row) is None, _time_mjd(row) or 0.0, str(row.get("id", ""))))
    missing = [record.get("id") for record in ordered_records if record.get("id") not in product_map]
    if missing:
        raise ValueError(f"FITS products are missing for observation IDs: {', '.join(map(str, missing))}")
    if len(ordered_records) == 1:
        return {
            "status": "insufficient_epochs",
            "tracks": [],
            "limitations": ["One epoch cannot establish motion. Select at least two observations; at least three are required for a motion track."],
        }
    if any(_time_mjd(record) is None for record in ordered_records):
        time_warning = "Some observation times are missing; those detections cannot contribute to a motion fit."
    else:
        time_warning = None

    ordered_products = {record["id"]: product_map[record["id"]] for record in ordered_records}
    data_modes = {record.get("data_mode") or ordered_products[record["id"]].data_mode for record in ordered_records}
    if len(data_modes) != 1:
        return {
            "status": "incompatible_data_modes",
            "tracks": [],
            "error": "Real archive data, local uploads, and synthetic demonstration data cannot be combined in one multi-epoch analysis.",
            "compatibility_matrix": [],
        }
    pair_rows = compatibility_matrix(ordered_records, ordered_products, cfg_merged)
    compatible_pairs = {
        (str(pair["observation_id_1"]), str(pair["observation_id_2"]))
        for pair in pair_rows if pair["astrometric_compatible"]
    }
    compatible_pairs |= {(right, left) for left, right in compatible_pairs}
    reference_product = ordered_products[ordered_records[0]["id"]]
    if reference_product.wcs is None:
        return {"status": "failed", "error": "Celestial WCS is required for multi-epoch sky association.", "compatibility_matrix": pair_rows}
    reference = reference_product.wcs.pixel_to_world((reference_product.shape[1] - 1) / 2, (reference_product.shape[0] - 1) / 2)
    epoch_detections: list[list[dict]] = []
    epoch_diagnostics = []
    for record in ordered_records:
        product = ordered_products[record["id"]]
        if product.wcs is None:
            epoch_detections.append([])
            epoch_diagnostics.append({"observation_id": record.get("id"), "status": "failed", "detail": "Celestial WCS unavailable."})
            continue
        detections, diagnostics = _detect_sources(product, record, cfg_merged, reference)
        epoch_detections.append(detections)
        epoch_diagnostics.append({"observation_id": record.get("id"), **diagnostics})
    raw_tracks = _associate_epochs(epoch_detections, cfg_merged, compatible_pairs)
    searchable_epochs = len(ordered_records)
    reports = [
        _track_report(index, track, searchable_epochs, cfg_merged)
        for index, track in enumerate(raw_tracks, start=1)
        if track["detections"]
    ]
    reports.sort(key=lambda item: (item["classification"] != "potential_moving_source", -item["detection_count"], item["track_id"]))
    # Keep IDs deterministic after sorting.
    for index, report in enumerate(reports, start=1):
        report["track_id"] = f"TRACK-{index:04d}"

    data_mode = next(iter(data_modes))
    limitations = [
        "Motion is fitted in a local ICRS tangent plane against sparse observations; it is not a continuous orbit or trajectory.",
        "The matching gate and peak-detection threshold affect which tracks are found; inspect their values and each epoch association.",
        "No calibrated false-positive probability or discovery claim is produced.",
    ]
    if len(ordered_records) < 3:
        limitations.insert(0, "2-EPOCH-ONLY: two positions provide a displacement, not a fitted motion track. No moving-source classification is allowed.")
    if time_warning:
        limitations.append(time_warning)
    run_inputs = []
    for record in ordered_records:
        product = ordered_products[record["id"]]
        run_inputs.append({
            "id": record.get("id"),
            "sha256": product.sha256,
            "access_url": record.get("cutout_url") or record.get("full_fits_url") or record.get("cached_path"),
            "origin": record.get("origin"),
            "date_utc": record.get("date_utc"),
            "mjd": _time_mjd(record),
            "band": record.get("band"),
            "coordinates": coordinates.coordinate_context(record),
        })
    run_cfg = {**cfg_merged, "analysis": "multi_epoch_peak_track_v1"}
    run_id = provenance.deterministic_run_id(run_inputs, run_cfg)
    warnings_out = [time_warning] if time_warning else []
    stages = [
        {"stage": "Pair compatibility", "status": "ok", "detail": f"Evaluated {len(pair_rows)} observation pair(s); astrometric and photometric compatibility are separate."},
        {"stage": "Per-epoch source detection", "status": "ok", "detail": f"Detected {sum(map(len, epoch_detections))} peak(s) across {len(ordered_records)} epoch(s)."},
        {"stage": "Source association and motion fit", "status": "ok", "detail": f"Built {len(reports)} candidate track(s); motion classifications require at least {cfg_merged['multi_epoch_min_epochs_for_motion']} detections."},
    ]
    provenance_record = provenance.build_provenance(
        inputs=run_inputs,
        processing_config=run_cfg,
        stages=stages,
        warnings=warnings_out + limitations,
        data_mode=data_mode,
        result_summary={
            "observations": len(ordered_records),
            "detections": sum(map(len, epoch_detections)),
            "tracks": len(reports),
            "potential_moving_source_tracks": sum(t["classification"] == "potential_moving_source" for t in reports),
            "compatibility_pairs": len(pair_rows),
            "disclaimer": "Multi-epoch screening only; no discovery or object identification.",
        },
        run_id=run_id,
    )
    provenance_path = config.OUTPUT_DIR / "multi_epoch" / f"{run_id}.json"
    provenance.write_json(provenance_path, provenance_record)
    return {
        "status": "ok",
        "run_id": run_id,
        "data_mode": data_mode,
        "reference_coordinate_icrs_deg": [float(reference.icrs.ra.deg), float(reference.icrs.dec.deg)],
        "compatibility_matrix": pair_rows,
        "epoch_diagnostics": epoch_diagnostics,
        "detections": [detection for epoch in epoch_detections for detection in epoch],
        "tracks": reports,
        "summary": {
            "observations": len(ordered_records),
            "detections": sum(map(len, epoch_detections)),
            "tracks": len(reports),
            "potential_moving_source_tracks": sum(t["classification"] == "potential_moving_source" for t in reports),
            "insufficient_evidence_tracks": sum(t["classification"] == "insufficient_evidence" for t in reports),
        },
        "configuration": cfg_merged,
        "stages": stages,
        "limitations": limitations,
        "warnings": warnings_out,
        "provenance": provenance_record,
        "provenance_path": str(provenance_path),
    }

