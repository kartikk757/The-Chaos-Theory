"""Observation compatibility checks and reprojection onto a common celestial grid."""
from __future__ import annotations

import copy
import logging
import warnings
from typing import Any

import astropy.units as u
import numpy as np
from astropy.coordinates import SkyCoord
from reproject import reproject_interp

import config
from modules import coordinates
from modules.models import AlignedPair, FitsProduct


class AlignmentError(RuntimeError):
    pass


def _bandwidth_estimate(product: FitsProduct) -> float | None:
    if product.center_bandwidth_um is not None and product.center_bandwidth_um > 0:
        return float(product.center_bandwidth_um)
    if product.center_wavelength_um is not None:
        band = config.DETECTOR_BANDS.get(product.metadata.get("observation", {}).get("band", ""), None)
        if band and band["R"]:
            return float(product.center_wavelength_um) / float(band["R"])
    return None


def check_observation_compatibility(
    record1: dict,
    record2: dict,
    product1: FitsProduct,
    product2: FitsProduct,
    cfg: dict,
) -> dict:
    checks: list[dict[str, str]] = []
    warnings_out: list[str] = []
    reasons: list[str] = []

    def add(name: str, ok: bool, detail: str, severity: str = "fail") -> None:
        checks.append({"name": name, "status": "pass" if ok else severity, "detail": detail})
        if not ok and severity == "fail":
            reasons.append(f"{name}: {detail}")

    same_band = record1.get("band") == record2.get("band")
    add(
        "Detector / LVF band",
        same_band,
        f"{record1.get('band')} vs {record2.get('band')}."
        + ("" if same_band else " Different detectors sample different wavelength ranges, so a difference would mix colour with variability."),
    )
    mode1 = record1.get("data_mode") or product1.data_mode
    mode2 = record2.get("data_mode") or product2.data_mode
    add(
        "Data mode",
        mode1 == mode2 and product1.data_mode == product2.data_mode,
        f"{mode1} vs {mode2}; real archive data, local uploads and synthetic demonstration products must never be mixed.",
    )
    unit_ok = (not product1.bunit or not product2.bunit) or product1.bunit == product2.bunit
    add("IMAGE BUNIT", unit_ok, f"{product1.bunit or 'not supplied'} vs {product2.bunit or 'not supplied'}")
    wcs_ok = product1.wcs is not None and product2.wcs is not None and product1.wcs.has_celestial and product2.wcs.has_celestial
    add(
        "Celestial WCS",
        wcs_ok,
        "Both products provide a celestial WCS." if wcs_ok else "At least one product has no celestial WCS, so reprojection is impossible.",
    )

    wl1, wl2 = product1.center_wavelength_um, product2.center_wavelength_um
    delta = None
    tolerance = None
    if wl1 is not None and wl2 is not None:
        delta = abs(float(wl1) - float(wl2))
        bw1, bw2 = _bandwidth_estimate(product1), _bandwidth_estimate(product2)
        reference_bandwidth = min([value for value in (bw1, bw2) if value is not None] or [0.0])
        tolerance = float(cfg["wavelength_tolerance_bandwidth_fraction"]) * reference_bandwidth if reference_bandwidth > 0 else None
        if tolerance is None:
            add("WCS-WAVE spectral compatibility", False, "Central wavelengths are available but bandwidths are not; compatibility cannot be quantified.")
        else:
            ok = delta <= tolerance
            add(
                "WCS-WAVE spectral compatibility",
                ok,
                f"Δλ = {delta:.5f} µm; allowed ≤ {tolerance:.5f} µm ({cfg['wavelength_tolerance_bandwidth_fraction']:.0%} of the narrower WCS-WAVE bandwidth).",
            )
            if not ok:
                warnings_out.append("A fixed sky position is observed at different LVF wavelengths in these two pointings; a direct brightness difference would be confounded by spectral response.")
    elif wl1 is None and wl2 is None and product1.data_mode == product2.data_mode and product1.data_mode in (
        config.DATA_MODES["LOCAL"], config.DATA_MODES["DEMO"]
    ):
        add("WCS-WAVE spectral compatibility", True, f"Both products are {product1.data_mode} and lack WCS-WAVE; spectral compatibility is explicitly unverified.")
        warnings_out.append(
            "Spectral compatibility is unverified because no WCS-WAVE lookup table is present in these files. "
            "This comparison is valid for pipeline demonstration only."
        )
    else:
        add("WCS-WAVE spectral compatibility", False, "Load both FITS cutouts so the WCS-WAVE central wavelength can be checked before comparison.")

    if record1.get("calibration_level") is not None and record2.get("calibration_level") is not None and record1["calibration_level"] != record2["calibration_level"]:
        warnings_out.append(f"Calibration levels differ ({record1['calibration_level']} vs {record2['calibration_level']}); check provenance before interpretation.")
    if record1.get("release") and record2.get("release") and record1["release"] != record2["release"]:
        warnings_out.append(f"Archive releases differ ({record1['release']} vs {record2['release']}); processing versions may differ even when pipeline metadata matches.")
    if record1.get("quality_flag") or record2.get("quality_flag"):
        warnings_out.append("At least one archive plane supplies a quality_flag; inspect it in Data Explorer.")

    record_position_separation = None
    footprint_separation = None
    coord1 = coordinates.coordinate_context(record1)
    coord2 = coordinates.coordinate_context(record2)
    record_sky1 = coord1.get("record_position_icrs_deg")
    record_sky2 = coord2.get("record_position_icrs_deg")
    if record_sky1 is not None and record_sky2 is not None:
        try:
            c1 = SkyCoord(float(record_sky1[0]) * u.deg, float(record_sky1[1]) * u.deg, frame="icrs")
            c2 = SkyCoord(float(record_sky2[0]) * u.deg, float(record_sky2[1]) * u.deg, frame="icrs")
            record_position_separation = float(c1.separation(c2).arcsec)
            if (
                coord1["record_position_role"] == "parent_observation_footprint_center"
                and coord2["record_position_role"] == "parent_observation_footprint_center"
            ):
                footprint_separation = record_position_separation
        except Exception:
            record_position_separation = None
            footprint_separation = None

    # Difference geometry always uses the delivered FITS WCS centres. Record and
    # parent-footprint positions are provenance only and never substitute for them.
    separation = None
    cutout_centres = [product1.center_sky, product2.center_sky]
    if all(centre is not None for centre in cutout_centres) and product1.pixel_scale_arcsec and product2.pixel_scale_arcsec:
        try:
            c1 = SkyCoord(cutout_centres[0][0] * u.deg, cutout_centres[0][1] * u.deg, frame="icrs")
            c2 = SkyCoord(cutout_centres[1][0] * u.deg, cutout_centres[1][1] * u.deg, frame="icrs")
            separation = float(c1.separation(c2).arcsec)
            fov1 = min(product1.shape) * product1.pixel_scale_arcsec
            fov2 = min(product2.shape) * product2.pixel_scale_arcsec
            allowed = 0.75 * min(fov1, fov2)
            add(
                "Cutout centre separation",
                separation <= allowed,
                f"Delivered cutout centres are {separation:.2f} arcsec apart; allowed ≤ {allowed:.1f} arcsec "
                f"(75% of the smaller 59-pixel cutout width). Exact pixel overlap is measured after reprojection.",
            )
        except Exception as exc:
            warnings_out.append(f"Cutout separation could not be computed ({exc}); overlap is measured after reprojection.")
    else:
        warnings_out.append("Cutout centre or pixel scale is unavailable; exact overlap will be measured from WCS footprints during alignment.")
    if footprint_separation is not None and footprint_separation > 600:
        warnings_out.append(
            f"The parent observation footprint centres are {footprint_separation / 3600:.2f}° apart. This is expected for "
            "SPHEREx strips and does not reduce comparability of the delivered cutouts, but it is recorded for transparency."
        )

    product_warnings = list(product1.warnings) + list(product2.warnings)
    return {
        "compatible": len(reasons) == 0,
        "checks": checks,
        "reasons": reasons,
        "warnings": warnings_out,
        "product_warnings": product_warnings,
        "wavelength_delta_um": delta,
        "wavelength_tolerance_um": tolerance,
        "center_separation_arcsec": separation,
        "record_position_separation_arcsec": record_position_separation,
        "footprint_center_separation_arcsec": footprint_separation,
        "time_separation_days": (
            abs(float(record2["mjd"]) - float(record1["mjd"]))
            if record1.get("mjd") is not None and record2.get("mjd") is not None else None
        ),
    }


def _simple_celestial_wcs(wcs, shape: tuple[int, int]):
    result = copy.deepcopy(wcs.celestial)
    result.array_shape = shape
    if hasattr(result, "sip"):
        result.sip = None
    try:
        result.wcs.ctype = [str(ctype).replace("-SIP", "") for ctype in result.wcs.ctype]
    except Exception:
        pass
    return result


def _reproject(data: np.ndarray, input_wcs, target_wcs, shape_out: tuple[int, int], order: str) -> tuple[np.ndarray, np.ndarray]:
    wcs = copy.deepcopy(input_wcs)
    wcs.array_shape = tuple(data.shape)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        previous = logging.getLogger("astropy").level
        logging.getLogger("astropy").setLevel(logging.ERROR)
        projected, footprint = reproject_interp(
            (np.asarray(data, dtype=float), wcs),
            target_wcs,
            shape_out=shape_out,
            order=order,
            return_footprint=True,
        )
        logging.getLogger("astropy").setLevel(previous)
    return np.asarray(projected, dtype=float), np.asarray(footprint, dtype=float)


def align_pair(product1: FitsProduct, product2: FitsProduct, cfg: dict) -> AlignedPair:
    if product1.wcs is None or product2.wcs is None:
        raise AlignmentError("One or both products have no valid celestial WCS; direct comparison is not possible.")
    shape = tuple(int(v) for v in product1.shape)
    target_wcs = _simple_celestial_wcs(product1.wcs, shape)
    warnings_out: list[str] = []
    fallback_used = False
    input_wcs1, input_wcs2 = product1.wcs, product2.wcs

    science1_raw = np.asarray(product1.science, dtype=float)
    science2_raw = np.asarray(product2.science, dtype=float)
    invalid1 = product1.bad_flag_mask | ~np.isfinite(science1_raw)
    invalid2 = product2.bad_flag_mask | ~np.isfinite(science2_raw)
    science1_filled = np.where(invalid1, 0.0, science1_raw)
    science2_filled = np.where(invalid2, 0.0, science2_raw)
    variance1_raw = np.asarray(product1.variance, dtype=float) if product1.variance is not None else None
    variance2_raw = np.asarray(product2.variance, dtype=float) if product2.variance is not None else None
    if variance1_raw is not None:
        variance1_raw = np.where(np.isfinite(variance1_raw) & (variance1_raw >= 0) & ~invalid1, variance1_raw, 0.0)
    if variance2_raw is not None:
        variance2_raw = np.where(np.isfinite(variance2_raw) & (variance2_raw >= 0) & ~invalid2, variance2_raw, 0.0)

    payloads = {
        "science1": (science1_filled, input_wcs1, "bilinear"),
        "science2": (science2_filled, input_wcs2, "bilinear"),
        "variance1": (variance1_raw if variance1_raw is not None else np.zeros(shape, dtype=float), input_wcs1, "bilinear"),
        "variance2": (variance2_raw if variance2_raw is not None else np.zeros(shape, dtype=float), input_wcs2, "bilinear"),
        "bad1": (invalid1.astype(float), input_wcs1, "nearest-neighbor"),
        "bad2": (invalid2.astype(float), input_wcs2, "nearest-neighbor"),
        "source1": (product1.source_flag_mask.astype(float), input_wcs1, "nearest-neighbor"),
        "source2": (product2.source_flag_mask.astype(float), input_wcs2, "nearest-neighbor"),
        "ones1": (np.ones(shape, dtype=float), input_wcs1, "nearest-neighbor"),
        "ones2": (np.ones(shape, dtype=float), input_wcs2, "nearest-neighbor"),
    }

    def run(wcs1, wcs2):
        out = {}
        for key, (data, wcs, order) in payloads.items():
            selected_wcs = wcs1 if key.endswith("1") else wcs2
            out[key] = _reproject(data, selected_wcs, target_wcs, shape, order)
        return out

    try:
        projected = run(input_wcs1, input_wcs2)
    except Exception as first_error:
        fallback_used = True
        warnings_out.append(f"Full SIP-aware reprojection failed ({first_error}); retried with SIP distortion removed from input WCS. Alignment is approximate.")
        simple1 = _simple_celestial_wcs(input_wcs1, product1.shape)
        simple2 = _simple_celestial_wcs(input_wcs2, product2.shape)
        try:
            projected = run(simple1, simple2)
        except Exception as second_error:
            raise AlignmentError(f"Reprojection onto a common WCS failed: {second_error}") from second_error

    science1 = projected["science1"][0]
    science2 = projected["science2"][0]
    variance1 = projected["variance1"][0] if product1.variance is not None else None
    variance2 = projected["variance2"][0] if product2.variance is not None else None
    bad1 = projected["bad1"][0] > 0.5
    bad2 = projected["bad2"][0] > 0.5
    source1 = projected["source1"][0] > 0.5
    source2 = projected["source2"][0] > 0.5
    coverage1 = projected["ones1"][1]
    coverage2 = projected["ones2"][1]
    overlap = coverage1 > 0.8
    overlap &= coverage2 > 0.8
    overlap_fraction = float(np.mean(overlap))
    if overlap_fraction < float(cfg["min_overlap_fraction"]):
        raise AlignmentError(
            f"Insufficient spatial overlap: {overlap_fraction:.1%} of the target grid is covered by both WCS footprints; "
            f"the configured minimum is {float(cfg['min_overlap_fraction']):.0%}."
        )
    valid = overlap & np.isfinite(science1) & np.isfinite(science2) & ~bad1 & ~bad2
    if variance1 is not None:
        valid &= np.isfinite(variance1) & (variance1 > 0)
    if variance2 is not None:
        valid &= np.isfinite(variance2) & (variance2 > 0)
    science1[~valid] = np.nan
    science2[~valid] = np.nan
    if variance1 is not None:
        variance1[~valid] = np.nan
    if variance2 is not None:
        variance2[~valid] = np.nan
    border = int(cfg.get("border_px", 3))
    if border > 0:
        valid[:border, :] = False
        valid[-border:, :] = False
        valid[:, :border] = False
        valid[:, -border:] = False
        science1[~valid] = np.nan
        science2[~valid] = np.nan
    if np.mean(valid) < 0.10:
        raise AlignmentError(f"Only {np.mean(valid):.1%} of pixels remain valid after FLAGS, WCS and finite-value filtering.")

    try:
        matrix = np.asarray(target_wcs.celestial.pixel_scale_matrix, dtype=float)
        values = np.sqrt(np.sum(matrix ** 2, axis=0)) * 3600.0
        pixel_scale = float(np.mean(np.abs(values)))
    except Exception:
        pixel_scale = float("nan")
    if product1.zodi is None or product2.zodi is None:
        warnings_out.append("At least one ZODI extension is absent; zodiacal-model subtraction could not be applied to that epoch.")
    if variance1 is None or variance2 is None:
        warnings_out.append("At least one VARIANCE extension is absent; SNR will fall back to a global robust difference-noise estimate.")
    if fallback_used:
        warnings_out.append("SIP distortion was removed from the input WCS fallback; residual registration artifacts are more likely.")
    diagnostics = {
        "target_reference": "epoch 1 celestial WCS, SIP removed from output grid",
        "target_shape": shape,
        "input_shapes": [list(product1.shape), list(product2.shape)],
        "reprojection": "reproject.reproject_interp; bilinear for science/variance, nearest for masks",
        "sip_fallback_used": fallback_used,
        "overlap_fraction": overlap_fraction,
        "valid_pixel_fraction": float(np.mean(valid)),
        "severe_flag_fraction_epoch1": float(np.mean(product1.bad_flag_mask)),
        "severe_flag_fraction_epoch2": float(np.mean(product2.bad_flag_mask)),
        "pixel_scale_arcsec": pixel_scale,
        "zodi_subtracted": [product1.zodi is not None, product2.zodi is not None],
        "variance_available": [variance1 is not None, variance2 is not None],
        "units": product1.bunit or "IMAGE BUNIT not supplied",
    }
    return AlignedPair(
        target_wcs=target_wcs,
        shape=shape,
        science1=science1,
        science2=science2,
        variance1=variance1,
        variance2=variance2,
        bad_mask1=bad1,
        bad_mask2=bad2,
        source_mask1=source1,
        source_mask2=source2,
        coverage1=coverage1,
        coverage2=coverage2,
        valid_mask=valid,
        overlap_fraction=overlap_fraction,
        pixel_scale_arcsec=pixel_scale,
        diagnostics=diagnostics,
        warnings=warnings_out,
    )
