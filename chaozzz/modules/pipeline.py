"""End-to-end scientific pipeline orchestration (Phases 3–5).

Order is fixed and recorded in provenance:
validation → alignment → background/normalization → difference → detection →
candidate measurement → catalogue cross-match.
"""
from __future__ import annotations

import json
from typing import Any

import numpy as np
import streamlit as st

import config
from modules import candidate_analysis, catalogue_matching, data_loader, difference_engine
from modules import image_normalization, provenance, source_detection, wcs_alignment


def processing_config(overrides: dict | None = None, crossmatch_enabled: bool | None = None) -> dict:
    cfg = dict(config.PROCESSING_DEFAULTS)
    if overrides:
        cfg.update(overrides)
    if crossmatch_enabled is not None:
        cfg["crossmatch_enabled"] = bool(crossmatch_enabled)
    return cfg


@st.cache_data(max_entries=12, show_spinner=False)
def _run_pipeline_cached(
    record1_json: str,
    record2_json: str,
    config_json: str,
    catalogue_key: str,
) -> dict:
    return _run_pipeline(
        json.loads(record1_json),
        json.loads(record2_json),
        json.loads(config_json),
        catalogue_key,
    )


def run_pipeline(
    record1: dict,
    record2: dict,
    cfg: dict,
    catalogue_key: str = "2MASS PSC",
    use_cache: bool = True,
) -> dict:
    payload1 = json.dumps(provenance.safe_json(record1), sort_keys=True)
    payload2 = json.dumps(provenance.safe_json(record2), sort_keys=True)
    payload_cfg = json.dumps(provenance.safe_json(cfg), sort_keys=True)
    if use_cache:
        return _run_pipeline_cached(payload1, payload2, payload_cfg, catalogue_key)
    return _run_pipeline(json.loads(payload1), json.loads(payload2), json.loads(payload_cfg), catalogue_key)


def _run_pipeline(record1: dict, record2: dict, cfg: dict, catalogue_key: str) -> dict:
    stages: list[dict] = []
    warnings_out: list[str] = []
    result: dict[str, Any] = {"status": "ok", "stages": stages, "warnings": warnings_out}

    def stage(name: str, status: str, detail: str, **extra) -> None:
        stages.append({"stage": name, "status": status, "detail": detail, **provenance.safe_json(extra)})

    # 1. Load + validate
    try:
        product1 = data_loader.load_observation_product(record1)
        product2 = data_loader.load_observation_product(record2)
    except Exception as exc:
        stage("Load observations", "failed", str(exc))
        return {"status": "failed", "error": f"Observation loading failed: {exc}", "stages": stages, "warnings": warnings_out}
    stage(
        "Load observations",
        "ok",
        f"Loaded {product1.shape[0]}×{product1.shape[1]} and {product2.shape[0]}×{product2.shape[1]} pixel IMAGE arrays.",
        epoch1=record1.get("label"), epoch2=record2.get("label"),
    )
    stage(
        "Validate metadata",
        "ok" if all(check["status"] != "fail" for check in product1.validation + product2.validation) else "warn",
        f"{sum(1 for c in product1.validation + product2.validation if c['status'] == 'pass')} checks passed, "
        f"{sum(1 for c in product1.validation + product2.validation if c['status'] == 'warn')} warnings.",
        epoch1_checks=product1.validation, epoch2_checks=product2.validation,
    )

    # 2. Compatibility
    compatibility = wcs_alignment.check_observation_compatibility(record1, record2, product1, product2, cfg)
    if not compatibility["compatible"]:
        stage("Compatibility check", "failed", " | ".join(compatibility["reasons"]))
        return {
            "status": "incompatible",
            "error": "Observations are not directly comparable.",
            "compatibility": compatibility,
            "stages": stages,
            "warnings": compatibility["warnings"],
            "product1": product1,
            "product2": product2,
        }
    delta_wavelength = compatibility.get("wavelength_delta_um")
    delta_time = compatibility.get("time_separation_days")
    stage(
        "Compatibility check",
        "ok",
        f"Δλ = {delta_wavelength:.5f} µm" if delta_wavelength is not None else "Δλ unverified (no WCS-WAVE in the products)"
        + (f", Δt = {delta_time:.2f} d." if delta_time is not None else "."),
        compatibility=compatibility,
    )

    # 3. Alignment
    try:
        aligned = wcs_alignment.align_pair(product1, product2, cfg)
    except Exception as exc:
        stage("WCS alignment", "failed", str(exc))
        return {
            "status": "alignment_failed",
            "error": f"WCS alignment failed: {exc}",
            "compatibility": compatibility,
            "stages": stages,
            "warnings": compatibility["warnings"],
            "product1": product1,
            "product2": product2,
        }
    warnings_out.extend(aligned.warnings)
    stage(
        "WCS alignment",
        "ok",
        f"Reprojected onto epoch-1 grid; overlap {aligned.overlap_fraction:.1%}, valid pixels {np.mean(aligned.valid_mask):.1%}.",
        diagnostics=aligned.diagnostics,
    )

    # 4. Background + normalization
    normalized = image_normalization.normalize_pair(aligned, cfg)
    warnings_out.extend(normalized.warnings)
    stage("Background correction", "ok", "Robust sky background removed from both epochs.", diagnostics=normalized.diagnostics)

    # 5. Difference
    units = product1.bunit or "IMAGE units"
    difference = difference_engine.compute_difference(normalized, cfg, units=units)
    warnings_out.extend(difference.warnings)
    stage(
        "Image comparison",
        "ok",
        f"D(x, y) = I₂(x, y) − I₁(x, y); global noise {difference.global_noise:.4g} {units}.",
        diagnostics=difference.diagnostics,
    )

    # 6. Detection
    sources, detection_diagnostics = source_detection.detect_residual_sources(aligned, normalized, difference, cfg)
    stage(
        "Source detection",
        "ok" if detection_diagnostics.get("status") == "ok" else "warn",
        f"{len(sources)} residual source(s) above {cfg['detection_sigma']}σ using Photutils DAOStarFinder.",
        diagnostics=detection_diagnostics,
    )

    run_id = provenance.deterministic_run_id(
        [
            {"id": record1.get("id"), "sha256": product1.sha256, "access_url": record1.get("cutout_url"), "origin": record1.get("origin")},
            {"id": record2.get("id"), "sha256": product2.sha256, "access_url": record2.get("cutout_url"), "origin": record2.get("origin")},
        ],
        cfg,
    )

    # 7. Candidate measurement
    data_mode = record1.get("data_mode") or product1.data_mode
    candidates, candidate_summary = candidate_analysis.analyse_candidates(
        aligned, normalized, difference, sources, product1, product2, record1, record2, cfg, run_id, data_mode
    )
    stage(
        "Candidate measurement",
        "ok",
        f"{len(candidates)} candidate record(s) measured; association radius {candidate_summary['association_radius_px']:.1f} px.",
        counts=candidate_summary["counts"],
    )

    # 8. Catalogue cross-match
    crossmatch: dict[str, dict] = {}
    catalogue_diagnostics: dict = {"status": "disabled", "detail": "Cross-match disabled in processing configuration."}
    if cfg.get("crossmatch_enabled") and candidates:
        centre = product1.center_sky or (record1.get("ra_deg"), record1.get("dec_deg"))
        if centre[0] is not None and centre[1] is not None:
            try:
                rows, catalogue_diagnostics = catalogue_matching.fetch_catalogue(
                    catalogue_key, float(centre[0]), float(centre[1]), float(cfg["catalogue_cone_radius_deg"])
                )
                crossmatch = catalogue_matching.cross_match_candidates(
                    candidates, rows, radius_arcsec=float(cfg["crossmatch_radius_arcsec"])
                )
                matched = sum(1 for value in crossmatch.values() if value["status"] == "matched")
                stage("Catalogue cross-match", "ok", f"{matched}/{len(candidates)} candidate(s) matched in {catalogue_diagnostics['catalogue']}.")
            except Exception as exc:
                catalogue_diagnostics = {"status": "failed", "detail": str(exc), "catalogue": catalogue_key}
                crossmatch = {
                    candidate.candidate_id: {"status": "unavailable", "matches": [], "detail": f"Catalogue service failed: {exc}"}
                    for candidate in candidates
                }
                stage("Catalogue cross-match", "warn", f"Catalogue cross-match unavailable: {exc}")
        else:
            catalogue_diagnostics = {"status": "skipped", "detail": "Sky centre unavailable for a cone search."}
            stage("Catalogue cross-match", "warn", "Sky centre unavailable; cross-match skipped.")
    else:
        stage("Catalogue cross-match", "ok" if not candidates else "warn", catalogue_diagnostics["detail"])

    inputs = [
        {
            "role": "epoch 1",
            "id": record1.get("id"),
            "label": record1.get("label"),
            "origin": record1.get("origin"),
            "access_url": record1.get("cutout_url") or record1.get("full_fits_url") or record1.get("cached_path"),
            "sha256": product1.sha256,
            "date_utc": record1.get("date_utc"),
            "band": record1.get("band"),
            "central_wavelength_um": product1.center_wavelength_um,
            "central_bandwidth_um": product1.center_bandwidth_um,
            "data_mode": record1.get("data_mode"),
        },
        {
            "role": "epoch 2",
            "id": record2.get("id"),
            "label": record2.get("label"),
            "origin": record2.get("origin"),
            "access_url": record2.get("cutout_url") or record2.get("full_fits_url") or record2.get("cached_path"),
            "sha256": product2.sha256,
            "date_utc": record2.get("date_utc"),
            "band": record2.get("band"),
            "central_wavelength_um": product2.center_wavelength_um,
            "central_bandwidth_um": product2.center_bandwidth_um,
            "data_mode": record2.get("data_mode"),
        },
    ]
    provenance_record = provenance.build_provenance(
        inputs=inputs,
        processing_config=cfg,
        stages=stages,
        warnings=warnings_out + [w for c in candidates for w in c.warnings],
        data_mode=data_mode,
        run_id=run_id,
        result_summary={
            "candidates": len(candidates),
            "classification_counts": candidate_summary["counts"],
            "overlap_fraction": aligned.overlap_fraction,
            "valid_pixel_fraction": float(np.mean(aligned.valid_mask)),
            "global_difference_noise": difference.global_noise,
            "units": units,
            "disclaimer": "Candidate records are screening results from two observations, not discoveries.",
        },
    )
    provenance.write_json(config.OUTPUT_DIR / "runs" / f"{run_id}.json", provenance_record)
    if candidates:
        export_rows = []
        for candidate in candidates:
            values = _candidate_row(candidate, crossmatch.get(candidate.candidate_id))
            export_rows.append(values)
        _write_candidates_csv(export_rows, run_id, record1, record2)
    return {
        "status": "ok",
        "run_id": run_id,
        "product1": product1,
        "product2": product2,
        "compatibility": compatibility,
        "aligned": aligned,
        "normalized": normalized,
        "difference": difference,
        "detection": detection_diagnostics,
        "candidates": candidates,
        "candidate_summary": candidate_summary,
        "crossmatch": crossmatch,
        "catalogue_diagnostics": catalogue_diagnostics,
        "provenance": provenance_record,
        "stages": stages,
        "warnings": warnings_out,
        "data_mode": data_mode,
        "units": units,
    }


def _candidate_row(candidate, crossmatch: dict | None) -> dict:
    match = (crossmatch or {}).get("matches") or []
    first = match[0] if match else {}
    return {
        "candidate_id": candidate.candidate_id,
        "classification": candidate.classification,
        "ra_deg": candidate.ra_deg,
        "dec_deg": candidate.dec_deg,
        "pixel_x": candidate.pixel_x,
        "pixel_y": candidate.pixel_y,
        "peak_snr": candidate.snr,
        "flux_epoch1_jy": candidate.flux1_jy,
        "flux_epoch2_jy": candidate.flux2_jy,
        "delta_flux_jy": candidate.delta_flux_jy,
        "surface_brightness_difference": candidate.surface_brightness_diff,
        "local_noise": candidate.local_noise,
        "angular_displacement_arcsec": candidate.metadata.get("angular_displacement_arcsec"),
        "position_angle_deg": candidate.metadata.get("position_angle_deg"),
        "delta_time_days": candidate.metadata.get("delta_time_days"),
        "epoch1_date_utc": candidate.metadata.get("epoch1", {}).get("date_utc"),
        "epoch2_date_utc": candidate.metadata.get("epoch2", {}).get("date_utc"),
        "epoch1_obs_id": candidate.metadata.get("epoch1", {}).get("id"),
        "epoch2_obs_id": candidate.metadata.get("epoch2", {}).get("id"),
        "crossmatch_status": (crossmatch or {}).get("status", "not attempted"),
        "crossmatch_source": first.get("designation", ""),
        "crossmatch_separation_arcsec": first.get("separation_arcsec", ""),
        "data_mode": candidate.data_mode,
        "provenance_run_id": candidate.provenance_run_id,
        "warnings": " | ".join(candidate.warnings),
        "reasons": " | ".join(candidate.reasons),
    }


def _write_candidates_csv(rows: list[dict], run_id: str, record1: dict, record2: dict) -> str:
    import csv

    directory = config.OUTPUT_DIR / "candidates"
    directory.mkdir(parents=True, exist_ok=True)
    path = directory / f"candidates_{run_id}.csv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)
        handle.write(
            f"# provenance_run_id={run_id} epoch1={record1.get('id')} epoch2={record2.get('id')} "
            f"data_mode={record1.get('data_mode')} disclaimer=candidate screening only, not a discovery\n"
        )
    return str(path)


def candidates_to_csv(rows: list[dict]) -> str:
    import csv
    import io

    buffer = io.StringIO()
    if not rows:
        return "candidate_id,classification\n"
    writer = csv.DictWriter(buffer, fieldnames=list(rows[0].keys()))
    writer.writeheader()
    writer.writerows(rows)
    return buffer.getvalue()
