"""Residual detection, dipole merging and classification guards."""
from __future__ import annotations

import numpy as np
import pytest

import config
from modules import candidate_analysis, difference_engine, image_normalization, source_detection, wcs_alignment
from modules.models import AlignedPair


def _pair_from_arrays(image1, image2, variance=1e-4, scale=6.14):
    shape = image1.shape
    return AlignedPair(
        target_wcs=None,
        shape=shape,
        science1=image1,
        science2=image2,
        variance1=np.full(shape, variance),
        variance2=np.full(shape, variance),
        bad_mask1=np.zeros(shape, dtype=bool),
        bad_mask2=np.zeros(shape, dtype=bool),
        source_mask1=np.zeros(shape, dtype=bool),
        source_mask2=np.zeros(shape, dtype=bool),
        coverage1=np.ones(shape),
        coverage2=np.ones(shape),
        valid_mask=np.ones(shape, dtype=bool),
        overlap_fraction=1.0,
        pixel_scale_arcsec=scale,
        diagnostics={},
        warnings=[],
    )


def _difference_for(aligned, cfg):
    normalized = image_normalization.normalize_pair(aligned, cfg)
    return normalized, difference_engine.compute_difference(normalized, cfg)


def test_detects_a_known_significant_residual():
    rng = np.random.default_rng(3)
    image1 = rng.normal(0.5, 0.01, (64, 64))
    image2 = image1.copy()
    image2[40, 22] += 0.08  # ~12 sigma with sigma = sqrt(2 * 1e-4) = 0.014
    aligned = _pair_from_arrays(image1, image2)
    cfg = dict(config.PROCESSING_DEFAULTS)
    normalized, difference = _difference_for(aligned, cfg)
    sources, diagnostics = source_detection.detect_residual_sources(aligned, normalized, difference, cfg)
    assert diagnostics["status"] == "ok"
    assert sources, diagnostics
    nearest = min(sources, key=lambda item: np.hypot(item["xcentroid"] - 22, item["ycentroid"] - 40))
    assert np.hypot(nearest["xcentroid"] - 22, nearest["ycentroid"] - 40) <= 1.5
    assert nearest["sign"] == 1
    assert nearest["peak_snr"] > cfg["detection_sigma"]


def test_pure_noise_yields_no_detections():
    rng = np.random.default_rng(5)
    image1 = rng.normal(0.5, 0.01, (64, 64))
    image2 = rng.normal(0.5, 0.01, (64, 64))
    aligned = _pair_from_arrays(image1, image2)
    cfg = dict(config.PROCESSING_DEFAULTS)
    cfg["detection_sigma"] = 6.0
    normalized, difference = _difference_for(aligned, cfg)
    sources, _ = source_detection.detect_residual_sources(aligned, normalized, difference, cfg)
    assert len(sources) <= 1  # a random 6-sigma fluctuation is tolerated but not a flood


def test_dipole_is_merged_into_one_candidate():
    rng = np.random.default_rng(9)
    image1 = rng.normal(0.5, 0.005, (64, 64))
    image2 = image1.copy()
    image1[30, 20] += 0.08
    image2[30, 22] += 0.08
    aligned = _pair_from_arrays(image1, image2, variance=1e-5)
    cfg = dict(config.PROCESSING_DEFAULTS)
    normalized, difference = _difference_for(aligned, cfg)
    sources, diagnostics = source_detection.detect_residual_sources(aligned, normalized, difference, cfg)
    positives = [item for item in sources if item["sign"] > 0]
    assert positives, sources
    merged = [item for item in sources if item["dipole_partner"] is not None]
    assert merged, "an opposite-sign neighbour within the merge radius must be attached as a dipole partner"
    assert diagnostics["positive_peaks"] >= 1


def test_threshold_is_respected_by_increasing_sigma():
    rng = np.random.default_rng(13)
    image1 = rng.normal(0.5, 0.01, (64, 64))
    image2 = image1.copy()
    image2[10, 10] += 0.05
    aligned = _pair_from_arrays(image1, image2)
    cfg = dict(config.PROCESSING_DEFAULTS)
    normalized, difference = _difference_for(aligned, cfg)
    low_cfg = dict(cfg, detection_sigma=3.0)
    high_cfg = dict(cfg, detection_sigma=12.0)
    low_sources, _ = source_detection.detect_residual_sources(aligned, normalized, difference, low_cfg)
    high_sources, _ = source_detection.detect_residual_sources(aligned, normalized, difference, high_cfg)
    assert len(low_sources) >= len(high_sources)


def test_classification_never_uses_discovery_language(demo_pair, demo_products):
    records, _ = demo_pair
    cfg = dict(config.PROCESSING_DEFAULTS)
    aligned = wcs_alignment.align_pair(demo_products[0], demo_products[1], cfg)
    normalized, difference = _difference_for(aligned, cfg)
    sources, _ = source_detection.detect_residual_sources(aligned, normalized, difference, cfg)
    candidates, summary = candidate_analysis.analyse_candidates(
        aligned, normalized, difference, sources, demo_products[0], demo_products[1],
        records[0], records[1], cfg, "CT-TESTRUN", config.DATA_MODES["DEMO"],
    )
    assert candidates, "the synthetic truth model contains changes that must be detected"
    forbidden = ("planet", "asteroid", "discovery", "confirmed")
    for candidate in candidates:
        assert candidate.classification in candidate_analysis.CATEGORY_DEFINITIONS
        text = " ".join(candidate.reasons + candidate.warnings).lower()
        for word in forbidden:
            assert word not in text
    assert summary["note"].lower().startswith("classifications are screening labels")
    moving = [c for c in candidates if c.classification == candidate_analysis.CLASS_MOVING]
    variable = [c for c in candidates if c.classification == candidate_analysis.CLASS_VARIABLE]
    assert moving, [c.classification for c in candidates]
    assert variable, [c.classification for c in candidates]


def test_candidate_measures_units_and_displacement(demo_pair, demo_products):
    records, _ = demo_pair
    cfg = dict(config.PROCESSING_DEFAULTS)
    aligned = wcs_alignment.align_pair(demo_products[0], demo_products[1], cfg)
    normalized, difference = _difference_for(aligned, cfg)
    sources, _ = source_detection.detect_residual_sources(aligned, normalized, difference, cfg)
    candidates, _ = candidate_analysis.analyse_candidates(
        aligned, normalized, difference, sources, demo_products[0], demo_products[1],
        records[0], records[1], cfg, "CT-TESTRUN2", config.DATA_MODES["DEMO"],
    )
    mover = next(c for c in candidates if c.classification == candidate_analysis.CLASS_MOVING)
    assert mover.metadata["angular_displacement_arcsec"] is not None
    assert mover.metadata["angular_displacement_arcsec"] > mover.quality["alignment_tolerance_arcsec"]
    assert mover.metadata["units"]["angular_displacement"] == "arcsec"
    assert mover.metadata["units"]["flux"] == "Jy"
    assert mover.metadata["delta_time_days"] == pytest.approx(abs(records[1]["mjd"] - records[0]["mjd"]))
    assert mover.quality["sign"] in (1, -1)
