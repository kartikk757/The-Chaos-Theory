"""Background normalization, difference imaging and noise propagation."""
from __future__ import annotations

import numpy as np
import pytest

import config
from modules import difference_engine, image_normalization, wcs_alignment


def _aligned_with_offsets(offset1: float = 0.8, offset2: float = 0.35):
    """Build a controlled AlignedPair: same sky, different backgrounds, one brighter source."""
    from modules.models import AlignedPair

    shape = (48, 48)
    rng = np.random.default_rng(11)
    base = rng.normal(0.0, 0.02, shape) + 0.5
    image1 = base + offset1
    image2 = base + offset2
    image2[20:22, 20:22] += 0.4  # a real brightness increase in epoch 2
    variance = np.full(shape, 0.02 ** 2)
    valid = np.ones(shape, dtype=bool)
    return AlignedPair(
        target_wcs=None,
        shape=shape,
        science1=image1,
        science2=image2,
        variance1=variance.copy(),
        variance2=variance.copy(),
        bad_mask1=np.zeros(shape, dtype=bool),
        bad_mask2=np.zeros(shape, dtype=bool),
        source_mask1=np.zeros(shape, dtype=bool),
        source_mask2=np.zeros(shape, dtype=bool),
        coverage1=np.ones(shape),
        coverage2=np.ones(shape),
        valid_mask=valid,
        overlap_fraction=1.0,
        pixel_scale_arcsec=6.14,
        diagnostics={},
        warnings=[],
    )


def test_background_offsets_are_removed_and_real_change_preserved():
    aligned = _aligned_with_offsets()
    cfg = dict(config.PROCESSING_DEFAULTS)
    normalized = image_normalization.normalize_pair(aligned, cfg)
    background_level = np.median(normalized.image1)
    assert abs(background_level) < 0.02
    assert np.median(normalized.image2) == pytest.approx(np.median(normalized.image1), abs=0.02)
    difference = difference_engine.compute_difference(normalized, cfg)
    patch = difference.difference[20:22, 20:22]
    assert np.nanmean(patch) > 0.25
    assert difference.global_noise > 0


def test_difference_obeys_equation_and_uses_variance():
    aligned = _aligned_with_offsets()
    cfg = dict(config.PROCESSING_DEFAULTS)
    normalized = image_normalization.normalize_pair(aligned, cfg)
    difference = difference_engine.compute_difference(normalized, cfg)
    manual = normalized.image2 - normalized.image1
    mask = difference.valid_mask
    assert np.allclose(difference.difference[mask], manual[mask], atol=1e-12)
    assert "VARIANCE" in difference.diagnostics["noise_model"]
    expected_sigma = np.sqrt(normalized.variance1 + normalized.variance2)
    assert np.nanmedian(difference.sigma_map) == pytest.approx(np.median(expected_sigma), rel=1e-6)
    assert np.nanmax(np.abs(difference.snr)) > 5


def test_invalid_pixels_become_nan_not_zero():
    aligned = _aligned_with_offsets()
    aligned.valid_mask[5:10, 5:10] = False
    cfg = dict(config.PROCESSING_DEFAULTS)
    normalized = image_normalization.normalize_pair(aligned, cfg)
    difference = difference_engine.compute_difference(normalized, cfg)
    assert np.isnan(difference.difference[5:10, 5:10]).all()
    assert difference.diagnostics["difference_pixels_nan"] >= 25


def test_noise_falls_back_when_variance_missing():
    aligned = _aligned_with_offsets()
    aligned.variance1 = None
    aligned.variance2 = None
    cfg = dict(config.PROCESSING_DEFAULTS)
    normalized = image_normalization.normalize_pair(aligned, cfg)
    difference = difference_engine.compute_difference(normalized, cfg)
    assert np.isfinite(difference.global_noise) and difference.global_noise > 0
    assert "robust" in difference.diagnostics["noise_model"].lower()
    assert any("VARIANCE" in warning for warning in difference.warnings)


def test_difference_declares_limitations():
    aligned = _aligned_with_offsets()
    cfg = dict(config.PROCESSING_DEFAULTS)
    normalized = image_normalization.normalize_pair(aligned, cfg)
    difference = difference_engine.compute_difference(normalized, cfg)
    assert any("does not establish" in warning for warning in difference.warnings)


def test_normalization_reports_background_model():
    aligned = _aligned_with_offsets()
    cfg = dict(config.PROCESSING_DEFAULTS)
    normalized = image_normalization.normalize_pair(aligned, cfg)
    diagnostics = normalized.diagnostics
    assert diagnostics["sky_pixels_used"] > 0
    assert "background_epoch1" in diagnostics and "background_epoch2" in diagnostics
    assert "units_after_normalization" in diagnostics
