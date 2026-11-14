"""Compatibility gating and WCS reprojection."""
from __future__ import annotations

import copy

import numpy as np
import pytest

import config
from modules import wcs_alignment


def test_compatible_pair_passes_all_checks(demo_pair, demo_products):
    records, _ = demo_pair
    cfg = dict(config.PROCESSING_DEFAULTS)
    result = wcs_alignment.check_observation_compatibility(records[0], records[1], demo_products[0], demo_products[1], cfg)
    assert result["compatible"] is True, result["reasons"]
    names = {check["name"] for check in result["checks"]}
    assert {"Detector / LVF band", "Data mode", "Celestial WCS", "WCS-WAVE spectral compatibility"} <= names


def test_different_band_is_rejected(demo_pair, demo_products):
    records, _ = demo_pair
    changed = copy.deepcopy(records[1])
    changed["band"] = "SPHEREx-D3"
    result = wcs_alignment.check_observation_compatibility(
        records[0], changed, demo_products[0], demo_products[1], dict(config.PROCESSING_DEFAULTS)
    )
    assert result["compatible"] is False
    assert any("band" in reason.lower() for reason in result["reasons"])


def test_mixed_data_modes_are_rejected(demo_pair, demo_products):
    records, _ = demo_pair
    changed = copy.deepcopy(records[1])
    changed["data_mode"] = config.DATA_MODES["LIVE"]
    result = wcs_alignment.check_observation_compatibility(
        records[0], changed, demo_products[0], demo_products[1], dict(config.PROCESSING_DEFAULTS)
    )
    assert result["compatible"] is False
    assert any("data mode" in reason.lower() for reason in result["reasons"])


def test_wavelength_mismatch_is_rejected(demo_pair, demo_products):
    records, _ = demo_pair
    product1 = copy.copy(demo_products[0])
    product2 = copy.copy(demo_products[1])
    product1.center_wavelength_um = 1.4200
    product1.center_bandwidth_um = 0.034
    product2.center_wavelength_um = 1.5200
    product2.center_bandwidth_um = 0.034
    result = wcs_alignment.check_observation_compatibility(records[0], records[1], product1, product2, dict(config.PROCESSING_DEFAULTS))
    assert result["compatible"] is False
    assert result["wavelength_delta_um"] == pytest.approx(0.1, abs=1e-9)
    assert any("WCS-WAVE" in reason for reason in result["reasons"])


def test_wavelength_within_tolerance_is_accepted(demo_pair, demo_products):
    records, _ = demo_pair
    product1 = copy.copy(demo_products[0])
    product2 = copy.copy(demo_products[1])
    product1.center_wavelength_um = 1.42020
    product1.center_bandwidth_um = 0.03422
    product2.center_wavelength_um = 1.42058
    product2.center_bandwidth_um = 0.03423
    result = wcs_alignment.check_observation_compatibility(records[0], records[1], product1, product2, dict(config.PROCESSING_DEFAULTS))
    assert result["compatible"] is True, result["reasons"]
    assert result["wavelength_tolerance_um"] == pytest.approx(0.25 * 0.03422, rel=1e-6)


def test_missing_wcs_is_rejected(demo_pair, demo_products):
    records, _ = demo_pair
    product2 = copy.copy(demo_products[1])
    product2.wcs = None
    result = wcs_alignment.check_observation_compatibility(records[0], records[1], demo_products[0], product2, dict(config.PROCESSING_DEFAULTS))
    assert result["compatible"] is False
    assert any("WCS" in reason for reason in result["reasons"])


def test_align_pair_reprojects_onto_common_grid(demo_products):
    cfg = dict(config.PROCESSING_DEFAULTS)
    aligned = wcs_alignment.align_pair(demo_products[0], demo_products[1], cfg)
    assert aligned.shape == demo_products[0].shape
    assert aligned.target_wcs is not None
    assert aligned.overlap_fraction > 0.9
    assert np.mean(aligned.valid_mask) > 0.5
    # NaN/NONFUNC patch must remain invalid after reprojection.
    assert not aligned.valid_mask[10, 84]
    # Non-finite values must never reach the science arrays used for arithmetic.
    assert np.isfinite(aligned.science1[aligned.valid_mask]).all()
    assert np.isfinite(aligned.science2[aligned.valid_mask]).all()


def test_align_pair_refuses_insufficient_overlap(demo_products):
    cfg = dict(config.PROCESSING_DEFAULTS)
    shifted = copy.copy(demo_products[1])
    shifted.wcs = copy.deepcopy(demo_products[1].wcs)
    shifted.wcs.wcs.crval = [demo_products[1].wcs.wcs.crval[0] + 5.0, demo_products[1].wcs.wcs.crval[1]]
    shifted.center_sky = (shifted.wcs.wcs.crval[0], shifted.wcs.wcs.crval[1])
    with pytest.raises(wcs_alignment.AlignmentError):
        wcs_alignment.align_pair(demo_products[0], shifted, cfg)


def test_pixel_scale_is_measured_from_wcs(demo_products):
    scale = demo_products[0].pixel_scale_arcsec
    assert scale is not None
    assert 5.0 < scale < 8.0


def test_parent_footprint_separation_does_not_replace_cutout_wcs_centres(demo_pair, demo_products):
    records, _ = demo_pair
    record1, record2 = copy.deepcopy(records)
    record1.update(ra_deg=5.0, dec_deg=0.0, coordinates_role="parent_observation_footprint_center")
    record2.update(ra_deg=15.0, dec_deg=0.0, coordinates_role="parent_observation_footprint_center")
    result = wcs_alignment.check_observation_compatibility(
        record1, record2, demo_products[0], demo_products[1], dict(config.PROCESSING_DEFAULTS)
    )

    assert result["compatible"] is True
    assert result["center_separation_arcsec"] < 1.0
    assert result["footprint_center_separation_arcsec"] > 30_000
