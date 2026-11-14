"""FITS loading, validation, FLAGS parsing and upload safety."""
from __future__ import annotations

import numpy as np
import pytest

import config
from modules import data_loader, fits_processor, provenance


def test_fits_product_has_required_extensions_and_units(demo_products):
    product = demo_products[0]
    assert product.shape == (96, 96)
    assert "IMAGE" in product.extension_names
    assert "FLAGS" in product.extension_names
    assert "VARIANCE" in product.extension_names
    assert "ZODI" in product.extension_names
    assert product.bunit.replace(" ", "").lower() == "mjy/sr"
    assert product.variance is not None and product.variance.shape == product.shape
    assert product.wcs is not None and product.wcs.has_celestial


def test_flag_definitions_are_read_from_header_and_source_is_informational(demo_products):
    product = demo_products[0]
    assert "SOURCE" in product.flag_definitions
    assert product.flag_definitions["SOURCE"]["bit_mask"] == 2097152
    assert "NONFUNC" in product.flag_definitions
    # The informational SOURCE bit must not mark every pixel as a defect.
    source_bit_only = (product.flags & 2097152) != 0
    assert source_bit_only.mean() > 0.9
    assert product.bad_flag_mask.mean() < 0.05
    assert np.array_equal(product.bad_flag_mask, product.flags.astype(np.int64) & 64 != 0)


def test_nan_pixels_are_preserved_and_flagged(demo_products):
    product = demo_products[0]
    assert np.isnan(product.image[10, 84])
    assert product.bad_flag_mask[10, 84]
    assert product.metadata["flag_diagnostics"]["available"] is True


def test_validation_records_state_every_check(demo_products):
    product = demo_products[0]
    names = {check["name"] for check in product.validation}
    assert {"2-D IMAGE array", "Celestial WCS", "FLAGS extension", "VARIANCE extension"} <= names
    assert all(check["status"] in {"pass", "warn", "fail"} for check in product.validation)
    assert "WCS-WAVE wavelength" in names
    # The synthetic product has no WCS-WAVE table, so the check must be a warning, not a pass.
    wavelength_check = next(check for check in product.validation if check["name"] == "WCS-WAVE wavelength")
    assert wavelength_check["status"] == "warn"


def test_checksum_matches_file_bytes(demo_pair):
    records, _ = demo_pair
    path = records[0]["cached_path"]
    product = data_loader.load_fits_product(path, records[0])
    assert product.sha256 == provenance.sha256_file(path)
    assert product.file_size_bytes > 0


def test_loader_rejects_non_fits_payload(fake_upload):
    with pytest.raises(data_loader.DataLoadError):
        data_loader.save_uploaded_fits(fake_upload("bad.fits", b"not a fits file at all"))


def test_loader_rejects_wrong_extension(fake_upload):
    with pytest.raises(data_loader.DataLoadError):
        data_loader.save_uploaded_fits(fake_upload("image.png", b"\x89PNG\r\n\x1a\n"))


def test_upload_accepts_a_valid_fits_file(fake_upload, toy_fits):
    payload = toy_fits.read_bytes()
    upload = data_loader.save_uploaded_fits(fake_upload("uploaded.fits", payload))
    assert upload["bytes"] == len(payload)
    assert upload["sha256"] == provenance.sha256_bytes(payload)
    product = data_loader.load_fits_product(upload["path"], {"id": "LOCAL-TEST", "data_mode": config.DATA_MODES["LOCAL"]})
    record = data_loader.make_local_observation(upload, product, 1)
    assert record["data_mode"] == config.DATA_MODES["LOCAL"]
    assert record["is_real_spherex"] is False
    assert record["date_utc"] == "2026-03-01T00:00:00"
    assert any("local" in warning.lower() or "WCS" in warning for warning in record["warnings"] + product.warnings)


def test_missing_file_raises_clear_error(tmp_path):
    with pytest.raises(data_loader.DataLoadError):
        data_loader.load_fits_product(tmp_path / "ghost.fits", {"id": "X"})


def test_flag_mask_keeps_informational_bits():
    flags = np.full((4, 4), 2097152, dtype=np.int64)
    definitions = {"SOURCE": {"bit_mask": 2097152}, "NONFUNC": {"bit_mask": 64}}
    bad, source, diagnostics = fits_processor.build_flag_masks(flags, definitions)
    assert bad.sum() == 0
    assert source.all()
    assert diagnostics["source_flag_fraction"] == 1.0
    flags[0, 0] = 64 | 2097152
    bad, _, _ = fits_processor.build_flag_masks(flags, definitions)
    assert bad[0, 0] and bad.sum() == 1


def test_wavelength_sampled_at_cutout_centre_not_archive_footprint(demo_products):
    product = demo_products[0]
    assert product.metadata["wavelength_sample_pixel"] == [47.5, 47.5]
    assert product.is_demo is True
    assert product.is_real_spherex is False
