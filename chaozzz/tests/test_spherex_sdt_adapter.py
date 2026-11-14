"""Offline tests for the optional Source Discovery Tool adapter."""
import builtins
import hashlib
import sys
from types import ModuleType, SimpleNamespace

import numpy as np
from astropy import units as u
from astropy.io import fits
from astropy.table import Table

from modules import spherex_sdt_adapter


def test_independent_catalogue_comparison_reports_agreement_and_disagreement():
    result = spherex_sdt_adapter.compare_with_photutils(
        [
            {"x_pixel": 10.0, "y_pixel": 20.0},
            {"x_pixel": 40.0, "y_pixel": 40.0},
        ],
        [
            {"xcentroid": 10.5, "ycentroid": 20.1},
            {"xcentroid": 70.0, "ycentroid": 70.0},
        ],
        match_radius_px=2.0,
    )
    assert result["matched_pairs"] == 1
    assert result["photutils_only_count"] == 1
    assert result["sdt_only_count"] == 1
    assert "not a calibrated" in result["interpretation"]


def test_optional_sdt_spectral_photometry_uses_each_validated_epoch_position(tmp_path, monkeypatch):
    image_header = fits.Header()
    image_header["BUNIT"] = "MJy/sr"
    image_header["DETECTOR"] = 2
    image_header["MJD-OBS"] = 61234.5
    primary_header = fits.Header()
    primary_header["VERSION"] = "7.0.5"
    flags_header = fits.Header()
    flags_header["MSKN0000"] = "HOT"
    science = np.ones((59, 59), dtype=np.float32)
    flags = np.zeros((59, 59), dtype=np.uint32)
    variance = np.full((59, 59), 0.25, dtype=np.float32)
    variance[0, 0] = np.nan
    hdus = fits.HDUList([
        fits.PrimaryHDU(header=primary_header),
        fits.ImageHDU(science, header=image_header, name="IMAGE"),
        fits.ImageHDU(flags, header=flags_header, name="FLAGS"),
        fits.ImageHDU(variance, name="VARIANCE"),
        fits.ImageHDU(science.copy(), name="ZODI"),
        fits.BinTableHDU.from_columns([], name="EPSF"),
        fits.BinTableHDU.from_columns([], name="WCS-WAVE"),
    ])
    path = tmp_path / "sdt_epoch.fits"
    hdus.writeto(path)
    checksum = hashlib.sha256(path.read_bytes()).hexdigest()
    product = SimpleNamespace(product_id="epoch-1", path=str(path), sha256=checksum, is_real_spherex=True)

    parent = ModuleType("spherex_source_discovery_tool")
    parent.__path__ = []
    photometry_module = ModuleType("spherex_source_discovery_tool.aperture_photometry")
    calls = {}

    def fake_extract_photometry(image, positions, aperture_radii, **kwargs):
        calls["position"] = positions[0]
        calls["mask_has_nan"] = bool(np.isnan(image).any())
        calls["annulus"] = (kwargs["annulus_rad_in"], kwargs["annulus_rad_out"])
        return Table({"aper_sum_bkgsub_r3.0": [12.0]})

    def fake_extract_dqf(flag_mask, positions, aperture_radii):
        return Table({"aperture_sum_r3.0": [float(np.count_nonzero(flag_mask))]})

    photometry_module.fetch_sapm = lambda paths: ({}, {})
    photometry_module.sa_mjy_scale_map = lambda band, ny, nx, **kwargs: np.ones((ny, nx))
    photometry_module.extract_photometry = fake_extract_photometry
    photometry_module.extract_dqf = fake_extract_dqf
    photometry_module.get_AB = lambda flux_jy: 18.0
    monkeypatch.setitem(sys.modules, "spherex_source_discovery_tool", parent)
    monkeypatch.setitem(sys.modules, "spherex_source_discovery_tool.aperture_photometry", photometry_module)

    class FakeWCS:
        def __init__(self, *args, **kwargs):
            self.sip = "unused"

        def world_to_pixel(self, coordinate):
            return 29.0, 30.0

        def pixel_to_world(self, x, y):
            return 1.24 * u.um, 0.03 * u.um

    monkeypatch.setattr("astropy.wcs.WCS", FakeWCS)
    result = spherex_sdt_adapter.extract_epoch_spectrum(
        [product], [(151.2, 3.4)], source_id="candidate-1", sapm_s3=["test-map"]
    )

    assert result["status"] == "ok"
    row = result["measurements"][0]
    assert row["observation_id"] == "epoch-1"
    assert row["wavelength_um"] == 1.24
    assert row["bandwidth_um"] == 0.03
    assert row["product_version"] == "7.0.5"
    assert row["flux"] == 12.0
    assert row["flux_unit"] == "mJy"
    assert row["flux_error"] > 0
    assert row["known_bad_flags_found"] == ["HOT"]
    assert calls["position"] == (29.0, 30.0)
    assert calls["mask_has_nan"]
    assert calls["annulus"] == (7, 16)


def test_optional_sdt_spectral_photometry_falls_back_when_package_is_missing(monkeypatch):
    original_import = builtins.__import__

    def import_without_sdt(name, *args, **kwargs):
        if name.startswith("spherex_source_discovery_tool"):
            raise ModuleNotFoundError("optional SDT is not installed")
        return original_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", import_without_sdt)
    result = spherex_sdt_adapter.extract_epoch_spectrum(
        [object()], [(151.2, 3.4)], source_id="candidate-1"
    )
    assert result["status"] == "unavailable"
    assert result["measurements"] == []
    assert "optional SDT is not installed" in result["detail"]

