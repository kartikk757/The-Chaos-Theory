"""Tests for documented SPHEREx ePSF zone selection and diagnostics."""
from __future__ import annotations

from types import SimpleNamespace

import numpy as np
from astropy.io import fits
import pytest
from photutils.psf import ImagePSF

from modules import epsf_processor


def _epsf_product(path, extensions=None):
    primary = fits.PrimaryHDU()
    primary.header["VERSION"] = "7.0.5"
    image = fits.ImageHDU(data=np.zeros((64, 64), dtype=np.float32), name="IMAGE")
    image.header["CRPIX1A"] = 1.0
    image.header["CRPIX2A"] = 1.0
    hdus = [primary, image, *(extensions or [])]
    fits.HDUList(hdus).writeto(path, overwrite=True)
    return SimpleNamespace(path=str(path), metadata={"crop_origin_pixel": [10, 20]}, center_wavelength_um=1.42)


def test_selects_nearest_epsf_zone_in_parent_detector_coordinates(tmp_path):
    model = np.ones((33, 33), dtype=np.float32)
    columns = [
        fits.Column(name="BINX", format="J", array=[0, 1]),
        fits.Column(name="BINY", format="J", array=[0, 1]),
        fits.Column(name="XCENTER", format="E", array=[15.0, 50.0]),
        fits.Column(name="YCENTER", format="E", array=[25.0, 50.0]),
        fits.Column(name="NSTAR", format="J", array=[20, 30]),
        fits.Column(name="NEFF_MEAN", format="E", array=[4.2, 5.1]),
        fits.Column(name="CWAVE", format="E", unit="um", array=[1.42, 1.50]),
        fits.Column(name="CBAND", format="E", unit="um", array=[0.03, 0.04]),
        fits.Column(name="EPSF", format="1089E", dim="(33,33)", array=[model, model * 2]),
    ]
    epsf = fits.BinTableHDU.from_columns(columns, name="EPSF")
    epsf.header["OVSMPX"] = 5
    epsf.header["OVSMPY"] = 5
    epsf.header["DETECTOR"] = "2"
    product = _epsf_product(tmp_path / "epsf.fits", [epsf])

    result = epsf_processor.epsf_for_position(product, 5.0, 5.0)

    assert result["status"] == "ok"
    assert result["source_pixel_original_cutout"] == [15.0, 25.0]
    assert result["parent_detector_pixel"] == [15.0, 25.0]
    assert result["zone_id"] == [0, 0]
    assert result["epsf_shape_oversampled"] == [33, 33]
    assert np.sum(result["epsf_normalized"]) == pytest.approx(1.0)
    assert result["epsf_effective_area"] == pytest.approx(4.2)


def test_legacy_psf_is_not_misreported_as_epsf(tmp_path):
    psf = fits.ImageHDU(data=np.ones((4, 33, 33), dtype=np.float32), name="PSF")
    product = _epsf_product(tmp_path / "legacy.fits", [psf])
    result = epsf_processor.epsf_for_position(product, 2.0, 3.0)
    assert result["status"] == "unsupported_legacy_psf"
    assert "does not" in result["detail"]


def test_photutils_epsf_fit_recovers_injected_subpixel_source(tmp_path):
    oversampling = 5
    grid = (np.arange(33, dtype=float) - 16.0) / oversampling
    yy_model, xx_model = np.meshgrid(grid, grid, indexing="ij")
    epsf_array = np.exp(-0.5 * ((xx_model / 0.58) ** 2 + (yy_model / 0.58) ** 2))
    epsf_array *= (oversampling**2) / float(epsf_array.sum())
    columns = [
        fits.Column(name="BINX", format="J", array=[0]),
        fits.Column(name="BINY", format="J", array=[0]),
        fits.Column(name="XCENTER", format="E", array=[40.0]),
        fits.Column(name="YCENTER", format="E", array=[52.0]),
        fits.Column(name="NSTAR", format="J", array=[30]),
        fits.Column(name="NEFF_MEAN", format="E", array=[4.2]),
        fits.Column(name="CWAVE", format="E", unit="um", array=[1.42]),
        fits.Column(name="CBAND", format="E", unit="um", array=[0.04]),
        fits.Column(name="EPSF", format="1089E", dim="(33,33)", array=[epsf_array]),
    ]
    epsf = fits.BinTableHDU.from_columns(columns, name="EPSF")
    epsf.header["OVSMPX"] = oversampling
    epsf.header["OVSMPY"] = oversampling
    epsf.header["DETECTOR"] = "2"
    path = tmp_path / "fitted_epsf.fits"
    _epsf_product(path, [epsf])

    true_x, true_y, true_flux = 30.35, 31.65, 500.0
    image_y, image_x = np.indices((64, 64), dtype=float)
    epsf_model = ImagePSF(
        epsf_array,
        oversampling=(oversampling, oversampling),
        flux=true_flux,
        x_0=true_x,
        y_0=true_y,
    )
    rng = np.random.default_rng(42)
    science = 100.0 + epsf_model(image_x, image_y)
    science += rng.normal(0.0, 0.02, size=science.shape)
    product = SimpleNamespace(
        path=str(path),
        metadata={"crop_origin_pixel": [10, 20]},
        center_wavelength_um=1.42,
        science=science,
        variance=np.full(science.shape, 0.02**2, dtype=float),
        bad_flag_mask=np.zeros(science.shape, dtype=bool),
        header={"DETECTOR": "2"},
        bunit="MJy / sr",
        wcs=None,
        sha256="a" * 64,
    )

    result = epsf_processor.fit_source(product, true_x, true_y)

    assert result["status"] == "ok"
    assert result["flux_fit_native"] == pytest.approx(true_flux, rel=0.03)
    assert result["x_fit"] == pytest.approx(true_x, abs=0.08)
    assert result["y_fit"] == pytest.approx(true_y, abs=0.08)
    assert result["epsf_oversampling_xy"] == [5, 5]
    assert result["epsf_normalization_ratio"] == pytest.approx(1.0, rel=0.02)
    assert "PSFPhotometry" in result["interpretation"]

