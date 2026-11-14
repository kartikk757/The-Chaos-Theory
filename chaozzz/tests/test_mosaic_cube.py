"""Offline tests for FITS cube validation and background-subtracted extraction."""
import io

import numpy as np
from astropy.io import fits

from modules import mosaic_cube


def _cube_bytes():
    data = np.ones((5, 32, 32), dtype=np.float32)
    for channel in range(data.shape[0]):
        data[channel, 16, 16] += 10.0 + channel
    header = fits.Header()
    header["CTYPE1"] = "RA---TAN"
    header["CTYPE2"] = "DEC--TAN"
    header["CTYPE3"] = "WAVE"
    header["CUNIT1"] = "deg"
    header["CUNIT2"] = "deg"
    header["CUNIT3"] = "um"
    header["CRPIX1"] = 16.5
    header["CRPIX2"] = 16.5
    header["CRPIX3"] = 1.0
    header["CRVAL1"] = 10.0
    header["CRVAL2"] = 20.0
    header["CRVAL3"] = 1.0
    header["CDELT1"] = -0.001
    header["CDELT2"] = 0.001
    header["CDELT3"] = 0.1
    header["BUNIT"] = "MJy / sr"
    image = fits.ImageHDU(data=data, header=header, name="IMAGE")
    channels = fits.BinTableHDU.from_columns([
        fits.Column(name="WAVELENGTH", format="E", unit="um", array=np.linspace(1.0, 1.4, 5)),
        fits.Column(name="BANDWIDTH", format="E", unit="um", array=np.full(5, 0.1)),
    ], name="SPECTRAL_CHANNELS")
    buffer = io.BytesIO()
    fits.HDUList([fits.PrimaryHDU(), image, channels]).writeto(buffer)
    return buffer.getvalue()


def test_mosaic_cube_validates_wavelengths_and_extracts_background_subtracted_aperture():
    cube = mosaic_cube.load_cube(_cube_bytes(), "controlled.fits")
    result = mosaic_cube.extract_aperture_spectrum(cube, 10.0, 20.0, aperture_radius_px=2.0, annulus_inner_px=5.0, annulus_outer_px=8.0)
    assert cube.metadata["wavelength_source"].startswith("SPECTRAL_CHANNELS")
    assert result["status"] == "ok"
    np.testing.assert_allclose(result["wavelength_um"], np.linspace(1.0, 1.4, 5), rtol=1e-6, atol=0.0)
    assert np.all(result["flux"] > 0)
    assert "mJy" in result["flux_unit"]
    assert len(result["limitations"]) >= 3


def test_cube_without_celestial_wcs_is_refused():
    header = fits.Header()
    data = np.ones((3, 8, 8), dtype=np.float32)
    buffer = io.BytesIO()
    fits.PrimaryHDU(data=data, header=header).writeto(buffer)
    try:
        mosaic_cube.load_cube(buffer.getvalue(), "no-wcs.fits")
    except ValueError as error:
        assert "celestial WCS" in str(error)
    else:
        raise AssertionError("Expected cube without WCS to be refused")

