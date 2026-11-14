"""Tests for truthful WCS axes and image-marker coordinates."""
import numpy as np
from astropy.wcs import WCS

from modules import visualization


def test_wcs_image_figure_uses_sky_coordinates_and_transforms_pixel_marker():
    wcs = WCS(naxis=2)
    wcs.wcs.crpix = [5.0, 5.0]
    wcs.wcs.crval = [10.0, 20.0]
    wcs.wcs.cdelt = [-0.001, 0.001]
    wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    image = np.zeros((10, 10))
    fig = visualization.image_figure(image, "test", wcs=wcs, marker={"x": 2.0, "y": 3.0})
    sky = wcs.pixel_to_world(2.0, 3.0).icrs
    assert np.allclose(fig.data[0].x, wcs.pixel_to_world(np.arange(10), np.full(10, 4.5)).icrs.ra.deg)
    assert fig.data[1].x[0] == sky.ra.deg
    assert fig.data[1].y[0] == sky.dec.deg
