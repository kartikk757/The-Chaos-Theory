"""Offline tests for HiPS properties and NESTED tile URL construction."""
import sys
import types

import pytest

from modules import hips_viewer


def test_properties_parser_and_nested_cube_url(monkeypatch):
    module = types.ModuleType("hpgeom")
    module.angle_to_pixel = lambda nside, ra, dec, nest, lonlat: 445
    monkeypatch.setitem(sys.modules, "hpgeom", module)
    properties = hips_viewer.parse_properties("hips_order=5\nhips_order_min=0\n# note\nhips_data_url=https://example.test/hips5_cube/version\n")
    result = hips_viewer.hips_cube_url(60.81, 51.34, properties=properties)
    assert result["ipix_nested"] == 445
    assert result["order"] == 5
    assert result["url"].endswith("/Norder5/Dir0/Npix445_cube.fits")


def test_invalid_coordinates_are_refused():
    with pytest.raises(ValueError):
        hips_viewer.healpix_nested_index(361, 0, 5)


def test_aladin_embed_uses_clamped_target_and_official_hi_ps_url():
    html = hips_viewer.aladin_lite_html(200.0, 10.0, 0.4)
    assert "200.0000000 10.0000000" in html
    assert "spherex_allsky_cube" in html
    assert "AladinLite/api/v3/latest/aladin.js" in html
