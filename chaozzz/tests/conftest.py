"""Shared fixtures.

Synthetic products are used only where a controlled truth model is required.
They are labelled DEMONSTRATION DATA in their headers and are never presented as
SPHEREx observations.
"""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parent.parent
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

import config  # noqa: E402
from modules import data_loader  # noqa: E402


def _load_generator():
    path = ROOT / "scripts" / "generate_demo_data.py"
    spec = importlib.util.spec_from_file_location("ct_generate_demo_data", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="session")
def demo_pair(tmp_path_factory):
    """Write the synthetic two-epoch pair into a temporary directory."""
    generator = _load_generator()
    directory = tmp_path_factory.mktemp("demo")
    records = []
    for epoch in (1, 2):
        hdul, truth = generator.build_epoch(epoch)
        path = directory / f"demo_epoch{epoch}.fits"
        hdul.writeto(path, overwrite=True)
        records.append({
            "id": f"TEST-DEMO-{epoch}",
            "label": f"synthetic epoch {epoch}",
            "origin": "pytest synthetic fixture (not a SPHEREx observation)",
            "data_mode": config.DATA_MODES["DEMO"],
            "ra_deg": generator.CENTER_RA,
            "dec_deg": generator.CENTER_DEC,
            "date_utc": truth["date_utc"],
            "mjd": truth["mjd"],
            "band": "SPHEREx-D2 (synthetic stand-in)",
            "release": "demo",
            "obs_id": f"TEST-DEMO-{epoch}",
            "cutout_url": "",
            "full_fits_url": "",
            "cached_path": str(path),
            "calibration_level": None,
            "quality_flag": "synthetic",
            "provenance_version": "DEMO-1.0",
            "product_type": "synthetic image",
            "content_type": "image/fits",
            "plane_id": "",
            "is_real_spherex": False,
            "is_demo": True,
            "central_wavelength_um": None,
            "central_bandwidth_um": None,
            "coordinates_source": "synthetic WCS image centre",
            "metadata": {"synthetic": True},
            "warnings": ["DEMONSTRATION DATA: synthetic product, not a SPHEREx observation."],
        })
    return records, generator


@pytest.fixture(scope="session")
def demo_products(demo_pair):
    records, _ = demo_pair
    return [data_loader.load_observation_product(record) for record in records]


@pytest.fixture()
def toy_fits(tmp_path):
    """A minimal but structurally valid single-extension FITS file."""
    from astropy.io import fits

    data = np.random.default_rng(7).normal(5.0, 0.2, (32, 32)).astype(np.float32)
    header = fits.Header()
    header["CTYPE1"] = "RA---TAN"
    header["CTYPE2"] = "DEC--TAN"
    header["CRVAL1"] = 150.0
    header["CRVAL2"] = 2.0
    header["CRPIX1"] = 16.0
    header["CRPIX2"] = 16.0
    header["CDELT1"] = -0.001
    header["CDELT2"] = 0.001
    header["BUNIT"] = "MJy / sr"
    header["DATE-OBS"] = "2026-03-01T00:00:00"
    path = tmp_path / "toy.fits"
    fits.PrimaryHDU(data=data, header=header).writeto(path)
    return path


class _FakeUpload:
    def __init__(self, name: str, payload: bytes):
        self.name = name
        self._payload = payload

    def getbuffer(self) -> bytes:
        return self._payload


@pytest.fixture()
def fake_upload():
    return _FakeUpload
