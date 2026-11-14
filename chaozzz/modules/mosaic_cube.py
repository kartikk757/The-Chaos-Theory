"""Validated SPHEREx mosaic/HiPS cube loading and quick-look extraction."""
from __future__ import annotations

import io
from dataclasses import dataclass
from typing import Any

import astropy.units as u
import numpy as np
import streamlit as st
from astropy.coordinates import SkyCoord
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_area
from photutils.aperture import CircularAnnulus, CircularAperture, aperture_photometry


@dataclass
class SpectralCube:
    data: np.ndarray
    wavelengths_um: np.ndarray
    celestial_wcs: WCS
    spectral_numpy_axis: int
    header: dict[str, Any]
    source_name: str
    sha256: str
    bunit: str
    variance: np.ndarray | None
    flags: np.ndarray | None
    nhits: np.ndarray | None
    metadata: dict[str, Any]

    @property
    def spectral_first(self) -> np.ndarray:
        return np.moveaxis(self.data, self.spectral_numpy_axis, 0)


def _hash_bytes(payload: bytes) -> str:
    import hashlib
    return hashlib.sha256(payload).hexdigest()


def _wave_from_table(hdul, channel_count: int) -> np.ndarray | None:
    hdu = next((item for item in hdul if str(item.name).upper() == "SPECTRAL_CHANNELS"), None)
    if hdu is None or hdu.data is None:
        return None
    column_names = getattr(hdu.data, "names", None) or getattr(getattr(hdu.data, "dtype", None), "names", None)
    if not column_names:
        return None
    names = {name.upper(): name for name in column_names}
    name = names.get("WAVELENGTH") or names.get("CWAVE") or names.get("WAVE")
    if not name:
        return None
    try:
        values = np.asarray(hdu.data[name], dtype=float).squeeze()
    except (TypeError, ValueError):
        return None
    if values.ndim != 1 or values.size != channel_count or not np.all(np.isfinite(values)):
        return None
    unit = getattr(hdu.columns[name], "unit", None)
    if unit:
        try:
            values = (values * u.Unit(unit)).to_value(u.um)
        except Exception:
            return None
    # The official SPHEREx SPECTRAL_CHANNELS WAVELENGTH field is in microns.
    return values


def _spectral_wcs_axis(wcs: WCS, shape: tuple[int, ...]) -> tuple[np.ndarray, int]:
    types = wcs.get_axis_types()
    world_axis = next((index for index, info in enumerate(types) if info.get("coordinate_type") == "spectral"), None)
    if world_axis is None:
        raise ValueError("Cube WCS has no spectral world axis.")
    correlation = np.asarray(wcs.axis_correlation_matrix, dtype=bool)
    pixel_axes = np.where(correlation[world_axis])[0]
    if len(pixel_axes) != 1:
        raise ValueError("Spectral WCS is coupled to multiple pixel axes; extraction is not supported for this cube layout.")
    pixel_axis = int(pixel_axes[0])
    if wcs.pixel_n_dim != len(shape):
        raise ValueError("FITS cube dimensionality does not match its WCS axes.")
    numpy_axis = len(shape) - 1 - pixel_axis
    count = shape[numpy_axis]
    pixels = np.zeros((count, wcs.pixel_n_dim), dtype=float)
    pixels[:, pixel_axis] = np.arange(count, dtype=float)
    world = np.asarray(wcs.all_pix2world(pixels, 0), dtype=float)
    values = world[:, world_axis]
    source_unit = wcs.world_axis_units[world_axis]
    try:
        wavelengths = (values * u.Unit(source_unit)).to_value(u.um)
    except Exception as exc:
        raise ValueError(f"Could not convert spectral WCS unit {source_unit!r} to microns: {exc}") from exc
    if not np.all(np.isfinite(wavelengths)) or np.any(wavelengths <= 0):
        raise ValueError("Spectral WCS produced invalid wavelength values.")
    return np.asarray(wavelengths, dtype=float), numpy_axis


def load_cube(payload: bytes | bytearray, source_name: str = "uploaded.fits", max_pixels: int = 80_000_000) -> SpectralCube:
    """Load a 3-D FITS cube and validate its spectral and celestial axes."""
    raw = bytes(payload)
    if len(raw) < 2880:
        raise ValueError("The upload is too small to be a FITS file.")
    with fits.open(io.BytesIO(raw), memmap=False, do_not_scale_image_data=False) as hdul:
        image_hdu = next((item for item in hdul if str(item.name).upper() == "IMAGE" and item.data is not None), None)
        if image_hdu is None:
            image_hdu = next((item for item in hdul if isinstance(item, (fits.ImageHDU, fits.PrimaryHDU)) and item.data is not None), None)
        if image_hdu is None or image_hdu.data is None:
            raise ValueError("No image data extension was found.")
        data = np.asarray(image_hdu.data, dtype=float)
        if data.ndim != 3:
            raise ValueError(f"Expected a 3-D spectral cube; found {data.ndim}-D image data.")
        if int(np.prod(data.shape)) > int(max_pixels):
            raise ValueError(f"Cube has {int(np.prod(data.shape)):,} pixels, above the {int(max_pixels):,}-pixel processing limit.")
        wcs = WCS(image_hdu.header, fobj=hdul)
        if not wcs.has_celestial:
            raise ValueError("Cube has no valid celestial WCS.")
        table_waves = _wave_from_table(hdul, data.shape[0])
        if table_waves is not None and data.shape[0] in (table_waves.size,):
            spectral_axis = 0
            wavelengths = table_waves
            wavelength_source = "SPECTRAL_CHANNELS.WAVELENGTH (micron)"
        else:
            wavelengths, spectral_axis = _spectral_wcs_axis(wcs, data.shape)
            if wavelengths.size != data.shape[spectral_axis]:
                raise ValueError("Wavelength vector length does not match the spectral dimension.")
            wavelength_source = "FITS spectral WCS"
        if wavelengths.size < 2 or len(np.unique(np.round(wavelengths, 10))) < 2:
            raise ValueError("Cube needs at least two distinct wavelength channels.")
        def extension(name: str):
            item = next((candidate for candidate in hdul if str(candidate.name).upper() == name), None)
            if item is None or item.data is None:
                return None
            arr = np.asarray(item.data)
            return arr.copy() if arr.shape == data.shape else None
        return SpectralCube(
            data=data.copy(),
            wavelengths_um=np.asarray(wavelengths, dtype=float),
            celestial_wcs=wcs.celestial,
            spectral_numpy_axis=spectral_axis,
            header=dict(image_hdu.header),
            source_name=str(source_name),
            sha256=_hash_bytes(raw),
            bunit=str(image_hdu.header.get("BUNIT", "not supplied")),
            variance=extension("VARIANCE"),
            flags=extension("FLAGS"),
            nhits=extension("NHITS"),
            metadata={
                "image_extension": str(image_hdu.name),
                "shape_numpy": list(data.shape),
                "spectral_numpy_axis": spectral_axis,
                "wavelength_source": wavelength_source,
                "wavelength_min_um": float(np.nanmin(wavelengths)),
                "wavelength_max_um": float(np.nanmax(wavelengths)),
                "wavelength_channels": int(wavelengths.size),
                "celestial_wcs_valid": True,
                "sha256": _hash_bytes(raw),
                "bunit": str(image_hdu.header.get("BUNIT", "not supplied")),
            },
        )


def extract_aperture_spectrum(
    cube: SpectralCube,
    ra_deg: float,
    dec_deg: float,
    aperture_radius_px: float = 2.0,
    annulus_inner_px: float = 5.0,
    annulus_outer_px: float = 8.0,
    max_channels: int = 120,
) -> dict[str, Any]:
    """Extract per-channel aperture-minus-annulus spectrum at an ICRS position."""
    if not 0.1 <= aperture_radius_px < annulus_inner_px < annulus_outer_px <= 100:
        raise ValueError("Require 0.1 <= aperture radius < inner annulus < outer annulus <= 100 pixels.")
    if cube.wavelengths_um.size > max_channels:
        raise ValueError(f"Cube has {cube.wavelengths_um.size} channels, above the {max_channels}-channel quick-look limit.")
    coordinate = SkyCoord(float(ra_deg) * u.deg, float(dec_deg) * u.deg, frame="icrs")
    x, y = cube.celestial_wcs.world_to_pixel(coordinate)
    planes = cube.spectral_first
    if planes.ndim != 3:
        raise ValueError("Spectral cube layout could not be reduced to [wavelength, y, x].")
    if not (0 <= x < planes.shape[2] and 0 <= y < planes.shape[1]):
        raise ValueError("Selected sky position falls outside the cube coverage.")
    aperture = CircularAperture([(float(x), float(y))], r=float(aperture_radius_px))
    annulus = CircularAnnulus([(float(x), float(y))], r_in=float(annulus_inner_px), r_out=float(annulus_outer_px))
    ap_area = float(aperture.area)
    ann_area = float(annulus.area)
    net_values, bkg_values, valid_fractions = [], [], []
    errors = []
    variance_planes = np.moveaxis(cube.variance, cube.spectral_numpy_axis, 0).astype(float) if cube.variance is not None else None
    flag_planes = np.moveaxis(cube.flags, cube.spectral_numpy_axis, 0) if cube.flags is not None else None
    hit_planes = np.moveaxis(cube.nhits, cube.spectral_numpy_axis, 0) if cube.nhits is not None else None
    for channel, plane in enumerate(planes):
        image = np.asarray(plane, dtype=float)
        mask = ~np.isfinite(image)
        # FLAGS can include informational bits, so retain it for diagnostics.
        # NHITS=0 is a safe, explicit no-coverage mask.
        if hit_planes is not None:
            mask |= np.asarray(hit_planes[channel]) <= 0
        ap_sum = float(aperture_photometry(np.nan_to_num(image, nan=0.0), aperture, mask=mask, method="exact")["aperture_sum"][0])
        ann_sum = float(aperture_photometry(np.nan_to_num(image, nan=0.0), annulus, mask=mask, method="exact")["aperture_sum"][0])
        ap_valid = float(aperture.area_overlap(image, mask=mask, method="exact"))
        ann_valid = float(annulus.area_overlap(image, mask=mask, method="exact"))
        if ann_valid <= 0 or ap_valid <= 0:
            bkg = float("nan")
            net = float("nan")
            valid_fraction = 0.0
        else:
            bkg = ann_sum / ann_valid
            net = ap_sum - bkg * ap_valid
            valid_fraction = min(ap_valid / ap_area, ann_valid / ann_area)
        net_values.append(net)
        bkg_values.append(bkg)
        valid_fractions.append(valid_fraction)
        if variance_planes is not None and np.isfinite(net):
            var_img = np.asarray(variance_planes[channel], dtype=float)
            var_mask = mask | ~np.isfinite(var_img) | (var_img < 0)
            var_ap = float(aperture_photometry(np.nan_to_num(var_img, nan=0.0), aperture, mask=var_mask, method="exact")["aperture_sum"][0])
            errors.append(float(np.sqrt(max(var_ap, 0.0))))
        else:
            errors.append(float("nan"))
    flux = np.asarray(net_values, dtype=float)
    unit = cube.bunit
    flux_unit = f"{unit} pixel-sum"
    try:
        pixel_area_sr = float((proj_plane_pixel_area(cube.celestial_wcs) * u.deg**2).to_value(u.sr))
    except Exception:
        pixel_area_sr = None
    if pixel_area_sr and "mjy" in unit.lower() and "sr" in unit.lower():
        flux = flux * pixel_area_sr * 1e9
        errors_array = np.asarray(errors, dtype=float) * pixel_area_sr * 1e9
        flux_unit = "mJy (aperture flux density; MJy/sr × pixel solid angle)"
    else:
        errors_array = np.asarray(errors, dtype=float)
    finite_channels = np.isfinite(flux) & np.isfinite(cube.wavelengths_um)
    return {
        "status": "ok",
        "ra_deg": float(ra_deg),
        "dec_deg": float(dec_deg),
        "pixel_x": float(x),
        "pixel_y": float(y),
        "wavelength_um": cube.wavelengths_um.copy(),
        "flux": flux,
        "background_per_pixel": np.asarray(bkg_values, dtype=float),
        "flux_error": errors_array,
        "valid_fraction": np.asarray(valid_fractions, dtype=float),
        "flagged_pixel_fraction_at_target": float(np.mean(flag_planes[:, int(round(y)), int(round(x))] != 0)) if flag_planes is not None else None,
        "flux_unit": flux_unit,
        "aperture_radius_px": float(aperture_radius_px),
        "background_annulus_px": [float(annulus_inner_px), float(annulus_outer_px)],
        "finite_channels": int(np.count_nonzero(finite_channels)),
        "cube_sha256": cube.sha256,
        "extraction": "CircularAperture flux sum minus annulus median estimate per channel; flagged pixels excluded; any variance-based errors omit background-estimation and inter-pixel covariance terms.",
        "limitations": [
            "This mosaic/HiPS spectrum is a quick-look measurement, not a substitute for the IRSA SPHEREx spectrophotometry product.",
            "The background annulus is local and can be biased by extended emission, crowding or gradients.",
            "HiPS all-sky cube products are for visualization; they do not preserve time-separated spectral-image epochs.",
        ],
        "metadata": cube.metadata,
    }
