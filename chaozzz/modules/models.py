"""Typed scientific data containers shared across processing layers."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np
from astropy.wcs import WCS


@dataclass
class Observation:
    id: str
    label: str
    origin: str
    data_mode: str
    ra_deg: float
    dec_deg: float
    date_utc: str | None = None
    mjd: float | None = None
    band: str = "unknown"
    detector: str = ""
    release: str = ""
    obs_id: str = ""
    cutout_url: str = ""
    full_fits_url: str = ""
    content_length: int | None = None
    calibration_level: int | None = None
    quality_flag: str = ""
    provenance_version: str = ""
    exposure_s: float | None = None
    energy_bounds_lower_m: float | None = None
    energy_bounds_upper_m: float | None = None
    product_type: str = ""
    content_type: str = ""
    plane_id: str = ""
    is_real_spherex: bool = False
    is_demo: bool = False
    central_wavelength_um: float | None = None
    central_bandwidth_um: float | None = None
    coordinates_source: str = "requested field centre"
    coordinates_role: str = "requested_field_center"
    metadata: dict[str, Any] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict[str, Any]:
        return {k: v for k, v in self.__dict__.items()}


@dataclass
class FitsProduct:
    product_id: str
    path: str
    origin: str
    data_mode: str
    image: np.ndarray
    flags: np.ndarray | None
    variance: np.ndarray | None
    zodi: np.ndarray | None
    wcs: WCS | None
    header: dict[str, Any]
    primary_header: dict[str, Any]
    extension_names: list[str]
    flag_definitions: dict[str, dict[str, Any]]
    bad_flag_mask: np.ndarray
    source_flag_mask: np.ndarray
    bunit: str
    shape: tuple[int, int]
    center_wavelength_um: float | None
    center_bandwidth_um: float | None
    center_pixel: tuple[float, float]  # Requested field centre pixel when available; geometric array centre otherwise.
    center_sky: tuple[float, float] | None  # Actual world coordinate at the geometric array centre.
    pixel_scale_arcsec: float | None
    file_size_bytes: int
    sha256: str
    metadata: dict[str, Any]
    warnings: list[str]
    validation: list[dict[str, str]]
    is_real_spherex: bool = False
    is_demo: bool = False

    @property
    def science(self) -> np.ndarray:
        """IMAGE minus documented ZODI model, in IMAGE units (usually MJy/sr)."""
        if self.zodi is None:
            return self.image.copy()
        return self.image - self.zodi


@dataclass
class AlignedPair:
    target_wcs: WCS
    shape: tuple[int, int]
    science1: np.ndarray
    science2: np.ndarray
    variance1: np.ndarray | None
    variance2: np.ndarray | None
    bad_mask1: np.ndarray
    bad_mask2: np.ndarray
    source_mask1: np.ndarray
    source_mask2: np.ndarray
    coverage1: np.ndarray
    coverage2: np.ndarray
    valid_mask: np.ndarray
    overlap_fraction: float
    pixel_scale_arcsec: float
    diagnostics: dict[str, Any]
    warnings: list[str]


@dataclass
class NormalizedPair:
    image1: np.ndarray
    image2: np.ndarray
    variance1: np.ndarray | None
    variance2: np.ndarray | None
    valid_mask: np.ndarray
    diagnostics: dict[str, Any]
    warnings: list[str]


@dataclass
class DifferenceResult:
    difference: np.ndarray
    snr: np.ndarray
    sigma_map: np.ndarray
    valid_mask: np.ndarray
    global_noise: float
    diagnostics: dict[str, Any]
    warnings: list[str]
    units: str = "MJy/sr"


@dataclass
class Candidate:
    candidate_id: str
    pixel_x: float
    pixel_y: float
    ra_deg: float
    dec_deg: float
    classification: str
    snr: float
    flux1_jy: float
    flux2_jy: float
    delta_flux_jy: float
    surface_brightness_diff: float
    local_noise: float
    quality: dict[str, Any]
    reasons: list[str]
    warnings: list[str]
    data_mode: str
    provenance_run_id: str = ""
    metadata: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        return {k: v for k, v in self.__dict__.items()}
