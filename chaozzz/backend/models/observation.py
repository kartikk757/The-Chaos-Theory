"""Structured observation metadata used by the service layer and the frontend contract."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import config


@dataclass
class ObservationRecord:
    """Normalised archive observation.

    Every field is taken from a real archive response, a verified local manifest or
    an explicitly labelled synthetic product. Nothing is synthesised here.
    """

    observation_id: str
    label: str
    collection: str
    data_mode: str
    ra: float | None
    dec: float | None
    mjd: float | None
    observation_date: str | None
    spectral_band: str
    detector: str
    release: str
    product_type: str
    calibration_level: int | None
    quality_flag: str
    provenance_version: str
    cutout_url: str = ""
    full_fits_url: str = ""
    cached_path: str = ""
    content_length: int | None = None
    exposure_seconds: float | None = None
    energy_bounds_um: tuple[float, float] | None = None
    measured_wavelength_um: float | None = None
    measured_bandwidth_um: float | None = None
    image_width: int | None = None
    image_height: int | None = None
    pixel_scale_arcsec: float | None = None
    is_real_spherex: bool = False
    is_demo: bool = False
    coordinates_source: str = ""
    coordinates_role: str = "recorded_observation_position"
    metadata: dict[str, Any] = field(default_factory=dict)
    warnings: list[str] = field(default_factory=list)

    # ------------------------------------------------------------------ factories
    @classmethod
    def from_archive_record(cls, record: dict) -> "ObservationRecord":
        bounds = record.get("wavelength_range_um")
        return cls(
            observation_id=str(record.get("obs_id") or record.get("id") or ""),
            label=str(record.get("label") or ""),
            collection="SPHEREx (NASA/IPAC IRSA)",
            data_mode=str(record.get("data_mode") or config.DATA_MODES["LIVE"]),
            ra=record.get("ra_deg"),
            dec=record.get("dec_deg"),
            mjd=record.get("mjd"),
            observation_date=record.get("date_utc"),
            spectral_band=str(record.get("band") or "unknown"),
            detector=str(record.get("detector") or ""),
            release=str(record.get("release") or ""),
            product_type=str(record.get("product_type") or ""),
            calibration_level=record.get("calibration_level"),
            quality_flag=str(record.get("quality_flag") or ""),
            provenance_version=str(record.get("provenance_version") or ""),
            cutout_url=str(record.get("cutout_url") or ""),
            full_fits_url=str(record.get("full_fits_url") or ""),
            cached_path=str(record.get("cached_path") or ""),
            content_length=record.get("content_length"),
            exposure_seconds=record.get("exposure_s"),
            energy_bounds_um=None if not bounds else (float(bounds[0]), float(bounds[1])),
            measured_wavelength_um=record.get("central_wavelength_um"),
            measured_bandwidth_um=record.get("central_bandwidth_um"),
            is_real_spherex=bool(record.get("is_real_spherex", False)),
            is_demo=bool(record.get("is_demo", False)),
            coordinates_source=str(record.get("coordinates_source") or ""),
            coordinates_role=str(record.get("coordinates_role") or "recorded_observation_position"),
            metadata=dict(record.get("metadata") or {}),
            warnings=list(record.get("warnings") or []),
        )

    # ------------------------------------------------------------------ contract
    def to_contract(self) -> dict[str, Any]:
        """JSON-compatible structure consumed by the frontend (Phase 9 schema)."""
        return {
            "observation_id": self.observation_id,
            "collection": self.collection,
            "observation_date": self.observation_date,
            "ra": self.ra,
            "dec": self.dec,
            "image_width": self.image_width,
            "image_height": self.image_height,
            "spectral_band": self.spectral_band,
            "data_status": "validated" if self.image_width else "metadata-only",
            "release": self.release,
            "calibration_level": self.calibration_level,
            "provenance_version": self.provenance_version,
            "quality_flag": self.quality_flag or None,
            "pixel_scale_arcsec": self.pixel_scale_arcsec,
            "measured_wavelength_um": self.measured_wavelength_um,
            "measured_bandwidth_um": self.measured_bandwidth_um,
            "exposure_seconds": self.exposure_seconds,
            "mjd": self.mjd,
            "data_mode": self.data_mode,
            "is_real_spherex": self.is_real_spherex,
            "is_demo": self.is_demo,
            "access_url": self.cutout_url or self.full_fits_url or None,
            "coordinates_source": self.coordinates_source,
            "coordinates_role": self.coordinates_role,
            "warnings": self.warnings,
        }

    def to_dict(self) -> dict[str, Any]:
        payload = dict(self.__dict__)
        payload["energy_bounds_um"] = list(self.energy_bounds_um) if self.energy_bounds_um else None
        return payload
