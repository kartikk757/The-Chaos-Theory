"""Candidate report contract for moving-source investigation."""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

# Classification vocabulary required by the investigation brief. These are
# observational categories, never physical identifications.
STATIONARY = "stationary_source"
POTENTIAL_MOVING = "potential_moving_source"
POSSIBLE_TRANSIENT = "possible_transient"
INSTRUMENTAL_ARTIFACT = "instrumental_artifact"
INSUFFICIENT = "insufficient_evidence"

CLASS_LABELS = {
    STATIONARY: "Stationary source (motion consistent with zero)",
    POTENTIAL_MOVING: "Potential moving source (investigative lead, not a discovery)",
    POSSIBLE_TRANSIENT: "Possible transient (seen in one epoch only)",
    INSTRUMENTAL_ARTIFACT: "Instrumental artifact (flagged or morphologically suspect)",
    INSUFFICIENT: "Insufficient evidence",
}


@dataclass
class EpochDetection:
    """One detection of one candidate in one epoch."""

    observation_id: str
    observation_date: str | None
    mjd: float
    x_pixel: float
    y_pixel: float
    ra: float
    dec: float
    xi_arcsec: float
    eta_arcsec: float
    peak_snr: float
    flux_jy: float | None
    flux_error_jy: float | None
    position_sigma_arcsec: float
    severe_flag_fraction: float
    pixels_above_threshold: int
    wavelength_um: float | None

    def to_contract(self) -> dict[str, Any]:
        return {
            "observation_id": self.observation_id,
            "observation_date": self.observation_date,
            "x_pixel": round(self.x_pixel, 3),
            "y_pixel": round(self.y_pixel, 3),
            "ra": round(self.ra, 7),
            "dec": round(self.dec, 7),
            "peak_snr": round(self.peak_snr, 2),
            "flux_jy": self.flux_jy,
            "flux_error_jy": self.flux_error_jy,
            "position_sigma_arcsec": round(self.position_sigma_arcsec, 3),
            "severe_flag_fraction": round(self.severe_flag_fraction, 4),
            "wavelength_um": self.wavelength_um,
        }


@dataclass
class CandidateReport:
    """Full scientific report for one multi-epoch candidate (Phase 7 / Phase 9)."""

    candidate_id: str
    status: str
    detections: list[EpochDetection]
    initial_position: dict[str, float]
    final_position: dict[str, float]
    angular_displacement_arcsec: float
    time_difference_hours: float
    angular_motion_arcsec_per_hour: float | None
    motion_uncertainty_arcsec_per_hour: float | None
    signal_to_noise: float
    measurement_uncertainty: float
    quality_flags: list[str]
    scientific_explanation: str
    evidence: list[str] = field(default_factory=list)
    limitations: list[str] = field(default_factory=list)
    motion_components: dict[str, float | None] = field(default_factory=dict)
    fit_statistics: dict[str, float | None] = field(default_factory=dict)
    brightness: dict[str, Any] = field(default_factory=dict)
    crossmatch: dict[str, Any] = field(default_factory=dict)
    ranking_score: dict[str, Any] = field(default_factory=dict)
    data_mode: str = ""
    provenance: dict[str, Any] = field(default_factory=dict)

    def to_contract(self) -> dict[str, Any]:
        """JSON-compatible payload following the Phase 9 candidate schema."""
        return {
            "candidate_id": self.candidate_id,
            "status": self.status,
            "status_label": CLASS_LABELS.get(self.status, self.status),
            "initial_position": self.initial_position,
            "final_position": self.final_position,
            "angular_displacement_arcsec": round(self.angular_displacement_arcsec, 4),
            "time_difference_hours": round(self.time_difference_hours, 4),
            "angular_motion_arcsec_per_hour": (
                None if self.angular_motion_arcsec_per_hour is None else round(self.angular_motion_arcsec_per_hour, 6)
            ),
            "motion_uncertainty_arcsec_per_hour": (
                None if self.motion_uncertainty_arcsec_per_hour is None else round(self.motion_uncertainty_arcsec_per_hour, 6)
            ),
            "signal_to_noise": round(self.signal_to_noise, 2),
            "measurement_uncertainty": round(self.measurement_uncertainty, 4),
            "quality_flags": self.quality_flags,
            "scientific_explanation": self.scientific_explanation,
            "motion_components": self.motion_components,
            "fit_statistics": self.fit_statistics,
            "brightness": self.brightness,
            "crossmatch": self.crossmatch,
            "ranking_score": self.ranking_score,
            "detections": [detection.to_contract() for detection in self.detections],
            "evidence": self.evidence,
            "limitations": self.limitations,
            "data_mode": self.data_mode,
            "provenance": self.provenance,
            "disclaimer": (
                "Investigative lead derived from archive measurements. This is not a discovery, not a confirmed "
                "solar-system object and not a Planet X identification."
            ),
        }

    def to_dict(self) -> dict[str, Any]:
        return self.to_contract()
