"""Controlled offline tests for multi-epoch association and 2-epoch limits."""
from __future__ import annotations

from modules import multi_epoch


def _detection(epoch_id: str, mjd: float, xi: float) -> dict:
    return {
        "observation_id": epoch_id,
        "observation_date": f"2026-01-0{int(mjd - 60000) + 1}T00:00:00",
        "mjd": mjd,
        "xi_arcsec": xi,
        "eta_arcsec": 0.0,
        "ra_deg": 10.0 + xi / 3600.0,
        "dec_deg": 20.0,
        "position_sigma_arcsec": 1.0,
        "peak_snr": 8.0,
        "severe_flag_fraction": 0.0,
        "support_pixels": 4,
        "artifact_flags": [],
    }


def test_three_compatible_epochs_associate_and_fit_linear_motion():
    detections = [
        [_detection("E1", 60000.0, 0.0)],
        [_detection("E2", 60001.0, 5.0)],
        [_detection("E3", 60002.0, 10.0)],
    ]
    pairs = {("E1", "E2"), ("E2", "E3")}
    tracks = multi_epoch._associate_epochs(detections, multi_epoch.DEFAULT_MULTI_EPOCH_CONFIG, pairs)

    assert len(tracks) == 1
    report = multi_epoch._track_report(1, tracks[0], 3, multi_epoch.DEFAULT_MULTI_EPOCH_CONFIG)
    assert report["detection_count"] == 3
    assert report["motion"]["speed_arcsec_per_day"] > 4.9
    assert report["classification"] == "potential_moving_source"
    assert "2-EPOCH-ONLY" not in report["quality_flags"]


def test_incompatible_pair_is_not_associated():
    detections = [
        [_detection("E1", 60000.0, 0.0)],
        [_detection("E2", 60001.0, 1.0)],
    ]
    tracks = multi_epoch._associate_epochs(detections, multi_epoch.DEFAULT_MULTI_EPOCH_CONFIG, set())
    assert len(tracks) == 2
    assert all(len(track["detections"]) == 1 for track in tracks)


def test_two_detection_track_keeps_explicit_two_epoch_only_limit():
    detections = [_detection("E1", 60000.0, 0.0), _detection("E2", 60001.0, 1.0)]
    report = multi_epoch._track_report(
        1,
        {"detections": detections, "quality_flags": [], "association_ambiguous": False},
        2,
        multi_epoch.DEFAULT_MULTI_EPOCH_CONFIG,
    )
    assert report["classification"] == "insufficient_evidence"
    assert "2-EPOCH-ONLY" in report["quality_flags"]
    assert "2-EPOCH-ONLY" in " ".join(report["evidence"])

