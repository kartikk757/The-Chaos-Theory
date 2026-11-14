"""Controlled tests for empirical-kernel injection/recovery and null control."""
import numpy as np

from modules import validation


def test_injection_recovery_records_completeness_checksums_and_null_control():
    rng = np.random.default_rng(17)
    control = rng.normal(0.0, 1.0, (128, 128))
    kernel = np.array([[0.05, 0.1, 0.05], [0.1, 0.4, 0.1], [0.05, 0.1, 0.05]])
    report = validation.run_injection_recovery(
        control,
        kernel,
        amplitudes_sigma=(8.0,),
        injections_per_level=8,
        threshold_sigma=5.0,
        random_seed=2026,
        control_assumed_source_free=True,
    )
    level = report["levels"][0]
    assert report["status"] == "ok"
    assert level["recovered_sources"] == 8
    assert level["empirical_recovery_fraction"] == 1.0
    assert report["no_injection_control"]["detected_peak_count"] < 5
    assert report["control_image"]["assumed_source_free_by_user"] is True
    assert len(report["control_image"]["sha256"]) == 64
    assert any("not a calibrated" in text for text in report["scientific_limitations"])


def test_even_sized_kernel_is_refused_for_unambiguous_centering():
    with np.testing.assert_raises(ValueError):
        validation.run_injection_recovery(np.random.default_rng(1).normal(size=(64, 64)), np.ones((4, 4)))
