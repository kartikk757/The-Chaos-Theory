"""End-to-end integration test of the processing pipeline.

Runs entirely on the labelled synthetic fixture (DEMONSTRATION DATA) so it needs
no network access. Results validate the software, not the sky.
"""
from __future__ import annotations

import copy
import json

import pandas as pd
import pytest

import config
from modules import pipeline, provenance


@pytest.fixture(scope="module")
def demo_result(demo_pair):
    records, _ = demo_pair
    cfg = pipeline.processing_config({"crossmatch_enabled": False})
    result = pipeline.run_pipeline(records[0], records[1], cfg, "2MASS PSC", use_cache=False)
    return records, cfg, result


def test_pipeline_completes_all_stages(demo_result):
    _, _, result = demo_result
    assert result["status"] == "ok", result.get("error")
    stage_names = [stage["stage"] for stage in result["stages"]]
    assert stage_names[:6] == [
        "Load observations", "Validate metadata", "Compatibility check",
        "WCS alignment", "Background correction", "Image comparison",
    ]
    assert "Source detection" in stage_names
    assert "Candidate measurement" in stage_names
    assert all(stage["status"] in {"ok", "warn"} for stage in result["stages"])


def test_pipeline_labels_data_mode_and_writes_provenance(demo_result):
    records, cfg, result = demo_result
    assert result["data_mode"] == config.DATA_MODES["DEMO"]
    record = result["provenance"]
    assert record["data_mode"] == config.DATA_MODES["DEMO"]
    assert len(record["inputs"]) == 2
    assert all(item["sha256"] for item in record["inputs"])
    coordinates = record["inputs"][0]["coordinate_provenance"]
    assert coordinates["coordinate_frame"] == "ICRS"
    assert coordinates["record_position_role"] == "synthetic_image_reference_center"
    assert coordinates["delivered_cutout_geometric_center_icrs_deg"] is not None
    assert coordinates["delivered_cutout_center_definition"].startswith("FITS WCS")
    assert record["processing_config"]["detection_sigma"] == cfg["detection_sigma"]
    assert "not an official nasa product" in record["attribution"]["notice"].lower()
    path = config.OUTPUT_DIR / "runs" / f"{record['run_id']}.json"
    assert path.exists()
    on_disk = json.loads(path.read_text(encoding="utf-8"))
    assert on_disk["run_id"] == record["run_id"]


def test_run_id_is_deterministic(demo_pair):
    records, _ = demo_pair
    cfg = pipeline.processing_config({"crossmatch_enabled": False})
    first = pipeline.run_pipeline(records[0], records[1], cfg, "2MASS PSC", use_cache=False)
    second = pipeline.run_pipeline(records[0], records[1], cfg, "2MASS PSC", use_cache=False)
    assert first["run_id"] == second["run_id"]


def test_run_id_changes_with_configuration(demo_pair):
    records, _ = demo_pair
    base = pipeline.run_pipeline(records[0], records[1], pipeline.processing_config({"crossmatch_enabled": False}), "2MASS PSC", use_cache=False)
    changed = pipeline.run_pipeline(records[0], records[1], pipeline.processing_config({"crossmatch_enabled": False, "detection_sigma": 6.0}), "2MASS PSC", use_cache=False)
    assert base["run_id"] != changed["run_id"]


def test_candidates_are_measured_and_exportable(demo_result):
    _, _, result = demo_result
    candidates = result["candidates"]
    assert candidates, "the synthetic truth model must produce residuals"
    candidate = candidates[0]
    assert candidate.ra_deg == candidate.ra_deg  # not NaN
    assert candidate.dec_deg == candidate.dec_deg
    assert candidate.snr > 0
    assert candidate.reasons
    assert candidate.metadata["epoch1"]["id"] != candidate.metadata["epoch2"]["id"]
    rows = [pipeline._candidate_row(item, result["crossmatch"].get(item.candidate_id)) for item in candidates]
    csv_text = pipeline.candidates_to_csv(rows)
    assert "candidate_id" in csv_text.splitlines()[0]
    assert "planetcandidate" not in csv_text.lower().replace(" ", "")
    exported = pd.read_csv(pd.io.common.StringIO(csv_text))
    assert len(exported) == len(candidates)
    written = sorted((config.OUTPUT_DIR / "candidates").glob(f"candidates_{result['run_id']}.csv"))
    assert written, "the pipeline must persist the candidate CSV for traceability"


def test_pipeline_refuses_an_incompatible_pair(demo_pair):
    records, _ = demo_pair
    incompatible = copy.deepcopy(records[1])
    incompatible["band"] = "SPHEREx-D4"
    cfg = pipeline.processing_config({"crossmatch_enabled": False})
    result = pipeline.run_pipeline(records[0], incompatible, cfg, "2MASS PSC", use_cache=False)
    assert result["status"] == "incompatible"
    assert any("band" in reason.lower() for reason in result["compatibility"]["reasons"])
    assert result["stages"][-1]["stage"] == "Compatibility check"


def test_pipeline_handles_unreadable_input_honestly(tmp_path):
    record = {
        "id": "BROKEN", "label": "missing file", "origin": "test",
        "data_mode": config.DATA_MODES["LOCAL"], "cached_path": str(tmp_path / "missing.fits"),
        "ra_deg": 10.0, "dec_deg": 10.0, "band": "unknown", "mjd": 61000.0,
    }
    result = pipeline.run_pipeline(record, record, pipeline.processing_config({"crossmatch_enabled": False}), "2MASS PSC", use_cache=False)
    assert result["status"] == "failed"
    assert "loading failed" in result["error"].lower()
    assert "does not exist" in result["error"] or "not exist" in result["error"]
