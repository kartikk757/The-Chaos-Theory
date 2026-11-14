"""Provenance, deterministic run identifiers, error logging and output helpers."""
from __future__ import annotations

import hashlib
import importlib.metadata
import json
import logging
import sys
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np

import config

_LOGGER: logging.Logger | None = None


def utc_now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def get_logger() -> logging.Logger:
    global _LOGGER
    if _LOGGER is not None:
        return _LOGGER
    logger = logging.getLogger("the_chaos_theory")
    logger.setLevel(logging.INFO)
    logger.propagate = False
    config.OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    path = config.OUTPUT_DIR / "application.log"
    handler_exists = any(
        isinstance(h, logging.FileHandler) and Path(h.baseFilename) == path.resolve()
        for h in logger.handlers
    )
    if not handler_exists:
        handler = logging.FileHandler(path, encoding="utf-8")
        handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
        logger.addHandler(handler)
    _LOGGER = logger
    return logger


def safe_json(value: Any) -> Any:
    if isinstance(value, dict):
        return {str(k): safe_json(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [safe_json(v) for v in value]
    if isinstance(value, np.ndarray):
        return value.tolist()
    if isinstance(value, np.generic):
        return value.item()
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, float) and (math.isnan(value) or math.isinf(value)):
        return None
    if hasattr(value, "to_string"):
        return value.to_string()
    return value


def write_json(path: str | Path, payload: Any) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(safe_json(payload), indent=2, sort_keys=True), encoding="utf-8")
    return path


def read_json(path: str | Path) -> Any:
    return json.loads(Path(path).read_text(encoding="utf-8"))


def log_event(message: str, level: str = "info", **details: Any) -> None:
    entry = {"utc": utc_now(), "level": level, "message": message, "details": safe_json(details)}
    config.OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    with (config.OUTPUT_DIR / "events.jsonl").open("a", encoding="utf-8") as handle:
        handle.write(json.dumps(entry, sort_keys=True) + "\n")
    logger = get_logger()
    logger.log(getattr(logging, level.upper(), logging.INFO), "%s %s", message, json.dumps(entry["details"]))


def sha256_file(path: str | Path) -> str:
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def sha256_bytes(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def deterministic_run_id(inputs: list[dict], processing_config: dict) -> str:
    material = {
        "inputs": [
            {
                "id": item.get("id"),
                "sha256": item.get("sha256"),
                "access_url": item.get("access_url"),
                "origin": item.get("origin"),
            }
            for item in inputs
        ],
        "config": safe_json(processing_config),
        "app_version": config.APP_VERSION,
    }
    digest = hashlib.sha256(json.dumps(material, sort_keys=True).encode("utf-8")).hexdigest()
    return "CT-" + digest[:12].upper()


def software_versions() -> dict[str, str]:
    packages = [
        "streamlit", "numpy", "scipy", "astropy", "astroquery", "pyvo",
        "photutils", "reproject", "pandas", "plotly", "matplotlib", "requests",
    ]
    result = {"python": sys.version.split()[0], "app": config.APP_VERSION}
    for package in packages:
        try:
            result[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            result[package] = "not installed"
    return result


def build_provenance(
    inputs: list[dict],
    processing_config: dict,
    stages: list[dict],
    warnings: list[str],
    data_mode: str,
    result_summary: dict,
    run_id: str | None = None,
) -> dict:
    """Assemble a provenance record.

    ``run_id`` should be supplied when the caller has already derived it from the
    exact inputs used for processing, so that filenames, CSV exports and the JSON
    record all carry the same identifier.
    """
    run_id = run_id or deterministic_run_id(inputs, processing_config)
    return {
        "run_id": run_id,
        "created_utc": utc_now(),
        "application": config.APP_TITLE,
        "data_mode": data_mode,
        "result_summary": safe_json(result_summary),
        "inputs": safe_json(inputs),
        "processing_config": safe_json(processing_config),
        "stages": safe_json(stages),
        "warnings": safe_json(warnings),
        "software_versions": software_versions(),
        "attribution": {
            "mission": "NASA SPHEREx",
            "archive": "NASA/IPAC Infrared Science Archive (IRSA), Caltech/IPAC",
            "documentation": [doc["url"] for doc in config.ARCHIVE_DOCS],
            "notice": "This prototype is not an official NASA product and makes no discovery claim.",
        },
    }
