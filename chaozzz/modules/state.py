"""Session-state helpers.

All user-specific state lives here and in ``st.session_state``. Heavy scientific
functions in ``modules/`` never read or write session state directly.
"""
from __future__ import annotations

import copy
import uuid
from datetime import datetime, timezone
from typing import Any

import config
from modules import provenance


def _defaults() -> dict[str, Any]:
    return {
        "session_id": "",
        "session_created_utc": "",
        "active_page": "Overview",
        "mode_key": "EMPTY",
        "mode_note": "",
        "query": copy.deepcopy(config.DEFAULT_QUERY),
        "observations": [],
        "coverage_rows": [],
        "coverage_status": {"status": "not queried", "detail": ""},
        "archive_status": {"status": "not checked", "detail": "Use Sky Explorer to query IRSA."},
        "selected_observation_ids": [],
        "selected_pair": [],
        "loaded_products": {},
        "fits_products": {},
        "pair_check": None,
        "processing_config": copy.deepcopy(config.PROCESSING_DEFAULTS),
        "comparison": None,
        "candidates": [],
        "selected_candidate_id": None,
        "provenance": None,
        "logs": [],
        "last_error": "",
        "viewer_controls": {},
    }


def init_state() -> None:
    import streamlit as st

    for key, value in _defaults().items():
        st.session_state.setdefault(key, copy.deepcopy(value))
    if not st.session_state.session_id:
        st.session_state.session_id = "CT-" + uuid.uuid4().hex[:10].upper()
        st.session_state.session_created_utc = datetime.now(timezone.utc).isoformat(timespec="seconds")
    if not st.session_state.logs:
        add_log("Session initialized", level="info", details={"app_version": config.APP_VERSION})


def get(key: str, default: Any = None) -> Any:
    import streamlit as st

    init_state()
    return st.session_state.get(key, default)


def set_value(key: str, value: Any) -> None:
    import streamlit as st

    init_state()
    st.session_state[key] = value


def data_mode_label() -> str:
    init_state()
    return config.DATA_MODES.get(get("mode_key", "EMPTY"), "UNKNOWN DATA MODE")


def set_data_mode(mode_key: str, note: str = "") -> None:
    if mode_key not in config.DATA_MODES:
        raise ValueError(f"Unknown data mode: {mode_key}")
    set_value("mode_key", mode_key)
    set_value("mode_note", note)
    add_log("Data mode changed", level="info", details={"mode": config.DATA_MODES[mode_key], "note": note})


def set_query(**updates: Any) -> None:
    query = copy.deepcopy(get("query", {}))
    query.update(updates)
    set_value("query", query)
    add_log("Sky query updated", level="info", details=updates)


def clear_results(keep_observations: bool = False) -> None:
    if not keep_observations:
        set_value("observations", [])
        set_value("coverage_rows", [])
        set_value("coverage_status", {"status": "not queried", "detail": ""})
    set_value("selected_observation_ids", [])
    set_value("selected_pair", [])
    set_value("loaded_products", {})
    set_value("fits_products", {})
    set_value("pair_check", None)
    set_value("comparison", None)
    set_value("candidates", [])
    set_value("selected_candidate_id", None)
    set_value("provenance", None)
    set_value("last_error", "")
    add_log("Processing results cleared", level="info")


def set_observations(records: list[dict], mode_key: str, note: str = "") -> None:
    set_data_mode(mode_key, note)
    set_value("observations", records)
    clear_results(keep_observations=True)
    add_log("Observation metadata loaded", level="info", details={"count": len(records), "mode": mode_key})


def set_coverage(rows: list[dict], status: dict) -> None:
    set_value("coverage_rows", rows)
    set_value("coverage_status", status)


def set_products(products: dict[str, dict]) -> None:
    current = copy.deepcopy(get("loaded_products", {}))
    current.update(products)
    set_value("loaded_products", current)
    add_log("FITS products loaded", level="info", details={"ids": list(products.keys())})


def set_selected_pair(ids: list[str]) -> None:
    if len(ids) != 2:
        raise ValueError("Exactly two observation IDs are required for a comparison pair")
    set_value("selected_pair", list(ids))
    set_value("selected_observation_ids", list(ids))
    add_log("Comparison pair selected", level="info", details={"ids": ids})


def set_processing_config(updates: dict) -> None:
    cfg = copy.deepcopy(get("processing_config", config.PROCESSING_DEFAULTS))
    cfg.update(updates)
    set_value("processing_config", cfg)
    add_log("Processing configuration updated", level="info", details=updates)


def _candidate_id(candidate: Any) -> str | None:
    if isinstance(candidate, dict):
        return candidate.get("candidate_id")
    return getattr(candidate, "candidate_id", None)


def set_comparison(result: dict) -> None:
    set_value("comparison", result)
    candidates = list(result.get("candidates", []) or [])
    set_value("candidates", candidates)
    set_value("provenance", result.get("provenance", None))
    if candidates:
        set_value("selected_candidate_id", _candidate_id(candidates[0]))
    else:
        set_value("selected_candidate_id", None)
    if result.get("status") not in (None, "ok"):
        set_value("last_error", str(result.get("error", "Processing failed")))
    add_log("Scientific comparison completed", level="info", details={
        "run_id": (result.get("provenance") or {}).get("run_id", ""),
        "status": result.get("status"),
        "candidates": len(candidates),
    })


def set_candidates(candidates: list[dict]) -> None:
    set_value("candidates", candidates)
    comparison = get("comparison")
    if comparison is not None:
        comparison = dict(comparison)
        comparison["candidates"] = candidates
        set_value("comparison", comparison)


def add_log(message: str, level: str = "info", details: dict | None = None) -> None:
    import streamlit as st

    # Do not call init_state() here: init_state() logs its own creation event,
    # which would otherwise recurse indefinitely.
    record = {
        "utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
        "level": level,
        "message": message,
        "details": details or {},
    }
    logs = list(st.session_state.get("logs", []))
    logs.append(record)
    st.session_state["logs"] = logs[-200:]
    provenance.log_event(message, level=level, **record["details"])


def reset_session() -> None:
    import streamlit as st

    preserve = {
        "session_id": st.session_state.get("session_id", ""),
        "session_created_utc": st.session_state.get("session_created_utc", ""),
    }
    for key in list(st.session_state.keys()):
        del st.session_state[key]
    for key, value in _defaults().items():
        st.session_state[key] = copy.deepcopy(value)
    st.session_state.update(preserve)
    add_log("Session reset", level="warning")


def observation_by_id(obs_id: str) -> dict | None:
    for record in get("observations", []):
        if record.get("id") == obs_id:
            return record
    return None


def loaded_product(obs_id: str) -> dict | None:
    return get("loaded_products", {}).get(obs_id)

def set_error(message: str, context: str = "") -> None:
    set_value("last_error", message)
    add_log("Error surfaced to user", level="error", details={"context": context, "message": message})


def candidate_list() -> list:
    return list(get("candidates", []) or [])


def selected_candidate():
    candidates = candidate_list()
    if not candidates:
        return None
    wanted = get("selected_candidate_id")
    for candidate in candidates:
        if _candidate_id(candidate) == wanted:
            return candidate
    return candidates[0]


def export_dir() -> str:
    return str(config.OUTPUT_DIR)



def state_summary() -> dict:
    observations = get("observations", [])
    products = get("loaded_products", {})
    candidates = get("candidates", [])
    comparison = get("comparison")
    query = get("query", {})
    return {
        "session_id": get("session_id"),
        "created_utc": get("session_created_utc"),
        "data_mode": data_mode_label(),
        "mode_key": get("mode_key"),
        "selected_region": query,
        "observation_count": len(observations),
        "product_count": len(products),
        "candidate_count": len(candidates),
        "comparison_ready": bool(comparison),
        "processing_state": (comparison or {}).get("status", "Idle") if comparison else "Idle",
        "run_id": (comparison or {}).get("run_id", ""),
        "last_error": get("last_error", ""),
        "archive_status": get("archive_status", {}),
    }
