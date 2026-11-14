"""Data Explorer: raw archive metadata, FITS headers, FLAGS tables and provenance."""
from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import streamlit as st

import config
from components import metric_cards, navbar, sidebar, status_panel
from modules import provenance, state

sidebar.bootstrap_page(
    "Data Explorer",
    "Raw metadata, FITS headers, FLAGS bit definitions and provenance records behind the interface.",
    "Data Explorer",
)

observations = state.get("observations", [])
products = state.get("fits_products", {})
comparison = state.get("comparison")

metric_cards.render_metric_grid([
    {"label": "Observation records", "value": str(len(observations)) if observations else "Awaiting data", "caption": "Archive rows in this session", "icon": "☷"},
    {"label": "FITS products loaded", "value": str(len(products)) if products else "Awaiting data", "caption": "Validated arrays in memory", "icon": "▤"},
    {"label": "Candidate records", "value": str(len((comparison or {}).get("candidates", []))) if comparison else "Awaiting data", "caption": "From the current processing run", "icon": "✦"},
    {"label": "Data mode", "value": state.data_mode_label(), "caption": state.get("mode_note", "")[:70], "icon": "◈"},
])

st.markdown("### Archive metadata as returned")
if observations:
    tap_rows = [obs.get("metadata", {}).get("tap_row", {}) for obs in observations]
    if any(tap_rows):
        st.dataframe(pd.DataFrame(tap_rows), use_container_width=True, hide_index=True)
        st.caption("These are the unmodified columns returned by the documented IRSA TAP query (spherex.artifact joined to spherex.plane).")
    with st.expander("Normalised observation records used by the interface", expanded=False):
        st.dataframe(
            pd.DataFrame([
                {
                    "id": obs.get("id"),
                    "obs_id": obs.get("obs_id"),
                    "date_utc": obs.get("date_utc"),
                    "mjd": obs.get("mjd"),
                    "band": obs.get("band"),
                    "release": obs.get("release"),
                    "calibration_level": obs.get("calibration_level"),
                    "quality_flag": obs.get("quality_flag"),
                    "provenance_version": obs.get("provenance_version"),
                    "coordinates_source": obs.get("coordinates_source"),
                    "cutout_url": obs.get("cutout_url"),
                    "cached_path": obs.get("cached_path"),
                }
                for obs in observations
            ]),
            use_container_width=True, hide_index=True,
        )
    st.download_button(
        "DOWNLOAD OBSERVATION RECORDS (JSON)",
        data=json.dumps(provenance.safe_json(observations), indent=2),
        file_name="observations.json",
        mime="application/json",
    )
else:
    st.info("No observation metadata in this session. Use Sky Explorer to query IRSA or load the cached/demonstration bundle.")

if products:
    st.markdown("### FITS header inspection")
    labels = {identifier: f"{identifier} · {product.metadata.get('observation', {}).get('obs_id', '')}" for identifier, product in products.items()}
    chosen = st.selectbox("Product", list(products), format_func=lambda value: labels[value])
    product = products[chosen]
    tabs = st.tabs(["IMAGE header", "Primary header", "Validation", "FLAGS definitions", "Extensions"])
    with tabs[0]:
        st.dataframe(pd.DataFrame([{"keyword": key, "value": str(value)} for key, value in product.header.items()]), use_container_width=True, hide_index=True, height=400)
    with tabs[1]:
        st.dataframe(pd.DataFrame([{"keyword": key, "value": str(value)} for key, value in product.primary_header.items()]), use_container_width=True, hide_index=True, height=300)
    with tabs[2]:
        st.dataframe(pd.DataFrame(product.validation), use_container_width=True, hide_index=True)
        st.caption(f"File: {product.path} · {product.file_size_bytes:,} bytes · sha256 {product.sha256}")
    with tabs[3]:
        if product.flag_definitions:
            st.dataframe(
                pd.DataFrame([
                    {
                        "flag": name,
                        "bit mask": definition["bit_mask"],
                        "treated as": "severe (masked)" if name not in config.INFO_FLAG_NAMES else "informational (kept)",
                        "description": definition.get("description", ""),
                    }
                    for name, definition in product.flag_definitions.items()
                ]),
                use_container_width=True, hide_index=True,
            )
            st.caption(
                "Bit names and masks are parsed from each FITS header at load time. SOURCE and FULLSAMPLE are informational: masking them would "
                "remove nearly every pixel of a real cutout, so they are never used as defects."
            )
            st.json(product.metadata.get("flag_diagnostics", {}))
        else:
            st.info("This product has no FLAGS extension or no named bit definitions in its header.")
    with tabs[4]:
        st.write(", ".join(product.extension_names))
        st.json({key: value for key, value in product.metadata.items() if key != "observation"})

st.markdown("### Processing provenance")
if comparison and comparison.get("provenance"):
    record = comparison["provenance"]
    st.caption(f"Run {record['run_id']} created {record['created_utc']} · data mode {record['data_mode']}")
    st.json(record)
    st.download_button(
        "DOWNLOAD PROVENANCE JSON",
        data=json.dumps(record, indent=2, default=str),
        file_name=f"provenance_{record['run_id']}.json",
        mime="application/json",
    )
else:
    runs = sorted((config.OUTPUT_DIR / "runs").glob("*.json")) if (config.OUTPUT_DIR / "runs").exists() else []
    if runs:
        st.caption(f"{len(runs)} provenance record(s) written to outputs/runs in previous sessions.")
        latest = provenance.read_json(runs[-1])
        st.json(latest)
    else:
        st.info("No processing run has been executed yet, so there is no provenance record to inspect.")

st.markdown("### Session event log")
logs = state.get("logs", [])
if logs:
    st.dataframe(pd.DataFrame(list(reversed(logs))), use_container_width=True, hide_index=True, height=320)
    st.caption("Events are also appended to outputs/events.jsonl with an application log at outputs/application.log.")
else:
    st.info("No events logged yet.")

st.markdown("### Local data directories")
directory_rows = []
for path, description in [
    (config.REAL_DATA_DIR, "Verified real SPHEREx cutouts + manifest"),
    (config.DEMO_DATA_DIR, "Synthetic demonstration products (never observations)"),
    (config.UPLOAD_DIR, "User-supplied FITS uploads"),
    (config.FITS_CACHE_DIR, "HTTP cutout cache keyed by URL hash"),
    (config.OUTPUT_DIR / "runs", "Provenance JSON per processing run"),
    (config.OUTPUT_DIR / "candidates", "Candidate CSV exports"),
]:
    exists = Path(path).exists()
    files = sorted(Path(path).glob("*")) if exists else []
    total = sum(item.stat().st_size for item in files if item.is_file())
    directory_rows.append({
        "path": str(Path(path).relative_to(config.ROOT)) if Path(path).is_relative_to(config.ROOT) else str(path),
        "purpose": description,
        "files": len(files),
        "size": f"{total / 1024 / 1024:.2f} MB",
        "status": "present" if exists else "missing",
    })
st.dataframe(pd.DataFrame(directory_rows), use_container_width=True, hide_index=True)

navbar.render_footer()
