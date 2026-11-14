"""Mission Control — real application state dashboard."""
from __future__ import annotations

from html import escape

import pandas as pd
import streamlit as st

import config
from components import metric_cards, navbar, status_panel
from components.sidebar import bootstrap_page
from modules import state

st.set_page_config(page_title="Mission Control · THE CHAOS THEORY", page_icon="◉", layout="wide")
bootstrap_page(
    "MISSION CONTROL",
    "Application state · archive status · processing provenance",
    "Overview",
)

summary = state.state_summary()
query = summary["selected_region"]
region = (
    f"RA {float(query.get('ra_deg', 0)):.5f}° / Dec {float(query.get('dec_deg', 0)):+.5f}°"
    if summary["observation_count"] else "Awaiting data"
)
band = query.get("band", "Awaiting data") if summary["observation_count"] else "Awaiting data"
archive = summary["archive_status"]

metric_cards.render_metric_grid([
    {"label": "Selected sky region", "value": region, "caption": "Requested ICRS centre; archive coverage is verified separately.", "icon": "⌖", "accent": True},
    {"label": "Available observations", "value": summary["observation_count"] if summary["observation_count"] else "Awaiting data", "caption": "Rows actually returned by IRSA TAP or a labelled local manifest.", "icon": "≡"},
    {"label": "Selected spectral band", "value": band, "caption": "Detector/LVF band; central wavelength is checked per FITS WCS-WAVE.", "icon": "◈"},
    {"label": "Detected candidates", "value": summary["candidate_count"] if summary["comparison_ready"] else "Awaiting data", "caption": "Residual sources from the currently loaded scientific comparison, if any.", "icon": "✧"},
    {"label": "Processing state", "value": summary["processing_state"], "caption": f"{summary['product_count']} FITS product(s) loaded in this session.", "icon": "⚙"},
])

st.markdown("<div style='height:1.2rem'></div>", unsafe_allow_html=True)
left, right = st.columns([1.35, 1])
with left:
    st.markdown("### Scientific workflow state")
    stages = [
        ("Sky region requested", bool(summary["observation_count"])),
        ("Observation metadata retrieved", bool(summary["observation_count"])),
        ("FITS products loaded", bool(summary["product_count"])),
        ("Images aligned and compared", summary["comparison_ready"]),
        ("Candidate records generated", bool(summary["comparison_ready"])),
    ]
    html = "<div class='ct-pipeline'>"
    for i, (label, done) in enumerate(stages):
        if i:
            html += "<span class='ct-step-arrow'>→</span>"
        html += f"<span class='ct-step {'done' if done else ''}'>{escape(label)}</span>"
    html += "</div>"
    st.html(html)
    st.markdown("### Quick start")
    q1, q2, q3 = st.columns(3)
    with q1:
        st.page_link(navbar.page("pages/02_Sky_Explorer.py"), label="QUERY IRSA SKY REGION", use_container_width=True)
    with q2:
        st.page_link(navbar.page("pages/03_Observation_Explorer.py"), label="BROWSE OBSERVATIONS", use_container_width=True)
    with q3:
        st.page_link(navbar.page("pages/07_Methodology.py"), label="REVIEW METHODOLOGY", use_container_width=True)
    st.markdown("### One-click verified demo")
    demo_left, demo_right = st.columns([1, 1])
    with demo_left:
        if st.button("LOAD VERIFIED REAL PAIR (M101, 2026-05-17 / 2026-06-22)", use_container_width=True):
            try:
                from modules import data_loader

                records, manifest = data_loader.load_cached_real_bundle()
                state.set_observations(records, "CACHED_REAL", note=f"manifest {manifest['generated_utc']}")
                state.set_selected_pair([records[0]["id"], records[1]["id"]])
                products = {}
                for record in records:
                    products[record["id"]] = data_loader.load_observation_product(record)
                state.set_value("fits_products", products)
                state.set_products({record["id"]: {"path": products[record["id"]].path} for record in records})
                from modules import wcs_alignment

                state.set_value(
                    "pair_check",
                    wcs_alignment.check_observation_compatibility(
                        records[0], records[1], products[records[0]["id"]], products[records[1]["id"]],
                        state.get("processing_config", config.PROCESSING_DEFAULTS),
                    ),
                )
                st.success("Verified real SPHEREx pair loaded. Open the Cosmic Image Lab to run the difference pipeline.")
            except Exception as error:
                status_panel.render_exception(error, "Verified real pair load")
    with demo_right:
        if st.button("LOAD DEMONSTRATION DATA (SYNTHETIC)", use_container_width=True):
            try:
                from modules import data_loader, provenance, wcs_alignment

                manifest = provenance.read_json(config.DEMO_DATA_DIR / "manifest.json")
                records = manifest["observations"]
                state.set_observations(records, "DEMO", note=manifest["notice"])
                state.set_selected_pair([records[0]["id"], records[1]["id"]])
                products = {record["id"]: data_loader.load_observation_product(record) for record in records}
                state.set_value("fits_products", products)
                state.set_products({record["id"]: {"path": products[record["id"]].path} for record in records})
                state.set_value(
                    "pair_check",
                    wcs_alignment.check_observation_compatibility(
                        records[0], records[1], products[records[0]["id"]], products[records[1]["id"]],
                        state.get("processing_config", config.PROCESSING_DEFAULTS),
                    ),
                )
                st.warning("DEMONSTRATION DATA loaded — synthetic products, not SPHEREx observations.")
            except Exception as error:
                status_panel.render_exception(error, "Demonstration data load")
    st.markdown("### Session log")
    logs = state.get("logs", [])
    if logs:
        frame = pd.DataFrame(logs)
        st.dataframe(frame[["utc", "level", "message"]], use_container_width=True, hide_index=True, height=230)
    else:
        st.caption("No log entries in this session.")
with right:
    st.markdown("### Data integrity")
    st.html(
        "<div class='ct-glass'>"
        f"<div class='ct-tag {'cyan' if summary['mode_key'] != 'EMPTY' else ''}'>{escape(summary['data_mode'])}</div>"
        f"<p class='ct-quiet' style='margin-top:.8rem'><strong>Session</strong><br/>{escape(str(summary['session_id']))}<br/>"
        f"Created {escape(str(summary['created_utc']))}</p>"
        f"<p class='ct-quiet'><strong>IRSA status</strong><br/>{escape(str(archive.get('status','Not checked')))}. "
        f"{escape(str(archive.get('detail','')))}</p>"
        "<p class='ct-quiet'>No value in this dashboard is synthesized. Empty states are explicitly labelled. "
        "Demonstration data, local FITS uploads, real observations and processed measurements never share an unlabelled status.</p>"
        "</div>"
    )
    st.markdown("### Current processing configuration")
    cfg = state.get("processing_config", config.PROCESSING_DEFAULTS)
    st.json({k: cfg[k] for k in ["detection_sigma", "normalization_mode", "wavelength_tolerance_bandwidth_fraction", "crossmatch_enabled"]})

navbar.render_footer()
