"""Detection engine: stage-by-stage transparency, candidate records and export."""
from __future__ import annotations

import numpy as np
import pandas as pd
import streamlit as st

import config
from components import candidate_table, metric_cards, navbar, sidebar, status_panel
import json

from modules import pipeline, state, visualization

sidebar.bootstrap_page(
    "Detection Engine",
    "Every stage of the residual-detection pipeline, with its configuration and outputs exposed.",
    "Detection Engine",
)

comparison = state.get("comparison")
if not comparison:
    status_panel.render_empty_state(
        "No processed comparison in this session",
        "Load two compatible observations and run the difference pipeline in the Cosmic Image Lab. Candidate records are never pre-filled.",
        "pages/04_Cosmic_Image_Lab.py",
        "OPEN COSMIC IMAGE LAB",
    )
    navbar.render_footer()
    st.stop()

if comparison.get("status") != "ok":
    st.error(comparison.get("error", "The comparison did not complete."))
    st.dataframe(pd.DataFrame(comparison.get("stages", [])), use_container_width=True, hide_index=True)
    if comparison.get("compatibility"):
        st.dataframe(pd.DataFrame(comparison["compatibility"].get("checks", [])), use_container_width=True, hide_index=True)
    st.page_link(navbar.page("pages/03_Observation_Explorer.py"), label="CHOOSE ANOTHER PAIR", use_container_width=True)
    navbar.render_footer()
    st.stop()

summary = comparison.get("candidate_summary", {})
counts = summary.get("counts", {})
candidates = comparison.get("candidates", [])
difference = comparison.get("difference")
aligned = comparison.get("aligned")

metric_cards.render_metric_grid([
    {"label": "Candidate records", "value": str(len(candidates)), "caption": f"run {comparison.get('run_id')}", "icon": "✦", "accent": True},
    {"label": "Moving-source candidates", "value": str(counts.get("Potential moving-source candidate", 0)), "caption": "Screening label, not a discovery", "icon": "⇄"},
    {"label": "Variable-source candidates", "value": str(counts.get("Potential variable-source candidate", 0)), "caption": "Brightness change at fixed position", "icon": "◐"},
    {"label": "Unmatched / artifact", "value": str(counts.get("Unmatched residual", 0) + counts.get("Likely image artifact (user review required)", 0)), "caption": "Require expert review before interpretation", "icon": "⚑"},
])

st.markdown("### Pipeline stages")
stage_rows = [
    {
        "stage": stage["stage"],
        "status": stage["status"],
        "detail": stage["detail"],
    }
    for stage in comparison.get("stages", [])
]
st.dataframe(pd.DataFrame(stage_rows), use_container_width=True, hide_index=True)

detection = comparison.get("detection", {})
with st.expander("Detection configuration and diagnostics", expanded=False):
    st.json(detection)
    st.caption(
        "Detection runs on the SNR map built from the official VARIANCE planes. Because SPHEREx full-frame pixels (~6.1″) undersample the PSF, "
        "a box-car peak finder is used instead of a multi-pixel DAO point-source model, and opposite-sign neighbouring peaks are merged as a "
        "dipole so that one shifted object is not counted twice."
    )

st.markdown("### Residual map with candidate positions")
if difference is not None and candidates:
    st.plotly_chart(
        visualization.candidate_overlay_figure(
            comparison["product1"].science,
            candidates,
            wcs=aligned.target_wcs if aligned else None,
            title="Epoch-1 image with residual candidates (circle = peak, label = candidate ID)",
        ),
        use_container_width=True, theme=None, config={"displaylogo": False, "scrollZoom": True},
    )
elif difference is not None:
    st.info("The pipeline found no residual above the configured threshold. A null result is a valid scientific outcome for this pair.")
    st.plotly_chart(
        visualization.difference_figure(difference.difference, "Difference image (no significant residuals)"),
        use_container_width=True, theme=None, config={"displaylogo": False},
    )

st.markdown("### Candidate records")
csv_payload = candidate_table.render_candidate_table(candidates, comparison.get("crossmatch"), "detect")

if candidates:
    legend = pd.DataFrame([
        {"classification": name, "definition": definition}
        for name, definition in summary.get("classification_definitions", {}).items()
    ])
    with st.expander("Classification definitions used here", expanded=False):
        st.dataframe(legend, use_container_width=True, hide_index=True)
        st.caption("Only these observational categories are used. The pipeline never assigns planet, asteroid or discovery labels.")

    columns = st.columns([1, 1])
    with columns[0]:
        if csv_payload:
            st.download_button(
                "EXPORT CANDIDATES CSV",
                data=csv_payload,
                file_name=f"the_chaos_theory_candidates_{comparison.get('run_id')}.csv",
                mime="text/csv",
                use_container_width=True,
            )
    with columns[1]:
        st.page_link(navbar.page("pages/06_Candidate_Investigation.py"), label="INVESTIGATE CANDIDATES →", use_container_width=True)

st.markdown("### Reproducibility")
provenance = comparison.get("provenance") or {}
st.caption(
    f"Run ID {comparison.get('run_id')} · deterministic hash of the input checksums, processing configuration and application version. "
    f"A JSON record was written to outputs/runs/{comparison.get('run_id')}.json."
)
st.download_button(
    "DOWNLOAD PROVENANCE RECORD",
    data=json.dumps(provenance, indent=2, default=str),
    file_name=f"provenance_{comparison.get('run_id')}.json",
    mime="application/json",
    use_container_width=False,
)

st.markdown("### What would strengthen or overturn these results")
st.markdown(
    "- **More epochs.** Two epochs cannot separate a slow trend from a one-off residual. Adding a third epoch tests whether a moving candidate "
    "continues along a consistent track and whether a variable candidate repeats.\n"
    "- **Independent imagery.** Cross-checking candidates in another survey or wavelength guards against SPHEREx-specific artifacts.\n"
    "- **Deeper catalogue matching.** A matched catalogue source usually explains a residual as a known object; an unmatched residual is not "
    "evidence of a new object, only of missing information.\n"
    "- **PSF-matched differencing.** Full PSF matching would reduce residuals caused by undersampling, which currently dominate at the few-percent level."
)

navbar.render_footer()
