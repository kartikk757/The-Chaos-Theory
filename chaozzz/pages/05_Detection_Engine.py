"""Detection engine: stage-by-stage transparency, candidate records and export."""
from __future__ import annotations

import numpy as np
import pandas as pd
import streamlit as st
import plotly.graph_objects as go
from dataclasses import replace

import config
from components import candidate_table, metric_cards, navbar, sidebar, status_panel
import json

from modules import multi_epoch, pipeline, spherex_sdt_adapter, state, validation, visualization

sidebar.bootstrap_page(
    "Detection Engine",
    "Every stage of the residual-detection pipeline, with its configuration and outputs exposed.",
    "Detection Engine",
)

comparison = state.get("comparison")
pair_run_available = bool(comparison and comparison.get("status") == "ok")
if not pair_run_available:
    if comparison:
        st.warning(comparison.get("error", "The last two-observation comparison did not complete."))
        st.dataframe(pd.DataFrame(comparison.get("stages", [])), use_container_width=True, hide_index=True)
        if comparison.get("compatibility"):
            st.dataframe(pd.DataFrame(comparison["compatibility"].get("checks", [])), use_container_width=True, hide_index=True)
    else:
        st.info("No two-observation residual run is available. The multi-epoch workflow below can run independently on validated epochs.")
    st.page_link(navbar.page("pages/03_Observation_Explorer.py"), label="LOAD OR CHOOSE OBSERVATION EPOCHS", use_container_width=True)
    comparison = {
        "status": "ok",
        "run_id": "no-two-epoch-run",
        "candidate_summary": {"counts": {}, "classification_definitions": {}},
        "candidates": [],
        "stages": [],
        "warnings": [],
    }

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
if pair_run_available:
    st.download_button(
        "DOWNLOAD PROVENANCE RECORD",
        data=json.dumps(provenance, indent=2, default=str),
        file_name=f"provenance_{comparison.get('run_id')}.json",
        mime="application/json",
        use_container_width=False,
    )

st.markdown("### Multi-epoch source association and motion screening")
epoch_ids = state.get("selected_epoch_ids", [])
observations = {row.get("id"): row for row in state.get("observations", [])}
products = state.get("fits_products", {})
ready_ids = [identifier for identifier in epoch_ids if identifier in products and identifier in observations]
if len(ready_ids) < 2:
    st.info("Load two or more validated epochs in Observation Explorer to run the pairwise coverage matrix and track screening.")
else:
    st.caption(f"Validated epoch set contains {len(ready_ids)} files. The processor sorts by UTC time, reports every compatibility pair, and requires three detections before a moving-source classification.")
    if st.button("RUN MULTI-EPOCH ASSOCIATION", type="primary", key="run_multi_epoch"):
        try:
            with st.spinner("Checking spectral compatibility and WCS footprints, then associating per-epoch detections…"):
                result = multi_epoch.analyse_multi_epoch(
                    [observations[identifier] for identifier in ready_ids],
                    {identifier: products[identifier] for identifier in ready_ids},
                    state.get("processing_config", config.PROCESSING_DEFAULTS),
                )
            state.set_value("multi_epoch_result", result)
            st.rerun()
        except Exception as error:
            st.error(f"Multi-epoch analysis failed: {error}")
    track_result = state.get("multi_epoch_result") or {}
    if track_result and track_result.get("status") == "ok":
        stats = track_result.get("summary", {})
        summary_cols = st.columns(4)
        summary_cols[0].metric("Epochs", str(stats.get("observations", 0)))
        summary_cols[1].metric("Peak detections", str(stats.get("detections", 0)))
        summary_cols[2].metric("Associated tracks", str(stats.get("tracks", 0)))
        summary_cols[3].metric("Potential moving leads", str(stats.get("potential_moving_source_tracks", 0)))
        if len(ready_ids) == 2:
            st.warning("2-EPOCH-ONLY: two positions provide a displacement; they cannot establish a fitted motion track or orbital trajectory.")
        if track_result.get("compatibility_matrix"):
            with st.expander("All epoch-pair compatibility and coverage checks", expanded=False):
                st.dataframe(pd.DataFrame(track_result["compatibility_matrix"]), use_container_width=True, hide_index=True)
        with st.expander("Per-epoch source detection diagnostics", expanded=False):
            st.dataframe(pd.DataFrame(track_result.get("epoch_diagnostics", [])), use_container_width=True, hide_index=True)
        tracks = track_result.get("tracks", [])
        track_rows = []
        for track in tracks:
            motion = track.get("motion", {})
            track_rows.append({
                "track": track.get("track_id"),
                "classification": track.get("status_label"),
                "detections": track.get("detection_count"),
                "speed (arcsec/day)": motion.get("speed_arcsec_per_day"),
                "speed uncertainty": motion.get("speed_sigma_arcsec_per_day"),
                "reduced chi-square": motion.get("reduced_chi2"),
                "quality flags": ", ".join(track.get("quality_flags", [])),
                "ambiguous association": track.get("association_ambiguous"),
            })
        st.dataframe(pd.DataFrame(track_rows), use_container_width=True, hide_index=True)
        if tracks:
            chosen_track = st.selectbox("Track evidence", tracks, format_func=lambda item: f"{item['track_id']} · {item['status_label']}", key="multi_epoch_track_detail")
            st.markdown("**Evidence**")
            for item in chosen_track.get("evidence", []):
                st.write(f"- {item}")
            st.markdown("**Limitations**")
            for item in chosen_track.get("limitations", []):
                st.write(f"- {item}")
            detections = chosen_track.get("detections", [])
            if detections:
                fig = go.Figure()
                fig.add_trace(go.Scatter(
                    x=[row["xi_arcsec"] for row in detections],
                    y=[row["eta_arcsec"] for row in detections],
                    mode="lines+markers+text",
                    text=[str(row.get("observation_date") or row.get("observation_id")) for row in detections],
                    textposition="top center",
                    customdata=[[row.get("mjd"), row.get("peak_snr"), row.get("observation_id")] for row in detections],
                    hovertemplate="ξ=%{x:.2f}″<br>η=%{y:.2f}″<br>MJD=%{customdata[0]}<br>SNR=%{customdata[1]:.1f}<br>%{customdata[2]}<extra></extra>",
                    name=chosen_track["track_id"],
                ))
                fig.update_layout(title=f"{chosen_track['track_id']} · local tangent-plane detections", xaxis_title="ξ east (arcsec)", yaxis_title="η north (arcsec)", template="plotly_dark", height=420)
                st.plotly_chart(fig, use_container_width=True, config={"displaylogo": False})
                st.dataframe(pd.DataFrame(detections), use_container_width=True, hide_index=True)
        with st.expander("Multi-epoch provenance and conservative interpretation", expanded=False):
            st.json({"run_id": track_result.get("run_id"), "limitations": track_result.get("limitations"), "provenance": track_result.get("provenance")})
            st.download_button(
                "DOWNLOAD MULTI-EPOCH REPORT",
                data=json.dumps(track_result, indent=2, default=str),
                file_name=f"multi_epoch_{track_result.get('run_id')}.json",
                mime="application/json",
                key="download_multi_epoch_report",
            )
    elif track_result and track_result.get("status") != "ok":
        st.warning(track_result.get("error", f"Multi-epoch result status: {track_result.get('status')}"))

st.markdown("### Optional IRSA Source Discovery Tool comparison")
sdt_state = spherex_sdt_adapter.availability()
st.caption(sdt_state["detail"] + " The built-in Photutils residual detector remains the primary path.")
if st.button("RUN OPTIONAL SExtractor COMPARISON", key="run_sdt", disabled=difference is None or aligned is None):
    try:
        comparison_product = replace(
            comparison["product1"],
            image=np.asarray(difference.difference, dtype=np.float32),
            zodi=None,
            wcs=aligned.target_wcs,
            shape=difference.difference.shape,
            variance=None,
            bad_flag_mask=np.zeros(difference.difference.shape, dtype=bool),
            header={"BUNIT": difference.units},
        )
        result = spherex_sdt_adapter.extract_sources(comparison_product)
        if result.get("status") == "ok":
            result["comparison"] = spherex_sdt_adapter.compare_with_photutils(
                result.get("sources", []),
                [{"xcentroid": item.pixel_x, "ycentroid": item.pixel_y} for item in candidates],
            )
        state.set_value("sdt_comparison", result)
        st.rerun()
    except Exception as error:
        st.error(f"Optional SExtractor comparison failed: {error}")
sdt_result = state.get("sdt_comparison") or {}
if sdt_result:
    if sdt_result.get("status") == "ok":
        st.success(f"SExtractor returned {sdt_result.get('count', 0)} catalogue source(s); {sdt_result.get('comparison', {}).get('matched_pairs', 0)} spatially coincide with built-in residual candidates within the configured pixel radius.")
        if sdt_result.get("sources"):
            st.dataframe(pd.DataFrame(sdt_result["sources"]), use_container_width=True, hide_index=True)
        st.json(sdt_result.get("comparison", {}))
        st.caption(sdt_result.get("limitation", ""))
    else:
        st.info(sdt_result.get("detail", sdt_result.get("status", "Optional extraction unavailable.")))

st.markdown("### Controlled injection/recovery and null-control evaluation")
st.caption(
    "This optional validation evaluates a local-maximum detector on a user-supplied control image and empirical pixel kernel. "
    "It is a controlled software experiment, not a measured SPHEREx sensitivity claim."
)
with st.expander("Configure an injection/recovery run", expanded=False):
    upload_cols = st.columns(2)
    control_file = upload_cols[0].file_uploader("Source-free control FITS image", type=["fits", "fit", "fts"], key="validation_control_fits")
    kernel_file = upload_cols[1].file_uploader("Empirical detector-sampled source kernel FITS", type=["fits", "fit", "fts"], key="validation_kernel_fits")
    st.warning("Use a control region you have inspected for sources. Use an observed detector-sampled source cutout or separately validated kernel; do not use the displayed oversampled ePSF directly.")
    source_free = st.checkbox("I have inspected the control and assume it is source-free", value=False, key="validation_source_free")
    settings_cols = st.columns(4)
    amplitude_levels = settings_cols[0].multiselect("Injected peak SNR levels", [3.0, 4.0, 5.0, 7.0, 10.0], default=[4.0, 5.0, 7.0], key="validation_amplitudes")
    injection_count = settings_cols[1].number_input("Injections per level", min_value=1, max_value=100, value=12, step=1, key="validation_count")
    detection_threshold = settings_cols[2].number_input("Detection threshold (σ)", min_value=1.0, max_value=20.0, value=5.0, step=0.5, key="validation_threshold")
    validation_seed = settings_cols[3].number_input("Random seed", min_value=0, max_value=2147483647, value=20261114, step=1, key="validation_seed")
    validate_clicked = st.button("RUN CONTROLLED VALIDATION", key="run_injection_recovery", disabled=control_file is None or kernel_file is None or not amplitude_levels)
    if validate_clicked:
        try:
            control_image, control_info = validation.load_2d_fits_image(control_file.getvalue())
            kernel_image, kernel_info = validation.load_2d_fits_image(kernel_file.getvalue())
            report = validation.run_injection_recovery(
                control_image,
                kernel_image,
                amplitudes_sigma=tuple(float(value) for value in amplitude_levels),
                injections_per_level=int(injection_count),
                threshold_sigma=float(detection_threshold),
                random_seed=int(validation_seed),
                control_assumed_source_free=source_free,
            )
            report["control_image"].update(control_info)
            report["empirical_kernel"].update(kernel_info)
            state.set_value("injection_recovery_report", report)
        except Exception as error:
            st.error(f"Controlled validation could not run: {error}")
    validation_report = state.get("injection_recovery_report") or {}
    if validation_report.get("status") == "ok":
        st.metric("No-injection control peak count", str(validation_report["no_injection_control"]["detected_peak_count"]), help="Only interpretable as a false-positive count if the source-free assumption is valid.")
        st.dataframe(pd.DataFrame(validation_report.get("levels", [])), use_container_width=True, hide_index=True)
        with st.expander("Validation configuration, checksums and limits", expanded=False):
            st.json({key: value for key, value in validation_report.items() if key != "no_injection_control"})
            st.json(validation_report.get("no_injection_control"))
            st.download_button("DOWNLOAD VALIDATION REPORT", data=json.dumps(validation_report, indent=2, default=str), file_name="injection_recovery_report.json", mime="application/json", key="download_validation_report")

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
