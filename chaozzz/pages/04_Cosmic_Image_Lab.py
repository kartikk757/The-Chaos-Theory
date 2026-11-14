"""Cosmic Image Lab: FITS rendering, blink/overlay comparison and difference imaging."""
from __future__ import annotations

import numpy as np
import streamlit as st
import plotly.graph_objects as go
import hashlib

import config
from components import image_viewer, metric_cards, navbar, sidebar, status_panel
from modules import mosaic_cube, pipeline, state, visualization

sidebar.bootstrap_page(
    "Cosmic Image Lab",
    "Render validated FITS arrays, compare epochs and compute a physically labelled difference image.",
    "Image Comparison",
)

products = state.get("fits_products", {})
pair = state.get("selected_pair", [])
observations = state.get("observations", [])
check = state.get("pair_check")


def _render_temporal_sequence() -> None:
    """Display loaded epochs with navigation, blink, overlays and registered comparisons."""
    epoch_ids = state.get("selected_epoch_ids", [])
    record_map = {record.get("id"): record for record in observations}
    ordered = [identifier for identifier in epoch_ids if identifier in products and identifier in record_map]
    ordered.sort(key=lambda identifier: (
        record_map[identifier].get("mjd") is None,
        record_map[identifier].get("mjd") or 0.0,
        str(record_map[identifier].get("date_utc") or ""),
    ))
    if len(ordered) < 2:
        return
    st.markdown("### Multi-epoch temporal viewer")
    st.caption("Epochs are ordered by observation UTC/MJD. Blink and overlay are for visual comparison; registered differences require a compatible pair.")
    labels = {identifier: f"{record_map[identifier].get('date_utc') or 'undated'} · {record_map[identifier].get('obs_id') or identifier} · {record_map[identifier].get('band')}" for identifier in ordered}
    ref_id = st.selectbox("Reference epoch", ordered, format_func=lambda value: labels[value], key="temporal_reference_epoch")
    view_col, speed_col, overlay_col, autoplay_col = st.columns([1.5, 1, 1.2, 1])
    with view_col:
        temporal_mode = st.selectbox("Viewer mode", ["Single epoch", "Opacity overlay", "Difference image", "Blink sequence"], key="temporal_view_mode")
    with speed_col:
        blink_interval = st.slider("Blink speed (seconds)", 0.3, 3.0, 1.0, 0.1, key="temporal_blink_speed")
    with overlay_col:
        show_track_overlay = st.toggle("Overlay track detections", value=True, key="temporal_track_overlay")
    with autoplay_col:
        autoplay = st.toggle("Play", value=False, key="temporal_autoplay")
    idx_key = "temporal_current_index"
    current_index = int(st.session_state.get(idx_key, min(1, len(ordered) - 1)))
    if current_index >= len(ordered):
        current_index = len(ordered) - 1
    if current_index < 0:
        current_index = 0
    nav_cols = st.columns([1, 3, 1])
    with nav_cols[0]:
        if st.button("PREVIOUS EPOCH", key="temporal_previous", use_container_width=True):
            current_index = (current_index - 1) % len(ordered)
    with nav_cols[2]:
        if st.button("NEXT EPOCH", key="temporal_next", use_container_width=True):
            current_index = (current_index + 1) % len(ordered)
    current_index = nav_cols[1].slider("Current epoch", 0, len(ordered) - 1, current_index, key="temporal_current_slider")
    st.session_state[idx_key] = current_index
    current_id = ordered[current_index]
    if current_id == ref_id and len(ordered) > 1:
        current_id = ordered[(current_index + 1) % len(ordered)]
    record = record_map[current_id]
    product = products[current_id]
    st.caption(
        f"Current: {labels[current_id]} · MJD {record.get('mjd') if record.get('mjd') is not None else 'unknown'} · "
        f"FITS WCS geometric centre {product.center_sky if product.center_sky else 'unavailable'} ICRS deg · "
        f"coordinate role {product.metadata.get('record_position_role', 'unspecified')}"
    )
    multi_result = state.get("multi_epoch_result") or {}
    overlay_detections = [row for row in multi_result.get("detections", []) if str(row.get("observation_id")) == str(current_id)]
    if temporal_mode in ("Opacity overlay", "Difference image"):
        comparison_state = state.get("temporal_pair_comparison") or {}
        pair_key = comparison_state.get("pair_ids")
        result = comparison_state.get("result") if pair_key == [ref_id, current_id] else None
        if result is None and st.button("REGISTER REFERENCE AND CURRENT EPOCH", key="register_temporal_pair", type="primary"):
            try:
                ref_record = record_map[ref_id]
                cfg = dict(state.get("processing_config", config.PROCESSING_DEFAULTS))
                cfg["crossmatch_enabled"] = False
                cfg["gaia_crossmatch_enabled"] = False
                with st.spinner("Checking compatibility and registering the selected pair…"):
                    result = pipeline.run_pipeline(ref_record, record, cfg, catalogue_key="2MASS PSC")
                state.set_value("temporal_pair_comparison", {"pair_ids": [ref_id, current_id], "result": result})
            except Exception as error:
                st.error(f"Temporal pair registration failed: {error}")
        if result is None:
            st.info("Run registration for this exact epoch pair to enable a WCS-aligned opacity comparison or difference image.")
        elif result.get("status") != "ok":
            st.warning(result.get("error", "The selected epochs are not compatible for direct comparison."))
            st.dataframe(result.get("compatibility", {}).get("checks", []), use_container_width=True, hide_index=True)
        elif temporal_mode == "Difference image":
            registered_difference = result.get("difference")
            if registered_difference is not None:
                st.plotly_chart(
                    visualization.difference_figure(registered_difference.difference, f"Registered difference · {labels[current_id]} − {labels[ref_id]}"),
                    use_container_width=True, config={"displaylogo": False, "scrollZoom": True},
                )
                st.caption("This registered residual is a pairwise difference. It is not a multi-epoch motion fit.")
        else:
            aligned_pair = result.get("aligned")
            if aligned_pair is not None:
                opacity = st.slider("Current epoch opacity", 0.0, 1.0, 0.5, 0.05, key="temporal_overlay_opacity")
                first = visualization._stretch(aligned_pair.science1, 0.5, 99.5)
                second = visualization._stretch(aligned_pair.science2, 0.5, 99.5)
                blend = (1.0 - opacity) * first + opacity * second
                fig = visualization.image_figure(blend, f"Registered display blend · {labels[ref_id]} / {labels[current_id]}", "Inferno", 0, 100, wcs=aligned_pair.target_wcs, height=460)
                st.plotly_chart(fig, use_container_width=True, config={"displaylogo": False, "scrollZoom": True})
                st.caption("The blend is a display-only percentile composite after WCS registration, not a physical flux composite.")
    else:
        def _render_temporal_frame():
            frame_index = int(st.session_state.get(idx_key, current_index))
            frame_id = ordered[frame_index]
            frame_record, frame_product = record_map[frame_id], products[frame_id]
            fig = visualization.image_figure(
                frame_product.science,
                f"{frame_record.get('date_utc') or 'undated'} · {frame_record.get('obs_id') or frame_id}",
                "Inferno", 0.5, 99.5, wcs=frame_product.wcs, height=460,
            )
            frame_detections = [row for row in multi_result.get("detections", []) if str(row.get("observation_id")) == str(frame_id)]
            if show_track_overlay and frame_detections and frame_product.wcs is not None:
                fig.add_trace(go.Scatter(
                    x=[row["ra_deg"] for row in frame_detections], y=[row["dec_deg"] for row in frame_detections],
                    mode="markers", marker={"size": 9, "symbol": "circle-open", "color": "cyan", "line": {"width": 2}},
                    text=[f"{row.get('peak_snr', 0):.1f}σ" for row in frame_detections],
                    name="Multi-epoch detections", hovertemplate="RA=%{x:.6f}°<br>Dec=%{y:.6f}°<br>%{text}<extra></extra>",
                ))
            st.plotly_chart(fig, use_container_width=True, config={"displaylogo": False, "scrollZoom": True})
            st.caption(f"{labels[frame_id]} · MJD {frame_record.get('mjd') if frame_record.get('mjd') is not None else 'unknown'} · WCS centre {frame_product.center_sky if frame_product.center_sky else 'unavailable'} ICRS deg")
        if temporal_mode == "Blink sequence" and autoplay:
            @st.fragment(run_every=float(blink_interval))
            def _auto_blink_frame():
                frame_index = int(st.session_state.get(idx_key, current_index))
                st.session_state[idx_key] = (frame_index + 1) % len(ordered)
                _render_temporal_frame()
            _auto_blink_frame()
        else:
            _render_temporal_frame()
        if show_track_overlay and not overlay_detections:
            st.caption("No multi-epoch detections are associated with this observation yet; source markers remain empty.")
    st.markdown("---")


_render_temporal_sequence()

st.markdown("### Spectral cube quick-look")
st.caption(
    "Upload a SPHEREx mosaic or HiPS cube FITS file to extract a background-subtracted aperture spectrum. "
    "This is a separate quick-look spectral product; it does not provide time-separated motion evidence."
)
cube_upload = st.file_uploader("SPHEREx 3-D FITS cube", type=["fits", "fit", "fts"], key="spectral_cube_upload")
if cube_upload is not None:
    cube_bytes = cube_upload.getvalue()
    cube_digest = hashlib.sha256(cube_bytes).hexdigest()
    if st.session_state.get("spectral_cube_sha256") != cube_digest:
        try:
            st.session_state["spectral_cube"] = mosaic_cube.load_cube(cube_bytes, cube_upload.name)
            st.session_state["spectral_cube_sha256"] = cube_digest
            st.session_state.pop("spectral_spectrum", None)
            st.session_state.pop("spectral_cube_error", None)
        except Exception as error:
            st.session_state.pop("spectral_cube", None)
            st.session_state["spectral_cube_error"] = str(error)
if st.session_state.get("spectral_cube_error"):
    st.error(f"Cube validation failed: {st.session_state.pop('spectral_cube_error')}")
cube = st.session_state.get("spectral_cube")
if cube is not None:
    st.success(
        f"Validated {cube.source_name}: {cube.data.shape} pixels, {cube.wavelengths_um.size} wavelength channels "
        f"({cube.wavelengths_um.min():.4f}–{cube.wavelengths_um.max():.4f} µm), BUNIT {cube.bunit}, SHA-256 {cube.sha256[:16]}…"
    )
    query_context = state.get("query", config.DEFAULT_QUERY) or config.DEFAULT_QUERY
    default_ra = float(query_context.get("ra_deg", config.DEFAULT_QUERY["ra_deg"]))
    default_dec = float(query_context.get("dec_deg", config.DEFAULT_QUERY["dec_deg"]))
    available_candidates = (state.get("comparison") or {}).get("candidates", [])
    coord_modes = ["Manual ICRS position", "Current residual candidate"] if available_candidates else ["Manual ICRS position"]
    coord_mode = st.radio("Spectrum position", coord_modes, horizontal=True, key="cube_position_mode")
    if coord_mode == "Current residual candidate":
        chosen = st.selectbox("Candidate", available_candidates, format_func=lambda item: item.candidate_id, key="cube_candidate")
        spectrum_ra, spectrum_dec = float(chosen.ra_deg), float(chosen.dec_deg)
        st.caption(f"Using candidate {chosen.candidate_id} at {spectrum_ra:.7f}°, {spectrum_dec:.7f}°. Candidate image coordinates are distinct from this cube measurement.")
    else:
        coord_cols = st.columns(2)
        spectrum_ra = coord_cols[0].number_input("RA (ICRS deg)", min_value=0.0, max_value=359.999999, value=default_ra, format="%.7f", key="cube_ra")
        spectrum_dec = coord_cols[1].number_input("Dec (ICRS deg)", min_value=-90.0, max_value=90.0, value=default_dec, format="%.7f", key="cube_dec")
    aperture_col, annulus_col, run_col = st.columns(3)
    aperture_radius = aperture_col.number_input("Aperture radius (px)", 0.5, 10.0, 2.0, 0.5, key="cube_aperture")
    annulus_inner = annulus_col.number_input("Background inner radius (px)", 1.0, 20.0, 5.0, 0.5, key="cube_annulus_in")
    run_cube = run_col.button("EXTRACT QUICK-LOOK SPECTRUM", type="primary", use_container_width=True, key="cube_extract")
    annulus_outer = st.number_input("Background outer radius (px)", min_value=float(annulus_inner + 0.5), max_value=30.0, value=max(8.0, float(annulus_inner + 1.0)), step=0.5, key="cube_annulus_out")
    if run_cube:
        try:
            st.session_state["spectral_spectrum"] = mosaic_cube.extract_aperture_spectrum(
                cube, float(spectrum_ra), float(spectrum_dec), float(aperture_radius), float(annulus_inner), float(annulus_outer)
            )
        except Exception as error:
            st.error(f"Spectrum extraction failed: {error}")
    spectrum = st.session_state.get("spectral_spectrum")
    if spectrum:
        figure = go.Figure()
        figure.add_trace(go.Scatter(x=spectrum["wavelength_um"], y=spectrum["flux"], mode="lines+markers", name="Aperture − background"))
        figure.update_layout(title="SPHEREx mosaic quick-look spectrum", xaxis_title="Wavelength (µm)", yaxis_title=spectrum["flux_unit"], template="plotly_dark", height=420)
        st.plotly_chart(figure, use_container_width=True, config={"displaylogo": False})
        st.caption(spectrum["extraction"])
        st.caption("This spectrum is not a time-domain measurement. For survey science, compare it with the official IRSA spectrophotometry products.")
        with st.expander("Spectral extraction provenance and limitations", expanded=False):
            st.json({key: value for key, value in spectrum.items() if key not in ("wavelength_um", "flux", "background_per_pixel", "flux_error", "valid_fraction")})
        csv_rows = ["wavelength_um,flux,flux_error,valid_fraction"]
        csv_rows.extend(f"{wave:.9g},{flux:.9g},{error:.9g},{valid:.6g}" for wave, flux, error, valid in zip(spectrum["wavelength_um"], spectrum["flux"], spectrum["flux_error"], spectrum["valid_fraction"]))
        st.download_button("DOWNLOAD SPECTRUM CSV", data="\n".join(csv_rows), file_name="spherex_quicklook_spectrum.csv", mime="text/csv", key="cube_spectrum_csv")
else:
    st.info("No spectral cube uploaded in this session. This optional tool is independent of the two-epoch difference pipeline.")
st.markdown("---")

if not pair or not all(identifier in products for identifier in pair):
    status_panel.render_empty_state(
        "No validated FITS pair loaded",
        "Load and validate two compatible observations in the Observation Explorer. Image arithmetic is never performed on unvalidated data.",
        "pages/03_Observation_Explorer.py",
        "OPEN OBSERVATION EXPLORER",
    )
    navbar.render_footer()
    st.stop()

if check and not check.get("compatible"):
    st.error("The selected pair failed the compatibility checks: " + " | ".join(check.get("reasons", [])))
    st.page_link(navbar.page("pages/03_Observation_Explorer.py"), label="CHOOSE ANOTHER PAIR", use_container_width=True)
    navbar.render_footer()
    st.stop()

record1 = next((obs for obs in observations if obs["id"] == pair[0]), {})
record2 = next((obs for obs in observations if obs["id"] == pair[1]), {})
product1, product2 = products[pair[0]], products[pair[1]]
comparison = state.get("comparison")
aligned = (comparison or {}).get("aligned")
difference = (comparison or {}).get("difference")

st.markdown("### Epochs under examination")
metric_cards.render_metric_grid([
    {"label": f"Epoch 1 · {record1.get('obs_id')}", "value": (record1.get("date_utc") or "undated")[:19], "caption": f"{record1.get('band')} · λ {product1.center_wavelength_um:.5f} µm" if product1.center_wavelength_um else f"{record1.get('band')}", "icon": "①", "accent": True},
    {"label": f"Epoch 2 · {record2.get('obs_id')}", "value": (record2.get("date_utc") or "undated")[:19], "caption": f"{record2.get('band')} · λ {product2.center_wavelength_um:.5f} µm" if product2.center_wavelength_um else f"{record2.get('band')}", "icon": "②", "accent": True},
    {"label": "Time baseline", "value": f"{(check or {}).get('time_separation_days', float('nan')):.2f} d" if (check or {}).get("time_separation_days") is not None else "unknown", "caption": "Measured from archive observation times", "icon": "⧗"},
    {"label": "Spectral separation", "value": f"Δλ {(check or {}).get('wavelength_delta_um', 0):.5f} µm" if (check or {}).get("wavelength_delta_um") is not None else "unverified", "caption": "Measured from each WCS-WAVE lookup table", "icon": "λ"},
])

st.markdown("### 01 · Single-epoch rendering")
controls = image_viewer.render_view_controls("lab")
reset_columns = st.columns([1, 3])
with reset_columns[0]:
    if st.button("RESET VIEW CONTROLS", use_container_width=True, help="Restore the default colormap and contrast percentiles."):
        for key in ("lab_cmap", "lab_low", "lab_high", "lab_invert"):
            st.session_state.pop(key, None)
        st.rerun()
with reset_columns[1]:
    st.caption(
        "Zoom with the mouse wheel and pan by dragging. Each panel is independent: synchronized zoom/pan across the two epochs is "
        "not implemented in this build, so use the blink comparator for pixel-by-pixel registration checks."
    )
tabs = st.tabs(["Epoch 1", "Epoch 2"])
with tabs[0]:
    image_viewer.render_image(product1.science, f"Epoch 1 · {record1.get('date_utc')} · IMAGE − ZODI", controls, wcs=product1.wcs)
    image_viewer.render_pixel_inspector(product1, "lab_p1")
with tabs[1]:
    image_viewer.render_image(product2.science, f"Epoch 2 · {record2.get('date_utc')} · IMAGE − ZODI", controls, wcs=product2.wcs)
    image_viewer.render_pixel_inspector(product2, "lab_p2")

st.caption(
    "Displayed arrays are IMAGE − ZODI in the file's native surface-brightness unit. The zodiacal model is subtracted as documented for "
    "SPHEREx Level-2 products; contrast limits affect display only and never the measurements."
)

st.markdown("### 02 · Comparison modes")
mode = st.segmented_control(
    "Comparison mode",
    ["Side by side", "Blink", "Opacity overlay", "Difference image"],
    default="Side by side",
    key="lab_mode",
)

if mode == "Side by side":
    columns = st.columns(2)
    with columns[0]:
        image_viewer.render_image(product1.science, f"Epoch 1 · {record1.get('date_utc')}", controls, wcs=product1.wcs)
    with columns[1]:
        image_viewer.render_image(product2.science, f"Epoch 2 · {record2.get('date_utc')}", controls, wcs=product2.wcs)
    st.caption("Both panels use the same contrast percentiles and colormap for a fair visual comparison.")

elif mode == "Blink":
    image_viewer.render_blink(product1, product2, controls, "lab")

elif mode == "Opacity overlay":
    image_viewer.render_overlay(product1, product2, controls, "lab")

else:
    if difference is None:
        st.info("Run the difference pipeline below to generate the difference image.")
    else:
        st.plotly_chart(
            visualization.difference_figure(
                difference.difference, f"D(x, y) = I₂ − I₁ · {record1.get('date_utc')} → {record2.get('date_utc')}"
            ),
            use_container_width=True, theme=None, config={"displaylogo": False, "scrollZoom": True},
        )
        st.plotly_chart(
            visualization.snr_figure(
                difference.snr, f"SNR of the difference · threshold {state.get('processing_config', {}).get('detection_sigma', 5.0):.1f}σ"
            ),
            use_container_width=True, theme=None, config={"displaylogo": False, "scrollZoom": True},
        )
        st.caption(
            "Positive (red) residuals mean the second epoch is brighter at that pixel. A difference image always contains residuals from "
            "resampling and background mismatch; it does not by itself establish that an object moved or changed."
        )

st.markdown("### 03 · Difference pipeline")
config_columns = st.columns(5)
processing = dict(state.get("processing_config", config.PROCESSING_DEFAULTS))
with config_columns[0]:
    processing["detection_sigma"] = st.slider("Detection threshold (σ)", 3.0, 12.0, float(processing.get("detection_sigma", 5.0)), 0.5)
with config_columns[1]:
    processing["fit_background_plane"] = st.toggle("Fit tilted background plane", value=bool(processing.get("fit_background_plane", True)))
with config_columns[2]:
    processing["aperture_radius_px"] = st.slider("Photometry aperture (px)", 1.0, 5.0, float(processing.get("aperture_radius_px", 2.0)), 0.5)
with config_columns[3]:
    processing["crossmatch_enabled"] = st.toggle("Catalogue cross-match", value=bool(processing.get("crossmatch_enabled", True)))
with config_columns[4]:
    processing["gaia_crossmatch_enabled"] = st.toggle("Gaia DR3 cross-match", value=bool(processing.get("gaia_crossmatch_enabled", True)))
st.caption("Gaia DR3 positions are propagated to the observation epoch using catalogue proper motions and matched with an uncertainty-scaled radius. The resulting match is a catalogue association, not an object identification.")

if st.button("RUN ALIGNMENT · NORMALIZATION · DIFFERENCE · DETECTION", type="primary", use_container_width=True):
    state.set_processing_config(processing)
    try:
        with st.spinner("Reprojecting onto a common WCS, correcting backgrounds, differencing and detecting residuals…"):
            result = pipeline.run_pipeline(record1, record2, processing, catalogue_key="2MASS PSC")
        if result["status"] == "ok":
            state.set_comparison(result)
            st.success(f"Pipeline complete · run {result['run_id']} · {len(result['candidates'])} candidate record(s).")
            st.rerun()
        else:
            state.set_comparison(result)
            st.error(result.get("error", "Processing did not complete."))
            st.dataframe(result.get("stages", []), use_container_width=True, hide_index=True)
    except Exception as error:
        status_panel.render_exception(error, "Difference pipeline")

if comparison:
    st.markdown("### 04 · Processing record")
    st.dataframe(
        [{"stage": stage["stage"], "status": stage["status"], "detail": stage["detail"]} for stage in comparison.get("stages", [])],
        use_container_width=True, hide_index=True,
    )
    if alignment := comparison.get("aligned"):
        st.caption(
            f"Overlap {alignment.overlap_fraction:.1%} · valid pixels {float(np.mean(alignment.valid_mask)):.1%} · "
            f"units {comparison.get('units')} · run {comparison.get('run_id')}"
        )
    if difference is not None:
        metrics = st.columns(4)
        metrics[0].metric("Global difference noise", f"{difference.global_noise:.4g}", help=f"Units: {difference.units}")
        metrics[1].metric("Median residual", f"{difference.diagnostics.get('difference_median', float('nan')):.4g}")
        metrics[2].metric("Residual max", f"{difference.diagnostics.get('difference_max', float('nan')):.4g}")
        metrics[3].metric("Noise model", "VARIANCE propagation" if "VARIANCE" in difference.diagnostics.get("noise_model", "") else "robust estimate")
        with st.expander("Normalization and noise diagnostics", expanded=False):
            st.json(comparison.get("normalized").diagnostics if comparison.get("normalized") else {})
            st.json(difference.diagnostics)
        for warning in comparison.get("warnings", []):
            st.warning(warning)
    image_viewer.render_download_buttons(product1, "lab_downloads")
    st.page_link(navbar.page("pages/05_Detection_Engine.py"), label="OPEN DETECTION ENGINE →", use_container_width=True)
else:
    st.info(
        "No processing has been run for this pair yet. Nothing is pre-computed: the difference image and every candidate record you see "
        "are produced from the two FITS files currently loaded, with the configuration shown above recorded in provenance."
    )

navbar.render_footer()
