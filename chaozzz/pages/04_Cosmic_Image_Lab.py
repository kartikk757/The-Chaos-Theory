"""Cosmic Image Lab: FITS rendering, blink/overlay comparison and difference imaging."""
from __future__ import annotations

import numpy as np
import streamlit as st

import config
from components import image_viewer, metric_cards, navbar, sidebar, status_panel
from modules import pipeline, state, visualization

sidebar.bootstrap_page(
    "Cosmic Image Lab",
    "Render validated FITS arrays, compare epochs and compute a physically labelled difference image.",
    "Image Comparison",
)

products = state.get("fits_products", {})
pair = state.get("selected_pair", [])
observations = state.get("observations", [])
check = state.get("pair_check")

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
config_columns = st.columns(4)
processing = dict(state.get("processing_config", config.PROCESSING_DEFAULTS))
with config_columns[0]:
    processing["detection_sigma"] = st.slider("Detection threshold (σ)", 3.0, 12.0, float(processing.get("detection_sigma", 5.0)), 0.5)
with config_columns[1]:
    processing["fit_background_plane"] = st.toggle("Fit tilted background plane", value=bool(processing.get("fit_background_plane", True)))
with config_columns[2]:
    processing["aperture_radius_px"] = st.slider("Photometry aperture (px)", 1.0, 5.0, float(processing.get("aperture_radius_px", 2.0)), 0.5)
with config_columns[3]:
    processing["crossmatch_enabled"] = st.toggle("Catalogue cross-match", value=bool(processing.get("crossmatch_enabled", True)))

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
