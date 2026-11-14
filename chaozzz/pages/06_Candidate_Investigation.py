"""Candidate investigation workspace: evidence, warnings, cross-match and export."""
from __future__ import annotations

import numpy as np
import streamlit as st

import config
from components import candidate_table, image_viewer, navbar, sidebar, status_panel
from modules import pipeline, state, visualization

sidebar.bootstrap_page(
    "Candidate Investigation",
    "Examine an individual residual candidate: position, brightness change, displacement and data quality.",
    "Candidate Analysis",
)

def _cutout(array, x, y, half):
    """Small pixel window around a candidate, clipped to the array bounds."""
    xi, yi = int(round(x)), int(round(y))
    return np.asarray(array)[max(yi - half, 0):yi + half + 1, max(xi - half, 0):xi + half + 1]


comparison = state.get("comparison")
candidates = (comparison or {}).get("candidates") or []
if not comparison or not candidates:
    status_panel.render_empty_state(
        "No candidates to investigate",
        "Run the detection pipeline on a compatible pair, or accept that this pair produced no residual above threshold. Nothing is invented here.",
        "pages/05_Detection_Engine.py",
        "OPEN DETECTION ENGINE",
    )
    navbar.render_footer()
    st.stop()

aligned = comparison.get("aligned")
difference = comparison.get("difference")
crossmatch = comparison.get("crossmatch", {})
product1 = comparison.get("product1")
product2 = comparison.get("product2")

csv_payload = candidate_table.render_candidate_table(candidates, crossmatch, "invest")
selected_id = st.session_state.get("selected_candidate_id")
candidate = next((item for item in candidates if item.candidate_id == selected_id), candidates[0])

st.markdown("---")
candidate_table.render_candidate_detail(candidate, crossmatch.get(candidate.candidate_id), "invest_detail")

st.markdown("### Where the candidate sits")
columns = st.columns(2)
with columns[0]:
    if product1 is not None:
        view = visualization.image_figure(
            product1.science,
            f"Epoch 1 · {candidate.metadata['epoch1'].get('date_utc')} · candidate marked",
            "Inferno",
            0.5,
            99.5,
            wcs=product1.wcs,
            height=420,
            marker={"x": candidate.pixel_x, "y": candidate.pixel_y, "label": candidate.candidate_id},
        )
        st.plotly_chart(view, use_container_width=True, theme=None, config={"displaylogo": False, "scrollZoom": True})
with columns[1]:
    if difference is not None:
        st.plotly_chart(
            visualization.cutout_figure(
                _cutout(difference.difference, candidate.pixel_x, candidate.pixel_y, 9),
                f"Difference cutout · ±4 pixels around {candidate.candidate_id}",
            ),
            use_container_width=True, theme=None, config={"displaylogo": False},
        )
        st.caption("The difference cutout shows the residual exactly as measured in D(x, y) = I₂ − I₁, in the FITS IMAGE units.")

st.markdown("### Measurements and quality indicators")
quality = candidate.quality
measurement_rows = [
    {"quantity": "Peak SNR of residual", "value": f"{candidate.snr:.2f} σ"},
    {"quantity": "Local difference noise", "value": f"{candidate.local_noise:.4g} {comparison.get('units')}"},
    {"quantity": "Surface-brightness difference at peak", "value": f"{candidate.surface_brightness_diff:.4g} {comparison.get('units')}"},
    {"quantity": "Aperture flux · epoch 1", "value": f"{candidate.flux1_jy * 1e6:.4f} µJy"},
    {"quantity": "Aperture flux · epoch 2", "value": f"{candidate.flux2_jy * 1e6:.4f} µJy"},
    {"quantity": "Flux difference (epoch 2 − epoch 1)", "value": f"{candidate.delta_flux_jy * 1e6:.4f} µJy"},
    {"quantity": "Angular displacement", "value": (f"{candidate.metadata['angular_displacement_arcsec']:.3f} arcsec" if candidate.metadata.get("angular_displacement_arcsec") is not None else "not measurable")},
    {"quantity": "Position angle of displacement", "value": (f"{candidate.metadata['position_angle_deg']:.2f}° east of north" if candidate.metadata.get("position_angle_deg") is not None else "not available")},
    {"quantity": "Time baseline", "value": (f"{candidate.metadata['delta_time_days']:.3f} d" if candidate.metadata.get("delta_time_days") is not None else "unknown")},
    {"quantity": "Pixels above threshold in 3×3 window", "value": str(quality.get("pixels_above_threshold_in_window"))},
    {"quantity": "Residual sign", "value": "epoch 2 brighter" if quality.get("sign", 1) > 0 else "epoch 1 brighter"},
    {"quantity": "Opposite-sign dipole separation", "value": (f"{quality['dipole_separation_px']:.2f} px" if quality.get("dipole_separation_px") is not None else "none detected")},
    {"quantity": "Alignment tolerance", "value": f"{quality.get('alignment_tolerance_arcsec', float('nan')):.2f} arcsec"},
    {"quantity": "Data mode", "value": candidate.data_mode},
]
st.dataframe(measurement_rows, use_container_width=True, hide_index=True)

st.markdown("### Interpretation boundaries")
st.markdown(
    "- A residual is a **measurement**, not a physical object identification. The categories above describe how the measurement relates "
    "to the two epochs; they never assert a planet, asteroid or discovery.\n"
    "- An angular displacement smaller than the alignment tolerance is **not** evidence of motion.\n"
    "- With only two epochs, a single-epoch residual can be a cosmic ray, a detector artifact or a resampling residue; three or more epochs "
    "are required to test persistence and track consistency."
)

if csv_payload:
    st.download_button(
        "EXPORT ALL CANDIDATES (CSV)",
        data=csv_payload,
        file_name=f"candidates_{comparison.get('run_id')}.csv",
        mime="text/csv",
    )
    st.caption(f"CSV provenance: run {comparison.get('run_id')} · epochs {candidate.metadata['epoch1'].get('id')} and {candidate.metadata['epoch2'].get('id')}.")

st.page_link(navbar.page("pages/05_Detection_Engine.py"), label="← BACK TO DETECTION ENGINE", use_container_width=True)
navbar.render_footer()
