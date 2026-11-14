"""Candidate table, inspection panel and CSV export widgets."""
from __future__ import annotations

from typing import Any

import pandas as pd
import streamlit as st

from modules import pipeline


def candidate_dataframe(candidates: list[Any], crossmatch: dict | None = None) -> pd.DataFrame:
    rows = []
    for candidate in candidates:
        match = (crossmatch or {}).get(candidate.candidate_id, {})
        rows.append({
            "Candidate": candidate.candidate_id,
            "Classification": candidate.classification,
            "SNR": round(candidate.snr, 2),
            "RA (deg)": round(candidate.ra_deg, 6) if candidate.ra_deg == candidate.ra_deg else None,
            "Dec (deg)": round(candidate.dec_deg, 6) if candidate.dec_deg == candidate.dec_deg else None,
            "x (px)": round(candidate.pixel_x, 2),
            "y (px)": round(candidate.pixel_y, 2),
            "Δ flux (µJy)": round(candidate.delta_flux_jy * 1e6, 4) if candidate.delta_flux_jy == candidate.delta_flux_jy else None,
            "Displacement (arcsec)": (
                round(candidate.metadata["angular_displacement_arcsec"], 3)
                if candidate.metadata.get("angular_displacement_arcsec") is not None else None
            ),
            "Cross-match": match.get("status", "not attempted"),
            "Warnings": len(candidate.warnings),
        })
    return pd.DataFrame(rows)


def render_candidate_table(candidates: list[Any], crossmatch: dict | None, key: str) -> str | None:
    if not candidates:
        st.info("No candidate records exist for this comparison. Adjust the detection threshold or choose a different epoch pair.")
        return None
    frame = candidate_dataframe(candidates, crossmatch)
    st.dataframe(frame, use_container_width=True, hide_index=True, key=f"{key}_table")
    identifiers = [candidate.candidate_id for candidate in candidates]
    labels = {candidate.candidate_id: f"{candidate.candidate_id} · {candidate.classification} · {candidate.snr:.1f}σ" for candidate in candidates}
    default = st.session_state.get("selected_candidate_id")
    index = identifiers.index(default) if default in identifiers else 0
    choice = st.selectbox("Select candidate to investigate", identifiers, index=index, format_func=lambda value: labels[value], key=f"{key}_select")
    st.session_state["selected_candidate_id"] = choice
    rows = [pipeline._candidate_row(candidate, (crossmatch or {}).get(candidate.candidate_id)) for candidate in candidates]
    return pipeline.candidates_to_csv(rows)


def render_candidate_detail(candidate: Any, crossmatch: dict | None, key: str) -> None:
    st.markdown(f"### {candidate.candidate_id} · {candidate.classification}")
    metrics = st.columns(4)
    metrics[0].metric("Peak SNR", f"{candidate.snr:.1f}σ")
    metrics[1].metric("Δ flux", f"{candidate.delta_flux_jy * 1e6:.4f} µJy", help="Aperture flux difference, epoch 2 minus epoch 1")
    metrics[2].metric(
        "Angular displacement",
        f"{candidate.metadata['angular_displacement_arcsec']:.2f} arcsec"
        if candidate.metadata.get("angular_displacement_arcsec") is not None else "not measurable",
    )
    metrics[3].metric(
        "Δt",
        f"{candidate.metadata['delta_time_days']:.2f} d" if candidate.metadata.get("delta_time_days") is not None else "unknown",
    )
    position = st.columns(4)
    position[0].metric("RA", f"{candidate.ra_deg:.6f}°" if candidate.ra_deg == candidate.ra_deg else "unavailable")
    position[1].metric("Dec", f"{candidate.dec_deg:.6f}°" if candidate.dec_deg == candidate.dec_deg else "unavailable")
    position[2].metric("Pixel", f"{candidate.pixel_x:.2f}, {candidate.pixel_y:.2f}")
    position[3].metric("Local Δ noise", f"{candidate.local_noise:.4g}", help="Robust local noise of the difference image")

    epochs = candidate.metadata
    st.markdown(
        f"<div class='ct-mono-note'>Epoch 1 · {epochs['epoch1'].get('date_utc')} · {epochs['epoch1'].get('band')} · "
        f"{epochs['epoch1'].get('origin')}<br/>Epoch 2 · {epochs['epoch2'].get('date_utc')} · {epochs['epoch2'].get('band')} · "
        f"{epochs['epoch2'].get('origin')}</div>",
        unsafe_allow_html=True,
    )
    with st.expander("WHY WAS THIS CANDIDATE FLAGGED?", expanded=True):
        for reason in candidate.reasons:
            st.markdown(f"- {reason}")
        st.caption("These statements are generated directly from the measurements listed above. They describe screening criteria, not a physical interpretation.")
    if candidate.warnings:
        st.warning("Processing warnings for this candidate:\n\n" + "\n".join(f"- {warning}" for warning in candidate.warnings))
    quality = st.columns(3)
    quality[0].metric("Severe flags at position", f"{candidate.quality['severe_flag_fraction_epoch1']:.0%} / {candidate.quality['severe_flag_fraction_epoch2']:.0%}", help="Epoch 1 / epoch 2 fraction of the 3×3 neighbourhood carrying severe FLAGS bits")
    quality[1].metric("Epoch detections searched", f"{candidate.quality['epoch_sources_found'][0]} / {candidate.quality['epoch_sources_found'][1]}", help="Compact sources detected independently in epoch 1 / epoch 2 images")
    quality[2].metric("Alignment tolerance", f"{candidate.quality['alignment_tolerance_arcsec']:.2f} arcsec", help="1.5 × the reprojected pixel scale; offsets below this are not evidence of motion")
    if crossmatch:
        status = crossmatch.get("status")
        if status == "matched":
            st.success(f"Catalogue cross-match: {crossmatch['detail']}")
            st.dataframe(pd.DataFrame(crossmatch["matches"]), use_container_width=True, hide_index=True)
        elif status == "no match":
            st.info(f"Catalogue cross-match: {crossmatch['detail']} An absent catalogue entry does not imply a new object; the catalogue may be shallower or at a different wavelength.")
        else:
            st.warning(f"Catalogue cross-match unavailable: {crossmatch.get('detail', '')}")
