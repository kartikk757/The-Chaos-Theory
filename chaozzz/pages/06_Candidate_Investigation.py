"""Candidate investigation workspace: evidence, warnings, cross-match and export."""
from __future__ import annotations

import json
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from astropy.coordinates import SkyCoord

import config
from components import candidate_table, image_viewer, navbar, sidebar, status_panel
from modules import epsf_processor, jpl_crossmatch, mosaic_cube, pipeline, provenance, spherex_sdt_adapter, state, visualization

sidebar.bootstrap_page(
    "Candidate Investigation",
    "Examine an individual residual candidate: position, brightness change, displacement and data quality.",
    "Candidate Analysis",
)

def _cutout(array, x, y, half):
    """Small pixel window around a candidate, clipped to the array bounds."""
    xi, yi = int(round(x)), int(round(y))
    return np.asarray(array)[max(yi - half, 0):yi + half + 1, max(xi - half, 0):xi + half + 1]


def _record_jpl_provenance(comparison: dict, query: str, *, jpl_result: dict | None = None, error: str | None = None) -> None:
    """Append the user-requested JPL hypothesis check to the reproducible run JSON."""
    run_id = comparison.get("run_id")
    if not run_id or run_id == "no-two-epoch-run":
        return
    record = comparison.get("provenance")
    if not isinstance(record, dict) or record.get("run_id") != run_id:
        path = config.OUTPUT_DIR / "runs" / f"{run_id}.json"
        record = provenance.read_json(path) if path.exists() else {"run_id": run_id}
    external = record.setdefault("external_catalogues", {})
    external["JPL SBDB / Horizons"] = {
        "query": str(query).strip(),
        "checked_utc": provenance.utc_now(),
        "result": jpl_result,
        "error": error,
        "scope": "User-specified known-object hypothesis only; not a blind discovery search.",
    }
    provenance.write_json(config.OUTPUT_DIR / "runs" / f"{run_id}.json", record)
    comparison["provenance"] = record
    state.set_value("comparison", comparison)


def _record_epsf_provenance(comparison: dict, candidate_id: str, result: dict) -> None:
    """Append supplemental ePSF fits to the reproducible candidate run JSON."""
    run_id = comparison.get("run_id")
    if not run_id or run_id == "no-two-epoch-run":
        return
    record = comparison.get("provenance")
    if not isinstance(record, dict) or record.get("run_id") != run_id:
        path = config.OUTPUT_DIR / "runs" / f"{run_id}.json"
        record = provenance.read_json(path) if path.exists() else {"run_id": run_id}
    measurements = record.setdefault("ancillary_measurements", {})
    key = "SPHEREx ePSF photometry"
    bucket = measurements.setdefault(key, {})
    if not isinstance(bucket, dict) or not isinstance(bucket.get("by_candidate"), dict):
        previous = bucket if isinstance(bucket, dict) and bucket.get("candidate_id") else None
        bucket = {"by_candidate": {}}
        if previous:
            bucket["by_candidate"][str(previous["candidate_id"])] = previous
        measurements[key] = bucket
    bucket["by_candidate"][str(candidate_id)] = {
        "candidate_id": candidate_id,
        "checked_utc": provenance.utc_now(),
        "result": result,
        "scope": "Supplemental isolated-source fit to the selected candidate; not used to change its classification or to claim a sensitivity improvement.",
    }
    provenance.write_json(config.OUTPUT_DIR / "runs" / f"{run_id}.json", record)
    comparison["provenance"] = record
    state.set_value("comparison", comparison)


def _record_sdt_provenance(comparison: dict, candidate_id: str, result: dict) -> None:
    """Preserve selected-epoch SDT aperture measurements in run provenance."""
    run_id = comparison.get("run_id")
    if not run_id or run_id == "no-two-epoch-run":
        return
    record = comparison.get("provenance")
    if not isinstance(record, dict) or record.get("run_id") != run_id:
        path = config.OUTPUT_DIR / "runs" / f"{run_id}.json"
        record = provenance.read_json(path) if path.exists() else {"run_id": run_id}
    measurements = record.setdefault("ancillary_measurements", {})
    key = "SPHEREx SDT spectral aperture extraction"
    bucket = measurements.setdefault(key, {})
    if not isinstance(bucket, dict) or not isinstance(bucket.get("by_candidate"), dict):
        previous = bucket if isinstance(bucket, dict) and bucket.get("candidate_id") else None
        bucket = {"by_candidate": {}}
        if previous:
            bucket["by_candidate"][str(previous["candidate_id"])] = previous
        measurements[key] = bucket
    bucket["by_candidate"][str(candidate_id)] = {
        "candidate_id": candidate_id,
        "checked_utc": provenance.utc_now(),
        "result": result,
    }
    provenance.write_json(config.OUTPUT_DIR / "runs" / f"{run_id}.json", record)
    comparison["provenance"] = record
    state.set_value("comparison", comparison)


comparison = state.get("comparison")
candidates = (comparison or {}).get("candidates") or []
multi_epoch_result = state.get("multi_epoch_result") or {}
multi_epoch_tracks = multi_epoch_result.get("tracks", []) if multi_epoch_result.get("status") == "ok" else []
if multi_epoch_tracks:
    st.markdown("### Multi-epoch track evidence")
    st.caption("Tracks are sparse-epoch screening results. The fitted local tangent-plane rate and uncertainty floor are not an orbit or a calibrated SPHEREx astrometric solution.")
    track = st.selectbox(
        "Investigate a multi-epoch track",
        multi_epoch_tracks,
        format_func=lambda item: f"{item['track_id']} · {item['status_label']} · {item['detection_count']} detections",
        key="candidate_investigation_track",
    )
    motion = track.get("motion", {})
    st.dataframe([
        {"measurement": "Classification", "value": track.get("status_label")},
        {"measurement": "Independent detections", "value": track.get("detection_count")},
        {"measurement": "Angular rate (arcsec/day)", "value": motion.get("speed_arcsec_per_day")},
        {"measurement": "Rate uncertainty (arcsec/day)", "value": motion.get("speed_sigma_arcsec_per_day")},
        {"measurement": "Reduced chi-square", "value": motion.get("reduced_chi2")},
        {"measurement": "Quality flags", "value": ", ".join(track.get("quality_flags", []))},
    ], use_container_width=True, hide_index=True)
    if track.get("detections"):
        st.dataframe(pd.DataFrame(track["detections"]), use_container_width=True, hide_index=True)
    for evidence in track.get("evidence", []):
        st.write(f"- {evidence}")
    for limitation in track.get("limitations", []):
        st.caption(f"Limit: {limitation}")
    if track.get("quality_flags") and "2-EPOCH-ONLY" in track["quality_flags"]:
        st.warning("2-EPOCH-ONLY: two positions give a displacement, not a fitted motion track.")
    st.markdown("---")

if not comparison or not candidates:
    status_panel.render_empty_state(
        "No candidates to investigate",
        "Run the detection pipeline on a compatible pair, or review the multi-epoch track evidence above. A null result is preserved as a valid outcome.",
        "pages/05_Detection_Engine.py",
        "OPEN DETECTION ENGINE",
    )
    navbar.render_footer()
    st.stop()

aligned = comparison.get("aligned")
difference = comparison.get("difference")
crossmatch = comparison.get("crossmatch", {})
gaia_crossmatch = comparison.get("gaia_crossmatch", {})
product1 = comparison.get("product1")
product2 = comparison.get("product2")

csv_payload = candidate_table.render_candidate_table(candidates, crossmatch, "invest")
selected_id = st.session_state.get("selected_candidate_id")
candidate = next((item for item in candidates if item.candidate_id == selected_id), candidates[0])

st.markdown("---")
candidate_table.render_candidate_detail(candidate, crossmatch.get(candidate.candidate_id), "invest_detail")

st.markdown("### Gaia DR3 epoch-aware cross-match")
gaia_result = gaia_crossmatch.get(candidate.candidate_id, {})
gaia_diagnostics = comparison.get("gaia_diagnostics", {})
st.caption(
    "Gaia DR3 reference-epoch positions are propagated to the candidate's epoch1 observation time before matching. "
    "The radius combines the Gaia propagated positional uncertainty with a conservative half-pixel SPHEREx position floor."
)
if gaia_result.get("status") == "matched":
    st.success(gaia_result.get("detail", "Gaia DR3 catalogue match found."))
    st.dataframe(pd.DataFrame(gaia_result.get("matches", [])), use_container_width=True, hide_index=True)
elif gaia_result.get("status") == "no match":
    st.info(gaia_result.get("detail", "No Gaia DR3 counterpart within the uncertainty-scaled match radius."))
elif gaia_result:
    st.warning(gaia_result.get("detail", "Gaia DR3 cross-match unavailable."))
else:
    st.info(gaia_diagnostics.get("detail", "Gaia DR3 matching was disabled or no residual candidate was available."))
if gaia_diagnostics:
    with st.expander("Gaia DR3 query and matching provenance", expanded=False):
        st.json(gaia_diagnostics)

st.markdown("### 2MASS PSC and AllWISE evidence")
irsa_result = crossmatch.get(candidate.candidate_id, {})
catalogue_diagnostics = comparison.get("catalogue_diagnostics", {})
if irsa_result.get("status") == "matched":
    st.success(irsa_result.get("detail", "At least one IRSA catalogue counterpart is within the uncertainty-scaled radius."))
    st.caption("Matches are ordered by separation relative to combined positional uncertainty. Catalogue positions have not been epoch-propagated; AllWISE apparent-motion fields are diagnostic only and may include parallax.")
    match_rows = irsa_result.get("matches", [])
    if match_rows:
        st.dataframe(pd.DataFrame(match_rows), use_container_width=True, hide_index=True)
        st.caption(f"Alternative counterparts retained: {irsa_result.get('all_match_count', len(match_rows))}. Multiple matches should be reviewed rather than collapsed to one identity.")
elif irsa_result:
    st.info(irsa_result.get("detail", "No catalogue counterpart was returned; this is not evidence of an unknown object."))
else:
    st.info(catalogue_diagnostics.get("detail", "IRSA catalogue matching was disabled or unavailable for this run."))
if catalogue_diagnostics:
    with st.expander("IRSA catalogue query and provenance", expanded=False):
        st.json(catalogue_diagnostics)

st.markdown("### JPL known-object hypothesis check")
st.caption(
    "Enter an already known asteroid or comet name, designation, or SPK ID. SBDB resolves the object and Horizons predicts its position "
    "at the measured candidate epochs. This is a known-object investigation only; Horizons is not a blind discovery search."
)
with st.form("jpl_known_object_form", clear_on_submit=False):
    jpl_cols = st.columns([2, 1, 1])
    with jpl_cols[0]:
        known_object = st.text_input("Known-object identifier", key="jpl_known_object_query", placeholder="e.g. a designation or SPK ID")
    with jpl_cols[1]:
        jpl_radius = st.number_input("Match radius (arcsec)", min_value=0.1, max_value=3600.0, value=60.0, step=5.0, key="jpl_match_radius")
    with jpl_cols[2]:
        st.markdown("<div style='height:1.55rem'></div>", unsafe_allow_html=True)
        jpl_submitted = st.form_submit_button("CHECK KNOWN OBJECT", use_container_width=True, type="primary")
if jpl_submitted:
    if not known_object.strip():
        st.warning("Enter a known-object identifier to run SBDB and Horizons.")
    else:
        try:
            with st.spinner("Resolving the object in JPL SBDB and requesting Horizons positions at the candidate epochs…"):
                jpl_result = jpl_crossmatch.investigate_known_object(known_object, [candidate], float(jpl_radius))
            state.set_value("jpl_candidate_result", {"candidate_id": candidate.candidate_id, "query": known_object, "result": jpl_result})
            _record_jpl_provenance(comparison, known_object, jpl_result=jpl_result)
        except Exception as error:
            state.set_value("jpl_candidate_result", {"candidate_id": candidate.candidate_id, "query": known_object, "error": str(error)})
            _record_jpl_provenance(comparison, known_object, error=str(error))
jpl_saved = state.get("jpl_candidate_result") or {}
if jpl_saved.get("candidate_id") == candidate.candidate_id and jpl_saved.get("query") == known_object:
    if jpl_saved.get("error"):
        st.warning(f"JPL investigation unavailable: {jpl_saved['error']}")
    else:
        jpl_result = jpl_saved.get("result", {})
        if jpl_result.get("status") == "ok":
            summary = jpl_result.get("summary", {})
            if summary.get("hypothesis_consistent_with_positions"):
                st.success(f"The known-object hypothesis is positionally consistent at {summary.get('epochs_within_radius')} of {summary.get('epochs_compared')} epoch(s). Minimum separation: {summary.get('minimum_separation_arcsec'):.2f} arcsec.")
            else:
                st.info(f"No candidate epoch fell within the selected match radius. Minimum separation: {summary.get('minimum_separation_arcsec'):.2f} arcsec. This does not imply an unknown object.")
            st.dataframe(pd.DataFrame(jpl_result.get("matches", [])), use_container_width=True, hide_index=True)
            with st.expander("SBDB object/orbit data and Horizons provenance", expanded=False):
                st.json({"sbdb": jpl_result.get("sbdb"), "ephemeris": jpl_result.get("ephemeris"), "limitations": jpl_result.get("limitations")})
        elif jpl_result.get("status") == "insufficient_candidate_epochs":
            st.info(jpl_result.get("detail", "Candidate position epochs are insufficient for a Horizons comparison."))

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

st.markdown("### Candidate-centered spectral quick-look")
cube = state.get("spectral_cube")
if cube is None:
    st.info("No validated spectral cube is loaded. Upload a mosaic or HiPS cube in Cosmic Image Lab to associate a quick-look spectrum with this candidate.")
else:
    st.caption(
        f"Using candidate ICRS position ({candidate.ra_deg:.7f}°, {candidate.dec_deg:.7f}°) in {cube.source_name}. "
        "This cube spectrum is a spatial/spectral measurement, separate from the time-domain candidate evidence."
    )
    aperture_cols = st.columns(3)
    candidate_aperture = aperture_cols[0].number_input("Candidate aperture radius (px)", 0.5, 10.0, 2.0, 0.5, key="candidate_cube_aperture")
    candidate_annulus_inner = aperture_cols[1].number_input("Background annulus inner radius (px)", 1.0, 20.0, 5.0, 0.5, key="candidate_cube_annulus_inner")
    candidate_annulus_outer = aperture_cols[2].number_input("Background annulus outer radius (px)", float(candidate_annulus_inner + 0.5), 30.0, max(8.0, float(candidate_annulus_inner + 1.0)), 0.5, key="candidate_cube_annulus_outer")
    if st.button("EXTRACT CANDIDATE SPECTRUM", key="extract_candidate_cube_spectrum", type="primary"):
        try:
            spectrum = mosaic_cube.extract_aperture_spectrum(
                cube, float(candidate.ra_deg), float(candidate.dec_deg),
                float(candidate_aperture), float(candidate_annulus_inner), float(candidate_annulus_outer),
            )
            state.set_value("candidate_cube_spectrum", {"candidate_id": candidate.candidate_id, "cube_sha256": cube.sha256, "spectrum": spectrum})
        except Exception as error:
            st.error(f"Candidate-centered spectral extraction failed: {error}")
    saved_spectrum = state.get("candidate_cube_spectrum") or {}
    if saved_spectrum.get("candidate_id") == candidate.candidate_id and saved_spectrum.get("cube_sha256") == cube.sha256:
        spectrum = saved_spectrum["spectrum"]
        spectrum_figure = go.Figure(go.Scatter(x=spectrum["wavelength_um"], y=spectrum["flux"], mode="lines+markers", name="Candidate aperture − background"))
        spectrum_figure.update_layout(title=f"{candidate.candidate_id} · cube quick-look spectrum", xaxis_title="Wavelength (µm)", yaxis_title=spectrum["flux_unit"], template="plotly_dark", height=380)
        st.plotly_chart(spectrum_figure, use_container_width=True, config={"displaylogo": False})
        st.caption(spectrum["extraction"])
        with st.expander("Candidate spectrum provenance and limitations", expanded=False):
            st.json({key: value for key, value in spectrum.items() if key not in ("wavelength_um", "flux", "background_per_pixel", "flux_error", "valid_fraction")})

with st.expander("Optional IRSA Source Discovery Tool spectral photometry", expanded=False):
    sdt_spectrum_state = spherex_sdt_adapter.spectrum_availability()
    if sdt_spectrum_state["status"] != "available":
        st.info(sdt_spectrum_state["detail"] + " The cube quick-look path above remains independent.")
    elif aligned is None or product1 is None or product2 is None:
        st.info("Load a compatible pair before extracting SDT measurements from selected epochs.")
    else:
        st.caption(
            "This uses the official SDT aperture, background-annulus, solid-angle-map, flag and AB-magnitude components on the exact loaded FITS products. "
            "Each epoch uses its own candidate position. Results are sparse, time-stamped flux measurements, not a variability-corrected continuous spectrum."
        )
        sdt_radius = st.selectbox("SDT aperture radius (pixels)", [3.0, 4.0, 5.0], index=0, key="candidate_sdt_aperture")
        if st.button("EXTRACT SDT SPECTRAL PHOTOMETRY FOR BOTH EPOCHS", key=f"extract_sdt_{candidate.candidate_id}"):
            products = [product1, product2]
            positions = []
            for index, product in enumerate(products, start=1):
                epoch = candidate.metadata.get(f"epoch{index}", {})
                ra_value = epoch.get("ra_deg", candidate.ra_deg)
                dec_value = epoch.get("dec_deg", candidate.dec_deg)
                try:
                    if product.wcs is None or not np.isfinite(float(ra_value)) or not np.isfinite(float(dec_value)):
                        raise ValueError("A finite ICRS candidate coordinate and product WCS are required.")
                    positions.append((float(ra_value), float(dec_value)))
                except Exception as exc:
                    positions.append((float("nan"), float("nan")))
            spectrum = spherex_sdt_adapter.extract_epoch_spectrum(
                products,
                positions,
                source_id=candidate.candidate_id,
                aperture_radius_px=float(sdt_radius),
            )
            candidate.metadata["sdt_spectral_photometry"] = spectrum
            comparison["candidates"] = candidates
            state.set_value("comparison", comparison)
            _record_sdt_provenance(comparison, candidate.candidate_id, spectrum)

        sdt_spectrum = candidate.metadata.get("sdt_spectral_photometry") or {}
        if sdt_spectrum:
            if sdt_spectrum.get("status") not in {"ok", "partial"}:
                st.warning(sdt_spectrum.get("detail", "The SDT extraction did not produce measurements."))
            measurements = sdt_spectrum.get("measurements", [])
            valid_measurements = [row for row in measurements if row.get("status") == "ok"]
            if valid_measurements:
                st.dataframe([
                    {
                        "observation": row.get("observation_id"),
                        "UTC date": row.get("date_utc"),
                        "wavelength (µm)": row.get("wavelength_um"),
                        "bandwidth (µm)": row.get("bandwidth_um"),
                        "flux": f"{row.get('flux'):.5g} {row.get('flux_unit')}",
                        "uncertainty": f"{row.get('flux_error'):.4g} {row.get('flux_unit')}",
                        "AB magnitude": row.get("ab_magnitude"),
                        "flagged aperture pixels": row.get("bad_flagged_aperture_pixels"),
                    }
                    for row in valid_measurements
                ], use_container_width=True, hide_index=True)
                units = {row.get("flux_unit") for row in valid_measurements}
                if len(valid_measurements) >= 2 and len(units) == 1:
                    ordered = sorted(valid_measurements, key=lambda row: (row["wavelength_um"], row["mjd_utc"]))
                    spectrum_figure = go.Figure(go.Scatter(
                        x=[row["wavelength_um"] for row in ordered],
                        y=[row["flux"] for row in ordered],
                        error_y={"type": "data", "array": [row["flux_error"] for row in ordered], "visible": True},
                        customdata=[[row["date_utc"], row["observation_id"]] for row in ordered],
                        hovertemplate="%{x:.4f} µm<br>%{y:.5g} " + next(iter(units)) + "<br>%{customdata[0]}<br>%{customdata[1]}<extra></extra>",
                        mode="markers+lines",
                        name="SDT epoch measurements",
                    ))
                    spectrum_figure.update_layout(
                        title=f"{candidate.candidate_id} · sparse SDT epoch measurements",
                        xaxis_title="Wavelength (µm)",
                        yaxis_title=next(iter(units)),
                        template="plotly_dark",
                        height=380,
                    )
                    st.plotly_chart(spectrum_figure, use_container_width=True, config={"displaylogo": False})
            for row in measurements:
                if row.get("status") != "ok":
                    st.warning(f"{row.get('observation_id', 'Epoch')}: {row.get('detail', row.get('status'))}")
            st.caption(sdt_spectrum.get("limitation", "SDT spectral measurement limitations are unavailable."))
            with st.expander("SDT extraction provenance", expanded=False):
                st.json({key: value for key, value in sdt_spectrum.items() if key != "measurements"})
                st.download_button(
                    "DOWNLOAD SDT SPECTRAL PHOTOMETRY (JSON)",
                    data=json.dumps(provenance.safe_json(sdt_spectrum), indent=2),
                    file_name=f"sdt_spectrum_{candidate.candidate_id}_{comparison.get('run_id', 'run')}.json",
                    mime="application/json",
                    key=f"download_sdt_spectrum_{candidate.candidate_id}",
                )

with st.expander("SPHEREx ePSF diagnostics and supplemental fit", expanded=False):
    st.caption(
        "The official EPSF table is selected by mapping the displayed sky coordinate through each FITS WCS and back to the parent detector zone. "
        "The original candidate measurements remain aperture-based. An optional single-source fit below uses Photutils ImagePSF with the FITS-provided normalization and is supplemental evidence only."
    )
    if aligned is None or product1 is None or product2 is None or aligned.target_wcs is None:
        st.info("Load a compatible pair before inspecting ePSF zones.")
    else:
        if st.button("FIT SELECTED CANDIDATE WITH THE ePSF IN BOTH EPOCHS", key=f"fit_epsf_{candidate.candidate_id}"):
            epoch_results = []
            for index, product in enumerate((product1, product2), start=1):
                epoch_key = f"epoch{index}"
                epoch_info = candidate.metadata.get(epoch_key, {})
                ra_value = epoch_info.get("ra_deg", candidate.ra_deg)
                dec_value = epoch_info.get("dec_deg", candidate.dec_deg)
                sky_position_deg = None
                try:
                    if product.wcs is None or not np.isfinite(float(ra_value)) or not np.isfinite(float(dec_value)):
                        raise ValueError("A finite sky position and celestial product WCS are required.")
                    ra_value, dec_value = float(ra_value), float(dec_value)
                    sky_position_deg = [ra_value, dec_value]
                    sky_position = SkyCoord(ra_value, dec_value, unit="deg", frame="icrs")
                    x_native, y_native = product.wcs.world_to_pixel(sky_position)
                    fit = epsf_processor.fit_source(product, float(x_native), float(y_native))
                except Exception as exc:
                    fit = {"status": "unavailable", "detail": str(exc)}
                epoch_results.append({
                    "epoch": index,
                    "observation_id": epoch_info.get("id"),
                    "date_utc": epoch_info.get("date_utc"),
                    "matched_source": epoch_info.get("matched_source"),
                    "sky_position_icrs_deg": sky_position_deg,
                    "fit": fit,
                })
            fit_report = {
                "status": "ok" if all(row["fit"].get("status") == "ok" for row in epoch_results) else "partial",
                "run_id": comparison.get("run_id"),
                "candidate_id": candidate.candidate_id,
                "epochs": epoch_results,
                "limitation": "Supplemental isolated-source fits do not update aperture classifications, deblend neighbors, calibrate completeness, or establish a sensitivity gain.",
            }
            candidate.metadata["epsf_photometry"] = fit_report
            comparison["candidates"] = candidates
            state.set_value("comparison", comparison)
            _record_epsf_provenance(comparison, candidate.candidate_id, fit_report)

        fit_report = candidate.metadata.get("epsf_photometry") or {}
        if fit_report.get("candidate_id") == candidate.candidate_id:
            for epoch_row in fit_report.get("epochs", []):
                fit = epoch_row.get("fit", {})
                st.markdown(f"**Epoch {epoch_row.get('epoch')} · {epoch_row.get('observation_id') or 'observation unavailable'}**")
                if fit.get("status") == "ok":
                    st.dataframe([
                        {"measurement": "Fitted ePSF amplitude", "value": f"{fit['flux_fit']:.6g} {fit['flux_unit']}"},
                        {"measurement": "Fit uncertainty", "value": (f"{fit['flux_error']:.4g} {fit['flux_unit']}" if fit.get("flux_error") is not None else "not available")},
                        {"measurement": "Fitted native pixel position", "value": [fit["x_fit"], fit["y_fit"]]},
                        {"measurement": "Reduced chi-square", "value": fit["reduced_chi_square"]},
                        {"measurement": "ePSF zone / oversampling", "value": [fit["epsf_zone_id"], fit["epsf_oversampling_xy"]]},
                        {"measurement": "Wavelength delta (µm)", "value": fit.get("wavelength_delta_um")},
                        {"measurement": "Normalization ratio", "value": fit.get("epsf_normalization_ratio")},
                    ], use_container_width=True, hide_index=True)
                    if fit.get("warnings"):
                        st.caption(" · ".join(fit["warnings"]))
                else:
                    st.info(f"ePSF fit {fit.get('status', 'unavailable')}: {fit.get('detail', 'No compatible measurement was produced.')}")
            st.download_button(
                "DOWNLOAD ePSF FIT REPORT (JSON)",
                data=json.dumps(provenance.safe_json(fit_report), indent=2),
                file_name=f"epsf_{candidate.candidate_id}_{comparison.get('run_id', 'run')}.json",
                mime="application/json",
                key=f"download_epsf_{candidate.candidate_id}",
            )

        candidate_sky = aligned.target_wcs.pixel_to_world(candidate.pixel_x, candidate.pixel_y)
        epsf_columns = st.columns(2)
        for index, (column, product) in enumerate(zip(epsf_columns, (product1, product2)), start=1):
            with column:
                st.markdown(f"**Epoch {index} · {product.metadata.get('observation', {}).get('obs_id', product.product_id)}**")
                if product.wcs is None:
                    st.info("Celestial WCS is unavailable, so no source-position ePSF zone can be selected.")
                    continue
                try:
                    x_native, y_native = product.wcs.world_to_pixel(candidate_sky)
                    epsf_info = epsf_processor.epsf_for_position(product, float(x_native), float(y_native))
                except Exception as exc:
                    st.warning(f"ePSF inspection unavailable: {exc}")
                    continue
                if epsf_info.get("status") != "ok":
                    st.info(epsf_info.get("detail", "No usable ePSF zone is available."))
                    continue
                st.plotly_chart(
                    visualization.cutout_figure(
                        epsf_info["epsf_normalized"],
                        f"Selected oversampled ePSF · zone {epsf_info['zone_id']}",
                    ),
                    use_container_width=True, theme=None, config={"displaylogo": False},
                )
                st.dataframe([
                    {"diagnostic": "Parent detector pixel", "value": epsf_info["parent_detector_pixel"]},
                    {"diagnostic": "ePSF zone [BINX, BINY]", "value": epsf_info["zone_id"]},
                    {"diagnostic": "ePSF model grid", "value": epsf_info["epsf_shape_oversampled"]},
                    {"diagnostic": "Oversampling [x, y]", "value": [epsf_info["oversampling_x"], epsf_info["oversampling_y"]]},
                    {"diagnostic": "Zone CWAVE / CBAND (µm)", "value": [epsf_info["epsf_center_wavelength_um"], epsf_info["epsf_bandwidth_um"]]},
                    {"diagnostic": "Image wavelength at position (µm)", "value": epsf_info["image_wavelength_at_source_um"]},
                    {"diagnostic": "ePSF effective area", "value": epsf_info["epsf_effective_area"]},
                    {"diagnostic": "ePSF finite-pixel sum", "value": epsf_info["epsf_sum_raw"]},
                ], use_container_width=True, hide_index=True)
                if epsf_info.get("warnings"):
                    st.warning("\n".join(epsf_info["warnings"]))
                st.caption(epsf_info["interpretation"])

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
