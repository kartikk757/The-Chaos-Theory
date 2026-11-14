"""Browse retrieved observations and select a comparable pair."""
from __future__ import annotations

import pandas as pd
import streamlit as st

import config
from components import metric_cards, navbar, sidebar, status_panel
from modules import data_loader, state, wcs_alignment

sidebar.bootstrap_page(
    "Observation Explorer",
    "Inspect retrieved SPHEREx observations and choose two that can be compared honestly.",
    "Observations",
)

observations = state.get("observations", [])
if not observations:
    status_panel.render_empty_state(
        "No observation metadata loaded",
        "Query IRSA in Sky Explorer, or load the verified cached pair / demonstration data. This page never fabricates observation rows.",
        "pages/02_Sky_Explorer.py",
        "OPEN SKY EXPLORER",
    )
    navbar.render_footer()
    st.stop()

bands = sorted({str(obs.get("band")) for obs in observations})
releases = sorted({str(obs.get("release")) for obs in observations})
filters = st.columns([1.3, 1, 1, 1])
with filters[0]:
    band_filter = st.multiselect("Band", bands, default=bands)
with filters[1]:
    release_filter = st.multiselect("Release", releases, default=releases)
with filters[2]:
    dated_only = st.toggle("Dated observations only", value=True)
with filters[3]:
    cached_only = st.toggle("Cached files only", value=False, help="Show only observations whose FITS file already exists locally.")

visible = [
    obs for obs in observations
    if str(obs.get("band")) in band_filter
    and str(obs.get("release")) in release_filter
    and (obs.get("date_utc") or not dated_only)
    and (obs.get("cached_path") or obs.get("cutout_url") or not cached_only)
]

metric_cards.render_metric_grid([
    {"label": "Observations (filtered)", "value": str(len(visible)), "caption": f"of {len(observations)} archive records in session", "icon": "☷", "accent": True},
    {"label": "Earliest epoch", "value": min((obs.get("date_utc") or "9999") for obs in visible)[:19] if visible else "Awaiting data", "caption": "UTC, from TAP time_bounds_lower", "icon": "◷"},
    {"label": "Latest epoch", "value": max((obs.get("date_utc") or "") for obs in visible)[:19] if visible else "Awaiting data", "caption": "UTC, from TAP time_bounds_lower", "icon": "◶"},
    {"label": "Data mode", "value": state.data_mode_label(), "caption": state.get("mode_note", "")[:70], "icon": "◈"},
])

table = pd.DataFrame([
    {
        "observation id": obs.get("obs_id"),
        "date (UTC)": obs.get("date_utc"),
        "band": obs.get("band"),
        "release": obs.get("release"),
        "Record RA (deg)": obs.get("ra_deg"),
        "Record Dec (deg)": obs.get("dec_deg"),
        "plane": obs.get("plane_id"),
        "quality": obs.get("quality_flag") or "not supplied",
        "file": "cached" if obs.get("cached_path") else ("remote cutout" if obs.get("cutout_url") else "no FITS URL"),
    }
    for obs in visible
])
st.dataframe(table, use_container_width=True, hide_index=True)
st.caption(
    "Record coordinates are the parent observation footprint centre when matched by SIA2, or the requested field centre otherwise. "
    "They are separate from each downloaded cutout's geometric centre, which is measured from that FITS file's WCS. "
    "The 'coordinates from' field in Sky Explorer identifies the record-coordinate source."
)

st.markdown("### Select two observations for comparison")
labels = {
    obs["id"]: f"{obs.get('date_utc') or 'undated'} · {obs.get('band')} · {obs.get('obs_id')} · {obs.get('release')}"
    for obs in visible
}
identifiers = list(labels)
if len(identifiers) < 2:
    st.warning("At least two observations are required. Adjust the filters or query a wider date range.")
    navbar.render_footer()
    st.stop()

default_pair = state.get("selected_pair") or identifiers[:2]
default_pair = [value for value in default_pair if value in identifiers][:2] or identifiers[:2]
selection = st.multiselect(
    "Comparison pair (select exactly two)",
    identifiers,
    default=default_pair,
    format_func=lambda value: labels[value],
    max_selections=2,
    help="Only two epochs are compared at a time in this prototype.",
)

default_epochs = state.get("selected_epoch_ids", [])
default_epochs = [value for value in default_epochs if value in identifiers]
default_epochs = list(dict.fromkeys([*default_pair, *default_epochs]))
epoch_selection = st.multiselect(
    "Epoch set for multi-epoch investigation (2–20 observations)",
    identifiers,
    default=default_epochs,
    format_func=lambda value: labels[value],
    max_selections=20,
    key="multi_epoch_selection",
    help="The engine sorts by UTC observation time and evaluates every pair. Motion tracks require at least three detections across compatible epochs.",
)

columns = st.columns([1, 2])
with columns[0]:
    load_and_check = st.button("LOAD FITS & CHECK COMPATIBILITY", type="primary", use_container_width=True, disabled=len(selection) != 2)
with columns[1]:
    st.caption(
        "Loading downloads the documented IRSA cutouts (a few MB each), validates the FITS headers and WCS, measures the WCS-WAVE central "
        "wavelength at the cutout centre, and reports incompatibilities before any image arithmetic is allowed."
    )

if load_and_check and len(selection) == 2:
    record1 = next(obs for obs in observations if obs["id"] == selection[0])
    record2 = next(obs for obs in observations if obs["id"] == selection[1])
    try:
        with st.spinner("Retrieving and validating FITS products from IRSA…"):
            product1 = data_loader.load_observation_product(record1)
        with st.spinner("Retrieving the second epoch…"):
            product2 = data_loader.load_observation_product(record2)
        products = dict(state.get("fits_products", {}))
        products[selection[0]] = product1
        products[selection[1]] = product2
        state.set_value("fits_products", products)
        state.set_products({selection[0]: {"path": product1.path, "sha256": product1.sha256}, selection[1]: {"path": product2.path, "sha256": product2.sha256}})
        state.set_selected_pair([selection[0], selection[1]])
        check = wcs_alignment.check_observation_compatibility(
            record1, record2, product1, product2, state.get("processing_config", config.PROCESSING_DEFAULTS)
        )
        state.set_value("pair_check", check)
        state.add_log("Compatibility evaluation completed", level="info" if check["compatible"] else "warning", details={
            "compatible": check["compatible"],
            "reasons": check["reasons"],
            "epochs": [record1.get("obs_id"), record2.get("obs_id")],
        })
    except Exception as error:
        status_panel.render_exception(error, "FITS retrieval / validation")

st.caption(
    "Loading an epoch set retrieves and validates each selected FITS file. The separate two-observation comparison pair remains available "
    "for the existing difference workflow; a 2-epoch analysis is always labelled 2-EPOCH-ONLY."
)
if st.button(
    "LOAD AND VALIDATE SELECTED EPOCH SET",
    key="load_multi_epoch_set",
    use_container_width=True,
    disabled=len(epoch_selection) < 2,
):
    products = dict(state.get("fits_products", {}))
    product_manifest = dict(state.get("products", {}))
    try:
        with st.spinner(f"Retrieving and validating {len(epoch_selection)} selected epoch(s)…"):
            for identifier in epoch_selection:
                if identifier not in products:
                    record = next(row for row in observations if row["id"] == identifier)
                    products[identifier] = data_loader.load_observation_product(record)
                product_manifest[identifier] = {"path": products[identifier].path, "sha256": products[identifier].sha256}
        state.set_value("fits_products", products)
        state.set_products(product_manifest)
        state.set_value("selected_epoch_ids", list(epoch_selection))
        if len(selection) == 2:
            pair_records = [next(row for row in observations if row["id"] == identifier) for identifier in selection]
            state.set_selected_pair(list(selection))
            pair_check = wcs_alignment.check_observation_compatibility(
                pair_records[0], pair_records[1], products[selection[0]], products[selection[1]],
                state.get("processing_config", config.PROCESSING_DEFAULTS),
            )
            state.set_value("pair_check", pair_check)
        st.success(f"Validated {len(epoch_selection)} epoch FITS product(s). Each file checksum is available in run provenance.")
        st.rerun()
    except Exception as error:
        status_panel.render_exception(error, "Multi-epoch FITS retrieval / validation")

products = state.get("fits_products", {})
pair = state.get("selected_pair", [])
check = state.get("pair_check")

if pair and all(identifier in products for identifier in pair):
    product1, product2 = products[pair[0]], products[pair[1]]
    st.markdown("### Loaded products")
    metric_cards.render_metric_grid([
        {"label": "Epoch 1 array", "value": f"{product1.shape[0]}×{product1.shape[1]}", "caption": f"{product1.bunit or 'BUNIT not supplied'} · {product1.center_wavelength_um:.5f} µm" if product1.center_wavelength_um else f"{product1.bunit or 'BUNIT not supplied'} · wavelength unavailable", "icon": "▤"},
        {"label": "Epoch 2 array", "value": f"{product2.shape[0]}×{product2.shape[1]}", "caption": f"{product2.bunit or 'BUNIT not supplied'} · {product2.center_wavelength_um:.5f} µm" if product2.center_wavelength_um else f"{product2.bunit or 'BUNIT not supplied'} · wavelength unavailable", "icon": "▥"},
        {"label": "Pixel scale", "value": f"{product1.pixel_scale_arcsec:.2f}″/px" if product1.pixel_scale_arcsec else "unavailable", "caption": "From the FITS celestial WCS", "icon": "⌗"},
        {"label": "Severe FLAGS", "value": f"{product1.metadata['flag_diagnostics'].get('bad_pixel_fraction', 0):.2%}", "caption": "Epoch 1 pixels with named defect bits (SOURCE excluded)", "icon": "⚑"},
    ])
    with st.expander("Epoch detail and FITS validation", expanded=False):
        for index, (identifier, product) in enumerate(((pair[0], product1), (pair[1], product2)), start=1):
            record = next(obs for obs in observations if obs["id"] == identifier)
            st.markdown(f"**Epoch {index} · {record.get('obs_id')} · {record.get('date_utc')}**")
            st.caption(f"{record.get('origin')} · release {record.get('release')} · plane {record.get('plane_id')} · calibration level {record.get('calibration_level')}")
            record_position = product.metadata.get("record_position_icrs_deg")
            requested_center = product.metadata.get("requested_field_center_icrs_deg")
            st.caption(
                f"Record position ({product.metadata.get('record_position_role', 'unspecified')}): "
                f"{record_position if record_position else 'unavailable'} ICRS deg · requested field centre: "
                f"{requested_center if requested_center else 'unavailable'} ICRS deg · FITS WCS cutout centre: "
                f"{product.center_sky if product.center_sky else 'unavailable'} ICRS deg"
            )
            st.dataframe(pd.DataFrame(product.validation), use_container_width=True, hide_index=True)
            st.caption(f"Extensions: {', '.join(product.extension_names)} · file {product.file_size_bytes:,} bytes · sha256 {product.sha256[:16]}…")
            if product.warnings:
                st.warning("\n".join(f"- {warning}" for warning in product.warnings))

if check:
    st.markdown("### Compatibility assessment")
    if check["compatible"]:
        st.success("These two observations passed every mandatory compatibility check and may be compared.")
    else:
        st.error("Direct comparison is refused: " + " | ".join(check["reasons"]))
    st.dataframe(pd.DataFrame(check["checks"]), use_container_width=True, hide_index=True)
    for warning in check["warnings"]:
        st.warning(warning)
    if check.get("time_separation_days") is not None:
        st.caption(f"Time separation between epochs: {check['time_separation_days']:.3f} days ({check['time_separation_days'] * 24:.1f} hours).")
    if check["compatible"]:
        st.page_link(navbar.page("pages/04_Cosmic_Image_Lab.py"), label="OPEN COSMIC IMAGE LAB →", use_container_width=True)

st.markdown("### Why a pair can be incompatible")
st.markdown(
    "- **Different detector / LVF band.** D1–D6 sample different wavelength ranges, so a brightness difference mixes colour with variability.\n"
    "- **WCS-WAVE mismatch inside one band.** The LVF central wavelength changes across the focal plane, so two D2 pointings at the same sky "
    "position can still be separated by tens of nanometres. The measured Δλ is compared with the narrower WCS-WAVE bandwidth.\n"
    "- **Missing or invalid WCS.** Without celestial WCS the arrays cannot be reprojected onto a common grid.\n"
    "- **Non-overlapping footprints.** A difference image requires pixels covered by both epochs after reprojection.\n"
    "- **Mixed data modes.** Real archive data, local uploads and synthetic demonstration products are never combined."
)

navbar.render_footer()
