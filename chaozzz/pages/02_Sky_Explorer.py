"""Sky region selection, IRSA discovery and coverage inspection."""
from __future__ import annotations

import pandas as pd
import streamlit as st

import config
from components import metric_cards, navbar, sidebar, status_panel
from modules import data_loader, irsa_client, state, visualization, wcs_alignment

sidebar.bootstrap_page(
    "Sky Explorer",
    "Define a sky region, query the NASA/IPAC IRSA SPHEREx archive and inspect real coverage.",
    "Explore Sky",
)

query = dict(state.get("query", config.DEFAULT_QUERY))

PRESETS = {
    "M101 field · verified SPHEREx region (RA 210.80227°, Dec +54.34895°)": {
        "ra_deg": 210.80227, "dec_deg": 54.34895, "radius_deg": 0.05, "size_deg": 0.10,
        "start_date": "2026-05-01", "end_date": "2026-06-30", "band": "SPHEREx-D2",
    },
    "M101 field · wide 2025–2026 date range": {
        "ra_deg": 210.80227, "dec_deg": 54.34895, "radius_deg": 0.06, "size_deg": 0.12,
        "start_date": "2025-05-01", "end_date": "2026-06-30", "band": "SPHEREx-D2",
    },
    "Ecliptic-region test field (RA 95.00000°, Dec +12.00000°)": {
        "ra_deg": 95.0, "dec_deg": 12.0, "radius_deg": 0.05, "size_deg": 0.10,
        "start_date": "2025-05-01", "end_date": "2026-06-30", "band": "SPHEREx-D1",
    },
}

st.markdown("### 01 · Target field")
preset = st.selectbox("Field preset", ["Custom / current session values", *PRESETS.keys()], index=0)
if preset != "Custom / current session values":
    selected = PRESETS[preset]
    # A preset must overwrite widget state, otherwise the widgets keep the values
    # they were first rendered with.
    st.session_state["sky_ra"] = str(selected["ra_deg"])
    st.session_state["sky_dec"] = str(selected["dec_deg"])
    st.session_state["sky_radius"] = float(selected["radius_deg"])
    st.session_state["sky_size"] = float(selected["size_deg"])
    st.session_state["sky_start"] = pd.to_datetime(selected["start_date"]).date()
    st.session_state["sky_end"] = pd.to_datetime(selected["end_date"]).date()
    st.session_state["sky_band"] = selected["band"]

with st.form("ct_sky_query", clear_on_submit=False):
    columns = st.columns([1, 1, 1, 1])
    with columns[0]:
        ra_text = st.text_input("Right Ascension (deg or hh:mm:ss)", value=str(query.get("ra_deg", "")), key="sky_ra", help="ICRS. Decimal degrees or sexagesimal hours, e.g. 14:03:12.6")
    with columns[1]:
        dec_text = st.text_input("Declination (deg or dd:mm:ss)", value=str(query.get("dec_deg", "")), key="sky_dec", help="ICRS. Decimal degrees or sexagesimal degrees, e.g. +54:20:56")
    with columns[2]:
        radius = st.number_input("Search radius (deg)", min_value=0.005, max_value=0.25, value=float(query.get("radius_deg", 0.05)), step=0.005, format="%.3f", key="sky_radius", help="Cone radius used for archive footprint searches. Max 0.25° keeps the query and cutouts small.")
    with columns[3]:
        size = st.number_input("Cutout width (deg)", min_value=0.01, max_value=0.5, value=float(query.get("size_deg", 0.10)), step=0.01, format="%.3f", key="sky_size", help="Width of the delivered FITS cutout, per the documented IRSA cutout service.")

    columns = st.columns([1, 1, 1, 1])
    with columns[0]:
        start_date = st.date_input("Start date (UTC)", value=pd.to_datetime(query.get("start_date", "2026-05-01")).date(), key="sky_start")
    with columns[1]:
        end_date = st.date_input("End date (UTC)", value=pd.to_datetime(query.get("end_date", "2026-06-30")).date(), key="sky_end")
    with columns[2]:
        band_options = ["ALL", *config.DETECTOR_BANDS.keys()]
        band = st.selectbox("SPHEREx detector band", band_options, index=band_options.index(query.get("band", "SPHEREx-D2")) if query.get("band", "SPHEREx-D2") in band_options else 1, key="sky_band", format_func=lambda value: "All bands" if value == "ALL" else f"{value} · {config.DETECTOR_BANDS[value]['range_um'][0]}–{config.DETECTOR_BANDS[value]['range_um'][1]} µm")
    with columns[3]:
        release = st.selectbox("Archive release", ["all", "dr1", "qr3"], index=["all", "dr1", "qr3"].index(query.get("release", "all")), key="sky_release", help="DR1 and QR3 are the releases searchable through SIA2; QR2 is superseded and not SIA2-searchable.")
    columns = st.columns([1, 1, 1])
    with columns[0]:
        survey = st.selectbox("Survey", ["wide", "deep", "all"], index=["wide", "deep", "all"].index(query.get("survey", "wide")), key="sky_survey")
    with columns[1]:
        max_rows = st.number_input("Maximum observations", min_value=1, max_value=300, value=int(query.get("max_rows", 120)), step=10, key="sky_max_rows")
    with columns[2]:
        st.markdown("<div style='height:1.55rem'></div>", unsafe_allow_html=True)
        submitted = st.form_submit_button("QUERY IRSA ARCHIVE", use_container_width=True, type="primary")

st.caption(
    "The requested position is the centre of the cone search and of each cutout. Archive coverage describes where SPHEREx actually "
    "pointed; a returned observation's footprint centre can differ from the requested position because SPHEREx images are long strips."
)

actions = st.columns(3)
with actions[0]:
    check = st.button("CHECK ARCHIVE CONNECTION", use_container_width=True, help="Query a single archive row from the TAP service and an SIA2 collection.")
with actions[1]:
    load_cached = st.button("LOAD VERIFIED CACHED PAIR", use_container_width=True, help="Load the two real SPHEREx cutouts already downloaded from IRSA for the M101 field.")
with actions[2]:
    load_demo = st.button("LOAD DEMONSTRATION DATA", use_container_width=True, help="Synthetic two-epoch products for offline software demonstration. Never a SPHEREx observation.")

if check:
    with st.spinner("Contacting NASA/IPAC IRSA…"):
        try:
            status = irsa_client.check_archive_status()
            state.set_value("archive_status", status)
            state.add_log("Archive status checked", level="info", details=status)
            st.success(f"IRSA: {status['status']} — {status['detail']} ({status['latency_seconds']} s)")
        except Exception as error:
            status_panel.render_exception(error, "Archive status check")
    with st.expander("VO client verification (PyVO / Astroquery)", expanded=False):
        if st.button("RUN VO CLIENT CHECK", key="ct_vo_check"):
            with st.spinner("Querying IRSA TAP through PyVO and Astroquery…"):
                result = irsa_client.verify_vo_clients()
            st.json(result)

if submitted:
    try:
        validated = irsa_client.validate_query(
            ra_deg=ra_text, dec_deg=dec_text, radius_deg=radius, size_deg=size,
            start_date=str(start_date), end_date=str(end_date), band=band,
            release=release, survey=survey, max_rows=max_rows,
        )
    except Exception as error:
        status_panel.render_exception(error, "Coordinate / query validation")
        validated = None
    if validated:
        with st.spinner("Querying spherex.plane + spherex.artifact through IRSA TAP and SIA2 footprints…"):
            try:
                result = irsa_client.query_archive(validated)
                state.set_query(**validated)
                state.set_observations(result["observations"], "LIVE", note=f"{len(result['observations'])} real archive rows")
                state.set_coverage(result["coverage_rows"], result["coverage_status"])
                state.set_value("archive_status", result["archive_status"])
                state.set_value("last_tap_diagnostics", result["tap_diagnostics"])
                st.success(f"IRSA returned {len(result['observations'])} real SPHEREx observation records.")
            except Exception as error:
                status_panel.render_exception(error, "IRSA archive query")
                st.info("If the archive is unreachable or this field has no coverage, use LOAD VERIFIED CACHED PAIR or LOAD DEMONSTRATION DATA above; both are labelled explicitly.")

if load_cached:
    try:
        with st.spinner("Loading cached real SPHEREx cutouts…"):
            records, manifest = data_loader.load_cached_real_bundle()
        state.set_observations(records, "CACHED_REAL", note=f"manifest {manifest['generated_utc']}")
        state.set_coverage([], {"status": "cached files", "detail": "Coverage footprints are not re-queried in cached mode; re-run a live IRSA query for SIA2 footprints."})
        state.set_query(ra_deg=config.DEFAULT_QUERY["ra_deg"], dec_deg=config.DEFAULT_QUERY["dec_deg"], radius_deg=config.DEFAULT_QUERY["radius_deg"], size_deg=config.DEFAULT_QUERY["size_deg"], band=config.DEFAULT_QUERY["band"])
        st.success(f"Loaded {len(records)} verified real SPHEREx cutouts; checksums are recorded in data/real/manifest.json (generated {manifest['generated_utc']}).")
    except Exception as error:
        status_panel.render_exception(error, "Cached real bundle load")

if load_demo:
    try:
        from modules import provenance

        manifest = provenance.read_json(config.DEMO_DATA_DIR / "manifest.json")
        state.set_observations(manifest["observations"], "DEMO", note=manifest["notice"])
        state.set_coverage([], {"status": "synthetic", "detail": "Synthetic products carry no archive footprint."})
        state.set_query(ra_deg=210.80227, dec_deg=54.34895, radius_deg=0.05, size_deg=0.10, band="SPHEREx-D2 (synthetic stand-in)")
        st.warning("DEMONSTRATION DATA loaded. These files are synthetic and are not NASA SPHEREx observations.")
    except Exception as error:
        status_panel.render_exception(error, "Demonstration data load")

with st.expander("LOCAL FITS UPLOAD — fallback when the archive cannot supply a comparable pair", expanded=False):
    st.caption(
        "Upload one or two FITS files (for example a SPHEREx cutout you downloaded yourself). Local files are processed normally, but their "
        "provenance cannot be verified, so they are labelled LOCAL FITS UPLOAD · UNVERIFIED PROVENANCE and can only be compared with other "
        "local files. Uploading synthetic or non-SPHEREx data is permitted for software testing and is shown as such."
    )
    uploaded = st.file_uploader(
        "FITS files (max two, 200 MB each)",
        type=["fits", "fit", "fts", "gz"],
        accept_multiple_files=True,
        key="local_upload_widget",
        help="Standard FITS with a 2-D IMAGE extension and a celestial WCS. Extra FLAGS/VARIANCE/ZODI extensions are used when present.",
    )
    if st.button("PROCESS LOCAL FITS", disabled=not uploaded, key="local_upload_process"):
        prepared_records = []
        prepared_products = {}
        try:
            for index, handle in enumerate(list(uploaded)[:2], start=1):
                if len(uploaded) > 2:
                    st.warning("Only the first two uploaded files are used.")
                upload = data_loader.save_uploaded_fits(handle)
                placeholder = {
                    "id": f"LOCAL-UPLOAD-{index}",
                    "origin": "User-supplied local FITS upload (provenance unverified)",
                    "data_mode": config.DATA_MODES["LOCAL"],
                }
                product = data_loader.load_fits_product(upload["path"], placeholder)
                record = data_loader.make_local_observation(upload, product, index)
                prepared_records.append(record)
                prepared_products[record["id"]] = product
            state.set_observations(prepared_records, "LOCAL", note=f"{len(prepared_records)} user-supplied FITS file(s)")
            state.set_value("fits_products", prepared_products)
            state.set_coverage([], {"status": "local upload", "detail": "Local FITS files carry no archive footprint."})
            if len(prepared_records) == 2:
                ids = [record["id"] for record in prepared_records]
                state.set_selected_pair(ids)
                state.set_value(
                    "pair_check",
                    wcs_alignment.check_observation_compatibility(
                        prepared_records[0], prepared_records[1], prepared_products[ids[0]], prepared_products[ids[1]],
                        state.get("processing_config", config.PROCESSING_DEFAULTS),
                    ),
                )
                st.success("Two local FITS files loaded and checked. Review the compatibility result in Observation Explorer.")
            else:
                st.success("Local FITS file loaded. Upload a second file of the same region to compare epochs.")
        except Exception as error:
            status_panel.render_exception(error, "Local FITS upload")

current = state.get("query", {})
observations = state.get("observations", [])
coverage_rows = state.get("coverage_rows", [])
coverage_status = state.get("coverage_status", {})

st.markdown("### 02 · Requested field and archive coverage")
metric_cards.render_metric_grid([
    {"label": "Requested RA / Dec", "value": f"{float(current.get('ra_deg', 0)):.5f}° / {float(current.get('dec_deg', 0)):.5f}°", "caption": "ICRS, verified through coordinate parsing", "icon": "⌖"},
    {"label": "Search radius", "value": f"{float(current.get('radius_deg', 0)):.3f}°", "caption": "Cone used for footprint queries", "icon": "◎"},
    {"label": "Cutout width", "value": f"{float(current.get('size_deg', 0)):.3f}°", "caption": "Documented IRSA cutout size parameter", "icon": "▣"},
    {"label": "Archive rows", "value": str(len(observations)) if observations else "Awaiting data", "caption": coverage_status.get("detail", "No archive query in this session")[:90], "icon": "≡", "accent": bool(observations)},
])

if observations:
    figure = visualization.sky_footprint_figure(
        float(current.get("ra_deg", 210.8)), float(current.get("dec_deg", 54.3)),
        float(current.get("radius_deg", 0.05)), coverage_rows, observations,
    )
    st.plotly_chart(figure, use_container_width=True, theme=None, config={"displaylogo": False})
    st.caption(
        "Cyan polygons are real SIA2 footprints returned by IRSA; the dotted violet circle is the requested search radius; dots mark "
        "observation centres. A requested position is not itself an observation — coverage is what the archive actually contains."
    )
    st.plotly_chart(visualization.timeline_figure(observations), use_container_width=True, theme=None, config={"displaylogo": False})
    frame = pd.DataFrame([
        {
            "date (UTC)": obs.get("date_utc"),
            "band": obs.get("band"),
            "release": obs.get("release"),
            "observation": obs.get("obs_id"),
            "RA (deg)": obs.get("ra_deg"),
            "Dec (deg)": obs.get("dec_deg"),
            "coordinates from": obs.get("coordinates_source"),
            "cutout size": f"{obs.get('content_length', 0):,} B" if obs.get("content_length") else None,
        }
        for obs in observations
    ])
    st.dataframe(frame, use_container_width=True, hide_index=True)
    st.page_link(navbar.page("pages/03_Observation_Explorer.py"), label="CONTINUE TO OBSERVATION EXPLORER →", use_container_width=True)
else:
    status_panel.render_empty_state(
        "No observations in this session",
        "Run a live IRSA query, load the verified cached pair, or load demonstration data. Nothing is filled in on your behalf.",
        "pages/03_Observation_Explorer.py",
        "OPEN OBSERVATION EXPLORER",
    )

with st.expander("Query provenance and archive diagnostics", expanded=False):
    diagnostics = state.get("last_tap_diagnostics") or {}
    if diagnostics:
        st.code(diagnostics.get("query", ""), language="sql")
        st.json({key: value for key, value in diagnostics.items() if key != "query"})
    else:
        st.caption("No live TAP query has been executed in this session.")
    st.json(coverage_status or {"status": "not queried"})
    st.markdown("Verified documentation: " + " · ".join(f"[{doc['title']}]({doc['url']})" for doc in config.ARCHIVE_DOCS[:4]))

navbar.render_footer()
