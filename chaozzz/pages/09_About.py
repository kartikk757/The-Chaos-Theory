"""About the project, the mission and the data sources."""
from __future__ import annotations

import streamlit as st

import config
from components import navbar, sidebar
from modules import provenance

sidebar.bootstrap_page(
    "About",
    "The project, NASA SPHEREx, the archive services used, and the boundaries of this prototype.",
    "About",
)

st.markdown("### THE CHAOS THEORY")
st.markdown(
    "**Explore the Universe. Detect the Unexpected.**\n\n"
    "THE CHAOS THEORY is a hackathon prototype for multi-epoch exploration of public NASA SPHEREx observations. It connects to the "
    "NASA/IPAC Infrared Science Archive (IRSA), retrieves real spectral-image cutouts, aligns them on a common celestial grid, and measures "
    "what changed between two visits to the same patch of sky. Every number displayed in the application derives from FITS data that was "
    "actually loaded; nothing is pre-computed or invented."
)

st.markdown("### The mission")
mission_columns = st.columns(2)
with mission_columns[0]:
    st.markdown(
        "**NASA SPHEREx** (Spectro-Photometer for the History of the Universe, Epoch of Reionization, and Ices Explorer) is a NASA "
        "astrophysics mission that maps the entire sky in near-infrared light. Its low-resolution spectrometers sample 0.75–5.00 µm in 102 "
        "wavelength channels across six detector bands (D1–D6), repeatedly revisiting the sky so that the same region is observed at many "
        "different epochs.\n\n"
        "SPHEREx is not a dedicated asteroid or planet-hunting telescope, and it was not designed for Planet X searches. Its repeated "
        "wide-field coverage is what makes multi-epoch change investigation possible in this prototype."
    )
with mission_columns[1]:
    st.markdown(
        "**Data used here**\n\n"
        "- Source: NASA/IPAC IRSA SPHEREx archive (`spherex.plane` and `spherex.artifact` via TAP; SIA2 collections for coverage).\n"
        "- Products: Level-2 spectral images delivered as cutouts through the documented cutout service.\n"
        "- Extensions read: `IMAGE`, `FLAGS`, `VARIANCE`, `ZODI`, `EPSF`, `WCS-WAVE`.\n"
        "- Units: surface brightness in MJy/sr for the image planes; fluxes converted to Jy using the pixel solid angle.\n"
        "- Attribution: NASA SPHEREx · NASA/IPAC Infrared Science Archive, Caltech/IPAC."
    )

st.markdown("### Verification trail")
st.markdown(
    "Before any retrieval code was written, the following official pages were read and the live services were probed: the SPHEREx Data "
    "Access guide, the data-product description (MEF extensions, FLAGS bitmap, LVF bands), the IRSA cutout tutorial (TAP join and cutout "
    "URL construction), the SPHEREx introduction tutorial (spectral WCS usage), and the IRSA TAP documentation. The archive endpoints in "
    "`config.py` are those documented values, and the queries executed by the application can be inspected in Data Explorer."
)
for doc in config.ARCHIVE_DOCS:
    st.markdown(f"- [{doc['title']}]({doc['url']})")

st.markdown("### What this prototype does not claim")
st.markdown(
    "- No discovery of any kind. Candidate records are screening labels produced from two observations.\n"
    "- No Planet X claim, and no claim that SPHEREx data can identify planets.\n"
    "- No use of demonstration data without an unmistakable label. Synthetic products are marked "
    "**DEMONSTRATION DATA — SYNTHETIC** in their FITS headers, in the manifest, and in the interface banner.\n"
    "- No fallback that silently substitutes one observation for another: if an archive request fails, the error is shown and alternatives "
    "(another observation, a local FITS upload, the cached real pair, or labelled demonstration data) are offered explicitly."
)

st.markdown("### Reproducibility")
versions = provenance.software_versions()
st.dataframe([{"component": key, "version": value} for key, value in versions.items()], use_container_width=True, hide_index=True)
st.caption(
    "Run `python3 scripts/build_demo_bundle.py --check` to verify the cached real cutouts, and `python3 -m pytest tests -q` to execute the "
    "unit and integration suite. Both commands are documented in README.md."
)

st.markdown("### Contact and licence")
st.markdown(
    "This is an independent hackathon project. SPHEREx data are public and are served by IRSA; please cite the mission and the archive in "
    "any downstream use. The application code is provided for demonstration and educational purposes, and is **not** an official NASA or "
    "Caltech/IPAC product."
)

navbar.render_footer()
