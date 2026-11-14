"""Methodology: how the science pipeline works and where its limits are."""
from __future__ import annotations

import streamlit as st

import config
from components import navbar, sidebar
from modules import state

sidebar.bootstrap_page(
    "Methodology",
    "What the pipeline does, why each step is necessary, and where its conclusions stop.",
    "Methodology",
)

STAGES = [
    ("RAW SPHEREx DATA", "Level-2 spectral image products retrieved from NASA/IPAC IRSA: IMAGE (MJy/sr), FLAGS, VARIANCE, ZODI, EPSF and the WCS-WAVE lookup table."),
    ("OBSERVATION VALIDATION", "FITS header inspection, 2-D array and unit checks, celestial WCS construction, extension shape consistency and WCS-WAVE evaluation at the cutout centre."),
    ("WCS IMAGE ALIGNMENT", "Both epochs are reprojected onto a common celestial grid (the epoch-1 WCS) with reproject. Footprint overlap is measured; if less than the configured minimum, the comparison is refused."),
    ("BACKGROUND CORRECTION", "The documented ZODI model is subtracted, then a robust tilted sky plane and a sigma-clipped constant are removed per epoch, excluding masked pixels."),
    ("IMAGE COMPARISON", "D(x, y) = I₂(x, y) − I₁(x, y), with a normalized SNR representation using the propagated VARIANCE planes where available."),
    ("SOURCE DETECTION", "Photutils peak detection on the SNR map (box-car, appropriate for an undersampled PSF), with ± dipoles merged so a shifted object is one record."),
    ("CATALOGUE CROSS-MATCH", "Cone search around each candidate in an IRSA-hosted survey table (2MASS PSC by default) to test whether the residual is a known object."),
    ("CANDIDATE INVESTIGATION", "Per-candidate measurements, quality flags, warnings and a plain-language explanation of why the record was flagged."),
]

st.markdown("### The pipeline")
for index, (title, description) in enumerate(STAGES):
    st.html(
        "<div class='ct-glass' style='padding:.85rem 1rem;margin-bottom:.5rem'>"
        f"<div class='ct-mono-note' style='letter-spacing:.22em;color:#46f0ff'>{index + 1:02d} · {title}</div>"
        f"<div class='ct-quiet' style='margin-top:.35rem'>{description}</div>"
        "</div>"
    )

st.markdown("### Why images must be aligned")
st.markdown(
    "SPHEREx pointings are not identical. The telescope scans, so the same sky position lands on different detector pixels in different "
    "visits, and each product carries its own tangent-plane projection. Subtracting two unaligned arrays produces a lattice of false "
    "positives wherever a gradient crosses a pixel edge. The pipeline therefore reprojects both epochs onto a single grid and refuses to "
    "proceed if the overlap is too small to be meaningful."
)

st.markdown("### Why backgrounds differ")
st.markdown(
    "The zodiacal background moves through the year, the detector response is not perfectly flat, and stray light changes between "
    "visits. SPHEREx Level-2 products include a ZODI model but no zodiacal subtraction in IMAGE; subtracting that model, plus a robust "
    "tilted plane, removes most of the additive difference. Residual multiplicative differences are handled with a robust linear sky "
    "scaling coefficient that is only applied when it is physically plausible."
)

st.markdown("### Why artifacts look like astronomical changes")
st.markdown(
    "Cosmic-ray hits, persistence, ghosts, streaks, snowballs and warm pixels can all appear in one epoch only, which is exactly the "
    "signature of a transient. Because of this, the pipeline reads the FLAGS bit names from each FITS header and masks severe bits, keeps "
    "informational bits such as SOURCE and FULLSAMPLE, records the flag fraction around every candidate and can label a record as a likely "
    "image artifact. It also reports single-pixel spikes as ambiguous rather than astrophysical."
)

st.markdown("### Why more than two epochs matter")
st.markdown(
    "Two epochs give one difference. A moving object should continue along a track and a variable object should repeat its behaviour; "
    "an artifact usually does not. With two epochs the pipeline can only report displacement between the two measurements and flag whether "
    "it exceeds the alignment tolerance. Multi-epoch tracking is the natural next step and is described in the README roadmap."
)

st.markdown("### Why a candidate is not a discovery")
st.markdown(
    "A candidate record means: a pixel group changed by more than the measured noise between two specific observations, with the stated "
    "quality caveats. Establishing an astrophysical object requires independent evidence — multiple epochs, consistent astrometry, "
    "cross-matching against known object catalogues and physical modelling. THE CHAOS THEORY never assigns planet, asteroid or discovery "
    "labels, and this prototype has no Planet X claim of any kind."
)

st.markdown("### Numbers you will see, and what they mean")
st.dataframe(
    [
        {"quantity": "D(x, y)", "units": "FITS IMAGE units (normally MJy/sr)", "meaning": "Difference of the two normalized epochs at one pixel."},
        {"quantity": "SNR", "units": "σ", "meaning": "Residual divided by sqrt(VARIANCE₁ + VARIANCE₂); a global robust estimate is used only when variance planes are missing."},
        {"quantity": "Aperture flux", "units": "Jy", "meaning": "Annulus-subtracted aperture sum converted with the pixel solid angle."},
        {"quantity": "Angular displacement", "units": "arcsec", "meaning": "Spherical separation between matched compact counterparts in the two epochs."},
        {"quantity": "Position angle", "units": "degrees east of north", "meaning": "Direction of that displacement on the sky."},
        {"quantity": "Δt", "units": "days", "meaning": "Separation between the archived observation times."},
        {"quantity": "Alignment tolerance", "units": "arcsec", "meaning": "1.5 × reprojected pixel scale; displacements below it are not treated as motion."},
    ],
    use_container_width=True, hide_index=True,
)

st.markdown("### Configuration actually used in this session")
st.json(state.get("processing_config", config.PROCESSING_DEFAULTS))

st.markdown("### Verified documentation")
for doc in config.ARCHIVE_DOCS:
    st.markdown(f"- [{doc['title']}]({doc['url']})")

st.caption(
    "Every statement on this page describes behaviour implemented in the modules/ layer. Where a step could not be implemented (for example "
    "PSF-matched differencing), it is described here as a limitation rather than presented as done."
)

navbar.render_footer()
