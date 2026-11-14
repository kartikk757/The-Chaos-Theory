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
    ("MULTI-EPOCH SCREENING", "A separate path checks every selected epoch pair, detects peaks, associates sky positions, fits local linear motion when at least three detections support it, and saves its own provenance."),
    ("OPTIONAL EXTERNAL EVIDENCE", "Gaia DR3 and IRSA catalogues add counterpart context; a user-supplied known-object query can test a JPL Horizons hypothesis. None of these calls makes a discovery."),
    ("SPECTRAL QUICK-LOOK", "A validated mosaic or HiPS FITS cube can provide an aperture-minus-annulus spectrum. It is kept separate from time-domain motion measurements."),
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
    "an artifact usually does not. With two epochs the pair pipeline can only report displacement between two measurements. The separate "
    "multi-epoch screen evaluates every selected pair for data mode, spectral compatibility and spatial overlap, then associates sky positions. "
    "A motion classification requires at least three detections. Sparse sampling, peak positions and an uncalibrated uncertainty floor mean "
    "that the fitted rate is a review lead, not an orbit or a confirmed moving object."
)

st.markdown("### Catalogue and known-object context")
st.markdown(
    "Gaia DR3 positions are propagated from the catalogue reference epoch using the reported linear proper motion, then compared at each "
    "candidate observation epoch with an uncertainty-scaled radius. The current model does not apply parallax or radial-velocity perspective "
    "effects. 2MASS PSC and AllWISE matches retain multiple counterparts and quality information; their positions are not propagated to each "
    "SPHEREx epoch. JPL SBDB/Horizons is queried only for a known identifier supplied by the user. A positional match tests that hypothesis "
    "and does not establish object identity; no match remains inconclusive."
)

st.markdown("### ePSF and spectral products")
st.markdown(
    "The official ePSF extension can be inspected by mapping an image position through its WCS to the parent detector zone. The model is "
    "currently diagnostic: candidate fluxes are not ePSF fits, and the app does not perform PSF-matched differencing. A mosaic or HiPS cube "
    "can be used for a quick-look aperture spectrum with annulus background subtraction and channel-quality diagnostics. This spectral "
    "measurement is spatial/spectral context and is not an additional time-separated epoch."
)

st.markdown("### Controlled injection and recovery")
st.markdown(
    "The optional validation panel injects an uploaded detector-sampled kernel into a user-selected control image and reports recovery "
    "fractions, unassociated peaks, a seed and input checksums. These results apply only to those settings. They are not a calibrated "
    "completeness curve, false-positive probability or sensitivity improvement."
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
        {"quantity": "Multi-epoch angular rate", "units": "arcsec/day", "meaning": "Weighted local tangent-plane linear fit; requires at least three detections and uses a stated, currently uncalibrated positional uncertainty floor."},
        {"quantity": "Gaia match radius", "units": "arcsec", "meaning": "Bounded radius based on Gaia propagated and candidate positional uncertainties at each measured epoch."},
        {"quantity": "JPL angular separation", "units": "arcsec", "meaning": "Separation between a candidate coordinate and the Horizons ephemeris for a user-supplied known object at the recorded UTC-derived time."},
        {"quantity": "Cube quick-look spectrum", "units": "FITS BUNIT by wavelength channel", "meaning": "Aperture flux minus annulus background at the cube WCS position; separate from time-domain evidence."},
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
