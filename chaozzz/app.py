"""THE CHAOS THEORY — cinematic landing page."""
from __future__ import annotations

import streamlit as st

import config
from components import navbar, theme
from modules import state

st.set_page_config(
    page_title="THE CHAOS THEORY · Explore the Universe",
    page_icon="◉",
    layout="wide",
    initial_sidebar_state="collapsed",
)
theme.load_styles(landing=True)
state.init_state()
state.set_value("active_page", "Landing")

# Top navigation
st.html(
    "<div class='ct-landing-topnav'>" + navbar.render_brand_html(compact=True) +
    "<div class='ct-eyebrow'>NASA SPHEREx · IRSA PUBLIC ARCHIVE · MULTI-EPOCH SCIENCE</div></div>"
)
c1, c2, c3, c4 = st.columns([1, 1, 1, 1])
with c1:
    st.page_link(navbar.page("pages/01_Mission_Control.py"), label="MISSION CONTROL")
with c2:
    st.page_link(navbar.page("pages/02_Sky_Explorer.py"), label="SKY EXPLORER")
with c3:
    st.page_link(navbar.page("pages/07_Methodology.py"), label="METHODOLOGY")
with c4:
    st.page_link(navbar.page("pages/09_About.py"), label="ABOUT")

st.markdown("<div style='height:4.2rem'></div>", unsafe_allow_html=True)
st.markdown("<div class='ct-kicker'>Advanced observatory interface · Public science data</div>", unsafe_allow_html=True)
st.markdown(
    "<h1 class='ct-hero-title'>The Universe Is<br/><span>Never Truly Still.</span></h1>",
    unsafe_allow_html=True,
)
st.markdown(
    "<div class='ct-hero-sub'>Compare genuine NASA SPHEREx observations of the same sky across different epochs. "
    "Align spectral images on a common celestial grid, correct background differences, inspect statistically "
    "significant residuals, and investigate change without mistaking a candidate for a discovery.</div>",
    unsafe_allow_html=True,
)
st.markdown(
    "<div class='ct-hero-meta'><span class='ct-tag cyan'>Live IRSA TAP + SIA2</span>"
    "<span class='ct-tag violet'>WCS-aware difference imaging</span>"
    "<span class='ct-tag'>Provenance on every measurement</span>"
    "<span class='ct-tag warn'>No fabricated discoveries</span></div>",
    unsafe_allow_html=True,
)

cta1, cta2, cta3 = st.columns([1.25, 1.15, 1.4])
with cta1:
    if st.button("INITIATE EXPLORATION", type="primary", use_container_width=True, key="hero_initiate"):
        st.switch_page(navbar.page("pages/02_Sky_Explorer.py"))
with cta2:
    if st.button("UNDERSTAND THE MISSION", use_container_width=True, key="hero_mission"):
        st.switch_page(navbar.page("pages/07_Methodology.py"))
with cta3:
    st.page_link(navbar.page("pages/01_Mission_Control.py"), label="OPEN MISSION CONTROL", use_container_width=True)

st.markdown("<div class='ct-scroll'>Scroll for mission overview ↓</div>", unsafe_allow_html=True)

# Project overview
st.markdown("<div class='ct-section'></div>", unsafe_allow_html=True)
st.html(
    "<div class='ct-section-head'><div><div class='ct-kicker'>01 / Project overview</div>"
    "<h2 class='ct-section-title'>An instrument for investigating change.</h2></div>"
    "<div class='ct-section-note'>The scientific objective is to investigate potential changes in observed sky brightness and position — "
    "not to claim the discovery of Planet X. SPHEREx is an all-sky spectral survey mission, not a dedicated planet or asteroid tracker.</div></div>"
)
st.html(
    "<div class='ct-grid'>"
    "<div class='ct-card'><div class='ct-card-index'>01 · SELECT</div><h3>Choose a verified sky region</h3>"
    "<p>Enter ICRS coordinates, a small cutout radius, an observation date range and a SPHEREx detector band. "
    "The interface distinguishes a requested position from actual archive coverage.</p></div>"
    "<div class='ct-card'><div class='ct-card-index'>02 · RETRIEVE</div><h3>Query the real archive</h3>"
    "<p>Observation metadata comes from IRSA TAP; SIA2 supplies footprint coverage; documented cutout URLs return compact "
    "multi-extension FITS files rather than 70 MB full spectral images.</p></div>"
    "<div class='ct-card'><div class='ct-card-index'>03 · COMPARE</div><h3>Align, normalize, subtract</h3>"
    "<p>WCS reprojection, zodiacal-model subtraction, robust background correction and a transparent normalization step "
    "produce D(x,y) = I₂(x,y) − I₁(x,y) and an approximate signal-to-noise map.</p></div>"
    "<div class='ct-card'><div class='ct-card-index'>04 · INVESTIGATE</div><h3>Read the evidence</h3>"
    "<p>Residual sources are measured, cross-matched where possible, classified conservatively, and accompanied by "
    "flags, warnings, quality indicators, processing configuration and exact input provenance.</p></div>"
    "</div>"
)

# Mission information
st.html(
    "<div class='ct-section-head' style='margin-top:3.5rem'><div><div class='ct-kicker'>02 / Mission</div>"
    "<h2 class='ct-section-title'>SPHEREx: an all-sky infrared spectral survey.</h2></div>"
    "<div class='ct-section-note'>Launched in March 2025, SPHEREx measures 0.75–5.0 µm spectra across the sky with deeper "
    "coverage near the ecliptic poles. Public DR1 and QR3 spectral images are served by NASA/IPAC IRSA.</div></div>"
)
st.html(
    "<div class='ct-grid'>"
    "<div class='ct-card'><div class='ct-card-index'>SPECTRAL RANGE</div><h3>0.75 – 5.0 µm</h3>"
    "<p>Six detector/LVF bands with wavelength varying across each detector. A fixed sky position can therefore be observed "
    "at a different central wavelength in another epoch; the app checks WCS-WAVE before comparison.</p></div>"
    "<div class='ct-card'><div class='ct-card-index'>LEVEL 2 PRODUCT</div><h3>Multi-extension FITS</h3>"
    "<p>IMAGE, FLAGS, VARIANCE, ZODI, ePSF and WCS-WAVE extensions. IMAGE is calibrated surface brightness in MJy/sr and "
    "has not had zodiacal light subtracted.</p></div>"
    "<div class='ct-card'><div class='ct-card-index'>PUBLIC ARCHIVE</div><h3>IRSA DR1 + QR3</h3>"
    "<p>SIA2 collections and CAOM TAP tables expose searchable metadata. The IRSA cutout service returns small spatial "
    "subsets suitable for limited hardware.</p></div>"
    "</div>"
)

# Methodology preview
st.html(
    "<div class='ct-section-head' style='margin-top:3.5rem'><div><div class='ct-kicker'>03 / Scientific pipeline</div>"
    "<h2 class='ct-section-title'>Every result has a traceable path.</h2></div>"
    "<div class='ct-section-note'>Image artifacts, imperfect registration, zodiacal background and LVF wavelength changes "
    "can all imitate astrophysical variability. The methodology panel explains each safeguard.</div></div>"
)
steps = [
    "RAW SPHEREx DATA", "OBSERVATION VALIDATION", "WCS IMAGE ALIGNMENT",
    "BACKGROUND CORRECTION", "IMAGE COMPARISON", "SOURCE DETECTION",
    "CATALOGUE CROSS-MATCH", "CANDIDATE INVESTIGATION",
]
pipeline_html = "<div class='ct-pipeline'>" + "<span class='ct-step'>" + "</span><span class='ct-step-arrow'>↓</span><span class='ct-step'>".join(steps) + "</span></div>"
st.html(pipeline_html)
m1, m2 = st.columns([1, 1])
with m1:
    st.page_link(navbar.page("pages/07_Methodology.py"), label="READ THE FULL METHODOLOGY", use_container_width=True)
with m2:
    st.page_link(navbar.page("pages/01_Mission_Control.py"), label="VIEW LIVE APPLICATION STATE", use_container_width=True)

navbar.render_footer()
