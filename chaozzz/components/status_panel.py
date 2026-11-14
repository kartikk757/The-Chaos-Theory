"""Data-mode banners, sidebar status, honest error reporting and empty states."""
from __future__ import annotations

import traceback
from html import escape

import streamlit as st

import config
from modules import state


def render_data_mode_banner() -> None:
    mode_key = state.get("mode_key", "EMPTY")
    label = config.DATA_MODES.get(mode_key, "UNKNOWN DATA MODE")
    note = state.get("mode_note", "")
    if mode_key == "EMPTY":
        st.html(
            "<div class='ct-banner'><span class='ct-tag'>Awaiting data</span> "
            "No observations or processed measurements are loaded. Values shown as "
            "<strong>Awaiting data</strong> are empty states, not scientific results.</div>"
        )
        return
    css = {"DEMO": "demo", "CACHED_REAL": "real", "LIVE": "real", "LOCAL": "local"}.get(mode_key, "")
    if mode_key == "DEMO":
        message = (
            "<strong>DEMONSTRATION DATA — SYNTHETIC.</strong> Every image, date, measurement and candidate in this "
            "session is generated for software validation. It is not a SPHEREx observation and must not be "
            "interpreted as astronomical evidence."
        )
    elif mode_key == "CACHED_REAL":
        message = (
            "<strong>REAL SPHEREx DATA — CACHED IRSA CUTOUTS.</strong> These are genuine, provenance-tracked "
            "SPHEREx FITS cutouts downloaded from NASA/IPAC IRSA. The cached files are the exact cutouts "
            "returned by the archive, including IMAGE, FLAGS, VARIANCE, ZODI, EPSF and WCS-WAVE."
        )
    elif mode_key == "LIVE":
        message = "<strong>LIVE IRSA SPHEREx.</strong> Observation metadata and FITS cutouts are retrieved from the official NASA/IPAC IRSA services."
    elif mode_key == "LOCAL":
        message = (
            "<strong>LOCAL FITS UPLOAD — UNVERIFIED PROVENANCE.</strong> The application can process the supplied "
            "arrays and WCS, but it cannot verify that the file is a SPHEREx observation."
        )
    else:
        message = "<strong>UNKNOWN DATA MODE.</strong> Provenance could not be established."
    suffix = f" <span class='ct-quiet'>{escape(note)}</span>" if note else ""
    st.html(
        f"<div class='ct-banner {css}'><span class='ct-tag {'danger' if mode_key == 'DEMO' else 'cyan'}'>"
        f"{escape(label)}</span> {message}{suffix}</div>"
    )


def render_status_panel() -> None:
    summary = state.state_summary()
    archive = summary["archive_status"] or {}
    st.sidebar.markdown("<div class='ct-nav-label'>Instrument state</div>", unsafe_allow_html=True)
    st.sidebar.markdown(
        "<div class='ct-glass' style='padding:.8rem;background:rgba(255,255,255,.018)'>"
        f"<div class='ct-quiet'><strong>Mode</strong><br/>{escape(summary['data_mode'])}</div>"
        f"<div class='ct-quiet' style='margin-top:.55rem'><strong>IRSA</strong><br/>{escape(str(archive.get('status','Not checked')))}</div>"
        f"<div class='ct-quiet' style='margin-top:.55rem'><strong>Observations</strong> {summary['observation_count']} · "
        f"<strong>Loaded FITS</strong> {summary['product_count']} · <strong>Candidates</strong> {summary['candidate_count']}</div>"
        f"<div class='ct-quiet' style='margin-top:.55rem'><strong>Session</strong><br/>{escape(str(summary['session_id']))}</div>"
        "</div>",
        unsafe_allow_html=True,
    )


def render_exception(error: Exception, context: str = "Operation") -> None:
    details = traceback.format_exc()
    state.set_value("last_error", f"{context}: {error}")
    state.add_log(f"{context} failed", level="error", details={"error": repr(error), "traceback": details[-4000:]})
    st.error(f"{context} failed: {error}")
    with st.expander("Technical details — full error is not hidden"):
        st.code(details or repr(error), language="python")


def render_empty_state(title: str, message: str, page: str, button_label: str) -> None:
    st.html(
        "<div class='ct-empty'><div class='ct-kicker'>Empty state</div>"
        f"<div class='ct-empty-title'>{escape(title)}</div>"
        f"<p class='ct-quiet'>{escape(message)}</p></div>"
    )
    st.page_link(st.Page(str(config.ROOT / page)), label=button_label, use_container_width=True)


def render_loading(label: str) -> None:
    st.html(
        "<div class='ct-glass' style='text-align:center;padding:1.4rem'>"
        f"<div class='ct-loading'></div><div class='ct-quiet'>{escape(label)}</div></div>"
    )
