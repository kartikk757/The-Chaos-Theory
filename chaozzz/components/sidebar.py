"""Persistent navigation and page bootstrap."""
from __future__ import annotations

import streamlit as st

import config
from components import navbar, status_panel, theme
from modules import state


def render_sidebar(active_label: str) -> None:
    from html import escape

    state.init_state()
    st.sidebar.markdown(
        f"<div class='ct-sidebar-head'>{navbar.render_brand_html(compact=True)}</div>",
        unsafe_allow_html=True,
    )
    st.sidebar.markdown("<div class='ct-nav-label'>Navigation</div>", unsafe_allow_html=True)
    for label, path, number, glyph in config.NAV_ITEMS:
        st.sidebar.page_link(
            navbar.page(path),
            label=f"{glyph}  {number} · {label}",
            use_container_width=True,
            help=f"Open {label}",
        )
    status_panel.render_status_panel()
    st.sidebar.markdown("<div style='height:.5rem'></div>", unsafe_allow_html=True)
    st.sidebar.button(
        "RESET SESSION",
        use_container_width=True,
        key="ct_reset_session",
        on_click=state.reset_session,
        help="Clear all session data and processing results",
    )
    st.sidebar.markdown(
        "<div class='ct-sidebar-foot'>NASA SPHEREx · NASA/IPAC IRSA<br/>"
        "Actual observations / processed measurements / candidates / demonstration data are kept distinct.<br/>"
        "<span style='color:#22d3ee'>Not an official NASA product.</span></div>",
        unsafe_allow_html=True,
    )


def bootstrap_page(
    title: str,
    subtitle: str,
    active_label: str,
    status_items: list[dict] | None = None,
) -> None:
    theme.load_styles(landing=False)
    state.init_state()
    state.set_value("active_page", active_label)
    render_sidebar(active_label)
    status_panel.render_data_mode_banner()
    navbar.render_top_bar(title, subtitle, status_items)
