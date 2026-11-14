"""Instrument-style top status bar and shared footer."""
from __future__ import annotations

from html import escape
from typing import Any

import streamlit as st

import config
from modules import state


def page(relative_path: str):
    """Return a Streamlit Page object that works in both MPA and AppTest contexts."""
    return st.Page(str(config.ROOT / relative_path))


def render_brand_html(compact: bool = False) -> str:
    from components.theme import logo_data_uri

    size = 38 if compact else 46
    return (
        f"<div class='ct-brand'><img src='{logo_data_uri()}' alt='THE CHAOS THEORY orbital logo' "
        f"style='width:{size}px;height:{size}px'/>"
        "<div><div class='ct-brand-name'>THE CHAOS THEORY</div>"
        "<div class='ct-brand-sub'>Explore the Universe · Detect the Unexpected</div></div></div>"
    )


def default_status_items() -> list[dict]:
    summary = state.state_summary()
    archive = summary["archive_status"] or {}
    return [
        {"label": "Data source", "value": config.DATA_MODES.get(summary["mode_key"], "NO DATA"), "dot": "off" if summary["mode_key"] == "EMPTY" else ""},
        {"label": "IRSA link", "value": archive.get("status", "Not checked"), "dot": "" if archive.get("status") == "online" else "warn"},
        {"label": "Processing", "value": summary["processing_state"], "dot": "" if summary["comparison_ready"] else "off"},
        {"label": "Session", "value": summary["session_id"], "dot": ""},
    ]


def render_top_bar(title: str, subtitle: str, status_items: list[dict] | None = None) -> None:
    items = status_items if status_items is not None else default_status_items()
    cells = []
    for item in items:
        dot = item.get("dot", "")
        cells.append(
            "<div class='ct-status-item'>"
            f"<div class='ct-status-label'>{escape(str(item.get('label','')))}</div>"
            f"<div class='ct-status-value'><span class='ct-dot {escape(dot)}'></span>{escape(str(item.get('value','')))}</div>"
            "</div>"
        )
    st.html(
        "<div class='ct-topbar ct-page-enter'>"
        f"<div><div class='ct-topbar-title'>{escape(title)}</div>"
        f"<div class='ct-topbar-sub'>{escape(subtitle)}</div></div>"
        + "".join(cells)
        + "</div>"
    )


def render_footer() -> None:
    links = " · ".join(
        f"<a href='{escape(doc['url'])}' target='_blank'>{escape(doc['title'])}</a>"
        for doc in config.ARCHIVE_DOCS[:4]
    )
    st.html(
        "<div class='ct-footer'>"
        "<strong>THE CHAOS THEORY</strong> — a hackathon prototype for multi-epoch exploration of NASA SPHEREx observations.<br/>"
        "Data: NASA SPHEREx / NASA-IPAC IRSA. Not an official NASA product. Candidate detections are not discoveries.<br/>"
        f"<span class='ct-quiet'>Verified service documentation: {links}</span>"
        "</div>"
    )
