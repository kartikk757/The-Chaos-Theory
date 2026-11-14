"""Reusable glass metric cards. Values are supplied by real application state."""
from __future__ import annotations

from html import escape

import streamlit as st


def render_metric_grid(cards: list[dict]) -> None:
    html = ["<div class='ct-grid'>"]
    for card in cards:
        label = escape(str(card.get("label", "")))
        value = escape(str(card.get("value", "Awaiting data")))
        caption = escape(str(card.get("caption", "")))
        icon = escape(str(card.get("icon", "")))
        accent = " accent" if card.get("accent") else ""
        html.append(
            "<div class='ct-metric'>"
            f"<div class='ct-metric-label'>{label}</div>"
            f"<div class='ct-metric-icon'>{icon}</div>"
            f"<div class='ct-metric-value{accent}'>{value}</div>"
            f"<div class='ct-metric-caption'>{caption}</div>"
            "</div>"
        )
    html.append("</div>")
    st.html("".join(html))
