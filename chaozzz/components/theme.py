"""Loads the fully self-contained design system (no external fonts/CDNs)."""
from __future__ import annotations

import base64
from functools import lru_cache
from pathlib import Path

import streamlit as st

import config


@lru_cache(maxsize=1)
def _styles() -> str:
    return (config.ASSETS_DIR / "styles.css").read_text(encoding="utf-8")


@lru_cache(maxsize=1)
def _image_data_uri(path: str) -> str:
    data = Path(path).read_bytes()
    suffix = Path(path).suffix.lower().lstrip(".")
    mime = "image/jpeg" if suffix in {"jpg", "jpeg"} else "image/png"
    return f"data:{mime};base64,{base64.b64encode(data).decode('ascii')}"


def logo_data_uri() -> str:
    return _image_data_uri(str(config.ASSETS_DIR / "logo.png"))


def background_data_uri() -> str:
    return _image_data_uri(str(config.ASSETS_DIR / "background.png"))


def load_styles(landing: bool = False) -> None:
    css = _styles()
    hero = background_data_uri()
    landing_css = """
    [data-testid="stSidebar"] { display: none !important; }
    .block-container { max-width: 1240px !important; padding-top: 1.2rem !important; }
    .ct-hero-visual { opacity: .96; }
    """
    st.markdown(
        f"<style>{css}\n:root {{ --ct-hero: url('{hero}'); }}\n{landing_css if landing else ''}</style>",
        unsafe_allow_html=True,
    )
    body = "<div class='ct-app-bg'></div><div class='ct-starfield'></div>"
    if landing:
        body += "<div class='ct-hero-visual'></div><div class='ct-orbital'></div>"
    st.html(body)
