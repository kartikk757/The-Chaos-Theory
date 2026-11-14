"""Reusable astronomical image viewer controls and panel rendering."""
from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np
import streamlit as st

from modules import visualization
from modules.models import FitsProduct


def render_view_controls(key: str, default_colormap: str = "Inferno (science default)") -> dict[str, Any]:
    columns = st.columns([1.4, 1, 1, 1])
    with columns[0]:
        colormap = st.selectbox("Colormap", list(visualization.COLORMAPS), index=list(visualization.COLORMAPS).index(default_colormap), key=f"{key}_cmap")
    with columns[1]:
        low = st.slider("Contrast low (%)", 0.0, 40.0, 0.5, 0.5, key=f"{key}_low")
    with columns[2]:
        high = st.slider("Contrast high (%)", 60.0, 100.0, 99.5, 0.5, key=f"{key}_high")
    with columns[3]:
        invert = st.toggle("Invert (dark features)", value=False, key=f"{key}_invert")
    return {"colormap": visualization.COLORMAPS[colormap], "colormap_label": colormap, "low": low, "high": high, "invert": invert}


def _marker_for(product: FitsProduct, candidates: list[Any], selected: Any) -> dict | None:
    if selected is None:
        return None
    if hasattr(selected, "pixel_x"):
        return {"x": selected.pixel_x, "y": selected.pixel_y, "label": selected.candidate_id}
    return None


def render_image(
    array: np.ndarray,
    title: str,
    controls: dict,
    wcs=None,
    candidates: list[Any] | None = None,
    selected: Any = None,
    height: int = 460,
    show_candidate_markers: bool = True,
) -> None:
    if candidates and show_candidate_markers:
        figure = visualization.candidate_overlay_figure(array, candidates, wcs=wcs, title=title, height=height)
    else:
        figure = visualization.image_figure(
            array,
            title,
            controls["colormap"],
            controls["low"],
            controls["high"],
            wcs=wcs,
            dark=controls.get("invert", False),
            height=height,
            marker=_marker_for(None, candidates or [], selected),
        )
    st.plotly_chart(figure, use_container_width=True, theme=None, config={"displaylogo": False, "scrollZoom": True})


def render_fullscreen(array: np.ndarray, title: str, controls: dict, wcs=None) -> None:
    with st.expander("FULL-SCREEN VIEWER", expanded=False):
        st.caption("Use your browser's full-screen shortcut (F11) for a distraction-free view. Zoom with the mouse wheel; pan by dragging.")
        figure = visualization.image_figure(array, title, controls["colormap"], controls["low"], controls["high"], wcs=wcs, dark=controls.get("invert", False), height=820)
        st.plotly_chart(figure, use_container_width=True, theme=None, config={"displaylogo": False, "scrollZoom": True})


def render_pixel_inspector(product: FitsProduct, key: str, candidates: list[Any] | None = None) -> None:
    rows, columns = product.shape
    st.markdown("<div class='ct-mini-label'>PIXEL INSPECTION</div>", unsafe_allow_html=True)
    default_x = int(product.center_pixel[0]) if 0 <= product.center_pixel[0] < columns else columns // 2
    default_y = int(product.center_pixel[1]) if 0 <= product.center_pixel[1] < rows else rows // 2
    col1, col2 = st.columns(2)
    with col1:
        x = st.number_input("x pixel", min_value=0, max_value=columns - 1, value=default_x, step=1, key=f"{key}_px_x")
    with col2:
        y = st.number_input("y pixel", min_value=0, max_value=rows - 1, value=default_y, step=1, key=f"{key}_px_y")
    xi, yi = int(x), int(y)
    value = float(product.image[yi, xi])
    science = float(product.science[yi, xi])
    variance = float(product.variance[yi, xi]) if product.variance is not None else None
    zodi = float(product.zodi[yi, xi]) if product.zodi is not None else None
    flag_value = int(product.flags[yi, xi]) if product.flags is not None else None
    active_flags = []
    if flag_value:
        for name, definition in product.flag_definitions.items():
            if flag_value & int(definition["bit_mask"]):
                active_flags.append(name)
    sky = None
    if product.wcs is not None:
        try:
            world = product.wcs.pixel_to_world(xi, yi)
            sky = (float(world.ra.deg), float(world.dec.deg))
        except Exception:
            sky = None
    metrics = st.columns(4)
    metrics[0].metric("IMAGE value", f"{value:.4g}", help=f"Units: {product.bunit or 'not supplied'}")
    metrics[1].metric("IMAGE − ZODI", f"{science:.4g}" if product.zodi is not None else "n/a", help="Zodiacal model subtracted (documented product convention)")
    metrics[2].metric("VARIANCE", f"{variance:.4g}" if variance is not None else "not supplied", help="Units: (IMAGE units)²")
    metrics[3].metric("ZODI", f"{zodi:.4g}" if zodi is not None else "not supplied")
    if sky:
        st.caption(f"ICRS sky position of pixel ({xi}, {yi}): RA {sky[0]:.6f}°, Dec {sky[1]:.6f}°")
    else:
        st.caption("Spatial WCS is unavailable for this pixel.")
    if flag_value is None:
        st.caption("FLAGS extension not present in this file.")
    else:
        st.caption(f"FLAGS raw value {flag_value} · named bits set: {', '.join(sorted(active_flags)) if active_flags else 'none'}")
        severe = [name for name in active_flags if name in set(product.flag_definitions) and name not in ("SOURCE", "FULLSAMPLE")]
        if severe:
            st.warning(f"Severe processing flags present at this pixel: {', '.join(severe)}. Treat measurements here with caution.")


def render_download_buttons(product: FitsProduct, key: str, candidates_csv: str | None = None) -> None:
    columns = st.columns(2 if candidates_csv else 1)
    path = Path(product.path)
    if path.exists():
        with columns[0]:
            st.download_button(
                "DOWNLOAD FITS CUTOUT",
                data=path.read_bytes(),
                file_name=path.name,
                mime="application/fits",
                key=f"{key}_dl_fits",
                help="Download the exact FITS file used for these measurements.",
                use_container_width=True,
            )
    if candidates_csv:
        with columns[1]:
            st.download_button(
                "EXPORT CANDIDATES CSV",
                data=candidates_csv,
                file_name="the_chaos_theory_candidates.csv",
                mime="text/csv",
                key=f"{key}_dl_csv",
                use_container_width=True,
            )


def render_blink(product1: FitsProduct, product2: FitsProduct, controls: dict, key: str, interval: float = 1.2) -> None:
    """Automatic blink comparator using a Streamlit fragment timer."""
    st.markdown("<div class='ct-mini-label'>BLINK COMPARATOR</div>", unsafe_allow_html=True)
    start = st.toggle("Start automatic blink", value=False, key=f"{key}_blink_start", help=f"Alternates epochs every {interval:.1f} s")
    if not start:
        st.caption("Blink is stopped. Enable it to alternate the two epochs automatically; the frame label shows which epoch is displayed.")
        return

    @st.fragment(run_every=interval)
    def _blink_frame() -> None:
        counter = st.session_state.get(f"{key}_blink_counter", 0)
        st.session_state[f"{key}_blink_counter"] = counter + 1
        first = counter % 2 == 0
        product = product1 if first else product2
        date = product.metadata.get("observation", {}).get("date_utc", "unknown")
        st.markdown(
            f"<div class='ct-blink-label'>{'EPOCH 1' if first else 'EPOCH 2'} · {date} · {product.metadata.get('observation', {}).get('band', '')}</div>",
            unsafe_allow_html=True,
        )
        figure = visualization.image_figure(
            product.science, "", controls["colormap"], controls["low"], controls["high"],
            wcs=product.wcs, dark=controls.get("invert", False), height=440,
        )
        st.plotly_chart(figure, use_container_width=True, theme=None, key=f"{key}_blink_{counter % 2}", config={"displaylogo": False, "staticPlot": True})

    _blink_frame()


def render_overlay(product1: FitsProduct, product2: FitsProduct, controls: dict, key: str) -> None:
    st.markdown("<div class='ct-mini-label'>OPACITY OVERLAY</div>", unsafe_allow_html=True)
    opacity = st.slider("Epoch 2 opacity over epoch 1", 0.0, 1.0, 0.5, 0.05, key=f"{key}_overlay_opacity")
    first = visualization._stretch(product1.science, controls["low"], controls["high"])
    second = visualization._stretch(product2.science, controls["low"], controls["high"])
    blended = (1 - opacity) * first + opacity * second
    figure = visualization.image_figure(blended, f"Blend · epoch 1 {1 - opacity:.0%} / epoch 2 {opacity:.0%}", "Inferno", 0, 100, wcs=product1.wcs, height=460)
    st.plotly_chart(figure, use_container_width=True, theme=None, config={"displaylogo": False, "scrollZoom": True})
    st.caption("The overlay is a display blend of two normalized percentiles, not a physical composite; it is intended for visual registration checks only.")
