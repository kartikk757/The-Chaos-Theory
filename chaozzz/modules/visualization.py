"""Plotly and Matplotlib rendering for scientific arrays and sky geometry.

Rendering only: no archive access, no Streamlit widgets, no state mutations.
"""
from __future__ import annotations

import io
from typing import Any

import numpy as np
import plotly.graph_objects as go
from astropy.visualization import AsinhStretch, ImageNormalize, ZScaleInterval
from matplotlib.figure import Figure

import config

PLOTLY_THEME = {
    "paper_bgcolor": "rgba(0,0,0,0)",
    "plot_bgcolor": "rgba(4,9,20,0.65)",
    "font": {"color": "#c8d6e5", "family": "IBM Plex Mono, monospace", "size": 11},
    "margin": {"l": 56, "r": 16, "t": 34, "b": 46},
}

COLORMAPS = {
    "Inferno (science default)": "Inferno",
    "Viridis": "Viridis",
    "Cividis": "Cividis",
    "Magma": "Magma",
    "Plasma": "Plasma",
    "Greys": "Greys",
    "Electric cyan": "ice",
    "Ultraviolet": "plotly3",
}


def _stretch(array: np.ndarray, low: float, high: float) -> np.ndarray:
    data = np.asarray(array, dtype=float)
    finite = np.isfinite(data)
    if not finite.any():
        return np.zeros_like(data)
    values = data[finite]
    vmin = np.nanpercentile(values, max(0.0, low))
    vmax = np.nanpercentile(values, min(99.99, high))
    if not np.isfinite(vmin) or not np.isfinite(vmax) or vmax <= vmin:
        vmin, vmax = float(np.nanmin(values)), float(np.nanmax(values))
    if vmax <= vmin:
        vmax = vmin + 1e-30
    return np.clip((data - vmin) / (vmax - vmin), 0.0, 1.0)


def image_figure(
    array: np.ndarray,
    title: str,
    colormap: str = "Inferno",
    low: float = 0.5,
    high: float = 99.5,
    wcs=None,
    dark: bool = False,
    height: int = 460,
    marker: dict | None = None,
) -> go.Figure:
    data = np.asarray(array, dtype=float)
    scaled = _stretch(data, low, high)
    if dark:
        scaled = 1.0 - scaled
    figure = go.Figure()
    if wcs is not None and getattr(wcs, "has_celestial", False):
        try:
            ny, nx = scaled.shape
            x_pixels = np.arange(nx, dtype=float)
            y_pixels = np.arange(ny, dtype=float)
            centre_y = (ny - 1) / 2.0
            centre_x = (nx - 1) / 2.0
            x_world = wcs.pixel_to_world(x_pixels, np.full(nx, centre_y)).icrs.ra.deg
            y_world = wcs.pixel_to_world(np.full(ny, centre_x), y_pixels).icrs.dec.deg
            x_axis = dict(title="RA (deg)", color="#8fa6c4")
            y_axis = dict(title="Dec (deg)", color="#8fa6c4")
            figure.add_trace(
                go.Heatmap(
                    z=scaled,
                    x=np.asarray(x_world, dtype=float),
                    y=np.asarray(y_world, dtype=float),
                    colorscale=colormap,
                    showscale=False,
                    customdata=data,
                    hovertemplate="RA %{x:.4f}°<br>Dec %{y:.4f}°<br>value %{customdata:.4g}<extra></extra>",
                    zmin=0,
                    zmax=1,
                )
            )
            figure.update_layout(xaxis=x_axis, yaxis=y_axis)
            figure.update_xaxes(autorange="reversed")
            figure.update_yaxes(scaleanchor="x", scaleratio=1)
        except Exception:
            wcs = None
    if wcs is None:
        figure.add_trace(
            go.Heatmap(
                z=scaled,
                colorscale=colormap,
                showscale=False,
                customdata=data,
                hovertemplate="x %{x} px<br>y %{y} px<br>value %{customdata:.4g}<extra></extra>",
                zmin=0,
                zmax=1,
            )
        )
        figure.update_xaxes(title="x (pixel)", color="#8fa6c4")
        figure.update_yaxes(title="y (pixel)", color="#8fa6c4", scaleanchor="x", scaleratio=1)
    if marker:
        marker_x, marker_y = marker.get("x"), marker.get("y")
        if wcs is not None and getattr(wcs, "has_celestial", False) and marker_x is not None and marker_y is not None:
            try:
                world = wcs.pixel_to_world(float(marker_x), float(marker_y)).icrs
                marker_x, marker_y = float(world.ra.deg), float(world.dec.deg)
            except Exception:
                pass
        figure.add_trace(
            go.Scatter(
                x=[marker_x],
                y=[marker_y],
                mode="markers",
                marker=dict(size=22, color="rgba(0,0,0,0)", line=dict(color="#46f0ff", width=2.5)),
                hoverinfo="skip",
                showlegend=False,
            )
        )
        figure.add_trace(
            go.Scatter(
                x=[marker.get("x")],
                y=[marker.get("y")],
                mode="markers+text",
                marker=dict(size=7, color="#46f0ff"),
                text=[marker.get("label", "")],
                textposition="top center",
                textfont=dict(color="#46f0ff", size=10),
                showlegend=False,
                hoverinfo="skip",
            )
        )
    figure.update_layout(title=dict(text=title, font=dict(size=12, color="#e6f7ff")), height=height, **PLOTLY_THEME)
    return figure


def difference_figure(difference: np.ndarray, title: str, scale_sigma: float = 6.0, wcs=None, height: int = 460) -> go.Figure:
    data = np.asarray(difference, dtype=float)
    finite = data[np.isfinite(data)]
    if finite.size == 0:
        return image_figure(data, title, "RdBu", height=height)
    limit = float(np.nanpercentile(np.abs(finite), 99.0)) or float(np.nanmax(np.abs(finite)))
    limit = max(limit, scale_sigma * float(np.nanstd(finite)) / 2.0, 1e-12)
    figure = go.Figure(
        go.Heatmap(
            z=data,
            colorscale="RdBu",
            zmid=0,
            zmin=-limit,
            zmax=limit,
            colorbar=dict(title=dict(text="Δ MJy/sr", font=dict(color="#c8d6e5")), tickfont=dict(color="#8fa6c4")),
            hovertemplate="value %{z:.4g} MJy/sr<extra></extra>",
        )
    )
    figure.update_layout(title=dict(text=title, font=dict(size=12, color="#e6f7ff")), height=height, **PLOTLY_THEME)
    figure.update_xaxes(title="x (pixel)", color="#8fa6c4")
    figure.update_yaxes(title="y (pixel)", color="#8fa6c4", scaleanchor="x", scaleratio=1)
    return figure


def snr_figure(snr: np.ndarray, title: str, threshold: float = 5.0, height: int = 420) -> go.Figure:
    data = np.asarray(snr, dtype=float)
    figure = go.Figure(
        go.Heatmap(
            z=data,
            colorscale="Electric" if False else "Turbo",
            zmin=-max(threshold * 2, 5),
            zmax=max(threshold * 2, 5),
            colorbar=dict(title=dict(text="SNR", font=dict(color="#c8d6e5")), tickfont=dict(color="#8fa6c4")),
            hovertemplate="SNR %{z:.2f}<extra></extra>",
        )
    )
    figure.update_layout(title=dict(text=title, font=dict(size=12, color="#e6f7ff")), height=height, **PLOTLY_THEME)
    figure.update_yaxes(scaleanchor="x", scaleratio=1)
    return figure


def cutout_figure(array: np.ndarray, title: str, low: float = 1.0, high: float = 99.0, height: int = 260) -> go.Figure:
    figure = go.Figure(
        go.Heatmap(
            z=_stretch(array, low, high),
            colorscale="Inferno",
            showscale=False,
            customdata=np.asarray(array, dtype=float),
            hovertemplate="value %{customdata:.4g}<extra></extra>",
            zmin=0,
            zmax=1,
        )
    )
    figure.update_layout(title=dict(text=title, font=dict(size=11, color="#e6f7ff")), height=height, **PLOTLY_THEME)
    figure.update_yaxes(scaleanchor="x", scaleratio=1)
    figure.update_xaxes(visible=False)
    figure.update_yaxes(visible=False)
    return figure


def sky_footprint_figure(
    center_ra: float,
    center_dec: float,
    radius_deg: float,
    coverage_rows: list[dict],
    observations: list[dict],
    height: int = 520,
) -> go.Figure:
    figure = go.Figure()
    for index, row in enumerate(coverage_rows):
        polygon = row.get("polygon") or []
        if len(polygon) < 3:
            continue
        ras = [point[0] for point in polygon] + [polygon[0][0]]
        decs = [point[1] for point in polygon] + [polygon[0][1]]
        figure.add_trace(
            go.Scatter(
                x=ras,
                y=decs,
                mode="lines",
                line=dict(color="rgba(70,240,255,0.55)", width=1.2),
                fill="toself",
                fillcolor="rgba(70,240,255,0.05)",
                name=str(row.get("collection", "footprint")),
                hovertemplate=f"{row.get('dataproduct_type','image')}<br>t_min {row.get('t_min','')}<br>{row.get('energy_bandpassname','')}<extra></extra>",
                showlegend=index < 6,
            )
        )
    theta = np.linspace(0, 2 * np.pi, 120)
    figure.add_trace(
        go.Scatter(
            x=center_ra + radius_deg * np.cos(theta) / max(np.cos(np.radians(center_dec)), 1e-3),
            y=center_dec + radius_deg * np.sin(theta),
            mode="lines",
            line=dict(color="#b07cff", width=2, dash="dot"),
            name="requested search radius",
        )
    )
    figure.add_trace(
        go.Scatter(
            x=[center_ra],
            y=[center_dec],
            mode="markers+text",
            marker=dict(size=11, color="#b07cff", symbol="x"),
            text=["requested centre"],
            textposition="top center",
            textfont=dict(color="#d9c2ff", size=10),
            name="requested centre",
        )
    )
    for obs in observations:
        if obs.get("ra_deg") is None or obs.get("dec_deg") is None:
            continue
        figure.add_trace(
            go.Scatter(
                x=[obs["ra_deg"]],
                y=[obs["dec_deg"]],
                mode="markers",
                marker=dict(size=6, color="#46f0ff", line=dict(color="#04101f", width=1)),
                name="selected observation centre",
                showlegend=False,
                hovertemplate=f"{obs.get('obs_id','')}<br>{obs.get('date_utc','')}<br>{obs.get('band','')}<extra></extra>",
            )
        )
    figure.update_layout(
        title=dict(text="Requested field vs. archive footprint coverage (ICRS)", font=dict(size=12, color="#e6f7ff")),
        height=height,
        showlegend=True,
        legend=dict(font=dict(size=9), bgcolor="rgba(0,0,0,0.3)"),
        **PLOTLY_THEME,
    )
    figure.update_xaxes(title="RA (deg, reversed)", autorange="reversed", color="#8fa6c4")
    figure.update_yaxes(title="Dec (deg)", color="#8fa6c4", scaleanchor="x", scaleratio=1 / max(np.cos(np.radians(center_dec)), 1e-3))
    return figure


def timeline_figure(observations: list[dict], height: int = 260) -> go.Figure:
    figure = go.Figure()
    points = [obs for obs in observations if obs.get("mjd") is not None]
    if not points:
        figure.update_layout(title=dict(text="No dated observations available", font=dict(size=12, color="#e6f7ff")), height=height, **PLOTLY_THEME)
        return figure
    import datetime as _dt

    dates = [_dt.datetime.utcfromtimestamp((float(obs["mjd"]) - 40587) * 86400) for obs in points]
    figure.add_trace(
        go.Scatter(
            x=dates,
            y=[obs.get("band", "") for obs in points],
            mode="markers",
            marker=dict(size=9, color="#46f0ff", line=dict(color="#0a1626", width=1)),
            text=[f"{obs.get('obs_id','')}<br>{obs.get('release','')}" for obs in points],
            hovertemplate="%{x|%Y-%m-%d %H:%M}<br>%{text}<extra></extra>",
            name="observation",
        )
    )
    figure.update_layout(
        title=dict(text="Observation epochs returned by the archive", font=dict(size=12, color="#e6f7ff")),
        height=height,
        **PLOTLY_THEME,
    )
    figure.update_xaxes(title="UTC", color="#8fa6c4")
    figure.update_yaxes(title="", color="#8fa6c4")
    return figure


def candidate_overlay_figure(
    base: np.ndarray,
    candidates: list[Any],
    wcs=None,
    title: str = "Residual candidates over epoch-1 image",
    height: int = 520,
) -> go.Figure:
    figure = image_figure(base, title, "Greys", wcs=wcs, height=height)
    colors = {
        "Potential moving-source candidate": "#46f0ff",
        "Potential variable-source candidate": "#b07cff",
        "Unmatched residual": "#ffd166",
        "Likely image artifact (user review required)": "#ff6b6b",
        "Unclassified detection": "#9bb1c9",
    }
    figure.add_trace(
        go.Scatter(
            x=[candidate.pixel_x for candidate in candidates],
            y=[candidate.pixel_y for candidate in candidates],
            mode="markers+text",
            marker=dict(
                size=13,
                color=[colors.get(candidate.classification, "#ffffff") for candidate in candidates],
                line=dict(color="#03101d", width=1.5),
                symbol="circle-open",
            ),
            text=[candidate.candidate_id for candidate in candidates],
            textposition="top center",
            textfont=dict(size=9, color="#dff6ff"),
            customdata=[[c.snr, c.classification, c.ra_deg, c.dec_deg] for c in candidates],
            hovertemplate="%{text}<br>SNR %{customdata[0]:.1f}<br>%{customdata[1]}<br>RA %{customdata[2]:.5f} Dec %{customdata[3]:.5f}<extra></extra>",
            name="candidates",
        )
    )
    return figure


def figure_to_png(figure: Figure, dpi: int = 130) -> bytes:
    buffer = io.BytesIO()
    figure.savefig(buffer, format="png", dpi=dpi, bbox_inches="tight", facecolor="#050b16")
    buffer.seek(0)
    return buffer.read()


def array_preview_png(array: np.ndarray, title: str = "", cmap: str = "inferno", width: float = 3.2) -> bytes:
    data = np.asarray(array, dtype=float)
    figure = Figure(figsize=(width, width * 0.92), facecolor="#050b16")
    axis = figure.add_subplot(111)
    axis.set_facecolor("#050b16")
    try:
        norm = ImageNormalize(data[np.isfinite(data)], interval=ZScaleInterval(), stretch=AsinhStretch())
        axis.imshow(data, cmap=cmap, norm=norm, origin="lower")
    except Exception:
        axis.imshow(data, cmap=cmap, origin="lower")
    axis.set_xticks([])
    axis.set_yticks([])
    if title:
        axis.set_title(title, color="#dff6ff", fontsize=8)
    return figure_to_png(figure)
