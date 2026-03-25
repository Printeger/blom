from __future__ import annotations

from typing import Any

import numpy as np
import pandas as pd
import plotly.graph_objects as go


PAPER_COLORS = {
    "blom": "#0B5FFF",
    "minco": "#D1495B",
    "bspline": "#2F9E44",
    "schemeA": "#9C36B5",
    "schemeB": "#E67700",
    "schemeC": "#0B7285",
    "diff": "#6C757D",
}


def trajectory_figure(curves: list[dict[str, Any]], title: str, x_label: str = "parameter", y_label: str = "value") -> go.Figure:
    fig = go.Figure()
    for curve in curves:
        fig.add_trace(
            go.Scatter(
                x=curve["x"],
                y=curve["y"],
                mode="lines",
                name=curve["name"],
                line=dict(color=curve.get("color"), width=curve.get("width", 2.2), dash=curve.get("dash", "solid")),
            )
        )
    fig.update_layout(
        title=title,
        template="plotly_white",
        xaxis_title=x_label,
        yaxis_title=y_label,
        legend_title_text="Curve",
        margin=dict(l=30, r=30, t=60, b=30),
    )
    return fig


def difference_heatmap(x: np.ndarray, base: np.ndarray, perturbed: np.ndarray, title: str) -> go.Figure:
    diff = np.asarray(perturbed, dtype=float) - np.asarray(base, dtype=float)
    fig = go.Figure(
        data=go.Heatmap(
            z=diff[None, :],
            x=x,
            y=["difference"],
            colorscale="RdBu",
            zmid=0.0,
            colorbar=dict(title="delta"),
        )
    )
    fig.update_layout(title=title, template="plotly_white", xaxis_title="parameter", yaxis_title="")
    return fig


def matrix_heatmap(J: np.ndarray, title: str, x_label: str, y_label: str) -> go.Figure:
    mat = np.asarray(J, dtype=float)
    fig = go.Figure(
        data=go.Heatmap(
            z=np.log10(np.abs(mat) + 1e-16),
            colorscale="Viridis",
            colorbar=dict(title="log10(|value|+1e-16)"),
        )
    )
    fig.update_layout(title=title, template="plotly_white", xaxis_title=x_label, yaxis_title=y_label)
    return fig


def binary_mask_heatmap(mask: np.ndarray, title: str, x_label: str, y_label: str) -> go.Figure:
    fig = go.Figure(data=go.Heatmap(z=np.asarray(mask, dtype=float), colorscale="Greys", showscale=False))
    fig.update_layout(title=title, template="plotly_white", xaxis_title=x_label, yaxis_title=y_label)
    return fig


def line_metric_figure(x: list[Any], series: list[dict[str, Any]], title: str, x_label: str, y_label: str) -> go.Figure:
    fig = go.Figure()
    for item in series:
        fig.add_trace(
            go.Scatter(
                x=x,
                y=item["y"],
                mode="lines+markers",
                name=item["name"],
                line=dict(color=item.get("color"), width=2.4),
            )
        )
    fig.update_layout(title=title, template="plotly_white", xaxis_title=x_label, yaxis_title=y_label)
    return fig


def bar_metric_figure(labels: list[str], values: list[float], title: str, y_label: str, color: str = "#0B5FFF") -> go.Figure:
    fig = go.Figure(data=go.Bar(x=labels, y=values, marker_color=color))
    fig.update_layout(title=title, template="plotly_white", xaxis_title="", yaxis_title=y_label)
    return fig


def dataframe_download(df: pd.DataFrame) -> bytes:
    return df.to_csv(index=False).encode("utf-8")
