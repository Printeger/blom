"""
Plotting helpers for Phase 4 scripts.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


def _prepare(save_path: str | Path) -> Path:
    path = Path(save_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def save_line_plot(
    x: np.ndarray,
    y: np.ndarray,
    save_path: str | Path,
    *,
    title: str,
    xlabel: str,
    ylabel: str,
    logx: bool = False,
    logy: bool = False,
    marker: str = "o",
) -> None:
    path = _prepare(save_path)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.plot(x, y, marker=marker, linewidth=1.8)
    if logx:
        ax.set_xscale("log")
    if logy:
        ax.set_yscale("log")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def save_scatter_with_diagonal(
    x: np.ndarray,
    y: np.ndarray,
    save_path: str | Path,
    *,
    title: str,
    xlabel: str,
    ylabel: str,
) -> None:
    path = _prepare(save_path)
    fig, ax = plt.subplots(figsize=(5.5, 5.5))
    ax.scatter(x, y, s=20, alpha=0.7)
    all_values = np.concatenate((np.asarray(x).reshape(-1), np.asarray(y).reshape(-1)))
    lower = float(np.min(all_values))
    upper = float(np.max(all_values))
    ax.plot([lower, upper], [lower, upper], linestyle="--", color="black", linewidth=1.0)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def save_histogram(
    values: np.ndarray,
    save_path: str | Path,
    *,
    title: str,
    xlabel: str,
    ylabel: str = "Count",
) -> None:
    path = _prepare(save_path)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    ax.hist(np.asarray(values, dtype=float).reshape(-1), bins=20, color="#4c78a8", edgecolor="black")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.25)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def save_heatmap(
    mat: np.ndarray,
    save_path: str | Path,
    *,
    title: str,
    xlabel: str,
    ylabel: str,
) -> None:
    path = _prepare(save_path)
    fig, ax = plt.subplots(figsize=(6, 5))
    image = ax.imshow(np.asarray(mat, dtype=float), cmap="viridis", aspect="auto")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    fig.colorbar(image, ax=ax)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def save_overlay_curves(
    x: np.ndarray,
    curves: list[tuple[str, np.ndarray]],
    save_path: str | Path,
    *,
    title: str,
    xlabel: str,
    ylabel: str,
) -> None:
    path = _prepare(save_path)
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    for label, values in curves:
        ax.plot(x, values, linewidth=1.8, label=label)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, alpha=0.3)
    ax.legend(loc="best")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def save_contour(
    X: np.ndarray,
    Y: np.ndarray,
    Z: np.ndarray,
    optimum: tuple[float, float],
    save_path: str | Path,
    *,
    title: str,
    xlabel: str,
    ylabel: str,
) -> None:
    path = _prepare(save_path)
    fig, ax = plt.subplots(figsize=(6.5, 5.0))
    contour = ax.contourf(X, Y, Z, levels=25, cmap="viridis")
    ax.scatter([optimum[0]], [optimum[1]], color="red", s=50, label="optimum")
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.legend(loc="best")
    fig.colorbar(contour, ax=ax)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def save_bar_chart(
    labels: list[str],
    values: np.ndarray,
    save_path: str | Path,
    *,
    title: str,
    xlabel: str,
    ylabel: str,
) -> None:
    path = _prepare(save_path)
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    positions = np.arange(len(labels))
    ax.bar(positions, np.asarray(values, dtype=float), color="#4c78a8", edgecolor="black")
    ax.set_xticks(positions, labels=labels)
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.grid(True, axis="y", alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)
