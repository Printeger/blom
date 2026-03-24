"""
Plotting helpers for Phase 3 BLOM-Strict validation artifacts.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np


def _prepare_path(save_path: str | Path) -> Path:
    path = Path(save_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def plot_window_layout(M: int, window: dict, save_path: str | Path) -> None:
    path = _prepare_path(save_path)
    fig, ax = plt.subplots(figsize=(8, 1.8))
    for segment in range(1, M + 1):
        color = "#d0d7de"
        if segment in window["segments"]:
            color = "#6aa84f"
        if segment == window["center"]:
            color = "#f6b26b"
        rect = plt.Rectangle((segment - 1, 0.0), 1.0, 1.0, facecolor=color, edgecolor="black")
        ax.add_patch(rect)
        ax.text(segment - 0.5, 0.5, f"{segment}", ha="center", va="center", fontsize=10)
    if window["touches_left_boundary"]:
        ax.text(-0.1, 1.1, "physical left", ha="left", va="bottom", fontsize=9, color="tab:blue")
    else:
        ax.text(window["L"] - 1.0, 1.1, "artificial left", ha="left", va="bottom", fontsize=9, color="tab:red")
    if window["touches_right_boundary"]:
        ax.text(M - 0.1, 1.1, "physical right", ha="right", va="bottom", fontsize=9, color="tab:blue")
    else:
        ax.text(window["R"], 1.1, "artificial right", ha="right", va="bottom", fontsize=9, color="tab:red")
    ax.set_xlim(0.0, float(M))
    ax.set_ylim(0.0, 1.35)
    ax.set_yticks([])
    ax.set_xticks(np.arange(0.5, M + 0.5, 1.0), labels=[str(j) for j in range(1, M + 1)])
    ax.set_xlabel("Global Segment Index")
    ax.set_title(f"Window W({window['center']},{window['k']}) = {window['segments']}")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_local_trajectory(
    samples: dict[str, np.ndarray],
    sigma_knots: np.ndarray,
    knot_values: np.ndarray,
    save_path: str | Path,
    title: str,
) -> None:
    path = _prepare_path(save_path)
    fig, axes = plt.subplots(5, 1, figsize=(9, 10), sharex=True)
    labels = ["p", "p'", "p''", "p'''", "p''''"]
    for order, ax in enumerate(axes):
        ax.plot(samples["sigma"], samples[f"order_{order}"], linewidth=1.8)
        for sigma in sigma_knots:
            ax.axvline(float(sigma), color="black", alpha=0.15, linewidth=1.0)
        if order == 0:
            ax.scatter(sigma_knots, knot_values, color="tab:red", zorder=3, label="interpolation knots")
            ax.legend(loc="best")
        ax.set_ylabel(labels[order])
        ax.grid(True, alpha=0.25)
    axes[0].set_title(title)
    axes[-1].set_xlabel("Local Time")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_continuity_jumps(jumps: np.ndarray, save_path: str | Path, title: str) -> None:
    path = _prepare_path(save_path)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    if jumps.size == 0:
        ax.text(0.5, 0.5, "No internal knots", ha="center", va="center")
        ax.set_axis_off()
    else:
        orders = np.arange(jumps.shape[1])
        for knot_idx in range(jumps.shape[0]):
            ax.plot(orders, np.abs(jumps[knot_idx]) + 1e-18, marker="o", label=f"knot {knot_idx + 1}")
        ax.set_yscale("log")
        ax.set_xlabel("Derivative Order r")
        ax.set_ylabel("|jump|")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="best")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_natural_bc_residuals(
    residuals: dict[str, np.ndarray],
    s: int,
    save_path: str | Path,
    title: str,
) -> None:
    path = _prepare_path(save_path)
    orders = np.arange(s, 2 * s - 1)
    fig, ax = plt.subplots(figsize=(7, 4.5))
    width = 0.35
    has_left = residuals["left"].size > 0
    has_right = residuals["right"].size > 0
    if has_left:
        ax.bar(orders - width / 2, residuals["left"] + 1e-18, width=width, label="left")
    else:
        ax.text(orders[0] - 0.6, 1e-16, "left: N/A", color="tab:blue")
    if has_right:
        ax.bar(orders + width / 2, residuals["right"] + 1e-18, width=width, label="right")
    else:
        ax.text(orders[-1] - 0.4, 1e-16, "right: N/A", color="tab:orange")
    ax.set_yscale("log")
    ax.set_xlabel("Derivative Order")
    ax.set_ylabel("Residual")
    ax.set_title(title)
    ax.grid(True, alpha=0.3)
    if has_left or has_right:
        ax.legend(loc="best")
    fig.subplots_adjust(left=0.12, right=0.88, bottom=0.15, top=0.88)
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_perturbation_response(
    epsilons: np.ndarray,
    coeff_diff_norms: np.ndarray,
    cond_numbers: np.ndarray,
    save_path: str | Path,
    title: str,
) -> None:
    path = _prepare_path(save_path)
    fig, ax1 = plt.subplots(figsize=(7, 4.5))
    ax1.loglog(epsilons, coeff_diff_norms, marker="o", linewidth=1.8, label="||delta c||_2")
    ax1.set_xlabel("Perturbation Magnitude eps")
    ax1.set_ylabel("Coefficient Difference Norm")
    ax1.grid(True, alpha=0.3)
    ax2 = ax1.twinx()
    ax2.semilogx(epsilons, cond_numbers, marker="s", linestyle="--", color="tab:red", label="cond(R_red)")
    ax2.set_ylabel("Condition Number")
    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(handles1 + handles2, labels1 + labels2, loc="best")
    ax1.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_uniqueness_multistart(distance_matrix: np.ndarray, save_path: str | Path, title: str) -> None:
    path = _prepare_path(save_path)
    fig, ax = plt.subplots(figsize=(5.5, 4.8))
    image = ax.imshow(np.log10(distance_matrix + 1e-18), cmap="viridis", origin="lower")
    ax.set_xlabel("Restart Index")
    ax.set_ylabel("Restart Index")
    ax.set_title(title)
    colorbar = fig.colorbar(image, ax=ax)
    colorbar.set_label("log10(pairwise coeff diff + 1e-18)")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)
