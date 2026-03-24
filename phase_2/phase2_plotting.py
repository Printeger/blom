"""
Plotting helpers for Phase 2 validation artifacts.
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np

from phase_2.phase2_validation import compute_segmentwise_influence_norms


def _prepare_path(save_path: str | Path) -> Path:
    path = Path(save_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    return path


def visualize_matrix_sparsity(mat: np.ndarray, save_path: str | Path, title: str) -> None:
    """Save a sparsity image for one matrix."""
    mat = np.asarray(mat, dtype=float)
    path = _prepare_path(save_path)
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.spy(mat, markersize=2, color="black")
    ax.set_title(title)
    ax.set_xlabel("Column")
    ax.set_ylabel("Row")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def visualize_heatmap(
    mat: np.ndarray,
    save_path: str | Path,
    title: str,
    log_scale: bool = True,
) -> None:
    """Save a heatmap of one matrix, optionally in log magnitude."""
    mat = np.asarray(mat, dtype=float)
    display = np.log10(np.abs(mat) + 1e-16) if log_scale else mat
    path = _prepare_path(save_path)
    fig, ax = plt.subplots(figsize=(7, 5))
    image = ax.imshow(display, aspect="auto", cmap="viridis", origin="lower")
    ax.set_title(title)
    ax.set_xlabel("Waypoint Index")
    ax.set_ylabel("Coefficient Row")
    colorbar = fig.colorbar(image, ax=ax)
    colorbar.set_label("log10(|value| + 1e-16)" if log_scale else "value")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def visualize_block_sparsity_Jq(J_q: np.ndarray, M_seg: int, save_path: str | Path) -> None:
    """Save the segment-by-waypoint block influence heatmap."""
    influence = compute_segmentwise_influence_norms(J_q, M_seg)
    path = _prepare_path(save_path)
    fig, ax = plt.subplots(figsize=(6, 5))
    image = ax.imshow(np.log10(influence + 1e-16), aspect="auto", cmap="magma", origin="lower")
    ax.set_title(f"Block Influence Heatmap, M={M_seg}")
    ax.set_xlabel("Interior Waypoint Index")
    ax.set_ylabel("Segment Index")
    colorbar = fig.colorbar(image, ax=ax)
    colorbar.set_label("log10(||dc_i / dq_j||_2 + 1e-16)")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def visualize_waypoint_influence_profile(
    profile: dict,
    save_path: str | Path,
    title: str,
) -> None:
    """Save the influence-decay curve for one interior waypoint."""
    path = _prepare_path(save_path)
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(profile["segment_indices"], profile["influence_norms"], marker="o", linewidth=1.8)
    ax.set_title(title)
    ax.set_xlabel("Segment Index")
    ax.set_ylabel("||dc_i / dq_j||_2")
    ax.grid(True, alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def visualize_jacobian_fd_compare(
    J_exact: np.ndarray,
    J_fd: np.ndarray,
    save_path: str | Path,
    title: str,
) -> None:
    """Save a comparison figure between exact and FD Jacobian entries."""
    J_exact = np.asarray(J_exact, dtype=float)
    J_fd = np.asarray(J_fd, dtype=float)
    path = _prepare_path(save_path)

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    axes[0].scatter(J_exact.reshape(-1), J_fd.reshape(-1), s=6, alpha=0.6)
    combined = np.concatenate((J_exact.reshape(-1), J_fd.reshape(-1))) if J_exact.size else np.array([0.0])
    lower = float(np.min(combined))
    upper = float(np.max(combined))
    axes[0].plot([lower, upper], [lower, upper], linestyle="--", color="black", linewidth=1.0)
    axes[0].set_title(f"{title}: Scatter")
    axes[0].set_xlabel("Exact Jacobian Entry")
    axes[0].set_ylabel("FD Jacobian Entry")

    diff = J_exact - J_fd
    image = axes[1].imshow(np.log10(np.abs(diff) + 1e-16), aspect="auto", cmap="viridis", origin="lower")
    axes[1].set_title(f"{title}: Error Heatmap")
    axes[1].set_xlabel("Waypoint Index")
    axes[1].set_ylabel("Coefficient Row")
    colorbar = fig.colorbar(image, ax=axes[1])
    colorbar.set_label("log10(|exact - fd| + 1e-16)")

    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def visualize_effective_bandwidth_vs_M(
    records: list[dict],
    save_path: str | Path,
    title: str,
) -> None:
    """Save the scaling plot of effective bandwidth statistics."""
    path = _prepare_path(save_path)
    M_values = np.asarray([record["M"] for record in records], dtype=float)
    max_bandwidth = np.asarray([record["max_effective_bandwidth"] for record in records], dtype=float)
    mean_bandwidth = np.asarray([record["mean_effective_bandwidth"] for record in records], dtype=float)
    far_ratio = np.asarray([record["far_nonzero_ratio"] for record in records], dtype=float)

    fig, ax1 = plt.subplots(figsize=(7, 4.5))
    ax1.plot(M_values, max_bandwidth, marker="o", linewidth=1.8, label="max effective bandwidth")
    ax1.plot(M_values, mean_bandwidth, marker="s", linewidth=1.8, label="mean effective bandwidth")
    ax1.set_title(title)
    ax1.set_xlabel("M")
    ax1.set_ylabel("Bandwidth")
    ax1.grid(True, alpha=0.3)

    ax2 = ax1.twinx()
    ax2.plot(M_values, far_ratio, marker="^", linestyle="--", color="tab:red", label="far-field ratio")
    ax2.set_ylabel("Far-field ratio")

    handles1, labels1 = ax1.get_legend_handles_labels()
    handles2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(handles1 + handles2, labels1 + labels2, loc="best")

    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)

