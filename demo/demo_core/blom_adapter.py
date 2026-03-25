from __future__ import annotations

import time
from typing import Any

import numpy as np
import streamlit as st

from phase_1.minco_scalar_baseline import evaluate_trajectory
from phase_5.blom_boundary_jump_check import assemble_scheme_A, assemble_scheme_B, assemble_scheme_C
from phase_6.blom_fd_jacobian_check import assembled_scheme_jacobians_fd
from phase_7.blom_convergence_vs_k import compute_actual_blom_k, compute_ideal_truncated_blom_k, compute_minco_reference
from phase_8.phase8_common import compute_phase8_matching_triplet, make_interior_sets, summarize_interior_matching
from phase_9.blom_backward_diff import compute_raw_schemeC_jacobians as compute_raw_schemeC_jacobians_k2
from phase_9.blom_space_time_opt_demo import default_obs_config as default_obs_minimal
from phase_9.blom_space_time_opt_demo import default_weights as default_weights_minimal
from phase_9.blom_space_time_opt_demo import run_minimal_optimization_demo
from phase_10.blom_full_backward_diff import (
    T_to_tau,
    compute_raw_schemeC_jacobians_general_k,
    default_bc_config,
    default_dyn_config,
    default_obs_config,
    default_reg_config,
    default_weights,
    evaluate_full_objective,
    full_backward_diff_dense,
    full_backward_diff_sparse,
)
from phase_10.blom_space_time_optimizer import run_space_time_optimization


def _sample_piecewise(coeffs: np.ndarray, T: np.ndarray, num_samples: int = 300) -> dict[str, np.ndarray]:
    T = np.asarray(T, dtype=float).reshape(-1)
    total_time = float(np.sum(T))
    t = np.linspace(0.0, total_time, num_samples)
    if t.size:
        t[-1] = np.nextafter(total_time, 0.0)
    y = np.asarray([evaluate_trajectory(coeffs, T, float(tt), order=0) for tt in t], dtype=float)
    return {"t": t, "y": y}


@st.cache_data(show_spinner=False)
def compute_actual(q_tuple: tuple[float, ...], T_tuple: tuple[float, ...], k: int, scheme: str = "C") -> dict[str, Any]:
    q = np.asarray(q_tuple, dtype=float)
    T = np.asarray(T_tuple, dtype=float)
    data = compute_actual_blom_k(q, T, k=k, scheme=scheme)
    return {
        "mode": "strict",
        "coeffs": np.asarray(data["coeffs"], dtype=float),
        "c_vec": np.asarray(data["c_vec"], dtype=float),
        "cost": float(data["cost"]),
        "meta": data.get("meta", {}),
        "stats": data.get("stats", {}),
        "stats_overview": data.get("stats_overview", {}),
        "raw": data,
    }


@st.cache_data(show_spinner=False)
def compute_ideal(q_tuple: tuple[float, ...], T_tuple: tuple[float, ...], k: int) -> dict[str, Any]:
    q = np.asarray(q_tuple, dtype=float)
    T = np.asarray(T_tuple, dtype=float)
    reference = compute_minco_reference(q, T)
    data = compute_ideal_truncated_blom_k(q, T, k=k, reference=reference)
    return {"mode": "strict", **data}


def sample_actual(coeffs: np.ndarray, T: np.ndarray, num_samples: int = 300) -> dict[str, np.ndarray]:
    return _sample_piecewise(coeffs, T, num_samples=num_samples)


@st.cache_data(show_spinner=False)
def compute_raw_jacobians(q_tuple: tuple[float, ...], T_tuple: tuple[float, ...], k: int) -> dict[str, Any]:
    q = np.asarray(q_tuple, dtype=float)
    T = np.asarray(T_tuple, dtype=float)
    if k == 2:
        data = compute_raw_schemeC_jacobians_k2(q, T, k=k, mode="analytic")
    else:
        data = compute_raw_schemeC_jacobians_general_k(q, T, k=k, mode="analytic")
    return data


@st.cache_data(show_spinner=False)
def compute_scheme_fd_jacobian(q_tuple: tuple[float, ...], T_tuple: tuple[float, ...], scheme: str) -> dict[str, Any]:
    q = np.asarray(q_tuple, dtype=float)
    T = np.asarray(T_tuple, dtype=float)
    data = assembled_scheme_jacobians_fd(q, T, scheme=scheme)
    return data


def assemble_scheme_compare(q: np.ndarray, T: np.ndarray, k: int) -> dict[str, Any]:
    outputs = {
        "A": assemble_scheme_A(q, T, k=k),
        "B": assemble_scheme_B(q, T, k=k),
        "C": assemble_scheme_C(q, T, k=k),
    }
    return outputs


def compute_matching_bridge(q: np.ndarray, T: np.ndarray, k: int, radius_mode: str = "default") -> dict[str, Any]:
    triplet = compute_phase8_matching_triplet(q, T, k=k)
    actual_vec = np.asarray(triplet["actual"]["c_vec"], dtype=float)
    ideal_vec = np.asarray(triplet["ideal"]["c_vec"], dtype=float)
    ref_vec = np.asarray(triplet["reference_window"]["c_vec"], dtype=float)
    full = summarize_interior_matching(actual_vec, ideal_vec, T.size, k, radius_mode=radius_mode)
    ref_gap = summarize_interior_matching(ref_vec, ideal_vec, T.size, k, radius_mode=radius_mode)
    raw_ref = summarize_interior_matching(actual_vec, ref_vec, T.size, k, radius_mode=radius_mode)
    return {
        "triplet": triplet,
        "actual_vs_ideal": full,
        "reference_vs_ideal": ref_gap,
        "actual_vs_reference": raw_ref,
        "interior_sets": make_interior_sets(T.size, k, radius_mode=radius_mode),
    }


@st.cache_data(show_spinner=True)
def run_optimizer_cached(
    q_tuple: tuple[float, ...],
    T_tuple: tuple[float, ...],
    objective_mode: str,
    n_steps: int,
    step_size: float,
    k: int,
) -> dict[str, Any]:
    q = np.asarray(q_tuple, dtype=float)
    T = np.asarray(T_tuple, dtype=float)
    if objective_mode == "minimal":
        result = run_minimal_optimization_demo(
            q,
            T,
            default_weights_minimal(),
            obs_config=default_obs_minimal(),
            n_steps=n_steps,
            step_size=step_size,
            T_min=0.2,
        )
        before = _sample_piecewise(np.asarray(result["coeffs_before"], dtype=float), np.asarray(T, dtype=float))
        after = _sample_piecewise(np.asarray(result["coeffs_after"], dtype=float), np.asarray(result["T_history"][-1], dtype=float))
        return {"mode": "minimal", "result": result, "before": before, "after": after}

    tau = T_to_tau(T, T_min=0.2)
    result = run_space_time_optimization(
        q,
        tau,
        default_weights(),
        obs_config=default_obs_config(),
        dyn_config=default_dyn_config(),
        bc_config=default_bc_config(q),
        reg_config=default_reg_config(),
        n_steps=n_steps,
        step_size=step_size,
        T_min=0.2,
        k=k,
    )
    before_coeffs = compute_actual(tuple(q.tolist()), tuple(T.tolist()), k=k, scheme="C")["coeffs"]
    before = _sample_piecewise(before_coeffs, T)
    after = _sample_piecewise(np.asarray(result["coeffs_final"], dtype=float), np.asarray(result["T_final"], dtype=float))
    dense = full_backward_diff_dense(q, T, default_weights(), obs_config=default_obs_config(), dyn_config=default_dyn_config(), bc_config=default_bc_config(q), reg_config=default_reg_config(), k=k)
    sparse = full_backward_diff_sparse(q, T, default_weights(), obs_config=default_obs_config(), dyn_config=default_dyn_config(), bc_config=default_bc_config(q), reg_config=default_reg_config(), k=k)
    return {"mode": "full", "result": result, "before": before, "after": after, "dense": dense, "sparse": sparse}


def runtime_actual(q: np.ndarray, T: np.ndarray, k: int, scheme: str = "C") -> tuple[dict[str, Any], float]:
    start = time.perf_counter()
    data = compute_actual(tuple(q.tolist()), tuple(T.tolist()), k=k, scheme=scheme)
    return data, time.perf_counter() - start
