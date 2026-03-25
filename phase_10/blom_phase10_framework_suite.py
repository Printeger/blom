"""
Unified Phase 10 framework suite.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

from phase_10.blom_benchmark_suite import run_benchmark_suite
from phase_10.blom_full_backward_diff import (
    DEFAULT_RESULTS_DIR,
    T_to_tau,
    default_bc_config,
    default_dyn_config,
    default_obs_config,
    default_reg_config,
    default_weights,
    ensure_results_dirs,
    representative_case,
    run_full_backward_diff_demo,
    save_csv,
    save_json,
)
from phase_10.blom_space_time_optimizer import run_space_time_optimization


def _write_interpretation_summary(backward: dict[str, Any], optimizer: dict[str, Any], benchmark: dict[str, Any], save_path: str | Path) -> None:
    backward_summary = backward["summary"]
    optimizer_summary = optimizer["summary"]
    benchmark_rows = benchmark["benchmark_rows"]
    implemented_k_rows = [row for row in benchmark["overview_rows"] if row["method"] == "raw_schemeC" and row["k"] == 2]
    lines = [
        "# Phase 10 Interpretation Summary",
        "",
        f"- Full objective layer working end-to-end: `{optimizer_summary['objective_drop'] > 0.0}`.",
        f"- Dense and sparse backward agree: `{backward_summary['dense_vs_sparse_max_abs'] < 1e-8}`.",
        f"- Reparameterized gradient matches FD: `{backward_summary['reparam_vs_fd_max_abs'] < 1e-4}`.",
        f"- Raw Scheme C optimizer beats heuristic on representative case: `{min((row['final_objective'] for row in benchmark_rows if row['method']=='raw_schemeC'), default=float('inf')) < min((row['final_objective'] for row in benchmark_rows if row['method']=='heuristic'), default=float('inf'))}`.",
        f"- Any explicit general-k placeholder rows remain: `{any(row['k'] != 2 and row['method']=='raw_schemeC' and row['final_objective'] != row['final_objective'] for row in benchmark['overview_rows'])}`.",
        "",
        "Interpretation:",
        "- Phase 10 now elevates the raw Scheme C differentiable loop into a full objective-layer + optimizer + benchmark stack.",
        "- The dense checker stays in the framework as a safety rail, while the sparse path is the intended implementation route.",
        "- The benchmark layer is already sufficient to discuss speed-quality trade-offs, even if some general-k rows are still explicitly marked as not yet implemented.",
    ]
    Path(save_path).write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_phase10_framework_suite(
    q0: Any = None,
    tau0: Any = None,
    weights: dict[str, float] | None = None,
    obs_config: dict[str, Any] | None = None,
    dyn_config: dict[str, Any] | None = None,
    bc_config: dict[str, Any] | None = None,
    reg_config: dict[str, Any] | None = None,
    M_values: list[int] | None = None,
    k_values: list[int] | None = None,
    benchmark_methods: list[str] | None = None,
    n_steps: int = 40,
    T_min: float = 0.2,
    save_dir: str | Path | None = None,
    seed: int = 42,
) -> dict[str, Any]:
    q_rep, T_rep = representative_case() if q0 is None else (q0, None)
    if q0 is None:
        q_rep, T_rep = representative_case()
    else:
        q_rep = q0
        T_rep = None
    if tau0 is None:
        if T_rep is None:
            raise ValueError("When q0 is provided explicitly, tau0 must also be provided or T must be recoverable.")
        tau_rep = T_to_tau(T_rep, T_min=T_min)
        T0 = T_rep
    else:
        tau_rep = tau0
        T0 = None
    if T0 is None:
        from phase_10.blom_full_backward_diff import tau_to_T

        T0 = tau_to_T(tau_rep, T_min=T_min)

    weights = {**default_weights(), **({} if weights is None else dict(weights))}
    obs_cfg = default_obs_config() if obs_config is None else {**default_obs_config(), **dict(obs_config)}
    dyn_cfg = default_dyn_config() if dyn_config is None else {**default_dyn_config(), **dict(dyn_config)}
    bc_cfg = default_bc_config(q_rep) if bc_config is None else {**default_bc_config(q_rep), **dict(bc_config)}
    reg_cfg = default_reg_config() if reg_config is None else {**default_reg_config(), **dict(reg_config)}
    results_dirs = ensure_results_dirs(DEFAULT_RESULTS_DIR if save_dir is None else save_dir)

    backward = run_full_backward_diff_demo(q=q_rep, T=T0, weights=weights, obs_config=obs_cfg, dyn_config=dyn_cfg, bc_config=bc_cfg, reg_config=reg_cfg, k_values=k_values, save_dir=results_dirs["base"], seed=seed)
    optimizer = run_space_time_optimization(q_rep, tau_rep, weights, obs_config=obs_cfg, dyn_config=dyn_cfg, bc_config=bc_cfg, reg_config=reg_cfg, n_steps=n_steps, T_min=T_min, save_dir=results_dirs["base"], seed=seed)
    benchmark = run_benchmark_suite(q0=q_rep, T0=T0, weights=weights, obs_config=obs_cfg, dyn_config=dyn_cfg, bc_config=bc_cfg, reg_config=reg_cfg, M_values=M_values, k_values=k_values, benchmark_methods=benchmark_methods, n_steps=max(20, n_steps // 2), T_min=T_min, save_dir=results_dirs["base"], seed=seed)

    overview_rows = [
        ["backward_diff", backward["summary"]["dense_vs_sparse_max_abs"], backward["summary"]["reparam_vs_fd_max_abs"], 0.0, 0.0],
        ["optimizer", 0.0, 0.0, optimizer["summary"]["objective_drop"], optimizer["summary"]["runtime_total"]],
        ["benchmark", 0.0, 0.0, min(row["final_objective"] for row in benchmark["benchmark_rows"]), sum(row["runtime_total"] for row in benchmark["benchmark_rows"])],
    ]
    save_csv(
        results_dirs["summary"] / "phase10_overview.csv",
        ["component", "dense_sparse_gap", "reparam_fd_gap", "objective_drop_or_best", "runtime_total"],
        overview_rows,
    )
    payload = {
        "backward": backward["summary"],
        "optimizer": optimizer["summary"],
        "benchmark": {
            "n_rows": len(benchmark["overview_rows"]),
            "n_ablations": len(benchmark["ablation_rows"]),
        },
        "config": {
            "n_steps": int(n_steps),
            "T_min": float(T_min),
            "seed": int(seed),
        },
    }
    save_json(results_dirs["summary"] / "phase10_suite_summary.json", payload)
    _write_interpretation_summary(backward, optimizer, benchmark, results_dirs["summary"] / "phase10_interpretation_summary.md")
    return {"backward": backward, "optimizer": optimizer, "benchmark": benchmark, "summary": payload}


def main() -> None:
    q0, T0 = representative_case()
    tau0 = T_to_tau(T0, T_min=0.2)
    run_phase10_framework_suite(q0=q0, tau0=tau0, n_steps=30, T_min=0.2)


if __name__ == "__main__":
    main()
