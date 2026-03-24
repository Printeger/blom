"""
Symbolic validation for the warm-up BLOM analytic model (s=2, k=2).
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np
import sympy as sp

from phase_4.utils.io_utils import ensure_dir, save_json, save_text
from phase_4.utils.plotting_utils import save_line_plot, save_scatter_with_diagonal


DEFAULT_RESULTS_DIR = Path("phase_4/results/s2_sympy")


def _symbols() -> dict[str, sp.Symbol]:
    q_L, q_R, q_im1, q_i, q_ip1 = sp.symbols("q_L q_R q_im1 q_i q_ip1", real=True)
    h, h_m, h_p = sp.symbols("h h_m h_p", positive=True, real=True)
    v_R, v_L, v_i = sp.symbols("v_R v_L v_i", real=True)
    tau = sp.symbols("tau", real=True)
    return {
        "q_L": q_L,
        "q_R": q_R,
        "q_im1": q_im1,
        "q_i": q_i,
        "q_ip1": q_ip1,
        "h": h,
        "h_m": h_m,
        "h_p": h_p,
        "v_R": v_R,
        "v_L": v_L,
        "v_i": v_i,
        "tau": tau,
    }


@lru_cache(maxsize=1)
def build_left_natural_cubic_symbolic() -> dict[str, Any]:
    sym = _symbols()
    a0, a1, a2, a3 = sp.symbols("a0 a1 a2 a3", real=True)
    p = a0 + a1 * sym["tau"] + a2 * sym["tau"] ** 2 + a3 * sym["tau"] ** 3
    equations = [
        sp.Eq(p.subs(sym["tau"], 0), sym["q_L"]),
        sp.Eq(p.subs(sym["tau"], 1), sym["q_R"]),
        sp.Eq(sp.diff(p, sym["tau"]).subs(sym["tau"], 1) / sym["h"], sym["v_R"]),
        sp.Eq(sp.diff(p, sym["tau"], 2).subs(sym["tau"], 0) / sym["h"] ** 2, 0),
    ]
    solution = sp.solve(equations, (a0, a1, a2, a3), dict=True)[0]
    p_exact = sp.expand(p.subs(solution))
    return {
        "polynomial": p_exact,
        "coefficients": {str(key): sp.simplify(value) for key, value in solution.items()},
    }


@lru_cache(maxsize=1)
def build_right_natural_cubic_symbolic() -> dict[str, Any]:
    sym = _symbols()
    a0, a1, a2, a3 = sp.symbols("b0 b1 b2 b3", real=True)
    p = a0 + a1 * sym["tau"] + a2 * sym["tau"] ** 2 + a3 * sym["tau"] ** 3
    equations = [
        sp.Eq(p.subs(sym["tau"], 0), sym["q_L"]),
        sp.Eq(p.subs(sym["tau"], 1), sym["q_R"]),
        sp.Eq(sp.diff(p, sym["tau"]).subs(sym["tau"], 0) / sym["h"], sym["v_L"]),
        sp.Eq(sp.diff(p, sym["tau"], 2).subs(sym["tau"], 1) / sym["h"] ** 2, 0),
    ]
    solution = sp.solve(equations, (a0, a1, a2, a3), dict=True)[0]
    p_exact = sp.expand(p.subs(solution))
    return {
        "polynomial": p_exact,
        "coefficients": {str(key): sp.simplify(value) for key, value in solution.items()},
    }


@lru_cache(maxsize=1)
def derive_s2_local_energy_symbolic() -> dict[str, Any]:
    sym = _symbols()
    left = build_left_natural_cubic_symbolic()["polynomial"]
    right = build_right_natural_cubic_symbolic()["polynomial"]

    left_second = sp.diff(left, sym["tau"], 2) / sym["h"] ** 2
    right_second = sp.diff(right, sym["tau"], 2) / sym["h"] ** 2
    J_left = sp.simplify(sp.integrate(left_second**2 * sym["h"], (sym["tau"], 0, 1)))
    J_right = sp.simplify(sp.integrate(right_second**2 * sym["h"], (sym["tau"], 0, 1)))

    s_minus = (sym["q_i"] - sym["q_im1"]) / sym["h_m"]
    s_plus = (sym["q_ip1"] - sym["q_i"]) / sym["h_p"]
    J_total = sp.expand(
        J_left.subs(
            {
                sym["q_L"]: sym["q_im1"],
                sym["q_R"]: sym["q_i"],
                sym["h"]: sym["h_m"],
                sym["v_R"]: sym["v_i"],
            }
        )
        + J_right.subs(
            {
                sym["q_L"]: sym["q_i"],
                sym["q_R"]: sym["q_ip1"],
                sym["h"]: sym["h_p"],
                sym["v_L"]: sym["v_i"],
            }
        )
    )
    quadratic_coeff = sp.expand(sp.diff(J_total, sym["v_i"], 2) / 2)
    return {
        "left_energy": sp.simplify(J_left),
        "right_energy": sp.simplify(J_right),
        "total_energy": sp.simplify(J_total),
        "s_minus": sp.simplify(s_minus),
        "s_plus": sp.simplify(s_plus),
        "quadratic_coeff": sp.simplify(quadratic_coeff),
    }


@lru_cache(maxsize=1)
def derive_s2_exact_velocity_symbolic() -> dict[str, Any]:
    sym = _symbols()
    energy = derive_s2_local_energy_symbolic()["total_energy"]
    gradient = sp.simplify(sp.diff(energy, sym["v_i"]))
    solution = sp.solve(sp.Eq(gradient, 0), sym["v_i"])[0]
    return {
        "gradient": gradient,
        "v_exact": sp.simplify(solution),
    }


@lru_cache(maxsize=1)
def verify_catmull_equivalence_symbolic() -> dict[str, Any]:
    sym = _symbols()
    derived = derive_s2_exact_velocity_symbolic()["v_exact"]
    s_minus = (sym["q_i"] - sym["q_im1"]) / sym["h_m"]
    s_plus = (sym["q_ip1"] - sym["q_i"]) / sym["h_p"]
    v_catmull = (sym["h_p"] * s_minus + sym["h_m"] * s_plus) / (sym["h_m"] + sym["h_p"])
    diff = sp.simplify(sp.factor(derived - v_catmull))
    return {
        "v_exact": derived,
        "v_catmull": sp.simplify(v_catmull),
        "difference": diff,
    }


def _energy_numeric(v_values: np.ndarray, q_im1: float, q_i: float, q_ip1: float, h_m: float, h_p: float) -> np.ndarray:
    energy = derive_s2_local_energy_symbolic()["total_energy"]
    sym = _symbols()
    fn = sp.lambdify(
        (sym["v_i"], sym["q_im1"], sym["q_i"], sym["q_ip1"], sym["h_m"], sym["h_p"]),
        energy,
        "numpy",
    )
    return np.asarray(fn(v_values, q_im1, q_i, q_ip1, h_m, h_p), dtype=float)


def _exact_velocity_numeric(q_im1: float, q_i: float, q_ip1: float, h_m: float, h_p: float) -> float:
    expr = derive_s2_exact_velocity_symbolic()["v_exact"]
    sym = _symbols()
    fn = sp.lambdify((sym["q_im1"], sym["q_i"], sym["q_ip1"], sym["h_m"], sym["h_p"]), expr, "numpy")
    return float(fn(q_im1, q_i, q_ip1, h_m, h_p))


def _catmull_velocity_numeric(q_im1: float, q_i: float, q_ip1: float, h_m: float, h_p: float) -> float:
    s_minus = (q_i - q_im1) / h_m
    s_plus = (q_ip1 - q_i) / h_p
    return float((h_p * s_minus + h_m * s_plus) / (h_m + h_p))


def main(results_dir: str | Path = DEFAULT_RESULTS_DIR, seed: int = 42) -> dict[str, Any]:
    results_dir = ensure_dir(results_dir)
    left = build_left_natural_cubic_symbolic()
    right = build_right_natural_cubic_symbolic()
    energy = derive_s2_local_energy_symbolic()
    exact_velocity = derive_s2_exact_velocity_symbolic()
    equivalence = verify_catmull_equivalence_symbolic()

    save_text(results_dir / "s2_left_segment_coeffs.txt", sp.pretty(left["coefficients"]))
    save_text(results_dir / "s2_right_segment_coeffs.txt", sp.pretty(right["coefficients"]))
    save_text(results_dir / "s2_energy_expression.txt", sp.pretty(energy["total_energy"]))
    save_text(results_dir / "s2_optimal_velocity_expression.txt", sp.pretty(exact_velocity["v_exact"]))
    save_text(results_dir / "s2_catmull_equivalence.txt", sp.pretty(equivalence["difference"]))

    q_im1, q_i, q_ip1 = 0.0, 1.3, 0.4
    h_m, h_p = 0.8, 1.6
    v_exact = _exact_velocity_numeric(q_im1, q_i, q_ip1, h_m, h_p)
    v_values = np.linspace(v_exact - 4.0, v_exact + 4.0, 400)
    energy_values = _energy_numeric(v_values, q_im1, q_i, q_ip1, h_m, h_p)
    save_line_plot(
        v_values,
        energy_values,
        results_dir / "s2_energy_vs_v.png",
        title=f"s=2 Local Energy vs v_i, h-= {h_m:.2f}, h+= {h_p:.2f}",
        xlabel="v_i",
        ylabel="J_i^(2)(v_i)",
    )

    rng = np.random.default_rng(seed)
    exact_samples = []
    catmull_samples = []
    for _ in range(128):
        q_vals = rng.normal(size=3)
        h_vals = rng.uniform(0.4, 1.8, size=2)
        exact_samples.append(_exact_velocity_numeric(q_vals[0], q_vals[1], q_vals[2], h_vals[0], h_vals[1]))
        catmull_samples.append(_catmull_velocity_numeric(q_vals[0], q_vals[1], q_vals[2], h_vals[0], h_vals[1]))
    exact_arr = np.asarray(exact_samples, dtype=float)
    catmull_arr = np.asarray(catmull_samples, dtype=float)
    save_scatter_with_diagonal(
        exact_arr,
        catmull_arr,
        results_dir / "s2_catmull_match_demo.png",
        title="s=2 Exact vs Catmull Velocity",
        xlabel="v_exact",
        ylabel="v_catmull",
    )

    metrics = {
        "quadratic_coeff": str(energy["quadratic_coeff"]),
        "quadratic_coeff_positive_demo": bool(float(sp.N(energy["quadratic_coeff"].subs({_symbols()["h_m"]: 1.2, _symbols()["h_p"]: 0.7}))) > 0.0),
        "symbolic_difference": str(equivalence["difference"]),
        "random_max_abs_error": float(np.max(np.abs(exact_arr - catmull_arr))),
        "random_mean_abs_error": float(np.mean(np.abs(exact_arr - catmull_arr))),
        "v_exact_demo": v_exact,
        "v_catmull_demo": _catmull_velocity_numeric(q_im1, q_i, q_ip1, h_m, h_p),
    }
    save_json(results_dir / "s2_metrics.json", metrics)
    return {
        "inputs": {
            "demo_q": [q_im1, q_i, q_ip1],
            "demo_h": [h_m, h_p],
        },
        "solution": {
            "v_exact_expression": str(exact_velocity["v_exact"]),
            "v_catmull_expression": str(equivalence["v_catmull"]),
        },
        "metrics": metrics,
        "fig_paths": {
            "energy_vs_v": str(results_dir / "s2_energy_vs_v.png"),
            "catmull_match_demo": str(results_dir / "s2_catmull_match_demo.png"),
        },
    }


if __name__ == "__main__":
    summary = main()
    print("Phase 4 s=2 symbolic validation")
    print(f"symbolic difference: {summary['metrics']['symbolic_difference']}")
    print(f"random max abs error: {summary['metrics']['random_max_abs_error']:.3e}")
    print(f"results dir: {DEFAULT_RESULTS_DIR}")

