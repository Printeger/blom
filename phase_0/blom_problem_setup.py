"""
Canonical Phase 0 problem setup for the scalar BLOM research baseline.

This module freezes the shared input object used by all later Phase 0 code:
global MINCO, local BLOM windows, trajectory evaluation, and validation.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
import json
from typing import Any

import numpy as np


DEFAULT_SEED = 42


@dataclass
class BLOMProblemSetup:
    """
    Canonical Phase 0 problem data.

    Shapes
    ------
    q : (M + 1,)
        Scalar waypoint sequence q_0, ..., q_M.
    T : (M,)
        Positive segment durations T_1, ..., T_M.
    """

    M: int
    q: np.ndarray
    T: np.ndarray
    s: int = 4
    k: int = 2
    d_i: int = 1
    dim: int = 1
    boundary_condition: str = "natural"
    seed: int | None = None

    def __post_init__(self) -> None:
        self.q = np.asarray(self.q, dtype=float).reshape(-1)
        self.T = np.asarray(self.T, dtype=float).reshape(-1)

    @classmethod
    def make_random(
        cls,
        M: int,
        *,
        seed: int = DEFAULT_SEED,
        q_mode: str = "random_walk",
        q_scale: float = 1.0,
        T_low: float = 0.5,
        T_high: float = 1.5,
    ) -> "BLOMProblemSetup":
        """Generate a deterministic Phase 0 random instance."""
        if M < 2:
            raise ValueError(f"Phase 0 requires M >= 2, got {M}.")
        if not (0.0 < T_low < T_high):
            raise ValueError(
                f"Require 0 < T_low < T_high, got T_low={T_low}, T_high={T_high}."
            )

        rng = np.random.default_rng(seed)
        if q_mode == "random_walk":
            steps = rng.normal(loc=0.0, scale=q_scale, size=M)
            q = np.concatenate(([0.0], np.cumsum(steps)))
        elif q_mode == "uniform":
            q = rng.uniform(low=-q_scale, high=q_scale, size=M + 1)
        else:
            raise ValueError(f"Unsupported q_mode: {q_mode}")

        T = rng.uniform(low=T_low, high=T_high, size=M)
        setup = cls(M=M, q=q, T=T, seed=seed)
        setup.validate()
        return setup

    def validate(self) -> None:
        """Validate the canonical Phase 0 assumptions."""
        if self.M < 2:
            raise ValueError(f"Phase 0 requires M >= 2, got {self.M}.")
        if self.s != 4:
            raise ValueError(f"Phase 0 requires s=4, got {self.s}.")
        if self.k != 2:
            raise ValueError(f"Phase 0 requires k=2, got {self.k}.")
        if self.d_i != 1:
            raise ValueError(f"Phase 0 requires d_i=1, got {self.d_i}.")
        if self.dim != 1:
            raise ValueError(f"Phase 0 requires dim=1, got {self.dim}.")
        if self.boundary_condition != "natural":
            raise ValueError(
                "Phase 0 only supports natural boundary conditions. "
                f"Got {self.boundary_condition!r}."
            )
        if self.q.shape != (self.M + 1,):
            raise ValueError(
                f"q must have shape ({self.M + 1},), got {self.q.shape}."
            )
        if self.T.shape != (self.M,):
            raise ValueError(f"T must have shape ({self.M},), got {self.T.shape}.")
        if not np.all(np.isfinite(self.q)):
            raise ValueError("q contains non-finite values.")
        if not np.all(np.isfinite(self.T)):
            raise ValueError("T contains non-finite values.")
        if not np.all(self.T > 0.0):
            raise ValueError("All segment durations T_i must be positive.")

    def theta(self) -> np.ndarray:
        """
        Return the canonical parameter vector
            theta = (q_1, ..., q_{M-1}, T_1, ..., T_M).
        """
        return np.concatenate((self.q[1:-1], self.T), axis=0)

    def window(self, i: int) -> tuple[int, ...]:
        """
        Return the 0-based active segment indices for W(i, 2).

        The input i is a 0-based segment index in [0, M-1].
        """
        if not 0 <= i < self.M:
            raise IndexError(f"center segment must be in [0, {self.M - 1}], got {i}.")
        left = max(0, i - self.k // 2)
        right = min(self.M - 1, i + self.k // 2)
        return tuple(range(left, right + 1))

    @property
    def total_time(self) -> float:
        """Return the total trajectory duration."""
        return float(np.sum(self.T))

    def to_dict(self) -> dict[str, Any]:
        """Serialize the setup into a JSON-friendly dictionary."""
        self.validate()
        return {
            "M": self.M,
            "q": self.q.tolist(),
            "T": self.T.tolist(),
            "s": self.s,
            "k": self.k,
            "d_i": self.d_i,
            "dim": self.dim,
            "boundary_condition": self.boundary_condition,
            "seed": self.seed,
            "theta": self.theta().tolist(),
            "total_time": self.total_time,
            "window_map": {str(i): list(self.window(i)) for i in range(self.M)},
        }

    @classmethod
    def from_dict(cls, data: dict[str, Any]) -> "BLOMProblemSetup":
        """Deserialize Phase 0 setup data from a dictionary."""
        d_i = data.get("d_i", data.get("d_intermediate", 1))
        setup = cls(
            M=int(data["M"]),
            q=np.asarray(data["q"], dtype=float).reshape(-1),
            T=np.asarray(data["T"], dtype=float).reshape(-1),
            s=int(data.get("s", 4)),
            k=int(data.get("k", 2)),
            d_i=int(d_i),
            dim=int(data.get("dim", 1)),
            boundary_condition=str(data.get("boundary_condition", "natural")),
            seed=None if data.get("seed") is None else int(data["seed"]),
        )
        setup.validate()
        return setup

    def save_json(self, path: str | Path) -> None:
        """Save the setup as UTF-8 JSON."""
        path = Path(path)
        path.write_text(
            json.dumps(self.to_dict(), indent=2, ensure_ascii=False),
            encoding="utf-8",
        )

    @classmethod
    def load_json(cls, path: str | Path) -> "BLOMProblemSetup":
        """Load a setup from UTF-8 JSON."""
        return cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))


def generate_random_problem(
    M: int,
    *,
    seed: int = DEFAULT_SEED,
    q_mode: str = "random_walk",
    q_scale: float = 1.0,
    T_low: float = 0.5,
    T_high: float = 1.5,
) -> BLOMProblemSetup:
    """Compatibility wrapper around BLOMProblemSetup.make_random()."""
    return BLOMProblemSetup.make_random(
        M,
        seed=seed,
        q_mode=q_mode,
        q_scale=q_scale,
        T_low=T_low,
        T_high=T_high,
    )


def print_problem_report(problem: BLOMProblemSetup) -> None:
    """Print a compact human-readable report for the canonical setup."""
    problem.validate()
    print("=" * 72)
    print("BLOM Phase 0 Problem Setup")
    print("=" * 72)
    print(f"Segments M              : {problem.M}")
    print(f"Output dimension        : {problem.dim}")
    print(f"Control order s         : {problem.s}")
    print(f"Window width k          : {problem.k}")
    print(f"Intermediate order d_i  : {problem.d_i}")
    print(f"Boundary condition      : {problem.boundary_condition}")
    print(f"Seed                    : {problem.seed}")
    print(f"Total time              : {problem.total_time:.6f}")
    print("-" * 72)
    print("Waypoints q:")
    print(problem.q)
    print("-" * 72)
    print("Durations T:")
    print(problem.T)
    print("-" * 72)
    print("theta = (q_1,...,q_{M-1}, T_1,...,T_M):")
    print(problem.theta())
    print("-" * 72)
    print("Window map in 0-based segment indices:")
    for i in range(problem.M):
        print(f"  segment {i}: {list(problem.window(i))}")
    print("=" * 72)


def main() -> None:
    """Generate the canonical sample JSON used by the repository."""
    problem = BLOMProblemSetup.make_random(M=6, seed=DEFAULT_SEED, q_scale=0.8)
    print_problem_report(problem)
    out_path = Path(__file__).with_name("blom_phase0_problem.json")
    problem.save_json(out_path)
    print(f"Saved JSON to: {out_path}")


if __name__ == "__main__":
    main()
