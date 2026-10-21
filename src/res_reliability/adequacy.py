"""Analytical adequacy indices from a capacity outage probability table (COPT).

For independent two-state conventional units and a load that is known hour by hour, the
expected loss-of-load hours and unserved energy follow exactly from the distribution of
available conventional capacity A:

    LOLE = sum_h P(A < L_h - R_h),     EENS = sum_h E[max(L_h - R_h - A, 0)]

R_h may be deterministic (historical weather) or multi-state with hour-specific state
probabilities (Markov model). Expected LOLE/EENS do not depend on chronology; sequential
Monte Carlo (``montecarlo``) is needed for frequency, duration and the spread of annual values.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class Copt:
    """Distribution of available capacity on a 1 MW grid (capacities are rounded)."""

    p: np.ndarray  # p[c] = P(available capacity == c MW)

    @classmethod
    def from_units(cls, capacities, forced_outage_rates) -> "Copt":
        p = np.array([1.0])
        for cap, q in zip(capacities, forced_outage_rates):
            c = int(round(cap))
            new = np.zeros(len(p) + c)
            new[: len(p)] += p * q
            new[c:] += p * (1.0 - q)
            p = new
        return cls(p)

    @property
    def installed(self) -> int:
        return len(self.p) - 1

    def _cum(self):
        if not hasattr(self, "_F"):
            c = np.arange(len(self.p))
            self._F = np.cumsum(self.p)
            self._G = np.cumsum(self.p * c)
        return self._F, self._G

    def shortfall(self, x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """P(A < x) and E[max(x - A, 0)] element-wise for net load x (MW)."""
        F, G = self._cum()
        x = np.asarray(x, dtype=float)
        k = np.ceil(x).astype(np.int64) - 1  # A < x  <=>  A <= ceil(x) - 1 on the integer grid
        k_clip = np.clip(k, 0, self.installed)
        prob = np.where(k < 0, 0.0, F[k_clip])
        expected = np.where(k < 0, 0.0, x * F[k_clip] - G[k_clip])
        above = k > self.installed  # net load above installed capacity: certain shortfall
        prob = np.where(above, 1.0, prob)
        expected = np.where(above, x - (self.p * np.arange(len(self.p))).sum(), expected)
        return prob, expected


def indices_deterministic(copt: Copt, load: np.ndarray, res: np.ndarray | float = 0.0) -> dict:
    """LOLE (h) and EENS (MWh) for one year with a known RES output per hour."""
    prob, exp = copt.shortfall(np.asarray(load) - res)
    return {"lole": float(prob.sum()), "eens": float(exp.sum())}


def indices_multistate(copt: Copt, load: np.ndarray, res_levels: np.ndarray, res_probs: np.ndarray) -> dict:
    """LOLE/EENS with a multi-state RES: ``res_levels`` and ``res_probs`` are (hours, states)."""
    x = np.asarray(load)[:, None] - res_levels
    prob, exp = copt.shortfall(x)
    return {"lole": float((prob * res_probs).sum()), "eens": float((exp * res_probs).sum())}
