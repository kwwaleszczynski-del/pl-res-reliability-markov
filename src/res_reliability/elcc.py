"""Effective load carrying capability (ELCC) by LOLE equality (Garver 1966; Keane et al. 2011).

ELCC is the constant load increase that the system can carry after adding the RES plant at
the same expected LOLE as before:  LOLE(conv + RES, L + ELCC) = LOLE(conv, L).

The analytical route represents RES output in every hour of every weather year as a discrete
distribution (levels, probabilities): a single level for a historical weather year, the
stationary distribution of the month for a Markov model. Structural reliability is either a
derating (many independent units) or a plant-level two-point mixture (thesis variant).
The Monte Carlo route uses stored margins of common-random-number runs.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from res_reliability import markov as mk
from res_reliability.adequacy import Copt


@dataclass
class ResDistribution:
    levels: np.ndarray  # (n_years, n_hours, n_states) MW
    probs: np.ndarray  # same shape, sums to 1 over the last axis
    capacity_mw: float

    @property
    def mean_mw(self) -> float:
        return float((self.levels * self.probs).sum(axis=-1).mean())

    def __add__(self, other: "ResDistribution") -> "ResDistribution":
        """Independent sum of two RES distributions (outer product of states)."""
        y = max(self.levels.shape[0], other.levels.shape[0])
        a_l = np.broadcast_to(self.levels, (y,) + self.levels.shape[1:])
        a_p = np.broadcast_to(self.probs, (y,) + self.probs.shape[1:])
        b_l = np.broadcast_to(other.levels, (y,) + other.levels.shape[1:])
        b_p = np.broadcast_to(other.probs, (y,) + other.probs.shape[1:])
        lv = (a_l[..., :, None] + b_l[..., None, :]).reshape(y, a_l.shape[1], -1)
        pr = (a_p[..., :, None] * b_p[..., None, :]).reshape(y, a_l.shape[1], -1)
        return ResDistribution(lv, pr, self.capacity_mw + other.capacity_mw)


def apply_structural(dist: ResDistribution, forced_outage_rate: float, mode: str = "units") -> ResDistribution:
    if forced_outage_rate <= 0 or mode == "none":
        return dist
    if mode == "units":  # many independent units: available share is ~deterministic
        return ResDistribution(dist.levels * (1 - forced_outage_rate), dist.probs, dist.capacity_mw)
    if mode == "plant_bernoulli":
        lv = np.concatenate([dist.levels, np.zeros_like(dist.levels)], axis=-1)
        pr = np.concatenate([dist.probs * (1 - forced_outage_rate), dist.probs * forced_outage_rate], axis=-1)
        return ResDistribution(lv, pr, dist.capacity_mw)
    raise ValueError(mode)


def historical(series: np.ndarray, capacity_mw: float, n_hours: int | None = None) -> ResDistribution:
    s = np.asarray(series, dtype=float)[:, :n_hours]
    return ResDistribution(capacity_mw * s[..., None], np.ones(s.shape + (1,)), capacity_mw)


def from_markov(
    model: mk.MarkovModel,
    months: np.ndarray,
    capacity_mw: float,
    envelope: np.ndarray | None = None,
    n_years: int = 1,
    n_hours: int | None = None,
) -> ResDistribution:
    """Marginal state distribution per hour: pi of the month (seasonal) or pi (homogeneous)."""
    months = np.asarray(months)[:n_hours]
    h = len(months)
    pis = np.stack([model.stationary(m) for m in range(1, 13)])  # (12, S)
    lvs = np.stack([model.level(m) for m in range(1, 13)])
    probs = pis[months - 1][None]  # (1, H, S)
    levels = capacity_mw * lvs[months - 1][None]
    if envelope is not None:
        env = np.asarray(envelope, dtype=float)[:, :h]
        levels = levels * env[..., None]
        probs = np.broadcast_to(probs, levels.shape)
    elif n_years > 1:
        levels = np.broadcast_to(levels, (n_years,) + levels.shape[1:])
        probs = np.broadcast_to(probs, levels.shape)
    return ResDistribution(np.clip(levels, 0, capacity_mw), np.asarray(probs), capacity_mw)


def lole(copt: Copt, load: np.ndarray, res: ResDistribution | None = None, shift: float = 0.0) -> float:
    """Expected LOLE (h/year) averaged over the weather/load years."""
    return indices(copt, load, res, shift)["lole"]


def indices(copt: Copt, load: np.ndarray, res: ResDistribution | None = None, shift: float = 0.0) -> dict:
    load = np.atleast_2d(np.asarray(load, dtype=float))
    if res is None:
        prob, exp = copt.shortfall(load + shift)
        return {"lole": float(prob.sum(axis=1).mean()), "eens": float(exp.sum(axis=1).mean())}
    h = res.levels.shape[1]
    y = max(load.shape[0], res.levels.shape[0])
    x = np.broadcast_to(load[:, :h], (y, h))[..., None] + shift - res.levels
    prob, exp = copt.shortfall(x)
    return {"lole": float((prob * res.probs).sum(axis=(1, 2)).mean()),
            "eens": float((exp * res.probs).sum(axis=(1, 2)).mean())}


def elcc(copt: Copt, load: np.ndarray, res: ResDistribution, target: float | None = None, xtol: float = 0.05) -> float:
    """ELCC in MW (load increase restoring the base LOLE)."""
    base = lole(copt, load) if target is None else target
    f = lambda s: lole(copt, load, res, s) - base  # noqa: E731
    lo, hi = 0.0, res.capacity_mw
    if f(lo) >= 0:
        return 0.0
    if f(hi) <= 0:
        return hi
    return float(brentq(f, lo, hi, xtol=xtol))


def mc_lole_at_shift(margin: np.ndarray, shift: float) -> float:
    return float((margin < shift).sum(axis=1).mean())


def mc_elcc(base_margin: np.ndarray, res_margin: np.ndarray, capacity_mw: float, tol: float = 0.1) -> float:
    """ELCC from Monte Carlo margins (available - load) of runs with common random numbers."""
    target = mc_lole_at_shift(base_margin, 0.0)
    lo, hi = 0.0, capacity_mw
    if mc_lole_at_shift(res_margin, lo) >= target:
        return 0.0
    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        if mc_lole_at_shift(res_margin, mid) < target:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)
