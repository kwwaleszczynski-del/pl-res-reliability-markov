"""Multi-state discrete-time Markov models of primary-energy availability.

A normalised output series x_t in [0, 1] (wind capacity factor, or the PV clear-sky index)
is discretised into N states. As in the thesis, the state thresholds are the mid-points
between equally spaced levels k/(N-1); the level attached to a state is either the bin
centre (thesis) or the conditional mean of x in the bin (default: preserves the energy).

The transition matrix is the maximum-likelihood estimate P_ij = n_ij / n_i counted over
consecutive valid hours. For PV only daylight hours are valid and the chain moves from the
last daylight hour of a day to the first daylight hour of the next one, which is also how it
is simulated. The model can be homogeneous (one matrix) or month-dependent (12 matrices).

The continuous-time intensities used in the thesis follow from lambda_ij = n_ij / T_i.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


def thresholds(n_states: int) -> np.ndarray:
    """Mid-points between the levels 0, 1/(N-1), ..., 1."""
    return (np.arange(n_states - 1) + 0.5) / (n_states - 1)


def discretize(x: np.ndarray, n_states: int) -> np.ndarray:
    """State index 0..N-1 for each value; NaN -> -1."""
    x = np.asarray(x, dtype=float)
    s = np.searchsorted(thresholds(n_states), x, side="right").astype(np.int64)
    s[np.isnan(x)] = -1
    return s


def state_levels(x: np.ndarray, states: np.ndarray, n_states: int, kind: str = "mean") -> np.ndarray:
    centres = np.arange(n_states) / (n_states - 1)
    if kind == "centre":
        return centres
    levels = centres.copy()
    for k in range(n_states):
        sel = states == k
        if sel.any():
            levels[k] = float(np.mean(x[sel]))
    return levels


def transition_counts(states: np.ndarray, n_states: int, valid: np.ndarray | None = None) -> np.ndarray:
    """Counts n_ij over consecutive valid observations (invalid ones are skipped, not bridged
    unless they are excluded by ``valid``, in which case the chain jumps over them)."""
    s = np.asarray(states)
    if valid is not None:
        s = s[np.asarray(valid, dtype=bool)]
    a, b = s[:-1], s[1:]
    ok = (a >= 0) & (b >= 0)
    counts = np.zeros((n_states, n_states))
    np.add.at(counts, (a[ok], b[ok]), 1.0)
    return counts


def normalise_rows(counts: np.ndarray, fallback: np.ndarray | None = None) -> np.ndarray:
    """Row-normalise counts; rows without data get ``fallback`` (default: stay put)."""
    counts = np.asarray(counts, dtype=float)
    n = counts.shape[0]
    tot = counts.sum(axis=1, keepdims=True)
    P = np.divide(counts, tot, out=np.zeros_like(counts), where=tot > 0)
    for i in np.flatnonzero(tot[:, 0] == 0):
        P[i] = fallback if fallback is not None else np.eye(n)[i]
    return P


def stationary(P: np.ndarray) -> np.ndarray:
    """Stationary distribution pi = pi P, sum(pi) = 1 (unique for an irreducible chain)."""
    P = np.asarray(P, dtype=float)
    n = P.shape[0]
    A = np.vstack([P.T - np.eye(n), np.ones(n)])
    b = np.zeros(n + 1)
    b[-1] = 1.0
    pi, *_ = np.linalg.lstsq(A, b, rcond=None)
    pi = np.clip(pi, 0.0, None)
    return pi / pi.sum()


def intensities(states: np.ndarray, n_states: int, dt_h: float = 1.0) -> np.ndarray:
    """Continuous-time intensity matrix Lambda (1/h): lambda_ij = n_ij / T_i, diagonal = -row sum."""
    counts = transition_counts(states, n_states)
    np.fill_diagonal(counts, 0.0)
    occupancy = np.bincount(np.asarray(states)[np.asarray(states) >= 0], minlength=n_states) * dt_h
    lam = np.divide(counts, occupancy[:, None], out=np.zeros_like(counts), where=occupancy[:, None] > 0)
    np.fill_diagonal(lam, -lam.sum(axis=1))
    return lam


@dataclass
class MarkovModel:
    """Fitted chain: ``P`` is (S, S) or (12, S, S); ``levels`` is (S,) or (12, S)."""

    P: np.ndarray
    levels: np.ndarray
    n_states: int

    @property
    def seasonal(self) -> bool:
        return self.P.ndim == 3

    def matrix(self, month: int) -> np.ndarray:
        return self.P[month - 1] if self.seasonal else self.P

    def level(self, month: int) -> np.ndarray:
        return self.levels[month - 1] if self.levels.ndim == 2 else self.levels

    def stationary(self, month: int | None = None) -> np.ndarray:
        if self.seasonal:
            return stationary(self.P[(month or 1) - 1])
        return stationary(self.P)

    def mean(self, month: int | None = None) -> float:
        if month is None and self.seasonal:
            return float(np.mean([self.mean(m) for m in range(1, 13)]))
        return float(self.stationary(month) @ self.level(month or 1))


def fit(
    x: np.ndarray,
    n_states: int = 5,
    months: np.ndarray | None = None,
    valid: np.ndarray | None = None,
    level_kind: str = "mean",
) -> MarkovModel:
    """Fit a homogeneous chain, or a month-dependent one when ``months`` (1..12) is given.

    For the seasonal chain a transition t -> t+1 belongs to the month of t.
    """
    x = np.asarray(x, dtype=float)
    valid = np.ones(len(x), dtype=bool) if valid is None else np.asarray(valid, dtype=bool)
    xv = np.where(valid, x, np.nan)
    states = discretize(xv, n_states)
    keep = valid & (states >= 0)
    if months is None:
        counts = transition_counts(states, n_states, keep)
        freq = np.bincount(states[keep], minlength=n_states) / keep.sum()
        P = normalise_rows(counts, fallback=freq)
        levels = state_levels(x[keep], states[keep], n_states, level_kind)
        return MarkovModel(P, levels, n_states)
    months = np.asarray(months)
    s_keep, m_keep, x_keep = states[keep], months[keep], x[keep]
    Ps, Ls = [], []
    for m in range(1, 13):
        sel = m_keep[:-1] == m
        counts = np.zeros((n_states, n_states))
        np.add.at(counts, (s_keep[:-1][sel], s_keep[1:][sel]), 1.0)
        in_m = m_keep == m
        freq = np.bincount(s_keep[in_m], minlength=n_states) / max(in_m.sum(), 1)
        Ps.append(normalise_rows(counts, fallback=freq))
        Ls.append(state_levels(x_keep[in_m], s_keep[in_m], n_states, level_kind))
    return MarkovModel(np.stack(Ps), np.stack(Ls), n_states)


def simulate_states(
    model: MarkovModel,
    uniforms: np.ndarray,
    months: np.ndarray,
    active: np.ndarray | None = None,
) -> np.ndarray:
    """Simulate many independent chains at once.

    ``uniforms`` is (n_sims, n_hours + 1): column 0 draws the initial state from the
    stationary distribution of the first month, column h+1 drives the step after hour h.
    The chain only moves on ``active`` hours (PV: daylight); on inactive hours it keeps its
    state, so an overnight gap is a single transition, as in the fitted data.
    Returns (n_sims, n_hours) int8 states.
    """
    u = np.asarray(uniforms)
    n_sims, n_cols = u.shape
    n_hours = n_cols - 1
    months = np.asarray(months)
    active = np.ones(n_hours, dtype=bool) if active is None else np.asarray(active, dtype=bool)
    if model.seasonal:
        cum = np.cumsum(model.P, axis=2)
        cum[:, :, -1] = 1.0
    else:
        cum = np.cumsum(model.P, axis=1)[None].repeat(12, axis=0)
        cum[:, :, -1] = 1.0
    pi0 = np.cumsum(model.stationary(int(months[0])))
    state = np.searchsorted(pi0, u[:, 0], side="right").clip(max=model.n_states - 1)
    out = np.empty((n_sims, n_hours), dtype=np.int8)
    last_active = None
    for h in range(n_hours):
        if active[h]:
            if last_active is not None:
                c = cum[months[last_active] - 1][state]
                state = (u[:, h + 1, None] > c).sum(axis=1).clip(max=model.n_states - 1)
            last_active = h
        out[:, h] = state
    return out
