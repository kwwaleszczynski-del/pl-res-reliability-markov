"""Uncertainty and goodness-of-fit tools for the weather models and the storage results.

* **Bootstrap over weather years.** Ten weather years are a small sample, so every quantity that depends on them
  gets a percentile interval from resampling the years with replacement (and, for Monte Carlo results, the
  simulated years inside each chosen weather year: two-stage bootstrap). The same draws are used for every
  model, so differences between models get intervals too.
* **Curves instead of re-running searches.** LOLE as a function of a load shift (for ELCC) or of the storage
  energy (for the energy needed) is computed once per simulated year on a grid; a bootstrap replicate is then a
  weighted mean of these curves, and its ELCC or energy is read off by linear interpolation.
* **Ten-year block test.** A model is simulated as many independent ten-year blocks (one synthetic version of every
  weather year per block). A statistic of the ERA5 decade is compared with its distribution over the blocks:
  percentile and two-sided p-value.
* **Probability integral transform (PIT).** For a held-out year, the share of the model's predictive sample
  below the observed value; well-calibrated predictions give PIT values spread evenly over [0, 1].
"""

from __future__ import annotations

import numpy as np

from res_reliability import storage as sto
from res_reliability.adequacy import Copt
from res_reliability.elcc import ResDistribution


def year_weights(years: np.ndarray, n_boot: int, rng: np.random.Generator, within: bool = True) -> np.ndarray:
    """(n_boot, n) weights summing to 1: resample weather years with replacement, and with ``within`` also the
    simulated years inside every chosen weather year."""
    years = np.asarray(years)
    groups = [np.flatnonzero(years == y) for y in np.unique(years)]
    w = np.zeros((n_boot, len(years)))
    for b in range(n_boot):
        for k in rng.integers(0, len(groups), len(groups)):
            g = groups[k]
            np.add.at(w[b], rng.choice(g, len(g)) if within else g, 1.0)
    return w / w.sum(axis=1, keepdims=True)


def interval(x: np.ndarray, level: float = 0.9) -> tuple[float, float]:
    a = (1 - level) / 2
    return float(np.nanquantile(x, a)), float(np.nanquantile(x, 1 - a))


def crossing(grid: np.ndarray, curves: np.ndarray, target: np.ndarray, increasing: bool = True,
             log: bool = False) -> np.ndarray:
    """Where each curve (rows of ``curves``, values on ``grid``) first reaches its ``target``.

    Increasing curves (LOLE against a load shift): the first grid point with curve >= target, interpolated
    with the point before it. Decreasing curves (LOLE against storage energy): the first point with
    curve <= target; ``log`` interpolates in log(grid). Returns the grid ends where the target is never or
    always reached.
    """
    grid = np.asarray(grid, dtype=float)
    c = np.atleast_2d(np.asarray(curves, dtype=float))
    t = np.broadcast_to(np.asarray(target, dtype=float), (c.shape[0],))
    hit = c >= t[:, None] if increasing else c <= t[:, None]
    g = np.log(grid) if log else grid
    out = np.empty(c.shape[0])
    for i in range(c.shape[0]):
        idx = np.flatnonzero(hit[i])
        if len(idx) == 0:
            out[i] = g[-1]
        elif idx[0] == 0:
            out[i] = g[0]
        else:
            j = idx[0]
            y0, y1 = c[i, j - 1], c[i, j]
            f = 1.0 if y1 == y0 else (t[i] - y0) / (y1 - y0)
            out[i] = g[j - 1] + f * (g[j] - g[j - 1])
    return np.exp(out) if log else out


def lole_by_year(copt: Copt, load: np.ndarray, res: ResDistribution | None, shift: float = 0.0) -> np.ndarray:
    """Expected LOLE of every weather/load year (exact COPT expectation)."""
    load = np.atleast_2d(np.asarray(load, dtype=float))
    if res is None:
        return copt.shortfall(load + shift)[0].sum(axis=1)
    h = res.levels.shape[1]
    y = max(load.shape[0], res.levels.shape[0])
    x = np.broadcast_to(load[:, :h], (y, h))[..., None] + shift - res.levels
    return (copt.shortfall(x)[0] * res.probs).sum(axis=(1, 2))


def storage_lole_curves(margin: np.ndarray, power: float, energy: float, shifts: np.ndarray, efficiency=1.0,
                        soc0=1.0) -> np.ndarray:
    """(n_shifts, n_sims) LOLE of every simulated year with one storage configuration, for every load shift."""
    margin = np.asarray(margin, dtype=float)
    shifts = np.asarray(shifts, dtype=float)
    out = np.zeros((len(shifts), margin.shape[0]))
    rel = margin.min(axis=1) < shifts.max()
    if rel.any():
        out[:, rel] = sto.dispatch(margin[rel], power, energy, efficiency, soc0, shifts)["lole"]
    return out


def energy_lole_curves(margin: np.ndarray, energies: np.ndarray, efficiency=1.0, soc0=1.0,
                       warmup_h: int = 0) -> np.ndarray:
    """(n_energies, n_sims) LOLE of every simulated year with storage of unlimited power and the given energies."""
    margin = np.asarray(margin, dtype=float)
    out = np.zeros((len(energies), margin.shape[0]))
    rel = margin.min(axis=1) < 0
    if rel.any():
        out[:, rel] = sto.dispatch(margin[rel], np.inf, np.asarray(energies, dtype=float), efficiency, soc0,
                                   warmup_h=warmup_h)["lole"]
    return out


def block_test(data_value: float, block_values: np.ndarray) -> dict:
    """Where the observed decade falls among the model's ten-year blocks."""
    v = np.asarray(block_values, dtype=float)
    v = v[np.isfinite(v)]
    below = float((v < data_value).mean() + 0.5 * (v == data_value).mean())
    lo, hi = interval(v, 0.9)
    return {"data": float(data_value), "model_median": float(np.median(v)), "model_p05": lo, "model_p95": hi,
            "percentile": 100 * below, "p_value": float(min(1.0, 2 * min(below, 1 - below)))}


def pit(observed: float, sample: np.ndarray) -> float:
    """P(model < observed) + half the probability of ties."""
    s = np.asarray(sample, dtype=float)
    return float((s < observed).mean() + 0.5 * (s == observed).mean())
