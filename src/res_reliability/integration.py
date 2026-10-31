"""Integration of the production (weather) model with the structural (failure) model.

The thesis composed the two with an independence formula and checked it with the Kronecker sum of the
intensity matrices,

    Lambda_int = Lambda_prod (x) I_2 + I_N (x) Lambda_struct,

a continuous-time chain on the pairs (production state i, unit up/down). With constant failure and repair
intensities this *is* the independence assumption: its stationary distribution is the product of the two
marginal ones (``kronecker_generator``; checked in the tests). Correlation needs intensities that depend on the
production state. The block generator (``block_generator``) keeps the production transitions of Lambda_prod in
both structural states and lets a unit fail with lambda_i and be repaired with mu_i while the weather is in
state i; with lambda_i = lambda and mu_i = mu it reduces to the Kronecker sum.

lambda_i and mu_i are estimated from turbine event logs as occurrence / exposure rates (the maximum-likelihood
estimate of a piecewise-constant hazard): failures that start in state i over the up-time spent in state i,
repairs that end in state i over the down-time spent in state i. The production state of an hour is the state
of the ERA5-driven model that the adequacy study uses, so the estimates plug straight into it.

Usage:
    exp = state_exposure(windows, outages, hour_ns, states, 5)   # up/down hours, failures, repairs per state
    tab = rate_table(exp)
    Q = block_generator(Lp, tab["lambda"], tab["mu"])
    A = availability_by_state(Lp, tab["lambda"], tab["mu"])     # P(unit up | production state i)
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats

from res_reliability import elcc, markov as mk
from res_reliability.farms import covered_until

HOUR_NS = 3600 * 10**9


# --- generators ---------------------------------------------------------------------------

def struct_generator(lam: float, mu: float) -> np.ndarray:
    """Two-state unit, states ordered (up, down)."""
    return np.array([[-lam, lam], [mu, -mu]], dtype=float)


def kronecker_generator(Lp: np.ndarray, lam: float, mu: float) -> np.ndarray:
    """Kronecker sum Lambda_prod (x) I_2 + I_N (x) Lambda_struct; state index = 2 i + (0 up, 1 down)."""
    Lp = np.asarray(Lp, dtype=float)
    n = Lp.shape[0]
    return np.kron(Lp, np.eye(2)) + np.kron(np.eye(n), struct_generator(lam, mu))


def block_generator(Lp: np.ndarray, lam: np.ndarray, mu: np.ndarray) -> np.ndarray:
    """Integrated generator with failure and repair intensities that depend on the production state."""
    Lp = np.asarray(Lp, dtype=float)
    n = Lp.shape[0]
    lam = np.broadcast_to(np.asarray(lam, dtype=float), (n,))
    mu = np.broadcast_to(np.asarray(mu, dtype=float), (n,))
    Q = np.kron(Lp, np.eye(2))
    for i in range(n):
        Q[2 * i:2 * i + 2, 2 * i:2 * i + 2] += struct_generator(lam[i], mu[i])
    return Q


def ctmc_stationary(Q: np.ndarray) -> np.ndarray:
    """pi Q = 0, sum(pi) = 1."""
    Q = np.asarray(Q, dtype=float)
    n = Q.shape[0]
    A = np.vstack([Q.T, np.ones(n)])
    b = np.zeros(n + 1)
    b[-1] = 1.0
    pi, *_ = np.linalg.lstsq(A, b, rcond=None)
    pi = np.clip(pi, 0.0, None)
    return pi / pi.sum()


def generator_from_P(P: np.ndarray) -> np.ndarray:
    """Hourly transition matrix -> intensities per hour, lambda_ij = n_ij / T_i = p_ij (i != j).

    This is the thesis estimate; its stationary distribution is that of P.
    """
    P = np.asarray(P, dtype=float)
    return P - np.eye(P.shape[-1])


def availability_by_state(Lp: np.ndarray, lam, mu) -> np.ndarray:
    """P(unit up | production state i) in the stationary integrated chain."""
    pi = ctmc_stationary(block_generator(Lp, lam, mu)).reshape(-1, 2)
    tot = pi.sum(axis=1)
    return np.divide(pi[:, 0], tot, out=np.ones_like(tot), where=tot > 0)


def joint_stationary(Lp: np.ndarray, lam, mu) -> np.ndarray:
    """(N, 2) stationary probabilities of (production state, up/down)."""
    return ctmc_stationary(block_generator(Lp, lam, mu)).reshape(-1, 2)


# --- estimation from event logs ---------------------------------------------------------------

def state_exposure(windows: pd.DataFrame, outages: pd.DataFrame, hour_ns: np.ndarray, states: np.ndarray,
                   n_states: int, min_duration_h: float = 1.0, cells: str | None = None) -> pd.DataFrame:
    """Up/down turbine-hours, failures and repairs per production state.

    ``windows`` (turbine, start, end): observation windows; ``outages`` (turbine, start, end, duration_h): merged
    forced outages; ``hour_ns``: UTC hour starts (int64 ns, consecutive) with their production ``states``
    (-1 = unknown). A failure is an outage of at least ``min_duration_h`` that starts inside a window; a repair
    is its end inside a window. Down-time is the part of those outages inside the windows; up-time the rest of
    the windows. With ``cells`` = "year" the table is split by turbine and calendar year as well (for the
    dispersion estimate).
    """
    hour_ns = np.asarray(hour_ns, dtype=np.int64)
    states = np.asarray(states)
    edges = np.r_[hour_ns, hour_ns[-1] + HOUR_NS].astype(np.float64)
    out = outages[outages["duration_h"] >= min_duration_h]
    years = pd.to_datetime(hour_ns, utc=True).year.to_numpy()
    rows = []
    for t, w in windows.groupby("turbine", sort=True):
        ws, we = _ns(w["start"]), _ns(w["end"])
        o = out[out["turbine"] == t]
        os_, oe = _ns(o["start"]), _ns(o["end"])
        W = np.diff(covered_until(ws, we, edges)) / HOUR_NS
        D = np.minimum(np.diff(covered_until(os_, oe, edges)) / HOUR_NS, W)
        U = W - D
        inside = lambda x: ((x[:, None] >= ws[None, :]) & (x[:, None] < we[None, :])).any(axis=1)  # noqa: E731
        f_idx = _hour_index(os_[inside(os_)], hour_ns) if len(os_) else np.array([], dtype=int)
        r_idx = _hour_index(oe[inside(oe)], hour_ns) if len(oe) else np.array([], dtype=int)
        key = years if cells == "year" else np.zeros(len(hour_ns), dtype=int)
        for k in np.unique(key):
            sel = key == k
            st = np.where(sel, states, -1)
            up = np.bincount(st[st >= 0], weights=U[st >= 0], minlength=n_states)
            dn = np.bincount(st[st >= 0], weights=D[st >= 0], minlength=n_states)
            fs = st[f_idx[(f_idx >= 0)]] if len(f_idx) else np.array([], dtype=int)
            rs = st[r_idx[(r_idx >= 0)]] if len(r_idx) else np.array([], dtype=int)
            nf = np.bincount(fs[fs >= 0], minlength=n_states)
            nr = np.bincount(rs[rs >= 0], minlength=n_states)
            for s in range(n_states):
                rows.append({"turbine": t, "cell": int(k), "state": s, "up_h": up[s], "down_h": dn[s],
                             "failures": int(nf[s]), "repairs": int(nr[s])})
    df = pd.DataFrame(rows)
    if cells is None:
        df = df.groupby("state", as_index=False)[["up_h", "down_h", "failures", "repairs"]].sum()
    return df


def _ns(t) -> np.ndarray:
    return pd.DatetimeIndex(pd.to_datetime(t, utc=True)).as_unit("ns").asi8


def _hour_index(t_ns: np.ndarray, hour_ns: np.ndarray) -> np.ndarray:
    """Index of the hour containing each time (-1 outside the grid)."""
    k = np.searchsorted(hour_ns, t_ns, side="right") - 1
    ok = (k >= 0) & (t_ns < hour_ns[-1] + HOUR_NS)
    return np.where(ok, k, -1)


def rate_table(exp: pd.DataFrame) -> pd.DataFrame:
    """lambda, mu (per hour) with exact Poisson 90 % intervals, and the state FOR = lambda / (lambda + mu)."""
    t = exp.copy()
    t["lambda"] = t["failures"] / t["up_h"].where(t["up_h"] > 0)
    t["mu"] = t["repairs"] / t["down_h"].where(t["down_h"] > 0)
    for k, n, T in (("lambda", "failures", "up_h"), ("mu", "repairs", "down_h")):
        lo, hi = poisson_ci(t[n].to_numpy(), 0.90)
        t[f"{k}_lo"] = lo / t[T].where(t[T] > 0)
        t[f"{k}_hi"] = hi / t[T].where(t[T] > 0)
    t["for_state"] = t["down_h"] / (t["up_h"] + t["down_h"]).where(t["up_h"] + t["down_h"] > 0)
    return t


def poisson_ci(n: np.ndarray, level: float = 0.90) -> tuple[np.ndarray, np.ndarray]:
    """Exact (Garwood) interval of a Poisson count."""
    n = np.asarray(n, dtype=float)
    a = (1 - level) / 2
    lo = np.where(n > 0, stats.chi2.ppf(a, 2 * n) / 2, 0.0)
    hi = stats.chi2.ppf(1 - a, 2 * n + 2) / 2
    return lo, hi


def poisson_lr(n: np.ndarray, T: np.ndarray) -> dict:
    """LR test H0: one rate for all states vs a rate per state (df = states with exposure - 1)."""
    n = np.asarray(n, dtype=float)
    T = np.asarray(T, dtype=float)
    ok = T > 0
    n, T = n[ok], T[ok]
    r0 = n.sum() / T.sum()
    r1 = n / T
    g2 = 2.0 * float(np.sum(np.where(n > 0, n * np.log(np.where(n > 0, r1, 1.0) / r0), 0.0)))
    df = int(ok.sum() - 1)
    # log-linear trend in the state index (df = 1): score test of b in log rate = a + b * state
    x = np.flatnonzero(ok).astype(float)
    e = r0 * T
    xb = (x * e).sum() / e.sum()
    u = float((x * (n - e)).sum())
    v = float((((x - xb) ** 2) * e).sum())
    z2 = u * u / v if v > 0 else float("nan")
    return {"g2": g2, "df": df, "p_value": float(stats.chi2.sf(g2, df)) if df > 0 else float("nan"),
            "trend_z2": z2, "trend_p": float(stats.chi2.sf(z2, 1)) if v > 0 else float("nan"),
            "rate_pooled": r0}


def dispersion(cells: pd.DataFrame, n_col: str, T_col: str) -> float:
    """Pearson chi2 / df of turbine x year x state cells around the state rates (>= 1: no shrinking)."""
    c = cells[cells[T_col] > 0]
    rate = c.groupby("state")[n_col].transform("sum") / c.groupby("state")[T_col].transform("sum")
    e = rate * c[T_col]
    ok = e > 0
    x2 = float((((c[n_col] - e) ** 2) / e)[ok].sum())
    df = int(ok.sum() - c["state"].nunique())
    return max(x2 / df, 1.0) if df > 0 else float("nan")


def common_state_effect(tables: dict[str, pd.DataFrame], n_col: str, T_col: str, iters: int = 500) -> dict:
    """Poisson model rate[f, s] = a_f * b_s fitted by iterative proportional fitting over several farms.

    Returns the relative state effects b_s (mean 1, weighted by exposure), LR tests of b_s = 1 and of the
    common effect against farm-specific state effects.
    """
    farms_ = list(tables)
    n = np.array([tables[f][n_col].to_numpy(dtype=float) for f in farms_])
    T = np.array([tables[f][T_col].to_numpy(dtype=float) for f in farms_])
    a = n.sum(axis=1) / T.sum(axis=1)
    b = np.ones(n.shape[1])
    for _ in range(iters):
        b_new = n.sum(axis=0) / np.maximum((a[:, None] * T).sum(axis=0), 1e-300)
        a = n.sum(axis=1) / np.maximum((T * b_new[None, :]).sum(axis=1), 1e-300)
        if np.allclose(b_new, b, rtol=1e-12, atol=0):
            b = b_new
            break
        b = b_new

    def ll(mu_):
        m = np.where(T > 0, mu_, 0.0)
        return float(np.sum(np.where(n > 0, n * np.log(np.where(m > 0, m, 1.0)), 0.0) - m))

    ll0 = ll(T * (n.sum(axis=1) / T.sum(axis=1))[:, None])
    ll1 = ll(T * a[:, None] * b[None, :])
    ll2 = ll(np.where(T > 0, n, 0.0))  # saturated: a rate per farm and state
    S, F = n.shape[1], n.shape[0]
    g_state = 2 * (ll1 - ll0)
    g_inter = 2 * (ll2 - ll1)
    w = T.sum(axis=0)
    b_norm = b / (np.sum(b * w) / w.sum())
    return {"farms": farms_, "state_effect": b_norm.tolist(),
            "lr_state": g_state, "df_state": S - 1, "p_state": float(stats.chi2.sf(g_state, S - 1)),
            "lr_interaction": g_inter, "df_interaction": (F - 1) * (S - 1),
            "p_interaction": float(stats.chi2.sf(g_inter, (F - 1) * (S - 1)))}


# --- the integrated model in the adequacy framework -------------------------------------------

def monthly_availability(model: mk.MarkovModel, lam, mu) -> np.ndarray:
    """(12, S) P(up | state) for every month of a (seasonal or homogeneous) production chain."""
    return np.stack([availability_by_state(generator_from_P(model.matrix(m)), lam, mu) for m in range(1, 13)])


def integrated_distribution(model: mk.MarkovModel, months: np.ndarray, capacity_mw: float, lam, mu,
                            n_years: int = 1, n_hours: int | None = None) -> elcc.ResDistribution:
    """Hourly output distribution of a plant of many units: level_i x P(up | i), with pi of the month.

    Many independent units given the weather: the available share in state i is P(up | i) (law of large
    numbers), the analogue of the deterministic derating ``elcc.apply_structural(..., "units")``.
    """
    d = elcc.from_markov(model, months, capacity_mw, None, n_years, n_hours)
    A = monthly_availability(model, lam, mu)
    m = np.asarray(months)[: d.levels.shape[1]]
    return elcc.ResDistribution(d.levels * A[m - 1][None], d.probs, capacity_mw)


@dataclass
class IntegratedProduction:
    """Monte Carlo production of a plant whose units fail and are repaired at state-dependent intensities.

    The production chain moves as in ``MarkovProduction``. Each of ``n_units`` units alternates up and down; it
    fails when the cumulative failure intensity since it came up, sum of lambda over the hours spent in each
    state, reaches an Exp(1) draw, and is repaired in the same way with mu (exact for intensities constant
    within an hour). Output = level(state) x share of units up at the start of the hour.
    """

    model: mk.MarkovModel
    months: np.ndarray
    lam: np.ndarray  # (S,) per hour
    mu: np.ndarray
    n_units: int = 100

    def sample(self, rngs, years, n_hours):
        lam = np.asarray(self.lam, dtype=float)
        mu = np.asarray(self.mu, dtype=float)
        A0 = availability_by_state(generator_from_P(self.model.matrix(int(self.months[0]))), lam, mu)
        pi0 = self.model.stationary(int(self.months[0]))
        p_up0 = float(pi0 @ A0)
        # the chain uses the same draws as MarkovProduction; the units continue each year's stream
        u = np.vstack([r.random(n_hours + 1) for r in rngs])
        months = np.asarray(self.months)[:n_hours]
        st = mk.simulate_states(self.model, u, months).astype(np.int64)
        lv = self.model.levels[months - 1][np.arange(n_hours)[None, :], st] if self.model.levels.ndim == 2 \
            else self.model.levels[st]
        out = np.empty((len(rngs), n_hours))
        for k, r in enumerate(rngs):
            out[k] = lv[k] * unit_availability(r, lam[st[k]], mu[st[k]], self.n_units, p_up0)
        return np.clip(out, 0.0, 1.0)


def unit_availability(rng: np.random.Generator, lam_h: np.ndarray, mu_h: np.ndarray, n_units: int,
                      p_up0: float) -> np.ndarray:
    """Share of ``n_units`` units up at each hour start, with hourly intensities ``lam_h``, ``mu_h``."""
    H = len(lam_h)
    cl = np.r_[0.0, np.cumsum(lam_h)]  # cumulative intensity at hour boundaries
    cm = np.r_[0.0, np.cumsum(mu_h)]
    t = np.zeros(n_units)
    up = rng.random(n_units) < p_up0
    diff = np.zeros(H + 1)
    start_up = int(up.sum())
    alive = np.ones(n_units, dtype=bool)
    while alive.any():
        idx = np.flatnonzero(alive)
        c = np.where(up[idx], np.interp(t[idx], np.arange(H + 1), cl), np.interp(t[idx], np.arange(H + 1), cm))
        target = c + rng.standard_exponential(len(idx))
        t_new = np.where(up[idx], _inverse(cl, target), _inverse(cm, target))
        done = t_new >= H
        moved = idx[~done]
        # a unit that changes state during hour h is in its new state from hour ceil(t)
        h = np.ceil(t_new[~done]).astype(np.int64)
        np.add.at(diff, h, np.where(up[moved], -1.0, 1.0))
        t[moved] = t_new[~done]
        up[moved] = ~up[moved]
        alive[idx[done]] = False
    return (start_up + np.cumsum(diff[:H])) / n_units


def _inverse(cum: np.ndarray, target: np.ndarray) -> np.ndarray:
    """Time (in hours, continuous) at which the piecewise-linear cumulative intensity reaches ``target``."""
    H = len(cum) - 1
    k = np.searchsorted(cum, target, side="left")
    out = np.full(len(target), float(H) + 1.0)
    ok = (k >= 1) & (k <= H)
    lo = cum[k[ok] - 1]
    hi = cum[k[ok]]
    step = np.where(hi > lo, (target[ok] - lo) / np.where(hi > lo, hi - lo, 1.0), 0.0)
    out[ok] = k[ok] - 1 + step
    return out
