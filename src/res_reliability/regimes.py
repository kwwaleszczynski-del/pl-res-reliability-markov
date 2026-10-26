"""Regime-switching weather model: daily wind regimes with explicit durations, hourly chains inside them.

A first-order hourly chain forgets its past within a day or two, so it cannot produce the multi-day calms that
matter for storage. This model adds a slow layer:

* **Days.** Every day of the national wind portfolio gets a regime (calm, ..., windy) from its mean capacity
  factor relative to the smoothed climatology of that day of the year. Regime runs are a semi-Markov process:
  the probability that a run ends after its a-th day (discrete hazard) depends on the regime, the season and
  the age of the run, so run lengths follow their empirical distribution instead of the geometric one of a
  first-order chain. When a run ends, the next regime is drawn from the observed switches.
* **Covariate.** Optionally the hazard and the switches also depend on an exogenous class of the day being
  entered (here the tercile of the temperature anomaly within the season). Simulated year i then uses the
  temperatures of weather year i, the same year as its load, so cold and calm days coincide as in the data.
* **Hours.** Within the regime of the day, wind (5 capacity-factor states) and PV (5 clear-sky-index states,
  daylight only) follow hourly chains with one matrix per (month, regime). The regime path is shared by wind
  and PV and by every grid point, which carries the spatial and wind-PV correlation at the daily scale.

Fitting is by counting. Thin cells are shrunk towards the pooled estimate: a hazard cell with n exposures and
e ends gets (e + kappa * h_pooled) / (n + kappa), first pooling over seasons, then over the covariate.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from res_reliability import markov as mk

SEASON_OF_MONTH = np.array([0, 0, 1, 1, 1, 2, 2, 2, 3, 3, 3, 0])  # month 1..12 -> DJF, MAM, JJA, SON
N_SEASONS = 4


def day_months(months_h: np.ndarray) -> np.ndarray:
    """Month of every day from the month of every hour (24 hours per day)."""
    return np.asarray(months_h).reshape(-1, 24)[:, 0]


def daily_mean(x: np.ndarray) -> np.ndarray:
    """(years, hours) -> (years, days)."""
    x = np.asarray(x, dtype=float)
    return x.reshape(x.shape[0], -1, 24).mean(axis=2)


def climatology(daily: np.ndarray, window: int = 31) -> np.ndarray:
    """Mean over years for every day of the year, smoothed by a centred circular moving average."""
    m = np.asarray(daily, dtype=float).mean(axis=0)
    ext = np.r_[m[-window:], m, m[:window]]
    return np.convolve(ext, np.ones(window) / window, mode="same")[window:-window]


def quantile_classes(values: np.ndarray, quantiles, groups: np.ndarray | None = None) -> np.ndarray:
    """Class 0..len(quantiles) of every value by quantile edges, computed within each group of days."""
    v = np.asarray(values, dtype=float)
    out = np.empty(v.shape, dtype=np.int64)
    if groups is None:
        return np.digitize(v, np.quantile(v, quantiles))
    groups = np.asarray(groups)
    for g in np.unique(groups):
        sel = np.broadcast_to(groups == g, v.shape)
        out[sel] = np.digitize(v[sel], np.quantile(v[sel], quantiles))
    return out


@dataclass
class RegimeProcess:
    """Semi-Markov process of daily regimes.

    ``hazard[r, s, c, a]``: probability that a run of regime r ends after its (a+1)-th day, when the next day
    is in season s with covariate class c (the last age bin is open-ended: ages >= max_age).
    ``switch[r, s, c, :]``: distribution of the next regime when a run ends (zero on the diagonal).
    ``init[s, :]``: regime frequencies per season, used to start a simulation (before a burn-in).
    """

    hazard: np.ndarray
    switch: np.ndarray
    init: np.ndarray

    @property
    def n_regimes(self) -> int:
        return self.hazard.shape[0]

    @property
    def n_cov(self) -> int:
        return self.hazard.shape[2]

    @property
    def max_age(self) -> int:
        return self.hazard.shape[3]


def _shrink(events: np.ndarray, exposure: np.ndarray, prior: np.ndarray, kappa: float) -> np.ndarray:
    return (events + kappa * prior) / (exposure + kappa)


def run_ages(regimes: np.ndarray) -> np.ndarray:
    """Age (1, 2, ...) of every day within its run of equal regimes."""
    r = np.asarray(regimes)
    age = np.ones(len(r), dtype=np.int64)
    for t in range(1, len(r)):
        if r[t] == r[t - 1]:
            age[t] = age[t - 1] + 1
    return age


def fit_process(
    regimes: np.ndarray,
    seasons: np.ndarray,
    n_regimes: int,
    cov: np.ndarray | None = None,
    n_cov: int = 1,
    max_age: int = 12,
    kappa: float = 10.0,
    valid_step: np.ndarray | None = None,
) -> RegimeProcess:
    """Fit hazards and switch probabilities to a sequence of daily regimes.

    ``seasons`` and ``cov`` are given per day; the step from day t to t+1 is attributed to the season and
    covariate class of day t+1. ``valid_step[t]`` False drops the step t -> t+1 (e.g. across a gap).
    """
    r = np.asarray(regimes, dtype=np.int64)
    s = np.asarray(seasons, dtype=np.int64)
    c = np.zeros(len(r), dtype=np.int64) if cov is None else np.asarray(cov, dtype=np.int64)
    age = np.minimum(run_ages(r), max_age) - 1
    ok = np.ones(len(r) - 1, dtype=bool) if valid_step is None else np.asarray(valid_step, bool)[: len(r) - 1]
    a, b = r[:-1][ok], r[1:][ok]
    ss, cc, ag = s[1:][ok], c[1:][ok], age[:-1][ok]
    end = (a != b).astype(float)
    shape = (n_regimes, N_SEASONS, n_cov, max_age)
    n = np.zeros(shape)
    e = np.zeros(shape)
    np.add.at(n, (a, ss, cc, ag), 1.0)
    np.add.at(e, (a, ss, cc, ag), end)
    # pooling: regime -> (regime, age) -> (regime, season, age) -> (regime, season, cov, age)
    h_r = e.sum(axis=(1, 2, 3)) / np.maximum(n.sum(axis=(1, 2, 3)), 1)
    h_ra = _shrink(e.sum(axis=(1, 2)), n.sum(axis=(1, 2)), h_r[:, None], kappa)
    h_rsa = _shrink(e.sum(axis=2), n.sum(axis=2), h_ra[:, None, :], kappa)
    hazard = _shrink(e, n, h_rsa[:, :, None, :], kappa)

    sw = np.zeros((n_regimes, N_SEASONS, n_cov, n_regimes))
    np.add.at(sw, (a[end > 0], ss[end > 0], cc[end > 0], b[end > 0]), 1.0)
    p_r = sw.sum(axis=(1, 2)) + 1e-9 * (1 - np.eye(n_regimes))
    p_r = p_r / p_r.sum(axis=1, keepdims=True)
    tot_rs = sw.sum(axis=2)
    p_rs = (tot_rs + kappa * p_r[:, None, :]) / (tot_rs.sum(axis=2, keepdims=True) + kappa)
    switch = (sw + kappa * p_rs[:, :, None, :]) / (sw.sum(axis=3, keepdims=True) + kappa)

    init = np.zeros((N_SEASONS, n_regimes))
    np.add.at(init, (s, r), 1.0)
    init = (init + 1e-9) / (init + 1e-9).sum(axis=1, keepdims=True)
    return RegimeProcess(hazard, switch, init)


def _draw(cum: np.ndarray, u: np.ndarray) -> np.ndarray:
    """Row-wise inverse-CDF draw: cum (n, K) cumulative probabilities, u (n,) uniforms."""
    return np.minimum((u[:, None] > cum).sum(axis=1), cum.shape[1] - 1)


def simulate_process(
    proc: RegimeProcess,
    u: np.ndarray,
    seasons: np.ndarray,
    cov: np.ndarray | None = None,
) -> np.ndarray:
    """Simulate daily regimes for many years at once.

    ``u`` is (n_sims, n_days, 2) uniforms (hazard draw, next-regime draw), ``seasons`` (n_days,), ``cov``
    (n_sims, n_days) covariate classes or None. Day 0 is drawn from the regime frequencies of its season with
    age 1; callers discard a burn-in so that the residual run lengths settle. Returns (n_sims, n_days) int64.
    """
    u = np.asarray(u)
    n_sims, n_days = u.shape[:2]
    seasons = np.asarray(seasons, dtype=np.int64)
    c = np.zeros((n_sims, n_days), dtype=np.int64) if cov is None else np.asarray(cov, dtype=np.int64)
    cum_init = np.cumsum(proc.init, axis=1)
    cum_sw = np.cumsum(proc.switch, axis=3)
    cum_sw[..., -1] = 1.0
    out = np.empty((n_sims, n_days), dtype=np.int64)
    r = _draw(np.broadcast_to(cum_init[seasons[0]], (n_sims, proc.n_regimes)), u[:, 0, 0])
    age = np.ones(n_sims, dtype=np.int64)
    out[:, 0] = r
    for d in range(1, n_days):
        s, cd = seasons[d], c[:, d]
        h = proc.hazard[r, s, cd, np.minimum(age, proc.max_age) - 1]
        end = u[:, d, 0] < h
        nxt = _draw(cum_sw[r, s, cd], u[:, d, 1])
        r = np.where(end, nxt, r)
        age = np.where(end, 1, age + 1)
        out[:, d] = r
    return out


@dataclass
class HourlyChains:
    """Hourly transition matrices ``P[month-1, regime]`` (S x S) and state levels ``levels[month-1, regime]``."""

    P: np.ndarray  # (12, K, S, S)
    levels: np.ndarray  # (12, K, S)

    @property
    def n_states(self) -> int:
        return self.P.shape[-1]


def fit_hourly(
    x: np.ndarray,
    months: np.ndarray,
    regimes_h: np.ndarray,
    n_regimes: int,
    n_states: int = 5,
    valid: np.ndarray | None = None,
    level_kind: str = "mean",
) -> HourlyChains:
    """One chain per (month, regime) on a flat hourly series.

    The transition into an hour belongs to that hour's month and regime (the day being entered decides where
    the chain goes). Invalid hours (PV at night) are skipped, so the chain jumps from the last valid hour to
    the next one. Rows without data fall back to the month's matrix pooled over regimes.
    """
    x = np.asarray(x, dtype=float)
    months = np.asarray(months, dtype=np.int64)
    rg = np.asarray(regimes_h, dtype=np.int64)
    valid = np.ones(len(x), dtype=bool) if valid is None else np.asarray(valid, dtype=bool)
    states = mk.discretize(np.where(valid, x, np.nan), n_states)
    keep = valid & (states >= 0)
    idx = np.flatnonzero(keep)
    a, b = states[idx[:-1]], states[idx[1:]]
    m_to, r_to = months[idx[1:]] - 1, rg[idx[1:]]
    counts = np.zeros((12, n_regimes, n_states, n_states))
    np.add.at(counts, (m_to, r_to, a, b), 1.0)
    pooled = counts.sum(axis=1)
    P = np.empty_like(counts)
    levels = np.empty((12, n_regimes, n_states))
    centres = np.arange(n_states) / (n_states - 1)
    for m in range(12):
        in_m = keep & (months == m + 1)
        freq = np.bincount(states[in_m], minlength=n_states) / max(in_m.sum(), 1)
        P_m = mk.normalise_rows(pooled[m], fallback=freq)
        lv_m = mk.state_levels(x[in_m], states[in_m], n_states, level_kind) if in_m.any() else centres
        for k in range(n_regimes):
            P[m, k] = mk.normalise_rows(counts[m, k])
            empty = counts[m, k].sum(axis=1) == 0
            P[m, k, empty] = P_m[empty]
            sel = in_m & (rg == k)
            lv = lv_m.copy()
            if level_kind == "mean":
                for j in range(n_states):
                    sj = sel & (states == j)
                    if sj.any():
                        lv[j] = float(x[sj].mean())
            levels[m, k] = lv
    return HourlyChains(P, levels)


def simulate_hourly(
    chains: HourlyChains,
    regimes_h: np.ndarray,
    months: np.ndarray,
    u: np.ndarray,
    active: np.ndarray | None = None,
) -> np.ndarray:
    """Hourly states of many chains, each following its own regime path.

    ``regimes_h`` (n_sims, n_hours) regime of every hour, ``months`` (n_hours,), ``u`` (n_sims, n_hours + 1):
    column 0 draws the initial state from the stationary distribution of the first hour's (month, regime),
    column h+1 drives the step into hour h. ``active`` (n_hours,) or (n_sims, n_hours): the chain only moves
    on active hours and keeps its state otherwise. Returns (n_sims, n_hours) int8 states.
    """
    rg = np.asarray(regimes_h, dtype=np.int64)
    n_sims, n_hours = rg.shape
    months = np.asarray(months, dtype=np.int64)
    S = chains.n_states
    act = np.ones((1, n_hours), dtype=bool) if active is None else np.atleast_2d(np.asarray(active, dtype=bool))
    cum = np.cumsum(chains.P, axis=3)
    cum[..., -1] = 1.0
    m0 = months[0] - 1
    pi0 = np.stack([np.cumsum(mk.stationary(chains.P[m0, k])) for k in range(chains.P.shape[1])])
    state = _draw(pi0[rg[:, 0]], u[:, 0])
    out = np.empty((n_sims, n_hours), dtype=np.int8)
    started = np.zeros(act.shape[0], dtype=bool)
    for h in range(n_hours):
        a = act[:, h]
        move = a & started
        if move.any():
            nxt = _draw(cum[months[h] - 1, rg[:, h], state], u[:, h + 1])
            state = np.where(np.broadcast_to(move, (n_sims,)), nxt, state)
        started |= a
        out[:, h] = state
    return out.clip(max=S - 1)


def levels_of(chains: HourlyChains, states: np.ndarray, regimes_h: np.ndarray, months: np.ndarray) -> np.ndarray:
    """Output level of every simulated hour (n_sims, n_hours)."""
    return chains.levels[np.asarray(months)[None, :] - 1, regimes_h, states]


@dataclass
class RegimeWeather:
    """National regime model of wind and PV, fitted on the portfolio series of the weather years.

    ``cov`` (n_years, 365) is the covariate class of every day of every weather year, or None for the
    unconditioned model. Simulated year i uses weather year i % n_years and is preceded by a burn-in of the
    last ``burn_in_days`` of the previous weather year.
    """

    process: RegimeProcess
    wind: HourlyChains
    pv: HourlyChains
    months: np.ndarray  # (8760,) month of every study hour
    edges: np.ndarray  # regime edges of the daily capacity-factor ratio to climatology
    cov: np.ndarray | None = None
    burn_in_days: int = 30
    hourly_cov: bool = False  # hourly chains per (month, regime, covariate class) instead of (month, regime)

    @property
    def n_regimes(self) -> int:
        return self.process.n_regimes

    @property
    def days(self) -> int:
        return len(self.months) // 24

    def extended_days(self) -> tuple[np.ndarray, np.ndarray]:
        """Seasons and months of the burn-in days followed by the study year."""
        md = day_months(self.months)
        md = np.r_[md[-self.burn_in_days:], md]
        return SEASON_OF_MONTH[md - 1], md

    def extended_cov(self, years: np.ndarray) -> np.ndarray | None:
        """(n_sims, burn_in_days + 365) covariate classes: end of the previous weather year, then the year."""
        if self.cov is None:
            return None
        ny = self.cov.shape[0]
        y = np.asarray(years) % ny
        return np.concatenate([self.cov[(y - 1) % ny, -self.burn_in_days:], self.cov[y]], axis=1)

    def regime_paths(self, rngs: list[np.random.Generator], years: np.ndarray) -> np.ndarray:
        """(n_sims, burn_in_days + 365) daily regimes; the first ``burn_in_days`` belong to the previous year."""
        seasons, _ = self.extended_days()
        u = np.stack([r.random((len(seasons), 2)) for r in rngs])
        return simulate_process(self.process, u, seasons, self.extended_cov(years))

    def hourly_classes(self, regimes_d: np.ndarray, years: np.ndarray) -> np.ndarray:
        """Index of the hourly chain of every simulated day (the regime, or regime x covariate class)."""
        if not self.hourly_cov or self.cov is None:
            return regimes_d
        return regimes_d * self.process.n_cov + self.extended_cov(years)

    def hourly(self, tech: str, rngs: list[np.random.Generator], classes_d: np.ndarray,
               active: np.ndarray | None = None) -> np.ndarray:
        """Normalised hourly output levels (wind CF, or PV clear-sky index) on the extended timeline.

        ``classes_d`` from ``hourly_classes``; ``active`` (n_sims, n_ext_hours) for PV (daylight), None for
        wind. Returns (n_sims, n_ext_hours).
        """
        chains = self.wind if tech == "wind" else self.pv
        _, md = self.extended_days()
        months_h = np.repeat(md, 24)
        rg_h = np.repeat(classes_d, 24, axis=1)
        u = np.stack([r.random(len(months_h) + 1) for r in rngs])
        states = simulate_hourly(chains, rg_h, months_h, u, active)
        return levels_of(chains, states, rg_h, months_h)


def fit_weather(
    wind: np.ndarray,
    pv_k: np.ndarray,
    daylight: np.ndarray,
    months: np.ndarray,
    quantiles=(0.2, 0.5, 0.8),
    temperature: np.ndarray | None = None,
    cov_quantiles=(1 / 3, 2 / 3),
    n_states: int = 5,
    max_age: int = 12,
    kappa: float = 10.0,
    burn_in_days: int = 30,
    hourly_cov: bool = False,
) -> RegimeWeather:
    """Fit the regime model to (years, 8760) portfolio series.

    ``temperature`` (years, 8760) adds the covariate: the class (quantiles within the season) of the daily
    temperature anomaly relative to the smoothed day-of-year climatology. With ``hourly_cov`` the hourly chains
    are also fitted per covariate class.
    """
    n_years = wind.shape[0]
    wd = daily_mean(wind)
    ratio = wd / climatology(wd)
    edges = np.quantile(ratio, quantiles)
    reg = np.digitize(ratio, edges)  # (years, days)
    md = day_months(months)
    seasons_d = SEASON_OF_MONTH[md - 1]
    cov = None
    n_cov = 1
    if temperature is not None:
        td = daily_mean(temperature)
        cov = quantile_classes(td - climatology(td), cov_quantiles, np.broadcast_to(seasons_d, td.shape))
        n_cov = len(cov_quantiles) + 1
    K = len(quantiles) + 1
    proc = fit_process(reg.reshape(-1), np.tile(seasons_d, n_years), K,
                       None if cov is None else cov.reshape(-1), n_cov, max_age, kappa)
    hc = hourly_cov and cov is not None
    cls = reg * n_cov + cov if hc else reg
    n_cls = K * n_cov if hc else K
    rg_h = np.repeat(cls, 24, axis=1).reshape(-1)
    m_h = np.tile(months, n_years)
    w = fit_hourly(wind.reshape(-1), m_h, rg_h, n_cls, n_states)
    p = fit_hourly(pv_k.reshape(-1), m_h, rg_h, n_cls, n_states, valid=daylight.reshape(-1))
    return RegimeWeather(proc, w, p, np.asarray(months), edges, cov, burn_in_days, hc)
