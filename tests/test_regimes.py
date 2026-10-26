import numpy as np
import pytest

from res_reliability import markov as mk, regimes as rgm, validation

MONTH_HOURS = np.repeat(np.arange(1, 13), np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]) * 24)
CALM_PMF = np.array([0.30, 0.05, 0.05, 0.10, 0.25, 0.15, 0.10])  # lull lengths 1..7 days, far from geometric


def alternating_regimes(n_days, rng, cov=None):
    """Calm (0) runs with lengths from CALM_PMF, windy (1) runs of 1-4 days. With ``cov`` (per day), calm runs
    that start on a class-1 day are drawn 3 days longer."""
    out = np.empty(n_days, dtype=np.int64)
    t, calm = 0, True
    while t < n_days:
        d = rng.choice(np.arange(1, 8), p=CALM_PMF) if calm else rng.integers(1, 5)
        if calm and cov is not None and cov[t] == 1:
            d += 3
        out[t:t + d] = 0 if calm else 1
        t += d
        calm = not calm
    return out


def run_pmf(runs, n=7):
    return np.bincount(np.minimum(runs, n), minlength=n + 1)[1:] / len(runs)


def test_semi_markov_process_reproduces_a_given_run_length_distribution():
    rng = np.random.default_rng(0)
    reg = alternating_regimes(60 * 365, rng)
    proc = rgm.fit_process(reg, np.zeros(len(reg), dtype=np.int64), 2, max_age=10)
    sim = rgm.simulate_process(proc, rng.random((20, 365 * 10, 2)), np.zeros(365 * 10, dtype=np.int64))
    runs = validation.low_runs((sim[:, 50:] != 0).astype(float), 0.5)
    assert np.abs(run_pmf(runs) - CALM_PMF).max() < 0.02
    # a first-order chain fitted to the same days only knows the mean length: geometric, too few 5-7 day runs
    P = mk.normalise_rows(mk.transition_counts(reg, 2))
    geo = P[0, 0] ** np.arange(7) * (1 - P[0, 0])
    assert geo[4:].sum() < 0.5 * CALM_PMF[4:].sum()


def synthetic_cf(n_years, rng):
    """Hourly wind CF: calm days around 5 %, windy days around 50 %, regime runs from ``alternating_regimes``."""
    reg = alternating_regimes(n_years * 365, rng).reshape(n_years, 365)
    level = np.where(reg == 0, 0.05, 0.5).repeat(24, axis=1)
    return np.clip(level + rng.normal(0, 0.04, level.shape), 0, 1), reg


def test_regime_model_reproduces_lull_lengths_and_the_monthly_chain_does_not():
    rng = np.random.default_rng(1)
    wind, reg = synthetic_cf(30, rng)
    k = np.clip(rng.uniform(0.2, 1.0, wind.shape), 0, 1)
    day = np.tile((np.arange(24) >= 7) & (np.arange(24) < 17), (30, 365))
    calm_share = float((reg == 0).mean())
    w = rgm.fit_weather(wind, k, day, MONTH_HOURS, quantiles=(calm_share,))
    years = np.arange(60) % 30
    regs = w.regime_paths([np.random.default_rng([2, i]) for i in range(60)], years)
    x = w.hourly("wind", [np.random.default_rng([3, i]) for i in range(60)], regs)[:, w.burn_in_days * 24:]
    sim_runs = validation.low_runs(rgm.daily_mean(x), 0.15)
    data_runs = validation.low_runs(rgm.daily_mean(wind), 0.15)
    assert np.abs(run_pmf(sim_runs) - CALM_PMF).max() < 0.04
    assert np.abs(run_pmf(data_runs) - CALM_PMF).max() < 0.04
    # hourly first-order monthly chain on the same series
    m = mk.fit(wind.reshape(-1), 5, months=np.tile(MONTH_HOURS, 30))
    s = mk.simulate_states(m, np.random.default_rng(4).random((60, 8761)), MONTH_HOURS)
    chain = m.levels[MONTH_HOURS - 1][np.arange(8760)[None, :], s]
    chain_runs = validation.low_runs(rgm.daily_mean(chain), 0.15)
    long_ = lambda r: (r >= 5).sum() / 60  # noqa: E731
    assert long_(sim_runs) == pytest.approx(long_(data_runs) * 2, rel=0.15)  # 60 synthetic vs 30 data years
    # the chain only keeps the mean lull length: geometric shape, about half of the 5-7 day lulls
    assert np.abs(run_pmf(chain_runs) - CALM_PMF).max() > 0.1
    assert long_(chain_runs) < 0.7 * long_(sim_runs)


def test_covariate_changes_the_hazard_in_the_right_direction():
    rng = np.random.default_rng(5)
    n = 80 * 365
    cov = np.repeat(rng.integers(0, 2, n // 5 + 1), 5)[:n]  # covariate classes persist for 5 days
    reg = alternating_regimes(n, rng, cov)
    proc = rgm.fit_process(reg, np.zeros(n, dtype=np.int64), 2, cov=cov, n_cov=2, max_age=12)
    h0, h1 = proc.hazard[0, 0, 0], proc.hazard[0, 0, 1]
    assert h1[:3].mean() < 0.5 * h0[:3].mean()  # calm runs in class 1 rarely end in their first days
    # simulated with the covariate, calm days pile up in class 1 as in the data
    c = np.tile(cov[: 365 * 4], (10, 1))
    sim = rgm.simulate_process(proc, rng.random((10, 365 * 4, 2)), np.zeros(365 * 4, dtype=np.int64), c)
    data_share = [(reg[cov == k] == 0).mean() for k in (0, 1)]
    sim_share = [(sim[c == k] == 0).mean() for k in (0, 1)]
    assert sim_share[1] - sim_share[0] > 0.5 * (data_share[1] - data_share[0]) > 0


def test_hazard_and_switch_tables_are_probabilities():
    rng = np.random.default_rng(6)
    reg = rng.integers(0, 3, 2000)
    seasons = rgm.SEASON_OF_MONTH[np.repeat(np.arange(12), 170)[:2000]]
    proc = rgm.fit_process(reg, seasons, 3, cov=rng.integers(0, 3, 2000), n_cov=3)
    assert proc.hazard.shape == (3, 4, 3, 12)
    assert np.all((proc.hazard >= 0) & (proc.hazard <= 1))
    assert np.allclose(proc.switch.sum(axis=3), 1.0)
    assert np.all(proc.switch[np.arange(3), :, :, np.arange(3)] == 0)  # a run ends in a different regime


def test_hourly_chains_follow_the_regime_of_the_day():
    P = np.zeros((12, 2, 2, 2))
    P[:, 0] = [[1.0, 0.0], [1.0, 0.0]]  # regime 0 pulls the chain into state 0
    P[:, 1] = [[0.0, 1.0], [0.0, 1.0]]  # regime 1 into state 1
    ch = rgm.HourlyChains(P, np.tile([0.0, 1.0], (12, 2, 1)))
    rg = np.repeat([[0, 1, 0]], 2, axis=0).repeat(24, axis=1)
    u = np.random.default_rng(0).random((2, 73))
    s = rgm.simulate_hourly(ch, rg, np.ones(72, dtype=int), u)
    assert np.all(s[:, 1:24] == 0) and np.all(s[:, 24:48] == 1) and np.all(s[:, 48:] == 0)
    active = np.zeros(72, dtype=bool)
    active[[0, 30, 60]] = True  # PV-like: moves only on active hours
    s = rgm.simulate_hourly(ch, rg, np.ones(72, dtype=int), u, active)
    assert np.all(s[:, 0:30] == s[:, [0]]) and np.all(s[:, 30:60] == 1) and np.all(s[:, 60:] == 0)


def test_fit_hourly_matrices_are_stochastic_and_levels_preserve_energy():
    rng = np.random.default_rng(7)
    x = np.clip(rng.beta(1.2, 3, 8760 * 3), 0, 1)
    reg = np.repeat(rng.integers(0, 3, 365 * 3), 24)
    ch = rgm.fit_hourly(x, np.tile(MONTH_HOURS, 3), reg, 3)
    assert ch.P.shape == (12, 3, 5, 5) and np.allclose(ch.P.sum(axis=3), 1.0)
    states = mk.discretize(x, 5)
    lv = ch.levels[np.tile(MONTH_HOURS, 3) - 1, reg, states]
    assert lv.mean() == pytest.approx(x.mean(), abs=1e-9)


def test_regime_paths_are_reproducible_and_seed_dependent():
    rng = np.random.default_rng(8)
    wind, _ = synthetic_cf(4, rng)
    day = np.tile((np.arange(24) >= 7) & (np.arange(24) < 17), (4, 365))
    temp = rng.normal(0, 5, wind.shape)
    w = rgm.fit_weather(wind, np.full(wind.shape, 0.7), day, MONTH_HOURS, temperature=temp)
    years = np.array([0, 1, 2, 3])
    a = w.regime_paths([np.random.default_rng([1, i]) for i in range(4)], years)
    b = w.regime_paths([np.random.default_rng([1, i]) for i in range(4)], years)
    c = w.regime_paths([np.random.default_rng([2, i]) for i in range(4)], years)
    assert a.shape == (4, w.burn_in_days + 365)
    assert np.array_equal(a, b) and not np.array_equal(a, c)
    assert w.extended_cov(years).shape == a.shape


def test_lull_helpers():
    daily = np.array([[0.05, 0.05, 0.3, 0.05, 0.3, 0.3, 0.05, 0.05, 0.05]])
    runs = validation.low_runs(daily, 0.1)
    assert sorted(runs.tolist()) == [1, 2, 3]  # runs at the start and the end of the year count
    assert validation.run_length_counts(runs, 1, 3).tolist() == [1, 1, 1]
    assert validation.rolling_min_mean(daily, 3)[0] == pytest.approx(0.05)
    assert rgm.climatology(np.ones((3, 365))).shape == (365,)
