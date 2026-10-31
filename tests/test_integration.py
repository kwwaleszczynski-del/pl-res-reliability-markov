import numpy as np
import pandas as pd
import pytest

from res_reliability import elcc, integration as ig, markov as mk
from res_reliability.farms import covered_until, interval_coverage


def random_generator(n=5, seed=1):
    rng = np.random.default_rng(seed)
    P = rng.random((n, n)) + 4 * np.eye(n)
    P /= P.sum(axis=1, keepdims=True)
    return ig.generator_from_P(P), P


def test_kronecker_sum_with_constant_rates_is_the_product_of_the_marginals():
    Lp, P = random_generator()
    lam, mu = 1 / 900.0, 1 / 15.0
    pi = ig.ctmc_stationary(ig.kronecker_generator(Lp, lam, mu))
    pi_prod = mk.stationary(P)
    pi_struct = np.array([mu, lam]) / (lam + mu)
    assert np.allclose(pi, np.kron(pi_prod, pi_struct), atol=1e-13)
    # the thesis independence formula follows: P(state i, available) = p_prod(i) * A
    assert np.allclose(pi.reshape(-1, 2)[:, 0], pi_prod * (1 - lam / (lam + mu)), atol=1e-13)


def test_block_generator_reduces_to_kronecker_sum():
    Lp, _ = random_generator(seed=2)
    Q1 = ig.kronecker_generator(Lp, 0.002, 0.08)
    Q2 = ig.block_generator(Lp, np.full(5, 0.002), np.full(5, 0.08))
    assert np.allclose(Q1, Q2)
    assert np.allclose(Q2.sum(axis=1), 0.0)


def test_state_dependent_rates_correlate_availability_with_production():
    Lp, P = random_generator(seed=3)
    lam = np.array([0.5, 1, 1, 2, 4]) * 1e-3
    mu = np.array([0.10, 0.10, 0.08, 0.04, 0.02])
    joint = ig.joint_stationary(Lp, lam, mu)
    # production transitions do not depend on the unit, so the production marginal is unchanged
    assert np.allclose(joint.sum(axis=1), mk.stationary(P), atol=1e-12)
    A = ig.availability_by_state(Lp, lam, mu)
    assert A[0] > A[-1]
    # the chain mixes the states, so the unit's lagged state pulls A towards the mean
    quasi_static = mu / (lam + mu)
    assert A[-1] > quasi_static[-1] and A[0] < quasi_static[0]


def test_covered_until_and_interval_coverage():
    s = np.array([1.0, 5.0])
    e = np.array([3.0, 6.0])
    assert covered_until(s, e, np.array([0.0, 2.0, 4.0, 5.5, 10.0])).tolist() == [0.0, 1.0, 2.0, 2.5, 3.0]
    assert interval_coverage(s, e, np.arange(0.0, 8.0)).tolist() == [0, 1, 1, 0, 0, 1, 0]
    assert covered_until(np.array([]), np.array([]), np.array([1.0])).tolist() == [0.0]


def simulate_log(lam, mu, n_turbines=40, n_hours=30000, seed=5):
    """Synthetic outages of turbines whose intensities follow an hourly production state."""
    rng = np.random.default_rng(seed)
    P = np.array([[0.95, 0.04, 0.01], [0.03, 0.94, 0.03], [0.01, 0.05, 0.94]])
    u = rng.random((1, n_hours + 1))
    states = mk.simulate_states(mk.MarkovModel(P, np.array([0.0, 0.5, 1.0]), 3), u,
                                np.ones(n_hours, dtype=int))[0].astype(int)
    hours = np.arange(n_hours + 1, dtype=float)
    cl, cm = np.r_[0.0, np.cumsum(lam[states])], np.r_[0.0, np.cumsum(mu[states])]
    t0 = pd.Timestamp("2020-01-01", tz="UTC")
    rows = []
    for k in range(n_turbines):
        t, up = 0.0, True
        while True:
            c = np.interp(t, hours, cl if up else cm) + rng.standard_exponential()
            t_new = float(np.interp(c, cl if up else cm, hours)) if c < (cl if up else cm)[-1] else np.inf
            if t_new >= n_hours:
                if not up:
                    rows.append((k, t, float(n_hours)))
                break
            if up:
                t_fail = t_new
            else:
                rows.append((k, t_fail, t_new))
            t, up = t_new, not up
    out = pd.DataFrame(rows, columns=["turbine", "a", "b"])
    out["start"] = t0 + pd.to_timedelta(out["a"], unit="h")
    out["end"] = t0 + pd.to_timedelta(out["b"], unit="h")
    out["duration_h"] = out["b"] - out["a"]
    win = pd.DataFrame({"turbine": range(n_turbines), "start": t0, "end": t0 + pd.Timedelta(hours=n_hours)})
    hour_ns = (t0 + pd.to_timedelta(np.arange(n_hours), unit="h")).as_unit("ns").asi8
    return out[["turbine", "start", "end", "duration_h"]], win, hour_ns, states


def test_rates_are_recovered_from_a_synthetic_log():
    lam = np.array([0.5e-3, 1e-3, 3e-3])
    mu = np.array([0.10, 0.05, 0.02])
    out, win, hour_ns, states = simulate_log(lam, mu)
    exp = ig.state_exposure(win, out, hour_ns, states, 3, min_duration_h=0.0)
    tab = ig.rate_table(exp)
    assert np.isclose(exp["up_h"].sum() + exp["down_h"].sum(), 40 * 30000)
    for k, true in (("lambda", lam), ("mu", mu)):
        assert np.all(tab[f"{k}_lo"] * 0.9 <= true) and np.all(true <= tab[f"{k}_hi"] * 1.1)
    assert ig.poisson_lr(tab["failures"], tab["up_h"])["p_value"] < 1e-6
    assert ig.poisson_lr(tab["repairs"], tab["down_h"])["p_value"] < 1e-6


def test_constant_rates_are_not_rejected():
    out, win, hour_ns, states = simulate_log(np.full(3, 1e-3), np.full(3, 0.05), seed=11)
    tab = ig.rate_table(ig.state_exposure(win, out, hour_ns, states, 3, min_duration_h=0.0))
    assert ig.poisson_lr(tab["failures"], tab["up_h"])["p_value"] > 0.01
    assert ig.poisson_lr(tab["repairs"], tab["down_h"])["p_value"] > 0.01
    cells = ig.state_exposure(win, out, hour_ns, states, 3, min_duration_h=0.0, cells="year")
    assert 0.7 < ig.dispersion(cells, "failures", "up_h") < 1.5


def test_common_state_effect_matches_the_generating_ratio():
    lam = np.array([1e-3, 2e-3, 4e-3])
    tabs = {}
    for f, scale in (("a", 1.0), ("b", 3.0)):
        out, win, hour_ns, states = simulate_log(lam * scale, np.full(3, 0.05), seed=int(scale * 7))
        tabs[f] = ig.rate_table(ig.state_exposure(win, out, hour_ns, states, 3, min_duration_h=0.0))
    r = ig.common_state_effect(tabs, "failures", "up_h")
    b = np.array(r["state_effect"])
    assert b[2] / b[0] == pytest.approx(4.0, rel=0.25)
    assert r["p_state"] < 1e-6 and r["p_interaction"] > 0.01


def test_integrated_distribution_with_constant_rates_is_the_units_derating():
    rng = np.random.default_rng(4)
    x = np.clip(np.cumsum(rng.normal(0, 0.05, 20000)) % 1.0, 0, 1)
    months = np.tile(np.repeat(np.arange(1, 13), 730)[:8760], 3)[:20000]
    model = mk.fit(x, 5, months=months)
    hour_months = np.repeat(np.arange(1, 13), 730)
    lam, mu = 1 / 1000.0, 1 / 20.0
    d1 = ig.integrated_distribution(model, hour_months, 100.0, np.full(5, lam), np.full(5, mu))
    d2 = elcc.apply_structural(elcc.from_markov(model, hour_months, 100.0), lam / (lam + mu), "units")
    assert np.allclose(d1.levels, d2.levels) and np.allclose(d1.probs, d2.probs)


def test_unit_availability_matches_the_stationary_value():
    rng = np.random.default_rng(9)
    lam, mu = 1 / 200.0, 1 / 10.0
    a = ig.unit_availability(rng, np.full(200000, lam), np.full(200000, mu), 50, mu / (lam + mu))
    assert a.mean() == pytest.approx(mu / (lam + mu), abs=0.003)
    # state-dependent: no repair while the intensity is zero
    a2 = ig.unit_availability(rng, np.r_[np.full(100, 0.5), np.full(100, 0.0)],
                              np.r_[np.full(100, 0.0), np.full(100, 0.0)], 20, 1.0)
    assert a2[0] == 1.0 and a2[-1] == 0.0


def test_integrated_production_mean_agrees_with_the_analytic_distribution():
    P = np.array([[0.9, 0.1, 0.0], [0.05, 0.9, 0.05], [0.0, 0.1, 0.9]])
    model = mk.MarkovModel(P, np.array([0.0, 0.5, 1.0]), 3)
    months = np.repeat(np.arange(1, 13), 730)
    lam = np.array([0.0, 0.01, 0.04])
    mu = np.array([0.2, 0.1, 0.02])
    prod = ig.IntegratedProduction(model, months, lam, mu, n_units=50)
    rngs = [np.random.default_rng(i) for i in range(40)]
    x = prod.sample(rngs, np.arange(40), 8760)
    d = ig.integrated_distribution(model, months, 1.0, lam, mu)
    assert x.mean() == pytest.approx(d.mean_mw, rel=0.03)
