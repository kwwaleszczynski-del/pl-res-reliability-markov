import numpy as np
import pytest

from res_reliability import markov as mk


def ar1_cf(n=20000, phi=0.95, seed=3):
    rng = np.random.default_rng(seed)
    z = np.empty(n)
    z[0] = 0.0
    e = rng.normal(size=n) * np.sqrt(1 - phi**2)
    for t in range(1, n):
        z[t] = phi * z[t - 1] + e[t]
    from scipy.stats import norm
    return norm.cdf(z) ** 1.5  # skewed towards low output, like wind


def test_thresholds_match_thesis():
    assert np.allclose(mk.thresholds(5), [0.125, 0.375, 0.625, 0.875])
    s = mk.discretize(np.array([0.0, 0.12, 0.13, 0.6, 0.7, 0.9, 1.0, np.nan]), 5)
    assert s.tolist() == [0, 0, 1, 2, 3, 4, 4, -1]


def test_transition_matrix_is_stochastic():
    x = ar1_cf()
    m = mk.fit(x, 5)
    assert m.P.shape == (5, 5)
    assert np.all(m.P >= 0)
    assert np.allclose(m.P.sum(axis=1), 1.0)
    # strong persistence on the diagonal for an hourly AR(1) with phi=0.95
    assert np.all(np.diag(m.P) > 0.5)


def test_stationary_distribution_solves_pi_P():
    P = np.array([[0.9, 0.1, 0.0], [0.2, 0.7, 0.1], [0.0, 0.3, 0.7]])
    pi = mk.stationary(P)
    assert np.isclose(pi.sum(), 1.0)
    assert np.allclose(pi @ P, pi)
    # birth-death chain: detailed balance pi_i p_i,i+1 = pi_i+1 p_i+1,i
    assert np.isclose(pi[0] * 0.1, pi[1] * 0.2)
    assert np.isclose(pi[1] * 0.1, pi[2] * 0.3)


def test_stationary_matches_empirical_frequencies():
    x = ar1_cf(50000)
    m = mk.fit(x, 5)
    freq = np.bincount(mk.discretize(x, 5), minlength=5) / len(x)
    assert np.abs(mk.stationary(m.P) - freq).max() < 2e-3


def test_mean_levels_preserve_energy_centre_levels_do_not_have_to():
    x = ar1_cf()
    mean_model = mk.fit(x, 5, level_kind="mean")
    assert np.isclose(mean_model.mean(), x.mean(), atol=2e-3)
    two = mk.fit(x, 2, level_kind="centre")
    assert two.levels.tolist() == [0.0, 1.0]


def test_rows_without_data_fall_back_to_frequencies():
    counts = np.array([[5.0, 5.0, 0.0], [2.0, 8.0, 0.0], [0.0, 0.0, 0.0]])
    P = mk.normalise_rows(counts, fallback=np.array([0.5, 0.5, 0.0]))
    assert np.allclose(P.sum(axis=1), 1)
    assert np.allclose(P[2], [0.5, 0.5, 0.0])


def test_seasonal_fit_has_twelve_stochastic_matrices():
    x = ar1_cf(12 * 1000)
    months = np.repeat(np.arange(1, 13), 1000)
    m = mk.fit(x, 4, months=months)
    assert m.P.shape == (12, 4, 4)
    assert np.allclose(m.P.sum(axis=2), 1.0)
    assert m.levels.shape == (12, 4)


def test_simulated_chain_reproduces_transition_probabilities():
    P = np.array([[0.8, 0.2, 0.0], [0.1, 0.8, 0.1], [0.0, 0.3, 0.7]])
    model = mk.MarkovModel(P, np.array([0.0, 0.5, 1.0]), 3)
    rng = np.random.default_rng(0)
    u = rng.random((40, 5001))
    s = mk.simulate_states(model, u, np.ones(5000, dtype=int))
    counts = sum(mk.transition_counts(row, 3) for row in s)
    assert np.abs(mk.normalise_rows(counts) - P).max() < 0.01
    freq = np.bincount(s.ravel(), minlength=3) / s.size
    assert np.abs(freq - mk.stationary(P)).max() < 0.01


def test_inactive_hours_freeze_the_chain():
    P = np.array([[0.0, 1.0], [1.0, 0.0]])  # alternates on every active step
    model = mk.MarkovModel(P, np.array([0.0, 1.0]), 2)
    active = np.array([1, 0, 0, 1, 1, 0, 1], dtype=bool)
    u = np.full((1, 8), 0.5)
    s = mk.simulate_states(model, u, np.ones(7, dtype=int), active)[0]
    assert s[1] == s[0] and s[2] == s[0]
    assert s[3] != s[0]
    assert s[4] != s[3]
    assert s[5] == s[4] and s[6] != s[4]


def test_intensities_rows_sum_to_zero():
    x = ar1_cf(5000)
    lam = mk.intensities(mk.discretize(x, 5), 5)
    assert np.allclose(lam.sum(axis=1), 0.0)
    assert np.all(lam - np.diag(np.diag(lam)) >= 0)


@pytest.mark.parametrize("n", [2, 3, 7, 10])
def test_any_number_of_states(n):
    m = mk.fit(ar1_cf(8000), n)
    assert m.P.shape == (n, n)
    assert np.allclose(m.P.sum(axis=1), 1)
