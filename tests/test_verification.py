import numpy as np
import pytest

from res_reliability import adequacy, elcc, regimes as rgm, storage as sto, verification as ver

MONTH_HOURS = np.repeat(np.arange(1, 13), np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]) * 24)


def test_year_weights_resample_whole_years():
    years = np.repeat(np.arange(5), 4)
    w = ver.year_weights(years, 200, np.random.default_rng(0), within=False)
    assert w.shape == (200, 20) and np.allclose(w.sum(axis=1), 1.0)
    # without the second stage every simulated year of a chosen weather year gets the same weight
    for row in w:
        assert np.allclose(row.reshape(5, 4), row.reshape(5, 4)[:, :1])
    w2 = ver.year_weights(years, 200, np.random.default_rng(0), within=True)
    assert np.allclose(w2.sum(axis=1), 1.0)
    assert np.allclose(w2.mean(axis=0), 1 / 20, atol=0.01)  # unbiased on average


def test_crossing_interpolates_both_ways():
    grid = np.array([0.0, 10.0, 20.0, 30.0])
    up = np.array([[1.0, 2.0, 4.0, 8.0]])
    assert ver.crossing(grid, up, np.array([3.0]))[0] == pytest.approx(15.0)
    assert ver.crossing(grid, up, np.array([0.5]))[0] == 0.0
    assert ver.crossing(grid, up, np.array([9.0]))[0] == 30.0
    e = np.array([1.0, 10.0, 100.0])
    down = np.array([[5.0, 3.0, 1.0]])
    assert ver.crossing(e, down, np.array([2.0]), increasing=False, log=True)[0] == pytest.approx(np.sqrt(10) * 10)


def test_crossing_recovers_storage_elcc_from_curves():
    rng = np.random.default_rng(1)
    t = np.arange(24 * 40)
    m = 60 + 100 * np.cos(2 * np.pi * t / 24)[None, :] + rng.normal(0, 20, (30, len(t)))
    shifts = np.linspace(0, 40, 161)
    curves = ver.storage_lole_curves(m, 40.0, 120.0, shifts)
    target = (m < 0).sum(axis=1).mean()
    assert curves.shape == (len(shifts), len(m))
    from_curves = ver.crossing(shifts, curves.mean(axis=1)[None], np.array([target]))[0]
    bisect = sto.elcc(m, np.array([40.0]), np.array([120.0]), iters=14)[0]
    assert from_curves == pytest.approx(bisect, abs=0.5)


def test_lole_by_year_averages_to_the_copt_expectation():
    copt = adequacy.Copt.from_units([100, 100, 50], [0.05, 0.05, 0.1])
    rng = np.random.default_rng(2)
    load = rng.uniform(80, 200, (3, 500))
    d = elcc.historical(rng.uniform(0, 1, (3, 500)), 60.0)
    by_year = ver.lole_by_year(copt, load, d, shift=5.0)
    assert by_year.shape == (3,)
    assert by_year.mean() == pytest.approx(elcc.lole(copt, load, d, 5.0))
    assert ver.lole_by_year(copt, load, None).mean() == pytest.approx(elcc.lole(copt, load))


def test_energy_curves_fall_with_energy():
    m = np.full((4, 300), 100.0)
    m[:, 100:110] = -20.0
    curves = ver.energy_lole_curves(m, np.array([50.0, 150.0, 250.0]))
    assert curves[:, 0].tolist() == [8.0, 3.0, 0.0]  # 200 MWh short in 10 h of 20 MW
    assert np.all(np.diff(curves, axis=0) <= 0)


def test_block_test_and_pit():
    rng = np.random.default_rng(3)
    blocks = rng.normal(0, 1, 5000)
    mid = ver.block_test(0.0, blocks)
    assert mid["percentile"] == pytest.approx(50, abs=2) and mid["p_value"] > 0.9
    far = ver.block_test(3.5, blocks)
    assert far["percentile"] > 99.9 and far["p_value"] < 0.01
    assert ver.pit(0.0, np.array([-1.0, 0.0, 1.0, 2.0])) == pytest.approx(0.375)
    # calibrated predictions give uniform PIT values
    pits = [ver.pit(rng.normal(), rng.normal(size=400)) for _ in range(300)]
    assert np.mean(pits) == pytest.approx(0.5, abs=0.05)


def test_regime_model_applies_its_temperature_classes_to_other_years():
    rng = np.random.default_rng(4)
    wind = np.clip(rng.beta(1.2, 3, (4, 8760)), 0, 1)
    k = rng.uniform(0.2, 1.0, (4, 8760))
    day = np.tile((np.arange(24) >= 7) & (np.arange(24) < 17), (4, 365))
    temp = rng.normal(0, 4, (4, 8760)) + 10 * np.sin(2 * np.pi * np.arange(8760) / 8760)
    w = rgm.fit_weather(wind, k, day, MONTH_HOURS, temperature=temp)
    assert np.array_equal(w.covariate_classes(temp), w.cov)  # same classes for the fitting years
    held = w.for_years(temp[[2]] - 30.0)  # a much colder year: every day in the coldest third
    assert held.cov.shape == (1, 365) and np.all(held.cov == 0)
    assert held.process is w.process and w.cov.shape == (4, 365)
