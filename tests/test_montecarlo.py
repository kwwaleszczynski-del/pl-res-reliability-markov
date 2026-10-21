import numpy as np
import pytest

from res_reliability import adequacy, elcc, markov as mk, montecarlo as mc, rts79


def rts_spec(**kw):
    u = rts79.units()
    return mc.SystemSpec(
        load=rts79.load_profile(),
        conv_capacity=np.array([a for a, _, _ in u]),
        conv_mttf=np.array([b for _, b, _ in u]),
        conv_mttr=np.array([c for _, _, c in u]),
        **kw,
    )


def test_two_state_marginal_availability():
    rng = np.random.default_rng(1)
    cap = np.ones(400)
    avail = mc.two_state_availability(rng, cap, np.full(400, 90.0), np.full(400, 10.0), 4000)
    assert avail.mean() / 400 == pytest.approx(0.9, abs=0.01)
    assert avail.min() >= 0 and avail.max() <= 400


def test_two_state_mean_down_time():
    rng = np.random.default_rng(2)
    avail = mc.two_state_availability(rng, np.ones(1), np.array([200.0]), np.array([20.0]), 400000)
    down = avail < 0.5
    starts = np.flatnonzero(down[1:] & ~down[:-1])
    ends = np.flatnonzero(~down[1:] & down[:-1])
    if ends[0] < starts[0]:
        ends = ends[1:]
    n = min(len(starts), len(ends))
    # sampled at hour starts: mean of ceil-discretised exponential durations ~ MTTR (+/- discretisation)
    assert np.mean(ends[:n] - starts[:n]) == pytest.approx(20.0, rel=0.08)


def test_same_seed_same_result_and_chunk_independence():
    spec = rts_spec()
    a = mc.run(spec, 60, seed=7, chunk=60)
    b = mc.run(spec, 60, seed=7, chunk=7)
    c = mc.run(spec, 60, seed=8)
    assert np.array_equal(a.lole, b.lole)
    assert np.array_equal(a.eens, b.eens)
    assert not np.array_equal(a.lole, c.lole)


def test_simulation_i_does_not_depend_on_number_of_simulations():
    spec = rts_spec()
    a = mc.run(spec, 30, seed=3)
    b = mc.run(spec, 10, seed=3)
    assert np.array_equal(a.lole[:10], b.lole)


def test_monte_carlo_matches_analytical_constant_load():
    """Single 100 MW unit, FOR 0.1, load 60 MW: LOLE = 0.1 * 8760 exactly in expectation."""
    spec = mc.SystemSpec(load=np.full(8760, 60.0), conv_capacity=np.array([100.0]),
                         conv_mttf=np.array([90.0]), conv_mttr=np.array([10.0]))
    r = mc.run(spec, 400, seed=11)
    s = r.summary()
    assert abs(s["lole"] - 876.0) < 3 * s["lole_std"] / np.sqrt(400)
    # one event per up->down cycle of 100 h, minus outages (or repairs) that fall between two
    # hourly samples: P(miss) = E[max(0, 1 - D)] = 1 - m (1 - exp(-1/m)) for D ~ Exp(mean m)
    miss_down = 1 - 10.0 * (1 - np.exp(-1 / 10.0))
    miss_up = 1 - 90.0 * (1 - np.exp(-1 / 90.0))
    assert s["lolf"] == pytest.approx(87.6 * (1 - miss_down) * (1 - miss_up), rel=0.03)


def test_monte_carlo_rts_within_ci_of_analytical():
    r = mc.run(rts_spec(), 1500, seed=0).summary()
    assert abs(r["lole"] - 9.394) < 2.6 * r["lole_std"] / np.sqrt(1500)


def test_markov_res_mc_matches_analytical_expectation():
    P = np.array([[0.9, 0.1, 0.0], [0.05, 0.9, 0.05], [0.0, 0.1, 0.9]])
    model = mk.MarkovModel(P, np.array([0.0, 0.5, 1.0]), 3)
    months = np.ones(2000, dtype=int)
    load = np.full(2000, 150.0)
    spec = mc.SystemSpec(load=load, conv_capacity=np.array([100.0, 100.0]),
                         conv_mttf=np.array([95.0, 95.0]), conv_mttr=np.array([5.0, 5.0]),
                         res=[mc.ResPlant("res", 100.0, mc.MarkovProduction(model, months))])
    r = mc.run(spec, 600, seed=5).summary()
    copt = adequacy.Copt.from_units([100, 100], [0.05, 0.05])
    dist = elcc.from_markov(model, months, 100.0)
    ana = elcc.indices(copt, load, dist)
    assert abs(r["lole"] - ana["lole"]) < 3 * r["lole_std"] / np.sqrt(600)


def test_mc_elcc_of_firm_capacity():
    base = mc.run(rts_spec(), 200, seed=4, keep_margin=True)
    firm = mc.HistoricalProduction(np.ones((1, 8736)))
    withres = mc.run(rts_spec(res=[mc.ResPlant("firm", 150.0, firm)]), 200, seed=4, keep_margin=True)
    assert elcc.mc_elcc(base.margin, withres.margin, 150.0) == pytest.approx(150.0, abs=0.2)


def test_structural_modes():
    rng = np.random.default_rng(0)
    plant = mc.Structural(0.05, 50.0, mode="plant_bernoulli").sample(rng, 20000)
    assert set(np.unique(plant)) <= {0.0, 1.0}
    assert 1 - plant.mean() == pytest.approx(0.05, abs=0.01)
    units = mc.Structural(0.05, 50.0, n_units=200).sample(rng, 20000)
    assert units.mean() == pytest.approx(0.95, abs=0.01)
    assert units.std() < 0.05
    assert np.all(mc.Structural(0.0, 50.0).sample(rng, 10) == 1.0)
    farms = mc.Structural(0.05, 50.0, n_units=29, mode="hourly_bernoulli").sample(rng, 20000)
    assert farms.mean() == pytest.approx(0.95, abs=0.005)
    assert np.allclose(farms * 29, np.round(farms * 29))
    # hourly draws are independent: no persistence, unlike the two-state process
    assert abs(np.corrcoef(farms[:-1], farms[1:])[0, 1]) < 0.05
    assert np.corrcoef(units[:-1], units[1:])[0, 1] > 0.9
