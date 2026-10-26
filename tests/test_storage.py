import numpy as np
import pytest

from res_reliability import montecarlo as mc, regimes as rgm, rts79, storage as sto

MONTH_HOURS = np.repeat(np.arange(1, 13), np.array([31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]) * 24)


def rare_deficits(n_sims=40, n_hours=2000, seed=0):
    """Margins with a daily cycle, plenty of surplus and a few deep dips of 1-6 h."""
    rng = np.random.default_rng(seed)
    t = np.arange(n_hours)
    m = 300 + 150 * np.sin(2 * np.pi * t / 24)[None, :] + rng.normal(0, 30, (n_sims, n_hours))
    for i in range(n_sims):
        for start in rng.choice(n_hours - 10, 3, replace=False):
            m[i, start:start + rng.integers(1, 7)] = -rng.uniform(0, 120)
    return m


def test_soc_balance_and_bounds():
    m = rare_deficits()
    r = sto.dispatch(m, 50.0, 150.0, efficiency=0.81, soc0=0.4, keep=True)
    soc, ch, dis = r["soc"][0], r["charge"][0], r["discharge"][0]
    e = np.sqrt(0.81)
    prev = np.c_[np.full(len(m), 0.4 * 150.0), soc[:, :-1]]
    assert np.allclose(soc, prev + e * ch - dis / e)  # hour by hour
    assert np.allclose(soc[:, -1], 0.4 * 150.0 + e * ch.sum(axis=1) - dis.sum(axis=1) / e)  # whole year
    assert soc.min() >= 0.0 and soc.max() <= 150.0 + 1e-9  # never below empty, never above full
    assert np.all(dis <= prev * e + 1e-9)  # cannot deliver more than it holds
    assert ch.max() <= 50.0 + 1e-9 and dis.max() <= 50.0 + 1e-9
    assert np.all(ch[m < 0] == 0) and np.all(dis[m > 0] == 0)  # no charging in a shortfall
    assert np.all(r["net"][0][m >= 0] >= -1e-9)  # charging never creates a shortfall


def test_single_event_is_covered_until_the_energy_runs_out():
    m = np.full((1, 100), 200.0)
    m[0, 50:53] = -50.0  # 3 h at 50 MW short, the storage holds 2 h of it
    r = sto.dispatch(m, 100.0, 100.0)
    assert r["lole"][0, 0] == 1 and r["eens"][0, 0] == pytest.approx(50.0)
    assert r["lolf"][0, 0] == 1
    r = sto.dispatch(m, 30.0, 1000.0)  # power-limited: 20 MW short in each of the 3 hours
    assert r["lole"][0, 0] == 3 and r["eens"][0, 0] == pytest.approx(60.0)
    r = sto.dispatch(m, 100.0, 100.0, soc0=0.0)  # starts empty but refills from the surplus before the event
    assert r["lole"][0, 0] == 1


def test_skipping_quiet_hours_is_exact():
    m = rare_deficits(seed=3)
    P = np.array([10.0, 40.0, 80.0, 80.0])
    for kw in ({}, {"efficiency": 0.85, "soc0": 0.5}, {"efficiency": 0.85, "soc0": 0.2, "warmup_h": 300}):
        for shift in (0.0, 35.0):
            fast = sto.dispatch(m, P, P * np.array([1, 2, 4, 8]), shift=shift, **kw)
            full = sto.dispatch(m, P, P * np.array([1, 2, 4, 8]), shift=shift, keep=True, **kw)
            for k in ("lole", "eens", "lolf", "charged", "discharged"):
                assert np.allclose(fast[k], full[k]), (kw, shift, k)


def test_colocated_storage_charges_only_from_res():
    m = rare_deficits(seed=5)
    res = np.full_like(m, 5.0)
    r = sto.dispatch(m, 50.0, 200.0, soc0=0.0, charge_cap=res, keep=True)
    assert r["charge"].max() <= 5.0 + 1e-9
    assert sto.dispatch(m, 50.0, 200.0, soc0=0.0)["lole"].sum() <= r["lole"].sum()


def test_storage_with_ample_energy_is_firm_capacity():
    """Unlimited energy and plenty of surplus between the dips: ELCC equals the power."""
    m = rare_deficits(seed=7)
    e = sto.elcc(m, np.array([40.0]), np.array([1e7]), iters=12)
    assert e[0] == pytest.approx(40.0, abs=0.05)
    # evening deficits of several hours every day: a longer duration earns more capacity credit
    t = np.arange(24 * 60)
    daily = np.tile(50 + 100 * np.cos(2 * np.pi * t / 24), (3, 1))
    e = sto.elcc(daily, np.full(3, 40.0), np.array([40.0, 160.0, 320.0]), iters=12)
    assert e[0] < e[1] < e[2] <= 40.0


def test_energy_needed_for_a_known_deficit():
    m = np.full((5, 300), 100.0)
    m[:, 150:160] = -20.0  # 200 MWh short in 10 h every year
    e = sto.energy_needed(m, 0.0, np.array([np.inf, 4.0, np.inf]), np.array([1.0, 1.0, 0.81]), iters=24)
    assert e[0] == pytest.approx(200.0, rel=1e-3)
    assert e[1] == pytest.approx(200.0, rel=1e-3)  # 4 h: power 50 MW is enough, energy binds
    assert e[2] == pytest.approx(200.0 / 0.9, rel=1e-3)  # discharge losses
    assert sto.energy_needed(m, 10.0, np.array([np.inf]))[0] == 0.0  # target already met


def rts_spec(**kw):
    u = rts79.units()
    return mc.SystemSpec(load=rts79.load_profile()[:8760], conv_capacity=np.array([a for a, _, _ in u]),
                         conv_mttf=np.array([b for _, b, _ in u]), conv_mttr=np.array([c for _, _, c in u]), **kw)


def test_storage_inside_monte_carlo_matches_dispatch_on_margins():
    st = sto.Storage(100.0, 400.0, efficiency=0.85, soc0=0.5)
    a = mc.run(rts_spec(), 40, seed=2, keep_margin=True)
    b = mc.run(rts_spec(storage=st), 40, seed=2)
    r = sto.dispatch(a.margin, 100.0, 400.0, 0.85, 0.5)
    assert np.array_equal(b.lole, r["lole"][0])
    assert b.lole.mean() < a.lole.mean()


def synthetic_regime_weather(n_days=365):
    rng = np.random.default_rng(1)
    n_years = 4
    wind = np.clip(np.where(rng.random((n_years, n_days)) < 0.3, 0.05, 0.5).repeat(24, axis=1)
                   + rng.normal(0, 0.03, (n_years, n_days * 24)), 0, 1)
    k = np.clip(rng.uniform(0.2, 1.0, (n_years, n_days * 24)), 0, 1)
    day = np.tile((np.arange(24) >= 7) & (np.arange(24) < 17), (n_years, n_days))
    return rgm.fit_weather(wind, k, day, MONTH_HOURS, quantiles=(0.3,)), day


def test_same_seed_same_years_with_regimes_and_storage():
    w, day = synthetic_regime_weather()
    env = day * 0.8
    res = [mc.ResPlant("wind", 300.0, mc.RegimeProduction(w, "wind")),
           mc.ResPlant("pv", 300.0, mc.RegimeProduction(w, "pv", envelope=env, active=day))]
    spec = rts_spec(res=res, weather=w, n_weather_years=4, storage=sto.Storage(100.0, 200.0, 0.85, 0.5))
    a = mc.run(spec, 24, seed=9, chunk=24, keep_margin=True)
    b = mc.run(spec, 24, seed=9, chunk=5, keep_margin=True)
    c = mc.run(spec, 24, seed=10, keep_margin=True)
    assert np.array_equal(a.margin, b.margin) and np.array_equal(a.lole, b.lole)
    assert not np.array_equal(a.margin, c.margin)
    # the shared regime path does not depend on how many plants there are
    wind_only = mc.run(rts_spec(res=res[:1], weather=w, n_weather_years=4), 6, seed=9, keep_res=True)
    zero_pv = [res[0], mc.ResPlant("pv", 0.0, res[1].production)]
    both = mc.run(rts_spec(res=zero_pv, weather=w, n_weather_years=4), 6, seed=9, keep_res=True)
    assert np.array_equal(wind_only.res, both.res)
    paths = [w.regime_paths([mc.weather_stream(9, i)], np.array([i % 4])) for i in range(3)]
    again = [w.regime_paths([mc.weather_stream(9, i)], np.array([i % 4])) for i in range(3)]
    assert all(np.array_equal(p, q) for p, q in zip(paths, again))


def test_weather_stream_is_a_fixed_child_of_the_year_seed():
    a = mc.weather_stream(5, 3).random(4)
    b = np.random.default_rng(np.random.SeedSequence(5, spawn_key=(3,)).spawn(mc.WEATHER_STREAM + 1)[-1]).random(4)
    assert np.array_equal(a, b)
    assert not np.array_equal(a, mc.year_streams(5, 3, 2)[1].random(4))


def test_pv_regime_output_is_zero_at_night_and_bounded():
    w, day = synthetic_regime_weather()
    prod = mc.RegimeProduction(w, "pv", envelope=day * 0.8, active=day)
    years = np.array([0, 1, 2, 3])
    regs = w.regime_paths([np.random.default_rng(i) for i in range(4)], years)
    x = prod.sample([np.random.default_rng(10 + i) for i in range(4)], years, 8760, regs)
    assert np.all(x[~day] == 0.0)
    assert x.max() <= 0.8 + 1e-9 and x.min() >= 0.0
