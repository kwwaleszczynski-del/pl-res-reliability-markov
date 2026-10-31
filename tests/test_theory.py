import numpy as np
import pandas as pd
import pytest

from res_reliability import farms, theory as th


def test_catalogue_curve_of_the_mm92_is_the_published_one():
    f = farms.FARMS["kelmarsh"]
    assert th.specific_power(f) == pytest.approx(th.MM92_SPECIFIC_POWER)
    assert th.catalogue_cf(np.array([8.0]), f)[0] == pytest.approx(991.2 / 2055.0)
    assert th.catalogue_cf(np.array([2.0, 13.0, 24.9, 25.0]), f).tolist() == [0.0, 1.0, 1.0, 0.0]


def test_specific_power_scaling_gives_the_rotor_area_in_partial_load():
    small = farms.FARMS["penmanshiel"]  # MM82: same rating, smaller rotor
    big = farms.FARMS["kelmarsh"]
    ws = np.array([5.0, 6.0, 7.0])
    p_small = th.catalogue_cf(ws, small) * small.rated_kw
    k = (th.MM92_SPECIFIC_POWER / th.specific_power(small)) ** (1 / 3)
    # P(v) = 0.5 rho A Cp(kv) v^3: the MM82 at v delivers the MM92 power at k v times (82/92)^2 / k^3
    p_big = th.catalogue_cf(ws * k, big) * big.rated_kw
    assert np.allclose(p_small, p_big)
    assert np.all(th.catalogue_cf(ws, small) < th.catalogue_cf(ws, big))
    # the higher specific power of the SWT-2.3-82 lowers it further
    assert np.all(th.catalogue_cf(ws, farms.FARMS["hill_of_towie"]) < th.catalogue_cf(ws, small))


def test_hub_speed_interpolates_with_the_hourly_shear():
    ws10, ws100 = np.array([5.0, 6.0, 0.0]), np.array([7.0, 6.0, 4.0])
    v = farms.hub_speed(ws10, ws100, 59.0)
    alpha = np.log(7.0 / 5.0) / np.log(10.0)
    assert v[0] == pytest.approx(7.0 * 0.59**alpha)
    assert v[1] == pytest.approx(6.0)  # no shear
    assert 0 < v[2] < 4.0  # calm at 10 m: shear capped
    assert farms.hub_speed(ws10, ws100, 100.0) == pytest.approx(ws100)


def test_quantile_map_maps_distributions_and_keeps_order():
    rng = np.random.default_rng(1)
    src = rng.beta(0.8, 1.6, 20000)
    dst = 0.8 * src
    f = th.quantile_map(src, dst)
    x = np.linspace(0, 1, 11)
    assert np.all(np.diff(f(x)) >= 0)
    assert np.allclose(f(np.quantile(src, [0.1, 0.5, 0.9])), np.quantile(dst, [0.1, 0.5, 0.9]), atol=0.01)
    # ties at zero do not break it
    src2 = np.r_[np.zeros(5000), src]
    assert np.isfinite(th.quantile_map(src2, np.r_[np.zeros(5000), dst])(np.array([0.0, 0.3]))).all()


def synthetic_site(n=24 * 365 * 2, seed=2):
    rng = np.random.default_rng(seed)
    idx = pd.date_range("2018-01-01", periods=n, freq="1h", tz="UTC")
    z = np.zeros(n)
    e = rng.normal(0, 2.5 * np.sqrt(1 - 0.97**2), n)
    for t in range(1, n):  # hourly AR(1) around a seasonal mean
        z[t] = 0.97 * z[t - 1] + e[t]
    ws = np.clip(8 + 2 * np.cos(2 * np.pi * np.arange(n) / 8760) + z, 0, 30)
    theory = th.catalogue_cf(ws, farms.FARMS["kelmarsh"])
    op = np.clip(0.8 * theory + rng.normal(0, 0.03, n), 0, 1)
    df = pd.DataFrame({"ws_hub": ws, "theory": theory, "model": theory, "operation": op}, index=idx)
    df["valid"] = True
    df.iloc[:100, df.columns.get_loc("valid")] = False
    df["month"] = idx.month
    return df


def test_series_metrics_and_wind_bins():
    df = synthetic_site()
    m = th.series_metrics(df, "theory", "operation")
    assert m["hours"] == len(df) - 100
    assert m["ratio"] == pytest.approx(0.8, abs=0.01)
    assert set(m["ratio_by_year"]) == {2018, 2019}
    t = th.wind_bin_table(df)
    assert t["hours"].sum() == len(df) - 100
    row = t[t["ws_bin"] == 8.0].iloc[0]
    assert row["operation"] == pytest.approx(0.8 * row["theory"], abs=0.01)


def test_markov_structure_of_a_series():
    df = synthetic_site()
    s = th.markov_structure(df["operation"].to_numpy(), df["valid"].to_numpy(), df["month"].to_numpy())
    assert np.isclose(sum(s["occupancy"]), 1.0)
    assert s["model"].P.shape == (12, 5, 5)
    assert 0 < s["acf"][24] < s["acf"][1] <= 1


def test_merge_and_subtract_intervals():
    t = pd.Timestamp("2020-01-01", tz="UTC")
    h = lambda x: t + pd.Timedelta(hours=x)  # noqa: E731
    ev = pd.DataFrame({"turbine": ["T01"] * 3, "start": [h(0), h(1.2), h(5)], "end": [h(1), h(2), h(6)]})
    m = farms.merge_intervals(ev, gap_h=0.5)
    assert m["duration_h"].tolist() == [2.0, 1.0]
    win = pd.DataFrame({"turbine": ["T01"], "start": [h(0)], "end": [h(10)]})
    cut = farms.subtract_intervals(win, pd.DataFrame({"turbine": ["T01"], "start": [h(2)], "end": [h(3)]}))
    assert [(r.start, r.end) for r in cut.itertuples()] == [(h(0), h(2)), (h(3), h(10))]


def test_observation_windows_cut_long_silences():
    t = pd.Timestamp("2020-01-01", tz="UTC")
    stamps = [t + pd.Timedelta(days=d) for d in (0, 1, 2, 12, 13)]
    rec = pd.DataFrame({"turbine": "T01", "start": stamps, "end": [pd.NaT] * 5})
    w = farms.observation_windows(rec, t, t + pd.Timedelta(days=14), ["T01"])
    # the 10-day silence between day 2 and day 12 is not observed
    assert w["hours"].sum() == pytest.approx(24 * 4)
    assert w["start"].iloc[1] == t + pd.Timedelta(days=12)
