import numpy as np
import pandas as pd
import pytest

from res_reliability import power, pse, structural, validation, weather
from res_reliability import markov as mk
from res_reliability.geo import GRID, fetch_order, poland_outline


def legacy_csv(rows):
    head = "Data;Godzina;Krajowe zapotrzebowanie na moc;Sumaryczna generacja JWCD\n"
    return head + "".join(f"{d};{h};{v};1\n" for d, h, v in rows)


def test_legacy_pse_autumn_dst_day_has_25_hours():
    hours = ["1", "2", "2A", *[str(h) for h in range(3, 25)]]
    rows = [("2015-10-25", h, f"{1000 + i},5") for i, h in enumerate(hours)]
    s = pse.parse_legacy_csv(legacy_csv(rows))
    assert len(s) == 25
    assert str(s.index[0]) == "2015-10-24 22:00:00+00:00"  # 00:00 CEST
    assert str(s.index[-1]) == "2015-10-25 22:00:00+00:00"  # 23:00 CET
    assert s.iloc[2] == pytest.approx(1002.5)
    assert s.index.is_monotonic_increasing and s.index.is_unique


def test_legacy_pse_spring_dst_day_has_23_hours():
    hours = ["1", *[str(h) for h in range(3, 25)]]
    rows = [("2015-03-29", h, "1000,0") for h in hours] + [("2015-03-30", "1", "999,0")]
    s = pse.parse_legacy_csv(legacy_csv(rows))
    assert len(s) == 24
    assert str(s.index[0]) == "2015-03-28 23:00:00+00:00"
    assert str(s.index[22]) == "2015-03-29 21:00:00+00:00"
    assert str(s.index[23]) == "2015-03-29 22:00:00+00:00"  # next day 00:00 CEST


def test_legacy_pse_missing_hour_does_not_shift_the_day():
    rows = [("2018-02-22", str(h), f"{h},0") for h in range(1, 25) if h not in (5, 6)]
    s = pse.parse_legacy_csv(legacy_csv(rows))
    assert len(s) == 22
    # hour 7 (06:00-07:00 CET) starts at 05:00 UTC whatever is missing before it
    assert s[pd.Timestamp("2018-02-22 05:00", tz="UTC")] == pytest.approx(7.0)
    assert pd.Timestamp("2018-02-22 03:00", tz="UTC") not in s.index


def test_api_quarter_hours_to_hourly_mean():
    rows = [{"dtime_utc": f"2024-06-14 0{h}:{m:02d}:00", "load_actual": v}
            for h, m, v in [(0, 15, 10.0), (0, 30, 20.0), (0, 45, 30.0), (1, 0, 40.0), (1, 15, 50.0)]]
    s = pse.parse_api_values(rows)
    assert str(s.index[0]) == "2024-06-14 00:00:00+00:00"
    assert s.iloc[0] == pytest.approx(25.0)
    assert s.iloc[1] == pytest.approx(50.0)


def test_splice_factor_recovers_level_shift():
    idx = pd.date_range("2024-05-20", "2024-07-10", freq="1h", tz="UTC")
    local_hour = idx.tz_convert("Europe/Warsaw").hour.to_numpy()
    base = pd.Series(15000 + 3000 * np.sin(2 * np.pi * local_hour / 24), index=idx)
    switch = pse.API_START.tz_localize("Europe/Warsaw")
    legacy, api = base[base.index < switch], base[base.index >= switch] * 0.96
    assert pse.splice_factor(legacy, api) == pytest.approx(1 / 0.96, rel=1e-6)


def test_year_matrix_shape_and_gap_rule():
    idx = pd.date_range("2016-01-01", "2017-12-31 23:00", freq="1h", tz="UTC")
    s = pd.Series(np.arange(len(idx), dtype=float), index=idx)
    m = pse.year_matrix(s, [2016, 2017])
    assert m.shape == (2, 8760)
    gappy = s.drop(s.index[100:110])  # 10 h: too long to interpolate
    with pytest.raises(ValueError):
        pse.year_matrix(gappy, [2016])
    small_gap = s.drop(s.index[100:104])
    assert np.allclose(pse.year_matrix(small_gap, [2016])[0, 100:104], [100, 101, 102, 103])


def test_weather_weight_and_parse():
    assert weather.request_weight(3, "2015-01-01", "2024-12-31") == pytest.approx(3 / 10 * 3653 / 14)
    assert weather.request_weight(1, "2020-01-01", "2020-01-02") == 1.0
    payload = {"hourly": {"time": [1420070400, 1420074000], "wind_speed_100m": [5.0, 6.0],
                          "shortwave_radiation": [0.0, 10.0], "temperature_2m": [-1.0, -2.0]}}
    df = weather.parse_response(payload)
    assert list(df.columns) == ["ws100", "ghi", "t2m"]
    assert str(df.index[0]) == "2015-01-01 00:00:00+00:00"


def test_grid_inside_outline_and_fetch_order_is_a_permutation():
    from matplotlib.path import Path
    path = Path(poland_outline())
    assert all(path.contains_point((p.lon, p.lat)) for p in GRID)
    assert sorted(p.id for p in fetch_order()) == sorted(p.id for p in GRID)
    assert 20 <= len(GRID) <= 40


def test_wind_power_curve_shape():
    v = np.array([0.0, 2.0, 3.5, 8.0, 13.0, 20.0, 25.0, 30.0])
    p = power.turbine_curve(v)
    assert p[0] == 0 and p[1] == 0 and 0 < p[2] < p[3] < p[4] <= 1
    assert p[5] == pytest.approx(1.0) and p[6] == 0.0 and p[7] == 0.0
    smooth = power.farm_curve(np.linspace(0, 30, 301))
    assert np.all(smooth >= 0) and np.all(smooth <= 1)
    assert power.wind_cf(np.array([13.0]), losses=0.1)[0] == pytest.approx(0.9 * power.farm_curve(np.array([13.0]))[0])


def test_pv_output_zero_at_night_and_bounded():
    idx = pd.date_range("2021-06-21", periods=48, freq="1h", tz="UTC")
    hour = idx.hour.to_numpy()
    ghi = np.where((hour >= 4) & (hour <= 18), 600.0, 0.0)
    p, p_cs = power.pv_output(idx, 52.0, 19.0, ghi, np.full(48, 20.0))
    assert np.all((p >= 0) & (p <= 1)) and np.all((p_cs >= 0) & (p_cs <= 1))
    assert p[0] == 0 and p_cs[0] == 0  # 00:00 UTC is night in Poland
    assert p_cs[10] > 0.6  # 10:00 UTC around noon, clear sky
    k, day = power.clear_sky_index(p, p_cs)
    assert np.all(np.isnan(k[~day]))


def test_merge_events_and_stats():
    t = pd.Timestamp("2020-01-01", tz="UTC")
    ev = pd.DataFrame({"turbine": [1, 1, 1, 2],
                       "start": [t, t + pd.Timedelta(hours=1, minutes=10), t + pd.Timedelta(hours=10), t],
                       "end": [t + pd.Timedelta(hours=1), t + pd.Timedelta(hours=2), t + pd.Timedelta(hours=13),
                               t + pd.Timedelta(hours=4)]})
    m = structural.merge_events(ev, gap_h=0.5)
    assert len(m) == 3
    assert m["duration_h"].tolist() == pytest.approx([2.0, 3.0, 4.0])
    st = structural.event_stats(m, exposure_h=1000.0, n_units=2)
    assert st["for"] == pytest.approx(9 / 1000)
    assert st["mttr_h"] == pytest.approx(3.0)
    assert st["mttf_h"] == pytest.approx((1000 - 9) / 3)
    assert structural.event_stats(m, 1000.0, 2, min_duration_h=3.0)["events"] == 2


def test_sojourns_and_acf():
    s = np.array([0, 0, 1, 1, 1, 0, 2, 2, 0])
    soj = validation.sojourns(s, 3)
    assert soj[1].tolist() == [3] and soj[0].tolist() == [1] and soj[2].tolist() == [2]
    x = np.sin(np.linspace(0, 20 * np.pi, 2000))
    a = validation.acf(x, 5)
    assert a[0] == pytest.approx(1.0) and a[1] > 0.99


def test_anderson_goodman_accepts_same_chain_rejects_different():
    P = np.array([[0.9, 0.1, 0.0], [0.1, 0.8, 0.1], [0.0, 0.2, 0.8]])
    model = mk.MarkovModel(P, np.array([0.0, 0.5, 1.0]), 3)
    rng = np.random.default_rng(4)
    s = mk.simulate_states(model, rng.random((1, 20001)), np.ones(20000, dtype=int))[0]
    same = validation.anderson_goodman(mk.transition_counts(s, 3), P)
    assert same["p_value"] > 0.01
    other = np.array([[0.8, 0.2, 0.0], [0.1, 0.8, 0.1], [0.0, 0.2, 0.8]])
    diff = validation.anderson_goodman(mk.transition_counts(s, 3), other)
    assert diff["p_value"] < 1e-6


def test_cross_validation_prefers_markov_over_iid():
    rng = np.random.default_rng(9)
    P = np.array([[0.95, 0.05, 0.0], [0.05, 0.9, 0.05], [0.0, 0.05, 0.95]])
    model = mk.MarkovModel(P, np.array([0.0, 0.5, 1.0]), 3)
    s = mk.simulate_states(model, rng.random((4, 3001)), np.ones(3000, dtype=int))
    x = np.array([0.0, 0.5, 1.0])[s]
    r = validation.cross_validate(x, np.ones(3000, dtype=int), np.arange(3), np.array([3]), n_states=3)
    assert r["ll_homogeneous"] > r["ll_iid"]
    assert r["ll_homogeneous"] == pytest.approx(r["ll_test_fit"], abs=0.02)
