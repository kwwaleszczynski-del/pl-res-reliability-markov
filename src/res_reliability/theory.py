"""Theoretical production models against operating data at real plants.

Three hourly series per wind farm, all on the same hours:

* **theory**: ERA5 interpolated to hub height with the hourly shear between 10 and 100 m, through the
  manufacturer's catalogue power curve of the installed turbine at 1.225 kg/m3, no losses. Only the Senvion MM92
  (Kelmarsh) has a public catalogue curve (Open Energy Platform turbine library, as distributed with
  windpowerlib); for the MM82 and the SWT-2.3-82 the MM92 curve is moved to their specific power (rated power
  over rotor area): f(v) = f_MM92(v (SP_MM92 / SP)^(1/3)), which keeps the MM92 power coefficient and gives
  the rotor area of the turbine in partial load;
* **this model**: what this repository uses for Poland: ERA5 at 100 m, the MM92 curve measured at Kelmarsh and
  smoothed for a farm, 10 % wake and electrical losses (``power.wind_cf``);
* **operation**: the 10-minute SCADA of the farm, turbines in a forced outage left out (technical
  availability is the structural model).

ERA5 wind is instantaneous at the top of the hour and SCADA hours are means, so the ERA5 speed of an hour is
the mean of its two ends. PV: ERA5 GHI and temperature through the PV model of this repository with the
nameplate data of PVDAQ system 34 (``power.pv_output``) against its measured AC output.

Usage:
    df = wind_site(FARMS["kelmarsh"], outages)
    pv = pv_site()
    m = series_metrics(df, "theory", "operation")
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from res_reliability import farms, markov as mk, power, validation
from res_reliability.config import DATA_DIR, WIND_LOSSES

# Senvion MM92 catalogue power curve at 1.225 kg/m3 (kW, 1 m/s steps), Open Energy Platform wind turbine
# library as distributed with windpowerlib (oedb/power_curves.csv), originally windenergie-im-binnenland.de.
# A published sales curve, not a warranted one.
MM92_CAT_WS = np.arange(0.0, 26.0)
MM92_CAT_KW = np.array([0.0, 0.0, 0.0, 22.0, 93.1, 207.2, 388.2, 642.7, 991.2, 1355.7, 1780.0, 2012.6, 2048.0]
                       + [2055.0] * 13)
MM92_SPECIFIC_POWER = 2050e3 / (np.pi * 46.0**2)
CUT_OUT_MS = 25.0
WS_BINS = np.arange(0.0, 26.0, 1.0)
SEASONS = {"DJF": (12, 1, 2), "MAM": (3, 4, 5), "JJA": (6, 7, 8), "SON": (9, 10, 11)}
N_STATES = 5


def specific_power(farm: farms.Farm) -> float:
    return farm.rated_kw * 1e3 / (np.pi * (farm.rotor_m / 2) ** 2)


def catalogue_cf(ws: np.ndarray, farm: farms.Farm) -> np.ndarray:
    """Normalised catalogue power of the farm's turbine (MM92 catalogue moved to the turbine's specific power)."""
    ws = np.asarray(ws, dtype=float)
    k = (MM92_SPECIFIC_POWER / specific_power(farm)) ** (1 / 3)
    p = np.interp(ws * k, MM92_CAT_WS, MM92_CAT_KW / MM92_CAT_KW.max(), left=0.0, right=1.0)
    return np.where(ws >= CUT_OUT_MS, 0.0, p)


def mid_hour(x: pd.Series) -> pd.Series:
    """Mean of the values at the start and the end of each hour (instantaneous -> hour mean)."""
    return 0.5 * (x + x.shift(-1))


def wind_site(farm: farms.Farm, outages: pd.DataFrame, data_dir: Path = DATA_DIR) -> pd.DataFrame:
    """Hourly theory, this model and operation at a farm (UTC hours with SCADA inside the log window)."""
    e = farms.farm_era5(farm, data_dir)
    h = farms.cached_scada_hourly(farm, outages, data_dir)
    h = h[(h.index >= pd.Timestamp(farm.log_start, tz="UTC")) & (h.index < pd.Timestamp(farm.log_end, tz="UTC"))]
    ws_hub = mid_hour(e["ws_hub"]).reindex(h.index)
    ws100 = mid_hour(e["ws100"]).reindex(h.index)
    df = pd.DataFrame({
        "ws_hub": ws_hub, "ws100": ws100,
        "theory": catalogue_cf(ws_hub.to_numpy(), farm),
        "model": power.wind_cf(ws100.to_numpy(), WIND_LOSSES),
        "operation": h["cf_avail"], "operation_all": h["cf_all"],
    }, index=h.index)
    df["valid"] = df["operation"].notna() & df["ws_hub"].notna()
    df["month"] = df.index.month
    return df


def pv_site(data_dir: Path = DATA_DIR) -> pd.DataFrame:
    """Hourly theory (ERA5 + nameplate PV model, per unit of inverter capacity) and measured PVDAQ 34 output."""
    s = farms.PV_SITE
    h = farms.pvdaq_hourly()
    y = s["years"]
    h = h[(h.index >= pd.Timestamp(f"{y[0]}-01-01", tz="UTC") - pd.Timedelta(hours=s["utc_offset_h"])) &
          (h.index < pd.Timestamp(f"{y[-1] + 1}-01-01", tz="UTC") - pd.Timedelta(hours=s["utc_offset_h"]))]
    e = farms.pv_era5(data_dir)
    ghi = e["ghi"].shift(-1).reindex(h.index)  # hour mean labelled by its end -> by its start
    t2m = e["t2m"].reindex(h.index).interpolate(limit_direction="both")
    p, p_cs = power.pv_output(h.index, s["lat"], s["lon"], ghi.fillna(0.0).to_numpy(), t2m.to_numpy(),
                              tilt=s["tilt"], azimuth=s["azimuth"], dc_ac=s["dc_kwp"] / s["ac_kw"], gamma=s["gamma"])
    df = pd.DataFrame({"theory": p, "clear_sky": p_cs,
                       "operation": (h["ac_kw"] / s["ac_kw"]).clip(0.0, 1.0),
                       "operation_all": (h["ac_kw_all"] / s["ac_kw"]).clip(0.0, 1.0),
                       "outage": h["outage"]}, index=h.index)
    df.loc[df["outage"], "operation"] = np.nan
    df["daylight"] = df["clear_sky"] > 0.02
    df["valid"] = df["operation"].notna() & ghi.notna().to_numpy()
    df["local"] = df.index.tz_convert(None) + pd.Timedelta(hours=s["utc_offset_h"])
    df["month"] = df["local"].dt.month
    return df


# --- comparison -------------------------------------------------------------------------------

def series_metrics(df: pd.DataFrame, ref: str, col: str) -> dict:
    """Energy, seasonal and per-wind-bin capacity factors of ``col`` against ``ref`` on the valid hours."""
    v = df[df["valid"]]
    out = {"hours": int(len(v)), "cf": float(v[col].mean()), "ratio": float(v[col].sum() / v[ref].sum())}
    for name, months in SEASONS.items():
        s = v[v["month"].isin(months)]
        out[f"cf_{name}"] = float(s[col].mean())
        out[f"ratio_{name}"] = float(s[col].sum() / s[ref].sum())
    years = v.index.year
    out["ratio_by_year"] = {int(y): float(g[col].sum() / g[ref].sum()) for y, g in v.groupby(years)}
    out["corr_hourly"] = float(np.corrcoef(v[col], v[ref])[0, 1])
    return out


def wind_bin_table(df: pd.DataFrame, cols=("theory", "model", "operation")) -> pd.DataFrame:
    """Mean capacity factor per 1 m/s bin of the ERA5 hub-height speed, with hours per bin."""
    v = df[df["valid"]]
    b = pd.cut(v["ws_hub"], WS_BINS, right=False, labels=WS_BINS[:-1])
    g = v.groupby(b, observed=True)
    t = g[list(cols)].mean()
    t["hours"] = g.size()
    t.index = t.index.astype(float)
    t.index.name = "ws_bin"
    return t.reset_index()


def markov_structure(x: np.ndarray, valid: np.ndarray, months: np.ndarray, lags=(1, 6, 24)) -> dict:
    """Five-state monthly chain of a series on its valid hours: occupancy, residence times, ACF."""
    xv = np.where(valid, x, np.nan)
    model = mk.fit(x, N_STATES, months=months, valid=valid)
    states = mk.discretize(xv, N_STATES)
    occ = np.bincount(states[states >= 0], minlength=N_STATES) / max((states >= 0).sum(), 1)
    soj = validation.sojourns(states, N_STATES)
    acf = validation.acf(xv[None, :], max(lags))
    hom = mk.fit(x, N_STATES, valid=valid)
    return {"model": model, "homogeneous": hom, "states": states,
            "occupancy": occ.tolist(),
            "residence_h": [float(s.mean()) if len(s) else float("nan") for s in soj],
            "acf": {int(k): float(acf[k]) for k in lags},
            "levels": hom.levels.tolist()}


def transition_test(states_a: np.ndarray, P_b: np.ndarray) -> dict:
    """Anderson-Goodman test of the transitions of series a against the homogeneous matrix of series b."""
    counts = mk.transition_counts(states_a, N_STATES)
    return validation.anderson_goodman(counts, P_b)


def quantile_map(src: np.ndarray, dst: np.ndarray, n_q: int = 201):
    """Monotone map x -> F_dst^-1(F_src(x)) estimated on paired samples (empirical quantiles)."""
    q = np.linspace(0.0, 1.0, n_q)
    xs = np.quantile(np.asarray(src, dtype=float), q)
    ys = np.quantile(np.asarray(dst, dtype=float), q)
    # ties in src (many hours at zero or at the plateau): average the target quantiles of a tie
    ux, inv = np.unique(xs, return_inverse=True)
    uy = np.bincount(inv, weights=ys) / np.bincount(inv)
    return lambda x: np.interp(np.asarray(x, dtype=float), ux, uy)
