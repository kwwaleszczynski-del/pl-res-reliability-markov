"""Open operating data of real plants: three UK wind farms and one US PV system, with ERA5 at each site.

=============  ==============================  ======================  ===================================
plant          turbines / inverter             10-minute SCADA used    event log used
=============  ==============================  ======================  ===================================
Kelmarsh       6 x Senvion MM92, 2.05 MW       2016-04-15 .. 2021      Greenbyte status log 2016-2021
Penmanshiel    14 x Senvion MM82, 2.05 MW      2018 and 2020           Greenbyte status log 2016-2024
Hill of Towie  21 x Siemens SWT-2.3-VS-82      2021                    controller alarm log 2016-04/2026
PVDAQ 34       146.6 kWp, one 135 kW inverter  15 min, 2016-2018       (inverter outages from the data)
=============  ==============================  ======================  ===================================

The archives (Zenodo, OEDI; CC BY 4.0) are downloaded and parsed by sibling projects; this module reads their
caches, so nothing large is downloaded here:

* ``KELMARSH_DATA_DIR`` (default ``../wind-power-curve-kelmarsh/data``): 10-minute SCADA in
  ``processed/scada_<year>.parquet``, status logs in the yearly SCADA zips;
* ``WIND_DOWNTIME_DATA_DIR`` (``../wind-downtime-availability/data``): ``processed/penmanshiel_<year>.parquet``,
  ``processed/hill_of_towie_2021.parquet``;
* ``WIND_OM_DATA_DIR`` (``../wind-om-reliability/data``): ``raw/penmanshiel/status_*.zip`` and
  ``raw/hill_of_towie/alarms_*.zip``;
* ``PVDAQ_LONG_DIR`` (``../pv-budget-bridge/data/pvdaq``): ``long_34.parquet`` (AC power in kW, local standard time).

Forced outages: the ``Stop`` records the operator files under the IEC 61400-26 category "Forced outage"
(Kelmarsh, Penmanshiel), and for Hill of Towie the alarm codes that ``wind-om-reliability`` classed as forced
outages from the 2021 SCADA (stopped and turbine-error timer on). Stops of one turbine within ``MERGE_GAP_H`` of
each other form one outage. A turbine is observed while its log is alive: silences of more than a week (missing
files) and communication losses of 6 h or more are cut out of the exposure.

ERA5 at each site comes from the Open-Meteo archive (one request per site, cached in ``data/farms``).

Usage:
    python -m res_reliability fetch-farms
    w = farm_era5(FARMS["kelmarsh"])           # hourly ws10, ws100 (UTC)
    ev, win = forced_outages(FARMS["penmanshiel"])
    h = scada_hourly(FARMS["kelmarsh"], ev)    # hourly farm CF with and without the forced outages
"""

from __future__ import annotations

import io
import json
import logging
import os
import re
import time
import zipfile
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd
import requests

from res_reliability import weather
from res_reliability.config import DATA_DIR, kelmarsh_dir

log = logging.getLogger(__name__)

MERGE_GAP_H = 0.5
MAX_SILENCE = pd.Timedelta(days=7)
MIN_COMM_GAP = pd.Timedelta(hours=6)

# Hill of Towie alarm codes classed as forced outages by wind-om-reliability (mappings/siemens_swt_hot.csv):
# stopped in at least 74 % of the 2021 SCADA intervals they cover, with the turbine-error timer running.
HOT_FORCED_CODES = frozenset({
    3123, 3125, 3126, 3171, 3172, 3173, 3219, 3220, 3223, 3411, 3412, 3413, 4111, 4112, 6102, 6160, 9303, 9306,
    10118, 10128, 13110, 13172, 13207, 13224, 13225, 13309, 13311, 13325, 13900, 13902, 14038, 14042, 14302,
    15006, 15033, 15037, 15038, 15042, 15214, 63005, 63016, 63025, 63026, 63027, 63036, 63038, 63207, 63304, 64066,
})
HOT_STATION_OFFSET = 2304509  # StationNr = turbine number + offset


@dataclass(frozen=True)
class Farm:
    key: str
    label: str
    model: str
    n_turbines: int
    rated_kw: float
    rotor_m: float
    hub_m: float  # mean hub height
    lat: float
    lon: float
    log_start: str  # observation window of the event log (UTC)
    log_end: str
    scada_years: tuple[int, ...]
    era5_end: str


FARMS = {
    "kelmarsh": Farm("kelmarsh", "Kelmarsh", "Senvion MM92", 6, 2050.0, 92.0, (4 * 78.5 + 2 * 68.5) / 6,
                     52.4015, -0.9431, "2016-04-15", "2022-01-01", (2016, 2017, 2018, 2019, 2020, 2021),
                     "2021-12-31"),
    "penmanshiel": Farm("penmanshiel", "Penmanshiel", "Senvion MM82", 14, 2050.0, 82.0, 59.0,
                        55.9036, -2.2913, "2016-09-01", "2025-01-01", (2018, 2020), "2024-12-31"),
    "hill_of_towie": Farm("hill_of_towie", "Hill of Towie", "Siemens SWT-2.3-VS-82", 21, 2300.0, 82.0, 59.0,
                          57.5056, -3.0663, "2016-01-02", "2026-04-01", (2021,), "2026-03-31"),
}

PV_SITE = {"key": "pvdaq34", "label": "PVDAQ 34 (Las Vegas)", "lat": 36.1952, "lon": -115.1582, "tilt": 11.2,
           "azimuth": 180.0, "dc_kwp": 146.64, "ac_kw": 135.0, "gamma": -0.0045, "utc_offset_h": -8,
           "years": (2016, 2017, 2018)}


def _ns(t) -> np.ndarray:
    """UTC datetimes (any resolution) -> int64 nanoseconds."""
    return pd.DatetimeIndex(pd.to_datetime(t, utc=True)).as_unit("ns").asi8


def _utc(t) -> pd.Timestamp:
    t = pd.Timestamp(t)
    return t.tz_localize("UTC") if t.tzinfo is None else t.tz_convert("UTC")


def _dir(env: str, default: str) -> Path:
    return Path(os.environ.get(env, default))


def wind_downtime_dir() -> Path:
    return _dir("WIND_DOWNTIME_DATA_DIR", "../wind-downtime-availability/data")


def wind_om_dir() -> Path:
    return _dir("WIND_OM_DATA_DIR", "../wind-om-reliability/data")


def pvdaq_long_dir() -> Path:
    return _dir("PVDAQ_LONG_DIR", "../pv-budget-bridge/data/pvdaq")


def farms_dir(data_dir: Path = DATA_DIR) -> Path:
    return Path(data_dir) / "farms"


# --- ERA5 at the sites ------------------------------------------------------------------

WIND_VARS = {"wind_speed_10m": "ws10", "wind_speed_100m": "ws100"}
PV_VARS = {"shortwave_radiation": "ghi", "temperature_2m": "t2m", "wind_speed_10m": "ws10"}
KELMARSH_ERA5_SEED = Path("../era5-check-pl/data/openmeteo/era5_kelmarsh_20160101_20211231.parquet")


def _fetch_point(lat: float, lon: float, start: str, end: str, variables: dict, out: Path, data_dir: Path,
                 budget: float) -> pd.DataFrame:
    ledger = weather.Ledger(Path(data_dir) / "openmeteo_ledger.json")
    weight = weather.request_weight(len(variables), start, end)
    if ledger.used() + weight > budget:
        raise RuntimeError(f"Open-Meteo budget: {ledger.used():.0f} + {weight:.0f} > {budget:.0f} in 24 h")
    params = {"latitude": lat, "longitude": lon, "start_date": start, "end_date": end,
              "hourly": ",".join(variables), "models": "era5", "timezone": "GMT", "timeformat": "unixtime",
              "wind_speed_unit": "ms"}
    r = requests.get(weather.URL, params=params, timeout=300, headers={"Accept-Encoding": "gzip"})
    if r.status_code == 429:
        ledger.add(1.0)
        raise weather.RateLimited(r.text[:200])
    r.raise_for_status()
    ledger.add(weight)
    payload = r.json()
    h = payload["hourly"]
    idx = pd.to_datetime(np.asarray(h["time"], dtype="int64"), unit="s", utc=True)
    df = pd.DataFrame({v: np.asarray(h[k], dtype="float64") for k, v in variables.items()}, index=idx)
    df.index.name = "time"
    out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(out)
    meta = {k: payload.get(k) for k in ("latitude", "longitude", "elevation")}
    out.with_suffix(".json").write_text(json.dumps({"requested": [lat, lon], "grid": meta, "weight": weight}))
    log.info("ERA5 %s: %d hours, weight %.1f", out.name, len(df), weight)
    return df


def fetch_all(data_dir: Path = DATA_DIR, budget: float = 2400.0, pause_s: float = 15.0) -> list[Path]:
    """ERA5 at the three farms (ws10, ws100) and the PV site (GHI, T2m, ws10); cached, one request each.

    Kelmarsh is copied from the cache of the sibling project era5-check-pl when it is there (same API and model).
    """
    got = []
    for f in FARMS.values():
        out = farms_dir(data_dir) / f"era5_{f.key}.parquet"
        if not out.exists():
            seed = KELMARSH_ERA5_SEED if f.key == "kelmarsh" else None
            if seed is not None and seed.exists():
                d = pd.read_parquet(seed).rename(columns=WIND_VARS)
                out.parent.mkdir(parents=True, exist_ok=True)
                d[["ws10", "ws100"]].to_parquet(out)
                log.info("ERA5 Kelmarsh copied from %s", seed)
            else:
                _fetch_point(f.lat, f.lon, "2016-01-01", f.era5_end, WIND_VARS, out, data_dir, budget)
                time.sleep(pause_s)
        got.append(out)
    out = farms_dir(data_dir) / "era5_pvdaq34.parquet"
    if not out.exists():
        y = PV_SITE["years"]
        _fetch_point(PV_SITE["lat"], PV_SITE["lon"], f"{y[0]}-01-01", f"{y[-1]}-12-31", PV_VARS, out, data_dir,
                     budget)
    got.append(out)
    return got


def farm_era5(farm: Farm, data_dir: Path = DATA_DIR) -> pd.DataFrame:
    """Hourly ERA5 ws10, ws100 (m/s, instantaneous at the top of the hour, UTC) and the hub-height speed."""
    d = pd.read_parquet(farms_dir(data_dir) / f"era5_{farm.key}.parquet")
    d = d[~d.index.duplicated()].asfreq("1h")
    d["ws_hub"] = hub_speed(d["ws10"].to_numpy(), d["ws100"].to_numpy(), farm.hub_m)
    return d


def hub_speed(ws10: np.ndarray, ws100: np.ndarray, hub_m: float, alpha_max: float = 0.6) -> np.ndarray:
    """Power-law interpolation between 10 and 100 m with the hourly shear exponent of the two ERA5 levels."""
    ws10 = np.asarray(ws10, dtype=float)
    ws100 = np.asarray(ws100, dtype=float)
    with np.errstate(divide="ignore", invalid="ignore"):
        alpha = np.log(ws100 / ws10) / np.log(10.0)
    alpha = np.clip(np.nan_to_num(alpha, nan=1 / 7, posinf=alpha_max, neginf=0.0), 0.0, alpha_max)
    return ws100 * (hub_m / 100.0) ** alpha


# --- event logs: forced outages and observation windows ---------------------------------

def _read_status_csv(raw: bytes) -> pd.DataFrame:
    df = pd.read_csv(io.BytesIO(raw), comment="#")
    return pd.DataFrame({
        "start": pd.to_datetime(df["Timestamp start"], format="ISO8601").dt.tz_localize("UTC"),
        "end": pd.to_datetime(df["Timestamp end"], errors="coerce", format="ISO8601").dt.tz_localize("UTC"),
        "status": df["Status"].astype(str),
        "category": df["IEC category"].fillna("").astype(str).str.strip(),
    })


def greenbyte_log(farm: Farm) -> pd.DataFrame:
    """All status records (turbine, start, end, status, category) of a Greenbyte farm."""
    if farm.key == "kelmarsh":
        root = kelmarsh_dir()
        zips = sorted((root / "raw").glob("Kelmarsh_SCADA_*.zip")) or sorted(root.glob("Kelmarsh_SCADA_*.zip"))
        pattern = r"Status_Kelmarsh_(\d+)_"
    else:
        zips = sorted((wind_om_dir() / "raw" / "penmanshiel").glob("status_*.zip"))
        pattern = r"Status_Penmanshiel_(\d+)_"
    if not zips:
        raise FileNotFoundError(f"no status archives for {farm.label}")
    parts = []
    for zp in zips:
        with zipfile.ZipFile(zp) as z:
            for name in z.namelist():
                m = re.match(pattern, Path(name).name)
                if m:
                    parts.append(_read_status_csv(z.read(name)).assign(turbine=f"T{int(m[1]):02d}"))
    d = pd.concat(parts, ignore_index=True)
    # the quarterly 2023 Penmanshiel exports overlap the yearly ones at their edges
    return d.drop_duplicates(["turbine", "start", "end", "status", "category"]).sort_values(["turbine", "start"])


def hot_log() -> pd.DataFrame:
    """Hill of Towie alarm log (turbine, start, end, code); open alarms keep a missing end."""
    zips = sorted((wind_om_dir() / "raw" / "hill_of_towie").glob("alarms_*.zip"))
    if not zips:
        raise FileNotFoundError("no Hill of Towie alarm archives (WIND_OM_DATA_DIR)")
    parts = []
    for zp in zips:
        with zipfile.ZipFile(zp) as z:
            for name in sorted(z.namelist()):
                parts.append(pd.read_csv(z.open(name), dtype=str))
    a = pd.concat(parts, ignore_index=True)
    out = pd.DataFrame({
        "turbine": "T" + (pd.to_numeric(a["StationNr"]) - HOT_STATION_OFFSET).astype(int).map("{:02d}".format),
        "start": pd.to_datetime(a["TimeOn"], format="ISO8601").dt.tz_localize("UTC"),
        "end": pd.to_datetime(a["TimeOff"], errors="coerce", format="ISO8601").dt.tz_localize("UTC"),
        "code": pd.to_numeric(a["Alarmcode"], errors="coerce").astype("Int64"),
    })
    out = out[out["turbine"].isin([f"T{i:02d}" for i in range(1, 22)])]
    return out.drop_duplicates(["turbine", "start", "code"]).sort_values(["turbine", "start"]).reset_index(drop=True)


def merge_intervals(ev: pd.DataFrame, gap_h: float = MERGE_GAP_H) -> pd.DataFrame:
    """Merge overlapping or nearly contiguous intervals per turbine -> turbine, start, end, duration_h."""
    rows = []
    gap = pd.Timedelta(hours=gap_h)
    for t, g in ev.sort_values("start").groupby("turbine", sort=True):
        cur_s = cur_e = None
        for s, e in zip(g["start"], g["end"]):
            if cur_s is None:
                cur_s, cur_e = s, e
            elif s <= cur_e + gap:
                cur_e = max(cur_e, e)
            else:
                rows.append((t, cur_s, cur_e))
                cur_s, cur_e = s, e
        if cur_s is not None:
            rows.append((t, cur_s, cur_e))
    out = pd.DataFrame(rows, columns=["turbine", "start", "end"])
    out["duration_h"] = (out["end"] - out["start"]).dt.total_seconds() / 3600
    return out


def subtract_intervals(a: pd.DataFrame, b: pd.DataFrame) -> pd.DataFrame:
    """Per turbine, intervals of ``a`` minus the union of ``b`` (both turbine, start, end)."""
    rows = []
    for t, ga in a.groupby("turbine", sort=True):
        gb = merge_intervals(b[b["turbine"] == t], 0.0) if len(b) else b
        cut = list(zip(gb["start"], gb["end"])) if len(gb) else []
        for s, e in zip(ga["start"], ga["end"]):
            pieces = [(s, e)]
            for cs, ce in cut:
                nxt = []
                for ps, pe in pieces:
                    if ce <= ps or cs >= pe:
                        nxt.append((ps, pe))
                        continue
                    if cs > ps:
                        nxt.append((ps, cs))
                    if ce < pe:
                        nxt.append((ce, pe))
                pieces = nxt
            rows += [(t, ps, pe) for ps, pe in pieces if pe > ps]
    return pd.DataFrame(rows, columns=["turbine", "start", "end"])


def observation_windows(records: pd.DataFrame, t0, t1, turbines: list[str],
                        comm_gaps: pd.DataFrame | None = None) -> pd.DataFrame:
    """Per turbine: [t0, t1) minus silences of the log longer than ``MAX_SILENCE`` and communication gaps."""
    t0, t1 = _utc(t0).as_unit("ns"), _utc(t1).as_unit("ns")
    rows = []
    for t in turbines:
        r = records[records["turbine"] == t]
        stamps = np.sort(np.r_[_ns(r["start"]), _ns(r["end"].dropna())])
        stamps = stamps[(stamps >= t0.value) & (stamps <= t1.value)]
        edges = np.r_[t0.value, stamps, t1.value]
        gaps = np.flatnonzero(np.diff(edges) > MAX_SILENCE.value)
        s = t0.value
        for k in gaps:
            if edges[k] > s:
                rows.append((t, s, edges[k]))
            s = edges[k + 1]
        if t1.value > s:
            rows.append((t, s, t1.value))
    win = pd.DataFrame(rows, columns=["turbine", "start", "end"])
    win["start"] = pd.to_datetime(win["start"].astype("int64"), utc=True)
    win["end"] = pd.to_datetime(win["end"].astype("int64"), utc=True)
    if comm_gaps is not None and len(comm_gaps):
        win = subtract_intervals(win, comm_gaps)
    win["hours"] = (win["end"] - win["start"]).dt.total_seconds() / 3600
    return win


def forced_outages(farm: Farm) -> tuple[pd.DataFrame, pd.DataFrame]:
    """(outages, windows): merged forced-outage intervals and the observation windows of every turbine."""
    t0, t1 = _utc(farm.log_start), _utc(farm.log_end)
    if farm.key == "hill_of_towie":
        lg = hot_log()
        turbines = [f"T{i:02d}" for i in range(1, farm.n_turbines + 1)]
        fo = lg[lg["code"].isin(HOT_FORCED_CODES)].copy()
        # an open alarm lasts until the turbine's next alarm (often a generator cut-in), at most a week
        nxt = lg.groupby("turbine")["start"].shift(-1)
        open_ = fo["end"].isna()
        fo.loc[open_, "end"] = nxt[fo.index[open_]].fillna(fo.loc[open_, "start"] + MAX_SILENCE).clip(
            upper=fo.loc[open_, "start"] + MAX_SILENCE)
        comm = None
        records = lg
    else:
        lg = greenbyte_log(farm)
        turbines = sorted(lg["turbine"].unique())
        fo = lg[(lg["status"] == "Stop") & (lg["category"].str.lower() == "forced outage")].dropna(subset=["end"])
        comm = lg[(lg["status"] == "Communication") & ((lg["end"] - lg["start"]) >= MIN_COMM_GAP)]
        records = lg
    fo = fo[(fo["end"] > fo["start"]) & (fo["end"] > t0) & (fo["start"] < t1)].copy()
    fo["start"] = fo["start"].clip(lower=t0)
    fo["end"] = fo["end"].clip(upper=t1)
    win = observation_windows(records, t0, t1, turbines, comm)
    out = merge_intervals(fo)
    return out, win


# --- SCADA: hourly farm capacity factor ---------------------------------------------------

def _scada_10min(farm: Farm) -> pd.DataFrame:
    """turbine, time (UTC, start of the 10-minute interval), power_kw."""
    if farm.key == "kelmarsh":
        parts = []
        for y in farm.scada_years:
            d = pd.read_parquet(kelmarsh_dir() / "processed" / f"scada_{y}.parquet", columns=["power_kw"])
            d = d.reset_index()
            parts.append(pd.DataFrame({"turbine": "T" + d["turbine"].astype(int).map("{:02d}".format),
                                       "time": d["timestamp"], "power_kw": d["power_kw"].astype(float)}))
        return pd.concat(parts, ignore_index=True)
    parts = []
    for y in farm.scada_years:
        d = pd.read_parquet(wind_downtime_dir() / "processed" / f"{farm.key}_{y}.parquet",
                            columns=["timestamp", "turbine_id", "power"])
        parts.append(pd.DataFrame({"turbine": d["turbine_id"].astype(str), "time": d["timestamp"],
                                   "power_kw": d["power"].astype(float)}))
    return pd.concat(parts, ignore_index=True)


def covered_until(starts: np.ndarray, ends: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Time covered by the intervals [starts, ends) before each x (same units; intervals must not overlap)."""
    s = np.asarray(starts, dtype=np.float64)
    e = np.asarray(ends, dtype=np.float64)
    x = np.asarray(x, dtype=np.float64)
    if not len(s):
        return np.zeros_like(x)
    order = np.argsort(s)
    s, e = s[order], e[order]
    cum_len = np.r_[0.0, np.cumsum(e - s)]
    # all intervals starting at or before x, minus the part of the last one that lies after x
    k = np.searchsorted(s, x, side="right")
    over = np.where(k > 0, np.maximum(e[np.clip(k - 1, 0, None)] - x, 0.0), 0.0)
    return cum_len[k] - over


def interval_coverage(starts: np.ndarray, ends: np.ndarray, grid: np.ndarray) -> np.ndarray:
    """Covered time inside every cell [grid[k], grid[k+1]) (non-overlapping intervals, same units)."""
    return np.diff(covered_until(starts, ends, grid))


def scada_hourly(farm: Farm, outages: pd.DataFrame, min_valid: float = 0.75) -> pd.DataFrame:
    """Hourly farm capacity factor from the 10-minute SCADA (UTC, hour start).

    ``cf_all``: mean of power / rated over the turbine-intervals with data (forced outages count as zero);
    ``cf_avail``: the same without turbine-intervals that are half or more in a forced outage, i.e. the output
    of the turbines that were technically available. An hour needs ``min_valid`` of the turbine-intervals with
    data (``cf_avail``: of the available ones), otherwise it is NaN.
    """
    d = _scada_10min(farm)
    d["time"] = pd.to_datetime(d["time"], utc=True).dt.floor("10min")
    d = d.drop_duplicates(["turbine", "time"])
    step = 600 * 10**9
    d["fo"] = 0.0
    for t, g in d.groupby("turbine"):
        o = outages[outages["turbine"] == t]
        if o.empty:
            continue
        ns = _ns(g["time"])
        a, b = _ns(o["start"]), _ns(o["end"])
        d.loc[g.index, "fo"] = (covered_until(a, b, ns + step) - covered_until(a, b, ns)) / step
    d["cf"] = (d["power_kw"] / farm.rated_kw).clip(0.0, 1.0)
    d["hour"] = d["time"].dt.floor("1h")
    valid = d["cf"].notna()
    avail = valid & (d["fo"] < 0.5)
    per = 6 * farm.n_turbines
    g = d.assign(cf_v=d["cf"].where(valid), cf_a=d["cf"].where(avail), v=valid, a=avail).groupby("hour")
    h = pd.DataFrame({"cf_all": g["cf_v"].mean(), "cf_avail": g["cf_a"].mean(),
                      "n_valid": g["v"].sum(), "n_avail": g["a"].sum()})
    h.loc[h["n_valid"] < min_valid * per, ["cf_all", "cf_avail"]] = np.nan
    h.loc[h["n_avail"] < min_valid * h["n_valid"].clip(lower=1), "cf_avail"] = np.nan
    h.loc[h["n_avail"] < 0.5 * per, "cf_avail"] = np.nan
    h.index.name = "time"
    return h


def cached_scada_hourly(farm: Farm, outages: pd.DataFrame, data_dir: Path = DATA_DIR) -> pd.DataFrame:
    path = farms_dir(data_dir) / f"scada_hourly_{farm.key}.parquet"
    if path.exists():
        return pd.read_parquet(path)
    h = scada_hourly(farm, outages)
    path.parent.mkdir(parents=True, exist_ok=True)
    h.to_parquet(path)
    return h


# --- PV: PVDAQ system 34 ---------------------------------------------------------------------

def pvdaq_hourly(poa_min: float = 200.0) -> pd.DataFrame:
    """Hourly mean AC power (kW, UTC hour start) of PVDAQ 34 and the hours in a forced outage.

    An hour is in a forced outage when a quarter hour of it has plane-of-array irradiance above ``poa_min``
    and no AC output (zero or no inverter record while the met station logs), as in ``structural``.
    """
    d = pd.read_parquet(pvdaq_long_dir() / "long_34.parquet")
    d = d[((d["device_id"] == "INV1") & (d["signal"] == "p_ac")) |
          ((d["device_id"] == "MET") & (d["signal"] == "poa_irradiance"))]
    w = d.pivot_table(index="timestamp", columns="signal", values="value", aggfunc="first")
    y = PV_SITE["years"]
    grid = pd.date_range(f"{y[0]}-01-01", f"{y[-1]}-12-31 23:45", freq="15min")
    w = w.reindex(grid)
    out_q = (w["poa_irradiance"] > poa_min) & ~(w["p_ac"] > 0)
    # local standard time -> UTC
    idx = grid - pd.Timedelta(hours=PV_SITE["utc_offset_h"])
    w.index = idx.tz_localize("UTC")
    out_q.index = w.index
    h = pd.DataFrame({"ac_kw": w["p_ac"].resample("1h").mean(),
                      "n_ac": w["p_ac"].notna().resample("1h").sum(),
                      "poa": w["poa_irradiance"].resample("1h").mean(),
                      "outage": out_q.resample("1h").max().astype(bool)})
    h.loc[h["n_ac"] < 3, "ac_kw"] = np.nan
    # inverter channels missing while the met station logs sun = no output (the logged outage pattern)
    h["ac_kw_all"] = h["ac_kw"].where(~(h["ac_kw"].isna() & (h["poa"] > poa_min)), 0.0)
    h.index.name = "time"
    return h


def pv_era5(data_dir: Path = DATA_DIR) -> pd.DataFrame:
    d = pd.read_parquet(farms_dir(data_dir) / "era5_pvdaq34.parquet")
    return d[~d.index.duplicated()].asfreq("1h")
