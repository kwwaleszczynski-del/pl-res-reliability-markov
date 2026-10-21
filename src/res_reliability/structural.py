"""Structural (technical) reliability: FOR, MTTF and MTTR from open event logs.

Wind: Kelmarsh wind farm status logs (6 x Senvion MM92, 2016-2021; Plumley 2022, Zenodo,
CC BY 4.0). Forced outages are the ``Stop`` records classified by the operator under the
IEC 61400-26 category "Forced outage"; overlapping or nearly contiguous stops of the same
turbine are merged into one event. Grid-side stops ("Out of Electrical Specification"),
scheduled maintenance and requested shutdowns are reported separately, not as forced outages.

PV: NREL PVDAQ system 34 (one inverter, 2016-2017, CC BY 4.0). An outage is a daylight
interval with plane-of-array irradiance above 200 W/m2 and no AC output (or no inverter
record while the irradiance sensor is logging); FOR is counted over daylight time only.
"""

from __future__ import annotations

import glob
import io
import logging
import zipfile
from pathlib import Path

import numpy as np
import pandas as pd

log = logging.getLogger(__name__)

KELMARSH_COD = pd.Timestamp("2016-04-15", tz="UTC")
KELMARSH_END = pd.Timestamp("2022-01-01", tz="UTC")


def load_kelmarsh_stops(data_dir: Path) -> pd.DataFrame:
    """All ``Stop`` records from the status files inside the yearly SCADA zips."""
    frames = []
    zips = sorted(glob.glob(str(Path(data_dir) / "raw" / "Kelmarsh_SCADA_*.zip"))) or \
        sorted(glob.glob(str(Path(data_dir) / "Kelmarsh_SCADA_*.zip")))
    if not zips:
        raise FileNotFoundError(f"no Kelmarsh_SCADA_*.zip under {data_dir} (set KELMARSH_DATA_DIR)")
    for zp in zips:
        with zipfile.ZipFile(zp) as z:
            for name in z.namelist():
                if not name.startswith("Status_"):
                    continue
                df = pd.read_csv(io.BytesIO(z.read(name)), comment="#")
                df = df[df["Status"] == "Stop"].copy()
                df["turbine"] = int(name.split("_")[2])
                frames.append(df)
    d = pd.concat(frames, ignore_index=True)
    out = pd.DataFrame({
        "turbine": d["turbine"],
        "start": pd.to_datetime(d["Timestamp start"]).dt.tz_localize("UTC"),
        "end": pd.to_datetime(d["Timestamp end"], errors="coerce").dt.tz_localize("UTC"),
        "category": d["IEC category"].fillna("").str.strip(),
        "message": d["Message"].fillna(""),
    })
    out = out.dropna(subset=["end"])
    out = out[out["end"] > out["start"]]
    return out.drop_duplicates(["turbine", "start", "end", "category"]).sort_values(["turbine", "start"])


def merge_events(ev: pd.DataFrame, gap_h: float = 0.5) -> pd.DataFrame:
    """Merge overlapping or nearly contiguous intervals per turbine."""
    rows = []
    for t, g in ev.sort_values("start").groupby("turbine"):
        cur_s, cur_e = None, None
        for s, e in zip(g["start"], g["end"]):
            if cur_s is None:
                cur_s, cur_e = s, e
            elif s <= cur_e + pd.Timedelta(hours=gap_h):
                cur_e = max(cur_e, e)
            else:
                rows.append((t, cur_s, cur_e))
                cur_s, cur_e = s, e
        if cur_s is not None:
            rows.append((t, cur_s, cur_e))
    out = pd.DataFrame(rows, columns=["turbine", "start", "end"])
    out["duration_h"] = (out["end"] - out["start"]).dt.total_seconds() / 3600
    return out


def event_stats(events: pd.DataFrame, exposure_h: float, n_units: int, min_duration_h: float = 0.0) -> dict:
    """FOR, MTTF, MTTR and failure rate from merged events and the total unit-hours observed."""
    e = events[events["duration_h"] >= min_duration_h]
    down = float(e["duration_h"].sum())
    n = int(len(e))
    up = exposure_h - down
    return {
        "events": n,
        "downtime_h": down,
        "exposure_unit_h": exposure_h,
        "for": down / exposure_h,
        "mttr_h": down / n if n else float("nan"),
        "mttf_h": up / n if n else float("nan"),
        "failures_per_unit_year": n / (exposure_h / 8760.0),
        "n_units": n_units,
    }


def kelmarsh_reliability(data_dir: Path) -> dict:
    stops = load_kelmarsh_stops(data_dir)
    stops = stops[(stops["end"] > KELMARSH_COD) & (stops["start"] < KELMARSH_END)].copy()
    stops["start"] = stops["start"].clip(lower=KELMARSH_COD)
    stops["end"] = stops["end"].clip(upper=KELMARSH_END)
    n_turb = int(stops["turbine"].nunique())
    exposure = n_turb * (KELMARSH_END - KELMARSH_COD).total_seconds() / 3600
    forced = merge_events(stops[stops["category"].str.lower() == "forced outage"])
    out = {
        "source": "Kelmarsh wind farm status logs 2016-04-15..2021-12-31, 6 x Senvion MM92",
        "all_forced": event_stats(forced, exposure, n_turb),
        "forced_ge_1h": event_stats(forced, exposure, n_turb, min_duration_h=1.0),
        "durations_h": forced["duration_h"].round(3).tolist(),
    }
    other = {}
    for cat in sorted(stops["category"].unique()):
        if not cat or cat.lower() == "forced outage":
            continue
        m = merge_events(stops[stops["category"] == cat])
        other[cat] = {"events": int(len(m)), "downtime_h": float(m["duration_h"].sum()),
                      "share_of_time": float(m["duration_h"].sum() / exposure)}
    out["other_categories"] = other
    yearly = []
    for y in range(2016, 2022):
        lo = max(pd.Timestamp(f"{y}-01-01", tz="UTC"), KELMARSH_COD)
        hi = pd.Timestamp(f"{y + 1}-01-01", tz="UTC")
        f = forced[(forced["start"] >= lo) & (forced["start"] < hi)]
        exp_y = n_turb * (hi - lo).total_seconds() / 3600
        yearly.append({"year": y, **event_stats(f, exp_y, n_turb)})
    out["yearly"] = yearly
    return out


def load_pvdaq(folder: Path) -> pd.DataFrame:
    files = sorted(glob.glob(str(Path(folder) / "*.csv")))
    if not files:
        raise FileNotFoundError(f"no PVDAQ csv files in {folder}")
    cols = ["measured_on", "ac_power_hw__2695", "poa_irradiance__2679"]
    # days without inverter channels (e.g. during an inverter outage) get NaN power
    d = pd.concat([pd.read_csv(f).reindex(columns=cols) for f in files], ignore_index=True)
    d.columns = ["time", "ac_w", "poa"]
    d["time"] = pd.to_datetime(d["time"]).dt.round("15min")
    return d.drop_duplicates("time").set_index("time").sort_index()


def pvdaq_reliability(folder: Path, poa_min: float = 200.0, step_min: int = 15) -> dict:
    """Outages of a single PV inverter, counted over daylight time (POA above ``poa_min``)."""
    d = load_pvdaq(folder)
    d = d.asfreq(f"{step_min}min")
    sunny = d["poa"] > poa_min  # irradiance logged, so the logger itself was running
    out_now = sunny & ~(d["ac_w"] > 0)
    obs = sunny
    # daylight-time clock: an outage continues across nights until AC output returns
    s = d.loc[obs | out_now, ["ac_w"]].copy()
    s["down"] = out_now[s.index]
    run_id = (s["down"] != s["down"].shift()).cumsum()
    runs = s.groupby(run_id).agg(down=("down", "first"), n=("down", "size"),
                                 start=("down", lambda x: x.index[0]), end=("down", lambda x: x.index[-1]))
    dt_h = step_min / 60.0
    down_runs = runs[runs["down"]]
    events = pd.DataFrame({
        "turbine": 0,
        "start": down_runs["start"].to_numpy(),
        "end": down_runs["end"].to_numpy(),
        "duration_h": down_runs["n"].to_numpy() * dt_h,
    })
    exposure = float(len(s)) * dt_h
    out = {
        "source": "NREL PVDAQ system 34, 2016-2017, single inverter; daylight time with POA > %g W/m2" % poa_min,
        "all": event_stats(events, exposure, 1),
        "ge_1h": event_stats(events, exposure, 1, min_duration_h=1.0),
        "events_ge_1h": [
            {"start": str(r.start), "end": str(r.end), "daylight_h": float(r.duration_h)}
            for r in events[events["duration_h"] >= 1.0].itertuples()
        ],
    }
    return out
