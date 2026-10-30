"""Hourly national power demand of the Polish power system (KSE) from PSE S.A.

Two public PSE sources are stitched together, both read with the pse-opendata client
(https://pypi.org/project/pse-opendata/), which caches whole months under ``<data_dir>/pse``:

* ``api.raporty.pse.pl`` report ``kse-load`` (field ``load_actual``, 15-minute values).
  The new reporting platform starts on 2024-06-14.
* the legacy PSE CSV export ``PL_WYK_KSE`` ("Krajowe zapotrzebowanie na moc", hourly, local
  time, hour 1 = 00:00-01:00, extra autumn hour labelled ``2A``) for everything before that.
  Every row is placed by its hour label, so a missing hour leaves a gap instead of shifting
  the rest of the day.

The two do not overlap; the API part is rescaled to the legacy level with a factor matched
over the two weeks either side of the switch (``splice_factor``; 1.007 for this data).
Everything is stored with a UTC hourly index labelled by the start of the hour.

Usage:
    python -m res_reliability fetch-load
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
import pse_opendata as pse

from res_reliability.config import DATA_DIR, HOURS, YEARS

log = logging.getLogger(__name__)

API_START = pd.Timestamp("2024-06-14")  # first business day on api.raporty.pse.pl


def make_client(data_dir: Path = DATA_DIR) -> pse.Client:
    """pse-opendata client with its monthly cache in ``<data_dir>/pse``."""
    return pse.Client(cache_dir=data_dir / "pse")


def fetch_legacy(start: pd.Timestamp, end: pd.Timestamp, client: pse.Client) -> pd.Series:
    """Hourly demand of the legacy export for local days start..end, labelled by hour start (UTC).

    Whole months are requested, so that they are cached; the server is slow and rate-limited,
    and pse-opendata spaces the requests.
    """
    frame = client.legacy("PL_WYK_KSE", start.replace(day=1), end, fields=["demand"])
    frame = frame[frame["business_date"] >= start.date()]
    return pd.Series(frame["demand"].to_numpy(), index=pd.DatetimeIndex(frame["start_utc"], name="time"),
                     name="demand_mw")


def fetch_api(start: pd.Timestamp, end: pd.Timestamp, client: pse.Client) -> pd.Series:
    """``kse-load`` for business days start..end -> hourly mean labelled by hour start (UTC)."""
    frame = client.get("kse-load", start, end, fields=["load_actual"])
    s = pd.Series(frame["load_actual"].to_numpy(), index=pd.DatetimeIndex(frame["start_utc"]))
    hourly = s.resample("1h").mean()
    hourly.index.name = "time"
    return hourly.rename("demand_mw")


def splice_factor(legacy: pd.Series, api: pd.Series, days: int = 14) -> float:
    """Level ratio legacy / API over the two weeks either side of the switch to the new platform.

    The two sources do not overlap and ``kse-load`` sits a few per cent below the legacy
    "Krajowe zapotrzebowanie na moc" (a definition change with the new market model), so
    the API part is rescaled. Hours are matched by weekday and hour of day (local time).
    """
    lo = legacy[API_START.tz_localize("Europe/Warsaw") - pd.Timedelta(days=days):]
    hi = api[:API_START.tz_localize("Europe/Warsaw") + pd.Timedelta(days=days)]
    key = lambda s: pd.Index([f"{t.weekday()}-{t.hour}" for t in s.index.tz_convert("Europe/Warsaw")])  # noqa: E731
    a = lo.groupby(key(lo)).mean()
    b = hi.groupby(key(hi)).mean()
    common = a.index.intersection(b.index)
    return float(a[common].mean() / b[common].mean())


def cache_path(data_dir: Path = DATA_DIR) -> Path:
    return data_dir / "load" / "kse_demand_hourly.parquet"


def fetch(years: list[int] = YEARS, data_dir: Path = DATA_DIR) -> pd.Series:
    out = cache_path(data_dir)
    if out.exists():
        return pd.read_parquet(out)["demand_mw"]
    client = make_client(data_dir)
    first, last = pd.Timestamp(f"{years[0] - 1}-12-31"), pd.Timestamp(f"{years[-1] + 1}-01-01")
    out.parent.mkdir(parents=True, exist_ok=True)
    legacy = fetch_legacy(first, min(last, API_START - pd.Timedelta(days=1)), client)
    api = fetch_api(API_START, last, client) if last >= API_START else pd.Series(dtype="float64")
    if len(api):
        k = splice_factor(legacy, api)
        log.info("api.raporty.pse.pl level rescaled by %.4f to the legacy definition", k)
        (out.parent / "splice_factor.txt").write_text(f"{k:.6f}")
        api = api * k
    s = pd.concat([legacy, api]).sort_index()
    s = s[~s.index.duplicated()]
    s = s.rename("demand_mw")
    out.parent.mkdir(parents=True, exist_ok=True)
    s.to_frame().to_parquet(out)
    return s


def year_matrix(s: pd.Series, years: list[int] = YEARS) -> np.ndarray:
    """(n_years, 8760) demand matrix, one row per UTC calendar year (leap day 31 Dec dropped).

    Gaps up to 6 hours (the source misses a few evening hours on single days) are
    interpolated linearly; longer gaps raise.
    """
    rows = []
    for y in years:
        idx = pd.date_range(f"{y}-01-01", periods=HOURS, freq="1h", tz="UTC")
        v = s.reindex(idx)
        miss = v.isna().to_numpy()
        if miss.any():
            runs = np.diff(np.r_[0, miss.astype(int), 0])
            longest = int((np.flatnonzero(runs == -1) - np.flatnonzero(runs == 1)).max())
            if longest > 6:
                raise ValueError(f"{y}: a gap of {longest} demand hours is too long to interpolate")
            v = v.interpolate(limit_direction="both")
            log.info("%d: interpolated %d missing demand hours", y, int(miss.sum()))
        rows.append(v.to_numpy())
    return np.vstack(rows)
