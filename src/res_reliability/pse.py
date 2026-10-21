"""Hourly national power demand of the Polish power system (KSE) from PSE S.A.

Two public PSE sources are stitched together:

* ``api.raporty.pse.pl`` report ``kse-load`` (field ``load_actual``, 15-minute values,
  ``dtime_utc`` = end of the interval). The new reporting platform starts on 2024-06-14.
* the legacy PSE CSV export ``PL_WYK_KSE`` ("Krajowe zapotrzebowanie na moc", hourly, local
  time, hour 1 = 00:00-01:00, extra autumn hour labelled ``2A``) for everything before that.

The two do not overlap; the API part is rescaled to the legacy level with a factor matched
over the two weeks either side of the switch (``splice_factor``; 1.007 for this data).
Everything is stored with a UTC hourly index labelled by the start of the hour.

Usage:
    python -m res_reliability fetch-load
"""

from __future__ import annotations

import io
import logging
import time
from pathlib import Path

import numpy as np
import pandas as pd
import requests

from res_reliability.config import DATA_DIR, HOURS, YEARS

log = logging.getLogger(__name__)

LEGACY_URL = "https://www.pse.pl/getcsv/-/export/csv/PL_WYK_KSE/data_od/{start}/data_do/{end}"
API_URL = "https://api.raporty.pse.pl/api/kse-load"
API_START = pd.Timestamp("2024-06-14")  # first business day on api.raporty.pse.pl
HEADERS = {"User-Agent": "Mozilla/5.0 (research script; public data)"}


def _day_labels(n_hours: int) -> list[str]:
    """Hour labels of a local day with 23 (spring), 24 or 25 (autumn) hours, in time order."""
    if n_hours == 25:
        return ["1", "2", "2A", *[str(h) for h in range(3, 25)]]
    if n_hours == 23:
        return ["1", *[str(h) for h in range(3, 25)]]
    return [str(h) for h in range(1, 25)]


def parse_legacy_csv(text: str) -> pd.Series:
    """Parse a legacy PSE CSV (``Data;Godzina;Krajowe zapotrzebowanie na moc;...``).

    Every row is placed by its hour label, so a missing hour leaves a gap instead of
    shifting the rest of the day.
    """
    df = pd.read_csv(io.StringIO(text), sep=";", decimal=",", usecols=[0, 1, 2], dtype=str)
    df.columns = ["date", "hour", "demand"]
    df["hour"] = df["hour"].str.strip().str.upper().str.lstrip("0")
    df["demand"] = pd.to_numeric(df["demand"].str.replace(",", ".", regex=False), errors="coerce")
    stamps, values = [], []
    for day, g in df.groupby("date", sort=True):
        start = pd.Timestamp(day).tz_localize("Europe/Warsaw")
        end = (pd.Timestamp(day) + pd.Timedelta(days=1)).tz_localize("Europe/Warsaw")
        n = int(round((end - start).total_seconds() / 3600))
        pos = {lab: i for i, lab in enumerate(_day_labels(n))}
        start_utc = start.tz_convert("UTC")
        for lab, v in zip(g["hour"], g["demand"]):
            if lab in pos:
                stamps.append(start_utc + pd.Timedelta(hours=pos[lab]))
                values.append(v)
    s = pd.Series(values, index=pd.DatetimeIndex(stamps, name="time"), name="demand_mw", dtype="float64")
    return s[~s.index.duplicated()].sort_index()


def parse_api_values(rows: list[dict]) -> pd.Series:
    """15-minute ``kse-load`` records -> hourly mean labelled by hour start (UTC)."""
    df = pd.DataFrame(rows)
    end = pd.to_datetime(df["dtime_utc"]).dt.tz_localize("UTC")
    start = end - pd.Timedelta(minutes=15)
    s = pd.Series(pd.to_numeric(df["load_actual"], errors="coerce").to_numpy(), index=start)
    s = s[~s.index.duplicated()].sort_index()
    hourly = s.resample("1h").mean()
    hourly.index.name = "time"
    return hourly.rename("demand_mw")


def _get(session: requests.Session, url: str, params: dict | None = None, tries: int = 10) -> requests.Response:
    for k in range(tries):
        r = session.get(url, params=params, headers=HEADERS, timeout=180)
        if r.status_code == 429:
            wait = 120 * (k + 1)
            log.warning("429 from %s; waiting %d s", url.split("/")[2], wait)
            time.sleep(wait)
            continue
        r.raise_for_status()
        return r
    raise RuntimeError(f"giving up on {url}")


def fetch_legacy(start: pd.Timestamp, end: pd.Timestamp, session: requests.Session, cache: Path,
                 pause_s: float = 5.0) -> pd.Series:
    """Legacy export one calendar month per request, cached per month.

    The server answers 429 to long ranges (a quarter or a year) but serves a month at a
    time, and it rejects one-day ranges.
    """
    parts = []
    for month_start in pd.date_range(start.normalize().replace(day=1), end, freq="MS"):
        lo = max(start, month_start)
        hi = min(end, month_start + pd.offsets.MonthEnd(0))
        if hi <= lo:
            lo = hi - pd.Timedelta(days=1)
        f = cache / f"legacy_{lo:%Y%m%d}_{hi:%Y%m%d}.parquet"
        if not f.exists():
            url = LEGACY_URL.format(start=lo.strftime("%Y%m%d"), end=hi.strftime("%Y%m%d"))
            part = parse_legacy_csv(_get(session, url).content.decode("cp1250"))
            part.to_frame().to_parquet(f)
            log.info("legacy PSE %s..%s: %d h", lo.date(), hi.date(), len(part))
            time.sleep(pause_s)
        parts.append(pd.read_parquet(f)["demand_mw"])
    return pd.concat(parts).sort_index()


def fetch_api(start: pd.Timestamp, end: pd.Timestamp, session: requests.Session, pause_s: float = 2.0) -> pd.Series:
    """``kse-load`` by business date, month by month, following ``nextLink`` paging."""
    parts = []
    lo = max(start, API_START)
    while lo <= end:
        hi = min(lo + pd.offsets.MonthEnd(0), end)
        flt = f"business_date ge '{lo.date()}' and business_date le '{hi.date()}'"
        rows: list[dict] = []
        url, params = API_URL, {"$filter": flt, "$first": 5000}
        while url:
            payload = _get(session, url, params).json()
            rows.extend(payload.get("value", []))
            url, params = payload.get("nextLink"), None
        if rows:
            parts.append(parse_api_values(rows))
        log.info("api.raporty.pse.pl %s..%s: %d records", lo.date(), hi.date(), len(rows))
        lo = hi + pd.Timedelta(days=1)
        time.sleep(pause_s)
    return pd.concat(parts).sort_index() if parts else pd.Series(dtype="float64")


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
    session = requests.Session()
    first, last = pd.Timestamp(f"{years[0] - 1}-12-31"), pd.Timestamp(f"{years[-1] + 1}-01-01")
    out.parent.mkdir(parents=True, exist_ok=True)
    legacy = fetch_legacy(first, min(last, API_START - pd.Timedelta(days=1)), session, out.parent)
    api = fetch_api(API_START, last, session) if last >= API_START else pd.Series(dtype="float64")
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
