"""Hourly ERA5 reanalysis per grid point from the Open-Meteo archive API (no key).

Open-Meteo counts a request as ``max(1, n_vars/10 * max(1, n_days/14))`` calls and allows
10 000 calls per day per IP, shared with every other process on the machine. Each point is
fetched once (one request for the whole period) and cached as Parquet; every request is
written to a ledger so that a daily budget can be respected across restarts.

Usage:
    python -m res_reliability fetch-weather --budget 2400
"""

from __future__ import annotations

import json
import logging
import time
from datetime import date
from pathlib import Path

import numpy as np
import pandas as pd
import requests

from res_reliability.config import DATA_DIR, YEARS
from res_reliability.geo import GRID, Point, fetch_order

log = logging.getLogger(__name__)

URL = "https://archive-api.open-meteo.com/v1/archive"
VARIABLES = ["wind_speed_100m", "shortwave_radiation", "temperature_2m"]
COLUMNS = {"wind_speed_100m": "ws100", "shortwave_radiation": "ghi", "temperature_2m": "t2m"}


def request_weight(n_vars: int, start: str, end: str) -> float:
    days = (date.fromisoformat(end) - date.fromisoformat(start)).days + 1
    return max(1.0, n_vars / 10 * max(1.0, days / 14))


def weather_dir(data_dir: Path = DATA_DIR) -> Path:
    return data_dir / "weather"


def cache_path(point: Point, data_dir: Path = DATA_DIR) -> Path:
    return weather_dir(data_dir) / f"era5_{point.id}.parquet"


class Ledger:
    """Append-only record of (unix time, weight) of the requests made by this project."""

    def __init__(self, path: Path):
        self.path = path
        self.entries = json.loads(path.read_text()) if path.exists() else []

    def used(self, window_s: float = 86400.0) -> float:
        now = time.time()
        return sum(w for t, w in self.entries if now - t < window_s)

    def add(self, weight: float) -> None:
        self.entries.append([time.time(), weight])
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self.path.write_text(json.dumps(self.entries))


def parse_response(payload: dict) -> pd.DataFrame:
    hourly = payload["hourly"]
    idx = pd.to_datetime(np.asarray(hourly["time"], dtype="int64"), unit="s", utc=True)
    df = pd.DataFrame({COLUMNS[v]: np.asarray(hourly[v], dtype="float64") for v in VARIABLES}, index=idx)
    df.index.name = "time"
    return df


def fetch_point(point: Point, start: str, end: str, session: requests.Session) -> tuple[pd.DataFrame, dict]:
    params = {
        "latitude": point.lat,
        "longitude": point.lon,
        "start_date": start,
        "end_date": end,
        "hourly": ",".join(VARIABLES),
        "models": "era5",
        "timezone": "GMT",
        "timeformat": "unixtime",
        "wind_speed_unit": "ms",
    }
    r = session.get(URL, params=params, timeout=300)
    if r.status_code == 429:
        raise RateLimited(r.text[:200])
    r.raise_for_status()
    payload = r.json()
    meta = {k: payload.get(k) for k in ("latitude", "longitude", "elevation")}
    return parse_response(payload), meta


class RateLimited(RuntimeError):
    pass


def fetch_all(
    points: list[Point] = GRID,
    start: str = f"{YEARS[0]}-01-01",
    end: str = f"{YEARS[-1]}-12-31",
    budget: float = 2400.0,
    pause_s: float = 15.0,
    backoff_s: float = 1800.0,
    max_backoffs: int = 24,
    data_dir: Path = DATA_DIR,
) -> list[Point]:
    """Fetch every missing point; returns the points that are cached afterwards."""
    ledger = Ledger(data_dir / "openmeteo_ledger.json")
    meta_path = weather_dir(data_dir) / "grid_meta.json"
    meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
    weight = request_weight(len(VARIABLES), start, end)
    session = requests.Session()
    session.headers["Accept-Encoding"] = "gzip"
    backoffs = 0
    for p in fetch_order(points):
        out = cache_path(p, data_dir)
        if out.exists():
            continue
        while True:
            used = ledger.used()
            if used + weight > budget:
                log.warning("daily budget reached (%.0f of %.0f used in 24 h); stopping", used, budget)
                return [q for q in points if cache_path(q, data_dir).exists()]
            try:
                df, m = fetch_point(p, start, end, session)
            except RateLimited as exc:
                ledger.add(1.0)
                backoffs += 1
                if backoffs > max_backoffs:
                    log.error("rate limited too many times; giving up")
                    return [q for q in points if cache_path(q, data_dir).exists()]
                log.warning("429 from Open-Meteo (%s); sleeping %.0f s", exc, backoff_s)
                time.sleep(backoff_s)
                continue
            ledger.add(weight)
            break
        out.parent.mkdir(parents=True, exist_ok=True)
        df.to_parquet(out)
        meta[p.id] = {"lat": p.lat, "lon": p.lon, "grid": m}
        meta_path.write_text(json.dumps(meta, indent=1))
        log.info("cached %s: %d hours, weight %.1f (24 h total %.0f)", p.id, len(df), weight, ledger.used())
        time.sleep(pause_s)
    return [q for q in points if cache_path(q, data_dir).exists()]


def load_point(point: Point, data_dir: Path = DATA_DIR) -> pd.DataFrame:
    return pd.read_parquet(cache_path(point, data_dir))


def available_points(points: list[Point] = GRID, data_dir: Path = DATA_DIR) -> list[Point]:
    return [p for p in points if cache_path(p, data_dir).exists()]
