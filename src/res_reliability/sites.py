"""Per-point hourly wind and PV series arranged as (year, hour-of-year) matrices.

All series are labelled by the start of the hour in UTC. Open-Meteo/ERA5 radiation is the
mean of the preceding hour, so it is shifted back by one hour; wind and temperature are
instantaneous values at the top of the hour.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from res_reliability import power, weather
from res_reliability.config import DATA_DIR, HOURS, WIND_LOSSES, YEARS
from res_reliability.geo import Point

log = logging.getLogger(__name__)


@dataclass
class SiteSeries:
    point: Point
    wind: np.ndarray  # (n_years, 8760) capacity factor after wake/electrical losses
    pv: np.ndarray  # (n_years, 8760) AC output per unit of inverter capacity
    pv_cs: np.ndarray  # (n_years, 8760) clear-sky AC output
    ws: np.ndarray  # (n_years, 8760) ERA5 100 m wind speed, m/s

    @property
    def pv_k(self) -> np.ndarray:
        k, _ = power.clear_sky_index(self.pv, self.pv_cs)
        return k

    @property
    def daylight(self) -> np.ndarray:
        return self.pv_cs > 0.02


def year_index(y: int) -> pd.DatetimeIndex:
    return pd.date_range(f"{y}-01-01", periods=HOURS, freq="1h", tz="UTC")


def study_months(year: int = 2015) -> np.ndarray:
    """Month (1..12) of each of the 8760 study hours (a non-leap reference year)."""
    return year_index(year).month.to_numpy()


def series_path(point: Point, data_dir: Path = DATA_DIR) -> Path:
    return data_dir / "series" / f"{point.id}.npz"


def build(point: Point, years: list[int] = YEARS, data_dir: Path = DATA_DIR) -> SiteSeries:
    path = series_path(point, data_dir)
    if path.exists():
        z = np.load(path)
        return SiteSeries(point, z["wind_gross"] * (1 - WIND_LOSSES), z["pv"], z["pv_cs"], z["ws"])
    w = weather.load_point(point, data_dir)
    w = w[~w.index.duplicated()]
    ghi = w["ghi"].shift(-1)  # hour-mean labelled by its end -> label by its start
    wind, pv, pv_cs, ws = [], [], [], []
    for y in years:
        idx = year_index(y)
        v = w["ws100"].reindex(idx).interpolate(limit=3, limit_direction="both").to_numpy()
        g = ghi.reindex(idx).interpolate(limit=3, limit_direction="both").fillna(0.0).to_numpy()
        t = w["t2m"].reindex(idx).interpolate(limit=3, limit_direction="both").to_numpy()
        if np.isnan(v).any() or np.isnan(t).any():
            raise ValueError(f"{point.id} {y}: gaps in ERA5 series")
        p, p_cs = power.pv_output(idx, point.lat, point.lon, g, t)
        wind.append(power.wind_cf(v))
        pv.append(p)
        pv_cs.append(p_cs)
        ws.append(v)
    gross = np.vstack(wind)
    s = SiteSeries(point, gross * (1 - WIND_LOSSES), np.vstack(pv), np.vstack(pv_cs), np.vstack(ws))
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, wind_gross=gross, pv=s.pv, pv_cs=s.pv_cs, ws=s.ws)
    log.info("%s: wind CF %.3f, PV CF %.3f", point.id, s.wind.mean(), s.pv.mean())
    return s


def portfolio(sites: list[SiteSeries]) -> SiteSeries:
    """Equal capacity at every point: the distributed national portfolio."""
    pt = Point(float(np.mean([s.point.lat for s in sites])), float(np.mean([s.point.lon for s in sites])))
    return SiteSeries(
        pt,
        np.mean([s.wind for s in sites], axis=0),
        np.mean([s.pv for s in sites], axis=0),
        np.mean([s.pv_cs for s in sites], axis=0),
        np.mean([s.ws for s in sites], axis=0),
    )
