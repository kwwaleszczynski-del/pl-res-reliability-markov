"""IEEE Reliability Test System 1979: generating units and the chronological load model.

Source: IEEE Subcommittee on the Application of Probability Methods, "IEEE Reliability Test
System", IEEE Trans. PAS-98(6), 1979, pp. 2047-2054 (Tables 1-3 and 6). With these data the
hourly-load LOLE of the generating system is 9.394 h/year and LOEE 1176 MWh/year
(Billinton & Allan, Reliability Evaluation of Power Systems, 2nd ed., 1996), which the tests
reproduce analytically.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

PEAK_LOAD_MW = 2850.0
HOURS = 8736  # 52 weeks x 7 days x 24 h


@dataclass(frozen=True)
class UnitType:
    capacity_mw: float
    count: int
    mttf_h: float
    mttr_h: float
    fuel: str

    @property
    def forced_outage_rate(self) -> float:
        return self.mttr_h / (self.mttf_h + self.mttr_h)


# Table 6: unit size, number, MTTF, MTTR (FOR = MTTR / (MTTF + MTTR)).
UNIT_TYPES = [
    UnitType(12, 5, 2940, 60, "oil/steam"),
    UnitType(20, 4, 450, 50, "oil/CT"),
    UnitType(50, 6, 1980, 20, "hydro"),
    UnitType(76, 4, 1960, 40, "coal/steam"),
    UnitType(100, 3, 1200, 50, "oil/steam"),
    UnitType(155, 4, 960, 40, "coal/steam"),
    UnitType(197, 3, 950, 50, "oil/steam"),
    UnitType(350, 1, 1150, 100, "coal/steam"),
    UnitType(400, 2, 1100, 150, "nuclear"),
]

# Table 1: weekly peak load, % of annual peak (week 1 = first week of January).
WEEKLY = np.array([
    86.2, 90.0, 87.8, 83.4, 88.0, 84.1, 83.2, 80.6, 74.0, 73.7, 71.5, 72.7, 70.4,
    75.0, 72.1, 80.0, 75.4, 83.7, 87.0, 88.0, 85.6, 81.1, 90.0, 88.7, 89.6, 86.1,
    75.5, 81.6, 80.1, 88.0, 72.2, 77.6, 80.0, 72.9, 72.6, 70.5, 78.0, 69.5, 72.4,
    72.4, 74.3, 74.4, 80.0, 88.1, 88.5, 90.9, 94.0, 89.0, 94.2, 97.0, 100.0, 95.2,
]) / 100.0

# Table 2: daily peak load, % of weekly peak, Monday..Sunday.
DAILY = np.array([93, 100, 98, 96, 94, 77, 75]) / 100.0

# Table 3: hourly load, % of daily peak, hour 0 = 00:00-01:00.
HOURLY = {
    ("winter", "weekday"): [67, 63, 60, 59, 59, 60, 74, 86, 95, 96, 96, 95,
                            95, 95, 93, 94, 99, 100, 100, 96, 91, 83, 73, 63],
    ("winter", "weekend"): [78, 72, 68, 66, 64, 65, 66, 70, 80, 88, 90, 91,
                            90, 88, 87, 87, 91, 100, 99, 97, 94, 92, 87, 81],
    ("summer", "weekday"): [64, 60, 58, 56, 56, 58, 64, 76, 87, 95, 99, 100,
                            99, 100, 100, 97, 96, 96, 93, 92, 92, 93, 87, 72],
    ("summer", "weekend"): [74, 70, 66, 65, 64, 62, 62, 66, 81, 86, 91, 93,
                            93, 92, 91, 91, 92, 94, 95, 95, 100, 93, 88, 80],
    ("spring_fall", "weekday"): [63, 62, 60, 58, 59, 65, 72, 85, 95, 99, 100, 99,
                                 93, 92, 90, 88, 90, 92, 96, 98, 96, 90, 80, 70],
    ("spring_fall", "weekend"): [75, 73, 69, 66, 65, 65, 68, 74, 83, 89, 92, 94,
                                 91, 90, 90, 86, 85, 88, 92, 100, 97, 95, 90, 85],
}


def season_of_week(week: int) -> str:
    """Season for a 1-based week number (winter 1-8 and 44-52, summer 18-30)."""
    if week <= 8 or week >= 44:
        return "winter"
    if 18 <= week <= 30:
        return "summer"
    return "spring_fall"


def load_profile(peak_mw: float = PEAK_LOAD_MW) -> np.ndarray:
    """Hourly load (MW) for the 8736-hour RTS year."""
    out = np.empty(HOURS)
    h = 0
    for w in range(52):
        season = season_of_week(w + 1)
        for d in range(7):
            kind = "weekend" if d >= 5 else "weekday"
            day_peak = peak_mw * WEEKLY[w] * DAILY[d]
            out[h:h + 24] = day_peak * np.asarray(HOURLY[(season, kind)]) / 100.0
            h += 24
    return out


def units() -> list[tuple[float, float, float]]:
    """One (capacity MW, MTTF h, MTTR h) tuple per generating unit (32 units, 3405 MW)."""
    return [(t.capacity_mw, t.mttf_h, t.mttr_h) for t in UNIT_TYPES for _ in range(t.count)]


def installed_capacity() -> float:
    return sum(t.capacity_mw * t.count for t in UNIT_TYPES)
