"""Study grid over Poland and the country outline (Natural Earth, public domain)."""

from __future__ import annotations

import json
from dataclasses import dataclass
from importlib import resources

import numpy as np


@dataclass(frozen=True)
class Point:
    lat: float
    lon: float

    @property
    def id(self) -> str:
        return f"N{self.lat:05.2f}E{self.lon:05.2f}"


# Staggered grid: rows every 0.9 deg latitude, 1.5 deg longitude spacing with a half-step
# offset on odd rows, kept where the point and a ~15 km neighbourhood fall inside Poland.
GRID = [
    Point(49.5, 22.1),
    Point(50.4, 18.35), Point(50.4, 19.85), Point(50.4, 21.35), Point(50.4, 22.85),
    Point(51.3, 16.1), Point(51.3, 17.6), Point(51.3, 19.1), Point(51.3, 20.6), Point(51.3, 22.1),
    Point(52.2, 15.35), Point(52.2, 16.85), Point(52.2, 18.35), Point(52.2, 19.85),
    Point(52.2, 21.35), Point(52.2, 22.85),
    Point(53.1, 14.6), Point(53.1, 16.1), Point(53.1, 17.6), Point(53.1, 19.1),
    Point(53.1, 20.6), Point(53.1, 22.1), Point(53.1, 23.6),
    Point(54.0, 15.35), Point(54.0, 16.85), Point(54.0, 18.35), Point(54.0, 19.85),
    Point(54.0, 21.35), Point(54.0, 22.85),
]


def fetch_order(points: list[Point] = GRID) -> list[Point]:
    """Order points so that any prefix covers the country evenly (farthest-point sampling)."""
    xy = np.array([[p.lon * np.cos(np.radians(52.0)), p.lat] for p in points])
    order = [int(np.argmin(np.abs(xy - xy.mean(axis=0)).sum(axis=1)))]
    dist = np.linalg.norm(xy - xy[order[0]], axis=1)
    while len(order) < len(points):
        nxt = int(np.argmax(dist))
        order.append(nxt)
        dist = np.minimum(dist, np.linalg.norm(xy - xy[nxt], axis=1))
    return [points[i] for i in order]


def poland_outline() -> np.ndarray:
    """(n, 2) array of lon/lat vertices of mainland Poland."""
    text = resources.files("res_reliability").joinpath("data/poland_outline.json").read_text("utf-8")
    return np.asarray(json.loads(text)["lonlat"], dtype=float)
