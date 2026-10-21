"""Paths and study-wide constants."""

from __future__ import annotations

import os
from pathlib import Path

DATA_DIR = Path(os.environ.get("RES_RELIABILITY_DATA", "data"))
RESULTS_DIR = Path("results")
DOCS_DIR = Path("docs")
IMG_DIR = DOCS_DIR / "img"

# Weather and load years (UTC calendar years). PSE demand is available for all of them;
# 2024-06-14 onwards comes from api.raporty.pse.pl, earlier hours from the legacy PSE export.
YEARS = list(range(2015, 2025))
HOURS = 8760  # one study year; the last day of leap years is dropped so every year has 8760 h

# Wake and electrical losses applied to the gross ERA5 + power-curve capacity factor
# (technical availability is modelled separately as structural reliability).
WIND_LOSSES = 0.10

# Polish reliability standard used for the capacity market (LOLE, h/year).
PL_RELIABILITY_STANDARD_LOLE = 3.0

# Literature structural-reliability parameters (forced outage rate, MTTR in hours) used in the
# thesis variant "OZE-LIT": wind FOR 5 % and PV FOR 3 % (Billinton & Allan 1996, Dhillon 1999).
LIT_WIND = {"for": 0.05, "mttr_h": 80.0}
LIT_PV = {"for": 0.03, "mttr_h": 24.0}


def kelmarsh_dir() -> Path:
    """Folder with the Kelmarsh SCADA zips (Zenodo record 16807551)."""
    default = Path("..") / "wind-power-curve-kelmarsh" / "data"
    return Path(os.environ.get("KELMARSH_DATA_DIR", default))
