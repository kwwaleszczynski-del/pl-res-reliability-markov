"""Weather -> normalised power: wind capacity factor and PV output with its clear-sky envelope.

Wind: measured IEC binned power curve of the Senvion MM92 (2.05 MW, 92 m rotor) from the
open Kelmarsh SCADA data (2017 reference year, mean of the six turbines; Plumley 2022,
CC BY 4.0), smoothed with a Gaussian kernel to represent a farm spread over an ERA5 cell.

PV: fixed-tilt, south-facing plant. ERA5 GHI (hour-average, labelled by hour start) is split
into beam and diffuse with the Erbs model, transposed with the isotropic sky model, derated
for cell temperature (Faiman, PVWatts temperature coefficient) and clipped at the inverter.
The same chain applied to Haurwitz clear-sky GHI gives the clear-sky envelope used by the
PV Markov model (state = clear-sky index of PV output).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pvlib

# Kelmarsh MM92 reference curve (bin centre m/s, power kW); rated 2050 kW above 16.5 m/s.
MM92_WS = np.array([0.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0,
                    9.5, 10.0, 10.5, 11.0, 11.5, 12.0, 12.5, 13.0, 13.5, 14.0])
MM92_KW = np.array([0.0, 0.0, 20.3, 59.3, 106.9, 163.6, 233.2, 320.0, 423.5, 548.2, 692.1, 855.4,
                    1033.0, 1225.3, 1421.0, 1610.4, 1771.6, 1887.1, 1963.6, 2006.5, 2027.7,
                    2038.8, 2044.1, 2048.1, 2050.0])
MM92_RATED_KW = 2050.0
CUT_OUT_MS = 24.0


def turbine_curve(ws: np.ndarray) -> np.ndarray:
    """Single-turbine normalised power curve (0..1)."""
    ws = np.asarray(ws, dtype=float)
    p = np.interp(ws, MM92_WS, MM92_KW / MM92_RATED_KW, left=0.0, right=1.0)
    return np.where(ws >= CUT_OUT_MS, 0.0, p)


def farm_curve(ws: np.ndarray, sigma0: float = 0.6, sigma1: float = 0.2) -> np.ndarray:
    """Turbine curve convolved with N(v, sigma0 + sigma1 * v): spatial spread inside a cell."""
    grid = np.arange(0.0, 35.0, 0.05)
    sig = sigma0 + sigma1 * grid
    z = np.linspace(-3, 3, 61)
    w = np.exp(-0.5 * z**2)
    w /= w.sum()
    smooth = (turbine_curve(grid[:, None] + sig[:, None] * z[None, :]) * w).sum(axis=1)
    return np.interp(np.asarray(ws, dtype=float), grid, smooth)


def wind_cf(ws100: np.ndarray, losses: float = 0.0) -> np.ndarray:
    return np.clip(farm_curve(ws100) * (1.0 - losses), 0.0, 1.0)


def pv_output(
    times: pd.DatetimeIndex,
    lat: float,
    lon: float,
    ghi: np.ndarray,
    t_air: np.ndarray,
    tilt: float = 35.0,
    azimuth: float = 180.0,
    dc_ac: float = 1.2,
    system_losses: float = 0.10,
    gamma: float = -0.0037,
) -> tuple[np.ndarray, np.ndarray]:
    """Normalised AC output (per unit of inverter capacity) and its clear-sky counterpart.

    ``times`` label the start of each hour; solar geometry is taken at mid-hour.
    """
    mid = times + pd.Timedelta(minutes=30)
    sp = pvlib.solarposition.get_solarposition(mid, lat, lon)
    zen = sp["apparent_zenith"].to_numpy()
    azi = sp["azimuth"].to_numpy()
    extra = pvlib.irradiance.get_extra_radiation(mid).to_numpy()
    cs_ghi = pvlib.clearsky.haurwitz(sp["apparent_zenith"])["ghi"].to_numpy()

    def chain(g: np.ndarray) -> np.ndarray:
        g = np.clip(np.nan_to_num(g), 0.0, None)
        erbs = pvlib.irradiance.erbs(g, zen, mid)
        dni = np.nan_to_num(np.asarray(erbs["dni"]))
        dhi = np.nan_to_num(np.asarray(erbs["dhi"]))
        poa = pvlib.irradiance.get_total_irradiance(
            tilt, azimuth, zen, azi, dni, g, dhi, dni_extra=extra, model="isotropic", albedo=0.2
        )["poa_global"]
        poa = np.clip(np.nan_to_num(np.asarray(poa)), 0.0, None)
        t_cell = pvlib.temperature.faiman(poa, t_air, wind_speed=1.0)
        dc = poa / 1000.0 * (1.0 + gamma * (t_cell - 25.0)) * (1.0 - system_losses) * dc_ac
        return np.clip(dc * 0.97, 0.0, 1.0)

    return chain(ghi), chain(cs_ghi)


def clear_sky_index(p: np.ndarray, p_cs: np.ndarray, threshold: float = 0.02) -> tuple[np.ndarray, np.ndarray]:
    """k = p / p_cs on daylight hours (p_cs above ``threshold``), clipped to [0, 1.5]."""
    day = p_cs > threshold
    k = np.where(day, np.divide(p, p_cs, out=np.zeros_like(p), where=day), np.nan)
    return np.clip(k, 0.0, 1.5), day
