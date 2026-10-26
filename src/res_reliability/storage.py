"""Energy storage in the sequential Monte Carlo: greedy operation on hourly margins.

The storage only reacts to the system: it discharges when available generation is below the load (as much as
its power and state of charge allow) and charges from any surplus (or, for a co-located plant, only from the
RES output) until it is full. It never charges during a shortfall, so it cannot cause one. With round-trip
efficiency eta split equally between charging and discharging,

    soc[t+1] = soc[t] + sqrt(eta) * charge[t] - discharge[t] / sqrt(eta),    0 <= soc <= E,
    charge, discharge <= P.

The margins of a Monte Carlo run do not depend on the storage, so one run serves any number of storage
configurations, and ELCC or the energy needed for a reliability target are found by bisection on the stored
margins (common random numbers). An optional warm-up replays the last ``warmup_h`` hours of the same simulated
year first, so the state of charge on 1 January comes from the system's own operation; loss of load is only
counted inside the year.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

TOL_MW = 1e-6


@dataclass
class Storage:
    power_mw: float
    energy_mwh: float
    efficiency: float = 1.0  # round trip
    soc0: float = 1.0  # state of charge at the start of the simulation, share of the energy
    charge_from: str = "any"  # "any" surplus, or "res" output only (co-located)
    warmup_h: int = 0

    @property
    def duration_h(self) -> float:
        return self.energy_mwh / self.power_mw if self.power_mw > 0 else float("inf")


def _cfg(*arrays):
    a = np.broadcast_arrays(*[np.atleast_1d(np.asarray(x, dtype=float)) for x in arrays])
    return [x[:, None] for x in a]


def dispatch(
    margin: np.ndarray,
    power,
    energy,
    efficiency=1.0,
    soc0=1.0,
    shift=0.0,
    charge_cap: np.ndarray | None = None,
    warmup_h: int = 0,
    keep: bool = False,
) -> dict:
    """Operate storage configurations on ``margin`` (n_sims, n_hours) = available generation - load, MW.

    ``power``, ``energy``, ``efficiency``, ``soc0`` and ``shift`` (a constant load added to every hour) are
    scalars or arrays of one value per configuration; ``power`` may be ``inf``. ``charge_cap`` (n_sims, n_hours)
    limits charging (co-located storage: the RES output). Returns per configuration and simulated year
    ``lole`` (h), ``eens`` (MWh), ``lolf`` (events), ``charged`` and ``discharged`` (MWh, grid side), each
    (n_cfg, n_sims); with ``keep`` also ``soc`` (state of charge at the end of every hour), ``charge``,
    ``discharge`` and ``net`` (margin after storage), each (n_cfg, n_sims, n_hours).
    """
    n_sims, n_hours = np.shape(margin)
    mT = np.ascontiguousarray(np.asarray(margin, dtype=float).T)  # hour rows: contiguous per-hour slices
    capT = None if charge_cap is None else np.ascontiguousarray(np.asarray(charge_cap, dtype=float).T)
    P, E, eta, s0, sh = _cfg(power, energy, efficiency, soc0, shift)
    n_cfg = P.shape[0]
    shape = (n_cfg, n_sims)
    ec = np.broadcast_to(np.sqrt(eta), shape).copy()
    P, E, sh = [np.broadcast_to(x, shape).copy() for x in (P, E, sh)]
    soc = np.broadcast_to(s0, shape) * E
    full_level = E - 1e-9 * np.maximum(E, 1.0)
    lole, eens, lolf, chg, dis_tot = (np.zeros(shape) for _ in range(5))
    m, sur, room, charge, deficit, discharge, tmp = (np.empty(shape) for _ in range(7))
    short, prev, starts = (np.zeros(shape, dtype=bool) for _ in range(3))
    if keep:
        rec = {k: np.empty((n_cfg, n_sims, n_hours)) for k in ("soc", "charge", "discharge", "net")}
    order = np.r_[np.arange(n_hours - warmup_h, n_hours), np.arange(n_hours)].astype(np.int64)
    # Hours where no simulated year can be short at any of the shifts. When every storage is full there,
    # nothing happens (no discharge, nothing to charge), so they are skipped; the result is exact.
    cand = (mT < sh.max()).any(axis=1)[order]
    nxt = np.full(len(order) + 1, len(order))
    idx = np.flatnonzero(cand)
    nxt[: len(order)] = np.r_[idx, len(order)][np.searchsorted(idx, np.arange(len(order)))]
    k = 0
    while k < len(order):
        t = order[k]
        np.subtract(mT[t], sh, out=m)
        np.maximum(m, 0.0, out=sur)
        if capT is not None:
            np.minimum(sur, capT[t], out=sur)
        np.minimum(sur, P, out=sur)
        np.subtract(E, soc, out=room)
        np.divide(room, ec, out=room)
        np.minimum(sur, room, out=charge)
        np.negative(m, out=deficit)
        np.maximum(deficit, 0.0, out=deficit)
        np.minimum(deficit, P, out=discharge)
        np.multiply(soc, ec, out=tmp)
        np.minimum(discharge, tmp, out=discharge)
        np.multiply(charge, ec, out=tmp)
        soc += tmp
        np.divide(discharge, ec, out=tmp)
        soc -= tmp
        np.clip(soc, 0.0, E, out=soc)
        if k >= warmup_h:
            np.subtract(deficit, discharge, out=tmp)
            np.greater(tmp, TOL_MW, out=short)
            lole += short
            tmp *= short
            eens += tmp
            np.logical_and(short, ~prev, out=starts)
            lolf += starts
            prev[...] = short
            chg += charge
            dis_tot += discharge
            if keep:
                rec["soc"][:, :, t] = soc
                rec["charge"][:, :, t] = charge
                rec["discharge"][:, :, t] = discharge
                rec["net"][:, :, t] = m - charge + discharge
        k += 1
        if not keep and k < len(order) and not cand[k] and (soc >= full_level).all():
            k = nxt[k]
            prev[...] = False  # no shortfall in the skipped hours
    out = {"lole": lole, "eens": eens, "lolf": lolf, "charged": chg, "discharged": dis_tot}
    if keep:
        out.update(rec)
    return out


def _relevant(margin: np.ndarray, max_shift: float) -> np.ndarray:
    """Simulated years that can have a shortfall at some shift <= max_shift (the others never use storage)."""
    return np.asarray(margin).min(axis=1) < max_shift


def lole(margin: np.ndarray, storage: Storage, shift: float = 0.0, res: np.ndarray | None = None) -> float:
    """Expected LOLE (h/year) of one storage configuration."""
    cap = res if storage.charge_from == "res" else None
    r = dispatch(margin, storage.power_mw, storage.energy_mwh, storage.efficiency, storage.soc0, shift, cap,
                 storage.warmup_h)
    return float(r["lole"].mean())


def elcc(
    margin: np.ndarray,
    power,
    energy,
    efficiency=1.0,
    soc0=1.0,
    charge_cap: np.ndarray | None = None,
    warmup_h: int = 0,
    target: float | None = None,
    iters: int = 10,
) -> np.ndarray:
    """ELCC (MW) of storage configurations: the constant load added to every hour at which LOLE with the
    storage equals ``target`` (default: LOLE of the margins without storage). Bisection on [0, power] for all
    configurations at once; precision power / 2**iters."""
    margin = np.asarray(margin, dtype=float)
    P, E, eta, s0 = [x[:, 0] for x in _cfg(power, energy, efficiency, soc0)]
    n_sims = margin.shape[0]
    target = float((margin < 0).sum(axis=1).mean()) if target is None else target
    rel = _relevant(margin, float(P.max()))
    m = margin[rel]
    cap = None if charge_cap is None else np.asarray(charge_cap)[rel]
    lo, hi = np.zeros_like(P), P.copy()
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        r = dispatch(m, P, E, eta, s0, mid, cap, warmup_h)
        ok = r["lole"].sum(axis=1) / n_sims < target
        lo = np.where(ok, mid, lo)
        hi = np.where(ok, hi, mid)
    return 0.5 * (lo + hi)


def energy_needed(
    margin: np.ndarray,
    target: float,
    duration_h,
    efficiency=1.0,
    soc0=1.0,
    warmup_h: int = 0,
    e_min: float = 1.0,
    e_max: float = 1e6,
    iters: int = 14,
) -> np.ndarray:
    """Smallest storage energy (MWh) that brings LOLE down to ``target``, for every configuration.

    ``duration_h`` sets the power as energy / duration (``inf`` = power not limited). Bisection on log(E)
    between ``e_min`` and ``e_max``; returns ``nan`` where even ``e_max`` is not enough and 0 where no storage
    is needed.
    """
    margin = np.asarray(margin, dtype=float)
    dur, eta, s0 = [x[:, 0] for x in _cfg(duration_h, efficiency, soc0)]
    n_sims = margin.shape[0]

    def lole_at(E):
        Pw = np.where(np.isinf(dur), np.inf, E / dur)
        rel = _relevant(margin, 0.0)
        r = dispatch(margin[rel], Pw, E, eta, s0, 0.0, None, warmup_h)
        return r["lole"].sum(axis=1) / n_sims

    none_needed = float((margin < 0).sum(axis=1).mean()) <= target
    lo = np.full(len(dur), np.log(e_min))
    hi = np.full(len(dur), np.log(e_max))
    reach = lole_at(np.exp(hi)) <= target
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        ok = lole_at(np.exp(mid)) <= target
        lo = np.where(ok, lo, mid)
        hi = np.where(ok, mid, hi)
    out = np.exp(hi)
    out = np.where(reach, out, np.nan)
    return np.where(none_needed, 0.0, out)
