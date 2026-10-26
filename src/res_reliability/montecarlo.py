"""Sequential (chronological) Monte Carlo simulation of a generating system, hour by hour.

Every simulated year ``i`` draws all of its random numbers from
``SeedSequence(seed, spawn_key=(i,))`` spawned into one child stream per component, so a
year is reproducible on its own, independent of chunking or of which other components are
present (common random numbers between scenarios), which keeps parallel runs reproducible
year by year.

Conventional units are two-state Markov processes with exponential up and down times
(MTTF, MTTR), initialised from their stationary probabilities. RES output comes from a
production model (Markov chain, regime model or a historical weather year) times the share of
RES units that are technically available (structural reliability). A regime model draws one
daily regime path per simulated year from its own stream, shared by all RES plants. Storage,
if present, operates on the resulting margins (``storage``).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Protocol

import numpy as np

from res_reliability import markov as mk
from res_reliability import regimes as rgm
from res_reliability import storage as sto

CONV_STREAM = 0
WEATHER_STREAM = 999  # fixed child index: the regime path does not depend on how many plants there are


def year_streams(seed: int, sim: int, n_streams: int) -> list[np.random.Generator]:
    ss = np.random.SeedSequence(seed, spawn_key=(sim,))
    return [np.random.default_rng(s) for s in ss.spawn(n_streams)]


def weather_stream(seed: int, sim: int) -> np.random.Generator:
    """Child ``WEATHER_STREAM`` of the year's seed sequence (same as spawning that many children)."""
    return np.random.default_rng(np.random.SeedSequence(seed, spawn_key=(sim, WEATHER_STREAM)))


def two_state_availability(
    rng: np.random.Generator,
    capacity: np.ndarray,
    mttf: np.ndarray,
    mttr: np.ndarray,
    n_hours: int,
) -> np.ndarray:
    """Available capacity per hour of independent two-state units (sampled at hour starts)."""
    capacity = np.asarray(capacity, dtype=float)
    mttf = np.asarray(mttf, dtype=float)
    mttr = np.asarray(mttr, dtype=float)
    n = len(capacity)
    q = mttr / (mttf + mttr)
    up0 = rng.random(n) >= q
    diff = np.zeros(n_hours + 1)
    start = np.where(up0, capacity, 0.0).sum()
    k = int(np.ceil(4 * n_hours / np.min(mttf + mttr))) + 8
    while True:
        e = rng.standard_exponential((n, 2 * k))
        # durations alternate starting with the initial state's mean (memoryless residual)
        first_up = up0[:, None]
        col_even = np.arange(2 * k) % 2 == 0
        mean_dur = np.where(col_even[None, :] == first_up, mttf[:, None], mttr[:, None])
        t = np.cumsum(e * mean_dur, axis=1)
        if np.all(t[:, -1] >= n_hours):
            break
        k *= 2
    hit = t < n_hours
    rows, cols = np.nonzero(hit)
    # transition number c moves a unit out of the state it started column c in
    leaving_up = (cols % 2 == 0) == up0[rows]
    delta = np.where(leaving_up, -capacity[rows], capacity[rows])
    np.add.at(diff, np.ceil(t[rows, cols]).astype(np.int64), delta)
    return start + np.cumsum(diff[:n_hours])


def plant_bernoulli_availability(rng: np.random.Generator, forced_outage_rate: float, n_hours: int) -> np.ndarray:
    """Thesis implementation: the whole plant is out with probability FOR, independently every hour."""
    return (rng.random(n_hours) >= forced_outage_rate).astype(float)


@dataclass
class Structural:
    """Technical availability of a RES plant made of ``n_units`` identical units."""

    forced_outage_rate: float
    mttr_h: float
    n_units: int = 100
    # "units": two-state process with MTTF/MTTR per unit; "hourly_bernoulli": each unit out with
    # probability FOR independently every hour (thesis, applied per farm); "plant_bernoulli": the
    # same for the whole plant at once; "none"
    mode: str = "units"

    def sample(self, rng: np.random.Generator, n_hours: int) -> np.ndarray:
        if self.mode == "none" or self.forced_outage_rate <= 0:
            return np.ones(n_hours)
        if self.mode == "plant_bernoulli":
            return plant_bernoulli_availability(rng, self.forced_outage_rate, n_hours)
        if self.mode == "hourly_bernoulli":
            return rng.binomial(self.n_units, 1.0 - self.forced_outage_rate, n_hours) / self.n_units
        mttf = self.mttr_h * (1 - self.forced_outage_rate) / self.forced_outage_rate
        cap = np.full(self.n_units, 1.0 / self.n_units)
        return two_state_availability(rng, cap, np.full(self.n_units, mttf), np.full(self.n_units, self.mttr_h), n_hours)


class Production(Protocol):
    def sample(self, rngs: list[np.random.Generator], years: np.ndarray, n_hours: int) -> np.ndarray:
        """(n_sims, n_hours) normalised output for the simulated years in ``years``."""


@dataclass
class MarkovProduction:
    """Output = level(state) x envelope(hour); the chain moves only on ``active`` hours."""

    model: mk.MarkovModel
    months: np.ndarray  # month (1..12) of every hour of the study year
    envelope: np.ndarray | None = None  # (n_hours,) or (n_years, n_hours); None = 1
    active: np.ndarray | None = None  # (n_hours,) or (n_years, n_hours)

    def sample(self, rngs, years, n_hours):
        u = np.vstack([r.random(n_hours + 1) for r in rngs])
        env = None if self.envelope is None else np.atleast_2d(self.envelope)[:, :n_hours]
        act = None if self.active is None else np.atleast_2d(self.active)[:, :n_hours]
        out = np.empty((len(rngs), n_hours))
        groups = [np.arange(len(rngs))] if act is None or act.shape[0] == 1 else \
            [np.flatnonzero(years % act.shape[0] == y) for y in range(act.shape[0])]
        for g in groups:
            if len(g) == 0:
                continue
            a = None if act is None else act[int(years[g[0]]) % act.shape[0]]
            states = mk.simulate_states(self.model, u[g], self.months, a)
            if self.model.levels.ndim == 2:
                lv = self.model.levels[self.months - 1]  # (n_hours, S)
                o = lv[np.arange(n_hours)[None, :], states]
            else:
                o = self.model.levels[states]
            if env is not None:
                o = o * env[years[g] % env.shape[0]]
            if a is not None:
                o = np.where(a[None, :], o, 0.0)
            out[g] = o
        return np.clip(out, 0.0, 1.0)


@dataclass
class HistoricalProduction:
    """Deterministic weather year: simulated year i uses weather year ``years[i]``."""

    series: np.ndarray  # (n_years, n_hours)

    def sample(self, rngs, years, n_hours):
        return np.asarray(self.series)[years % len(self.series), :n_hours]


@dataclass
class RegimeProduction:
    """Output of one technology of a regime model, driven by the year's shared regime path.

    Wind: capacity factor. PV: clear-sky index x ``envelope`` (clear-sky output of the weather year), moving
    only on ``active`` (daylight) hours. ``envelope`` and ``active`` are (n_years, 8760); the burn-in uses the
    end of the previous weather year.
    """

    weather: rgm.RegimeWeather
    tech: str  # "wind" or "pv"
    envelope: np.ndarray | None = None
    active: np.ndarray | None = None
    needs_regimes = True

    def sample(self, rngs, years, n_hours, regimes=None):
        if regimes is None:
            raise ValueError("RegimeProduction needs the shared regime paths (SystemSpec.weather)")
        b = self.weather.burn_in_days * 24
        classes = self.weather.hourly_classes(regimes, years)
        if self.active is None:
            x = self.weather.hourly(self.tech, rngs, classes)
            return np.clip(x[:, b:b + n_hours], 0.0, 1.0)
        out = np.empty((len(rngs), n_hours))
        ny = self.active.shape[0]
        for y in np.unique(years % ny):
            g = np.flatnonzero(years % ny == y)
            act = np.r_[self.active[(y - 1) % ny, -b:], self.active[y]]
            x = self.weather.hourly(self.tech, [rngs[i] for i in g], classes[g], act)[:, b:b + n_hours]
            out[g] = np.where(self.active[y, :n_hours], x * self.envelope[y, :n_hours], 0.0)
        return np.clip(out, 0.0, 1.0)


@dataclass
class ResPlant:
    name: str
    capacity_mw: float
    production: Production
    structural: Structural | None = None


@dataclass
class SystemSpec:
    load: np.ndarray  # (n_hours,) or (n_years, n_hours)
    conv_capacity: np.ndarray
    conv_mttf: np.ndarray
    conv_mttr: np.ndarray
    res: list[ResPlant] = field(default_factory=list)
    n_weather_years: int = 1
    weather: rgm.RegimeWeather | None = None  # shared regime process of RegimeProduction plants
    storage: sto.Storage | None = None


@dataclass
class McResult:
    lole: np.ndarray  # per simulated year
    eens: np.ndarray
    lolf: np.ndarray  # loss-of-load events per year
    margin: np.ndarray | None = None  # (n_sims, n_hours) float32 available - load (before storage)
    res: np.ndarray | None = None  # (n_sims, n_hours) float32 total RES output, MW

    def summary(self) -> dict:
        n = len(self.lole)
        out = {"n_sims": n}
        for k in ("lole", "eens", "lolf"):
            v = getattr(self, k)
            out[k] = float(v.mean())
            out[f"{k}_std"] = float(v.std(ddof=1))
            out[f"{k}_ci95"] = float(1.96 * v.std(ddof=1) / np.sqrt(n))
        dur = self.lole.sum() / max(self.lolf.sum(), 1)
        out["lold_h"] = float(dur)
        return out


def run(spec: SystemSpec, n_sims: int, seed: int = 0, chunk: int = 250, keep_margin: bool = False,
        keep_res: bool = False) -> McResult:
    """Simulate ``n_sims`` years. Year ``i`` uses weather/load year ``i % n_weather_years``.

    With ``spec.storage`` the indices are those after storage; ``margin`` is always the margin before it.
    """
    load = np.atleast_2d(spec.load)
    n_hours = load.shape[1]
    n_streams = 1 + 2 * len(spec.res)
    st = spec.storage
    lole, eens, lolf, margins, res_out = [], [], [], [], []
    for lo in range(0, n_sims, chunk):
        sims = np.arange(lo, min(lo + chunk, n_sims))
        streams = [year_streams(seed, int(i), n_streams) for i in sims]
        years = sims % max(spec.n_weather_years, 1)
        regimes = None
        if spec.weather is not None:
            regimes = spec.weather.regime_paths([weather_stream(seed, int(i)) for i in sims], years)
        gen = np.vstack([
            two_state_availability(s[CONV_STREAM], spec.conv_capacity, spec.conv_mttf, spec.conv_mttr, n_hours)
            for s in streams
        ])
        need_res = keep_res or (st is not None and st.charge_from == "res")
        res_gen = np.zeros_like(gen) if need_res else None
        for j, plant in enumerate(spec.res):
            rngs = [s[1 + 2 * j] for s in streams]
            if getattr(plant.production, "needs_regimes", False):
                prod = plant.production.sample(rngs, years, n_hours, regimes)
            else:
                prod = plant.production.sample(rngs, years, n_hours)
            if plant.structural is not None:
                avail = np.vstack([plant.structural.sample(s[2 + 2 * j], n_hours) for s in streams])
                prod = prod * avail
            gen += plant.capacity_mw * prod
            if need_res:
                res_gen += plant.capacity_mw * prod
        m = gen - load[years % load.shape[0]]
        if st is None:
            short = m < 0
            lole.append(short.sum(axis=1))
            eens.append(np.where(short, -m, 0.0).sum(axis=1))
            starts = short[:, 1:] & ~short[:, :-1]
            lolf.append(starts.sum(axis=1) + short[:, 0])
        else:
            r = sto.dispatch(m, st.power_mw, st.energy_mwh, st.efficiency, st.soc0,
                             charge_cap=res_gen if st.charge_from == "res" else None, warmup_h=st.warmup_h)
            lole.append(r["lole"][0])
            eens.append(r["eens"][0])
            lolf.append(r["lolf"][0])
        if keep_margin:
            margins.append(m.astype(np.float32))
        if keep_res:
            res_out.append(res_gen.astype(np.float32))
    return McResult(
        np.concatenate(lole).astype(float),
        np.concatenate(eens),
        np.concatenate(lolf).astype(float),
        np.vstack(margins) if keep_margin else None,
        np.vstack(res_out) if keep_res else None,
    )
