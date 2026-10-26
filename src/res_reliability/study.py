"""Study context: grid series, portfolio, load profiles, conventional system, RES models.

Two load profiles are used, as in the thesis:

* IEEE RTS-79 (8736 h, peak 2850 MW) with its own 32-unit generating system;
* the real KSE demand of each weather year (PSE), normalised to its annual peak and scaled so
  that the RTS-79 generating system just meets the Polish reliability standard
  (LOLE = 3 h/year) on average over the ten years. The thesis instead scaled the 2023
  ENTSO-E profile to 2850 MW; that scaling is reported as well.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import brentq

from res_reliability import adequacy, elcc, markov as mk, montecarlo as mc, pse, regimes as rgm, rts79, sites, weather
from res_reliability.config import DATA_DIR, HOURS, PL_RELIABILITY_STANDARD_LOLE, YEARS

log = logging.getLogger(__name__)


@dataclass
class LoadCase:
    name: str
    matrix: np.ndarray  # (1 or n_years, n_hours) MW
    note: str

    @property
    def n_hours(self) -> int:
        return self.matrix.shape[1]

    @property
    def peak(self) -> float:
        return float(self.matrix.max())


@dataclass
class Context:
    sites: list[sites.SiteSeries]
    portfolio: sites.SiteSeries
    months: np.ndarray
    copt: adequacy.Copt
    loads: dict[str, LoadCase]
    kse_shape: np.ndarray | None = None
    t2m: np.ndarray | None = None  # (n_years, 8760) mean ERA5 temperature of the grid points
    extra: dict = field(default_factory=dict)


def rts_copt() -> adequacy.Copt:
    t = rts79.UNIT_TYPES
    return adequacy.Copt.from_units([u.capacity_mw for u in t for _ in range(u.count)],
                                    [u.forced_outage_rate for u in t for _ in range(u.count)])


def conv_arrays() -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    u = rts79.units()
    return (np.array([a for a, _, _ in u]), np.array([b for _, b, _ in u]), np.array([c for _, _, c in u]))


def kse_shape(data_dir=DATA_DIR) -> np.ndarray:
    s = pse.fetch(data_dir=data_dir)
    m = pse.year_matrix(s, YEARS)
    return m / m.max(axis=1, keepdims=True)


def scale_to_lole(copt: adequacy.Copt, shape: np.ndarray, target: float) -> float:
    f = lambda peak: elcc.lole(copt, shape * peak) - target  # noqa: E731
    return float(brentq(f, 1500.0, 3405.0, xtol=0.01))


def build(data_dir=DATA_DIR, with_kse: bool = True) -> Context:
    pts = weather.available_points(data_dir=data_dir)
    if not pts:
        raise RuntimeError("no weather cached; run fetch-weather first")
    ss = [sites.build(p, data_dir=data_dir) for p in pts]
    pf = sites.portfolio(ss)
    copt = rts_copt()
    loads = {"RTS": LoadCase("RTS", rts79.load_profile()[None], "IEEE RTS-79, 8736 h, peak 2850 MW")}
    shape = None
    if with_kse:
        shape = kse_shape(data_dir)
        peak = scale_to_lole(copt, shape, PL_RELIABILITY_STANDARD_LOLE)
        loads["KSE"] = LoadCase("KSE", shape * peak,
                                f"PSE KSE demand {YEARS[0]}-{YEARS[-1]}, each year scaled to peak {peak:.0f} MW "
                                f"(conventional-only LOLE = {PL_RELIABILITY_STANDARD_LOLE:g} h/yr)")
        loads["KSE2850"] = LoadCase("KSE2850", shape * rts79.PEAK_LOAD_MW,
                                    "PSE KSE demand scaled to the RTS-79 peak of 2850 MW (thesis scaling)")
    log.info("%d grid points, portfolio wind CF %.3f, PV CF %.3f", len(ss), pf.wind.mean(), pf.pv.mean())
    t2m = sites.mean_temperature(pts, data_dir=data_dir)
    return Context(ss, pf, sites.study_months(), copt, loads, shape, t2m)


# --- RES models -------------------------------------------------------------------------

def fit_wind(s: sites.SiteSeries, months: np.ndarray, n_states: int = 5, seasonal: bool = True,
             level_kind: str = "mean", rows: np.ndarray | None = None) -> mk.MarkovModel:
    x = s.wind if rows is None else s.wind[rows]
    m = np.tile(months, x.shape[0]) if seasonal else None
    return mk.fit(x.reshape(-1), n_states, months=m, level_kind=level_kind)


def fit_pv(s: sites.SiteSeries, months: np.ndarray, n_states: int = 5, seasonal: bool = True,
           level_kind: str = "mean", rows: np.ndarray | None = None) -> mk.MarkovModel:
    k, day = s.pv_k, s.daylight
    if rows is not None:
        k, day = k[rows], day[rows]
    m = np.tile(months, k.shape[0]) if seasonal else None
    return mk.fit(k.reshape(-1), n_states, months=m, valid=day.reshape(-1), level_kind=level_kind)


def pv_envelope(s: sites.SiteSeries) -> np.ndarray:
    return s.pv_cs * s.daylight


def dist_wind(ctx: Context, s: sites.SiteSeries, cap: float, model: str, n_hours: int, **fit_kw) -> elcc.ResDistribution:
    if model == "hist":
        return elcc.historical(s.wind, cap, n_hours)
    mdl = fit_wind(s, ctx.months, seasonal=(model == "seasonal"), **fit_kw)
    return elcc.from_markov(mdl, ctx.months, cap, None, s.wind.shape[0], n_hours)


def dist_pv(ctx: Context, s: sites.SiteSeries, cap: float, model: str, n_hours: int, **fit_kw) -> elcc.ResDistribution:
    if model == "hist":
        return elcc.historical(s.pv, cap, n_hours)
    mdl = fit_pv(s, ctx.months, seasonal=(model == "seasonal"), **fit_kw)
    return elcc.from_markov(mdl, ctx.months, cap, pv_envelope(s), s.pv.shape[0], n_hours)


REGIME_MODELS = ("regime", "regime_t")  # regime switching; "_t": conditioned on the temperature anomaly


def regime_weather(ctx: Context, model: str, **fit_kw) -> rgm.RegimeWeather:
    """Regime model of the national portfolio (cached on the context for the default settings)."""
    key = f"weather:{model}"
    if fit_kw or key not in ctx.extra:
        w = rgm.fit_weather(ctx.portfolio.wind, ctx.portfolio.pv_k, ctx.portfolio.daylight, ctx.months,
                            temperature=ctx.t2m if model == "regime_t" else None, **fit_kw)
        if fit_kw:
            return w
        ctx.extra[key] = w
    return ctx.extra[key]


def production(ctx: Context, s: sites.SiteSeries, tech: str, model: str, n_states: int = 5,
               level_kind: str = "mean") -> mc.Production:
    """Monte Carlo production model of a site for ``tech`` in {"wind", "pv"}.

    Regime models are fitted on the national portfolio; their plants need the shared regime process in the
    system specification (``system_spec(..., weather=regime_weather(ctx, model))``).
    """
    if model == "hist":
        return mc.HistoricalProduction(s.wind if tech == "wind" else s.pv)
    if model in REGIME_MODELS:
        w = regime_weather(ctx, model)
        if tech == "wind":
            return mc.RegimeProduction(w, "wind")
        return mc.RegimeProduction(w, "pv", envelope=pv_envelope(s), active=s.daylight)
    seasonal = model == "seasonal"
    if tech == "wind":
        mdl = fit_wind(s, ctx.months, n_states, seasonal, level_kind)
        return mc.MarkovProduction(mdl, ctx.months)
    mdl = fit_pv(s, ctx.months, n_states, seasonal, level_kind)
    return mc.MarkovProduction(mdl, ctx.months, envelope=pv_envelope(s), active=s.daylight)


def system_spec(ctx: Context, load: LoadCase, res: list[mc.ResPlant], weather: rgm.RegimeWeather | None = None,
                storage=None) -> mc.SystemSpec:
    cap, mttf, mttr = conv_arrays()
    if weather is None:
        weather = next((p.production.weather for p in res if isinstance(p.production, mc.RegimeProduction)), None)
    return mc.SystemSpec(load=load.matrix, conv_capacity=cap, conv_mttf=mttf, conv_mttr=mttr,
                         res=res, n_weather_years=len(YEARS), weather=weather, storage=storage)


def hours_of(load: LoadCase) -> int:
    return min(load.n_hours, HOURS)
