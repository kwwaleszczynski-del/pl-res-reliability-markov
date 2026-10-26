"""Study stages. Each stage writes JSON/CSV into ``results/`` and returns what it computed.

Usage:
    python -m res_reliability run --stages markov structural variants elcc sensitivity
    python -m res_reliability run --stages weather storage --jobs 6
"""

from __future__ import annotations

import json
import logging
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.optimize import brentq

from res_reliability import config, elcc, markov as mk, montecarlo as mc, regimes as rgm, sites, storage as sto
from res_reliability import structural, study, validation
from res_reliability.config import HOURS, LIT_PV, LIT_WIND, RESULTS_DIR, YEARS

log = logging.getLogger(__name__)

SEED = 20251216
RES_SHARE = 0.10  # RES capacity per technology in the variant table, share of the peak load
N_UNITS = 100  # turbines / inverters per RES plant in the structural model
KWD = [  # capacity-market availability (derating) factors, onshore wind and solar, %
    {"auction": "main auction, delivery 2029", "source": "Dz.U. 2024 poz. 1154", "wind": 13.79, "pv": 2.80},
    {"auction": "main auction, delivery 2030", "source": "Dz.U. 2025 poz. 1066", "wind": 12.52, "pv": 2.41},
    {"auction": "additional auctions, delivery 2028", "source": "Dz.U. 2026 poz. 1359", "wind": 12.62, "pv": 1.56},
]


def save_json(name: str, obj, results_dir: Path = RESULTS_DIR) -> None:
    results_dir.mkdir(parents=True, exist_ok=True)
    (results_dir / name).write_text(json.dumps(obj, indent=1, default=float), encoding="utf-8")


def save_csv(name: str, df: pd.DataFrame, results_dir: Path = RESULTS_DIR) -> None:
    results_dir.mkdir(parents=True, exist_ok=True)
    df.to_csv(results_dir / name, index=False, float_format="%.6g")


def _r(x, n=4):
    return np.round(np.asarray(x, dtype=float), n).tolist()


# --- 0. Load profiles ---------------------------------------------------------------------

def load_stage(ctx: study.Context, data_dir: Path = config.DATA_DIR) -> dict:
    """Yearly KSE statistics, the scaling used, and a cross-check of PSE 2023 against ENTSO-E."""
    from res_reliability import pse, rts79

    s = pse.fetch(data_dir=data_dir)
    m = pse.year_matrix(s, YEARS)
    out = {"years": [{"year": y, "peak_mw": float(r.max()), "mean_mw": float(r.mean()),
                      "load_factor": float(r.mean() / r.max()), "energy_twh": float(r.sum() / 1e6)}
                     for y, r in zip(YEARS, m)],
           "kse_scaled_peak_mw": ctx.loads["KSE"].peak, "kse_note": ctx.loads["KSE"].note,
           "rts_load_factor": float(rts79.load_profile().mean() / rts79.PEAK_LOAD_MW),
           "api_start": str(pse.API_START.date())}
    sf = data_dir / "load" / "splice_factor.txt"
    if sf.exists():
        out["api_splice_factor"] = float(sf.read_text().strip())
    f = data_dir / "load" / "entsoe_actual_load_PL_2023.csv"
    if f.exists():
        e = pd.read_csv(f, parse_dates=["timestamp"]).set_index("timestamp")["demand_mw"]
        best = None
        for tz in ("UTC", "Europe/Warsaw"):
            idx = e.index.tz_localize(tz, ambiguous="NaT", nonexistent="NaT").tz_convert("UTC")
            ee = pd.Series(e.to_numpy(), index=idx)
            ee = ee[ee.index.notna()]
            ee = ee[~ee.index.duplicated()]
            j = pd.concat([s.rename("pse"), ee.rename("entsoe")], axis=1, join="inner").dropna()
            r = float(np.corrcoef(j["pse"], j["entsoe"])[0, 1])
            if best is None or r > best["r"]:
                best = {"assumed_timezone": tz, "hours": int(len(j)), "r": r,
                        "mean_pse_mw": float(j["pse"].mean()), "mean_entsoe_mw": float(j["entsoe"].mean()),
                        "mean_rel_diff": float((j["entsoe"].mean() - j["pse"].mean()) / j["pse"].mean())}
        out["entsoe_2023_check"] = best
    save_json("load_profiles.json", out)
    return out


# --- 1. Markov models and their validation ----------------------------------------------

def simulate_years(ctx: study.Context, s: sites.SiteSeries, tech: str, model: str, n: int, seed: int,
                   n_states: int = 5) -> np.ndarray:
    prod = study.production(ctx, s, tech, model, n_states)
    rngs = [mc.year_streams(seed, i, 2)[1] for i in range(n)]
    return prod.sample(rngs, np.arange(n), config.HOURS)


def markov_stage(ctx: study.Context) -> dict:
    pf, months = ctx.portfolio, ctx.months
    out: dict = {"portfolio": {}, "validation": {}}
    for tech in ("wind", "pv"):
        fit = study.fit_wind if tech == "wind" else study.fit_pv
        hom, sea = fit(pf, months, seasonal=False), fit(pf, months, seasonal=True)
        centre = fit(pf, months, seasonal=False, level_kind="centre")
        data = pf.wind if tech == "wind" else pf.pv
        k = pf.wind if tech == "wind" else pf.pv_k
        valid = None if tech == "wind" else pf.daylight
        states = mk.discretize(np.where(valid, k, np.nan) if valid is not None else k, 5)
        out["portfolio"][tech] = {
            "P_homogeneous": _r(hom.P), "pi_homogeneous": _r(hom.stationary()),
            "levels_mean": _r(hom.levels), "levels_centre": _r(centre.levels),
            "P_monthly": _r(sea.P), "levels_monthly": _r(sea.levels),
            "pi_monthly": [_r(sea.stationary(m)) for m in range(1, 13)],
            "intensities_per_h": _r(mk.intensities(states.reshape(-1) if valid is None else
                                                   states.reshape(-1)[valid.reshape(-1)], 5), 5),
            "mean_data": float(data.mean()),
        }
        # simulated years for chronology checks (one per weather year)
        sims = {m: simulate_years(ctx, pf, tech, m, len(YEARS), SEED + i)
                for i, m in enumerate(("homogeneous", "seasonal"))}
        v: dict = {}
        v["mean"] = {"data": float(data.mean()), **{m: float(x.mean()) for m, x in sims.items()}}
        v["monthly_mean"] = {"data": _r([data[:, months == mo].mean() for mo in range(1, 13)]),
                             **{m: _r([x[:, months == mo].mean() for mo in range(1, 13)]) for m, x in sims.items()}}
        hours = np.arange(config.HOURS) % 24
        v["diurnal_mean_utc"] = {"data": _r([data[:, hours == h].mean() for h in range(24)]),
                                 **{m: _r([x[:, hours == h].mean() for h in range(24)]) for m, x in sims.items()}}
        v["acf"] = {"data": _r(validation.acf(data, 72)), **{m: _r(validation.acf(x, 72)) for m, x in sims.items()}}
        if tech == "wind":
            s_data = mk.discretize(data.reshape(-1), 5)
            soj = {"data": validation.sojourns(s_data, 5)}
            for m, x in sims.items():
                soj[m] = validation.sojourns(mk.discretize(x.reshape(-1), 5), 5)
            v["sojourn_mean_h"] = {key: [float(a.mean()) if len(a) else None for a in val] for key, val in soj.items()}
            v["sojourn_geometric_h"] = _r(1.0 / (1.0 - np.diag(hom.P)), 2)
            # long low-output spells (state 0, CF below 12.5 %): exceedance counts per year
            spells = {key: val[0] for key, val in soj.items()}
            grid = [6, 12, 24, 36, 48, 72, 96, 120]
            v["low_spells_per_year_longer_than"] = {
                "hours": grid, **{key: [float((a >= g).sum() / len(YEARS)) for g in grid] for key, a in spells.items()}}
            v["longest_low_spell_h"] = {key: int(a.max()) for key, a in spells.items()}
        out["validation"][tech] = v

    # cross-validation: chronological split and leave-one-year-out, every grid point and the portfolio
    rows = []
    train, test = np.arange(len(YEARS) - 2), np.arange(len(YEARS) - 2, len(YEARS))
    for s in ctx.sites + [pf]:
        name = "portfolio" if s is pf else s.point.id
        for tech in ("wind", "pv"):
            x = s.wind if tech == "wind" else s.pv_k
            valid = None if tech == "wind" else s.daylight
            r = validation.cross_validate(x, months, train, test, 5, valid)
            rows.append({"site": name, "tech": tech, "split": f"{YEARS[0]}-{YEARS[-3]} / {YEARS[-2]}-{YEARS[-1]}",
                         **{k: val for k, val in r.items() if not isinstance(val, list)}})
    cv = pd.DataFrame(rows)
    save_csv("markov_cross_validation.csv", cv)
    loyo = []
    for tech in ("wind", "pv"):
        x = pf.wind if tech == "wind" else pf.pv_k
        valid = None if tech == "wind" else pf.daylight
        for i, y in enumerate(YEARS):
            tr = np.array([j for j in range(len(YEARS)) if j != i])
            r = validation.cross_validate(x, months, tr, np.array([i]), 5, valid)
            loyo.append({"tech": tech, "test_year": y, **{k: val for k, val in r.items() if not isinstance(val, list)}})
    save_csv("markov_leave_one_year_out.csv", pd.DataFrame(loyo))

    # per-point model summary for the map
    pts = []
    for s in ctx.sites:
        w_sea, p_sea = study.fit_wind(s, months), study.fit_pv(s, months)
        w_hom = study.fit_wind(s, months, seasonal=False)
        pts.append({"site": s.point.id, "lat": s.point.lat, "lon": s.point.lon,
                    "wind_cf": float(s.wind.mean()), "pv_cf": float(s.pv.mean()), "ws100_mean": float(s.ws.mean()),
                    "wind_cf_winter": float(s.wind[:, np.isin(months, [12, 1, 2])].mean()),
                    "wind_p00": float(w_hom.P[0, 0]), "wind_pi0": float(w_hom.stationary()[0]),
                    "wind_model_cf_seasonal": w_sea.mean(), "pv_k_mean_model": p_sea.mean()})
    save_csv("markov_points.csv", pd.DataFrame(pts))
    mats = {s.point.id: {"wind": _r(study.fit_wind(s, months, seasonal=False).P, 3),
                         "pv": _r(study.fit_pv(s, months, seasonal=False).P, 3),
                         "wind_levels": _r(study.fit_wind(s, months, seasonal=False).levels, 3),
                         "pv_levels": _r(study.fit_pv(s, months, seasonal=False).levels, 3)} for s in ctx.sites}
    save_json("markov_point_matrices.json", mats)
    save_json("markov_models.json", out)
    return out


# --- 2. Structural reliability ------------------------------------------------------------

def structural_stage(kelmarsh_dir: Path | None = None, pvdaq_dir: Path | None = None) -> dict:
    kel = structural.kelmarsh_reliability(kelmarsh_dir or config.kelmarsh_dir())
    pvdaq_dir = pvdaq_dir or Path("..") / "pv-anomaly-autoencoder" / "data" / "pvdaq" / "system_34"
    pv = structural.pvdaq_reliability(pvdaq_dir)
    w, p = kel["forced_ge_1h"], pv["ge_1h"]
    params = {
        "literature": {"wind": LIT_WIND, "pv": LIT_PV},
        "data": {"wind": {"for": w["for"], "mttr_h": w["mttr_h"], "source": kel["source"]},
                 "pv": {"for": p["for"], "mttr_h": p["mttr_h"], "source": pv["source"]}},
    }
    out = {"kelmarsh": {k: v for k, v in kel.items() if k != "durations_h"}, "pvdaq": pv, "parameters": params}
    save_json("structural.json", out)
    durations = pd.DataFrame({"duration_h": kel["durations_h"]})
    save_csv("kelmarsh_forced_outage_durations.csv", durations)
    return out


def structural_params(results_dir: Path = RESULTS_DIR) -> dict:
    return json.loads((results_dir / "structural.json").read_text(encoding="utf-8"))["parameters"]


# --- 3. System variants (sequential Monte Carlo) ------------------------------------------

def _plants(ctx, load, variant: dict, params: dict, techs=("wind", "pv")) -> list[mc.ResPlant]:
    cap = RES_SHARE * load.peak
    plants = []
    for tech in techs:
        prod = study.production(ctx, ctx.portfolio, tech, variant["model"], level_kind=variant.get("levels", "mean"))
        src = variant["structural"]
        if src == "none":
            st = None
        else:
            p = params[src][tech]
            mode = variant.get("structural_mode", "units")
            n = len(ctx.sites) if mode == "hourly_bernoulli" else N_UNITS  # thesis: FOR per farm
            st = mc.Structural(p["for"], p["mttr_h"], n, mode)
        plants.append(mc.ResPlant(tech, cap, prod, st))
    return plants


VARIANTS = [
    {"id": "BAU", "label": "conventional only", "techs": ()},
    {"id": "LIT", "label": "Markov monthly, literature FOR", "model": "seasonal", "structural": "literature"},
    {"id": "DATA", "label": "Markov monthly, FOR from open data", "model": "seasonal", "structural": "data"},
    {"id": "NOFOR", "label": "Markov monthly, no structural outages", "model": "seasonal", "structural": "none"},
    {"id": "THESIS", "label": "thesis implementation: homogeneous chain, bin-centre levels, hourly FOR per farm (literature)",
     "model": "homogeneous", "levels": "centre", "structural": "literature", "structural_mode": "hourly_bernoulli"},
    {"id": "HOM", "label": "Markov homogeneous, FOR from open data", "model": "homogeneous", "structural": "data"},
    {"id": "HIST", "label": "historical weather years, FOR from open data", "model": "hist", "structural": "data"},
    {"id": "DATA-W", "label": "wind only, Markov monthly, FOR from data", "model": "seasonal", "structural": "data",
     "techs": ("wind",)},
    {"id": "DATA-PV", "label": "PV only, Markov monthly, FOR from data", "model": "seasonal", "structural": "data",
     "techs": ("pv",)},
]


def _analytic(ctx, load, variant, params, techs):
    h = study.hours_of(load)
    cap = RES_SHARE * load.peak
    dist = None
    for tech in techs:
        if variant["model"] == "hist":
            d = elcc.historical(ctx.portfolio.wind if tech == "wind" else ctx.portfolio.pv, cap, h)
        else:
            fn = study.dist_wind if tech == "wind" else study.dist_pv
            d = fn(ctx, ctx.portfolio, cap, variant["model"], h, level_kind=variant.get("levels", "mean"))
        src = variant["structural"]
        if src != "none":
            # many independent units (turbines, or 29 farms in the thesis variant): ~deterministic derating
            d = elcc.apply_structural(d, params[src][tech]["for"], "units")
        dist = d if dist is None else dist + d
    return dist


def variants_stage(ctx: study.Context, n_sims: int = 2000, loads=("RTS", "KSE")) -> pd.DataFrame:
    params = structural_params()
    rows = []
    for lname in loads:
        load = ctx.loads[lname]
        base = mc.run(study.system_spec(ctx, load, []), n_sims, SEED, keep_margin=True)
        for v in VARIANTS:
            techs = v.get("techs", ("wind", "pv"))
            if v["id"] == "BAU":
                res, ana = base, elcc.indices(ctx.copt, load.matrix)
                e_mc = e_an = float("nan")
            else:
                plants = _plants(ctx, load, v, params, techs)
                res = mc.run(study.system_spec(ctx, load, plants), n_sims, SEED, keep_margin=True)
                dist = _analytic(ctx, load, v, params, techs)
                ana = elcc.indices(ctx.copt, load.matrix, dist)
                cap = RES_SHARE * load.peak * len(techs)
                e_mc = elcc.mc_elcc(base.margin, res.margin, cap) / cap
                e_an = elcc.elcc(ctx.copt, load.matrix, dist) / cap
            s = res.summary()
            rows.append({"load": lname, "variant": v["id"], "label": v["label"],
                         "res_mw_per_tech": 0.0 if v["id"] == "BAU" else RES_SHARE * load.peak,
                         "techs": "+".join(techs) or "-", **s,
                         "lole_analytic": ana["lole"], "eens_analytic": ana["eens"],
                         "elcc_mc_pct": 100 * e_mc, "elcc_analytic_pct": 100 * e_an})
            log.info("%s %-7s LOLE %.2f h (an. %.2f)  EENS %.0f  LOLF %.2f  ELCC %.1f %%", lname, v["id"],
                     s["lole"], ana["lole"], s["eens"], s["lolf"], 100 * e_mc if e_mc == e_mc else float("nan"))
    df = pd.DataFrame(rows)
    save_csv("variants.csv", df)
    # thesis-comparable: KSE profile scaled to 2850 MW, conventional only
    if "KSE2850" in ctx.loads:
        k = elcc.indices(ctx.copt, ctx.loads["KSE2850"].matrix)
        save_json("kse2850_bau.json", {"lole": k["lole"], "eens": k["eens"],
                                       "note": ctx.loads["KSE2850"].note})
    return df


# --- 4. ELCC: map, penetration curves, comparison with the capacity market ----------------

SHARES = [0.02, 0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0]


def elcc_stage(ctx: study.Context, loads=("KSE", "RTS")) -> dict:
    rows = []
    for lname in loads:
        load = ctx.loads[lname]
        h = study.hours_of(load)
        cap = RES_SHARE * load.peak
        for s in ctx.sites:
            for tech, fn in (("wind", study.dist_wind), ("pv", study.dist_pv)):
                for model in ("hist", "seasonal"):
                    d = fn(ctx, s, cap, model, h)
                    e = elcc.elcc(ctx.copt, load.matrix, d)
                    rows.append({"load": lname, "site": s.point.id, "lat": s.point.lat, "lon": s.point.lon,
                                 "tech": tech, "model": model, "capacity_mw": cap, "elcc_mw": e,
                                 "elcc_pct": 100 * e / cap, "cf_pct": 100 * d.mean_mw / cap})
        log.info("ELCC per point done for %s", lname)
    pts = pd.DataFrame(rows)
    save_csv("elcc_points.csv", pts)

    curves, lole_rows = [], []
    for lname in loads:
        load = ctx.loads[lname]
        h = study.hours_of(load)
        base = elcc.lole(ctx.copt, load.matrix)
        lole_rows.append({"load": lname, "mix": "none", "model": "-", "share": 0.0, "lole": base})
        for model in ("hist", "seasonal", "homogeneous"):
            for share in SHARES:
                cap = share * load.peak
                dw = study.dist_wind(ctx, ctx.portfolio, cap, model, h)
                dp = study.dist_pv(ctx, ctx.portfolio, cap, model, h)
                half_w = study.dist_wind(ctx, ctx.portfolio, cap / 2, model, h)
                half_p = study.dist_pv(ctx, ctx.portfolio, cap / 2, model, h)
                for mix, d in (("wind", dw), ("pv", dp), ("wind+pv", half_w + half_p)):
                    e = elcc.elcc(ctx.copt, load.matrix, d)
                    curves.append({"load": lname, "mix": mix, "model": model, "share": share,
                                   "capacity_mw": cap, "elcc_mw": e, "elcc_pct": 100 * e / cap,
                                   "cf_pct": 100 * d.mean_mw / cap})
                    lole_rows.append({"load": lname, "mix": mix, "model": model, "share": share,
                                      "lole": elcc.lole(ctx.copt, load.matrix, d)})
            log.info("ELCC curves %s %s done", lname, model)
    cdf = pd.DataFrame(curves)
    save_csv("elcc_curves.csv", cdf)
    save_csv("lole_vs_share.csv", pd.DataFrame(lole_rows))
    save_json("kwd.json", {"values": KWD, "note": "Korekcyjne współczynniki dyspozycyjności (capacity market "
                           "derating factors) set yearly by regulation of the minister for energy; verify the "
                           "current regulation before use."})
    return {"points": pts, "curves": cdf}


# --- headline numbers ------------------------------------------------------------------

def summary(results_dir: Path = RESULTS_DIR) -> dict:
    """Collect the headline numbers quoted in the README into results/summary.json."""
    j = lambda n: json.loads((results_dir / n).read_text(encoding="utf-8"))  # noqa: E731
    out: dict = {}
    v = pd.read_csv(results_dir / "variants.csv")
    out["variants"] = {f"{r.load}:{r.variant}": {"lole": round(r.lole, 3), "lole_ci95": round(r.lole_ci95, 3),
                                                  "lole_analytic": round(r.lole_analytic, 3), "eens": round(r.eens, 1),
                                                  "lolf": round(r.lolf, 3), "lold_h": round(r.lold_h, 2),
                                                  "elcc_mc_pct": round(r.elcc_mc_pct, 2),
                                                  "elcc_analytic_pct": round(r.elcc_analytic_pct, 2)}
                       for r in v.itertuples()}
    c = pd.read_csv(results_dir / "elcc_curves.csv")
    out["elcc_curves"] = {f"{r.load}:{r.mix}:{r.model}:{r.share:g}": round(r.elcc_pct, 2) for r in c.itertuples()}
    p = pd.read_csv(results_dir / "elcc_points.csv")
    out["elcc_points"] = {f"{k[0]}:{k[1]}:{k[2]}": {"min": round(g.elcc_pct.min(), 2), "median": round(g.elcc_pct.median(), 2),
                                                    "max": round(g.elcc_pct.max(), 2),
                                                    "cf_min": round(g.cf_pct.min(), 2), "cf_max": round(g.cf_pct.max(), 2),
                                                    "best_site": g.loc[g.elcc_pct.idxmax(), "site"],
                                                    "corr_elcc_cf": round(float(np.corrcoef(g.elcc_pct, g.cf_pct)[0, 1]), 3)}
                          for k, g in p.groupby(["load", "tech", "model"])}
    out["structural"] = j("structural.json")["parameters"]
    m = j("markov_models.json")["validation"]
    out["markov"] = {t: {"mean": m[t]["mean"], "acf24": {k: round(a[24], 3) for k, a in m[t]["acf"].items()}}
                     for t in ("wind", "pv")}
    out["markov"]["wind"]["spells"] = m["wind"]["low_spells_per_year_longer_than"]
    out["markov"]["wind"]["longest"] = m["wind"]["longest_low_spell_h"]
    cv = pd.read_csv(results_dir / "markov_cross_validation.csv")
    out["cross_validation"] = {t: {"ll_iid": round(g.ll_iid.median(), 4), "ll_homogeneous": round(g.ll_homogeneous.median(), 4),
                                   "ll_seasonal": round(g.ll_seasonal.median(), 4), "ll_test_fit": round(g.ll_test_fit.median(), 4),
                                   "ag_p_median": round(g.ag_p_value.median(), 4),
                                   "ag_p_below_005": int((g.ag_p_value < 0.05).sum()), "n": int(len(g)),
                                   "max_abs_dP_median": round(g.ag_max_abs_diff.median(), 4)}
                               for t, g in cv[cv.site != "portfolio"].groupby("tech")}
    for name in ("load_profiles.json", "kse2850_bau.json"):
        if (results_dir / name).exists():
            out[name.removesuffix(".json")] = j(name)
    s = pd.read_csv(results_dir / "sensitivity_states.csv")
    out["sensitivity_states"] = s[["n_states", "levels", "model", "elcc_pct", "lole_analytic"]].round(3).to_dict("records")
    ln = pd.read_csv(results_dir / "sensitivity_length.csv").groupby("years")["elcc_pct"]
    out["sensitivity_length"] = {int(k): {"mean": round(g.mean(), 2), "min": round(g.min(), 2), "max": round(g.max(), 2)}
                                 for k, g in ln}
    yr = pd.read_csv(results_dir / "sensitivity_years.csv")
    out["sensitivity_years"] = {t: {"min": round(g.elcc_pct.min(), 2), "max": round(g.elcc_pct.max(), 2),
                                    "min_year": int(g.loc[g.elcc_pct.idxmin(), "year"]),
                                    "max_year": int(g.loc[g.elcc_pct.idxmax(), "year"])} for t, g in yr.groupby("tech")}
    pr = pd.read_csv(results_dir / "sensitivity_weather_load_pairing.csv")
    out["weather_load_pairing"] = {t: {"matched": round(float(g[g.matched].elcc_pct.iloc[0]), 2),
                                       "shifted_mean": round(float(g[~g.matched].elcc_pct.mean()), 2),
                                       "shifted_min": round(float(g[~g.matched].elcc_pct.min()), 2),
                                       "shifted_max": round(float(g[~g.matched].elcc_pct.max()), 2)}
                                   for t, g in pr.groupby("tech")}
    ws = pd.read_csv(results_dir / "sensitivity_wind_speed.csv")
    out["sensitivity_wind_speed"] = ws.round(3).to_dict("records")
    fr = pd.read_csv(results_dir / "sensitivity_for.csv")
    out["sensitivity_for"] = fr[fr["for"].isin([0.0, 0.05, 0.1])].round(3).to_dict("records")
    if (results_dir / "weather_models.json").exists():
        w = j("weather_models.json")
        out["weather_models"] = {
            m: {"lulls_ge5d_0.15_per_year": v["lulls"]["0.15:all"]["at_least_per_year"][4],
                "winter_lulls_ge5d_0.15_per_year": v["lulls"]["0.15:winter"]["at_least_per_year"][4],
                "lulls_ge3d_0.1_per_year": v["lulls"]["0.1:all"]["at_least_per_year"][2],
                "wind_pv_lulls_ge3d_0.1_per_year": v["lulls"]["wind+pv 0.1:all"]["at_least_per_year"][2],
                "longest_lull_0.15_mean_d": round(v["longest_lull_0.15"]["mean"], 2),
                "spells_ge48h_per_year": round(v["hourly_low_spells_per_year"][">=48h"], 2),
                "acf24": v["acf"][24], "acf72": v["acf"][72], "acf168": v["acf"][168],
                "min_7d_wind": round(v["min_7d_wind"], 4), "min_7d_wind_pv": round(v["min_7d_wind_pv"], 4),
                "wind_top100_load_hours": round(v.get("wind_top100_load_hours", float("nan")), 4),
                "djf_below_0.15_by_tercile": [round(x, 3) for x in v["djf_by_temperature_tercile"]["share_days_below_0.15"]]}
            for m, v in w["metrics"].items()}
        out["elcc_no_storage"] = {f"{r['model']}:{r['mix']}": round(r["elcc_pct"], 2) for r in w.get("elcc_no_storage", [])}
        out["spatial"] = w["spatial"]
        out["regime_model"] = w["regime_model"]
    if (results_dir / "storage_elcc.csv").exists():
        s = pd.read_csv(results_dir / "storage_elcc.csv")
        out["storage_elcc_losses"] = {f"{r.model}:{r.res_share:g}:{r.duration_h:g}h:{r.power_pct_peak:g}%":
                                      [round(r.elcc_pct, 1), round(r.marginal_elcc_pct, 1)]
                                      for r in s[s["policy"] == "losses"].itertuples()}
        out["storage_base_lole"] = {f"{k[0]}:{k[1]:g}": round(g["base_lole"].iloc[0], 3)
                                    for k, g in s.groupby(["model", "res_share"])}
    if (results_dir / "storage_fleet_sweep.csv").exists():
        s = pd.read_csv(results_dir / "storage_fleet_sweep.csv")
        out["storage_fleet_sweep"] = {f"{r.model}:{r.duration_h:g}h:{r.power_pct_peak:g}%":
                                      [round(r.elcc_pct, 1), round(r.marginal_elcc_pct, 1)] for r in s.itertuples()}
    if (results_dir / "storage_lulls.csv").exists():
        s = pd.read_csv(results_dir / "storage_lulls.csv")
        out["storage_lulls_gwh_kse"] = {f"{r.model}:{r.energy_share:g}:{r.policy}:{r.duration_h:g}h":
                                        round(r.energy_gwh_kse, 1) for r in s.itertuples()}
        out["lole_without_storage"] = {f"{r.model}:{r.energy_share:g}": round(r.lole_no_storage, 2)
                                       for r in s.drop_duplicates(["model", "energy_share"]).itertuples()}
    save_json("summary.json", out, results_dir)
    return out


# --- 5. Sensitivity ---------------------------------------------------------------------

def sensitivity_stage(ctx: study.Context, lname: str = "KSE", n_sims: int = 500) -> dict:
    load = ctx.loads[lname]
    h = study.hours_of(load)
    pf, months = ctx.portfolio, ctx.months
    cap = RES_SHARE * load.peak
    out = {}

    # (a) number of states, level definition, homogeneous vs monthly; LOLF needs chronology -> MC
    rows = []
    base = mc.run(study.system_spec(ctx, load, []), n_sims, SEED, keep_margin=True)
    hist = elcc.historical(pf.wind, cap, h)
    ref = elcc.elcc(ctx.copt, load.matrix, hist)
    for n in (2, 3, 4, 5, 7, 10, 15):
        for levels in ("centre", "mean"):
            for model in ("homogeneous", "seasonal"):
                d = study.dist_wind(ctx, pf, cap, model, h, n_states=n, level_kind=levels)
                e = elcc.elcc(ctx.copt, load.matrix, d)
                row = {"n_states": n, "levels": levels, "model": model,
                       "lole_analytic": elcc.lole(ctx.copt, load.matrix, d), "elcc_pct": 100 * e / cap,
                       "elcc_hist_pct": 100 * ref / cap, "cf_pct": 100 * d.mean_mw / cap}
                if levels == "mean":
                    prod = study.production(ctx, pf, "wind", model, n_states=n, level_kind=levels)
                    r = mc.run(study.system_spec(ctx, load, [mc.ResPlant("wind", cap, prod)]), n_sims, SEED)
                    row.update({f"mc_{k}": v for k, v in r.summary().items() if k != "n_sims"})
                rows.append(row)
    r = mc.run(study.system_spec(ctx, load, [mc.ResPlant("wind", cap, mc.HistoricalProduction(pf.wind))]), n_sims, SEED)
    rows.append({"n_states": 0, "levels": "-", "model": "hist", "lole_analytic": elcc.lole(ctx.copt, load.matrix, hist),
                 "elcc_pct": 100 * ref / cap, "elcc_hist_pct": 100 * ref / cap, "cf_pct": 100 * hist.mean_mw / cap,
                 **{f"mc_{k}": v for k, v in r.summary().items() if k != "n_sims"}})
    out["states"] = pd.DataFrame(rows)
    save_csv("sensitivity_states.csv", out["states"])
    save_json("sensitivity_base.json", {"bau_mc": base.summary(), "load": lname, "wind_mw": cap})

    # (b) length of the series used to fit the (monthly) chain: every window of k consecutive years
    rows = []
    for k in (1, 2, 3, 5, 10):
        for start in range(0, len(YEARS) - k + 1):
            sel = np.arange(start, start + k)
            mdl = study.fit_wind(pf, months, rows=sel)
            d = elcc.from_markov(mdl, months, cap, None, len(YEARS), h)
            e = elcc.elcc(ctx.copt, load.matrix, d)
            rows.append({"years": k, "first_year": YEARS[start], "elcc_pct": 100 * e / cap,
                         "lole": elcc.lole(ctx.copt, load.matrix, d), "cf_model_pct": 100 * mdl.mean()})
    out["length"] = pd.DataFrame(rows)
    save_csv("sensitivity_length.csv", out["length"])

    # (c) climate: single weather years (historical, with the matching KSE year) and +/-5 % wind speed
    rows = []
    for i, y in enumerate(YEARS):
        lm = load.matrix[i % load.matrix.shape[0]][None]
        for tech, series in (("wind", pf.wind), ("pv", pf.pv)):
            d = elcc.historical(series[i][None], cap, h)
            e = elcc.elcc(ctx.copt, lm, d)
            rows.append({"year": y, "tech": tech, "elcc_pct": 100 * e / cap, "cf_pct": 100 * series[i].mean(),
                         "lole_bau": elcc.lole(ctx.copt, lm), "lole_with_res": elcc.lole(ctx.copt, lm, d)})
    out["years"] = pd.DataFrame(rows)
    save_csv("sensitivity_years.csv", out["years"])
    # weather-demand correlation: pair every load year with a different weather year
    rows = []
    for tech, series in (("wind", pf.wind), ("pv", pf.pv)):
        for k in range(len(YEARS)):
            d = elcc.historical(np.roll(series, -k, axis=0), cap, h)
            rows.append({"tech": tech, "year_shift": k, "matched": k == 0,
                         "elcc_pct": 100 * elcc.elcc(ctx.copt, load.matrix, d) / cap})
    out["decorrelated"] = pd.DataFrame(rows)
    save_csv("sensitivity_weather_load_pairing.csv", out["decorrelated"])
    rows = []
    from res_reliability import power
    for f in (0.90, 0.95, 1.0, 1.05, 1.10):
        w = np.mean([power.wind_cf(s.ws * f) * (1 - config.WIND_LOSSES) for s in ctx.sites], axis=0)
        sc = sites.SiteSeries(pf.point, w, pf.pv, pf.pv_cs, pf.ws * f)
        for model in ("hist", "seasonal"):
            d = study.dist_wind(ctx, sc, cap, model, h)
            e = elcc.elcc(ctx.copt, load.matrix, d)
            rows.append({"wind_speed_factor": f, "model": model, "cf_pct": 100 * d.mean_mw / cap,
                         "elcc_pct": 100 * e / cap, "lole": elcc.lole(ctx.copt, load.matrix, d)})
    out["climate"] = pd.DataFrame(rows)
    save_csv("sensitivity_wind_speed.csv", out["climate"])

    # (d) forced outage rate of wind: many independent turbines vs whole-plant hourly outages (thesis)
    rows = []
    d0 = study.dist_wind(ctx, pf, cap, "seasonal", h)
    for q in np.round(np.arange(0.0, 0.101, 0.01), 2):
        for mode in ("units", "plant_bernoulli"):
            d = elcc.apply_structural(d0, float(q), mode)
            e = elcc.elcc(ctx.copt, load.matrix, d)
            rows.append({"for": q, "mode": mode, "lole": elcc.lole(ctx.copt, load.matrix, d),
                         "eens": elcc.indices(ctx.copt, load.matrix, d)["eens"], "elcc_pct": 100 * e / cap})
    out["for"] = pd.DataFrame(rows)
    save_csv("sensitivity_for.csv", out["for"])
    return out


# --- 6. Weather models with memory: regimes, temperature, multi-day lulls ------------------

WEATHER_MODELS = ("seasonal", "regime", "regime_t")
LULL_THRESHOLDS = (0.10, 0.15)  # daily mean wind capacity factor of the national portfolio
LULL_MAX_DAYS = 10
N_WEATHER_PER_YEAR = 20  # synthetic years per weather year in the validation
WINTER_MONTHS = (10, 11, 12, 1, 2, 3)


def simulate_weather(ctx: study.Context, model: str, n_per_year: int = N_WEATHER_PER_YEAR, seed: int = SEED) -> dict:
    """Synthetic wind and PV years of the national portfolio with the random streams of ``mc.run``.

    Synthetic year i belongs to weather (and load) year i % n_years. Returns ``wind``, ``pv`` (n, 8760),
    ``years`` and, for regime models, the daily ``regimes`` (n, 365).
    """
    n = n_per_year * len(YEARS)
    sims = np.arange(n)
    years = sims % len(YEARS)
    streams = [mc.year_streams(seed, int(i), 5) for i in sims]
    regimes = None
    if model in study.REGIME_MODELS:
        w = study.regime_weather(ctx, model)
        regimes = w.regime_paths([mc.weather_stream(seed, int(i)) for i in sims], years)
    out = {"years": years}
    for j, tech in enumerate(("wind", "pv")):
        prod = study.production(ctx, ctx.portfolio, tech, model)
        rngs = [s[1 + 2 * j] for s in streams]
        out[tech] = prod.sample(rngs, years, HOURS, regimes) if regimes is not None else prod.sample(rngs, years, HOURS)
    if regimes is not None:
        out["regimes"] = regimes[:, study.regime_weather(ctx, model).burn_in_days:]
    return out


def _longest(daily: np.ndarray, threshold: float) -> np.ndarray:
    return np.array([max(validation.low_runs(row[None], threshold), default=0) for row in daily])


def weather_metrics(ctx: study.Context, wind: np.ndarray, pv: np.ndarray, years: np.ndarray) -> dict:
    """Chronology, lulls and demand coupling of (n, 8760) wind and PV series; ``years`` = load year of each."""
    n = len(wind)
    md = rgm.day_months(ctx.months)
    wd, pdl = rgm.daily_mean(wind), rgm.daily_mean(pv)
    comb = 0.5 * (wd + pdl)  # equal wind and PV capacity
    winter = np.isin(md, WINTER_MONTHS)
    out: dict = {"n_years": n, "lulls": {}}
    for th in LULL_THRESHOLDS:
        for season, x in (("all", wd), ("winter", np.where(winter, wd, 1.0))):
            runs = validation.low_runs(x, th)
            counts = validation.run_length_counts(runs, n, LULL_MAX_DAYS)
            out["lulls"][f"{th:g}:{season}"] = {"per_year": _r(counts),
                                                "at_least_per_year": _r(counts[::-1].cumsum()[::-1])}
        long_ = _longest(wd, th)
        out[f"longest_lull_{th:g}"] = {"mean": float(long_.mean()), "p90": float(np.quantile(long_, 0.9)),
                                       "max": int(long_.max())}
    runs = validation.low_runs(comb, 0.10)
    counts = validation.run_length_counts(runs, n, LULL_MAX_DAYS)
    out["lulls"]["wind+pv 0.1:all"] = {"per_year": _r(counts), "at_least_per_year": _r(counts[::-1].cumsum()[::-1])}
    soj = validation.sojourns(mk.discretize(wind.reshape(-1), 5), 5)[0]
    out["hourly_low_spells_per_year"] = {f">={h}h": float((soj >= h).sum() / n) for h in (24, 48, 96)}
    out["acf"] = _r(validation.acf(wind, 168), 3)
    out["acf_pv_daily"] = _r(validation.acf(pdl, 10), 3)
    out["monthly_mean"] = {"wind": _r([wind[:, ctx.months == m].mean() for m in range(1, 13)]),
                           "pv": _r([pv[:, ctx.months == m].mean() for m in range(1, 13)])}
    for k in (3, 7):
        out[f"min_{k}d_wind"] = float(validation.rolling_min_mean(wd, k).mean())
        out[f"min_{k}d_wind_pv"] = float(validation.rolling_min_mean(comb, k).mean())
    if "KSE" in ctx.loads:
        load = ctx.loads["KSE"].matrix
        top = np.argsort(load, axis=1)[:, -100:]
        out["wind_top100_load_hours"] = float(np.mean([wind[i, top[y]].mean() for i, y in enumerate(years)]))
        out["wind_annual_mean"] = float(wind.mean())
        ld = rgm.daily_mean(load)[years]
        nov_feb = np.isin(md, (11, 12, 1, 2))
        out["r_daily_wind_load_nov_feb"] = float(np.corrcoef(wd[:, nov_feb].ravel(), ld[:, nov_feb].ravel())[0, 1])
    w_t = study.regime_weather(ctx, "regime_t")
    cov = w_t.cov[years % w_t.cov.shape[0]]
    djf = np.broadcast_to(rgm.SEASON_OF_MONTH[md - 1] == 0, wd.shape)
    out["djf_by_temperature_tercile"] = {
        "mean_wind": [float(wd[djf & (cov == c)].mean()) for c in range(w_t.process.n_cov)],
        "share_days_below_0.15": [float((wd[djf & (cov == c)] < 0.15).mean()) for c in range(w_t.process.n_cov)]}
    return out


def _sampled_dist(x: np.ndarray, capacity_mw: float, n_years: int) -> elcc.ResDistribution:
    """Equally likely simulated years per weather year as an hourly RES distribution for the COPT."""
    n = len(x) // n_years
    lv = capacity_mw * x[: n * n_years].reshape(n, n_years, -1).transpose(1, 2, 0)
    return elcc.ResDistribution(lv, np.full(lv.shape, 1.0 / n), capacity_mw)


def spatial_check(ctx: study.Context, n_per_year: int = 5, seed: int = SEED) -> dict:
    """Does a regime path shared by all grid points carry their spatial correlation?

    Every point gets hourly wind chains conditioned on the national regime of the day; the points then move
    independently given one shared regime path. Compared with independent monthly chains per point and with
    ERA5: mean correlation of deseasonalised daily capacity factors between points, and the day-to-day spread
    and lulls of the equal-weight portfolio rebuilt from the points.
    """
    w = study.regime_weather(ctx, "regime_t")
    wd = rgm.daily_mean(ctx.portfolio.wind)
    reg = np.digitize(wd / rgm.climatology(wd), w.edges)
    rg_h = np.repeat(reg, 24, axis=1).reshape(-1)
    m_h = np.tile(ctx.months, len(YEARS))
    n = n_per_year * len(YEARS)
    sims = np.arange(n)
    years = sims % len(YEARS)
    regs = w.regime_paths([mc.weather_stream(seed, int(i)) for i in sims], years)
    b = w.burn_in_days * 24
    series = {"data": [], "regime": [], "independent": []}
    for j, s in enumerate(ctx.sites):
        ch = rgm.fit_hourly(s.wind.reshape(-1), m_h, rg_h, w.n_regimes)
        pw = rgm.RegimeWeather(w.process, ch, w.pv, w.months, w.edges, w.cov, w.burn_in_days)
        x = pw.hourly("wind", [np.random.default_rng([seed, int(i), j]) for i in sims], regs)[:, b:b + HOURS]
        ind = study.production(ctx, s, "wind", "seasonal").sample(
            [np.random.default_rng([seed, int(i), 1000 + j]) for i in sims], years, HOURS)
        series["data"].append(rgm.daily_mean(s.wind))
        series["regime"].append(rgm.daily_mean(x))
        series["independent"].append(rgm.daily_mean(ind))
    out = {}
    for key, daily in series.items():
        a = np.stack([d - rgm.climatology(d) for d in daily])  # (points, years, days)
        c = np.corrcoef(a.reshape(len(daily), -1))
        iu = np.triu_indices(len(daily), 1)
        port = np.mean(daily, axis=0)
        anom = port - rgm.climatology(port)
        runs = validation.low_runs(port, 0.15)
        out[key] = {"mean_pair_corr": float(c[iu].mean()), "portfolio_daily_std": float(anom.std()),
                    "portfolio_lulls_ge5d_per_year": float((runs >= 5).sum() / len(port))}
    return out


def weather_stage(ctx: study.Context, n_per_year: int = N_WEATHER_PER_YEAR) -> dict:
    pf = ctx.portfolio
    w0, w_t = study.regime_weather(ctx, "regime"), study.regime_weather(ctx, "regime_t")
    wd = rgm.daily_mean(pf.wind)
    reg = np.digitize(wd / rgm.climatology(wd), w_t.edges)
    md = rgm.day_months(ctx.months)
    djf = np.broadcast_to(rgm.SEASON_OF_MONTH[md - 1] == 0, reg.shape)
    run_len = [validation.low_runs((reg != k).astype(float), 0.5) for k in range(w_t.n_regimes)]
    out: dict = {
        "regime_model": {
            "regime_quantiles": [0.2, 0.5, 0.8], "edges_ratio_to_climatology": _r(w_t.edges),
            "regimes": ["calm", "below normal", "above normal", "windy"],
            "covariate": "tercile of the daily mean temperature anomaly within the season (DJF, MAM, JJA, SON)",
            "max_age_days": w_t.process.max_age, "kappa": 10.0, "burn_in_days": w_t.burn_in_days,
            "mean_run_days_data": _r([r.mean() for r in run_len], 2),
            "hazard_calm_djf_by_tercile": _r(w_t.process.hazard[0, 0], 3),
            "hazard_calm_djf_pooled": _r(w0.process.hazard[0, 0, 0], 3),
            "hazard_windy_djf_by_tercile": _r(w_t.process.hazard[-1, 0], 3),
            "regime_freq_djf_by_tercile_data": [
                _r(np.bincount(reg[djf & (w_t.cov == c)], minlength=w_t.n_regimes) / (djf & (w_t.cov == c)).sum(), 3)
                for c in range(w_t.process.n_cov)],
        },
        "metrics": {"data": weather_metrics(ctx, pf.wind, pf.pv, np.arange(len(YEARS)))},
    }
    sims = {}
    for m in WEATHER_MODELS:
        sims[m] = simulate_weather(ctx, m, n_per_year)
        out["metrics"][m] = weather_metrics(ctx, sims[m]["wind"], sims[m]["pv"], sims[m]["years"])
        log.info("weather %s: lulls >= 5 d below 0.15: %.2f /yr", m,
                 out["metrics"][m]["lulls"]["0.15:all"]["at_least_per_year"][4])
    rs = sims["regime_t"]["regimes"]
    cov = w_t.cov[sims["regime_t"]["years"]]
    djf_s = np.broadcast_to(rgm.SEASON_OF_MONTH[md - 1] == 0, rs.shape)
    out["regime_model"]["regime_freq_djf_by_tercile_regime_t"] = [
        _r(np.bincount(rs[djf_s & (cov == c)], minlength=w_t.n_regimes) / (djf_s & (cov == c)).sum(), 3)
        for c in range(w_t.process.n_cov)]

    rows = []
    for key, mt in out["metrics"].items():
        for name, v in mt["lulls"].items():
            th, season = name.rsplit(":", 1)
            for k, val in enumerate(v["per_year"]):
                rows.append({"model": key, "series": th, "season": season, "days": k + 1, "per_year": val,
                             "at_least_per_year": v["at_least_per_year"][k]})
    save_csv("weather_lulls.csv", pd.DataFrame(rows))

    # ELCC without storage (KSE load, 10 % of peak per technology): only the hour-by-hour distribution of RES
    # output and its alignment with the load matter, so simulated years enter the exact COPT expectation
    if "KSE" in ctx.loads:
        load = ctx.loads["KSE"]
        h, cap = study.hours_of(load), RES_SHARE * load.peak
        erows = []
        for m in ("hist",) + WEATHER_MODELS:
            if m == "hist":
                dw, dp = elcc.historical(pf.wind, cap, h), elcc.historical(pf.pv, cap, h)
                dwp = elcc.historical(pf.wind + pf.pv, 2 * cap, h)
            else:
                x = sims[m]
                dw, dp = _sampled_dist(x["wind"], cap, len(YEARS)), _sampled_dist(x["pv"], cap, len(YEARS))
                dwp = _sampled_dist(x["wind"] + x["pv"], 2 * cap, len(YEARS))
            for mix, d in (("wind", dw), ("pv", dp), ("wind+pv", dwp)):
                e = elcc.elcc(ctx.copt, load.matrix, d)
                erows.append({"model": m, "mix": mix, "capacity_mw": d.capacity_mw,
                              "elcc_pct": 100 * e / d.capacity_mw, "lole": elcc.lole(ctx.copt, load.matrix, d)})
            log.info("ELCC without storage, %s: wind %.1f %%", m, erows[-3]["elcc_pct"])
        out["elcc_no_storage"] = erows
    out["spatial"] = spatial_check(ctx)
    save_json("weather_models.json", out)
    return out


# --- 7. Storage: ELCC by duration, fleet size and RES share; storage needed for lulls ---------

STORAGE_MODELS = ("seasonal", "regime", "regime_t", "hist")
STORAGE_DURATIONS = (1, 2, 4, 8)  # h
STORAGE_FLEET = (0.01, 0.05, 0.10, 0.20)  # storage power, share of the peak load
STORAGE_RES = (0.0, 0.10, 0.25, 0.50)  # wind and PV capacity each, share of the peak load
POLICIES = {
    "ideal": {"efficiency": 1.0, "soc0": 1.0, "charge_from": "any"},
    "losses": {"efficiency": 0.85, "soc0": 0.5, "charge_from": "any"},
    "colocated": {"efficiency": 0.85, "soc0": 0.5, "charge_from": "res"},
}
LULL_SHARES = (0.0, 0.1, 0.2, 0.3, 0.4)  # wind + PV share of annual demand
LULL_DURATIONS = (4.0, 8.0, 24.0, float("inf"))  # inf: power not limited
LULL_WARMUP_H = 720
KWD_STORAGE = [  # capacity-market derating factor of batteries (one value for every duration), %
    {"auction": "main auction, delivery 2028", "source": "Dz.U. 2023 poz. 1561, par. 10 pkt 9", "kwd": 95.00},
    {"auction": "main auction, delivery 2029", "source": "Dz.U. 2024 poz. 1154, par. 10 pkt 9", "kwd": 61.30},
    {"auction": "main auction, delivery 2030",
     "source": "Dz.U. 2025 poz. 1066, par. 10 pkt 9 (one group with pumped storage)", "kwd": 13.39},
]

_WORKER: dict = {}


def _ctx(data_dir) -> study.Context:
    if "ctx" not in _WORKER:
        _WORKER["ctx"] = study.build(data_dir)
    return _WORKER["ctx"]


def _map(fn, tasks: list[tuple], jobs: int) -> list:
    if jobs <= 1:
        return [fn(*t) for t in tasks]
    with ProcessPoolExecutor(jobs) as ex:
        return list(ex.map(fn, *zip(*tasks)))


def _res_plants(ctx: study.Context, model: str, wind_mw: float, pv_mw: float, params: dict) -> list[mc.ResPlant]:
    plants = []
    for tech, cap in (("wind", wind_mw), ("pv", pv_mw)):
        if cap > 0:
            p = params["data"][tech]
            plants.append(mc.ResPlant(tech, cap, study.production(ctx, ctx.portfolio, tech, model),
                                      mc.Structural(p["for"], p["mttr_h"], N_UNITS)))
    return plants


def _load_factor(ctx: study.Context, share: float) -> float:
    """Factor on the KSE load that brings LOLE with ``share`` of wind and PV (historical weather) back to 3 h."""
    load = ctx.loads["KSE"]
    if share <= 0:
        return 1.0
    cap, h = share * load.peak, study.hours_of(load)
    d = elcc.historical(ctx.portfolio.wind, cap, h) + elcc.historical(ctx.portfolio.pv, cap, h)
    f = lambda k: elcc.lole(ctx.copt, k * load.matrix, d) - config.PL_RELIABILITY_STANDARD_LOLE  # noqa: E731
    return float(brentq(f, 1.0, 2.0, xtol=1e-4))


def _storage_elcc_task(model: str, share: float, n_sims: int, data_dir) -> list[dict]:
    ctx = _ctx(data_dir)
    load = ctx.loads["KSE"]
    params = structural_params()
    k = _load_factor(ctx, share)
    plants = _res_plants(ctx, model, share * load.peak, share * load.peak, params)
    lc = study.LoadCase("KSE", k * load.matrix, "")
    r = mc.run(study.system_spec(ctx, lc, plants), n_sims, SEED, keep_margin=True, keep_res=share > 0)
    base = r.summary()
    P = np.repeat(STORAGE_FLEET, len(STORAGE_DURATIONS)) * load.peak
    D = np.tile(STORAGE_DURATIONS, len(STORAGE_FLEET)).astype(float)
    rows = []
    for pol, pp in POLICIES.items():
        if pp["charge_from"] == "res" and share <= 0:
            continue
        cap = r.res if pp["charge_from"] == "res" else None
        e = np.empty(len(P))
        for p in np.unique(P):  # one power level at a time: fewer relevant years, more skipped hours
            sel = P == p
            e[sel] = sto.elcc(r.margin, P[sel], P[sel] * D[sel], pp["efficiency"], pp["soc0"], cap)
        for pw, dur, val in zip(P, D, e):
            rows.append({"model": model if share > 0 else "none", "res_share": share, "load_factor": k,
                         "policy": pol, "duration_h": dur, "power_pct_peak": 100 * pw / load.peak, "power_mw": pw,
                         "elcc_mw": val, "elcc_pct": 100 * val / pw, "base_lole": base["lole"],
                         "base_lole_ci95": base["lole_ci95"], "base_lold_h": base["lold_h"]})
    log.info("storage ELCC %s, RES %.0f %% x2: done (base LOLE %.2f h)", model, 100 * share, base["lole"])
    return rows


SWEEP_FLEET = (0.01, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30, 0.40, 0.50)
SWEEP_DURATIONS = (2, 4, 8)
SWEEP_RES = 0.25


def _fleet_sweep_task(model: str, n_sims: int, data_dir) -> list[dict]:
    """ELCC of ever larger storage fleets (2 x 25 % wind and PV, round trip 85 %), for the comparison with KWD."""
    ctx = _ctx(data_dir)
    load = ctx.loads["KSE"]
    k = _load_factor(ctx, SWEEP_RES)
    plants = _res_plants(ctx, model, SWEEP_RES * load.peak, SWEEP_RES * load.peak, structural_params())
    lc = study.LoadCase("KSE", k * load.matrix, "")
    r = mc.run(study.system_spec(ctx, lc, plants), n_sims, SEED, keep_margin=True)
    pp = POLICIES["losses"]
    rows = []
    for f in SWEEP_FLEET:
        p = f * load.peak
        d = np.array(SWEEP_DURATIONS, dtype=float)
        e = sto.elcc(r.margin, np.full(len(d), p), p * d, pp["efficiency"], pp["soc0"])
        rows += [{"model": model, "res_share": SWEEP_RES, "policy": "losses", "duration_h": dur,
                  "power_pct_peak": 100 * f, "power_mw": p, "elcc_mw": val, "elcc_pct": 100 * val / p}
                 for dur, val in zip(d, e)]
    log.info("fleet sweep %s done", model)
    return rows


def fleet_sweep(ctx: study.Context, n_sims: int = 2000, jobs: int = 1, data_dir=config.DATA_DIR) -> pd.DataFrame:
    _WORKER["ctx"] = ctx
    rows = [r for rr in _map(_fleet_sweep_task, [(m, n_sims, data_dir) for m in ("hist", "regime_t", "seasonal")],
                             jobs) for r in rr]
    df = pd.DataFrame(rows).sort_values(["model", "duration_h", "power_mw"])
    prev = df.groupby(["model", "duration_h"])[["elcc_mw", "power_mw"]].shift(1).fillna(0.0)
    df["marginal_elcc_pct"] = 100 * (df["elcc_mw"] - prev["elcc_mw"]) / (df["power_mw"] - prev["power_mw"])
    lp = json.loads((RESULTS_DIR / "load_profiles.json").read_text(encoding="utf-8"))
    df["power_gw_kse"] = df["power_pct_peak"] / 100 * float(lp["years"][-1]["peak_mw"]) / 1000
    save_csv("storage_fleet_sweep.csv", df)
    return df


def _lull_task(model: str, share: float, target: float, n_sims: int, data_dir) -> list[dict]:
    """Storage energy that keeps LOLE at ``target`` when demand grows and wind + PV cover the growth."""
    ctx = _ctx(data_dir)
    load = ctx.loads["KSE"]
    params = structural_params()
    growth = share / (1 - share)
    cap = growth * load.matrix.mean() / (ctx.portfolio.wind.mean() + ctx.portfolio.pv.mean())
    plants = _res_plants(ctx, model, cap, cap, params)
    lc = study.LoadCase("KSE", (1 + growth) * load.matrix, "")
    r = mc.run(study.system_spec(ctx, lc, plants), n_sims, SEED, keep_margin=True)
    s = r.summary()
    durs = np.array(LULL_DURATIONS * 2)
    pols = ["ideal"] * len(LULL_DURATIONS) + ["losses"] * len(LULL_DURATIONS)
    eff = np.array([POLICIES[p]["efficiency"] for p in pols])
    e = sto.energy_needed(r.margin, target, durs, eff, 1.0, LULL_WARMUP_H, e_min=1.0, e_max=5e5, iters=13)
    rows = []
    for dur, pol, val in zip(durs, pols, e):
        rows.append({"model": model if share > 0 else "none", "energy_share": share, "growth": growth,
                     "wind_mw": cap, "pv_mw": cap, "load_peak_mw": (1 + growth) * load.peak,
                     "lole_no_storage": s["lole"], "eens_no_storage": s["eens"], "lole_target": target,
                     "duration_h": dur, "policy": pol, "energy_mwh": val,
                     "power_mw": val / dur if np.isfinite(dur) else np.nan})
    log.info("lull storage %s, share %.0f %%: LOLE %.1f h without storage", model, 100 * share, s["lole"])
    return rows


def storage_stage(ctx: study.Context, n_sims: int = 2000, jobs: int = 1, data_dir=config.DATA_DIR) -> dict:
    if "KSE" not in ctx.loads:
        raise RuntimeError("the storage stage needs the KSE load")
    _WORKER["ctx"] = ctx
    tasks = [("-", 0.0, n_sims, data_dir)] + [(m, s, n_sims, data_dir) for s in STORAGE_RES if s > 0
                                               for m in STORAGE_MODELS]
    rows = [r for rr in _map(_storage_elcc_task, tasks, jobs) for r in rr]
    df = pd.DataFrame(rows).sort_values(["model", "res_share", "policy", "duration_h", "power_mw"])
    # incremental ELCC of each block of storage over the fleet below it
    prev = df.groupby(["model", "res_share", "policy", "duration_h"])[["elcc_mw", "power_mw"]].shift(1).fillna(0.0)
    df["marginal_elcc_pct"] = 100 * (df["elcc_mw"] - prev["elcc_mw"]) / (df["power_mw"] - prev["power_mw"])
    save_csv("storage_elcc.csv", df)
    fleet_sweep(ctx, n_sims, jobs, data_dir)

    base = mc.run(study.system_spec(ctx, ctx.loads["KSE"], []), n_sims, SEED)
    target = float(base.lole.mean())
    tasks = [("-", 0.0, target, n_sims, data_dir)] + [(m, s, target, n_sims, data_dir) for s in LULL_SHARES
                                                       if s > 0 for m in STORAGE_MODELS]
    ldf = pd.DataFrame([r for rr in _map(_lull_task, tasks, jobs) for r in rr])
    lp = json.loads((RESULTS_DIR / "load_profiles.json").read_text(encoding="utf-8"))
    kse_peak = float(lp["years"][-1]["peak_mw"])
    ldf["energy_h_of_peak"] = ldf["energy_mwh"] / ctx.loads["KSE"].peak
    ldf["energy_gwh_kse"] = ldf["energy_h_of_peak"] * kse_peak / 1000
    save_csv("storage_lulls.csv", ldf)
    save_json("kwd_storage.json", {"values": KWD_STORAGE, "kse_peak_mw": kse_peak, "kse_peak_year": YEARS[-1],
                                   "model_peak_mw": ctx.loads["KSE"].peak,
                                   "note": "Capacity-market derating factors (KWD) of batteries, flywheels and "
                                   "supercapacitors in the main auctions; one value for every storage duration. "
                                   "Set yearly by regulation; verify the current regulation before use."})
    return {"elcc": df, "lulls": ldf}
