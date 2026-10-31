"""Stages on operating data: theory against SCADA (``theory``) and weather-dependent failures (``failures``).

Usage:
    python -m res_reliability fetch-farms
    python -m res_reliability run --stages theory failures --n-sims 2000
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

from res_reliability import config, elcc, farms, integration as ig, markov as mk, montecarlo as mc, power, study
from res_reliability import theory as th
from res_reliability import sites as site_series
from res_reliability.analyses import RES_SHARE, SEED, save_csv, save_json, structural_params
from res_reliability.config import HOURS, LIT_PV, LIT_WIND, WIND_LOSSES

log = logging.getLogger(__name__)

MIN_OUTAGE_H = 1.0  # forced outages of an hour or more, as for the structural parameters of the main study
WS_EDGES = (3.0, 6.0, 9.0, 12.0, 15.0)  # hub-height wind-speed bins of the farm-level integrated model
WS_LABELS = ("<3", "3-6", "6-9", "9-12", "12-15", ">=15")
N_UNITS = 100


def cached_outages(farm: farms.Farm, data_dir: Path = config.DATA_DIR) -> tuple[pd.DataFrame, pd.DataFrame]:
    d = farms.farms_dir(data_dir)
    fo, win = d / f"fo_{farm.key}.parquet", d / f"win_{farm.key}.parquet"
    if fo.exists() and win.exists():
        return pd.read_parquet(fo), pd.read_parquet(win)
    out, w = farms.forced_outages(farm)
    d.mkdir(parents=True, exist_ok=True)
    out.to_parquet(fo)
    w.to_parquet(win)
    return out, w


def site_reliability(out: pd.DataFrame, win: pd.DataFrame) -> dict:
    e = out[out["duration_h"] >= MIN_OUTAGE_H]
    exposure = float(win["hours"].sum())
    return {"for": float(e["duration_h"].sum() / exposure), "mttr_h": float(e["duration_h"].mean()),
            "outages": int(len(e)), "turbine_years": exposure / 8760.0,
            "outages_per_turbine_year": float(len(e) / (exposure / 8760.0))}


def _adequacy(ctx: study.Context, lname: str, model: mk.MarkovModel, for_: float, envelope=None) -> dict:
    load = ctx.loads[lname]
    h = study.hours_of(load)
    cap = RES_SHARE * load.peak
    d = elcc.from_markov(model, ctx.months, cap, envelope, 1, h)
    d = elcc.apply_structural(d, for_, "units")
    ind = elcc.indices(ctx.copt, load.matrix, d)
    return {"load": lname, "cf": d.mean_mw / cap, "lole": ind["lole"], "eens": ind["eens"],
            "elcc_pct": 100 * elcc.elcc(ctx.copt, load.matrix, d) / cap}


# --- A. theory against operating data -------------------------------------------------------

def theory_stage(ctx: study.Context, n_sims: int = 2000, data_dir: Path = config.DATA_DIR) -> dict:
    sites, bins, structure, adequacy, tests = [], [], [], [], {}
    maps = {}
    lit = LIT_WIND["for"]
    for key, f in farms.FARMS.items():
        out, win = cached_outages(f, data_dir)
        rel = site_reliability(out, win)
        df = th.wind_site(f, out, data_dir)
        v = df["valid"].to_numpy()
        for col in ("theory", "model", "operation_all", "operation"):
            m = th.series_metrics(df, "theory", col)
            sites.append({"site": key, "tech": "wind", "series": col, **{k: m[k] for k in m if k != "ratio_by_year"},
                          "ratio_year_min": min(m["ratio_by_year"].values()),
                          "ratio_year_max": max(m["ratio_by_year"].values()), "n_years": len(m["ratio_by_year"]),
                          "for_site": rel["for"], "mttr_site_h": rel["mttr_h"]})
        bins.append(th.wind_bin_table(df).assign(site=key))
        pr = float(df.loc[v, "operation"].sum() / df.loc[v, "theory"].sum())
        series = {"theory": df["theory"], "theory_scaled": df["theory"] * pr, "model": df["model"],
                  "operation": df["operation"]}
        months = df["month"].to_numpy()
        fitted, st = {}, {}
        for name, x in series.items():
            s = th.markov_structure(x.to_numpy(), v, months)
            fitted[name], st[name] = s["model"], s
            for k in range(th.N_STATES):
                structure.append({"site": key, "series": name, "state": k, "occupancy": s["occupancy"][k],
                                  "residence_h": s["residence_h"][k], "level": s["levels"][k],
                                  "acf1": s["acf"][1], "acf6": s["acf"][6], "acf24": s["acf"][24]})
        tests[key] = {"operation_vs_theory": th.transition_test(st["operation"]["states"],
                                                                st["theory"]["homogeneous"].P),
                      "production_ratio": pr, "reliability": rel}
        maps[key] = th.quantile_map(df.loc[v, "model"], df.loc[v, "operation"])
        cases = [("theory", "literature", lit), ("theory", "site", rel["for"]),
                 ("theory_scaled", "site", rel["for"]), ("model", "site", rel["for"]),
                 ("operation", "site", rel["for"])]
        for lname in ("KSE", "RTS"):
            for name, src, q in cases:
                adequacy.append({"site": key, "tech": "wind", "series": name, "structural": src, "for": q,
                                 **_adequacy(ctx, lname, fitted[name], q)})
        # chronology: frequency and duration of loss-of-load events (KSE load, sequential Monte Carlo)
        load = ctx.loads["KSE"]
        cap = RES_SHARE * load.peak
        for name in ("theory", "theory_scaled", "operation"):
            plant = mc.ResPlant("wind", cap, mc.MarkovProduction(fitted[name], ctx.months),
                                mc.Structural(rel["for"], rel["mttr_h"], N_UNITS))
            r = mc.run(study.system_spec(ctx, load, [plant]), n_sims, SEED).summary()
            for a in adequacy:
                if (a["site"], a["series"], a["structural"], a["load"]) == (key, name, "site", "KSE"):
                    a.update({"mc_lole": r["lole"], "mc_lole_ci95": r["lole_ci95"], "mc_lolf": r["lolf"],
                              "mc_lolf_ci95": r["lolf_ci95"],
                              "mc_lold_h": r["lold_h"], "n_sims": r["n_sims"]})
        log.info("theory %s: PR %.3f, FOR %.4f", key, pr, rel["for"])

    pv = pv_theory(ctx)
    sites += pv["sites"]
    structure += pv["structure"]
    adequacy += pv["adequacy"]
    tests["pvdaq34"] = pv["tests"]

    transfer = transfer_to_poland(ctx, maps)
    save_csv("theory_sites.csv", pd.DataFrame(sites))
    save_csv("theory_wind_bins.csv", pd.concat(bins, ignore_index=True))
    save_csv("theory_markov.csv", pd.DataFrame(structure))
    save_csv("theory_adequacy.csv", pd.DataFrame(adequacy))
    save_csv("theory_transfer.csv", transfer)
    save_json("theory.json", {"tests": tests, "farms": {k: vars(f) for k, f in farms.FARMS.items()},
                              "pv_site": farms.PV_SITE,
                              "specific_power_w_m2": {k: th.specific_power(f) for k, f in farms.FARMS.items()}})
    return {"sites": sites, "adequacy": adequacy, "transfer": transfer}


def _local_matrix(x: pd.Series, local: pd.Series, years) -> np.ndarray:
    """(n_years, 8760) matrix of an hourly series by local standard time (leap days: last day dropped)."""
    s = pd.Series(np.asarray(x, dtype=float), index=pd.DatetimeIndex(local))
    rows = []
    for y in years:
        idx = pd.date_range(f"{y}-01-01", periods=HOURS, freq="1h")
        rows.append(s.reindex(idx).to_numpy())
    return np.vstack(rows)


def pv_theory(ctx: study.Context) -> dict:
    """PVDAQ 34: ERA5 + nameplate model against measured output; adequacy on the RTS-79 load (local time)."""
    df = th.pv_site()
    sites, structure, adequacy = [], [], []
    for_site = structural_params()["data"]["pv"]["for"]  # PVDAQ 34 forced-outage rate of the main study
    for col in ("theory", "operation_all", "operation"):
        m = th.series_metrics(df, "theory", col)
        sites.append({"site": "pvdaq34", "tech": "pv", "series": col, **{k: m[k] for k in m if k != "ratio_by_year"},
                      "ratio_year_min": min(m["ratio_by_year"].values()),
                      "ratio_year_max": max(m["ratio_by_year"].values()), "n_years": len(m["ratio_by_year"]),
                      "for_site": for_site})
    years = farms.PV_SITE["years"]
    cs = _local_matrix(df["clear_sky"], df["local"], years)
    day = cs > 0.02
    valid = _local_matrix(df["valid"].astype(float), df["local"], years) > 0.5
    pr = float(df.loc[df["valid"], "operation"].sum() / df.loc[df["valid"], "theory"].sum())
    env = np.where(day, cs, 0.0)
    fitted = {}
    for name, col, scale in (("theory", "theory", 1.0), ("theory_scaled", "theory", pr), ("operation", "operation", 1.0)):
        x = _local_matrix(df[col] * scale, df["local"], years)
        k, _ = power.clear_sky_index(np.nan_to_num(x), cs)
        ok = (day & valid).reshape(-1)
        months = np.tile(ctx.months, len(years))
        fitted[name] = mk.fit(k.reshape(-1), th.N_STATES, months=months, valid=ok)
        hom = mk.fit(k.reshape(-1), th.N_STATES, valid=ok)
        states = mk.discretize(np.where(ok, k.reshape(-1), np.nan), th.N_STATES)
        occ = np.bincount(states[states >= 0], minlength=th.N_STATES) / (states >= 0).sum()
        for s in range(th.N_STATES):
            structure.append({"site": "pvdaq34", "series": name, "state": s, "occupancy": occ[s],
                              "residence_h": float("nan"), "level": hom.levels[s]})
    for name, src, q in (("theory", "literature", LIT_PV["for"]), ("theory", "site", for_site),
                         ("theory_scaled", "site", for_site), ("operation", "site", for_site)):
        adequacy.append({"site": "pvdaq34", "tech": "pv", "series": name, "structural": src, "for": q,
                         **_adequacy(ctx, "RTS", fitted[name], q, env)})
    return {"sites": sites, "structure": structure, "adequacy": adequacy,
            "tests": {"production_ratio": pr, "for_site": for_site}}


ENERGY_RATIOS = (0.7, 0.8, 0.9, 1.0, 1.1, 1.2)


def transfer_to_poland(ctx: study.Context, maps: dict, ratios=ENERGY_RATIOS) -> pd.DataFrame:
    """Sensitivity of the Polish portfolio's wind ELCC to an error in the energy of this model.

    Every grid point's wind series is multiplied by a ratio (capped at full output), and, to see whether the
    shape of a real farm matters beyond its energy, put through each UK farm's quantile map (this model ->
    operation). ELCC of 10 % of peak, KSE load, historical weather and the monthly chain. Hypothetical: nothing
    here was measured in Poland, and the UK ratios are not an estimate of the Polish one.
    """
    load = ctx.loads["KSE"]
    h = study.hours_of(load)
    cap = RES_SHARE * load.peak
    rows = []

    def ev(case, wind, ratio):
        s = site_series.SiteSeries(ctx.portfolio.point, wind, ctx.portfolio.pv, ctx.portfolio.pv_cs, ctx.portfolio.ws)
        for model in ("hist", "seasonal"):
            d = study.dist_wind(ctx, s, cap, model, h)
            rows.append({"case": case, "energy_ratio": ratio, "model": model, "cf": float(wind.mean()),
                         "elcc_pct": 100 * elcc.elcc(ctx.copt, load.matrix, d) / cap,
                         "lole": elcc.lole(ctx.copt, load.matrix, d)})

    base = ctx.portfolio.wind
    for r in ratios:
        ev("energy ratio", np.clip(np.mean([np.clip(s.wind * r, 0, 1) for s in ctx.sites], axis=0), 0, 1), r)
    for key, fmap in maps.items():
        mapped = np.mean([np.clip(fmap(s.wind), 0, 1) for s in ctx.sites], axis=0)
        ev(f"quantile map of {key}", mapped, float(mapped.mean() / base.mean()))
    return pd.DataFrame(rows)


# --- B. weather-dependent failures and repairs ---------------------------------------------

def _farm_states(farm: farms.Farm, data_dir: Path) -> pd.DataFrame:
    e = farms.farm_era5(farm, data_dir)
    e = e[(e.index >= pd.Timestamp(farm.log_start, tz="UTC")) & (e.index < pd.Timestamp(farm.log_end, tz="UTC"))]
    cf = power.wind_cf(e["ws100"].to_numpy(), WIND_LOSSES)
    return pd.DataFrame({"cf_state": mk.discretize(cf, 5), "ws_bin": np.digitize(e["ws_hub"].to_numpy(), WS_EDGES),
                         "theory": th.catalogue_cf(e["ws_hub"].to_numpy(), farm)}, index=e.index)


def _corr_joint(pi: np.ndarray, levels: np.ndarray) -> float:
    """Correlation between the production level and the unit being up under a joint (state, up/down) law."""
    p_up = pi[:, 0].sum()
    m_x = float(pi.sum(axis=1) @ levels)
    cov = float(pi[:, 0] @ levels) - m_x * p_up
    var_x = float(pi.sum(axis=1) @ levels**2) - m_x**2
    var_u = p_up * (1 - p_up)
    return cov / np.sqrt(var_x * var_u) if var_x > 0 and var_u > 0 else float("nan")


def _fill(rates: np.ndarray, pooled: float) -> np.ndarray:
    r = np.asarray(rates, dtype=float)
    return np.where(np.isfinite(r), r, pooled)


def failures_stage(ctx: study.Context, n_sims: int = 2000, data_dir: Path = config.DATA_DIR) -> dict:
    rate_rows, tests, farm_rows, pl_rows = [], {}, [], []
    tables = {"cf_state": {}, "ws_bin": {}}
    for key, f in farms.FARMS.items():
        out, win = cached_outages(f, data_dir)
        sd = _farm_states(f, data_dir)
        hns = sd.index.as_unit("ns").asi8
        tests[key] = {}
        for cov, n_s in (("cf_state", 5), ("ws_bin", 6)):
            states = sd[cov].to_numpy()
            exp = ig.state_exposure(win, out, hns, states, n_s, MIN_OUTAGE_H)
            tab = ig.rate_table(exp)
            tables[cov][key] = tab
            cells = ig.state_exposure(win, out, hns, states, n_s, MIN_OUTAGE_H, cells="year")
            t = {}
            for what, n, T in (("failure", "failures", "up_h"), ("repair", "repairs", "down_h")):
                r = ig.poisson_lr(tab[n], tab[T])
                phi = ig.dispersion(cells, n, T)
                r.update({"dispersion": phi, "quasi_p": float(stats.chi2.sf(r["g2"] / phi, r["df"])) if r["df"] else None,
                          "quasi_trend_p": float(stats.chi2.sf(r["trend_z2"] / phi, 1))})
                t[what] = r
            tests[key][cov] = t
            for row in tab.to_dict("records"):
                rate_rows.append({"site": key, "covariate": cov,
                                  "label": WS_LABELS[int(row["state"])] if cov == "ws_bin" else str(int(row["state"])),
                                  **row})
        farm_rows += farm_integration(ctx, key, sd, tables["ws_bin"][key])
    for cov in ("cf_state", "ws_bin"):
        for n, T, what in (("failures", "up_h", "failure"), ("repairs", "down_h", "repair")):
            ok = sorted(set.intersection(*[set(v.loc[v[T] > 0, "state"]) for v in tables[cov].values()]))
            common = {k: v.set_index("state").loc[ok].reset_index() for k, v in tables[cov].items()}
            tests.setdefault("common", {}).setdefault(cov, {})[what] = ig.common_state_effect(common, n, T)
    pl_rows, mc_rows = pl_integration(ctx, tables["cf_state"], n_sims)
    save_csv("failure_rates.csv", pd.DataFrame(rate_rows))
    save_csv("integration_farms.csv", pd.DataFrame(farm_rows))
    save_csv("integration_pl.csv", pd.DataFrame(pl_rows + mc_rows))
    save_json("failure_tests.json", {"tests": tests, "min_outage_h": MIN_OUTAGE_H, "ws_bins": list(WS_LABELS),
                                     "cf_state_thresholds": mk.thresholds(5).tolist()})
    return {"rates": rate_rows, "tests": tests, "farms": farm_rows, "pl": pl_rows + mc_rows}


def farm_integration(ctx: study.Context, key: str, sd: pd.DataFrame, tab: pd.DataFrame) -> list[dict]:
    """At the farm: hub-wind-bin chain of the log period with its own rates, Kronecker vs block."""
    n_s = len(WS_LABELS)
    states = sd["ws_bin"].to_numpy()
    P = mk.normalise_rows(mk.transition_counts(states, n_s))
    Lp = ig.generator_from_P(P)
    levels = np.array([sd.loc[states == s, "theory"].mean() if (states == s).any() else 0.0 for s in range(n_s)])
    lam0 = tab["failures"].sum() / tab["up_h"].sum()
    mu0 = tab["repairs"].sum() / tab["down_h"].sum()
    lam = _fill(tab["lambda"], lam0)
    mu = _fill(tab["mu"], mu0)
    pi_k = ig.ctmc_stationary(ig.kronecker_generator(Lp, lam0, mu0)).reshape(-1, 2)
    pi_b = ig.joint_stationary(Lp, lam, mu)
    pi_p = mk.stationary(P)
    A_b = ig.availability_by_state(Lp, lam, mu)
    A_k = np.full(n_s, mu0 / (lam0 + mu0))
    A_match = np.full(n_s, float(pi_p @ A_b))
    mu_h = mu.copy()
    mu_h[[WS_LABELS.index("12-15"), WS_LABELS.index(">=15")]] = 0.0
    A_h = ig.availability_by_state(Lp, lam, mu_h)
    A_h_match = np.full(n_s, float(pi_p @ A_h))
    rows = []
    for name, A in (("independent", A_k), ("state-dependent", A_b), ("independent, same FOR", A_match),
                    ("hypothetical: no repairs at 12 m/s or more", A_h),
                    ("hypothetical: independent, same FOR", A_h_match)):
        r = {"site": key, "model": name, "mean_availability": float(pi_p @ A),
             "mean_output": float(pi_p @ (levels * A)), "output_loss_pct": 100 * float(pi_p @ (levels * (1 - A)) /
                                                                                     (pi_p @ levels))}
        for s in range(n_s):
            r[f"A_{WS_LABELS[s]}"] = float(A[s])
            r[f"pi_{WS_LABELS[s]}"] = float(pi_p[s])
        for lname in ("KSE", "RTS"):
            load = ctx.loads[lname]
            h = study.hours_of(load)
            cap = RES_SHARE * load.peak
            d = elcc.ResDistribution(cap * (levels * A)[None, None, :] * np.ones((1, h, 1)),
                                     np.broadcast_to(pi_p, (1, h, n_s)).copy(), cap)
            r[f"lole_{lname}"] = elcc.lole(ctx.copt, load.matrix, d)
            r[f"elcc_pct_{lname}"] = 100 * elcc.elcc(ctx.copt, load.matrix, d) / cap
        rows.append(r)
    # distance of the joint stationary distribution from the product of its marginals: zero for the Kronecker sum
    rows[0]["max_abs_diff_joint_vs_product"] = float(np.abs(pi_k - np.outer(pi_p, [mu0, lam0]) / (lam0 + mu0)).max())
    a = float(pi_p @ A_b)
    rows[1]["max_abs_diff_joint_vs_product"] = float(np.abs(pi_b - np.outer(pi_p, [a, 1 - a])).max())
    rows[1]["corr_output_availability"] = _corr_joint(pi_b, levels)
    return rows


STRESS = {  # hypothetical cases on top of the pooled rates of the three farms (state = portfolio CF state)
    "H1 (hypothetical) no repairs in the highest state": {"mu_factor": {3: 0.0}},
    "H2 (hypothetical) H1 + 5x failures in the highest state": {"mu_factor": {3: 0.0}, "lam_factor": {3: 5.0}},
}


def pl_integration(ctx: study.Context, tabs: dict, n_sims: int) -> tuple[list[dict], list[dict]]:
    """National portfolio (monthly chain, wind 10 % of peak): independence vs state-dependent rates."""
    model = study.fit_wind(ctx.portfolio, ctx.months)
    pooled = pd.concat(tabs.values()).groupby("state")[["failures", "repairs", "up_h", "down_h"]].sum()
    cases = {}
    for key, tab in list(tabs.items()) + [("pooled", ig.rate_table(pooled.reset_index()))]:
        lam0 = tab["failures"].sum() / tab["up_h"].sum()
        mu0 = tab["repairs"].sum() / tab["down_h"].sum()
        lam, mu = _fill(tab["lambda"], lam0), _fill(tab["mu"], mu0)
        cases[key] = (lam0, mu0, lam, mu)
    rows, mc_rows = [], []
    occ = np.stack([model.stationary(m) for m in range(1, 13)])
    for key, (lam0, mu0, lam, mu) in cases.items():
        variants = [("independent", np.full(5, lam0), np.full(5, mu0)), ("state-dependent", lam, mu)]
        if key == "pooled":
            for name, s in STRESS.items():
                lam_s, mu_s = lam.copy(), mu.copy()
                for i, fct in s.get("lam_factor", {}).items():
                    lam_s[i] *= fct
                for i, fct in s.get("mu_factor", {}).items():
                    mu_s[i] *= fct
                variants.append((name, lam_s, mu_s))
        A_sd = None
        for name, la, mu_ in variants:
            A = ig.monthly_availability(model, la, mu_)
            mean_a = float(np.mean((occ * A).sum(axis=1)))
            if name == "state-dependent":
                A_sd = mean_a
            var = [(name, la, mu_)]
            if name != "independent":
                # same mean availability, no correlation
                q = 1 - mean_a
                var.append((f"{name}: independent, same FOR", np.full(5, q * mu0 / (1 - q)), np.full(5, mu0)))
            for vname, la2, mu2 in var:
                r = {"rates": key, "variant": vname, **{f"lambda_{i}": la2[i] for i in range(5)},
                     **{f"mu_{i}": mu2[i] for i in range(5)}}
                A2 = ig.monthly_availability(model, la2, mu2)
                r["mean_availability"] = float(np.mean((occ * A2).sum(axis=1)))
                r["output_weighted_availability"] = float(np.mean((occ * A2 * model.levels).sum(axis=1)) /
                                                          np.mean((occ * model.levels).sum(axis=1)))
                for lname in ("KSE", "RTS"):
                    load = ctx.loads[lname]
                    h = study.hours_of(load)
                    cap = RES_SHARE * load.peak
                    d = ig.integrated_distribution(model, ctx.months, cap, la2, mu2, 1, h)
                    ind = elcc.indices(ctx.copt, load.matrix, d)
                    r[f"lole_{lname}"] = ind["lole"]
                    r[f"eens_{lname}"] = ind["eens"]
                    r[f"elcc_pct_{lname}"] = 100 * elcc.elcc(ctx.copt, load.matrix, d, xtol=0.01) / cap
                rows.append(r)
        log.info("integration %s: mean availability %.4f (state-dependent)", key, A_sd)
    # Monte Carlo check with common random numbers (KSE load): pooled rates
    load = ctx.loads["KSE"]
    cap = RES_SHARE * load.peak
    lam0, mu0, lam, mu = cases["pooled"]
    mc_cases = [("independent", np.full(5, lam0), np.full(5, mu0)), ("state-dependent", lam, mu)]
    for name, s in STRESS.items():
        la, m_ = lam.copy(), mu.copy()
        for i, fct in s.get("lam_factor", {}).items():
            la[i] *= fct
        for i, fct in s.get("mu_factor", {}).items():
            m_[i] *= fct
        mc_cases.append((name, la, m_))
    base_margin = mc.run(study.system_spec(ctx, load, []), n_sims, SEED, keep_margin=True).margin
    for name, la, m_ in mc_cases:
        prod = ig.IntegratedProduction(model, ctx.months, la, m_, N_UNITS)
        res = mc.run(study.system_spec(ctx, load, [mc.ResPlant("wind", cap, prod)]), n_sims, SEED, keep_margin=True)
        s = res.summary()
        mc_rows.append({"rates": "pooled", "variant": f"MC: {name}", "lole_KSE": s["lole"],
                        "lole_ci95_KSE": s["lole_ci95"], "lolf_KSE": s["lolf"], "lold_h_KSE": s["lold_h"],
                        "elcc_pct_KSE": 100 * elcc.mc_elcc(base_margin, res.margin, cap) / cap, "n_sims": n_sims})
    return rows, mc_rows
