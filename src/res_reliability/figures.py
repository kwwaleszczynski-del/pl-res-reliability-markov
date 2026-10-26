"""README charts (PNG) rendered from ``results/``."""

from __future__ import annotations

import json
import logging
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

from res_reliability.config import IMG_DIR, RESULTS_DIR  # noqa: E402
from res_reliability.geo import poland_outline  # noqa: E402

log = logging.getLogger(__name__)

SURFACE = "#fcfcfb"
INK, INK2, MUTED = "#0b0b0b", "#52514e", "#898781"
GRID, BASE = "#e1e0d9", "#c3c2b7"
C1, C2, C3, C4, C5 = "#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"
BLUES = LinearSegmentedColormap.from_list(
    "blues", ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#2a78d6", "#1c5cab", "#104281", "#0d366b"])
RAMP4 = ["#6da7ec", "#2a78d6", "#1c5cab", "#0d366b"]  # ordered magnitudes (durations, fleet sizes, shares)
MODEL_COLORS = {"data": C1, "hist": C1, "homogeneous": C2, "seasonal": C3, "regime": C4, "regime_t": C5}
MODEL_NAMES = {"data": "ERA5 series", "hist": "historical weather", "homogeneous": "Markov, homogeneous",
               "seasonal": "Markov, monthly", "regime": "regime switching", "regime_t": "regimes + temperature"}
MONTHS = list("JFMAMJJASOND")


def style() -> None:
    plt.rcParams.update({
        "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
        "axes.edgecolor": BASE, "axes.labelcolor": INK2, "axes.titlecolor": INK, "axes.titlesize": 11,
        "axes.titleweight": "semibold", "axes.titlelocation": "left", "axes.labelsize": 9.5,
        "xtick.color": MUTED, "ytick.color": MUTED, "xtick.labelsize": 8.5, "ytick.labelsize": 8.5,
        "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8, "axes.axisbelow": True,
        "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
        "legend.fontsize": 8.5, "legend.labelcolor": INK2, "font.family": "sans-serif",
        "font.sans-serif": ["Segoe UI", "DejaVu Sans", "Arial"], "lines.linewidth": 2,
        "lines.solid_capstyle": "round", "figure.dpi": 110, "savefig.dpi": 150,
    })


def save(fig, name: str, img_dir: Path = IMG_DIR) -> Path:
    img_dir.mkdir(parents=True, exist_ok=True)
    p = img_dir / name
    fig.savefig(p, bbox_inches="tight")
    plt.close(fig)
    log.info("wrote %s", p)
    return p


def _load_json(name, results_dir=RESULTS_DIR):
    return json.loads((results_dir / name).read_text(encoding="utf-8"))


def _cell_text_color(v: float, vmax: float) -> str:
    return "#ffffff" if v / vmax > 0.55 else INK


def transition_matrices(results_dir=RESULTS_DIR):
    m = _load_json("markov_models.json", results_dir)["portfolio"]
    fig, axes = plt.subplots(1, 2, figsize=(9.4, 4.2))
    for ax, tech, title in zip(axes, ("wind", "pv"), ("Wind, capacity factor states", "PV, clear-sky index states")):
        P = np.array(m[tech]["P_homogeneous"])
        pi = np.array(m[tech]["pi_homogeneous"])
        keep = pi > 0
        P = P[np.ix_(keep, keep)]
        lv = np.array(m[tech]["levels_mean"])[keep]
        ax.imshow(P, cmap=BLUES, vmin=0, vmax=1)
        ax.grid(False)
        n = P.shape[0]
        for i in range(n):
            for j in range(n):
                ax.text(j, i, f"{P[i, j]:.2f}", ha="center", va="center", fontsize=8.5,
                        color=_cell_text_color(P[i, j], 1.0))
        labels = [f"{100 * v:.0f}%" for v in lv]
        ax.set_xticks(range(n), labels)
        ax.set_yticks(range(n), [f"{a}\nπ={p:.2f}" for a, p in zip(labels, pi[keep])])
        ax.set_xlabel("to state (mean level)")
        ax.set_ylabel("from state")
        ax.set_title(title)
        for s in ax.spines.values():
            s.set_visible(False)
    fig.suptitle("Hourly transition matrices of the national portfolio (homogeneous, 5 states; empty states hidden)",
                 x=0.01, ha="left", fontsize=11.5, color=INK, fontweight="semibold")
    fig.tight_layout()
    return save(fig, "transition_matrices.png")


def seasonality_and_persistence(results_dir=RESULTS_DIR):
    v = _load_json("markov_models.json", results_dir)["validation"]
    fig, axes = plt.subplots(2, 2, figsize=(10, 6.6))
    for ax, tech, title in ((axes[0, 0], "wind", "Wind: monthly mean capacity factor"),
                            (axes[0, 1], "pv", "PV: monthly mean capacity factor")):
        for key in ("data", "homogeneous", "seasonal"):
            ax.plot(range(12), 100 * np.array(v[tech]["monthly_mean"][key]), color=MODEL_COLORS[key],
                    marker="o", markersize=4, label=MODEL_NAMES[key])
        ax.set_xticks(range(12), MONTHS)
        ax.set_ylabel("%")
        ax.set_ylim(bottom=0)
        ax.set_title(title)
    axes[0, 0].legend(loc="lower left")
    ax = axes[1, 0]
    lags = np.arange(len(v["wind"]["acf"]["data"]))
    for key in ("data", "homogeneous", "seasonal"):
        ax.plot(lags, v["wind"]["acf"][key], color=MODEL_COLORS[key], label=MODEL_NAMES[key])
    ax.set_xticks([0, 12, 24, 36, 48, 60, 72])
    ax.set_xlabel("lag, h")
    ax.set_title("Wind: autocorrelation of hourly output")
    ax.set_ylim(0, 1)
    ax = axes[1, 1]
    sp = v["wind"]["low_spells_per_year_longer_than"]
    x = np.arange(len(sp["hours"]))
    w = 0.26
    for k, key in enumerate(("data", "homogeneous", "seasonal")):
        ax.bar(x + (k - 1) * w, sp[key], width=w * 0.9, color=MODEL_COLORS[key], label=MODEL_NAMES[key])
    ax.set_yscale("log")
    ax.set_xticks(x, [f"≥{h}" for h in sp["hours"]])
    ax.set_xlabel("spell length, h")
    ax.set_ylabel("spells per year (log)")
    ax.set_title("Wind: low-output spells (portfolio CF < 12.5 %)")
    longest = v["wind"]["longest_low_spell_h"]
    ax.text(0.98, 0.95, "longest: data {data} h, homogeneous {homogeneous} h, monthly {seasonal} h".format(**longest),
            transform=ax.transAxes, ha="right", va="top", fontsize=8, color=INK2)
    fig.tight_layout()
    return save(fig, "markov_validation.png")


def variants(results_dir=RESULTS_DIR):
    df = pd.read_csv(results_dir / "variants.csv")
    order = ["BAU", "LIT", "DATA", "NOFOR", "HOM", "HIST", "THESIS", "DATA-W", "DATA-PV"]
    names = {"BAU": "conventional only", "LIT": "Markov monthly, FOR literature", "DATA": "Markov monthly, FOR data",
             "NOFOR": "Markov monthly, no outages", "HOM": "Markov homogeneous, FOR data",
             "HIST": "historical weather, FOR data", "THESIS": "thesis implementation",
             "DATA-W": "wind only (FOR data)", "DATA-PV": "PV only (FOR data)"}
    loads = [l for l in ("RTS", "KSE") if l in set(df["load"])]
    fig, axes = plt.subplots(1, len(loads), figsize=(5.2 * len(loads), 4.6), sharey=True, squeeze=False)
    for ax, load in zip(axes[0], loads):
        d = df[df["load"] == load].set_index("variant").reindex(order).dropna(subset=["lole"])
        y = np.arange(len(d))[::-1]
        ax.barh(y, d["lole"], height=0.55, color=[MUTED if i == "BAU" else C1 for i in d.index])
        ax.errorbar(d["lole"], y, xerr=d["lole_ci95"], fmt="none", ecolor=INK2, elinewidth=1, capsize=2)
        ax.scatter(d["lole_analytic"], y, marker="|", s=120, color=INK, zorder=3, label="analytical expectation")
        right = np.maximum(d["lole"] + d["lole_ci95"], d["lole_analytic"])
        for yy, val, r in zip(y, d["lole"], right):
            ax.text(r + d["lole"].max() * 0.03, yy, f"{val:.2f}", va="center", fontsize=8, color=INK2)
        ax.set_xlim(0, right.max() * 1.15)
        ax.set_yticks(y, [names[i] for i in d.index])
        ax.set_xlabel("LOLE, h/year (Monte Carlo mean, 95 % CI)")
        title = "IEEE RTS-79 load" if load == "RTS" else "KSE load (PSE), scaled to LOLE 3 h"
        ax.set_title(title)
        ax.grid(axis="y", visible=False)
    axes[0][0].legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=1)
    fig.suptitle("Adequacy of the RTS-79 generating system with 10 % wind + 10 % PV (of peak)",
                 x=0.01, ha="left", fontsize=11.5, color=INK, fontweight="semibold")
    fig.tight_layout()
    return save(fig, "variants_lole.png")


def _outline(ax):
    o = poland_outline()
    ax.plot(o[:, 0], o[:, 1], color=BASE, linewidth=1.2)
    ax.set_aspect(1 / np.cos(np.radians(52)))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    for s in ax.spines.values():
        s.set_visible(False)


def elcc_map(results_dir=RESULTS_DIR, load="KSE", model="hist"):
    df = pd.read_csv(results_dir / "elcc_points.csv")
    if load not in set(df["load"]):
        load = df["load"].iloc[0]
    d = df[(df["load"] == load) & (df["model"] == model)]
    w, p = d[d["tech"] == "wind"], d[d["tech"] == "pv"]
    panels = [(w, "cf_pct", "Wind capacity factor, %"), (w, "elcc_pct", "Wind ELCC, % of capacity"),
              (p, "elcc_pct", "PV ELCC, % of capacity")]
    fig, axes = plt.subplots(1, 3, figsize=(12.5, 4.3))
    for ax, (data, col, title) in zip(axes, panels):
        _outline(ax)
        sc = ax.scatter(data["lon"], data["lat"], c=data[col], cmap=BLUES, s=170, edgecolors=SURFACE,
                        linewidths=2, zorder=3)
        for _, r in data.iterrows():
            ax.text(r["lon"], r["lat"] - 0.33, f"{r[col]:.1f}", ha="center", va="top", fontsize=7.5, color=INK2)
        cb = fig.colorbar(sc, ax=ax, shrink=0.7, pad=0.02)
        cb.outline.set_visible(False)
        cb.ax.tick_params(labelsize=8, colors=MUTED)
        ax.set_title(title)
    lname = "KSE load" if load == "KSE" else "RTS-79 load"
    fig.suptitle(f"ELCC of a plant of 10 % of peak load at each ERA5 grid point ({lname}, {MODEL_NAMES[model]})",
                 x=0.01, ha="left", fontsize=11.5, color=INK, fontweight="semibold")
    fig.tight_layout()
    return save(fig, "elcc_map.png")


def elcc_curves(results_dir=RESULTS_DIR):
    c = pd.read_csv(results_dir / "elcc_curves.csv")
    kwd = _load_json("kwd.json", results_dir)["values"]
    loads = [l for l in ("KSE", "RTS") if l in set(c["load"])]
    fig, axes = plt.subplots(1, len(loads), figsize=(5.6 * len(loads), 4.4), sharey=True, squeeze=False)
    colors = {"wind": C1, "pv": C2, "wind+pv": C3}
    for ax, load in zip(axes[0], loads):
        d = c[c["load"] == load]
        for mix in ("wind", "pv", "wind+pv"):
            for model, ls in (("hist", "-"), ("seasonal", "--")):
                s = d[(d["mix"] == mix) & (d["model"] == model)]
                ax.plot(100 * s["share"], s["elcc_pct"], color=colors[mix], linestyle=ls,
                        label=f"{mix.replace('pv', 'PV')}, {MODEL_NAMES[model]}")
        if load == "KSE":
            wv = [k["wind"] for k in kwd]
            pv = [k["pv"] for k in kwd]
            ax.axhspan(min(wv), max(wv), color=C1, alpha=0.10, linewidth=0)
            ax.axhspan(min(pv), max(pv), color=C2, alpha=0.12, linewidth=0)
            ax.text(1, min(wv) - 0.3, f"KWD onshore wind {min(wv):.1f}-{max(wv):.1f} %", ha="left", va="top",
                    fontsize=8, color=INK2)
            ax.text(1, min(pv) - 0.3, f"KWD solar {min(pv):.1f}-{max(pv):.1f} %", ha="left", va="top",
                    fontsize=8, color=INK2)
        ax.set_xlabel("installed RES capacity, % of peak load")
        ax.set_ylabel("ELCC, % of installed capacity")
        ax.set_xlim(0, 102)
        ax.set_ylim(bottom=0)
        ax.set_title("KSE load (PSE), conventional LOLE 3 h/yr" if load == "KSE" else "IEEE RTS-79 load")
    axes[0][-1].legend(loc="upper right", ncol=1)
    fig.tight_layout()
    return save(fig, "elcc_curves.png")


def lole_vs_share(results_dir=RESULTS_DIR):
    df = pd.read_csv(results_dir / "lole_vs_share.csv")
    loads = [l for l in ("KSE", "RTS") if l in set(df["load"])]
    fig, axes = plt.subplots(1, len(loads), figsize=(5.4 * len(loads), 4.0), squeeze=False)
    colors = {"wind": C1, "pv": C2, "wind+pv": C3}
    for ax, load in zip(axes[0], loads):
        d = df[df["load"] == load]
        base = float(d[d["mix"] == "none"]["lole"].iloc[0])
        for mix in ("wind", "pv", "wind+pv"):
            s = d[(d["mix"] == mix) & (d["model"] == "hist")]
            ax.plot([0] + list(100 * s["share"]), [base] + list(s["lole"]), color=colors[mix],
                    label=mix.replace("pv", "PV"))
        ax.set_yscale("log")
        ax.set_xlabel("RES capacity added, % of peak load")
        ax.set_ylabel("LOLE, h/year (log)")
        ax.set_title("KSE load" if load == "KSE" else "IEEE RTS-79 load")
    axes[0][0].legend()
    fig.tight_layout()
    return save(fig, "lole_vs_share.png")


def sensitivity(results_dir=RESULTS_DIR):
    st = pd.read_csv(results_dir / "sensitivity_states.csv")
    ln = pd.read_csv(results_dir / "sensitivity_length.csv")
    fr = pd.read_csv(results_dir / "sensitivity_for.csv")
    yr = pd.read_csv(results_dir / "sensitivity_years.csv")
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7))
    ax = axes[0, 0]
    ref = float(st[st["model"] == "hist"]["elcc_pct"].iloc[0])
    for (levels, model), color, ls in ((("centre", "homogeneous"), C2, "--"), (("mean", "homogeneous"), C2, "-"),
                                       (("centre", "seasonal"), C3, "--"), (("mean", "seasonal"), C3, "-")):
        s = st[(st["levels"] == levels) & (st["model"] == model)]
        ax.plot(s["n_states"], s["elcc_pct"], color=color, linestyle=ls, marker="o", markersize=4,
                label=f"{model}, {'bin-centre' if levels == 'centre' else 'conditional-mean'} levels")
    ax.axhline(ref, color=C1, linewidth=1.5, label="historical weather")
    ax.set_xlabel("number of Markov states")
    ax.set_ylabel("wind ELCC, %")
    ax.set_title("Number of states and state levels")
    ax.legend(fontsize=7.5)
    ax = axes[0, 1]
    g = ln.groupby("years")["elcc_pct"]
    ax.errorbar(g.mean().index, g.mean(), yerr=[g.mean() - g.min(), g.max() - g.mean()], fmt="o",
                color=C3, ecolor=C3, elinewidth=1.5, capsize=3, markersize=5)
    ax.set_xticks(g.mean().index)
    ax.set_xlabel("years of data used to fit the monthly chain")
    ax.set_ylabel("wind ELCC, % (range over windows)")
    ax.set_title("Length of the weather series")
    ax = axes[1, 0]
    for mode, color, label in (("units", C1, "independent turbines (sequential)"),
                               ("plant_bernoulli", C2, "portfolio as one plant, hourly outages")):
        s = fr[fr["mode"] == mode]
        ax.plot(100 * s["for"], s["elcc_pct"], color=color, marker="o", markersize=4, label=label)
    ax.set_xlabel("wind forced outage rate, %")
    ax.set_ylabel("wind ELCC, %")
    ax.set_title("Structural reliability (FOR)")
    ax.legend()
    ax = axes[1, 1]
    for tech, color in (("wind", C1), ("pv", C2)):
        s = yr[yr["tech"] == tech]
        ax.plot(s["year"], s["elcc_pct"], color=color, marker="o", markersize=4, label=f"{'PV' if tech == 'pv' else 'wind'} ELCC")
    ax.set_xticks(yr["year"].unique()[::2])
    ax.set_ylabel("ELCC, %")
    ax.set_title("Single weather years (historical, matching load year)")
    ax.legend(loc="center right")
    pair = results_dir / "sensitivity_weather_load_pairing.csv"
    if pair.exists():
        pr = pd.read_csv(pair)
        w = pr[pr["tech"] == "wind"]
        ax.text(0.02, 0.04, "all ten years, wind: load and weather of the same year {:.1f} %;\n"
                "load paired with other weather years {:.1f} % ({:.1f}-{:.1f})".format(
                    w[w["matched"]]["elcc_pct"].iloc[0], w[~w["matched"]]["elcc_pct"].mean(),
                    w[~w["matched"]]["elcc_pct"].min(), w[~w["matched"]]["elcc_pct"].max()),
                transform=ax.transAxes, fontsize=8, color=INK2, va="bottom")
    fig.tight_layout()
    return save(fig, "sensitivity.png")


def structural_fig(results_dir=RESULTS_DIR):
    s = _load_json("structural.json", results_dir)
    dur = pd.read_csv(results_dir / "kelmarsh_forced_outage_durations.csv")["duration_h"].to_numpy()
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.9))
    ax = axes[0]
    p = s["parameters"]
    labels = ["wind, literature", "wind, Kelmarsh", "PV, literature", "PV, PVDAQ 34"]
    vals = [p["literature"]["wind"]["for"], p["data"]["wind"]["for"], p["literature"]["pv"]["for"], p["data"]["pv"]["for"]]
    y = np.arange(4)[::-1]
    ax.barh(y, 100 * np.array(vals), height=0.55, color=[MUTED, C1, MUTED, C2])
    for yy, v in zip(y, vals):
        ax.text(100 * v + 0.1, yy, f"{100 * v:.2f} %", va="center", fontsize=8.5, color=INK2)
    ax.set_yticks(y, labels)
    ax.set_xlabel("forced outage rate, %")
    ax.set_title("Literature vs open event logs")
    ax.grid(axis="y", visible=False)
    ax = axes[1]
    bins = np.logspace(-2, np.log10(max(dur.max(), 200)), 30)
    ax.hist(dur, bins=bins, color=C1, rwidth=0.85)
    ax.set_xscale("log")
    ax.set_xlabel("forced outage duration, h (log)")
    ax.set_ylabel("events")
    k = s["kelmarsh"]["forced_ge_1h"]
    ax.set_title(f"Kelmarsh forced outages ({s['kelmarsh']['all_forced']['events']} events, 6 turbines)")
    ax.text(0.98, 0.95, f"≥1 h: {k['events']} events, MTTR {k['mttr_h']:.1f} h, MTTF {k['mttf_h']:.0f} h",
            transform=ax.transAxes, ha="right", va="top", fontsize=8, color=INK2)
    fig.tight_layout()
    return save(fig, "structural.png")


WEATHER_ORDER = ("data", "seasonal", "regime", "regime_t")
STORAGE_ORDER = ("hist", "seasonal", "regime", "regime_t")


def weather_models(results_dir=RESULTS_DIR):
    m = _load_json("weather_models.json", results_dir)["metrics"]
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.2))
    k = np.arange(1, 11)
    for ax, key, title in ((axes[0, 0], "0.15:all", "Wind lulls, whole year (daily CF < 15 %)"),
                           (axes[0, 1], "0.15:winter", "Wind lulls, October-March (daily CF < 15 %)")):
        for model in WEATHER_ORDER:
            y = np.array(m[model]["lulls"][key]["at_least_per_year"], dtype=float)
            ax.plot(k, np.where(y > 0, y, np.nan), color=MODEL_COLORS[model], marker="o", markersize=4,
                    label=MODEL_NAMES[model])
        ax.set_yscale("log")
        ax.set_xticks(k, [str(i) for i in k])
        ax.set_xlabel("lulls lasting at least ... days")
        ax.set_ylabel("lulls per year (log)")
        ax.set_title(title)
    axes[0, 0].legend(loc="lower left")
    ax = axes[1, 0]
    lags = np.arange(len(m["data"]["acf"]))
    for model in WEATHER_ORDER:
        ax.plot(lags, m[model]["acf"], color=MODEL_COLORS[model], label=MODEL_NAMES[model])
    ax.set_xticks([0, 24, 48, 72, 96, 120, 144, 168])
    ax.set_xlabel("lag, h")
    ax.set_ylim(0, 1)
    ax.set_title("Wind: autocorrelation of hourly output over a week")
    ax = axes[1, 1]
    x = np.arange(3)
    wdt = 0.2
    for i, model in enumerate(WEATHER_ORDER):
        v = 100 * np.array(m[model]["djf_by_temperature_tercile"]["share_days_below_0.15"])
        ax.bar(x + (i - 1.5) * wdt, v, width=wdt * 0.9, color=MODEL_COLORS[model], label=MODEL_NAMES[model])
    ax.set_xticks(x, ["coldest third", "middle third", "warmest third"])
    ax.set_xlabel("winter days (DJF) by temperature anomaly")
    ax.set_ylabel("days with wind CF < 15 %, %")
    ax.set_title("Cold and calm: winter lulls by temperature")
    ax.grid(axis="x", visible=False)
    fig.suptitle("Weather models against the ERA5 national wind portfolio (200 synthetic years per model)",
                 x=0.01, ha="left", fontsize=11.5, color=INK, fontweight="semibold")
    fig.tight_layout()
    return save(fig, "weather_models.png")


def _kwd_lines(ax, kwd, side="right", x=None, below=False):
    x0, x1 = ax.get_xlim()
    x = (x1 if side == "right" else x0) if x is None else x
    for k in kwd:
        ax.axhline(k["kwd"], color=MUTED, linewidth=1, linestyle=":")
        ax.text(x, k["kwd"] + (-1.2 if below else 1.2), f"KWD {k['auction'][-4:]}: {k['kwd']:.2f} %", ha=side,
                va="top" if below else "bottom", fontsize=7.5, color=INK2)


def storage_elcc_fig(results_dir=RESULTS_DIR):
    df = pd.read_csv(results_dir / "storage_elcc.csv")
    kwd = _load_json("kwd_storage.json", results_dir)["values"]
    los = df[df["policy"] == "losses"]
    fig, axes = plt.subplots(1, 3, figsize=(13.8, 4.6), sharey=True)
    ax = axes[0]
    d = los[(los["model"] == "hist") & (los["res_share"] == 0.10)]
    for color, (p, g) in zip(RAMP4, d.groupby("power_pct_peak")):
        ax.plot(g["duration_h"], g["elcc_pct"], color=color, marker="o", markersize=5, label=f"fleet {p:g} % of peak")
    ax.set_xscale("log", base=2)
    ax.set_xticks([1, 2, 4, 8], ["1", "2", "4", "8"])
    ax.set_xlim(0.85, 9.5)
    ax.set_ylim(0, 105)
    ax.set_xlabel("storage duration, h")
    ax.set_ylabel("ELCC, % of storage power")
    ax.set_title("Duration and fleet size (10 % wind + 10 % PV)")
    _kwd_lines(ax, kwd, "left")
    ax.legend(loc="lower right", bbox_to_anchor=(1.0, 0.17), fontsize=7.5)
    ax = axes[1]
    sw = pd.read_csv(results_dir / "storage_fleet_sweep.csv")
    sw = sw[sw["model"] == "hist"]
    for color, (dur, g) in zip(RAMP4[1:], sw.groupby("duration_h")):
        ax.plot(g["power_pct_peak"], g["marginal_elcc_pct"], color=color, marker="o", markersize=5,
                label=f"{dur:g} h")
    gw = _load_json("kwd_storage.json", results_dir)["kse_peak_mw"] / 1000
    ticks = [0, 10, 20, 30, 40, 50]
    ax.set_xticks(ticks, [f"{t}\n{t * gw / 100:.1f} GW" for t in ticks])
    ax.set_xlim(0, 52)
    ax.set_xlabel("storage fleet, % of peak load (GW at KSE scale)")
    ax.set_title("Marginal ELCC of a growing fleet (2 x 25 % RES)")
    _kwd_lines(ax, kwd, "left", x=24, below=True)
    ax.legend(loc="upper right", fontsize=7.5)
    ax = axes[2]
    d = los[(los["duration_h"] == 4) & (los["res_share"] == 0.25)]
    for model in STORAGE_ORDER:
        g = d[d["model"] == model]
        ax.plot(g["power_pct_peak"], g["elcc_pct"], color=MODEL_COLORS[model], marker="o", markersize=5,
                label=MODEL_NAMES[model])
    ax.set_xlim(0, 21)
    ax.set_xticks([1, 5, 10, 20])
    ax.set_xlabel("4 h storage fleet, % of peak load")
    ax.set_title("Weather model (25 % wind + 25 % PV)")
    _kwd_lines(ax, kwd)
    ax.legend(loc="lower left", fontsize=7.5)
    fig.suptitle("Capacity value of storage: KSE load, system at LOLE 3 h/yr, round trip 85 %, historical weather "
                 "unless stated",
                 x=0.01, ha="left", fontsize=11.5, color=INK, fontweight="semibold")
    fig.tight_layout()
    return save(fig, "storage_elcc.png")


def storage_lulls_fig(results_dir=RESULTS_DIR):
    df = pd.read_csv(results_dir / "storage_lulls.csv")
    none = df[df["model"] == "none"].iloc[0]
    sel = (df["policy"] == "losses") & np.isinf(df["duration_h"])
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.3))
    for model in STORAGE_ORDER:
        g = df[sel & (df["model"] == model)]
        axes[0].plot([0] + list(100 * g["energy_share"]), [none["lole_no_storage"]] + list(g["lole_no_storage"]),
                     color=MODEL_COLORS[model], marker="o", markersize=5, label=MODEL_NAMES[model])
        axes[1].plot(100 * g["energy_share"], g["energy_gwh_kse"], color=MODEL_COLORS[model], marker="o",
                     markersize=5, label=MODEL_NAMES[model])
        if model in ("hist", "seasonal"):  # the two ends; the table in the README has every value
            last = g.iloc[-1]
            axes[1].text(100 * last["energy_share"] + 0.8, last["energy_gwh_kse"], f"{last['energy_gwh_kse']:.0f}",
                         va="center", fontsize=8, color=INK2)
    for ax in axes:
        ax.set_yscale("log")
        ax.set_xlabel("wind + PV share of annual demand, %")
    axes[0].set_ylabel("LOLE without storage, h/yr (log)")
    axes[0].set_title("Demand growth covered by wind and PV energy")
    axes[0].legend(loc="upper left")
    axes[1].set_xlim(8, 44)
    axes[1].set_ylabel("GWh at the scale of the 2024 KSE peak (log)")
    axes[1].set_title("Storage energy that keeps LOLE at today's level")
    fig.suptitle("How much storage the lulls need (power not limited, round trip 85 %)",
                 x=0.01, ha="left", fontsize=11.5, color=INK, fontweight="semibold")
    fig.tight_layout()
    return save(fig, "storage_lulls.png")


VERIFY_METRICS = {  # short labels of the block-test metrics, top to bottom
    "lulls_ge3_015": "lulls 3+ d < 15 %", "lulls_ge5_015": "lulls 5+ d < 15 %", "winter_ge3_015": "winter lulls 3+ d",
    "lulls_ge3_010": "lulls 3+ d < 10 %", "wind_pv_ge3_010": "wind + PV lulls 3+ d", "longest_015": "longest lull",
    "min7_wind": "weakest week (CF)", "spells_ge48h": "hourly spells 48 h+", "acf24": "autocorrelation 24 h",
    "lull_lengths_chi2": "lull lengths (chi-square)", "cold_calm_share": "cold winter days calm",
    "top100_wind": "wind in top-100 load hours", "wind_elcc": "wind ELCC",
}


def verification_fig(results_dir=RESULTS_DIR):
    v = _load_json("verification.json", results_dir)
    boot = pd.DataFrame(v["bootstrap"])
    blocks = pd.read_csv(results_dir / "verification_blocks.csv")
    loyo = pd.read_csv(results_dir / "verification_loyo.csv")
    fig, axes = plt.subplots(2, 2, figsize=(11, 8.2))
    for ax, q, xlabel, title in (
            (axes[0, 0], "elcc_no_storage:wind", "wind ELCC without storage, %", "Wind ELCC (10 % of peak)"),
            (axes[0, 1], "lull_energy_gwh:30%", "storage energy, GWh at KSE scale",
             "Storage that keeps LOLE (30 % wind + PV)")):
        d = boot[boot["quantity"] == q].set_index("model")
        for k, m in enumerate(STORAGE_ORDER[::-1]):
            lo, hi = d.loc[m, "ci90"]
            ax.plot([lo, hi], [k, k], color=MODEL_COLORS[m], linewidth=3, solid_capstyle="round")
            ax.scatter(d.loc[m, "point"], k, color=MODEL_COLORS[m], s=55, zorder=3, edgecolors=SURFACE, linewidths=1.5)
            ax.text(hi, k + 0.18, f"  {d.loc[m, 'point']:.1f}  ({lo:.1f}-{hi:.1f})", fontsize=8, color=INK2, va="bottom")
        ax.set_yticks(range(4), [MODEL_NAMES[m] for m in STORAGE_ORDER[::-1]])
        ax.set_ylim(-0.6, 3.8)
        ax.set_xlabel(xlabel)
        ax.set_title(title + ", 90 % bootstrap interval")
        ax.grid(axis="y", visible=False)
    ax = axes[1, 0]
    ax.axvspan(5, 95, color=GRID, alpha=0.6, linewidth=0)
    names = list(VERIFY_METRICS)[::-1]
    for i, m in enumerate(WEATHER_ORDER[1:]):
        d = blocks[blocks["model"] == m].set_index("metric")
        ax.scatter(d.loc[names, "percentile"], np.arange(len(names)) + (i - 1) * 0.22, color=MODEL_COLORS[m], s=28,
                   zorder=3, label=MODEL_NAMES[m])
    ax.set_yticks(range(len(names)), [VERIFY_METRICS[n] for n in names])
    ax.set_xlim(-3, 103)
    ax.set_xlabel("percentile of ERA5 2015-2024 among 100 simulated decades")
    ax.set_title("Ten-year block test (shaded: 5-95 %)")
    ax.grid(axis="y", visible=False)
    ax.legend(loc="upper left", fontsize=7.5)
    ax = axes[1, 1]
    e = loyo[loyo["metric"] == "wind_elcc"]
    lims = [e[["observed", "model_mean"]].min().min() - 1, e[["observed", "model_mean"]].max().max() + 1]
    ax.plot(lims, lims, color=BASE, linewidth=1)
    for m in WEATHER_ORDER[1:]:
        g = e[e["model"] == m]
        r = np.corrcoef(g["observed"], g["model_mean"])[0, 1]
        ax.scatter(g["observed"], g["model_mean"], color=MODEL_COLORS[m], s=36, label=f"{MODEL_NAMES[m]} (r = {r:.2f})",
                   edgecolors=SURFACE, linewidths=1)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_xlabel("wind ELCC of the held-out year, historical weather, %")
    ax.set_ylabel("predicted by the model fitted on the other 9 years, %")
    ax.set_title("Leave one year out: wind ELCC of every year")
    ax.legend(loc="lower right", fontsize=7.5)
    fig.suptitle("Verification: uncertainty from ten weather years, simulated decades, held-out years",
                 x=0.01, ha="left", fontsize=11.5, color=INK, fontweight="semibold")
    fig.tight_layout()
    return save(fig, "verification.png")


def make_all(results_dir=RESULTS_DIR) -> list[Path]:
    style()
    out = []
    for fn in (transition_matrices, seasonality_and_persistence, variants, elcc_map, elcc_curves, lole_vs_share,
               sensitivity, structural_fig, weather_models, storage_elcc_fig, storage_lulls_fig, verification_fig):
        try:
            out.append(fn(results_dir))
        except FileNotFoundError as exc:
            log.warning("%s skipped: %s", fn.__name__, exc)
    return out
