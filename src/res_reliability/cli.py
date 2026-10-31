"""Reliability of wind and PV in Poland: Markov models, sequential Monte Carlo, ELCC.

Usage:
    python -m res_reliability fetch-weather --budget 2400   # ERA5 grid from Open-Meteo (cached)
    python -m res_reliability fetch-load                    # KSE demand from PSE (cached)
    python -m res_reliability fetch-farms                   # ERA5 at the open wind farms and the PV site (cached)
    python -m res_reliability run                           # all stages -> results/
    python -m res_reliability run --stages markov elcc --n-sims 500
    python -m res_reliability run --stages weather storage --jobs 6   # regime models, storage
    python -m res_reliability --data-dir data/long fetch-weather --start-year 1980 --end-year 2014 --budget 8000
    python -m res_reliability --data-dir data/long fetch-load --start-year 2010 --end-year 2014
    python -m res_reliability run --stages long             # models fitted on 2015-2024 against 1980-2014
    python -m res_reliability run --stages theory failures  # open farm data: theory vs SCADA, failure rates
    python -m res_reliability figures                       # README charts -> docs/img/
    python -m res_reliability dashboard                     # docs/index.html
"""

from __future__ import annotations

import argparse
import logging
import time
from pathlib import Path

from res_reliability import config

STAGES = ["structural", "load", "markov", "variants", "elcc", "sensitivity", "weather", "storage", "verify",
          "theory", "failures"]


def _fetch_weather(args: argparse.Namespace) -> int:
    from res_reliability import weather

    got = weather.fetch_all(start=f"{args.start_year}-01-01", end=f"{args.end_year}-12-31", budget=args.budget,
                            pause_s=args.pause, data_dir=args.data_dir)
    logging.getLogger(__name__).info("%d of %d grid points cached", len(got), len(weather.GRID))
    return 0


def _fetch_load(args: argparse.Namespace) -> int:
    from res_reliability import pse

    s = pse.fetch(list(range(args.start_year, args.end_year + 1)), data_dir=args.data_dir)
    logging.getLogger(__name__).info("KSE demand: %d hours, %s .. %s", len(s), s.index[0], s.index[-1])
    return 0


def _fetch_farms(args: argparse.Namespace) -> int:
    from res_reliability import farms

    got = farms.fetch_all(data_dir=args.data_dir, budget=args.budget)
    logging.getLogger(__name__).info("ERA5 cached for %d sites", len(got))
    return 0


def _run(args: argparse.Namespace) -> int:
    from res_reliability import analyses, study

    log = logging.getLogger(__name__)
    stages = args.stages or STAGES
    ctx = None
    for st in stages:
        t0 = time.time()
        if st == "structural":
            analyses.structural_stage(Path(args.kelmarsh_dir) if args.kelmarsh_dir else None)
        else:
            if ctx is None:
                ctx = study.build(args.data_dir, with_kse=not args.no_kse)
            loads = [k for k in ("RTS", "KSE") if k in ctx.loads]
            if st == "load":
                if "KSE" in ctx.loads:
                    analyses.load_stage(ctx, args.data_dir)
            elif st == "markov":
                analyses.markov_stage(ctx)
            elif st == "variants":
                analyses.variants_stage(ctx, args.n_sims, loads)
            elif st == "elcc":
                analyses.elcc_stage(ctx, loads[::-1])
            elif st == "sensitivity":
                analyses.sensitivity_stage(ctx, "KSE" if "KSE" in loads else "RTS", max(args.n_sims // 4, 100))
            elif st == "weather":
                analyses.weather_stage(ctx)
            elif st == "storage":
                analyses.storage_stage(ctx, args.n_sims, args.jobs, args.data_dir)
            elif st == "verify":
                analyses.verify_stage(ctx, args.n_sims, args.jobs, args.data_dir)
            elif st == "theory":
                from res_reliability import operating
                operating.theory_stage(ctx, args.n_sims, args.data_dir)
            elif st == "failures":
                from res_reliability import operating
                operating.failures_stage(ctx, args.n_sims, args.data_dir)
            elif st == "long":
                analyses.long_stage(ctx, args.long_dir, list(range(args.long_years[0], args.long_years[1] + 1)))
        log.info("stage %s done in %.0f s", st, time.time() - t0)
    try:
        analyses.summary()
    except FileNotFoundError as exc:
        log.info("summary.json not written yet (%s)", exc.filename)
    return 0


def _figures(args: argparse.Namespace) -> int:
    from res_reliability import figures

    figures.make_all()
    return 0


def _dashboard(args: argparse.Namespace) -> int:
    from res_reliability import dashboard

    dashboard.build()
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="res_reliability", description=__doc__.splitlines()[0])
    p.add_argument("-v", "--verbose", action="store_true")
    p.add_argument("--data-dir", type=Path, default=config.DATA_DIR)
    sub = p.add_subparsers(dest="cmd", required=True)

    s = sub.add_parser("fetch-weather", help="download hourly ERA5 for the grid (Open-Meteo)")
    s.add_argument("--budget", type=float, default=2400.0, help="max Open-Meteo call weight per 24 h")
    s.add_argument("--pause", type=float, default=15.0, help="seconds between requests")
    s.add_argument("--start-year", type=int, default=config.YEARS[0])
    s.add_argument("--end-year", type=int, default=config.YEARS[-1])
    s.set_defaults(func=_fetch_weather)

    s = sub.add_parser("fetch-load", help="download hourly KSE demand from PSE")
    s.add_argument("--start-year", type=int, default=config.YEARS[0])
    s.add_argument("--end-year", type=int, default=config.YEARS[-1])
    s.set_defaults(func=_fetch_load)

    s = sub.add_parser("fetch-farms", help="download ERA5 at the open wind farms and the PV site (Open-Meteo)")
    s.add_argument("--budget", type=float, default=2400.0, help="max Open-Meteo call weight per 24 h")
    s.set_defaults(func=_fetch_farms)

    s = sub.add_parser("run", help="run study stages and write results/")
    s.add_argument("--stages", nargs="+", choices=STAGES + ["long"], help="default: all except long")
    s.add_argument("--n-sims", type=int, default=2000, help="Monte Carlo years per variant")
    s.add_argument("--no-kse", action="store_true", help="skip the KSE load profile (RTS-79 only)")
    s.add_argument("--kelmarsh-dir", help="Kelmarsh data folder (default: KELMARSH_DATA_DIR)")
    s.add_argument("--jobs", type=int, default=1, help="worker processes for the storage and verify stages")
    s.add_argument("--long-dir", type=Path, default=config.DATA_DIR / "long",
                   help="data folder of the earlier weather years (stage long)")
    s.add_argument("--long-years", type=int, nargs=2, default=[1980, 2014], metavar=("FIRST", "LAST"))
    s.set_defaults(func=_run)

    s = sub.add_parser("figures", help="render README charts from results/")
    s.set_defaults(func=_figures)

    s = sub.add_parser("dashboard", help="build the static dashboard docs/index.html")
    s.set_defaults(func=_dashboard)
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
    return args.func(args)
