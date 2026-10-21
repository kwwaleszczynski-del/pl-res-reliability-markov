"""Static interactive dashboard (Plotly) in ``docs/`` for GitHub Pages.

The page embeds the result tables as JSON and draws every chart client-side; plotly.min.js
is copied next to it, so the page has no external dependencies.
"""

from __future__ import annotations

import json
import logging
import shutil
from importlib import resources
from pathlib import Path

import pandas as pd

from res_reliability.config import DOCS_DIR, RESULTS_DIR, YEARS
from res_reliability.geo import poland_outline

log = logging.getLogger(__name__)


def _records(path: Path) -> list[dict]:
    return json.loads(pd.read_csv(path).round(4).to_json(orient="records"))


def collect(results_dir: Path = RESULTS_DIR) -> dict:
    j = lambda n: json.loads((results_dir / n).read_text(encoding="utf-8"))  # noqa: E731
    models = j("markov_models.json")
    port = {t: {k: models["portfolio"][t][k] for k in ("P_homogeneous", "P_monthly", "levels_mean",
                                                        "levels_monthly", "pi_homogeneous", "pi_monthly")}
            for t in ("wind", "pv")}
    val = {t: {k: models["validation"][t][k] for k in ("monthly_mean", "acf")} for t in ("wind", "pv")}
    st = j("structural.json")
    return {
        "years": [YEARS[0], YEARS[-1]],
        "outline": poland_outline().round(3).tolist(),
        "points": _records(results_dir / "elcc_points.csv"),
        "curves": _records(results_dir / "elcc_curves.csv"),
        "lole": _records(results_dir / "lole_vs_share.csv"),
        "variants": _records(results_dir / "variants.csv"),
        "matrices": j("markov_point_matrices.json"),
        "portfolio": port,
        "validation": val,
        "kwd": j("kwd.json")["values"],
        "structural": st["parameters"],
    }


def build(results_dir: Path = RESULTS_DIR, docs_dir: Path = DOCS_DIR) -> Path:
    data = collect(results_dir)
    docs_dir.mkdir(parents=True, exist_ok=True)
    js = resources.files("plotly").joinpath("package_data/plotly.min.js")
    with resources.as_file(js) as src:
        shutil.copyfile(src, docs_dir / "plotly.min.js")
    html = TEMPLATE.replace("__DATA__", json.dumps(data, separators=(",", ":")))
    out = docs_dir / "index.html"
    out.write_text(html, encoding="utf-8")
    (docs_dir / ".nojekyll").write_text("", encoding="utf-8")
    log.info("wrote %s (%.0f kB)", out, out.stat().st_size / 1024)
    return out


TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Wind and PV adequacy in Poland</title>
<meta name="description" content="Multi-state Markov models of wind and PV, sequential Monte Carlo LOLE/EENS and ELCC on open data (ERA5, PSE, IEEE RTS-79).">
<script src="plotly.min.js"></script>
<style>
:root {
  color-scheme: light;
  --page: #f9f9f7; --surface: #fcfcfb; --ink: #0b0b0b; --ink2: #52514e; --muted: #898781;
  --grid: #e1e0d9; --base: #c3c2b7; --border: rgba(11,11,11,0.10);
  --s1: #2a78d6; --s2: #eb6834; --s3: #1baf7a;
  --seq: #cde2fb,#9ec5f4,#6da7ec,#3987e5,#2a78d6,#1c5cab,#104281,#0d366b;
}
@media (prefers-color-scheme: dark) {
  :root:not([data-theme="light"]) {
    color-scheme: dark;
    --page: #0d0d0d; --surface: #1a1a19; --ink: #f0efec; --ink2: #c3c2b7; --muted: #898781;
    --grid: #2c2c2a; --base: #383835; --border: rgba(255,255,255,0.10);
    --s1: #3987e5; --s2: #d95926; --s3: #199e70;
    --seq: #0d366b,#104281,#184f95,#1c5cab,#256abf,#2a78d6,#5598e7,#9ec5f4;
  }
}
:root[data-theme="dark"] {
  color-scheme: dark;
  --page: #0d0d0d; --surface: #1a1a19; --ink: #f0efec; --ink2: #c3c2b7; --muted: #898781;
  --grid: #2c2c2a; --base: #383835; --border: rgba(255,255,255,0.10);
  --s1: #3987e5; --s2: #d95926; --s3: #199e70;
  --seq: #0d366b,#104281,#184f95,#1c5cab,#256abf,#2a78d6,#5598e7,#9ec5f4;
}
* { box-sizing: border-box; }
body { margin: 0; background: var(--page); color: var(--ink);
  font: 15px/1.5 system-ui, -apple-system, "Segoe UI", sans-serif; }
main { max-width: 1180px; margin: 0 auto; padding: 24px 16px 48px; }
header { display: flex; justify-content: space-between; align-items: flex-start; gap: 16px; }
h1 { font-size: 26px; margin: 0 0 6px; font-weight: 650; }
h2 { font-size: 17px; margin: 0 0 4px; font-weight: 620; }
p.lead { color: var(--ink2); margin: 0 0 18px; max-width: 860px; }
.note { color: var(--muted); font-size: 13px; margin: 4px 0 0; }
.grid { display: grid; grid-template-columns: repeat(12, 1fr); gap: 16px; }
.card { background: var(--surface); border: 1px solid var(--border); border-radius: 10px; padding: 16px; min-width: 0; }
.span12 { grid-column: span 12; } .span7 { grid-column: span 7; } .span5 { grid-column: span 5; } .span6 { grid-column: span 6; }
.span3 { grid-column: span 3; }
@media (max-width: 900px) { .span7, .span5, .span6, .span3 { grid-column: span 12; } }
.tile .label { color: var(--ink2); font-size: 13px; }
.tile .value { font-size: 30px; font-weight: 620; margin: 2px 0; }
.tile .sub { color: var(--muted); font-size: 12.5px; }
.controls { display: flex; flex-wrap: wrap; gap: 8px 14px; margin: 8px 0 4px; align-items: center; font-size: 13px; color: var(--ink2); }
select, button { font: inherit; font-size: 13px; color: var(--ink); background: var(--surface);
  border: 1px solid var(--base); border-radius: 6px; padding: 4px 8px; }
button { cursor: pointer; }
.chart { width: 100%; height: 380px; }
.chart.tall { height: 470px; }
table { border-collapse: collapse; width: 100%; font-size: 13px; }
th, td { text-align: right; padding: 6px 8px; border-bottom: 1px solid var(--grid); font-variant-numeric: tabular-nums; }
th:first-child, td:first-child { text-align: left; }
th { color: var(--ink2); font-weight: 600; }
.tablewrap { overflow-x: auto; }
a { color: var(--s1); }
footer { color: var(--muted); font-size: 12.5px; margin-top: 24px; }
</style>
</head>
<body>
<main>
<header>
  <div>
    <h1>Wind and PV adequacy in Poland</h1>
    <p class="lead">Multi-state Markov models of wind and PV availability fitted to ERA5 at 29 grid points
    (<span id="yrs"></span>), sequential Monte Carlo LOLE/EENS on the IEEE RTS-79 generating system with the
    RTS-79 and real KSE load (PSE), and ELCC by LOLE equality compared with the capacity-market derating
    factors (KWD).</p>
  </div>
  <button id="theme" title="Switch light/dark">Theme</button>
</header>

<section class="grid" id="tiles"></section>

<section class="grid" style="margin-top:16px">
  <div class="card span7">
    <h2>ELCC and capacity factor per grid point</h2>
    <div class="controls">
      <label>Metric <select id="m-metric">
        <option value="wind:elcc_pct">Wind ELCC</option><option value="wind:cf_pct">Wind capacity factor</option>
        <option value="pv:elcc_pct">PV ELCC</option><option value="pv:cf_pct">PV capacity factor</option></select></label>
      <label>Model <select id="m-model"><option value="hist">historical weather</option>
        <option value="seasonal">Markov, monthly</option></select></label>
      <label>Load <select id="m-load"><option value="KSE">KSE (PSE)</option><option value="RTS">IEEE RTS-79</option></select></label>
    </div>
    <div id="map" class="chart tall"></div>
    <p class="note">Plant of 10 % of peak load at one point. Click a point to show its transition matrix.</p>
  </div>
  <div class="card span5">
    <h2 id="hm-title">Transition matrix</h2>
    <div class="controls">
      <label>Source <select id="h-tech"><option value="wind">wind</option><option value="pv">PV</option></select></label>
      <label>Month <select id="h-month"></select></label>
      <button id="h-reset">National portfolio</button>
    </div>
    <div id="heat" class="chart tall"></div>
    <p class="note">P(state at h+1 | state at h). States are 5 equal-width bins of the capacity factor (wind) or of
    the clear-sky index of PV output (daylight hours), labelled by the mean level in each bin.</p>
  </div>
</section>

<section class="grid" style="margin-top:16px">
  <div class="card span6">
    <h2>LOLE as RES is added</h2>
    <div class="controls"><label>Load <select id="l-load"><option value="KSE">KSE (PSE)</option>
      <option value="RTS">IEEE RTS-79</option></select></label></div>
    <div id="lole" class="chart"></div>
    <p class="note">Expected LOLE of the RTS-79 generating system with the national portfolio added; solid =
    historical weather years, dashed = monthly Markov chains. Dotted line: Polish reliability standard 3 h/yr.</p>
  </div>
  <div class="card span6">
    <h2>ELCC versus penetration</h2>
    <div class="controls"><label>Load <select id="e-load"><option value="KSE">KSE (PSE)</option>
      <option value="RTS">IEEE RTS-79</option></select></label></div>
    <div id="elcc" class="chart"></div>
    <p class="note">Shaded bands: KWD for onshore wind and solar in recent auctions (regulations of the minister
    for energy, Dz.U. 2024 poz. 1154, 2025 poz. 1066, 2026 poz. 1359; they change every year).</p>
  </div>
</section>

<section class="grid" style="margin-top:16px">
  <div class="card span6">
    <h2>Markov model vs ERA5: seasonality</h2>
    <div class="controls"><label>Source <select id="v-tech"><option value="wind">wind</option>
      <option value="pv">PV</option></select></label></div>
    <div id="season" class="chart"></div>
  </div>
  <div class="card span6">
    <h2>Markov model vs ERA5: persistence</h2>
    <div id="acf" class="chart" style="margin-top:38px"></div>
  </div>
</section>

<section class="card" style="margin-top:16px">
  <h2>System variants (sequential Monte Carlo)</h2>
  <div class="controls"><label>Load <select id="t-load"><option value="KSE">KSE (PSE)</option>
    <option value="RTS">IEEE RTS-79</option></select></label></div>
  <div class="tablewrap"><table id="vtable"></table></div>
  <p class="note">10 % wind + 10 % PV of peak load. ± is the 95 % confidence interval of the Monte Carlo mean;
  "analytical" is the exact expectation from the capacity outage probability table.</p>
</section>

<footer>
  Data: ERA5 (Copernicus C3S, via Open-Meteo), PSE S.A. (KSE demand), IEEE RTS-79, Kelmarsh wind farm SCADA
  (Plumley 2022, CC BY 4.0), NREL PVDAQ (CC BY 4.0), Natural Earth. Method after the author's master's thesis and
  P. Marchel's doctoral dissertation (Warsaw University of Technology, 2021).
</footer>
</main>
<script>
const D = __DATA__;
const css = n => getComputedStyle(document.documentElement).getPropertyValue(n).trim();
const MONTHS = ['all year','Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'];
document.getElementById('yrs').textContent = D.years[0] + '-' + D.years[1];
let hmPoint = null;

function seq() { const s = css('--seq').split(',').map(x => x.trim()); return s.map((c, i) => [i / (s.length - 1), c]); }
function base(extra) {
  return Object.assign({
    paper_bgcolor: css('--surface'), plot_bgcolor: css('--surface'),
    font: {family: 'system-ui, -apple-system, "Segoe UI", sans-serif', size: 12, color: css('--ink2')},
    margin: {l: 52, r: 16, t: 10, b: 44},
    xaxis: {gridcolor: css('--grid'), linecolor: css('--base'), zerolinecolor: css('--base'), tickfont: {color: css('--muted')}},
    yaxis: {gridcolor: css('--grid'), linecolor: css('--base'), zerolinecolor: css('--base'), tickfont: {color: css('--muted')}},
    legend: {orientation: 'h', y: -0.2, font: {color: css('--ink2')}},
    hoverlabel: {bgcolor: css('--surface'), bordercolor: css('--base'), font: {color: css('--ink')}},
  }, extra || {});
}
const cfg = {displaylogo: false, responsive: true, modeBarButtonsToRemove: ['lasso2d', 'select2d']};
const fmt = (v, d = 1) => (v === null || v === undefined || Number.isNaN(v)) ? '-' : Number(v).toFixed(d);

function tiles() {
  const pick = (load, mix, model, share) => D.curves.find(r => r.load === load && r.mix === mix && r.model === model && Math.abs(r.share - share) < 1e-9);
  const w = pick('KSE', 'wind', 'hist', 0.1), p = pick('KSE', 'pv', 'hist', 0.1);
  const bau = D.variants.find(r => r.load === 'KSE' && r.variant === 'BAU');
  const dat = D.variants.find(r => r.load === 'KSE' && r.variant === 'DATA');
  const s = D.structural;
  const items = [
    ['Wind ELCC, national portfolio', fmt(w && w.elcc_pct) + ' %', 'KSE load, 10 % of peak, historical weather'],
    ['PV ELCC, national portfolio', fmt(p && p.elcc_pct) + ' %', 'same system; KWD solar 1.6-2.8 %'],
    ['LOLE with 10 % wind + 10 % PV', fmt(dat && dat.lole, 2) + ' h/yr', 'vs ' + fmt(bau && bau.lole, 2) + ' h/yr conventional only'],
    ['Wind FOR: literature vs data', fmt(100 * s.literature.wind.for, 1) + ' % / ' + fmt(100 * s.data.wind.for, 2) + ' %', 'Kelmarsh forced outages (≥ 1 h)'],
  ];
  document.getElementById('tiles').innerHTML = items.map(([l, v, s2]) =>
    `<div class="card span3 tile"><div class="label">${l}</div><div class="value">${v}</div><div class="sub">${s2}</div></div>`).join('');
}

function drawMap() {
  const [tech, col] = document.getElementById('m-metric').value.split(':');
  const model = document.getElementById('m-model').value, load = document.getElementById('m-load').value;
  const rows = D.points.filter(r => r.tech === tech && r.model === model && r.load === load);
  const label = (col === 'cf_pct' ? 'capacity factor' : 'ELCC') + ', %';
  const k = Math.cos(52 * Math.PI / 180);
  const traces = [
    {x: D.outline.map(p => p[0]), y: D.outline.map(p => p[1]), mode: 'lines', line: {color: css('--base'), width: 1.5},
     hoverinfo: 'skip', showlegend: false},
    {x: rows.map(r => r.lon), y: rows.map(r => r.lat), customdata: rows.map(r => r.site), mode: 'markers+text',
     text: rows.map(r => fmt(r[col])), textposition: 'bottom center', textfont: {size: 10, color: css('--ink2')},
     marker: {size: 18, color: rows.map(r => r[col]), colorscale: seq(), line: {color: css('--surface'), width: 2},
              colorbar: {title: {text: label, side: 'right'}, thickness: 12, outlinewidth: 0, tickfont: {color: css('--muted')}}},
     hovertemplate: '%{customdata}<br>' + label + ': %{marker.color:.1f}<extra></extra>', showlegend: false},
  ];
  const lay = base({margin: {l: 8, r: 8, t: 8, b: 8}});
  lay.xaxis = {visible: false, range: [13.9, 24.4], constrain: 'domain'};
  lay.yaxis = {visible: false, scaleanchor: 'x', scaleratio: 1 / k, constrain: 'domain'};
  Plotly.react('map', traces, lay, cfg);
  const el = document.getElementById('map');
  el.removeAllListeners && el.removeAllListeners('plotly_click');
  el.on('plotly_click', ev => { const pt = ev.points.find(q => q.customdata); if (pt) { hmPoint = pt.customdata; drawHeat(); } });
}

function drawHeat() {
  const tech = document.getElementById('h-tech').value;
  const msel = document.getElementById('h-month');
  let P, lv, title;
  if (hmPoint) {
    msel.value = '0'; msel.disabled = true;
    P = D.matrices[hmPoint][tech]; lv = D.matrices[hmPoint][tech + '_levels'];
    title = (tech === 'pv' ? 'PV' : 'Wind') + ' transition matrix, point ' + hmPoint;
  } else {
    msel.disabled = false;
    const m = +msel.value;
    P = m ? D.portfolio[tech].P_monthly[m - 1] : D.portfolio[tech].P_homogeneous;
    lv = m ? D.portfolio[tech].levels_monthly[m - 1] : D.portfolio[tech].levels_mean;
    title = (tech === 'pv' ? 'PV' : 'Wind') + ' transition matrix, national portfolio, ' + MONTHS[m];
  }
  document.getElementById('hm-title').textContent = title;
  // hide states that are never entered (e.g. CF above 87.5 % for the national wind portfolio)
  const keep = P.map((_, j) => P.some(r => r[j] > 0));
  const idx = keep.map((k, i) => k ? i : -1).filter(i => i >= 0);
  P = idx.map(i => idx.map(j => P[i][j]));
  const labels = idx.map(i => 'S' + i + ' ' + (100 * lv[i]).toFixed(0) + '%');
  const text = P.map(r => r.map(v => v.toFixed(2)));
  const trace = {z: P, x: labels, y: labels, type: 'heatmap', colorscale: seq(), zmin: 0, zmax: 1, xgap: 2, ygap: 2,
    hovertemplate: 'from %{y}<br>to %{x}<br>p = %{z:.3f}<extra></extra>', colorbar: {thickness: 12, outlinewidth: 0}};
  const lay = base({margin: {l: 80, r: 10, t: 10, b: 60}});
  // cell labels in ink or white depending on the fill (the ramp runs the other way in dark mode)
  const darkMode = css('--surface') === '#1a1a19';
  lay.annotations = P.flatMap((r, i) => r.map((v, j) => ({x: labels[j], y: labels[i], text: text[i][j], showarrow: false,
    font: {size: 11, color: v > 0.55 ? (darkMode ? '#0b0b0b' : '#ffffff') : (darkMode ? '#f0efec' : '#0b0b0b')}})));
  lay.xaxis = Object.assign(lay.xaxis, {title: {text: 'to state'}, side: 'bottom', showgrid: false});
  lay.yaxis = Object.assign(lay.yaxis, {title: {text: 'from state'}, autorange: 'reversed', showgrid: false});
  Plotly.react('heat', [trace], lay, cfg);
}

const MIX = {'wind': '--s1', 'pv': '--s2', 'wind+pv': '--s3'};
function drawLole() {
  const load = document.getElementById('l-load').value;
  const rows = D.lole.filter(r => r.load === load);
  const b = rows.find(r => r.mix === 'none').lole;
  const tr = [];
  for (const mix of ['wind', 'pv', 'wind+pv']) for (const model of ['hist', 'seasonal']) {
    const s = rows.filter(r => r.mix === mix && r.model === model);
    tr.push({x: [0, ...s.map(r => 100 * r.share)], y: [b, ...s.map(r => r.lole)], mode: 'lines',
      name: mix + (model === 'hist' ? '' : ' (Markov)'), line: {color: css(MIX[mix]), width: 2, dash: model === 'hist' ? 'solid' : 'dash'},
      hovertemplate: mix + ', ' + (model === 'hist' ? 'historical' : 'Markov monthly') + '<br>%{x:.0f} % of peak: LOLE %{y:.2f} h/yr<extra></extra>'});
  }
  const lay = base();
  lay.yaxis = Object.assign(lay.yaxis, {type: 'log', title: {text: 'LOLE, h/yr'}});
  lay.xaxis = Object.assign(lay.xaxis, {title: {text: 'RES capacity, % of peak load'}});
  lay.shapes = [{type: 'line', xref: 'paper', x0: 0, x1: 1, y0: 3, y1: 3, line: {color: css('--muted'), width: 1, dash: 'dot'}}];
  Plotly.react('lole', tr, lay, cfg);
}

function drawElcc() {
  const load = document.getElementById('e-load').value;
  const rows = D.curves.filter(r => r.load === load);
  const tr = [];
  for (const mix of ['wind', 'pv', 'wind+pv']) for (const model of ['hist', 'seasonal']) {
    const s = rows.filter(r => r.mix === mix && r.model === model);
    tr.push({x: s.map(r => 100 * r.share), y: s.map(r => r.elcc_pct), mode: 'lines+markers', marker: {size: 6},
      name: mix + (model === 'hist' ? '' : ' (Markov)'), line: {color: css(MIX[mix]), width: 2, dash: model === 'hist' ? 'solid' : 'dash'},
      hovertemplate: mix + ', ' + (model === 'hist' ? 'historical' : 'Markov monthly') + '<br>%{x:.0f} % of peak: ELCC %{y:.1f} %<extra></extra>'});
  }
  const lay = base();
  lay.yaxis = Object.assign(lay.yaxis, {title: {text: 'ELCC, % of capacity'}, rangemode: 'tozero'});
  lay.xaxis = Object.assign(lay.xaxis, {title: {text: 'installed RES capacity, % of peak load'}});
  if (load === 'KSE') {
    const w = D.kwd.map(k => k.wind), p = D.kwd.map(k => k.pv);
    lay.shapes = [
      {type: 'rect', xref: 'paper', x0: 0, x1: 1, y0: Math.min(...w), y1: Math.max(...w), fillcolor: css('--s1'), opacity: 0.12, line: {width: 0}},
      {type: 'rect', xref: 'paper', x0: 0, x1: 1, y0: Math.min(...p), y1: Math.max(...p), fillcolor: css('--s2'), opacity: 0.15, line: {width: 0}}];
    lay.annotations = [
      {xref: 'paper', x: 1, y: Math.max(...w), yanchor: 'bottom', xanchor: 'right', showarrow: false, text: 'KWD onshore wind', font: {size: 11, color: css('--ink2')}},
      {xref: 'paper', x: 1, y: Math.max(...p), yanchor: 'bottom', xanchor: 'right', showarrow: false, text: 'KWD solar', font: {size: 11, color: css('--ink2')}}];
  }
  Plotly.react('elcc', tr, lay, cfg);
}

const MODELC = {data: '--s1', homogeneous: '--s2', seasonal: '--s3'};
const MODELN = {data: 'ERA5 series', homogeneous: 'Markov, homogeneous', seasonal: 'Markov, monthly'};
function drawValidation() {
  const tech = document.getElementById('v-tech').value;
  const v = D.validation[tech];
  const tr = Object.keys(MODELC).map(k => ({x: MONTHS.slice(1), y: v.monthly_mean[k].map(x => 100 * x), mode: 'lines+markers',
    name: MODELN[k], line: {color: css(MODELC[k]), width: 2}, marker: {size: 7},
    hovertemplate: MODELN[k] + '<br>%{x}: %{y:.1f} %<extra></extra>'}));
  const lay = base();
  lay.yaxis = Object.assign(lay.yaxis, {title: {text: 'mean capacity factor, %'}, rangemode: 'tozero'});
  Plotly.react('season', tr, lay, cfg);
  const a = D.validation.wind.acf;
  const tr2 = Object.keys(MODELC).map(k => ({x: a[k].map((_, i) => i), y: a[k], mode: 'lines', name: MODELN[k],
    line: {color: css(MODELC[k]), width: 2}, hovertemplate: MODELN[k] + '<br>lag %{x} h: %{y:.3f}<extra></extra>'}));
  const lay2 = base();
  lay2.yaxis = Object.assign(lay2.yaxis, {title: {text: 'wind autocorrelation'}, range: [0, 1]});
  lay2.xaxis = Object.assign(lay2.xaxis, {title: {text: 'lag, h'}, dtick: 12});
  Plotly.react('acf', tr2, lay2, cfg);
}

function drawTable() {
  const load = document.getElementById('t-load').value;
  const rows = D.variants.filter(r => r.load === load);
  const head = '<tr><th>Variant</th><th>LOLE h/yr</th><th>analytical</th><th>EENS MWh/yr</th><th>LOLF occ/yr</th>' +
    '<th>mean duration h</th><th>ELCC % (MC)</th><th>ELCC % (analytical)</th></tr>';
  const body = rows.map(r => `<tr><td>${r.variant}: ${r.label}</td><td>${fmt(r.lole, 2)} ± ${fmt(r.lole_ci95, 2)}</td>` +
    `<td>${fmt(r.lole_analytic, 2)}</td><td>${fmt(r.eens, 0)}</td><td>${fmt(r.lolf, 2)}</td><td>${fmt(r.lold_h, 1)}</td>` +
    `<td>${fmt(r.elcc_mc_pct)}</td><td>${fmt(r.elcc_analytic_pct)}</td></tr>`).join('');
  document.getElementById('vtable').innerHTML = head + body;
}

function drawAll() { tiles(); drawMap(); drawHeat(); drawLole(); drawElcc(); drawValidation(); drawTable(); }
document.getElementById('h-month').innerHTML = MONTHS.map((m, i) => `<option value="${i}">${m}</option>`).join('');
for (const [id, fn] of [['m-metric', drawMap], ['m-model', drawMap], ['m-load', drawMap], ['h-tech', drawHeat],
  ['h-month', drawHeat], ['l-load', drawLole], ['e-load', drawElcc], ['v-tech', drawValidation], ['t-load', drawTable]])
  document.getElementById(id).addEventListener('change', fn);
document.getElementById('h-reset').addEventListener('click', () => { hmPoint = null; drawHeat(); });
document.getElementById('theme').addEventListener('click', () => {
  const dark = matchMedia('(prefers-color-scheme: dark)').matches;
  const cur = document.documentElement.dataset.theme || (dark ? 'dark' : 'light');
  document.documentElement.dataset.theme = cur === 'dark' ? 'light' : 'dark';
  drawAll();
});
matchMedia('(prefers-color-scheme: dark)').addEventListener('change', drawAll);
if (!D.variants.some(r => r.load === 'KSE')) for (const id of ['m-load', 'l-load', 'e-load', 't-load'])
  document.getElementById(id).value = 'RTS';
drawAll();
</script>
</body>
</html>
"""
