# Reliability of wind and PV in Poland: Markov models, Monte Carlo adequacy and ELCC

[![tests](https://github.com/kwwaleszczynski-del/pl-res-reliability-markov/actions/workflows/ci.yml/badge.svg)](https://github.com/kwwaleszczynski-del/pl-res-reliability-markov/actions/workflows/ci.yml)

How much firm capacity is a wind or PV plant in Poland worth, and how much does the answer depend on how its
output is modelled? This project builds **multi-state Markov models of primary-energy availability** for wind and
PV at 29 ERA5 grid points over Poland (2015-2024), adds **structural reliability** estimated from open event logs,
runs **sequential Monte Carlo** adequacy assessment (LOLE, EENS, LOLF) on the IEEE Reliability Test System with
its own load and with the real hourly demand of the Polish power system (PSE), and computes the **effective load
carrying capability (ELCC)** of each point and of a distributed national portfolio, compared with the derating
factors (KWD) of the Polish capacity market.

It reproduces, on open data only, the method of the author's master's thesis *"Verification of reliability models
in photovoltaic and wind power systems"* (2025), and extends and refines it.

**Interactive dashboard:** <https://kwwaleszczynski-del.github.io/pl-res-reliability-markov/>

## Key results

National portfolio = equal capacity at the 29 grid points; weather and demand 2015-2024; generating system =
IEEE RTS-79 (3405 MW). Unless stated otherwise RES is 10 % of the peak load per technology, and Monte Carlo uses
5000 simulated years per variant (± = 95 % confidence interval of the mean).

**Test system check.** The analytical engine reproduces the published RTS-79 indices exactly (LOLE 9.394 h/yr,
LOEE 1176 MWh/yr); sequential Monte Carlo gives 9.22 ± 0.45 h/yr and 1153 ± 81 MWh/yr with 1.87 loss-of-load
events per year. With the real KSE load shape the same fleet meets the 3 h/yr standard at a peak of 2505 MW;
scaled to 2850 MW, as in the thesis, it would have LOLE 35 h/yr.

| Variant (wind 10 % + PV 10 % of peak) | KSE load: LOLE h/yr | EENS MWh/yr | events/yr | ELCC % | RTS-79 load: LOLE h/yr | ELCC % |
|---|---|---|---|---|---|---|
| conventional only | 3.06 ± 0.24 | 301 | 0.63 | | 9.22 ± 0.45 | |
| Markov monthly, FOR from literature (wind 5 %, PV 3 %) | 1.14 ± 0.13 | 102 | 0.28 | 20.7 | 3.91 ± 0.27 | 20.0 |
| Markov monthly, FOR from open data | 1.11 ± 0.13 | 99 | 0.27 | 21.2 | 3.82 ± 0.26 | 20.5 |
| Markov monthly, no structural outages | 1.09 ± 0.12 | 98 | 0.27 | 21.6 | 3.76 ± 0.26 | 20.8 |
| Markov homogeneous (one matrix), FOR from data | 1.24 ± 0.13 | 112 | 0.31 | 19.1 | 4.26 ± 0.28 | 18.3 |
| **historical weather years, FOR from data** | **1.27 ± 0.14** | **115** | **0.29** | **18.6** | **3.94 ± 0.27** | **19.8** |
| thesis implementation¹ | 1.32 ± 0.14 | 121 | 0.33 | 18.0 | 4.48 ± 0.29 | 17.1 |
| wind only (Markov monthly, FOR data) | 1.45 ± 0.15 | 131 | 0.33 | 31.3 | 4.90 ± 0.30 | 29.3 |
| PV only (Markov monthly, FOR data) | 2.43 ± 0.20 | 235 | 0.55 | 10.2 | 7.47 ± 0.40 | 10.1 |

¹ homogeneous chain, bin-centre state levels, every farm out with probability FOR (5 % / 3 %) independently in
each hour. ELCC from the Monte Carlo margins with common random numbers; the exact COPT values agree within
0.7 percentage points, and every Monte Carlo LOLE lies within its confidence interval of the exact expectation.

- **Model structure matters more than the reliability parameters.** Literature FOR (5 %) against the open-data
  value (0.42 %) moves the portfolio ELCC by 0.5 points and LOLE by 3 %. Replacing the monthly chain by a single
  homogeneous matrix costs 2.1 points: one matrix spreads the winter wind maximum over the whole year, which for
  wind alone means 25.4 % instead of 30.8 %.
- **Real weather-demand correlation costs wind ELCC.** With the KSE load and weather of the same year the wind
  portfolio's ELCC is 26.5 %; pairing the same loads with the *other* nine weather years gives 29.1 % (28.6-30.0 %
  for every shift). Cold, calm winter spells coincide with high demand. A Markov chain is independent of the load
  and cannot see this, so the monthly chain overstates wind ELCC (30.8 %) and only historical weather is a safe
  reference. PV shows no such effect (10.2 % in both pairings). On the synthetic RTS-79 load the chain and
  historical weather agree within about a point.
- **The number of states matters only through the state levels.** With the thesis's bin-centre levels (0, 25,
  50, 75, 100 %), a two-state model gives wind ELCC 7.6 % and 47 % more LOLE than 15 states, and the results settle
  only from about 7 states. With energy-preserving (conditional-mean) levels every model from 2 to 15 states is
  within ±1 point of ELCC (25.1-26.6 % homogeneous, 30.4-32.0 % monthly) and 4 % of LOLE. Event frequency
  (0.33-0.40 /yr) and duration (4.2-4.8 h) hardly move with the number of states either.
- **Where wind is worth most.** Wind ELCC per point (a plant of 10 % of peak, historical weather, KSE load)
  ranges 21.7-29.4 %, highest on the West Pomeranian coast, and follows the capacity factor (r = 0.94). PV ranges
  only 9.1-10.9 %, highest in the south and south-east.
- **Capacity market.** ELCC falls with penetration: onshore wind from 29.6 % (2 % of peak) to 19.1 % at 40 % of
  peak and 13.0 % at 100 %; PV from 12.4 % to 3.9 % at 80 % and 3.3 % at 100 % (historical weather, KSE load).
  The capacity-market factors for 2028-2030 are 12.5-13.8 % for onshore wind and 1.6-2.8 % for solar. With today's
  fleet (about 11 GW of onshore wind and over 20 GW of PV against a peak of about 28 GW) the model lands in the
  same order of magnitude, above the official factors. The official values come from the real system's stress
  hours, which this stylised fleet does not reproduce, so the comparison is indicative.
- **Structural reliability from open data.** Kelmarsh: 332 forced outages in 300 528 turbine-hours, FOR 0.44 %
  (98 stops ≥ 1 h: FOR 0.42 %, MTTR 12.9 h, MTTF 3054 h). That is about a tenth of the literature value of 5 %.
  Grid-side stops (0.74 % of the time), scheduled maintenance (0.69 %) and requested shutdowns (0.40 %) are
  separate categories. PVDAQ 34: one inverter failure of 34 days (389 daylight hours) in two years, FOR 5.7 % of
  daylight time against 3 % in the literature, from a single event. Raising wind FOR from 0 to 10 % lowers wind
  ELCC from 30.8 % to 28.0 % with independent turbines, or 26.7 % if the whole portfolio fails as one plant.
- **How much data.** A monthly chain fitted on one year gives wind ELCC between 27.4 % and 33.5 % depending on
  the year; five-year windows narrow this to 30.6-31.5 %. Single historical years range from 23.3 % (2023) to
  31.7 % (2020) for wind and 8.4-12.3 % for PV. Wind speeds ±5 % move wind ELCC by about ±2.4 points.

![ELCC map](docs/img/elcc_map.png)

![ELCC vs penetration](docs/img/elcc_curves.png)

*ELCC of the national portfolio versus installed capacity; bands are the capacity-market factors (KWD) for onshore
wind (blue) and solar (orange) in the 2028-2030 auctions.*

![Variants](docs/img/variants_lole.png)

### Markov models against the ERA5 series

![Transition matrices](docs/img/transition_matrices.png)

- **Energy and seasonality.** Mean wind CF is 28.7 % in the data, 28.6 % (homogeneous) and 28.9 % (monthly) in the
  chains. Only the monthly chain reproduces the swing between 39 % in winter and 18 % in August. PV: 14.3 % in the data,
  14.4 % and 14.8 % in the chains (the monthly PV chain overstates June-July by up to 2.4 points).
- **Persistence.** Mean residence times match by construction (12 h in the lowest wind state). The autocorrelation
  decays faster than in the data: 0.31 (homogeneous) and 0.40 (monthly) against 0.52 at 24 h. Low-output spells
  of 48 h or more occur 4.5 times a year in the data but 2.8-2.9 times in the chains, and the longest calm lasts
  165 h against 92-102 h. A first-order chain is good enough for expected LOLE/EENS, but not for studies where
  energy matters, such as storage sizing.
- **Out of sample.** Fitted on 2015-2022 and tested on 2023-2024 at all 29 points, the chains predict the test
  transitions as well as matrices fitted on the test years themselves. The median log-likelihood per transition
  for wind is -0.474 for both, -0.470 for the monthly chain, and -1.31 without persistence (i.i.d.). The
  Anderson-Goodman test still detects that the transition probabilities differ between the periods at 16 of 29
  points for wind and 13 of 29 for PV (p < 0.05). With up to 17 500 test transitions it detects small changes:
  the median largest difference in a single probability is 0.012 for wind and 0.029 for PV.

![Markov validation](docs/img/markov_validation.png)

![Sensitivity](docs/img/sensitivity.png)

![Structural reliability](docs/img/structural.png)

### From the thesis implementation to this repository

| Thesis implementation | This repository |
|---|---|
| simplified RTS-79 load built from daily, weekday and monthly factors (peak above 2850 MW) | published RTS-79 load model, reference LOLE and LOEE reproduced exactly |
| conventional units drawn independently every hour from FOR | two-state units with exponential MTTF / MTTR (sequential) |
| random numbers drawn in parallel workers without a fixed seed per year | `SeedSequence(seed, spawn_key=(year,))`: every simulated year reproducible on its own |
| one homogeneous transition matrix, bin-centre levels | also month-dependent chains, energy-preserving levels, and historical weather as the reference |
| historical weather and Markov simulation discussed together | Markov chain vs historical weather years compared explicitly |
| ELCC: removing whole conventional units, linear interpolation, fresh random numbers at each step | load shift solved exactly on the COPT, checked with common-random-number Monte Carlo |
| chi-square / KS tests on state frequencies | out-of-sample likelihood and the Anderson-Goodman test on transition counts |
| two farms and a PV aggregate | 29-point grid, map of ELCC and CF, distributed portfolio, comparison with KWD |
| synthetic multi-year KSE profile built from annual reports | real hourly PSE demand 2015-2024 aligned with the weather years |

## Method

### Pipeline

```
ERA5 (Open-Meteo), 29 points, 2015-2024 --> wind CF (MM92 power curve) / PV output + clear-sky envelope (pvlib)
                                                 |
                                    discretise into N states --> transition matrix P, stationary pi
                                    (homogeneous, or one matrix per month)
                                                 |
Kelmarsh status logs, PVDAQ --> FOR, MTTF, MTTR -+
                                                 |
IEEE RTS-79 units (MTTF/MTTR) + load ------------+--> sequential Monte Carlo (one seed per simulated year)
PSE KSE demand 2015-2024 (scaled) ---------------+        LOLE, EENS, LOLF, event duration
                                                 +--> exact expectation from the COPT --> ELCC (LOLE equality)
```

### 1. Multi-state Markov model of primary-energy availability

Each grid point gets an hourly wind capacity factor (ERA5 100 m wind speed, the measured IEC power curve of the
2.05 MW Senvion MM92 from the open Kelmarsh data, smoothed for a farm, 10 % wake and electrical losses) and an
hourly PV output (fixed tilt 35°, south, DC/AC 1.2: Erbs decomposition, isotropic transposition, Faiman cell
temperature). As in the thesis, the normalised output $`x_t`$ is split into $`N`$ states with thresholds halfway
between the levels $`0, \tfrac{1}{N-1}, \dots, 1`$ (for $`N=5`$: 12.5, 37.5, 62.5 and 87.5 %), and the hourly
transition probabilities are the maximum-likelihood estimate

```math
p_{ij} = \frac{n_{ij}}{\sum_k n_{ik}}, \qquad \boldsymbol{\pi} = \boldsymbol{\pi}\mathbf{P}, \quad \sum_i \pi_i = 1,
\qquad \lambda_{ij} = \frac{n_{ij}}{T_i}
```

where $`n_{ij}`$ counts transitions and $`T_i`$ is the time spent in state $`i`$ (the continuous-time intensities
$`\lambda_{ij}`$ of the thesis follow directly). For PV the chain runs on the clear-sky index of PV output during
daylight hours and the output is $`P_t = \ell(s_t)\,P^{cs}_t`$, so the daily and seasonal envelope is
deterministic and only cloudiness is stochastic. Two extensions of the thesis model:

- **state levels:** the level of a state is the conditional mean of $`x_t`$ in the bin (energy-preserving), not
  the bin centre (0, 25, 50, 75, 100 %);
- **seasonality:** besides the single homogeneous matrix of the thesis, a month-dependent chain (12 matrices).

### 2. Structural reliability

Every RES plant consists of 100 identical units (turbines or inverters), each a two-state process with
exponential up and down times:

```math
\mathrm{FOR} = \frac{\mathrm{MTTR}}{\mathrm{MTTF} + \mathrm{MTTR}}, \qquad
\mathbf{\Lambda}_{struct} = \begin{bmatrix} -1/\mathrm{MTTF} & 1/\mathrm{MTTF} \\ 1/\mathrm{MTTR} & -1/\mathrm{MTTR} \end{bmatrix}
```

The thesis's simplification (a farm is out with probability FOR independently in every hour, no outage
duration) is kept as a variant. Parameters come from the literature (wind FOR 5 %, PV FOR 3 %, as in the thesis)
or from open event logs: the forced-outage stops of the six Kelmarsh turbines (2016-2021) and the inverter of
NREL PVDAQ system 34.

### 3. Sequential Monte Carlo

Every simulated year runs hour by hour: the 32 RTS-79 units follow their MTTF/MTTR, the wind and PV chains move
from hour to hour (or a historical weather year is replayed), RES units fail and are repaired, and

```math
\mathrm{LOLE} = \mathbb{E}\Big[\sum_h \mathbf{1}(G_h < L_h)\Big], \qquad
\mathrm{EENS} = \mathbb{E}\Big[\sum_h \max(L_h - G_h, 0)\Big]
```

together with the frequency (LOLF) and mean duration of loss-of-load events. Year $`i`$ draws all its random
numbers from `SeedSequence(seed, spawn_key=(i,))`, so every year is reproducible on its own, independent of
chunking and parallel execution, and all system variants share the same conventional-unit histories (common
random numbers). Without storage the *expected* LOLE and EENS depend only on the hour-by-hour distribution of
available capacity, so they are also computed exactly from the capacity outage probability table (COPT); the
Monte Carlo estimate has to agree with it within its confidence interval, and adds frequency, duration and the
spread between years.

### 4. Test systems and variants

- **IEEE RTS-79:** 32 units, 3405 MW, chronological load model with a 2850 MW peak (8736 h). The implementation
  reproduces the reference adequacy indices exactly: LOLE **9.394 h/yr**, LOEE **1176 MWh/yr**, daily-peak LOLE
  1.3689 d/yr (Billinton & Allan).
- **KSE:** the actual hourly demand of the Polish power system in each of the ten weather years (PSE), every year
  normalised to its own peak and scaled so that the RTS-79 generating system alone just meets the Polish
  reliability standard of LOLE = 3 h/yr. Load and weather of the same calendar year stay aligned, so the
  historical-weather variant keeps the real weather-demand correlation.

RES is the distributed national portfolio (equal capacity at the 29 points): 10 % of peak load of wind and 10 %
of PV. The variants follow the thesis (conventional only, literature parameters, parameters from data, with and
without structural outages), plus the model-structure comparison (homogeneous vs monthly chain vs historical
weather) and the thesis implementation itself.

### 5. ELCC

ELCC is the constant load increase that the system can carry after the plant is added, at unchanged LOLE
(Garver 1966; Keane et al. 2011):

```math
\mathrm{LOLE}\big(C + R,\; L + \mathrm{ELCC}\big) = \mathrm{LOLE}\big(C,\; L\big)
```

It is solved by root finding on the exact COPT expectation, and checked by bisection on Monte Carlo margins with
common random numbers. The thesis instead removed whole conventional units and interpolated between the discrete
steps, with fresh random numbers at every step.

### 6. Validation of the Markov models

State residence times, autocorrelation, monthly means, long low-output spells, and cross-validation: fit on
2015-2022 and test on 2023-2024 at every point, plus leave-one-year-out for the portfolio. The tests are the
out-of-sample log-likelihood per transition (against an i.i.d. model and against the matrix fitted on the test
years themselves) and the Anderson-Goodman likelihood-ratio test of equal transition probabilities. They complement
the chi-square and Kolmogorov-Smirnov tests on state frequencies used in the thesis with tests that account for
the dependence between consecutive hours.

## Limitations

- **A stylised system, not the real KSE.** The generating side is the IEEE RTS-79 fleet (no imports, planned
  maintenance, storage or demand response); only the *shape* of the load is Polish. ELCC values are a property of
  this system at a given penetration, indicative for Poland, not a forecast of KWD.
- **KWD are not ELCC.** The capacity-market derating factors are set every year by regulation from the TSO's own
  analysis of the real system. They are an external reference here, and they change from auction to auction
  (check the current regulation before using them).
- **First-order chains forget too fast.** Even the monthly chain underestimates multi-day lulls (see validation).
  Expected LOLE/EENS are not affected in a system without storage, but a study of storage, demand response or
  energy-limited units would need a higher-order or regime-switching model, or historical weather.
- **Independence.** The wind and PV chains are simulated independently of each other and of the load; the
  historical-weather variant keeps both correlations and is the reference.
- **Weather and power models.** ERA5 at 0.25° smooths local wind; one turbine type, a 10 % loss factor and one PV
  configuration are used at all points. The national portfolio is an equal-weight grid, not the real fleet.
- **Thin structural data.** Kelmarsh is one UK farm with six turbines, and its operator files a stop as "Forced
  outage" more narrowly than the thesis did (every stop of an hour or more while wind was available). PVDAQ system 34 is a
  single US inverter with one long failure in two years; its FOR is an illustration, not an estimate for Polish PV.

## Data sources and attribution

- **ERA5:** Hersbach, H. et al. (2020), The ERA5 global reanalysis, *Q. J. R. Meteorol. Soc.* 146, 1999-2049;
  Copernicus Climate Change Service. Accessed through the Open-Meteo historical weather API (Zippenfenig, P.,
  Open-Meteo.com Weather API, doi:10.5281/zenodo.7970649), CC BY 4.0.
- **Polish power system demand:** PSE S.A., hourly "Krajowe zapotrzebowanie na moc" from the legacy CSV export
  on pse.pl (2015 to 2024-06-13) and report `kse-load` of <https://api.raporty.pse.pl> (15-minute actual load,
  averaged to hours, from 2024-06-14, the first day on the new platform). Public data of the Polish TSO; no
  explicit open-data licence is stated. Cross-checked with the actual total load for Poland in 2023 from the
  ENTSO-E Transparency Platform.
- **IEEE RTS-79:** IEEE Subcommittee on the Application of Probability Methods (1979), IEEE Reliability Test
  System, *IEEE Trans. Power Apparatus and Systems* PAS-98(6), 2047-2054; reference indices from Billinton, R. &
  Allan, R. N. (1996), *Reliability Evaluation of Power Systems*, 2nd ed., Plenum.
- **Kelmarsh wind farm data:** C. Plumley and R. Takeuchi, Cubico Sustainable Investments,
  doi:10.5281/zenodo.16807551, CC BY 4.0 (status logs, and the MM92 power curve derived from the 10-minute SCADA).
- **PVDAQ:** NREL, Photovoltaic Data Acquisition (PVDAQ) Public Datasets, Open Energy Data Initiative,
  doi:10.25984/1846021, CC BY 4.0 (system 34).
- **Capacity market:** regulations of the minister responsible for energy on auction parameters: Dz.U. 2024
  poz. 1154 (main auction, delivery 2029), Dz.U. 2025 poz. 1066 (main auction, 2030), Dz.U. 2026 poz. 1359
  (additional auctions, 2028).
- **Country outline:** Natural Earth, public domain.
- **Method:** Marchel, P., *Modele niezawodnościowe jednostek wytwórczych z ograniczeniami zdolności
  generacyjnej* (Reliability models of generating units with limited generating capability), doctoral
  dissertation, Warsaw University of Technology, supervisor J. Paska, defended 2022 (abstract and reviews in the
  university's public information bulletin, full text in the WUT Main Library); Garver, L. L. (1966), Effective
  load carrying capability of generating units, *IEEE Trans. PAS* 85(8); Keane, A. et al. (2011), Capacity value
  of wind power, *IEEE Trans. Power Systems* 26(2); Anderson, T. W. & Goodman, L. A. (1957), Statistical inference
  about Markov chains, *Ann. Math. Statist.* 28(1); Billinton, R., Chen, H. & Ghajar, R. (1996), A sequential
  simulation technique for adequacy evaluation of generating systems including wind energy, *IEEE Trans. Energy
  Conversion* 11(4); pvlib python (Anderson, K. S. et al. 2023, *JOSS* 8(92), 5994).

Code: MIT. Only derived results (the tables in `results/` and the charts) are in the repository; raw data are
downloaded into `data/` on the first run and are not redistributed.

## Running it

```bash
python -m venv .venv
.venv/Scripts/python -m pip install -e ".[test]"      # Linux/macOS: .venv/bin/python
python -m pytest                                      # synthetic data, no network, ~15 s
python -m res_reliability fetch-weather               # ERA5, one request per point (~2300 Open-Meteo call units)
python -m res_reliability fetch-load                  # PSE demand 2015-2024, one month per request
python -m res_reliability run                         # all stages -> results/ (~20 min; --n-sims to change)
python -m res_reliability figures                     # docs/img/*.png
python -m res_reliability dashboard                   # docs/index.html
```

The structural stage reads the Kelmarsh SCADA zips from `../wind-power-curve-kelmarsh/data` (or
`KELMARSH_DATA_DIR`) and PVDAQ system 34 from `../pv-anomaly-autoencoder/data/pvdaq/system_34`; both sibling
projects download them.

```
src/res_reliability/
  weather.py     Open-Meteo ERA5 download with a call-weight ledger
  pse.py         KSE demand: legacy PSE export + api.raporty.pse.pl, DST-safe UTC index
  rts79.py       IEEE RTS-79 units and load model
  power.py       wind power curve, PV model and clear-sky envelope
  markov.py      discretisation, transition matrices, stationary distribution, simulation
  structural.py  FOR / MTTF / MTTR from Kelmarsh and PVDAQ event data
  adequacy.py    capacity outage probability table, exact LOLE/EENS
  montecarlo.py  sequential Monte Carlo with a seed per simulated year
  elcc.py        ELCC by LOLE equality (analytical and Monte Carlo)
  validation.py  residence times, autocorrelation, cross-validation, Anderson-Goodman test
  analyses.py    study stages -> results/
  figures.py, dashboard.py
```
