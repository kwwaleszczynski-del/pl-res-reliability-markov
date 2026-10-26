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
in photovoltaic and wind power systems"* (2025), and extends and refines it. A second part adds what first-order
chains cannot do: a **regime-switching weather model** with multi-day lulls and the cold-calm link to demand, and
**storage** in the sequential Monte Carlo, with the ELCC of 1-8 h storage against the capacity-market factors and the
storage energy needed to keep LOLE unchanged as wind and PV grow ([Storage and multi-day lulls](#storage-and-multi-day-lulls)).

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
  energy matters, such as storage sizing (see [Storage and multi-day lulls](#storage-and-multi-day-lulls)).
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
| first-order chains only | regime-switching model with multi-day memory, conditioned on temperature; storage in the Monte Carlo |

## Storage and multi-day lulls

The first-order chains above are good enough for adequacy without storage, but they forget a calm within a day or
two, and the monthly chain overstates wind ELCC because it cannot see cold, calm days. Storage turns both into an
energy question. This part adds a **regime-switching weather model** with explicit lull durations, shared by wind,
PV and all grid points and optionally conditioned on temperature (method 7), and **storage inside the sequential
Monte Carlo** (method 8). The same batteries in the market: [pl-bess-revenue-stack](https://github.com/kwwaleszczynski-del/pl-bess-revenue-stack)
(arbitrage and balancing reserves), [pl-bess-saturation](https://github.com/kwwaleszczynski-del/pl-bess-saturation)
(spreads as the fleet grows) and [bess-capacity-augmentation](https://github.com/kwwaleszczynski-del/bess-capacity-augmentation)
(ageing against the capacity-market obligation).

### Weather models with memory

200 synthetic years per model against the ten ERA5 years of the national portfolio. "Regimes + temperature" is
the regime model conditioned on the temperature of the weather year; the lull thresholds are absolute and the
model is not fitted to them.

| National portfolio | ERA5 | Markov, monthly | regime switching | regimes + temperature |
|---|---|---|---|---|
| wind lulls of 5+ days (daily CF < 15 %), per year | **4.2** | 1.06 | 3.60 | 3.62 |
| same, October-March only | **0.7** | 0.06 | 0.70 | 0.80 |
| wind lulls of 5+ days below 10 %, per year | **0.6** | 0.06 | 0.55 | 0.60 |
| wind + PV lulls of 3+ days below 10 % (equal capacity), per year | **2.2** | 0.55 | 2.36 | 2.45 |
| longest lull of the year (CF < 15 %), mean, days | **9.0** | 5.4 | 8.1 | 8.1 |
| lowest 7-day mean wind CF of the year, % | **7.7** | 11.3 | 8.9 | 8.9 |
| hourly spells below 12.5 % of 48 h or more, per year | **4.5** | 3.4 | 4.7 | 4.8 |
| wind autocorrelation at 24 / 72 / 168 h | **0.52 / 0.26 / 0.17** | 0.39 / 0.18 / 0.15 | 0.56 / 0.27 / 0.16 | 0.55 / 0.27 / 0.16 |
| winter days with CF < 15 %: coldest / warmest third, % | **28 / 3** | 8 / 8 | 17 / 17 | 25 / 10 |
| wind CF in the 100 highest-load hours (annual mean 28.7 %), % | **26.2** | 38.4 | 38.3 | 31.6 |
| wind ELCC without storage (10 % of peak, KSE load), % | **26.5** | 30.8 | 30.7 | 28.6 |
| wind + PV ELCC without storage (10 % + 10 %), % | **33.5** | 37.2 | 37.0 | 35.0 |

![Weather models](docs/img/weather_models.png)

- **Memory is restored.** The monthly chain produces a quarter of the 5-day lulls and almost no lull longer than a
  week (0.06 a year in winter against 0.7). The regime model follows the ERA5 distribution of lull lengths from 1
  to 10 days, the autocorrelation over a whole week and the frequency of dark, calm spells of wind and PV
  together. It still makes the deepest week of a year a little milder than the data (8.9 % against 7.7 %).
- **Cold calm is half restored.** Conditioning on temperature makes calm winter days two and a half times more
  common in the coldest third of days than in the warmest (25 % against 10 %; ERA5 28 % against 3 %; the
  first-order chain is flat at 8 %). Wind output in the 100 highest-load hours drops from 38 % to 32 % of
  capacity (ERA5 26 %) and wind ELCC from 30.7 % to 28.6 %, halfway to the 26.5 % of historical weather. Without
  the covariate, memory alone changes nothing for ELCC without storage: only the hour-by-hour distribution and its
  alignment with the load matter there.
- **One regime for all points carries most of the spatial correlation.** Simulating the 29 grid points
  separately with a shared daily regime gives a mean correlation of daily wind anomalies between points of 0.63
  (ERA5 0.77, independent chains 0.01), and the portfolio rebuilt from the points has 4.3 lulls of 5+ days a year
  (ERA5 4.2); independent chains average out to a portfolio without a single one.

### Capacity value of storage

ELCC of a storage fleet in % of its power (KSE load scaled so that the system with its wind and PV is at LOLE
3 h/yr; round trip 85 %, half full on 1 January; 2000 simulated years with common random numbers; historical
weather, wind and PV 10 % of peak each):

| Storage fleet (power, % of peak) | 1 h | 2 h | 4 h | 8 h | marginal ELCC of the last block, 4 h |
|---|---|---|---|---|---|
| 1 % | 49.6 | 66.0 | 80.9 | 92.2 | 80.9 |
| 5 % | 35.1 | 53.5 | 73.1 | 89.6 | 71.1 |
| 10 % | 28.4 | 42.9 | 62.9 | 84.0 | 52.8 |
| 20 % | 21.7 | 33.0 | 49.1 | 64.3 | 35.2 |

![Storage ELCC](docs/img/storage_elcc.png)

- **Duration decides.** Loss-of-load events in this system last 4.6 h on average, and ELCC is set by the longest
  ones, so even the first small block is worth only half of its power at 1 h, two thirds at 2 h, 81 % at 4 h and
  92 % at 8 h.
- **The fleet saturates its own value.** Each block flattens the peak that the next one would have to cover:
  with a fleet of 20 % of the peak the average ELCC of 4 h storage is 49 % and the next block adds only 35 %.
- **Wind and PV help storage a little.** With 2 x 10, 25 and 50 % of the peak in wind and PV, 4 h storage at a 5 %
  fleet is worth 73, 76 and 79 % (68 % without any): PV removes the daytime from the stress hours, so the events
  left are shorter (4.6 h, 4.2 h and 3.9 h on average).
- **Co-located storage loses at scale.** Charging only from wind and PV, as behind a shared grid connection,
  costs little for a small fleet but a lot for a large one: 4 h at 20 % of the peak is worth 37.5 % instead of
  49.1 % (2 x 10 % RES), because a battery that cannot charge from the grid stays empty through a cold calm.
  Losses and the start half full cost 2 points on average (at most 5) against an ideal storage.
- **Against the capacity market.** The KWD of batteries fell from 95.00 % (main auction for 2028) to 61.30 %
  (2029) and 13.39 % (2030), one value for every duration. Following a growing fleet of storage with 2 x 25 % of
  the peak in wind and PV (historical weather), the marginal ELCC of the next block of 4 h storage is 85 % for
  the first one, falls to 61 % for a block added to a fleet of about 5-10 % of the peak (1.4-2.8 GW at the scale
  of the KSE) and to 13 % for one added at 30-40 % of the peak (8.5-11.3 GW). 95 % is only reached by the first
  block of 8 h storage (94.5 %). The 2030 value thus corresponds to a fleet of storage of the order of 10 GW in
  this stylised system, not to a short duration: from about 30 % of the peak the marginal value of 2, 4 and 8 h
  storage is the same (10-20 %), because on consecutive high-demand days there is too little surplus left to
  recharge, whatever the duration. This saturation comes from the load and the conventional fleet: the regime and
  monthly weather models give the same marginal values within 4 points.
- **The weather model hardly matters for 1-8 h storage.** The four models differ by 2 points in the median and 6
  at most (the monthly chain is the most optimistic at high RES, e.g. 77 % against 72 % for 4 h at a 10 % fleet
  with 2 x 50 % RES).
  Storage of a few hours is limited by the shape of the evening peak and cycles every day; the multi-day lulls
  that separate the models matter for energy, below.

### Storage for the lulls

How much storage keeps LOLE at today's level when demand grows and wind and PV supply exactly the extra energy
(method 8)? The conventional fleet stays, so without storage LOLE climbs steeply: the added wind and PV
(468-2808 MW each in the model) cover the energy, not the cold evening peaks. Storage energy is the smallest that
brings LOLE back to 3.0 h/yr, with power not limited and a round trip of 85 %, in GWh scaled to the 2024 KSE peak
(28.3 GW):

| Wind + PV share of annual demand | LOLE without storage, h/yr | historical weather | regimes + temperature | regime switching | Markov, monthly | hours of peak load (historical) |
|---|---|---|---|---|---|---|
| 10 % | 7.1 | **4.6** | 3.7 | 2.8 | 2.7 | 0.16 |
| 20 % | 24 | **24.5** | 20.5 | 17.5 | 16.3 | 0.87 |
| 30 % | 101 | **117** | 101 | 83 | 65 | 4.2 |
| 40 % | 370 | **441** | 414 | 366 | 252 | 15.6 |

![Storage for lulls](docs/img/storage_lulls.png)

- **The answer grows far faster than the share.** From a sixth of an hour of peak load at 10 % to four hours at
  30 % and more than half a day at 40 %. Up to about 20 % the storage has to cover a few evening hours and power
  binds: with a 4 h duration it needs 6.2 and 24.7 GWh, with 24 h 28 and 83 GWh. From 30 % the cold, calm weeks
  set the energy and any duration of 4 h or more needs the same.
- **The weather model changes the answer by up to 45 %.** The first-order monthly chain asks for 34-45 % less
  storage than historical weather (65 against 117 GWh at 30 %). Memory alone (regime switching) closes part of
  the gap at high shares (17 % short at 40 %) but none at 10 %, where the coincidence of cold, calm and high
  demand decides; regimes + temperature come within 6-20 %. A model that cannot produce long or cold lulls makes
  storage look cheaper than it is.
- **Losses matter more as storage gets larger.** A round trip of 85 % instead of 100 % adds 7 % at 10 % and 13-15 %
  at 30-40 %, because energy for long lulls is charged long before it is needed.
- **Scale.** These are energies for a stylised 2.5 GW system scaled up linearly, with no imports, demand response,
  pumped hydro or new firm capacity; in the real KSE all of them reduce the need. The robust result is the
  ratio between the weather models and the shape of the curve, not the GWh.

## Method

### Pipeline

```
ERA5 (Open-Meteo), 29 points, 2015-2024 --> wind CF (MM92 power curve) / PV output + clear-sky envelope (pvlib)
                                                 |
                                    discretise into N states --> transition matrix P, stationary pi
                                    (homogeneous, or one matrix per month)
ERA5 2 m temperature --> daily anomaly terciles -+-> daily regimes with explicit durations (semi-Markov),
                                                 |   hourly chains per (month, regime), shared by wind and PV
Kelmarsh status logs, PVDAQ --> FOR, MTTF, MTTR -+
                                                 |
IEEE RTS-79 units (MTTF/MTTR) + load ------------+--> sequential Monte Carlo (one seed per simulated year)
PSE KSE demand 2015-2024 (scaled) ---------------+        LOLE, EENS, LOLF, event duration
                                                 +--> exact expectation from the COPT --> ELCC (LOLE equality)
                                                 +--> storage on the hourly margins --> storage ELCC,
                                                      storage energy needed for a LOLE target
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

### 7. Regime-switching weather model

A slow daily layer is added on top of the hourly chains, fitted to the national portfolio:

- **Daily regimes.** Day $`d`$ gets a regime $`r_d \in \{0, 1, 2, 3\}`$ (calm, below normal, above normal, windy)
  from the ratio of its mean wind capacity factor to the smoothed climatology of that day of the year
  (31-day circular moving average over the ten years), with edges at the 20, 50 and 80 % quantiles of the ratio
  (0.49, 0.91 and 1.47).
- **Explicit durations.** Regime runs are a semi-Markov process. The probability that a run of regime $`r`$
  ends after its $`a`$-th day, when the next day is in season $`s`$ with temperature class $`c`$, is the
  discrete hazard
  ```math
  h_{r,s,c}(a) = P\big(r_{d+1} \neq r \,\big|\, r_d = r,\; a_d = a,\; s_{d+1} = s,\; c_{d+1} = c\big),
  \qquad P(D = a) = h(a) \prod_{b < a} \big(1 - h(b)\big)
  ```
  for ages 1-11 days and one open-ended bin from 12 days. The run lengths $`D`$ then follow their empirical
  distribution rather than the geometric distribution that any first-order chain implies. When a run ends, the
  next regime is drawn from the observed switches $`P(r' \mid r, s, c)`$.
- **Temperature covariate.** $`c_d`$ is the tercile (within the season DJF, MAM, JJA, SON) of the daily mean ERA5
  2 m temperature anomaly of the grid points. Simulated year $`i`$ uses the temperatures of weather year
  $`i \bmod 10`$, the same year as its load, so cold days, calm days and high demand can coincide. The model
  without the covariate ("regime switching") has $`c \equiv 0`$.
- **Shrinkage.** Thin cells borrow from pooled estimates: with $`n`$ exposures and $`e`$ ended runs a cell gets
  $`\hat h = (e + \kappa \hat h_{pooled}) / (n + \kappa)`$, $`\kappa = 10`$, pooling first over seasons and then over
  temperature classes (the same for the switch probabilities).
- **Hours.** Within the regime of the day, wind (5 capacity-factor states) and PV (5 clear-sky-index states,
  daylight hours) follow hourly chains with one matrix and one set of energy-preserving levels per (month,
  regime). The transition into an hour belongs to that hour's month and regime.
- **One regime path for everything.** Each simulated year draws one regime path from its own random stream
  (`SeedSequence(seed, spawn_key=(i, 999))`, independent of how many plants there are), shared by wind, PV and,
  in the spatial check, by all 29 grid points; every series then moves independently given the regimes. A
  30-day burn-in at the end of the previous weather year sets the residual run length on 1 January.

Validation against ERA5 uses quantities the model is not fitted to: lulls defined by an absolute threshold of the
daily capacity factor (10 % and 15 %, whole year and October-March), wind + PV lulls, hourly spells below
12.5 %, the autocorrelation over a week, the lowest 7-day mean of every year, the wind capacity factor in the 100
highest-load hours, the winter share of calm days by temperature tercile, wind ELCC without storage and the
correlation between grid points. Every model gets 200 synthetic years (20 per weather year).

### 8. Storage in the sequential Monte Carlo

A storage fleet of power $`P`$ and energy $`E`$ (duration $`E/P`$ = 1, 2, 4 or 8 h) works greedily on the
hourly margin $`m_t = G_t - L_t`$ of each simulated year: it discharges in a shortfall and charges from a
surplus,

```math
c_t = \min\!\Big(m_t^+,\; P,\; \frac{E - s_t}{\sqrt{\eta}}\Big), \qquad
d_t = \min\!\big(m_t^-,\; P,\; \sqrt{\eta}\, s_t\big), \qquad
s_{t+1} = s_t + \sqrt{\eta}\, c_t - \frac{d_t}{\sqrt{\eta}}
```

with $`0 \le s_t \le E`$ and round-trip efficiency $`\eta`$. It never charges during a shortfall, so it cannot cause
one. Three operating variants: **ideal** ($`\eta = 1`$, full on 1 January), **losses** ($`\eta = 0.85`$, half full
on 1 January; the default) and **co-located** (as losses, but charging only from the wind and PV output, like
storage behind one grid connection). The margins of a Monte Carlo run do not depend on the storage, so one run
with common random numbers serves every storage configuration. Hours in which no simulated year can be short and
every storage is full are skipped; the result is exactly the same as stepping through every hour (tested).

- **ELCC of storage** is found as for wind and PV, by bisection on a constant load added to every hour:
  ```math
  \mathrm{LOLE}\big(C + R + S,\; L + \mathrm{ELCC}_S\big) = \mathrm{LOLE}\big(C + R,\; L\big)
  ```
  for fleets of 1, 5, 10 and 20 % of the peak load, with 0, 2 x 10, 2 x 25 and 2 x 50 % of wind and PV. For every
  RES level the KSE load is scaled up so that LOLE with historical weather is back at 3 h/yr, so storage is
  always valued in a system at the reliability standard. The marginal ELCC is the increment over the next
  smaller fleet.
- **Storage for lulls.** Demand grows by $`g`$ and wind and PV (equal capacity) are sized to deliver exactly the
  extra energy, $`\bar R_g \cdot 8760 = g \sum_h L_h`$, so their share of demand is $`g / (1 + g)`$ = 10-40 %; the
  conventional fleet stays. The storage energy needed to keep LOLE at the level of today's system (conventional
  only, LOLE about 3 h/yr, same random numbers) is
  ```math
  E^* = \min\big\{E : \mathrm{LOLE}\big(C + R_g + S_E,\; (1 + g)\,L\big) \le \mathrm{LOLE}(C, L)\big\}
  ```
  found by bisection on $`\log E`$, with power not limited, or set by a duration of 4, 8 or 24 h. Here a 30-day
  warm-up (the end of the same simulated year) precedes 1 January, so a large storage does not start the year
  with free energy. $`E^*`$ is reported in hours of peak load and, scaled by the ratio of the 2024 KSE peak
  (28.3 GW) to the model peak (2505 MW), in GWh.

## Limitations

- **A stylised system, not the real KSE.** The generating side is the IEEE RTS-79 fleet (no imports, planned
  maintenance, pumped storage, demand response or reserve requirements); only the *shape* of the load is Polish.
  ELCC values are a property of this system at a given penetration, indicative for Poland, not a forecast of KWD.
  The same holds for storage: the GWh figures scale a 2.5 GW model to the 28 GW KSE and ignore interconnectors,
  flexible demand and hydro, all of which shrink the storage needed.
- **Storage operation.** The dispatch is greedy and reacts only to shortfalls: no price-driven arbitrage that
  could leave the storage empty before a peak, no reserve provision, no degradation and no forced outages of the
  storage itself. It is the best case for adequacy that a non-anticipating operator can achieve; real fleets that
  also trade energy and reserves (see the related projects) would earn less capacity credit.
- **KWD are not ELCC.** The capacity-market derating factors are set every year by regulation from the TSO's own
  analysis of the real system. They are an external reference here, and they change from auction to auction
  (check the current regulation before using them). For batteries one value covers every duration, and the
  drop from 95 % to 13.39 % reflects the TSO's view of the whole future fleet (including the large volume of
  storage already contracted), which this model only mimics through the fleet size.
- **First-order chains forget too fast.** Even the monthly chain underestimates multi-day lulls (see validation).
  Expected LOLE/EENS are not affected in a system without storage, but for storage the regime model or
  historical weather is needed.
- **The regime model is still an approximation.** It is fitted on the same ten years it is validated on (only
  quantities it is not fitted to are compared, but there is no out-of-sample test of lulls: ten years hold too
  few long ones). The temperature covariate recovers about half of the cold-calm coupling, not all of it (wind
  CF in the 100 highest-load hours and wind ELCC sit between the monthly chain and the data). Points share the
  daily regime but move independently within it, which leaves the inter-point correlation below the data.
- **Independence.** The monthly wind and PV chains are simulated independently of each other and of the load;
  the regime models share one daily regime between wind and PV and condition it on temperature; the
  historical-weather variant keeps every correlation and is the reference.
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
  (additional auctions, 2028; regulation of 9 October 2026, with the same factors as the main auction for 2028 in
  Dz.U. 2023 poz. 1561: onshore wind 12.62 %, solar 1.56 %, batteries 95.00 %). KWD of batteries, kinetic storage
  and supercapacitors, par. 10 pkt 9 of each:
  95.00 % in Dz.U. 2023 poz. 1561 (regulation of 4 August 2023, main auction, delivery 2028), 61.30 % in Dz.U.
  2024 poz. 1154 (26 July 2024, delivery 2029) and 13.39 % in Dz.U. 2025 poz. 1066 (1 August 2025, delivery 2030,
  one group with pumped-storage plants), read from the published texts on eli.gov.pl in October 2026; to be
  verified against the current regulation before any use.
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
python -m pytest                                      # synthetic data, no network, ~1 min
python -m res_reliability fetch-weather               # ERA5, one request per point (~2300 Open-Meteo call units)
python -m res_reliability fetch-load                  # PSE demand 2015-2024, one month per request
python -m res_reliability run                         # all stages -> results/ (~20 min + storage; --n-sims)
python -m res_reliability run --stages weather storage --jobs 8   # regime models and storage only
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
  regimes.py     daily regimes with explicit durations, temperature covariate, hourly chains per regime
  structural.py  FOR / MTTF / MTTR from Kelmarsh and PVDAQ event data
  adequacy.py    capacity outage probability table, exact LOLE/EENS
  montecarlo.py  sequential Monte Carlo with a seed per simulated year
  storage.py     greedy storage on hourly margins, storage ELCC, energy needed for a LOLE target
  elcc.py        ELCC by LOLE equality (analytical and Monte Carlo)
  validation.py  residence times, autocorrelation, cross-validation, Anderson-Goodman test, lulls
  analyses.py    study stages -> results/
  figures.py, dashboard.py
```
