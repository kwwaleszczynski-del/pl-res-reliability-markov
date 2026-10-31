# Theory against operating data: method and details

Part of [pl-res-reliability-markov](../README.md) (method 10; the summary is in the README section
[Theory against operating data](../README.md#theory-against-operating-data)). All numbers come from
`python -m res_reliability run --stages theory` and are in `results/theory_*.csv`, `results/theory.json` and
`results/summary.json`.

The question of the thesis was whether a reliability model built from reanalysis weather and a catalogue
power curve is good enough for adequacy studies, or whether it has to be calibrated on operating data (the
thesis variants OZE-LIT and OZE-EMP). Its SCADA cannot be published, so the comparison is repeated here on the
open operating data of three wind farms and one PV system.

## 10.1 Plants and data

| Plant | Units | Hub, rotor, specific power | SCADA hours used | Event log |
|---|---|---|---|---|
| Kelmarsh (UK) | 6 x Senvion MM92, 2.05 MW | 68.5-78.5 m, 92 m, 308 W/m² | 49 276 (15.04.2016-2021) | Greenbyte status log 2016-2021 |
| Penmanshiel (UK) | 14 x Senvion MM82, 2.05 MW | 59 m, 82 m, 388 W/m² | 16 984 (2018 and 2020) | Greenbyte status log 2016-2024 |
| Hill of Towie (UK) | 21 x Siemens SWT-2.3-VS-82, 2.3 MW | 59 m, 82 m, 436 W/m² | 8 743 (2021) | controller alarm log 2016-04/2026 |
| PVDAQ 34 (Las Vegas, US) | 146.6 kWp, one 135 kW inverter, 11.2° south | - | 25 585 (2016-2018) | outages from the data |

The 10-minute SCADA is read from the caches of the sibling projects that download the Zenodo and OEDI archives
([wind-power-curve-kelmarsh](https://github.com/kwwaleszczynski-del/wind-power-curve-kelmarsh),
[wind-downtime-availability](https://github.com/kwwaleszczynski-del/wind-downtime-availability),
[wind-om-reliability](https://github.com/kwwaleszczynski-del/wind-om-reliability),
[pv-budget-bridge](https://github.com/kwwaleszczynski-del/pv-budget-bridge)). ERA5 at each site comes from the
Open-Meteo archive, one request per site (`python -m res_reliability fetch-farms`, 124 call units in total; the
Kelmarsh series is reused from [era5-check-pl](https://github.com/kwwaleszczynski-del/era5-check-pl)).

## 10.2 Three series per site

**Theory** is what a planner without operating data would use. ERA5 wind at 10 and 100 m is taken to hub height
with the shear of the hour,

```math
v_{hub} = v_{100}\left(\frac{h}{100}\right)^{\alpha}, \qquad \alpha = \frac{\ln(v_{100}/v_{10})}{\ln 10} \in [0, 0.6],
```

averaged over the two ends of the hour (ERA5 wind is instantaneous, SCADA hours are means), and put through
the manufacturer's catalogue power curve at 1.225 kg/m³ without any losses. Only the MM92 has a public catalogue
curve (Open Energy Platform turbine library, distributed with windpowerlib). For the MM82 and the SWT-2.3-82 the
MM92 curve is moved to their specific power $`SP`$ (rated power over rotor area),

```math
f(v) = f_{MM92}\big(k v\big), \qquad k = \left(\frac{SP_{MM92}}{SP}\right)^{1/3},
```

which keeps the power coefficient of the MM92 at the scaled speed and gives the turbine's own rotor area in
partial load, $`P(v) = \tfrac{1}{2}\rho A\,C_p(kv)\,v^3`$ (checked in `tests/test_theory.py`).

**This model** is the series this repository uses for every Polish grid point: ERA5 at 100 m, the MM92 curve
measured at Kelmarsh and smoothed for a farm, 10 % wake and electrical losses.

**Operation** is the hourly mean of power over rated power of the turbines with data, leaving out the
turbine-intervals that are half or more in a logged forced outage: technical availability is the structural
model's job, so it must not be counted twice. An hour needs three quarters of the turbine-intervals with data.
The series with forced outages counted as zeros is reported as well; it is 0.6-2.4 % lower at the wind farms.

**PV.** Theory is the PV model of this repository (Erbs, isotropic transposition, Faiman, 10 % system losses,
inverter clipping) with the nameplate data of PVDAQ 34 (DC/AC 1.086, tilt 11.2°, temperature coefficient
-0.45 %/K) and ERA5 GHI and 2 m temperature. Operation is the measured AC power; the 420 hours of the 2016
inverter outage are left out.

## 10.3 Energy, seasons and wind speed

Production ratio = operation / theory over the common valid hours (`results/theory_sites.csv`):

| Site | CF theory | CF this model | CF operation | ratio | by year | DJF | MAM | JJA | SON |
|---|---|---|---|---|---|---|---|---|---|
| Kelmarsh | 36.2 % | 37.4 % | 29.1 % | 0.80 | 0.78-0.81 | 0.83 | 0.84 | 0.77 | 0.77 |
| Penmanshiel | 35.7 % | 42.3 % | 32.6 % | 0.91 | 0.90-0.93 | 0.87 | 0.95 | 0.98 | 0.90 |
| Hill of Towie | 22.9 % | 33.7 % | 24.6 % | 1.08 | (one year) | 1.05 | 1.09 | 1.17 | 1.06 |
| PVDAQ 34 | 21.6 % | (= theory) | 19.2 % | 0.89 | 0.89 | 0.85 | 0.90 | 0.91 | 0.89 |

Mean capacity factor by ERA5 hub-height wind (`results/theory_wind_bins.csv`), selected bins:

| Bin, m/s | Kelmarsh theory / operation | Penmanshiel theory / operation | Hill of Towie theory / operation |
|---|---|---|---|
| 5-6 | 0.15 / 0.13 | 0.11 / 0.12 | 0.10 / 0.13 |
| 8-9 | 0.57 / 0.42 | 0.46 / 0.42 | 0.41 / 0.45 |
| 11-12 | 0.99 / 0.80 | 0.94 / 0.79 | 0.89 / 0.78 |
| 14-15 | 1.00 / 0.96 | 1.00 / 0.92 | 1.00 / 0.93 |
| 17-18 | 1.00 / 0.96 | 1.00 / 0.68 | (4 h) |

- The sign of the energy error changes from site to site: ERA5 is too windy at Kelmarsh (era5-check-pl found
  +10 % speed at hub height there) and not windy enough on the hills of Hill of Towie. At one site it is stable:
  0.78-0.81 in each of six years at Kelmarsh, 0.89 in each of three years for PV.
- The catalogue curve climbs faster and reaches rated power at 12 m/s; the farms reach 92-96 % only at 14-15 m/s
  (wakes, turbulence, the spread of wind over a farm). Penmanshiel's output falls above 16 m/s (high-wind
  derating and cut-out), which no catalogue curve shows.
- This model, built for an average Polish site, gives 29-37 % more energy than the three farms delivered. Part of
  it is the turbine (the MM92 has 308 W/m², the MM82 388 W/m² and the SWT-2.3-82 436 W/m²), part the hub height
  (100 m against 59-79 m) and the rest losses beyond 10 %.

## 10.4 Markov structure

Monthly five-state chains, thesis thresholds, energy-preserving levels, fitted on the same valid hours
(`results/theory_markov.csv`). Occupancy and mean residence time per state, theory / operation:

| Site | state 0 (< 12.5 %) | 1 | 2 | 3 | 4 (> 87.5 %) | ACF 1 h | ACF 24 h |
|---|---|---|---|---|---|---|---|
| Kelmarsh, occupancy | 0.33 / 0.38 | 0.28 / 0.33 | 0.15 / 0.14 | 0.09 / 0.09 | 0.14 / 0.06 | 0.984 / 0.951 | 0.40 / 0.39 |
| Kelmarsh, residence, h | 13.8 / 8.3 | 6.4 / 4.1 | 4.1 / 2.5 | 3.2 / 2.6 | 10.4 / 4.8 | | |
| Penmanshiel, occupancy | 0.36 / 0.40 | 0.25 / 0.25 | 0.14 / 0.13 | 0.09 / 0.11 | 0.15 / 0.11 | 0.987 / 0.947 | 0.40 / 0.32 |
| Penmanshiel, residence, h | 19.2 / 9.7 | 6.9 / 3.2 | 4.3 / 2.0 | 3.4 / 2.2 | 11.9 / 5.1 | | |
| Hill of Towie, occupancy | 0.52 / 0.50 | 0.26 / 0.24 | 0.11 / 0.12 | 0.06 / 0.09 | 0.06 / 0.05 | 0.985 / 0.952 | 0.39 / 0.40 |
| Hill of Towie, residence, h | 24.5 / 14.0 | 7.1 / 3.6 | 4.3 / 2.3 | 3.3 / 2.7 | 8.0 / 3.4 | | |

The operating series stays half as long in every state: hour-to-hour variability of real farms that a 31 km
reanalysis cell smooths away. The full-output state is much rarer in operation (5-11 % of hours against 6-15 %),
because the farms do not sit at rated power. The largest difference in a single transition probability between
the homogeneous matrices of operation and theory is 0.16-0.25. PV (clear-sky-index chain): theory puts 82 % of
daylight hours in the clearest state, operation 56 %.

## 10.5 Adequacy

Each series becomes a monthly chain; the plant is 10 % of the peak load; structural reliability is the
derating of many independent units with the literature FOR (wind 5 %, PV 3 %, variant OZE-LIT) or the plant's
own FOR from its log (forced outages of an hour or more: 0.44, 1.49 and 2.17 %; PVDAQ 34: 5.7 %, variant OZE-EMP).
LOLE, EENS and ELCC are exact (COPT); sequential Monte Carlo with 2000 years gives frequency and duration.
Wind uses the KSE load and RTS-79; PV only RTS-79, aligned by the local clock, because the KSE load belongs to
Polish weather years. `results/theory_adequacy.csv`:

| Site, load | theory, literature FOR | theory, own FOR | theory x ratio | this model | operation | LOLE operation, h/yr |
|---|---|---|---|---|---|---|
| Kelmarsh, KSE | 32.5 % | 33.6 % | 28.5 % | 37.1 % | 28.3 % | 1.60 |
| Penmanshiel, KSE | 32.7 % | 33.6 % | 31.5 % | 42.1 % | 29.7 % | 1.55 |
| Hill of Towie, KSE | 21.6 % | 22.1 % | 23.1 % | 33.3 % | 23.5 % | 1.78 |
| Kelmarsh, RTS-79 | 29.9 % | 31.0 % | 26.1 % | 34.7 % | 25.9 % | 5.39 |
| Penmanshiel, RTS-79 | 27.1 % | 27.8 % | 26.0 % | 36.2 % | 25.1 % | 5.52 |
| Hill of Towie, RTS-79 | 18.5 % | 18.9 % | 20.0 % | 29.4 % | 21.0 % | 6.03 |
| PVDAQ 34, RTS-79 | 21.3 % | 20.8 % | 19.2 % | | 17.4 % | 6.55 |

(Without the plant LOLE is 3.00 h/yr with the KSE load and 9.39 h/yr with RTS-79.) Sequential Monte Carlo, KSE load, own FOR:

| Site | LOLE theory / operation, h/yr | events per year, theory / operation (± 95 %) | mean event, h |
|---|---|---|---|
| Kelmarsh | 1.36 / 1.56 | 0.31 / 0.38 (± 0.05) | 4.4 / 4.1 |
| Penmanshiel | 1.33 / 1.47 | 0.31 / 0.37 (± 0.05) | 4.3 / 3.9 |
| Hill of Towie | 1.77 / 1.71 | 0.41 / 0.41 (± 0.05) | 4.3 / 4.1 |

- **Uncalibrated theory misstates ELCC by -1.4 to +5.3 points** at the wind farms and +3.4 points for PV, and the
  error follows the error in energy.
- **One production ratio fixes most of it.** Scaled by the farm's own ratio, the theoretical chain lands within
  0.2 (Kelmarsh), 1.8 (Penmanshiel) and 0.4 points (Hill of Towie) of the operating one, and within 1.7 points for
  PV, whose loss is larger in winter (ratio 0.85) than in summer (0.91), when RTS-79 has its peak.
- **The literature FOR matters less.** Replacing 5 % by the plant's own FOR moves ELCC by 0.5-1.2 points.
- **Chronology:** the operating chains give slightly more frequent and slightly shorter loss-of-load events, by
  about one Monte Carlo interval or less; expected LOLE and ELCC depend only on the hourly distribution.

## 10.6 What transfers to Poland, and what does not

Transfers: the ranking (energy level first, then shape, then the FOR source), the fact that a farm's ratio
is stable from year to year, and that operating series are less persistent than reanalysis ones. Does not
transfer: the ratios themselves. ERA5's wind bias belongs to the site, with opposite signs at the three UK farms,
and in Poland era5-check-pl found a median 10 m bias of +8 % with a range of -14 % to +54 % at 17 lowland and
coastal stations. One US desert rooftop says nothing about the losses of Polish PV. There is no open SCADA of a
Polish wind farm to measure the Polish ratio.

What an error in energy would do to the Polish results is a sensitivity, not a correction
(`results/theory_transfer.csv`; KSE load, wind 10 % of peak, every grid point scaled before averaging):

| Energy of the Polish series / this model | 0.7 | 0.8 | 0.9 | 1.0 | 1.1 | 1.2 |
|---|---|---|---|---|---|---|
| wind ELCC, historical weather | 19.3 % | 21.8 % | 24.2 % | 26.5 % | 28.8 % | 31.0 % |
| wind ELCC, monthly chain | 22.4 % | 25.3 % | 28.1 % | 30.8 % | 33.4 % | 35.9 % |

Every 10 % of energy error moves the portfolio's ELCC by 2.2-2.5 points (2.5-2.9 for the chain). Hypothetically,
if the Polish fleet related to this model as each UK farm does (the farm's quantile map from this model to
operation applied to every grid point), the energy would be 0.59-0.70 of the model's and the ELCC 15.4-18.4 %
with historical weather; at the same energy a pure ratio gives about 1 point more, so the shape of a real farm
adds a little to the energy effect. These numbers say how large the question is, not what the answer is in
Poland.
