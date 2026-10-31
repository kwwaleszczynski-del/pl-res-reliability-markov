# Verification of the weather models

Part of [pl-res-reliability-markov](../README.md). Ten weather years are a small sample, and the regime model was
fitted on them. Three checks on 2015-2024 (bootstrap, ten-year block test, leave one year out) and a test
against 35 earlier ERA5 years.

## Method

Ten weather years are a small sample and the regime model is fitted on them, so three checks are added
(`run --stages verify`):

- **Bootstrap over weather years.** The ten years are resampled with replacement (2000 times for the exact COPT
  results, 1000 times for Monte Carlo, where the simulated years inside every chosen weather year are resampled
  too). LOLE is computed once per year on a grid of load shifts (for ELCC) or storage energies (for $`E^*`$); a
  replicate is a weighted mean of these curves, read off by linear interpolation. The same draws serve every
  model, so differences between models get intervals as well. Intervals are the 5 and 95 % percentiles.
- **Ten-year block test.** Every model is simulated as 100 independent decades, each with one synthetic version
  of every weather year (with its temperatures and its load). A statistic of the ERA5 decade (lull counts, the
  longest lull, the weakest week, autocorrelation, cold-calm days, wind in the highest-load hours, wind ELCC) is
  placed among the 100 simulated decades: its percentile and a two-sided p-value. For the whole distribution of
  lull lengths (1, 2, 3, 4, 5-6 and 7+ days) a chi-square distance to the model mean is compared the same way.
- **Leave one year out.** Each year in turn is removed, all three models are fitted on the other nine and 100
  synthetic versions of the held-out year are drawn (the temperature model with that year's temperatures,
  classified with the climatology of the nine fitting years). Calibration of every annual statistic is the
  probability integral transform (PIT): the share of synthetic versions below the observed value, uniform on
  [0, 1] for a calibrated model. The wind ELCC of each held-out year (with its own load) is compared with the
  model's prediction.
- **Earlier years.** `run --stages long` tests the models fitted on 2015-2024 against ERA5 1980-2014 (and PSE
  demand 2010-2014): distributions, PIT and return levels of the longest lull and of the weakest week (results
  below).

## Results on 2015-2024

Ten weather years are a small sample, and the regime model was fitted on them. Three checks (see Method above): a
bootstrap over the weather years, the ERA5 decade against 100 simulated decades of every model, and every year
predicted by models fitted on the other nine.

**Bootstrap, 90 % intervals** (same draws for every model):

| | historical weather | Markov, monthly | regime switching | regimes + temperature |
|---|---|---|---|---|
| wind ELCC without storage (10 % of peak), % | **26.5** (25.3-27.9) | 30.7 (30.3-31.2) | 30.7 (30.3-31.2) | 28.6 (27.8-29.5) |
| difference to historical weather, points | | +2.8 to +5.7 | +2.7 to +5.6 | +0.9 to +3.4 |
| 4 h storage ELCC, fleet 5 % of peak (2 x 10 % RES), % | **73.0** (71.3-74.8) | 74.0 (72.1-75.9) | 74.4 (72.4-76.3) | 74.6 (72.2-76.5) |
| 4 h storage ELCC, fleet 20 % of peak, % | **49.0** (48.2-49.7) | 48.9 (47.9-49.9) | 49.1 (48.2-50.0) | 49.6 (48.5-50.5) |
| storage for 30 % wind + PV, GWh (KSE scale) | **117** (100-135) | 65 (60-71) | 83 (74-95) | 101 (91-112) |
| ratio to historical weather | | 0.48-0.66 | 0.60-0.86 | 0.76-0.99 |

**Ten-year block test**: percentile of the ERA5 decade among 100 simulated decades of each model (outside 5-95 %
marked *):

| Statistic of the decade | Markov, monthly | regime switching | regimes + temperature |
|---|---|---|---|
| wind lulls of 3+ / 5+ days below 15 % | 100* / 100* | 94 / 83 | 92 / 86 |
| October-March lulls of 3+ days | 100* | 80 | 83 |
| lulls of 3+ days below 10 %; wind + PV lulls of 3+ days | 100* ; 100* | 68 ; 47 | 77 ; 39 |
| longest lull of the year | 100* | 83 | 90 |
| distribution of lull lengths (chi-square) | 100* | 57 | 50 |
| hourly spells below 12.5 % of 48 h+ | 98* | 34 | 35 |
| weakest week of the year (lowest 7-day mean CF) | 0* | 0* | 0* |
| autocorrelation at 24 h | 100* | 2* | 1* |
| coldest-third winter days with CF below 15 % | 100* | 100* | 70 |
| wind CF in the 100 highest-load hours | 0* | 0* | 1* |
| wind ELCC without storage | 0* | 0* | 0* |

![Verification](img/verification.png)

- **The conclusions survive the small sample.** The monthly chain overstates wind ELCC by 2.8-5.7 points and
  understates the storage needed for lulls by 34-52 %; both differences exclude zero. Storage of 1-8 h is valued
  the same by every model within the intervals (the differences to historical weather include zero in 14 of 15
  comparisons), so for short storage the choice of weather model really is secondary.
- **The monthly chain fails as a model of lulls.** The ERA5 decade has more lulls of every kind than all 100
  simulated decades, and in leave-one-year-out every held-out year has more lulls than the chain predicts (mean
  PIT 0.93-0.96; its 80 % intervals contain the observed value in 2-3 of 10 years).
- **The regime models pass on lulls**, including the whole distribution of lull lengths (chi-square at the 50-57th
  percentile), winter lulls, dark calms of wind and PV and the longest lull of the year. Leave-one-year-out agrees:
  PIT 0.53-0.65, 80 % intervals covering 7-10 of 10 years. Only the temperature model also gets the cold, calm
  winter days right (70th percentile, PIT 0.65; without temperature the 100th).
- **Where they still fail.** The weakest week of the year is too mild (8.9 % against 7.7 % CF; every simulated
  decade is milder), the hourly autocorrelation at 24 h is slightly too high (0.55 against 0.52), and the coupling
  with demand is only half there: wind output in the highest-load hours and wind ELCC stay above the data in
  almost every simulated decade. This is why the temperature model still asks for 1-24 % less storage than
  historical weather.
- **It predicts the bad years.** Fitted on nine years and driven by the temperatures of the tenth, the temperature
  model predicts the wind ELCC of the held-out year with a correlation of 0.74 (RMSE 3.2 points, bias +2.6). The
  monthly chain has no skill (r = 0.06, bias +4.1), memory without temperature little (r = 0.31).

## Earlier weather years: ERA5 1980-2014

The three weather models keep their 2015-2024 parameters and are tested against 35 earlier ERA5 years
(`run --stages long`, data fetched with `--data-dir data/long fetch-weather --start-year 1980 --end-year 2014`;
results in [`results/verification_long.json`](../results/verification_long.json) and
[`results/verification_long.csv`](../results/verification_long.csv), one row per model, year and statistic). Every
test year gets 40 synthetic versions; the temperature model is driven by that year's temperatures,
classified with the 2015-2024 climatology. PSE demand exists for 2010-2014 only, so wind in the highest-load hours
and wind ELCC are tested on those five years.

Model columns: mean of the synthetic years, then (mean PIT; share of the 35 years inside the model's 80 % interval).
A calibrated model has PIT about 0.5 and 80 % coverage.

| Statistic | ERA5 1980-2014 (90 % interval) | ERA5 2015-2024 | Markov, monthly | regime switching | regimes + temperature |
|---|---|---|---|---|---|
| wind lulls of 3+ days below 15 %, per year | **11.7** (10.9-12.5) | 12.5 | 6.1 (0.95; 20 %) | 11.4 (0.52; 80 %) | 12.1 (0.46; 80 %) |
| wind lulls of 5+ days below 15 %, per year | **3.7** (3.1-4.3) | 4.2 | 1.1 (0.86; 40 %) | 3.6 (0.50; 80 %) | 3.8 (0.49; 77 %) |
| October-March lulls of 3+ days, per year | **3.4** (2.9-3.9) | 3.6 | 0.9 (0.90; 37 %) | 3.1 (0.53; 86 %) | 3.7 (0.44; 89 %) |
| lulls of 3+ days below 10 %, per year | **3.9** (3.4-4.3) | 3.7 | 0.9 (0.94; 29 %) | 3.4 (0.56; 80 %) | 3.6 (0.53; 89 %) |
| wind + PV lulls of 3+ days below 10 %, per year | **2.3** (1.9-2.7) | 2.2 | 0.7 (0.82; 49 %) | 2.2 (0.52; 74 %) | 2.6 (0.45; 86 %) |
| longest lull of the year, days | **8.0** (7.1-9.0) | 9.0 | 5.4 (0.79; 60 %) | 8.2 (0.44; 80 %) | 8.1 (0.43; 77 %) |
| weakest week of the year (lowest 7-day mean CF) | **0.087** (0.082-0.091) | 0.077 | 0.113 (0.09; 20 %) | 0.088 (0.46; 69 %) | 0.088 (0.47; 66 %) |
| hourly spells below 12.5 % of 48 h+, per year | **5.7** (5.1-6.3) | 4.5 | 3.3 (0.81; 63 %) | 4.8 (0.62; 77 %) | 5.1 (0.58; 83 %) |
| wind CF in the 100 highest-load hours (2010-2014) | **0.312** (0.255-0.377) | 0.262 | 0.379 (0.20; 20 %) | 0.384 (0.20; 20 %) | 0.307 (0.45; 60 %) |
| coldest-third winter days with CF below 15 %, share | **0.238** (0.213-0.263) | 0.248 | 0.083 (0.94; 23 %) | 0.164 (0.74; 66 %) | 0.243 (0.49; 80 %) |

- **The monthly chain fails again, in every decade.** It produces about half of the observed 3-day lulls (6.1
  against 11.7 a year) and a third of the 5-day lulls (1.1 against 3.7); its 80 % interval contains the observed
  count of 3-day lulls in 20 % of the years. In 1980, for example, it gives 6.5 lulls against 9 observed, and the
  gap is larger in 29 of the 35 years.
- **The regime models pass out of sample on lulls**, fitted on one decade and tested on three and a half others:
  lull counts, 5-day lulls, winter lulls, dark calms of wind and PV, the longest lull of the year and the weakest
  week, with PIT 0.43-0.62 and 80 % intervals containing 66-89 % of the years. The weakest week, too mild against
  2015-2024 (above), is right against 1980-2014: that decade had an unusually weak week (CF 0.077 against 0.087 on
  average in 1980-2014).
- **Only the temperature model gets the cold calms right** (0.243 against 0.238; without temperature 0.164), and
  only it gets wind in the highest-load hours of 2010-2014 (0.307 against 0.312; the other two 0.379-0.384). Wind
  ELCC of 2010-2014 (10 % of peak, demand of the same year) is 22.9-29.8 % in the data; the temperature model
  overstates it by +1.6 points on average (-1.1 to +4.2), the monthly chain and the regime model without temperature
  by +5.2 and +5.3.
- **The rarest lulls are a little short.** The longest lull exceeded once in 20 years is 15.3 days in the 35 ERA5
  years and 13 days in both regime models (monthly chain 8); at 2, 5 and 10 years they agree within a day (7.0 vs 8,
  9.0 vs 10, 10.6 vs 11). The 20-year value rests on one or two observed lulls, so the gap is uncertain, but it
  points the same way as the weakest-week result of 2015-2024: the regime models are slightly too mild in the
  extreme tail.
- 1980-2009 had fewer 3-day lulls than 2015-2024 (11.3 against 12.5 a year) and 2010-2014 more (13.8).
