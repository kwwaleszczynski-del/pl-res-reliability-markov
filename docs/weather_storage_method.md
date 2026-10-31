# Regime-switching weather and storage: method

Part of [pl-res-reliability-markov](../README.md) (methods 7 and 8; methods 1-6 are in the README, the
verification is in [verification.md](verification.md)).

## 7. Regime-switching weather model

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

## 8. Storage in the sequential Monte Carlo

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
