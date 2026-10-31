# Weather-dependent failures and repairs: method and details

Part of [pl-res-reliability-markov](../README.md) (method 11; the summary is in the README section
[Weather-dependent failures and repairs](../README.md#weather-dependent-failures-and-repairs)). All numbers come
from `python -m res_reliability run --stages failures` and are in `results/failure_rates.csv`,
`results/failure_tests.json`, `results/integration_farms.csv`, `results/integration_pl.csv` and
`results/summary.json`.

## 11.1 Independence, and why the Kronecker sum does not remove it

The thesis composed the production model (weather) and the structural model (failures) by assuming them
independent: a plant is in production state $`i`$ and available with probability $`p_{prod}(i)\,A`$. It checked
this with the Kronecker sum of the two intensity matrices,

```math
\mathbf{\Lambda}_{int} = \mathbf{\Lambda}_{prod} \otimes \mathbf{I}_2 + \mathbf{I}_N \otimes \mathbf{\Lambda}_{struct},
\qquad
\mathbf{\Lambda}_{struct} = \begin{bmatrix} -\lambda & \lambda \\ \mu & -\mu \end{bmatrix},
```

and found the same LOLE. That is not a coincidence. With constant $`\lambda`$ and $`\mu`$ the product of the
stationary marginals solves the integrated chain,

```math
(\boldsymbol{\pi}_{prod} \otimes \boldsymbol{\pi}_{struct})\,\mathbf{\Lambda}_{int}
= (\boldsymbol{\pi}_{prod}\mathbf{\Lambda}_{prod}) \otimes \boldsymbol{\pi}_{struct}
+ \boldsymbol{\pi}_{prod} \otimes (\boldsymbol{\pi}_{struct}\mathbf{\Lambda}_{struct}) = \mathbf{0},
```

so the Kronecker sum *is* the independence assumption written as a chain (checked to 1e-13 in
`tests/test_integration.py`). A correlation between weather and failures needs intensities that depend on the
production state. The block generator keeps the production transitions in both structural states and lets a unit
fail with $`\lambda_i`$ and be repaired with $`\mu_i`$ while the weather is in state $`i`$ ($`\mathbf{e}_i`$ is
the $`i`$-th unit vector):

```math
\mathbf{\Lambda}_{int} = \mathbf{\Lambda}_{prod} \otimes \mathbf{I}_2 + \sum_{i} \mathbf{e}_i\mathbf{e}_i^{\top} \otimes
\begin{bmatrix} -\lambda_i & \lambda_i \\ \mu_i & -\mu_i \end{bmatrix}
```

With $`\lambda_i = \lambda`$, $`\mu_i = \mu`$ it is the Kronecker sum again. Its stationary distribution gives the
availability in each production state, $`A_i = \pi(i, \mathrm{up}) / \pi_{prod}(i)`$; the production marginal is
unchanged because the weather does not depend on the turbines.

## 11.2 Events, exposure and rates

- **Forced outages.** Kelmarsh and Penmanshiel: the `Stop` records the operator files under the IEC 61400-26
  category "Forced outage"; Hill of Towie: the 49 alarm codes that
  [wind-om-reliability](https://github.com/kwwaleszczynski-del/wind-om-reliability) classed as forced outages from
  the 2021 SCADA (stopped, turbine-error timer on); an alarm without an end lasts until the turbine's next alarm,
  at most a week. Stops of one turbine within 30 minutes form one outage; outages of an hour or more count, as for
  the structural parameters of the main study.
- **Exposure.** A turbine is observed from commissioning (or the first log) to the end of the log, minus silences
  of more than a week (missing files) and communication losses of 6 h or more. This gives 33.2, 104.2 and 213.6
  turbine-years (wind-om-reliability: 32.7, 104.3 and 214.8 with its own rules).
- **Production state of an hour.** ERA5 at the farm through the same model as the Polish grid points (100 m, MM92,
  -10 %), discretised with the thesis thresholds, so the rates plug straight into the Polish chain. The highest
  state (> 87.5 %) is never reached with 10 % losses, here or in Poland. For the physics, the same is done with
  ERA5 hub-height wind in 3 m/s bins.
- **Rates.** Up-time $`U_i`$ and down-time $`D_i`$ in state $`i`$ are summed over turbines and hours; a failure
  belongs to the state of the hour it starts in, a repair to the state of the hour it ends in:

```math
\hat\lambda_i = \frac{n_i}{U_i}, \qquad \hat\mu_i = \frac{m_i}{D_i}
```

  the maximum-likelihood estimate of a piecewise-constant intensity, with exact Poisson intervals.
- **Tests.** Likelihood ratio of one rate against a rate per state,
  $`G^2 = 2\sum_i n_i \ln(\hat\lambda_i / \hat\lambda)`$ with $`S - 1`$ degrees of freedom, and a score test for a
  log-linear trend in the state index (one degree of freedom). Failures cluster (repeat faults, bad-actor turbines),
  so the Poisson tests are too confident; the dispersion $`\phi`$ (Pearson $`\chi^2`$ per degree of freedom over
  turbine x year x state cells) divides $`G^2`$ for the quasi-Poisson p-value. Across farms,
  $`\lambda_{f,i} = a_f\, b_i`$ (iterative proportional fitting) gives a common relative state effect $`b_i`$.

## 11.3 Results

Rates by production state (`results/failure_rates.csv`, failures per turbine-year of exposure, mean time to
repair):

| State (CF of this model) | Kelmarsh | Penmanshiel | Hill of Towie |
|---|---|---|---|
| 0 (< 12.5 %) | 3.9 /yr, 13.9 h | 6.8 /yr, 14.6 h | 8.1 /yr, 19.2 h |
| 1 (12.5-37.5 %) | 1.0 /yr, 28.0 h | 5.3 /yr, 18.1 h | 7.6 /yr, 18.1 h |
| 2 (37.5-62.5 %) | 2.3 /yr, 8.4 h | 5.8 /yr, 21.1 h | 7.7 /yr, 18.4 h |
| 3 (62.5-87.5 %) | 5.4 /yr, 11.5 h | 8.1 /yr, 11.4 h | 11.8 /yr, 20.8 h |

By ERA5 hub-height wind:

| Hub wind | Kelmarsh | Penmanshiel | Hill of Towie |
|---|---|---|---|
| < 12 m/s | 1.2-4.4 /yr, 12-26 h | 5.5-6.3 /yr, 15-19 h | 7.3-10.4 /yr, 18-20 h |
| 12-15 m/s | 4.0 /yr, 20.5 h | 5.4 /yr, 12.4 h | 12.3 /yr, 24.5 h |
| ≥ 15 m/s | 98 /yr, 4.4 h | 70 /yr, 5.1 h | 42 /yr, 27.2 h |
| share of outages starting at ≥ 15 m/s | 33 % | 15 % | 2.6 % |

Tests (`results/failure_tests.json`):

| | Kelmarsh | Penmanshiel | Hill of Towie |
|---|---|---|---|
| failures, 4 production states: Poisson p / dispersion / quasi-Poisson p | 7e-7 / 1.2 / 1e-5 | 3e-4 / 1.6 / 0.01 | 2e-13 / 2.9 / 7e-5 |
| repairs, 4 production states | 0.004 / 3.7 / 0.30 | 1e-8 / 4.8 / 0.04 | 0.14 / 7.7 / 0.87 |
| failures, 6 wind bins, quasi-Poisson p | 7e-30 | 3e-37 | 2.5e-9 |
| repairs, 6 wind bins, quasi-Poisson p | 0.07 | 0.004 | 0.82 |

- **Failures depend on the weather.** In the highest occupied production state the intensity is 1.32 times the
  average of all states (common effect of the three farms; 0.83-0.96 in the others); the size of the effect
  differs between farms (p = 0.002). Above 15 m/s at hub height the rate is 4-6 times (Hill of Towie), 11-13 times
  (Penmanshiel) and 20-80 times (Kelmarsh) the rate below 12 m/s. At Kelmarsh these are mostly tower-oscillation and pitch run-away trips, short protective stops
  (MTTR 4-5 h at the two Senvion farms).
- **Repairs do not.** No common effect of the production state on repair intensity (p = 0.54); the farm-level
  differences vanish once the clustering of repair times is accounted for, except a weak one at Penmanshiel.
  Repairs after storms are faster at the Senvion farms (the trips reset), slower at Hill of Towie.

## 11.4 The integrated chain at the farms

Hub-wind-bin chain of the log period (homogeneous), levels from the catalogue curve, the farm's own rates
(`results/integration_farms.csv`). Unavailable share of the turbines:

| | Kelmarsh | Penmanshiel | Hill of Towie |
|---|---|---|---|
| independence (Kronecker sum), every state | 0.43 % | 1.11 % | 1.88 % |
| state-dependent, < 3 m/s / 6-9 m/s / ≥ 15 m/s | 0.58 / 0.35 / 2.91 % | 1.09 / 1.13 / 2.60 % | 1.79 / 1.85 / 3.74 % |
| output lost to forced outages, independence / state-dependent | 0.43 / 0.57 % | 1.11 / 1.20 % | 1.88 / 2.06 % |
| largest distance of the joint distribution from the product of the marginals | 0.0005 | 0.0002 | 0.0004 |
| correlation of output level and unit availability | -0.016 | -0.006 | -0.013 |
| ELCC of 10 % of peak (KSE), independence / state-dependent | 27.69 / 27.68 % | 25.10 / 25.09 % | 19.25 / 19.24 % |
| hypothetical: no repairs at 12 m/s or more | 27.63 % | 24.99 % | 19.22 % |

The correlation is real but shifts unavailability into the windiest hours, which are the hours the system needs
the plant least; the energy lost grows by 8-34 % of a small number.

## 11.5 The Polish portfolio

Monthly chain of the national portfolio, wind 10 % of peak, rates per production state from each farm and pooled
(`results/integration_pl.csv`). Exact (COPT) values; Monte Carlo with 2000 years and common random numbers for
the pooled rates:

| Rates | Independent | State-dependent | Independent, same mean FOR |
|---|---|---|---|
| Kelmarsh, KSE | 30.68 % | 30.69 % | 30.68 % |
| Penmanshiel, KSE | 30.49 % | 30.48 % | 30.48 % |
| Hill of Towie, KSE | 30.28 % | 30.31 % | 30.33 % |
| pooled, KSE | 30.38 % (LOLE 1.520 h) | 30.40 % (1.519 h) | 30.41 % |
| pooled, RTS-79 | 29.01 % | 29.03 % | 29.02 % |
| pooled, KSE, Monte Carlo | 30.73 % (LOLE 1.37 ± 0.24 h) | 30.73 % (1.37 ± 0.24 h) | |

The national portfolio spends less time in its top state than a single farm, so the pooled constant rate, an
average over the farms' own weather, is slightly pessimistic for it: the state-dependent model gives 98.6 % mean
availability against 98.5 %.

**Hypothetical stress cases** (not estimated from data; pooled rates, top production state = portfolio CF above
62.5 %, 8 % of the hours):

| Case | mean availability | wind ELCC, KSE | of which correlation |
|---|---|---|---|
| H1: no repairs while the portfolio is in its top state (access restrictions in strong wind) | 98.3 % | -0.15 points | -0.11 |
| H2: H1 and five times the failure intensity in the top state | 97.1 % | -0.91 points (Monte Carlo -0.85) | -0.54 |

"Of which correlation" is the difference to an independent model with the same mean availability. Even a fleet
that fails five times more often in strong wind and cannot be repaired until it calms down loses less than a
point of ELCC, about half of it through the correlation.

## 11.6 Limits

- The rates come from three UK farms of 2.0-2.3 MW turbines built 2012-2016, with UK access rules and service
  contracts; the Polish fleet is younger, larger and serviced differently. What transfers is the size of the
  effect relative to the independence model, not the rates.
- The many-unit limit treats the turbines as independent given the weather path. Fleet-wide outages (grid, park
  controller) and serial faults are not separated from single-turbine faults; they are part of the dispersion.
- Production states come from ERA5 through this model, which overstates the farms' output by 29-37 %
  ([theory against operation](theory_vs_operation.md)); a state is a weather class, not the farm's own output.
- Hill of Towie's forced outages rest on alarm codes inferred from one year of SCADA.
