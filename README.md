# cbwsdid

`cbwsdid` implements weighted stacked difference-in-differences and
covariate-balanced weighted stacked DID for panel data.

The methodology implemented in the package is described in my preprint. If you use `cbwsdid`, please cite it:

Ustyuzhanin, Vadim (2026). Covariate-Balanced Weighted Stacked
Difference-in-Differences. arXiv preprint.
[https://doi.org/10.48550/arXiv.2604.02293](https://doi.org/10.48550/arXiv.2604.02293)

The package currently supports:

- absorbing staggered adoption designs;
- repeated `0 -> 1` switch-on designs;
- repeated `1 -> 0` switch-off designs via treatment recoding;
- within-subexperiment refinement via `MatchIt` or `WeightIt`;
- post-estimation helpers for balance, quantities of interest, and
  subexperiment summaries.

## Installation

```r
# install.packages("devtools")
devtools::install_github("vadvu/cbwsdid")
```

## Main functions

- `cbwsdid()` fits the estimator and returns a `fixest` model object.
- `cbwsdid_balance()` computes standardized mean difference tables based on
  refinement covariates.
- `cbwsdid_qoi()` extracts simple, dynamic, pooled, or cohort-level quantities
  of interest.
- `cbwsdid_subexperiments()` summarizes the raw and refined donor pools across
  subexperiments.

## Minimal example

```r
library(cbwsdid)

data(cbwsdid_sim)
cbwsdid_sim_small <- subset(cbwsdid_sim, id <= 150)

# --- 1. Weighted stacked DiD ----

fit.ws <- cbwsdid(
  data = cbwsdid_sim_small,
  y = "outcome",
  d = "D",
  id = c("id", "year"),
  kappa = c(-2, 2), # event-time window
  refinement.method = "none" # no refinement = Weighted stacked DiD
)

fit.s <- update(fit.ws, weights = NULL) # to get stacked DiD

# --- 2. Covariate-balanced Weighted stacked DiD ----

fit.cbws <- cbwsdid(
  data = cbwsdid_sim_small,
  y = "outcome",
  d = "D",
  id = c("id", "year"),
  kappa = c(-2, 1),
  refinement.method = "weightit",
  covs.formula = ~ lag(outcome, 1:2) + lag(x1, 1:2),
  refinement.args = list(method = "glm")
)

model <- attr(fit.cbws, "cbwsdid") # all information about the model

cbwsdid_balance(fit.cbws) # balance statistics by variables in `covs.formula`

cbwsdid_subexperiments(fit.cbws) # subexperiment diagnostics

cbwsdid_qoi(fit.cbws, type = "simple") # for ATT
cbwsdid_qoi(fit.cbws, type = "dynamic") # for event study estimates 
cohort_qoi <- cbwsdid_qoi(fit.cbws, type = "cohort") # subexperiment-level effects
cohort_qoi$dynamic
cohort_qoi$simple

iplot(fit.cbws) # event study plot using fixest function
```
