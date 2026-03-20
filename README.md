# cbwsdid

`cbwsdid` implements weighted stacked difference-in-differences and
covariate-balanced weighted stacked DID for panel data.

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
- `cbwsdid_qoi()` extracts simple or dynamic quantities of interest.
- `cbwsdid_subexperiments()` summarizes the raw and refined donor pools across
  subexperiments.

## Minimal example

```r
library(cbwsdid)

fit <- cbwsdid(
  data = panel_df,
  y = "y",
  d = "d",
  id = c("unit", "time"),
  kappa = c(-4, 4),
  refinement.method = "none"
)

cbwsdid_qoi(fit, type = "simple")
cbwsdid_qoi(fit, type = "dynamic")
```
