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
- `cbwsdid_qoi()` extracts simple or dynamic quantities of interest.
- `cbwsdid_subexperiments()` summarizes the raw and refined donor pools across
  subexperiments.

## Minimal example

```r
library(cbwsdid)

# --- 1. Weighted stacked DiD ----

fit.ws <- cbwsdid(
  data = panel_df, # dataset with id-time as unit of analysis
  y = "y", # dependent variable
  d = "d", # treatment variable
  id = c("unit", "time"),
  kappa = c(-4, 4), # event-time window
  refinement.method = "none" # no refinement = Weighted stacked DiD
)

fit.s <- update(fit.ws, weights = NULL) # to get stacked DiD

# --- 2. Covariate-balanced Weighted stacked DiD ----

fit.cbws <- cbwsdid(
  data = panel_df,
  y = "y",
  d = "d",
  id = c("unit", "time"),
  kappa = c(-4, 4),
  refinement.method = "matchit" # for matching
  covs.formula = ~ lag(y, 1:3) + lag(x1, 1:3), 
  refinement.args = list( # arguements for MatchIt
                      method = "nearest",
                      distance = "mahalanobis",
                      ratio = 4,
                      replace = T
                    )
)

model <- attr(fit.cbws, "cbwsdid") # all information about the model

cbwsdid_balance(fit.cbws) # balance statistics by variables in `covs.formula`

cbwsdid_qoi(fit, type = "simple") # for att
cbwsdid_qoi(fit, type = "dynamic") # for event study estimates 

iplot(fit) # event study plot using fixest function
```
