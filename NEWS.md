# cbwsdid 0.2.0

- Units with missing refinement features are no longer dropped silently:
  treated exclusions now raise an error unless `allow_treated_drop = TRUE`
  (then a warning), control exclusions emit a message, and doubly robust fits
  warn when rows cannot be residualized because outcome-model features or the
  base-period outcome are missing.
- Cohort-level analytic simple standard errors now zero-fill units with
  partial post-period coverage.
- For switch designs, `history.length` is now optional: it defaults to
  `abs(kappa[1])`, cannot exceed it, and shorter values warn that control
  treatment paths deeper than the matched history are unrestricted inside the
  event window.
- `cbwsdid_balance()` now computes `Diff.Un` on the full design sample,
  including units that refinement discarded. 
- Corrected moderated second-stage models fitted via `moderation.formula`.
  They now include moderator-by-event-time main effects (and post-by-moderator
  main effects in the standalone pooled moderated model).
- Added `cbwsdid_qoi_moderation_slopes()` for extracting the event-time path of
  moderation slopes, including pre-treatment placebo slopes.
- Added pipeline-bootstrap inference for
  `cbwsdid_qoi(model_type = "moderated", se.method = "bootstrap")` at fixed
  moderator profiles.
- Redefined the overall `type = "pooled"` quantity as one equal-weighted
  pre-/post contrast of the dynamic coefficients. Conditional, analytic, and
  bootstrap inference now share exactly the same point estimate, and overall
  pooled QOIs no longer require `pooled = TRUE`.
- Added analytic influence-function standard errors for cohort-level dynamic
  effects and their simple post-treatment averages.
- Added unified inference selection through
  `cbwsdid_qoi(se.method = c("conditional", "analytic", "bootstrap"))`.
  Analytic inference uses unit-clustered influence-function contributions and
  exposes the full event-time covariance matrix.
- `cbwsdid()` now always stores compact analytic-variance components and, with
  the default `keep_data = TRUE`, stores the raw panel so pipeline-bootstrap
  QOIs no longer require users to pass `data=` again.
- Added optional doubly robust CBWSDID estimation via `outcome.formula`, which
  residualizes within-stack outcome changes using pre-treatment covariates
  before fitting the usual weighted stacked DID regression.
- Added `cbwsdid_bootstrap()`, a post-estimation pipeline cluster bootstrap
  that reruns subexperiment construction, refinement, stacked weighting, and
  QOI aggregation on resampled clusters.

# cbwsdid 0.1.1

- Covariate/exact/moderation lags requested beyond the pre-treatment window are
  now clamped to `|kappa_pre|` with a single informative warning, instead of
  silently emptying the design sample and dropping the cohort.
- The cohort quantity of interest now builds its cohort-by-event-time design with
  a single sparse `fixest::i()` interaction instead of one dense dummy column per
  cell. Results are identical, but memory/time no longer blow up for switch
  designs with many episode types.
- The cohort dynamic output now reports per-cell `n_treated` (number of treated
  units contributing to each cohort/episode-by-event-time cell).
- Added `testthat` unit tests for lag clamping, the stacked-weight identities,
  and the cohort quantity of interest.

# cbwsdid 0.1.0

- New QOI -- cohort
- Simulated dataset for examples is added
- Vignette on usage

# cbwsdid 0.0.0.1

- Initial development version.
- Added `cbwsdid()` for weighted stacked DID and covariate-balanced weighted
  stacked DID.
- Added post-estimation helpers for balance, quantities of interest, and
  subexperiment summaries.
