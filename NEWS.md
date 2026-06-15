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
