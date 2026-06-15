# Tests for clamping of covariate/exact/moderation lags to the event window.
# A sub-experiment only spans event times in [kappa_pre, kappa_post], so a lag k
# can be constructed only if k <= |kappa_pre|. Deeper lags are clamped, not dropped.

va <- stats::setNames(c("id", "time", "y", "d"),
                      c("id", "year", "outcome", "D"))

test_that("deep lags are clamped to |kappa_pre|", {
  spec <- cbwsdid:::parse_feature_formula_cbwsdid(
    ~ lag(outcome, 1:10), var_aliases = va, kappa = c(-2, 2)
  )
  expect_setequal(spec$lag_k, c(1L, 2L))
  expect_setequal(spec$feature_name, c("outcome_l1", "outcome_l2"))
  expect_setequal(spec$et_lookup, c(-1L, -2L))
  expect_true(all(spec$data_var == "y"))
})

test_that("a fully out-of-window range collapses to a single boundary lag", {
  spec <- cbwsdid:::parse_feature_formula_cbwsdid(
    ~ lag(outcome, 5:10), var_aliases = va, kappa = c(-2, 2)
  )
  expect_equal(nrow(spec), 1L)
  expect_equal(spec$feature_name, "outcome_l2")
  expect_equal(spec$lag_k, 2L)
})

test_that("lags inside the window are left unchanged", {
  spec <- cbwsdid:::parse_feature_formula_cbwsdid(
    ~ lag(outcome, 1:2) + lag(x1, 1:2), var_aliases = va, kappa = c(-3, 2)
  )
  expect_setequal(
    spec$feature_name,
    c("outcome_l1", "outcome_l2", "x1_l1", "x1_l2")
  )
})

test_that("kappa = NULL disables clamping (backward compatible)", {
  spec <- cbwsdid:::parse_feature_formula_cbwsdid(
    ~ lag(outcome, 1:10), var_aliases = va, kappa = NULL
  )
  expect_equal(max(spec$lag_k), 10L)
})

test_that("a single informative warning is emitted when lags exceed the window", {
  expect_warning(
    cbwsdid:::warn_clamped_lags_cbwsdid(
      list(~ lag(outcome, 1:10), ~ lag(x1, 1:8), NULL),
      var_aliases = va, kappa = c(-2, 2)
    ),
    "clamped to lag 2"
  )
  expect_silent(
    cbwsdid:::warn_clamped_lags_cbwsdid(
      list(~ lag(outcome, 1:2)),
      var_aliases = va, kappa = c(-3, 2)
    )
  )
})

test_that("clamping retains cohorts that the old behaviour would have dropped", {
  set.seed(2)
  n <- 80L; years <- 2000:2012
  u <- tibble::tibble(
    id = seq_len(n), x1 = stats::rnorm(n),
    adopt = sample(c(2006, 2008, Inf), n, replace = TRUE, prob = c(.3, .3, .4))
  )
  dat <- tidyr::expand_grid(id = u$id, year = years) |>
    dplyr::left_join(u, by = "id") |>
    dplyr::group_by(id) |>
    dplyr::arrange(year, .by_group = TRUE) |>
    dplyr::mutate(
      outcome = cumsum(stats::rnorm(dplyr::n(), 0.1)) + 0.4 * x1,
      D = ifelse(is.infinite(adopt), 0L, as.integer(year >= adopt))
    ) |>
    dplyr::ungroup()

  db <- dat |>
    dplyr::select(id, year, outcome, D) |>
    dplyr::rename(time = year, y = outcome, d = D) |>
    dplyr::arrange(id, time) |>
    dplyr::group_by(id) |>
    dplyr::mutate(
      treated = ifelse(sum(d) == 0, 0, 1),
      A = ifelse(treated == 0, Inf, min(ifelse(d == 1, time, Inf)))
    ) |>
    dplyr::ungroup()

  cands <- sort(unique(db$A[is.finite(db$A)]))
  kap <- c(-2, 2)

  retained <- function(covf) {
    res <- lapply(cands, function(a) {
      s <- cbwsdid:::subexp.row.construct(db, a, kap)
      cbwsdid:::subexp.refine(
        s, var_aliases = va, kappa = kap,
        refinement.method = "none", covs.formula = covf
      )
    })
    sa <- purrr::list_rbind(purrr::map(res, "subexp.refined"))
    if (nrow(sa) == 0) 0L else dplyr::n_distinct(sa$subexp.id)
  }

  deep <- suppressWarnings(retained(~ lag(outcome, 1:5)))   # needs et = -3,-4,-5 (absent)
  shallow <- retained(~ lag(outcome, 1:2))

  expect_gt(deep, 0L)
  expect_equal(deep, shallow)
})
