test_that("outcome.formula = NULL nests the weighting-only estimator", {
  skip_if_not_installed("fixest")
  data(cbwsdid_sim)
  dat <- subset(cbwsdid_sim, id <= 140)
  
  fit_default <- cbwsdid(
    data = dat,
    y = "outcome",
    d = "D",
    id = c("id", "year"),
    kappa = c(-2, 1),
    refinement.method = "none",
    pooled = FALSE
  )
  fit_null <- cbwsdid(
    data = dat,
    y = "outcome",
    d = "D",
    id = c("id", "year"),
    kappa = c(-2, 1),
    refinement.method = "none",
    pooled = FALSE,
    outcome.formula = NULL
  )
  
  expect_equal(
    cbwsdid_qoi(fit_null, type = "dynamic"),
    cbwsdid_qoi(fit_default, type = "dynamic"),
    tolerance = 1e-12
  )
  expect_false(attr(fit_null, "cbwsdid")$dr)
})


test_that("doubly robust fit runs and returns well-formed QOIs", {
  skip_if_not_installed("fixest")
  skip_if_not_installed("WeightIt")
  data(cbwsdid_sim)
  dat <- subset(cbwsdid_sim, id <= 150)
  
  fit <- suppressWarnings(cbwsdid(
    data = dat,
    y = "outcome",
    d = "D",
    id = c("id", "year"),
    kappa = c(-2, 1),
    refinement.method = "weightit",
    covs.formula = ~ lag(outcome, 1:2) + lag(x1, 1:2),
    refinement.args = list(method = "glm"),
    outcome.formula = ~ lag(outcome, 1:2) + lag(x1, 1:2),
    keep_data = TRUE
  ))
  
  meta <- attr(fit, "cbwsdid")
  expect_s3_class(fit, "fixest")
  expect_true(meta$dr)
  expect_true("y_raw" %in% names(meta$stacked.data))
  expect_true(nrow(meta$outcome.model.summary) > 0)
  
  expect_s3_class(cbwsdid_qoi(fit, type = "dynamic"), "tbl_df")
  expect_s3_class(cbwsdid_qoi(fit, type = "simple"), "tbl_df")
  expect_type(suppressWarnings(cbwsdid_qoi(fit, type = "cohort")), "list")
})


test_that("DR regression equals manually augmented within-stack contrasts", {
  skip_if_not_installed("fixest")
  data(cbwsdid_sim)
  dat <- subset(cbwsdid_sim, id <= 160)
  
  fit <- cbwsdid(
    data = dat,
    y = "outcome",
    d = "D",
    id = c("id", "year"),
    kappa = c(-2, 1),
    refinement.method = "none",
    outcome.formula = ~ lag(outcome, 1) + lag(x1, 1),
    keep_data = TRUE,
    pooled = FALSE
  )
  
  dyn <- cbwsdid_qoi(fit, type = "dynamic")
  meta <- attr(fit, "cbwsdid")
  sd0 <- meta$stacked.data
  outcome.spec <- meta$outcome.spec
  feature.cols <- outcome.spec$feature_name
  
  manual <- purrr::map_dfr(unique(sd0$subexp.id), function(subexp.i){
    subexp.data <- sd0 |>
      dplyr::filter(.data$subexp.id == subexp.i)
    feature.sample <- cbwsdid:::build_feature_sample_cbwsdid(
      subexp = dplyr::mutate(subexp.data, y = .data$y_raw),
      feature.spec = outcome.spec,
      include_treated = TRUE
    )
    base.y <- subexp.data |>
      dplyr::filter(.data$et == -1) |>
      dplyr::transmute(id = .data$id, subexp.id = .data$subexp.id, y_base = .data$y_raw)
    
    purrr::map_dfr(sort(unique(subexp.data$et[subexp.data$et != -1])), function(et.i){
      cell <- subexp.data |>
        dplyr::filter(.data$et == et.i) |>
        dplyr::left_join(base.y, by = c("id", "subexp.id")) |>
        dplyr::left_join(feature.sample, by = c("id", "subexp.id", "treated_sa")) |>
        dplyr::mutate(delta_y = .data$y_raw - .data$y_base)
      
      complete.rows <- stats::complete.cases(cell[, feature.cols, drop = FALSE]) &
        !is.na(cell$delta_y)
      controls <- cell[cell$treated_sa == 0 & complete.rows, , drop = FALSE]
      bt <- function(x) paste0("`", x, "`")
      outcome.formula <- stats::as.formula(
        paste("delta_y ~", paste(bt(feature.cols), collapse = " + "))
      )
      outcome.fit <- stats::lm(outcome.formula, data = controls)
      cell$mu_hat <- 0
      cell$mu_hat[complete.rows] <- stats::predict(
        outcome.fit,
        newdata = cell[complete.rows, , drop = FALSE]
      )
      cell <- cell |>
        dplyr::mutate(delta_dr = .data$delta_y - .data$mu_hat)
      
      contrast <- cell |>
        dplyr::group_by(.data$treated_sa) |>
        dplyr::summarise(
          mean_dr = stats::weighted.mean(.data$delta_dr, .data$Qsa),
          .groups = "drop"
        ) |>
        tidyr::pivot_wider(names_from = "treated_sa", values_from = "mean_dr",
                           names_prefix = "g")
      
      tibble::tibble(
        subexp.id = subexp.i,
        et = et.i,
        DID = contrast$g1 - contrast$g0
      )
    })
  })
  
  shares <- sd0 |>
    dplyr::distinct(.data$subexp.id, .data$Nda) |>
    dplyr::mutate(w = .data$Nda / sum(.data$Nda))
  
  agg <- manual |>
    dplyr::left_join(shares, by = "subexp.id") |>
    dplyr::group_by(.data$et) |>
    dplyr::summarise(theta = sum(.data$w * .data$DID), .groups = "drop")
  
  cmp <- dplyr::inner_join(
    dplyr::select(dyn, "et", beta = "estimate"),
    agg,
    by = "et"
  )
  expect_equal(cmp$beta, cmp$theta, tolerance = 1e-6)
})


test_that("DR bias correction works with nearest-neighbour matching", {
  skip_if_not_installed("fixest")
  skip_if_not_installed("MatchIt")
  data(cbwsdid_sim)
  dat <- subset(cbwsdid_sim, id <= 120)
  
  fit <- suppressWarnings(cbwsdid(
    data = dat,
    y = "outcome",
    d = "D",
    id = c("id", "year"),
    kappa = c(-2, 1),
    refinement.method = "matchit",
    covs.formula = ~ lag(outcome, 1) + lag(x1, 1),
    outcome.formula = ~ lag(outcome, 1) + lag(x1, 1),
    keep_data = TRUE,
    pooled = FALSE
  ))
  
  cohort <- cbwsdid_qoi(fit, type = "cohort")
  expect_true(attr(fit, "cbwsdid")$dr)
  expect_true(all(cohort$dynamic$n_treated > 0))
  expect_true("n_treated" %in% names(cohort$simple))
})


test_that("bootstrap reuses outcome.formula for DR fits", {
  skip_if_not_installed("fixest")
  data(cbwsdid_sim)
  dat <- subset(cbwsdid_sim, id <= 100)
  
  fit <- cbwsdid(
    data = dat,
    y = "outcome",
    d = "D",
    id = c("id", "year"),
    kappa = c(-1, 1),
    refinement.method = "none",
    outcome.formula = ~ lag(outcome, 1),
    keep_data = TRUE,
    pooled = FALSE
  )
  
  boot <- cbwsdid_bootstrap(
    fit,
    data = dat,
    B = 5,
    qoi = "dynamic",
    seed = 1,
    progress = FALSE
  )
  
  expect_true(attr(fit, "cbwsdid")$dr)
  expect_equal(attr(fit, "cbwsdid")$outcome.formula, ~ lag(outcome, 1))
  expect_s3_class(boot, "cbwsdid_bootstrap")
  expect_equal(boot$n_failed, 0L)
  expect_true(all(boot$summary$n_eff == boot$B_effective))
})

test_that("bare covariates in outcome.formula run, warn, and match lag(.,1)", {
  skip_if_not_installed("fixest")
  data(cbwsdid_sim)
  dat <- subset(cbwsdid_sim, id <= 160)

  # A bare term used to error in cbwsdid_outcome_residualize() because its
  # materialized column name collides with the raw same-named column. It must now
  # run, and emit the "undated" notice that bare terms are read at the base period.
  expect_warning(
    fit_bare <- cbwsdid(
      data = dat, y = "outcome", d = "D", id = c("id", "year"),
      kappa = c(-2, 1), refinement.method = "none",
      outcome.formula = ~ x1, keep_data = TRUE, pooled = FALSE
    ),
    "undated"
  )

  # `lag(x1, 1)` looks up x1 at the same base period (et = -1) as the bare term,
  # so the doubly robust dynamic estimates must be identical.
  fit_lag <- cbwsdid(
    data = dat, y = "outcome", d = "D", id = c("id", "year"),
    kappa = c(-2, 1), refinement.method = "none",
    outcome.formula = ~ lag(x1, 1), keep_data = TRUE, pooled = FALSE
  )

  expect_true(attr(fit_bare, "cbwsdid")$dr)
  expect_equal(
    cbwsdid_qoi(fit_bare, type = "dynamic"),
    cbwsdid_qoi(fit_lag, type = "dynamic"),
    tolerance = 1e-9
  )
})


# TODO: add a simulation demonstrating double robustness when either the
# refinement weights or the outcome model is correctly specified.
