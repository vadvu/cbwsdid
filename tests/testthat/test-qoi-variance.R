qoi_variance_fit <- function(dr = FALSE,
                             pooled = FALSE,
                             kappa = c(-2, 2),
                             n = 180L){
  data(cbwsdid_sim)
  dat <- subset(cbwsdid_sim, id <= n)
  args <- list(
    data = dat,
    y = "outcome",
    d = "D",
    id = c("id", "year"),
    kappa = kappa,
    refinement.method = "weightit",
    covs.formula = ~ lag(x1, 1:2),
    refinement.args = list(method = "glm"),
    pooled = pooled
  )
  if(isTRUE(dr)){
    args$outcome.formula <- ~ lag(x1, 1:2)
  }

  list(
    fit = suppressWarnings(do.call(cbwsdid, args)),
    data = dat
  )
}


test_that("pooled estimate is one functional under every SE method", {
  skip_if_not_installed("fixest")
  fitted <- qoi_variance_fit(pooled = FALSE, n = 120L)

  conditional <- cbwsdid_qoi(
    fitted$fit,
    "pooled",
    se.method = "conditional"
  )
  analytic <- cbwsdid_qoi(
    fitted$fit,
    "pooled",
    se.method = "analytic"
  )
  bootstrap <- cbwsdid_qoi(
    fitted$fit,
    "pooled",
    se.method = "bootstrap",
    B = 5,
    seed = 11
  )

  expect_equal(conditional$estimate, analytic$estimate, tolerance = 1e-8)
  expect_equal(conditional$estimate, bootstrap$estimate, tolerance = 1e-8)
  expect_null(attr(fitted$fit, "cbwsdid")$pooled_model)
})


test_that("pooled equals the equal-weighted dynamic pre/post contrast", {
  skip_if_not_installed("fixest")
  fitted <- qoi_variance_fit(pooled = FALSE)
  dynamic <- cbwsdid_qoi(fitted$fit, "dynamic")
  post <- dynamic$et >= 0
  pre <- dynamic$et <= -2
  manual <- mean(dynamic$estimate[post]) - mean(dynamic$estimate[pre])

  pooled <- cbwsdid_qoi(fitted$fit, "pooled")
  expect_equal(pooled$estimate, manual, tolerance = 1e-8)
})


test_that("new pooled estimate is close to treated_post on balanced data", {
  skip_if_not_installed("fixest")
  fitted <- qoi_variance_fit(pooled = TRUE)
  new.pooled <- cbwsdid_qoi(fitted$fit, "pooled")
  old.pooled <- cbwsdid:::cbwsdid_qoi_linear_combo(
    attr(fitted$fit, "cbwsdid")$pooled_model,
    weights = stats::setNames(1, "treated_post")
  )

  expect_equal(new.pooled$estimate, old.pooled$estimate, tolerance = 0.10)
})


test_that("DR pooled analytic and conditional SEs agree to first order", {
  skip_if_not_installed("fixest")
  fitted <- qoi_variance_fit(dr = TRUE, pooled = FALSE)
  conditional <- cbwsdid_qoi(
    fitted$fit,
    "pooled",
    se.method = "conditional"
  )
  analytic <- cbwsdid_qoi(
    fitted$fit,
    "pooled",
    se.method = "analytic"
  )

  expect_equal(conditional$estimate, analytic$estimate, tolerance = 1e-8)
  expect_lt(
    abs(analytic$std.error / conditional$std.error - 1),
    0.15
  )
})


test_that("pooled requires a pre-treatment placebo period", {
  skip_if_not_installed("fixest")
  fitted <- qoi_variance_fit(
    pooled = FALSE,
    kappa = c(-1, 1),
    n = 120L
  )

  expect_error(
    cbwsdid_qoi(fitted$fit, "pooled"),
    "requires at least one pre-treatment period"
  )
})


test_that("cohort influence is mean zero within every cell", {
  skip_if_not_installed("fixest")
  fitted <- qoi_variance_fit(dr = TRUE)
  influence <- cbwsdid:::cbwsdid_influence_cohort(fitted$fit)
  cell.sums <- influence %>%
    dplyr::group_by(.data$subexp.id, .data$et) %>%
    dplyr::summarise(
      contribution = sum(.data$contribution),
      .groups = "drop"
    )

  expect_lt(max(abs(cell.sums$contribution)), 1e-8)
})


test_that("cohort analytic SEs are well formed and estimates are unchanged", {
  skip_if_not_installed("fixest")
  fitted <- qoi_variance_fit(dr = TRUE)
  conditional <- cbwsdid_qoi(
    fitted$fit,
    "cohort",
    se.method = "conditional"
  )
  expect_no_warning(
    analytic <- cbwsdid_qoi(
      fitted$fit,
      "cohort",
      se.method = "analytic"
    )
  )

  expect_equal(
    conditional$dynamic$estimate,
    analytic$dynamic$estimate,
    tolerance = 1e-12
  )
  expect_equal(
    conditional$simple$estimate,
    analytic$simple$estimate,
    tolerance = 1e-12
  )
  expect_true(all(is.finite(analytic$dynamic$std.error)))
  expect_true(all(is.finite(analytic$simple$std.error)))
  expect_true(all(analytic$dynamic$conf.low < analytic$dynamic$estimate))
  expect_true(all(analytic$dynamic$estimate < analytic$dynamic$conf.high))
  expect_true(all(analytic$simple$conf.low < analytic$simple$estimate))
  expect_true(all(analytic$simple$estimate < analytic$simple$conf.high))
})


test_that("dynamic and simple conditional paths remain unchanged", {
  skip_if_not_installed("fixest")
  fitted <- qoi_variance_fit(pooled = FALSE)

  expect_identical(
    cbwsdid_qoi(fitted$fit, "dynamic"),
    cbwsdid:::cbwsdid_qoi_dynamic_overall(fitted$fit)
  )
  expect_identical(
    cbwsdid_qoi(fitted$fit, "simple"),
    cbwsdid:::cbwsdid_qoi_simple_overall(fitted$fit)
  )
})


test_that("cohort analytic simple SE zero-fills units with partial post coverage", {
  skip_if_not_installed("fixest")
  data(cbwsdid_sim, envir = environment())
  dat <- subset(cbwsdid_sim, id <= 150)
  dat$outcome[dat$id == 120 & dat$year == 2005] <- NA

  fit <- suppressWarnings(cbwsdid(
    data = dat,
    y = "outcome",
    d = "D",
    id = c("id", "year"),
    kappa = c(-2, 1),
    refinement.method = "none",
    pooled = FALSE
  ))

  analytic <- cbwsdid_qoi(fit, "cohort", se.method = "analytic")
  influence <- attr(analytic$dynamic, "influence")
  post <- sort(unique(analytic$dynamic$et[analytic$dynamic$et >= 0]))

  presence <- influence %>%
    dplyr::filter(.data$et %in% post) %>%
    dplyr::count(.data$id, .data$subexp.key)
  expect_true(any(presence$n < length(post)))

  zero_fill <- influence %>%
    dplyr::filter(.data$et %in% post) %>%
    dplyr::group_by(.data$id, .data$subexp.key) %>%
    dplyr::summarise(
      psi = sum(.data$contribution) / length(post),
      .groups = "drop"
    ) %>%
    dplyr::group_by(.data$subexp.key) %>%
    dplyr::summarise(se = sqrt(sum(.data$psi^2)), .groups = "drop")

  drop_convention <- influence %>%
    dplyr::filter(.data$et %in% post) %>%
    dplyr::group_by(.data$id, .data$subexp.key) %>%
    dplyr::summarise(
      n_periods = dplyr::n_distinct(.data$et),
      psi = mean(.data$contribution),
      .groups = "drop"
    ) %>%
    dplyr::filter(.data$n_periods == length(post)) %>%
    dplyr::group_by(.data$subexp.key) %>%
    dplyr::summarise(se = sqrt(sum(.data$psi^2)), .groups = "drop")

  got <- analytic$simple$std.error[
    match(zero_fill$subexp.key, as.character(analytic$simple$subexp.id))
  ]
  expect_equal(got, zero_fill$se, tolerance = 1e-12)

  drop_se <- drop_convention$se[
    match(zero_fill$subexp.key, drop_convention$subexp.key)
  ]
  expect_false(isTRUE(all.equal(zero_fill$se, drop_se)))
})
