moderation_dgp <- function(n_per_group = 35L,
                           rho = 0.5,
                           tau0 = 1,
                           psi = 0,
                           shift = 0,
                           seed = 123){
  set.seed(seed)
  years <- 2000:2006
  n_per_group <- as.integer(n_per_group)
  adoption <- rep(c(2003, 2005, Inf), each = n_per_group)
  id <- seq_along(adoption)
  m0 <- rep(seq(-1.5, 1.5, length.out = n_per_group), times = 3) +
    stats::rnorm(length(adoption), sd = 0.03)
  
  unit <- tibble::tibble(id = id, adoption = adoption, m0 = m0)
  dat <- tidyr::expand_grid(id = id, year = years) %>%
    dplyr::left_join(unit, by = "id") %>%
    dplyr::mutate(
      untreated = rho^(.data$year - min(years)) * .data$m0 +
        0.15 * (.data$year - min(years)),
      base_m = dplyr::if_else(
        is.finite(.data$adoption),
        rho^(.data$adoption - 1 - min(years)) * .data$m0 +
          0.15 * (.data$adoption - 1 - min(years)),
        NA_real_
      ),
      D = as.integer(is.finite(.data$adoption) & .data$year >= .data$adoption),
      treatment_effect = dplyr::if_else(
        .data$D == 1L,
        tau0 + psi * (.data$base_m - shift),
        0
      ),
      outcome = .data$untreated + .data$treatment_effect,
      outcome_shift = .data$outcome + 3,
      treated_x_foo = .data$outcome
    ) %>%
    dplyr::select("id", "year", "outcome", "outcome_shift", "treated_x_foo", "D")
  
  attr(dat, "truth") <- list(tau0 = tau0, psi = psi, shift = shift)
  dat
}


fit_moderation_dgp <- function(dat, moderation.formula = ~ lag(outcome, 1),
                               pooled = TRUE){
  suppressWarnings(cbwsdid(
    data = dat,
    y = "outcome",
    d = "D",
    id = c("id", "year"),
    kappa = c(-2, 1),
    refinement.method = "none",
    moderation.formula = moderation.formula,
    pooled = pooled,
    keep_data = TRUE
  ))
}


test_that("moderator main effects remove untreated trend-gradient confounding", {
  skip_if_not_installed("fixest")
  dat <- moderation_dgp(psi = 0)
  fit <- fit_moderation_dgp(dat, pooled = TRUE)
  
  slopes <- cbwsdid_qoi_moderation_slopes(fit) %>%
    dplyr::filter(.data$moderator == "outcome_l1", .data$et >= 0)
  
  expect_true(nrow(slopes) > 0)
  expect_lt(max(abs(slopes$estimate)), 1e-8)
})


test_that("moderated slopes and profile effects recover a true interaction", {
  skip_if_not_installed("fixest")
  dat <- moderation_dgp(tau0 = 1.2, psi = 0.7, shift = 0.25)
  fit <- fit_moderation_dgp(dat)
  truth <- attr(dat, "truth")
  
  slopes <- cbwsdid_qoi_moderation_slopes(fit) %>%
    dplyr::filter(.data$moderator == "outcome_l1", .data$et >= 0)
  expect_equal(slopes$estimate, rep(truth$psi, nrow(slopes)), tolerance = 1e-8)
  
  ref <- attr(fit, "cbwsdid")$moderation.reference
  profiles <- tibble::tibble(
    outcome_l1 = as.numeric(stats::quantile(ref$outcome_l1, c(0.25, 0.75)))
  )
  dynamic <- cbwsdid_qoi(
    fit,
    type = "dynamic",
    model_type = "moderated",
    newdata = profiles
  ) %>%
    dplyr::filter(.data$et >= 0)
  expected <- truth$tau0 + truth$psi * (dynamic$outcome_l1 - truth$shift)
  
  expect_equal(dynamic$estimate, expected, tolerance = 1e-8)
})


test_that("moderated effects are invariant to shifting the moderator", {
  skip_if_not_installed("fixest")
  dat <- moderation_dgp(tau0 = 0.8, psi = -0.4, shift = 0.1)
  fit <- fit_moderation_dgp(dat, moderation.formula = ~ lag(outcome, 1))
  fit_shift <- fit_moderation_dgp(dat, moderation.formula = ~ lag(outcome_shift, 1))
  
  slopes <- cbwsdid_qoi_moderation_slopes(fit) %>%
    dplyr::filter(.data$moderator == "outcome_l1") %>%
    dplyr::arrange(.data$et)
  slopes_shift <- cbwsdid_qoi_moderation_slopes(fit_shift) %>%
    dplyr::filter(.data$moderator == "outcome_shift_l1") %>%
    dplyr::arrange(.data$et)
  expect_equal(slopes$estimate, slopes_shift$estimate, tolerance = 1e-8)
  
  ref <- attr(fit, "cbwsdid")$moderation.reference
  profile <- tibble::tibble(outcome_l1 = mean(ref$outcome_l1))
  profile_shift <- tibble::tibble(outcome_shift_l1 = profile$outcome_l1 + 3)
  
  dynamic <- cbwsdid_qoi(
    fit,
    type = "dynamic",
    model_type = "moderated",
    newdata = profile
  )
  dynamic_shift <- cbwsdid_qoi(
    fit_shift,
    type = "dynamic",
    model_type = "moderated",
    newdata = profile_shift
  )
  
  expect_equal(dynamic$estimate, dynamic_shift$estimate, tolerance = 1e-8)
})


test_that("pre-period moderation slopes are placebo checks", {
  skip_if_not_installed("fixest")
  dat <- moderation_dgp(tau0 = 1.2, psi = 0.7, shift = 0.25)
  fit <- fit_moderation_dgp(dat)
  
  pre_slopes <- cbwsdid_qoi_moderation_slopes(fit) %>%
    dplyr::filter(.data$moderator == "outcome_l1", .data$et < 0)
  
  expect_true(nrow(pre_slopes) > 0)
  expect_lt(max(abs(pre_slopes$estimate)), 1e-8)
})


test_that("pooled moderated model removes untreated trend-gradient confounding", {
  skip_if_not_installed("fixest")
  dat <- moderation_dgp(psi = 0)
  fit <- fit_moderation_dgp(dat, pooled = TRUE)
  pooled_model <- attr(fit, "cbwsdid")$pooled_moderation_model
  pooled_coef <- stats::coef(pooled_model)
  
  expect_true("treated_post_x_outcome_l1" %in% names(pooled_coef))
  expect_lt(abs(unname(pooled_coef["treated_post_x_outcome_l1"])), 1e-8)
})


test_that("primary non-moderated path is unchanged when moderation is fitted", {
  skip_if_not_installed("fixest")
  dat <- moderation_dgp(psi = 0.5)
  plain <- suppressWarnings(cbwsdid(
    data = dat,
    y = "outcome",
    d = "D",
    id = c("id", "year"),
    kappa = c(-2, 1),
    refinement.method = "none",
    pooled = FALSE,
    keep_data = TRUE
  ))
  moderated <- fit_moderation_dgp(dat, pooled = FALSE)
  
  expect_equal(
    cbwsdid_qoi(plain, type = "dynamic"),
    cbwsdid_qoi(moderated, type = "dynamic"),
    tolerance = 1e-12
  )
})


test_that("moderator names cannot collide with treatment interaction names", {
  skip_if_not_installed("fixest")
  dat <- moderation_dgp(psi = 0.2)
  fit <- fit_moderation_dgp(dat, moderation.formula = ~ lag(treated_x_foo, 1))
  
  dynamic <- cbwsdid_qoi(
    fit,
    type = "dynamic",
    model_type = "moderated",
    newdata = tibble::tibble(treated_x_foo_l1 = 0)
  )
  coef_names <- names(stats::coef(attr(fit, "cbwsdid")$moderation_model))
  treated_terms <- coef_names[grepl("^et::(-?[0-9]+):treated_sa$", coef_names)]

  expect_s3_class(dynamic, "tbl_df")
  expect_equal(length(treated_terms), length(unique(dynamic$et)))

  # The slopes helper must return only true interaction terms
  # (`et::<h>:treated_x_treated_x_foo_l1`), never the main-effect terms
  # `et::<h>:treated_x_foo_l1` misread as slopes of a phantom moderator.
  slopes <- cbwsdid_qoi_moderation_slopes(fit)
  expect_equal(unique(slopes$moderator), "treated_x_foo_l1")
  slope_terms <- paste0("et::", slopes$et, ":treated_x_", slopes$moderator)
  expect_true(all(
    slope_terms %in% coef_names[grepl(":treated_x_treated_x_foo_l1$", coef_names)]
  ))
  post_slopes <- slopes$estimate[slopes$et >= 0]
  expect_equal(post_slopes, rep(0.2, length(post_slopes)), tolerance = 1e-8)
})


test_that("factor moderator levels remain available for explicit profiles", {
  skip_if_not_installed("fixest")
  dat <- moderation_dgp(n_per_group = 20, psi = 0.2)
  dat$control_only <- factor(
    ifelse(dat$id > 40L, "control_only", "treated_support")
  )
  fit <- suppressMessages(
    fit_moderation_dgp(dat, moderation.formula = ~ control_only)
  )
  
  expect_no_error(
    cbwsdid_qoi(
      fit,
      type = "dynamic",
      model_type = "moderated",
      newdata = tibble::tibble(control_only = "control_only")
    )
  )
})


moderation_pkg_fit <- function(n = 90L, pooled = TRUE){
  data(cbwsdid_sim)
  dat <- subset(cbwsdid_sim, id <= n)
  fit <- suppressMessages(suppressWarnings(cbwsdid(
    data = dat,
    y = "outcome",
    d = "D",
    id = c("id", "year"),
    kappa = c(-2, 1),
    refinement.method = "none",
    moderation.formula = ~ lag(outcome, 1),
    pooled = pooled,
    keep_data = TRUE
  )))
  
  list(fit = fit, data = dat)
}


test_that("moderated bootstrap returns the conditional schema and is reproducible", {
  skip_if_not_installed("fixest")
  fitted <- moderation_pkg_fit(n = 90, pooled = TRUE)
  
  conditional_dynamic <- cbwsdid_qoi(
    fitted$fit,
    type = "dynamic",
    model_type = "moderated"
  )
  boot_dynamic_1 <- cbwsdid_qoi(
    fitted$fit,
    type = "dynamic",
    model_type = "moderated",
    se.method = "bootstrap",
    B = 7,
    seed = 31
  )
  boot_dynamic_2 <- cbwsdid_qoi(
    fitted$fit,
    type = "dynamic",
    model_type = "moderated",
    se.method = "bootstrap",
    B = 7,
    seed = 31
  )
  boot_simple <- cbwsdid_qoi(
    fitted$fit,
    type = "simple",
    model_type = "moderated",
    se.method = "bootstrap",
    B = 7,
    seed = 32
  )
  boot_pooled <- cbwsdid_qoi(
    fitted$fit,
    type = "pooled",
    model_type = "moderated",
    se.method = "bootstrap",
    B = 7,
    seed = 33
  )
  
  expect_equal(names(boot_dynamic_1), names(conditional_dynamic))
  expect_equal(boot_dynamic_1$estimate, conditional_dynamic$estimate)
  expect_true(any(boot_dynamic_1$std.error > 0, na.rm = TRUE))
  expect_true(any(boot_simple$std.error > 0, na.rm = TRUE))
  expect_true(any(boot_pooled$std.error > 0, na.rm = TRUE))
  expect_equal(boot_dynamic_1, boot_dynamic_2)
  expect_s3_class(attr(boot_dynamic_1, "bootstrap"), "cbwsdid_bootstrap")
})


test_that("moderated bootstrap freezes the default profile", {
  skip_if_not_installed("fixest")
  fitted <- moderation_pkg_fit(n = 80, pooled = TRUE)
  default_profile <- mean(
    attr(fitted$fit, "cbwsdid")$moderation.reference$outcome_l1,
    na.rm = TRUE
  )
  
  boot_default <- cbwsdid_qoi(
    fitted$fit,
    type = "dynamic",
    model_type = "moderated",
    se.method = "bootstrap",
    B = 5,
    seed = 41
  )
  boot_explicit <- cbwsdid_qoi(
    fitted$fit,
    type = "dynamic",
    model_type = "moderated",
    newdata = tibble::tibble(outcome_l1 = default_profile),
    se.method = "bootstrap",
    B = 5,
    seed = 41
  )
  
  draws_default <- attr(boot_default, "bootstrap")$draws %>%
    dplyr::arrange(.data$rep, .data$et) %>%
    dplyr::pull(.data$estimate)
  draws_explicit <- attr(boot_explicit, "bootstrap")$draws %>%
    dplyr::arrange(.data$rep, .data$et) %>%
    dplyr::pull(.data$estimate)
  
  expect_equal(boot_default$estimate, boot_explicit$estimate)
  expect_equal(draws_default, draws_explicit)
})


test_that("moderated bootstrap ignores term-level NA draws", {
  skip_if_not_installed("fixest")
  data(cbwsdid_sim)
  dat <- subset(cbwsdid_sim, id <= 36)
  dat$rare_group <- factor(ifelse(dat$id %in% c(1L, 2L), "rare", "common"))
  
  fit <- suppressMessages(suppressWarnings(cbwsdid(
    data = dat,
    y = "outcome",
    d = "D",
    id = c("id", "year"),
    kappa = c(-2, 1),
    refinement.method = "none",
    moderation.formula = ~ rare_group,
    pooled = FALSE,
    keep_data = TRUE
  )))
  
  boot <- suppressMessages(cbwsdid_bootstrap(
    fit,
    data = dat,
    B = 15,
    qoi = "dynamic",
    model_type = "moderated",
    newdata = tibble::tibble(rare_group = "rare"),
    seed = 51,
    progress = FALSE
  ))
  
  expect_s3_class(boot, "cbwsdid_bootstrap")
  expect_true(all(boot$summary$n_eff <= boot$B_effective))
  expect_true(any(boot$summary$n_eff < boot$B_effective))
})


test_that("moderated bootstrap guard rails are informative", {
  skip_if_not_installed("fixest")
  fitted <- moderation_pkg_fit(n = 70, pooled = FALSE)
  
  expect_message(
    boot <- cbwsdid_bootstrap(
      fitted$fit,
      data = fitted$data,
      B = 3,
      model_type = "moderated",
      seed = 61,
      progress = FALSE
    ),
    "Dropping `qoi = \"pooled\"`"
  )
  expect_false("pooled" %in% boot$qoi)
  
  expect_error(
    cbwsdid_qoi(
      fitted$fit,
      type = "dynamic",
      model_type = "moderated",
      se.method = "analytic"
    ),
    "conditional.*bootstrap"
  )
})
