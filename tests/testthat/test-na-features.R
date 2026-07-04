na_feature_data <- function(n_ids = 150L){
  data(cbwsdid_sim, envir = environment())
  dat <- subset(cbwsdid_sim, id <= n_ids)

  adoption <- stats::aggregate(
    year ~ id,
    data = subset(dat, D == 1),
    FUN = min
  )
  names(adoption) <- c("id", "adoption")

  list(
    dat = dat,
    treated_ids = adoption$id,
    adoption = adoption,
    control_ids = setdiff(unique(dat$id), adoption$id)
  )
}


test_that("treated units with missing refinement features hit the allow_treated_drop gate", {
  skip_if_not_installed("fixest")
  info <- na_feature_data()

  victim <- info$treated_ids[1]
  base_year <- info$adoption$adoption[info$adoption$id == victim] - 1

  dat_na <- info$dat
  dat_na$x1[dat_na$id == victim & dat_na$year == base_year] <- NA

  expect_error(
    suppressWarnings(cbwsdid(
      data = dat_na,
      y = "outcome",
      d = "D",
      id = c("id", "year"),
      kappa = c(-2, 1),
      refinement.method = "weightit",
      covs.formula = ~ lag(x1, 1),
      pooled = FALSE
    )),
    "missing \\(NA\\).*allow_treated_drop"
  )

  expect_warning(
    fit <- suppressMessages(cbwsdid(
      data = dat_na,
      y = "outcome",
      d = "D",
      id = c("id", "year"),
      kappa = c(-2, 1),
      refinement.method = "weightit",
      covs.formula = ~ lag(x1, 1),
      allow_treated_drop = TRUE,
      pooled = FALSE
    )),
    "no longer the full treated-cohort ATT"
  )

  design.weights <- attr(fit, "cbwsdid")$design.weights
  expect_false(victim %in% design.weights$id[design.weights$treated_sa == 1])
})


test_that("control units with missing refinement features are dropped with a message", {
  skip_if_not_installed("fixest")
  info <- na_feature_data()

  victim <- info$control_ids[1]
  dat_na <- info$dat
  dat_na$x1[dat_na$id == victim] <- NA

  expect_message(
    fit <- suppressWarnings(cbwsdid(
      data = dat_na,
      y = "outcome",
      d = "D",
      id = c("id", "year"),
      kappa = c(-2, 1),
      refinement.method = "weightit",
      covs.formula = ~ lag(x1, 1),
      pooled = FALSE
    )),
    "control unit\\(s\\) with missing refinement features"
  )

  design.weights <- attr(fit, "cbwsdid")$design.weights
  expect_false(victim %in% design.weights$id)
})


test_that("doubly robust fits warn about unresidualized rows from missing features", {
  skip_if_not_installed("fixest")
  info <- na_feature_data()

  victim <- info$treated_ids[1]
  base_year <- info$adoption$adoption[info$adoption$id == victim] - 1

  dat_na <- info$dat
  dat_na$x1[dat_na$id == victim & dat_na$year == base_year] <- NA

  expect_warning(
    fit <- cbwsdid(
      data = dat_na,
      y = "outcome",
      d = "D",
      id = c("id", "year"),
      kappa = c(-2, 1),
      refinement.method = "none",
      outcome.formula = ~ lag(x1, 1),
      pooled = FALSE
    ),
    "unresidualized"
  )

  stacked <- attr(fit, "cbwsdid")$stacked.data
  victim_rows <- stacked[stacked$id == victim & stacked$et != -1, ]
  expect_true(nrow(victim_rows) > 0)
  expect_equal(victim_rows$y, victim_rows$y_raw)
})
