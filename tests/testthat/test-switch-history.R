switch_history_panel <- function(){
  years <- 1:10
  d_paths <- list(
    t1 = c(0, 0, 0, 0, 0, 1, 1, 0, 0, 0),  # switches on at year 6
    c1 = c(0, 0, 1, 0, 0, 0, 0, 0, 0, 0),  # treated at year 3 = tau - 3
    c2 = rep(0, 10),
    c3 = rep(0, 10)
  )

  purrr::map_dfr(names(d_paths), function(u){
    tibble::tibble(
      id = u,
      year = years,
      d = d_paths[[u]],
      y = 0.2 * years + (u == "t1") * 0.5 * (years >= 6) + as.numeric(factor(u))
    )
  })
}


switch_history_fit <- function(dat, ...){
  cbwsdid(
    data = dat,
    y = "y",
    d = "d",
    id = c("id", "year"),
    kappa = c(-3, 1),
    design = "switch01",
    refinement.method = "none",
    pooled = FALSE,
    ...
  )
}


test_that("history.length defaults to abs(kappa[1]) and excludes in-window-treated controls", {
  skip_if_not_installed("fixest")
  dat <- switch_history_panel()

  fit <- suppressWarnings(switch_history_fit(dat))
  meta <- attr(fit, "cbwsdid")

  expect_identical(meta$history.length, 3L)
  controls <- meta$design.weights$id[meta$design.weights$treated_sa == 0]
  expect_false("c1" %in% controls)
  expect_true(all(c("c2", "c3") %in% controls))
})


test_that("history.length below abs(kappa[1]) warns and admits window-contaminated controls", {
  skip_if_not_installed("fixest")
  dat <- switch_history_panel()

  warnings_seen <- capture_warnings(
    fit <- switch_history_fit(dat, history.length = 1)
  )
  expect_true(any(grepl("history-matched only over the most recent", warnings_seen)))

  controls <- attr(fit, "cbwsdid")$design.weights
  expect_true("c1" %in% controls$id[controls$treated_sa == 0])
})


test_that("history.length above abs(kappa[1]) and invalid kappa error early and informatively", {
  skip_if_not_installed("fixest")
  dat <- switch_history_panel()

  expect_error(
    switch_history_fit(dat, history.length = 4),
    "cannot exceed"
  )

  expect_error(
    cbwsdid(
      data = dat,
      y = "y",
      d = "d",
      id = c("id", "year"),
      kappa = c(-3, -1),
      design = "switch01",
      refinement.method = "none",
      pooled = FALSE
    ),
    "kappa_post"
  )
})
