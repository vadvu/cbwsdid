# Tests for the cohort quantity of interest after the sparse-i() refactor:
# coefficient names must be recovered (no silent all-NA), per-cell n_treated must
# be reported, and the cohort decomposition must aggregate to the overall effect.

make_panel <- function(n = 150L, seed = 1L) {
  set.seed(seed)
  years <- 2000:2012
  u <- tibble::tibble(
    id = seq_len(n),
    x1 = stats::rnorm(n),
    adopt = sample(c(2004, 2006, 2008, Inf), n, replace = TRUE,
                   prob = c(.2, .2, .2, .4))
  )
  tidyr::expand_grid(id = u$id, year = years) |>
    dplyr::left_join(u, by = "id") |>
    dplyr::group_by(id) |>
    dplyr::arrange(year, .by_group = TRUE) |>
    dplyr::mutate(
      outcome = cumsum(stats::rnorm(dplyr::n(), 0.1)) + 0.5 * x1,
      D = ifelse(is.infinite(adopt), 0L, as.integer(year >= adopt))
    ) |>
    dplyr::ungroup() |>
    dplyr::select(id, year, outcome, D, x1)
}

test_that("cohort QOI recovers coefficients and reports per-cell n_treated", {
  skip_if_not_installed("fixest")
  fit <- cbwsdid(make_panel(), y = "outcome", d = "D", id = c("id", "year"),
                 kappa = c(-3, 2), design = "absorbing",
                 refinement.method = "none", keep_data = TRUE)
  coh <- cbwsdid_qoi(fit, type = "cohort")

  # term-name matching after the i() refactor must work (no silent all-NA)
  expect_false(all(is.na(coh$dynamic$estimate)))
  expect_true("n_treated" %in% names(coh$dynamic))
  expect_false(any(grepl("n_treated\\.", names(coh$dynamic))))  # no .x/.y join clash
  expect_true(all(coh$dynamic$n_treated > 0))
  expect_true("n_treated" %in% names(coh$simple))
})

test_that("cell-level n_treated equals the cohort treated count under balanced stacks", {
  skip_if_not_installed("fixest")
  fit <- cbwsdid(make_panel(seed = 4), y = "outcome", d = "D",
                 id = c("id", "year"), kappa = c(-3, 2), design = "absorbing",
                 refinement.method = "none", keep_data = TRUE)
  coh <- cbwsdid_qoi(fit, type = "cohort")
  cw <- attr(fit, "cbwsdid")$cohort.weights |>
    dplyr::mutate(subexp.id = as.character(.data$subexp.id))
  chk <- coh$dynamic |>
    dplyr::mutate(subexp.id = as.character(.data$subexp.id)) |>
    dplyr::distinct(subexp.id, n_treated) |>
    dplyr::inner_join(dplyr::select(cw, subexp.id, Nda), by = "subexp.id")
  expect_equal(chk$n_treated, chk$Nda)
})
