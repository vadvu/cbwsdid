# Tests for the stacked-weight construction and the algebraic identity between the
# weighted stacked regression and the treated-share average of within-stack DID
# contrasts (Proposition 1, step 3).

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

build_stack <- function(dat, kappa = c(-3, 2)) {
  va <- stats::setNames(c("id", "time", "y", "d"),
                        c("id", "year", "outcome", "D"))
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
  res <- lapply(cands, function(a) {
    s <- cbwsdid:::subexp.row.construct(db, a, kappa)
    cbwsdid:::subexp.refine(s, var_aliases = va, kappa = kappa,
                            refinement.method = "none")
  })
  purrr::list_rbind(purrr::map(res, "subexp.refined"))
}

test_that("treated stacked weights are 1 and control mass matches the Wing identity", {
  subexp.all <- build_stack(make_panel())
  W <- cbwsdid:::compute_stack_weights(subexp.all)
  dw <- W$design.weights

  expect_true(all(abs(dw$stack_weight[dw$treated_sa == 1] - 1) < 1e-12))

  NtC <- sum(dw$bsa[dw$treated_sa == 0])
  NdT <- dw |>
    dplyr::distinct(subexp.id, Nda) |>
    dplyr::summarise(s = sum(Nda)) |>
    dplyr::pull(s)

  chk <- dw |>
    dplyr::filter(treated_sa == 0) |>
    dplyr::group_by(subexp.id) |>
    dplyr::summarise(sumW = sum(stack_weight), Nda = dplyr::first(Nda),
                     .groups = "drop") |>
    dplyr::mutate(target = NtC * Nda / NdT)

  expect_equal(chk$sumW, chk$target, tolerance = 1e-9)
})

test_that("the weight identity holds for arbitrary within-stack control weights", {
  subexp.all <- build_stack(make_panel(seed = 11))
  set.seed(99)
  bsa_tbl <- subexp.all |>
    dplyr::distinct(id, subexp.id, treated_sa) |>
    dplyr::mutate(bsa = ifelse(treated_sa == 1, 1, stats::runif(dplyr::n(), 0.1, 3)))
  subexp.all2 <- subexp.all |>
    dplyr::select(-bsa) |>
    dplyr::left_join(bsa_tbl, by = c("id", "subexp.id", "treated_sa"))

  W <- cbwsdid:::compute_stack_weights(subexp.all2)
  dw <- W$design.weights
  NtC <- sum(dw$bsa[dw$treated_sa == 0])
  NdT <- dw |> dplyr::distinct(subexp.id, Nda) |>
    dplyr::summarise(s = sum(Nda)) |> dplyr::pull(s)
  chk <- dw |>
    dplyr::filter(treated_sa == 0) |>
    dplyr::group_by(subexp.id) |>
    dplyr::summarise(sumW = sum(stack_weight), Nda = dplyr::first(Nda),
                     .groups = "drop") |>
    dplyr::mutate(target = NtC * Nda / NdT)
  expect_equal(chk$sumW, chk$target, tolerance = 1e-9)
})

test_that("weighted stacked regression equals the treated-share average of within-stack DIDs", {
  skip_if_not_installed("fixest")
  dat <- make_panel(seed = 7)
  kappa <- c(-3, 2)
  fit <- cbwsdid(dat, y = "outcome", d = "D", id = c("id", "year"),
                 kappa = kappa, design = "absorbing",
                 refinement.method = "none", keep_data = TRUE)
  dyn <- cbwsdid_qoi(fit, type = "dynamic")

  sd0 <- attr(fit, "cbwsdid")$stacked.data
  ref <- sd0 |> dplyr::filter(et == -1) |>
    dplyr::select(id, subexp.id, treated_sa, yref = y)
  mm <- sd0 |> dplyr::left_join(ref, by = c("id", "subexp.id", "treated_sa")) |>
    dplyr::mutate(dY = y - yref)
  shares <- mm |> dplyr::distinct(subexp.id, Nda) |>
    dplyr::mutate(w = Nda / sum(Nda))
  manual <- mm |> dplyr::filter(et != -1) |>
    dplyr::group_by(subexp.id, et, treated_sa) |>
    dplyr::summarise(m = stats::weighted.mean(dY, Qsa), .groups = "drop") |>
    tidyr::pivot_wider(names_from = treated_sa, values_from = m, names_prefix = "g") |>
    dplyr::mutate(DID = g1 - g0) |>
    dplyr::left_join(shares, by = "subexp.id") |>
    dplyr::group_by(et) |>
    dplyr::summarise(theta = sum(w * DID), .groups = "drop")

  cmp <- dplyr::inner_join(
    dplyr::select(dyn, et, beta = estimate), manual, by = "et"
  )
  expect_equal(cmp$beta, cmp$theta, tolerance = 1e-6)
})

test_that("cohort decomposition aggregates back to the overall dynamic effect", {
  skip_if_not_installed("fixest")
  dat <- make_panel(seed = 3)
  fit <- cbwsdid(dat, y = "outcome", d = "D", id = c("id", "year"),
                 kappa = c(-3, 2), design = "absorbing",
                 refinement.method = "none", keep_data = TRUE)
  dyn <- cbwsdid_qoi(fit, type = "dynamic")
  coh <- cbwsdid_qoi(fit, type = "cohort")
  agg <- coh$dynamic |>
    dplyr::group_by(et) |>
    dplyr::summarise(agg = sum(cohort_weight * estimate), .groups = "drop")
  cmp <- dplyr::inner_join(dplyr::select(dyn, et, beta = estimate), agg, by = "et")
  expect_equal(cmp$beta, cmp$agg, tolerance = 1e-6)
})
