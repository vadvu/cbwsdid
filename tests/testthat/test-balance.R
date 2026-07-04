test_that("Diff.Un is computed on the full design sample under matchit", {
  skip_if_not_installed("fixest")
  skip_if_not_installed("MatchIt")
  skip_if_not_installed("cobalt")
  data(cbwsdid_sim, envir = environment())
  dat <- subset(cbwsdid_sim, id <= 200)

  fit <- suppressWarnings(cbwsdid(
    data = dat,
    y = "outcome",
    d = "D",
    id = c("id", "year"),
    kappa = c(-2, 1),
    refinement.method = "matchit",
    covs.formula = ~ lag(x1, 1) + lag(outcome, 1),
    pooled = FALSE
  ))

  bal <- cbwsdid_balance(fit)
  detail <- attr(fit, "cbwsdid")$refinement.details[[1]]
  subexp.id <- as.character(detail$design.sample$subexp.id[1])

  design.sample <- detail$design.sample
  x <- design.sample$x1_l1
  treated <- design.sample$treated_sa == 1
  smd_un_hand <- (mean(x[treated]) - mean(x[!treated])) / stats::sd(x[treated])

  smd.un <- bal$by_subexp[[subexp.id]]$smd.un
  smd.adj <- bal$by_subexp[[subexp.id]]$smd.adj
  expect_equal(smd.un$x1[smd.un$et == -1], smd_un_hand, tolerance = 1e-6)

  # 1:1 matching without replacement must show adjustment: with unmatched
  # controls back in the unadjusted column, the two tables can no longer be
  # identical (the old bug made them coincide exactly).
  expect_false(isTRUE(all.equal(smd.un, smd.adj)))
})
