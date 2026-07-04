analytic_test_fit <- function(dr = FALSE,
                              refinement.method = "weightit",
                              n = 220L,
                              pooled = TRUE,
                              keep_data = TRUE){
  data(cbwsdid_sim)
  dat <- subset(cbwsdid_sim, id <= n)

  args <- list(
    data = dat,
    y = "outcome",
    d = "D",
    id = c("id", "year"),
    kappa = c(-2, 2),
    refinement.method = refinement.method,
    covs.formula = ~ lag(x1, 1:2),
    pooled = pooled,
    keep_data = keep_data
  )
  if(identical(refinement.method, "weightit")){
    args$refinement.args <- list(method = "glm")
  }
  if(isTRUE(dr)){
    args$outcome.formula <- ~ lag(x1, 1:2)
  }

  suppressWarnings(do.call(cbwsdid, args))
}


test_that("analytic influence contributions are mean zero", {
  skip_if_not_installed("fixest")

  weighting.fit <- analytic_test_fit(dr = FALSE)
  dr.fit <- analytic_test_fit(dr = TRUE)

  weighting.if <- cbwsdid:::cbwsdid_influence(weighting.fit)
  dr.if <- cbwsdid:::cbwsdid_influence(dr.fit)

  expect_lt(max(abs(colSums(weighting.if))), 1e-8)
  expect_lt(max(abs(colSums(dr.if))), 1e-8)
})


test_that("analytic and conditional SEs agree to first order", {
  skip_if_not_installed("fixest")

  for(fit in list(
    analytic_test_fit(dr = FALSE),
    analytic_test_fit(dr = TRUE)
  )){
    conditional <- cbwsdid_qoi(fit, "dynamic", se.method = "conditional")
    analytic <- cbwsdid_qoi(fit, "dynamic", se.method = "analytic")
    relative.gap <- abs(analytic$std.error / conditional$std.error - 1)

    expect_true(all(is.finite(analytic$std.error)))
    expect_lt(max(relative.gap), 0.15)
  }
})


test_that("analytic dynamic output has intervals and a joint covariance", {
  skip_if_not_installed("fixest")
  fit <- analytic_test_fit(dr = TRUE)

  out <- cbwsdid_qoi(fit, "dynamic", se.method = "analytic")
  variance <- attr(out, "vcov")

  expect_equal(length(out$std.error), nrow(out))
  expect_true(all(is.finite(out$std.error)))
  expect_true(all(out$conf.low < out$estimate))
  expect_true(all(out$estimate < out$conf.high))
  expect_equal(variance, t(variance), tolerance = 1e-12)
  expect_equal(nrow(variance), nrow(out))
  expect_gte(min(eigen(variance, symmetric = TRUE, only.values = TRUE)$values), -1e-10)
})


test_that("simple analytic SE is a linear combination of dynamic influence", {
  skip_if_not_installed("fixest")
  fit <- analytic_test_fit(dr = TRUE)

  influence <- cbwsdid:::cbwsdid_influence(fit)
  event.times <- as.numeric(colnames(influence))
  post.columns <- which(event.times >= 0)
  manual.influence <- rowMeans(influence[, post.columns, drop = FALSE])
  manual.se <- sqrt(sum(manual.influence^2))

  simple <- cbwsdid_qoi(fit, "simple", se.method = "analytic")
  expect_equal(simple$std.error, manual.se, tolerance = 1e-8)
})


test_that("pooled analytic inference uses the dynamic influence matrix", {
  skip_if_not_installed("fixest")
  fit <- analytic_test_fit(dr = FALSE, pooled = TRUE)

  influence <- cbwsdid:::cbwsdid_influence(fit)
  event.times <- as.numeric(colnames(influence))
  weights <- ifelse(
    event.times >= 0,
    1 / sum(event.times >= 0),
    -1 / sum(event.times <= -2)
  )
  manual.se <- sqrt(sum(as.numeric(influence %*% weights)^2))

  pooled <- cbwsdid_qoi(fit, "pooled", se.method = "analytic")
  expect_equal(pooled$std.error, manual.se, tolerance = 1e-8)
})


test_that("matching guard falls back only for pure matching", {
  skip_if_not_installed("fixest")
  skip_if_not_installed("MatchIt")

  matching.fit <- analytic_test_fit(
    dr = FALSE,
    refinement.method = "matchit",
    n = 140L,
    pooled = FALSE
  )
  conditional <- cbwsdid_qoi(matching.fit, "dynamic")
  expect_warning(
    analytic <- cbwsdid_qoi(
      matching.fit,
      "dynamic",
      se.method = "analytic"
    ),
    "not valid for pure fixed-k matching"
  )
  expect_equal(analytic, conditional)

  dr.matching.fit <- analytic_test_fit(
    dr = TRUE,
    refinement.method = "matchit",
    n = 140L,
    pooled = FALSE
  )
  expect_no_warning(
    dr.analytic <- cbwsdid_qoi(
      dr.matching.fit,
      "dynamic",
      se.method = "analytic"
    )
  )
  expect_true(all(is.finite(dr.analytic$std.error)))
})


test_that("bootstrap QOI auto-picks stored data and is reproducible", {
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
    pooled = FALSE
  )
  boot1 <- cbwsdid_qoi(
    fit,
    "dynamic",
    se.method = "bootstrap",
    B = 5,
    seed = 7
  )
  boot2 <- cbwsdid_qoi(
    fit,
    "dynamic",
    se.method = "bootstrap",
    B = 5,
    seed = 7
  )

  expect_identical(attr(fit, "cbwsdid")$data, dat)
  expect_equal(boot1, boot2)
  expect_s3_class(attr(boot1, "bootstrap"), "cbwsdid_bootstrap")

  lean.fit <- cbwsdid(
    data = dat,
    y = "outcome",
    d = "D",
    id = c("id", "year"),
    kappa = c(-1, 1),
    refinement.method = "none",
    pooled = FALSE,
    keep_data = FALSE
  )
  expect_false(is.null(attr(lean.fit, "cbwsdid")$if_components))
  expect_error(
    cbwsdid_qoi(
      lean.fit,
      "dynamic",
      se.method = "bootstrap",
      B = 2
    ),
    "Refit with `keep_data = TRUE`"
  )
})


test_that("conditional remains the no-op default", {
  skip_if_not_installed("fixest")
  fit <- analytic_test_fit(dr = FALSE, pooled = FALSE)

  expect_identical(
    cbwsdid_qoi(fit, "dynamic"),
    cbwsdid_qoi(fit, "dynamic", se.method = "conditional")
  )
})


# TODO: compare analytic SEs with Monte Carlo SD in simulation_plan.md,
# Exhibit C.
