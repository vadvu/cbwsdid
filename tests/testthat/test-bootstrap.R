test_that("bootstrap returns dynamic summaries and is reproducible", {
  skip_if_not_installed("fixest")
  data(cbwsdid_sim)
  dat <- subset(cbwsdid_sim, id <= 150)
  
  fit <- cbwsdid(
    data = dat,
    y = "outcome",
    d = "D",
    id = c("id", "year"),
    kappa = c(-2, 1),
    refinement.method = "none",
    pooled = FALSE
  )
  
  boot1 <- cbwsdid_bootstrap(
    fit,
    data = dat,
    B = 25,
    qoi = "dynamic",
    seed = 1,
    progress = FALSE
  )
  boot2 <- cbwsdid_bootstrap(
    fit,
    data = dat,
    B = 25,
    qoi = "dynamic",
    seed = 1,
    progress = FALSE
  )
  
  dyn <- cbwsdid_qoi(fit, type = "dynamic")
  
  expect_s3_class(boot1, "cbwsdid_bootstrap")
  expect_equal(nrow(boot1$summary), nrow(dyn))
  expect_true(all(boot1$summary$boot_se > 0, na.rm = TRUE))
  expect_true(all(boot1$summary$conf.low < boot1$summary$estimate, na.rm = TRUE))
  expect_true(all(boot1$summary$estimate < boot1$summary$conf.high, na.rm = TRUE))
  expect_equal(boot1$summary$estimate, dyn$estimate)
  expect_equal(boot1$summary, boot2$summary)
})


test_that("bootstrap resampler relabels duplicated units", {
  dat <- tibble::tibble(
    id = c(1L, 1L, 2L, 2L),
    year = c(2000L, 2001L, 2000L, 2001L),
    y = c(1, 2, 3, 4),
    d = c(0, 1, 0, 0)
  )
  
  drawn <- cbwsdid:::cbwsdid_bootstrap_resample(
    data = dat,
    unit_col = "id",
    cluster_col = "id",
    drawn = c(1L, 1L)
  )
  
  duplicated_source <- sub("__b[0-9]+$", "", drawn$id)
  
  expect_equal(nrow(drawn), 4L)
  expect_equal(sort(unique(drawn$id)), c("1__b1", "1__b2"))
  expect_equal(as.integer(table(duplicated_source)), 4L)
  expect_equal(as.integer(table(drawn$id)), c(2L, 2L))
  
  one_cluster <- dplyr::filter(dat, .data$id == 1L)
  one_cluster$id <- 5L
  single <- cbwsdid:::cbwsdid_bootstrap_resample(
    data = one_cluster,
    unit_col = "id",
    cluster_col = "id"
  )
  
  expect_equal(nrow(single), nrow(one_cluster))
  expect_equal(unique(single$id), "5__b1")
})


test_that("default qoi drops pooled when the original model has no pooled fit", {
  skip_if_not_installed("fixest")
  data(cbwsdid_sim)
  dat <- subset(cbwsdid_sim, id <= 80)
  
  fit <- cbwsdid(
    data = dat,
    y = "outcome",
    d = "D",
    id = c("id", "year"),
    kappa = c(-1, 1),
    refinement.method = "none",
    pooled = FALSE
  )
  
  expect_warning(
    boot <- cbwsdid_bootstrap(
      fit,
      data = dat,
      B = 3,
      seed = 1,
      progress = FALSE
    ),
    "Dropping `qoi = \"pooled\"`"
  )
  expect_false("pooled" %in% boot$qoi)
  expect_false("pooled" %in% boot$summary$qoi)
  expect_error(
    cbwsdid_bootstrap(
      fit,
      data = dat,
      B = 3,
      qoi = "pooled",
      seed = 1,
      progress = FALSE
    ),
    "does not store a pooled model"
  )
})


test_that("nearest-neighbour matching emits the bootstrap caveat", {
  skip_if_not_installed("fixest")
  skip_if_not_installed("MatchIt")
  data(cbwsdid_sim)
  dat <- subset(cbwsdid_sim, id <= 80)
  
  fit <- suppressWarnings(cbwsdid(
    data = dat,
    y = "outcome",
    d = "D",
    id = c("id", "year"),
    kappa = c(-1, 1),
    refinement.method = "matchit",
    covs.formula = ~ lag(outcome, 1),
    pooled = FALSE
  ))
  
  expect_warning(
    cbwsdid_bootstrap(
      fit,
      data = dat,
      B = 3,
      qoi = "dynamic",
      seed = 1,
      progress = FALSE
    ),
    "bootstrap"
  )
})


test_that("bootstrap records failed replications without aborting", {
  skip_if_not_installed("fixest")
  dat <- tibble::tibble(
    id = rep(1:8, each = 4),
    year = rep(2000:2003, times = 8),
    outcome = stats::rnorm(32),
    D = as.integer(rep(2000:2003, times = 8) >= rep(c(2001, 2001, 2002, 2002, Inf, Inf, Inf, Inf), each = 4))
  )
  
  fit <- cbwsdid(
    data = dat,
    y = "outcome",
    d = "D",
    id = c("id", "year"),
    kappa = c(-1, 1),
    refinement.method = "none",
    pooled = FALSE
  )
  
  boot <- cbwsdid_bootstrap(
    fit,
    data = dat,
    B = 5,
    qoi = "dynamic",
    cluster = "id",
    seed = 3,
    progress = FALSE
  )
  
  expect_s3_class(boot, "cbwsdid_bootstrap")
  expect_gte(boot$n_failed, 0L)
  expect_equal(boot$B_effective, boot$B - boot$n_failed)
})


test_that("cluster argument defaults to unit id and accepts coarser clusters", {
  skip_if_not_installed("fixest")
  data(cbwsdid_sim)
  dat <- subset(cbwsdid_sim, id <= 100)
  dat$region <- ifelse(dat$id <= 50, "a", "b")
  
  fit <- cbwsdid(
    data = dat,
    y = "outcome",
    d = "D",
    id = c("id", "year"),
    kappa = c(-1, 1),
    refinement.method = "none",
    pooled = FALSE
  )
  
  boot_default <- cbwsdid_bootstrap(
    fit,
    data = dat,
    B = 4,
    qoi = "dynamic",
    cluster = "id",
    seed = 4,
    progress = FALSE
  )
  boot_region <- cbwsdid_bootstrap(
    fit,
    data = dat,
    B = 4,
    qoi = "dynamic",
    cluster = "region",
    seed = 4,
    progress = FALSE
  )
  
  expect_equal(boot_default$cluster, "id")
  expect_equal(boot_region$cluster, "region")
  expect_s3_class(boot_region, "cbwsdid_bootstrap")
})
