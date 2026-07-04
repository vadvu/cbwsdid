#' Pipeline Cluster Bootstrap for `cbwsdid`
#'
#' Re-runs the full [cbwsdid()] pipeline on cluster bootstrap samples and
#' summarizes bootstrap uncertainty for overall or moderated quantities of
#' interest.
#'
#' @param model A `fixest` model object returned by [cbwsdid()].
#' @param data Original raw panel dataframe passed to [cbwsdid()].
#' @param B Integer. Number of bootstrap replications.
#' @param cluster Optional character. Name of the resampling cluster column in
#'   `data`. Defaults to the unit identifier supplied to [cbwsdid()].
#' @param qoi Character vector. Quantities of interest to bootstrap. Supported
#'   values are `"dynamic"`, `"simple"`, and `"pooled"`.
#' @param post Optional integer vector of event times passed to
#'   `cbwsdid_qoi(type = "simple")`.
#' @param model_type Character. `"overall"` bootstraps the baseline QOIs.
#'   `"moderated"` bootstraps moderated QOIs at fixed moderator profiles.
#' @param newdata Optional data frame of moderator profiles used only when
#'   `model_type = "moderated"`. If `NULL`, the treated mean profile from the
#'   original model is frozen across bootstrap draws.
#' @param ci_type Confidence interval type. `"percentile"` uses bootstrap
#'   quantiles; `"normal"` uses the original estimate plus or minus a normal
#'   critical value times the bootstrap standard error.
#' @param level Confidence level.
#' @param parallel Logical. If `TRUE`, run bootstrap replications in parallel
#'   using `future` and `future.apply`.
#' @param workers Optional integer. Number of parallel workers. If `NULL`, a
#'   conservative default based on available cores is used.
#' @param seed Optional integer seed for reproducible bootstrap draws.
#' @param progress Logical. If `TRUE`, print lightweight progress messages.
#'
#' @details
#' This is a pipeline bootstrap: each replication rebuilds sub-experiments,
#' re-estimates refinement weights, recomputes stacked weights, refits the
#' second-stage model, and then extracts the requested overall quantity of
#' interest. The original `model` is used to recover estimation arguments,
#' point estimates, and, for moderated QOIs, fixed moderator profiles.
#'
#' The pooled quantity is recomputed in every replication as the equally
#' weighted post-treatment dynamic average minus the equally weighted
#' pre-treatment placebo average. It does not use the standalone
#' `treated_post` regression. For moderated QOIs, dynamic and simple
#' quantities are fixed linear combinations of the moderated dynamic model,
#' while `qoi = "pooled"` uses the standalone pooled moderated model and
#' therefore requires `pooled = TRUE` in the original fit.
#'
#' Moderated bootstrap profiles are resolved once from the original model and
#' kept fixed in every draw. This propagates uncertainty in the effect at those
#' evaluation points, not uncertainty in a re-estimated average profile.
#'
#' Duplicated clusters are relabelled with fresh unit identifiers before each
#' refit. This is required so that fixed effects and stacked-treatment logic
#' treat repeated bootstrap copies as distinct histories.
#'
#' The nonparametric bootstrap is not generally valid for nearest-neighbour
#' matching estimators. When `model` was fitted with nearest-neighbour
#' `MatchIt` refinement, this function warns but still runs.
#'
#' @return An object of class `"cbwsdid_bootstrap"`, a list with `summary`,
#'   `draws`, `B`, `B_effective`, `n_failed`, `cluster`, `qoi`, `ci_type`,
#'   `level`, and `call`.
#' @export
#'
#' @examples
#' data(cbwsdid_sim)
#' cbwsdid_sim_small <- subset(cbwsdid_sim, id <= 150)
#'
#' fit <- cbwsdid(
#'   data = cbwsdid_sim_small,
#'   y = "outcome",
#'   d = "D",
#'   id = c("id", "year"),
#'   kappa = c(-2, 1),
#'   refinement.method = "none"
#' )
#'
#' boot <- cbwsdid_bootstrap(
#'   fit,
#'   data = cbwsdid_sim_small,
#'   B = 19,
#'   seed = 1,
#'   progress = FALSE
#' )
#' boot
cbwsdid_bootstrap <- function(model,
                              data,
                              B = 999L,
                              cluster = NULL,
                              qoi = c("dynamic", "simple", "pooled"),
                              post = NULL,
                              model_type = c("overall", "moderated"),
                              newdata = NULL,
                              ci_type = c("percentile", "normal"),
                              level = 0.95,
                              parallel = FALSE,
                              workers = NULL,
                              seed = NULL,
                              progress = FALSE){
  call <- match.call()
  qoi_missing <- missing(qoi)
  qoi <- match.arg(qoi, several.ok = TRUE)
  model_type <- match.arg(model_type)
  ci_type <- match.arg(ci_type)
  meta <- cbwsdid_bootstrap_validate_inputs(
    model = model,
    data = data,
    B = B,
    cluster = cluster,
    level = level
  )
  B <- meta$B
  cluster_col <- meta$cluster
  meta <- meta$meta
  qoi <- cbwsdid_bootstrap_resolve_qoi(
    qoi = qoi,
    meta = meta,
    qoi_missing = qoi_missing,
    model_type = model_type
  )
  cluster_index <- cbwsdid_bootstrap_cluster_index(
    data = data,
    cluster_col = cluster_col
  )
  profiles <- NULL
  original_qoi <- NULL
  if(identical(model_type, "moderated")){
    profiles <- cbwsdid_prepare_moderator_profiles(meta = meta, newdata = newdata)
    if("simple" %in% qoi){
      post <- cbwsdid_bootstrap_moderated_resolve_post(
        model = model,
        post = post
      )
    }
    original_qoi <- cbwsdid_bootstrap_extract_qoi_moderated(
      model = model,
      qoi = qoi,
      post = post,
      level = level,
      profiles = profiles
    )
  } else {
    if(!is.null(newdata)){
      stop("`newdata` is used only when `model_type = \"moderated\"`.")
    }
    original_qoi <- cbwsdid_bootstrap_extract_qoi(
      model = model,
      qoi = qoi,
      post = post,
      level = level
    )
  }

  cbwsdid_bootstrap_warn_matching(meta)

  resolved.workers <- cbwsdid_bootstrap_resolve_workers(
    parallel = parallel,
    workers = workers
  )

  if(!is.null(seed)){
    seed <- as.integer(seed[1])
    if(is.na(seed)){
      stop("`seed` must be an integer or `NULL`.")
    }
    set.seed(seed)
  }

  if(isTRUE(progress)){
    message(
      "Running ", B, " bootstrap replication",
      if(B == 1L) "" else "s",
      " clustered by `", cluster_col, "`."
    )
  }

  worker <- function(b){
    cbwsdid_bootstrap_worker(
      b = b,
      data = data,
      meta = meta,
      cluster_col = cluster_col,
      cluster_index = cluster_index,
      qoi = qoi,
      post = post,
      model_type = model_type,
      profiles = profiles,
      level = level
    )
  }

  reps <- seq_len(B)
  if(isTRUE(parallel) && resolved.workers > 1L){
    old.plan <- future::plan()
    on.exit(future::plan(old.plan), add = TRUE)
    future::plan(future::multisession, workers = resolved.workers)
    results <- future.apply::future_lapply(
      reps,
      worker,
      future.seed = seed %||% TRUE
    )
  } else {
    progress_every <- max(1L, floor(B / 10L))
    results <- lapply(reps, function(b){
      if(isTRUE(progress) && (b == 1L || b == B || b %% progress_every == 0L)){
        message("Bootstrap replication ", b, "/", B)
      }
      worker(b)
    })
  }

  failed <- purrr::map_lgl(results, is.null)
  n_failed <- sum(failed)
  B_effective <- B - n_failed

  if(n_failed > 0L && n_failed / B > 0.10){
    warning(
      n_failed,
      " of ",
      B,
      " bootstrap replications failed. Interpret bootstrap summaries with care."
    )
  }

  draws <- cbwsdid_bootstrap_collect_draws(
    results = results,
    canonical = original_qoi
  )

  summary <- cbwsdid_bootstrap_summarise(
    original_qoi = original_qoi,
    draws = draws,
    ci_type = ci_type,
    level = level
  )

  out <- list(
    summary = summary,
    draws = draws,
    B = B,
    B_effective = B_effective,
    n_failed = n_failed,
    cluster = cluster_col,
    qoi = qoi,
    model_type = model_type,
    ci_type = ci_type,
    level = level,
    call = call
  )
  class(out) <- "cbwsdid_bootstrap"
  out
}


cbwsdid_bootstrap_validate_inputs <- function(model,
                                              data,
                                              B,
                                              cluster = NULL,
                                              level = 0.95){
  if(!inherits(model, "fixest")){
    stop("`model` must be a fixest object returned by `cbwsdid()`.")
  }

  meta <- attr(model, "cbwsdid")
  if(is.null(meta)){
    stop("`model` does not contain `cbwsdid` metadata.")
  }

  required_meta <- c("y", "d", "id")
  missing_meta <- required_meta[!purrr::map_lgl(required_meta, ~ !is.null(meta[[.x]]))]
  if(length(missing_meta) > 0){
    stop(
      "`model` metadata is missing required field(s): ",
      paste(missing_meta, collapse = ", "),
      ". Refit the model with the current version of `cbwsdid()`."
    )
  }

  if(length(meta$id) != 2L){
    stop("`model` metadata field `id` must contain unit and time column names.")
  }

  if(!is.data.frame(data)){
    stop("`data` must be the original raw panel dataframe passed to `cbwsdid()`.")
  }

  B <- as.integer(B[1])
  if(is.na(B) || B < 1L){
    stop("`B` must be a positive integer.")
  }

  if(!is.numeric(level) || length(level) != 1L || is.na(level) ||
     level <= 0 || level >= 1){
    stop("`level` must be a single number between 0 and 1.")
  }

  unit_col <- meta$id[1]
  cluster_col <- cluster %||% unit_col
  if(!is.character(cluster_col) || length(cluster_col) != 1L || is.na(cluster_col)){
    stop("`cluster` must be a single column name or `NULL`.")
  }

  required_vars <- cbwsdid_bootstrap_required_vars(meta, cluster_col)
  missing_vars <- setdiff(required_vars, names(data))
  if(length(missing_vars) > 0){
    stop(
      "`data` is missing required column(s): ",
      paste(missing_vars, collapse = ", "),
      ". Pass the original raw panel used to fit `model`."
    )
  }

  if(anyNA(data[[cluster_col]])){
    stop("`cluster` column `", cluster_col, "` contains missing values.")
  }

  list(
    meta = meta,
    B = B,
    cluster = cluster_col
  )
}


cbwsdid_bootstrap_required_vars <- function(meta, cluster_col){
  alias_map <- stats::setNames(
    c(meta$id[1], meta$id[2], meta$y, meta$d),
    c("id", "time", "y", "d")
  )

  formula_vars <- unique(c(
    all.vars(meta$covs.formula),
    all.vars(meta$exact.formula),
    all.vars(meta$moderation.formula),
    all.vars(meta$outcome.formula)
  ))

  formula_vars <- purrr::map_chr(formula_vars, function(x){
    if(x %in% names(alias_map)){
      return(unname(alias_map[[x]]))
    }
    x
  })

  unique(c(meta$id, meta$y, meta$d, cluster_col, formula_vars))
}


cbwsdid_bootstrap_resolve_qoi <- function(qoi,
                                          meta,
                                          qoi_missing = FALSE,
                                          model_type = "overall"){
  if(!"pooled" %in% qoi){
    return(qoi)
  }

  if(identical(model_type, "moderated") &&
     is.null(meta$pooled_moderation_model)){
    message(
      "Dropping `qoi = \"pooled\"` because no moderated pooled model is ",
      "stored in this `cbwsdid` object. Re-estimate with `pooled = TRUE` and ",
      "`moderation.formula` to bootstrap moderated pooled effects."
    )
    qoi <- setdiff(qoi, "pooled")
    if(length(qoi) == 0L){
      stop("No bootstrap QOI remain after dropping unavailable moderated pooled effects.")
    }
    return(qoi)
  }

  has.pre.placebo <- !is.null(meta$kappa) && meta$kappa[1] <= -2
  if(!has.pre.placebo && isTRUE(qoi_missing)){
    warning(
      "Dropping `qoi = \"pooled\"` because the event window has no ",
      "pre-treatment placebo period at or before event time -2."
    )
    return(setdiff(qoi, "pooled"))
  }

  qoi
}


cbwsdid_bootstrap_warn_matching <- function(meta){
  method <- meta$refinement.args$method
  nearest <- is.null(method) || identical(method, "nearest")

  if(identical(meta$refinement.method, "matchit") && nearest){
    warning(
      "The nonparametric bootstrap is not generally valid for nearest-neighbour ",
      "matching (Abadie & Imbens 2008). Prefer smooth balancing weights ",
      "(`refinement.method = \"weightit\"` with entropy balancing or CBPS), ",
      "or use bias-corrected matching, when bootstrap inference is the goal."
    )
  }
}


cbwsdid_bootstrap_resolve_workers <- function(parallel = FALSE, workers = NULL){
  resolved.workers <- 1L
  if(!isTRUE(parallel)){
    return(resolved.workers)
  }

  if(!requireNamespace("future", quietly = TRUE) ||
     !requireNamespace("future.apply", quietly = TRUE)){
    stop(
      "Parallel execution requires the optional packages `future` and ",
      "`future.apply`. Please install them or set `parallel = FALSE`."
    )
  }

  if(is.null(workers)){
    detected.cores <- parallel::detectCores()
    if(is.na(detected.cores) || detected.cores < 2L){
      resolved.workers <- 1L
    } else {
      resolved.workers <- max(1L, as.integer(detected.cores - 1L))
    }
  } else {
    resolved.workers <- as.integer(workers[1])
    if(is.na(resolved.workers) || resolved.workers < 1L){
      stop("`workers` must be a positive integer or `NULL`.")
    }
  }

  resolved.workers
}


cbwsdid_bootstrap_worker <- function(b,
                                     data,
                                     meta,
                                     cluster_col,
                                     cluster_index,
                                     qoi,
                                     post,
                                     model_type = "overall",
                                     profiles = NULL,
                                     level = 0.95){
  tryCatch(
    {
      boot_data <- cbwsdid_bootstrap_resample(
        data = data,
        unit_col = meta$id[1],
        cluster_col = cluster_col,
        cluster_index = cluster_index
      )

      args <- list(
        data = boot_data,
        y = meta$y,
        d = meta$d,
        id = meta$id,
        kappa = meta$kappa,
        design = meta$design,
        post_path = meta$post_path,
        history.length = meta$history.length,
        first_switch_only = meta$first_switch_only,
        refinement.method = meta$refinement.method,
        covs.formula = meta$covs.formula,
        exact.formula = meta$exact.formula,
        refinement.args = meta$refinement.args,
        allow_treated_drop = meta$allow_treated_drop,
        moderation.formula = if(identical(model_type, "moderated")) meta$moderation.formula else NULL,
        outcome.formula = meta$outcome.formula,
        pooled = identical(model_type, "moderated") && "pooled" %in% qoi,
        keep_data = FALSE,
        parallel = FALSE
      )

      boot_fit <- do.call(cbwsdid, args)
      if(identical(model_type, "moderated")){
        return(cbwsdid_bootstrap_extract_qoi_moderated(
          model = boot_fit,
          qoi = qoi,
          post = post,
          level = level,
          profiles = profiles,
          strict = TRUE
        ))
      }

      cbwsdid_bootstrap_extract_qoi(
        model = boot_fit,
        qoi = qoi,
        post = post,
        level = level
      )
    },
    error = function(e){
      NULL
    }
  )
}


cbwsdid_bootstrap_resample <- function(data,
                                       unit_col,
                                       cluster_col,
                                       cluster_index = NULL,
                                       drawn = NULL){
  if(is.null(cluster_index)){
    cluster_index <- cbwsdid_bootstrap_cluster_index(
      data = data,
      cluster_col = cluster_col
    )
  }

  clusters <- cluster_index$clusters
  cluster_rows <- cluster_index$rows

  if(is.null(drawn)){
    drawn_pos <- sample.int(length(clusters), length(clusters), replace = TRUE)
  } else {
    drawn_pos <- match(drawn, clusters)
    if(anyNA(drawn_pos)){
      stop("`drawn` contains cluster values not present in `data`.")
    }
  }

  purrr::imap_dfr(drawn_pos, function(cluster_pos, k){
    rows <- data[cluster_rows[[as.character(cluster_pos)]], , drop = FALSE]
    rows[[unit_col]] <- paste0(rows[[unit_col]], "__b", k)
    rows
  })
}


cbwsdid_bootstrap_cluster_index <- function(data, cluster_col){
  clusters <- unique(data[[cluster_col]])
  if(length(clusters) == 0L){
    stop("`data` contains no bootstrap clusters.")
  }

  cluster_id <- match(data[[cluster_col]], clusters)
  rows <- split(seq_len(nrow(data)), cluster_id)

  list(
    clusters = clusters,
    rows = rows
  )
}


cbwsdid_bootstrap_extract_qoi <- function(model,
                                          qoi,
                                          post = NULL,
                                          level = 0.95){
  purrr::map_dfr(qoi, function(qoi_i){
    out_i <- cbwsdid_qoi(
      model = model,
      type = qoi_i,
      model_type = "overall",
      post = post,
      level = level
    )

    if(identical(qoi_i, "dynamic")){
      return(
        out_i %>%
          dplyr::transmute(
            qoi = .data$type,
            term = paste0("et=", .data$et),
            et = .data$et,
            estimate = .data$estimate
          )
      )
    }

    out_i %>%
      dplyr::transmute(
        qoi = .data$type,
        term = .data$type,
        et = NA_real_,
        estimate = .data$estimate
      )
  })
}


cbwsdid_bootstrap_extract_qoi_moderated <- function(model,
                                                    qoi,
                                                    post = NULL,
                                                    level = 0.95,
                                                    profiles = NULL,
                                                    strict = FALSE){
  meta <- attr(model, "cbwsdid")
  if(is.null(meta) || is.null(meta$moderation.formula) ||
     nrow(meta$moderation.spec) == 0L){
    stop("`model` does not contain moderated second-stage metadata.")
  }
  if(is.null(profiles)){
    profiles <- cbwsdid_prepare_moderator_profiles(meta = meta, newdata = NULL)
  }

  purrr::map_dfr(qoi, function(qoi_i){
    if(identical(qoi_i, "dynamic")){
      moderation_model <- meta$moderation_model
      if(is.null(moderation_model)){
        stop("No moderated dynamic model is stored in this `cbwsdid` object.")
      }
      dynamic_ets <- cbwsdid_bootstrap_moderated_event_times(moderation_model)

      return(purrr::map_dfr(seq_len(nrow(profiles$expanded)), function(i){
        profile_i <- profiles$expanded[i, , drop = FALSE]
        profile_label <- profiles$profile$profile[i]

        rows_i <- purrr::map_dfr(dynamic_ets, function(et_i){
          combo_i <- tryCatch(
            cbwsdid_qoi_moderated_dynamic_combo(
              model = moderation_model,
              et = et_i,
              moderator_values = profile_i,
              level = level,
              strict = strict
            ),
            error = function(e) cbwsdid_bootstrap_na_combo()
          )

          tibble::tibble(
            qoi = "dynamic",
            term = paste0("dynamic:", profile_label, ":et=", et_i),
            et = et_i,
            estimate = combo_i$estimate
          )
        })

        dplyr::bind_cols(profiles$profile[i, , drop = FALSE], rows_i)
      }))
    }

    if(identical(qoi_i, "simple")){
      moderation_model <- meta$moderation_model
      if(is.null(moderation_model)){
        stop("No moderated dynamic model is stored in this `cbwsdid` object.")
      }
      post_i <- cbwsdid_bootstrap_moderated_resolve_post(model, post = post)

      return(purrr::map_dfr(seq_len(nrow(profiles$expanded)), function(i){
        profile_i <- profiles$expanded[i, , drop = FALSE]
        profile_label <- profiles$profile$profile[i]
        combo_i <- tryCatch(
          cbwsdid_qoi_moderated_simple_combo(
            model = moderation_model,
            post = post_i,
            moderator_values = profile_i,
            level = level,
            strict = strict
          ),
          error = function(e) cbwsdid_bootstrap_na_combo()
        )

        dplyr::bind_cols(
          profiles$profile[i, , drop = FALSE],
          tibble::tibble(
            qoi = "simple",
            term = paste0("simple:", profile_label),
            et = NA_real_,
            estimate = combo_i$estimate
          )
        )
      }))
    }

    pooled_moderation_model <- meta$pooled_moderation_model
    purrr::map_dfr(seq_len(nrow(profiles$expanded)), function(i){
      profile_i <- profiles$expanded[i, , drop = FALSE]
      profile_label <- profiles$profile$profile[i]
      combo_i <- if(is.null(pooled_moderation_model)){
        cbwsdid_bootstrap_na_combo()
      } else {
        tryCatch(
          cbwsdid_qoi_moderated_pooled_combo(
            model = pooled_moderation_model,
            moderator_values = profile_i,
            level = level,
            strict = strict
          ),
          error = function(e) cbwsdid_bootstrap_na_combo()
        )
      }

      dplyr::bind_cols(
        profiles$profile[i, , drop = FALSE],
        tibble::tibble(
          qoi = "pooled",
          term = paste0("pooled:", profile_label),
          et = NA_real_,
          estimate = combo_i$estimate
        )
      )
    })
  })
}


cbwsdid_bootstrap_moderated_event_times <- function(model){
  coef_names <- names(stats::coef(model))
  term_pattern <- "^et::(-?[0-9]+):treated_sa$"
  terms <- coef_names[grepl(term_pattern, coef_names)]
  if(length(terms) == 0L){
    stop("No moderated dynamic treatment coefficients were found.")
  }

  sort(as.numeric(stringr::str_match(terms, term_pattern)[, 2]))
}


cbwsdid_bootstrap_moderated_resolve_post <- function(model, post = NULL){
  meta <- attr(model, "cbwsdid")
  moderation_model <- meta$moderation_model
  if(is.null(moderation_model)){
    stop("No moderated dynamic model is stored in this `cbwsdid` object.")
  }

  dynamic_ets <- cbwsdid_bootstrap_moderated_event_times(moderation_model)
  if(is.null(post)){
    post <- dynamic_ets[dynamic_ets >= 0]
  }

  post <- sort(unique(post))
  if(length(post) == 0L){
    stop("No post-treatment event times were selected.")
  }

  missing_post <- setdiff(post, dynamic_ets)
  if(length(missing_post) > 0L){
    stop(
      "Some requested event times are not present in the model: ",
      paste(missing_post, collapse = ", ")
    )
  }

  post
}


cbwsdid_bootstrap_na_combo <- function(){
  list(
    estimate = NA_real_,
    std.error = NA_real_,
    conf.low = NA_real_,
    conf.high = NA_real_
  )
}


cbwsdid_bootstrap_collect_draws <- function(results, canonical){
  successful_reps <- which(!purrr::map_lgl(results, is.null))
  key_cols <- intersect(c("qoi", "term", "profile", "et"), names(canonical))

  if(length(successful_reps) == 0L){
    empty <- canonical[0, key_cols, drop = FALSE]
    empty$rep <- integer()
    empty$estimate <- numeric()
    return(empty %>% dplyr::select("rep", dplyr::all_of(key_cols), "estimate"))
  }

  draws_raw <- purrr::map_dfr(successful_reps, function(rep_i){
    results[[rep_i]] %>%
      dplyr::mutate(rep = rep_i, .before = 1)
  })

  tidyr::expand_grid(
    rep = successful_reps,
    canonical %>% dplyr::select(dplyr::all_of(key_cols))
  ) %>%
    dplyr::left_join(
      draws_raw %>% dplyr::select("rep", dplyr::all_of(key_cols), "estimate"),
      by = c("rep", key_cols)
    ) %>%
    dplyr::arrange(
      .data$rep,
      .data$qoi,
      dplyr::across(dplyr::any_of("profile")),
      .data$et,
      .data$term
    )
}


cbwsdid_bootstrap_summarise <- function(original_qoi,
                                        draws,
                                        ci_type = "percentile",
                                        level = 0.95){
  alpha <- 1 - level
  crit <- stats::qnorm(1 - alpha / 2)
  key_cols <- intersect(c("qoi", "term", "profile", "et"), names(original_qoi))

  if(nrow(draws) == 0L){
    stats <- original_qoi %>%
      dplyr::transmute(
        dplyr::across(dplyr::all_of(key_cols)),
        boot_se = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_,
        boot_mean = NA_real_,
        n_eff = 0L
      )
  } else {
    stats <- draws %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(key_cols))) %>%
      dplyr::summarise(
        n_eff = sum(!is.na(.data$estimate)),
        boot_mean = cbwsdid_bootstrap_mean(.data$estimate),
        boot_se = cbwsdid_bootstrap_sd(.data$estimate),
        conf.low = cbwsdid_bootstrap_ci_low(
          .data$estimate,
          point = NA_real_,
          boot_se = NA_real_,
          ci_type = ci_type,
          alpha = alpha,
          crit = crit
        ),
        conf.high = cbwsdid_bootstrap_ci_high(
          .data$estimate,
          point = NA_real_,
          boot_se = NA_real_,
          ci_type = ci_type,
          alpha = alpha,
          crit = crit
        ),
        .groups = "drop"
      )
  }

  out <- original_qoi %>%
    dplyr::left_join(
      stats %>%
        dplyr::select(dplyr::all_of(key_cols), "boot_se", "conf.low",
                      "conf.high", "boot_mean", "n_eff"),
      by = key_cols
    )

  if(identical(ci_type, "normal")){
    out <- out %>%
      dplyr::mutate(
        conf.low = .data$estimate - crit * .data$boot_se,
        conf.high = .data$estimate + crit * .data$boot_se
      )
  }

  known_cols <- c(
    "qoi", "term", "profile", "et", "estimate", "boot_se",
    "conf.low", "conf.high", "boot_mean", "boot_bias", "n_eff"
  )
  extra_cols <- setdiff(names(out), known_cols)

  out %>%
    dplyr::mutate(boot_bias = .data$boot_mean - .data$estimate) %>%
    dplyr::select(
      "qoi",
      dplyr::any_of("profile"),
      dplyr::all_of(extra_cols),
      "term",
      "et",
      "estimate",
      "boot_se",
      "conf.low",
      "conf.high",
      "boot_mean",
      "boot_bias",
      "n_eff"
    )
}


cbwsdid_bootstrap_mean <- function(x){
  x <- x[!is.na(x)]
  if(length(x) == 0L){
    return(NA_real_)
  }
  mean(x)
}


cbwsdid_bootstrap_sd <- function(x){
  x <- x[!is.na(x)]
  if(length(x) < 2L){
    return(NA_real_)
  }
  stats::sd(x)
}


cbwsdid_bootstrap_ci_low <- function(x,
                                     point,
                                     boot_se,
                                     ci_type,
                                     alpha,
                                     crit){
  x <- x[!is.na(x)]
  if(length(x) == 0L){
    return(NA_real_)
  }

  if(identical(ci_type, "percentile")){
    return(unname(stats::quantile(x, probs = alpha / 2, na.rm = TRUE)))
  }

  point - crit * boot_se
}


cbwsdid_bootstrap_ci_high <- function(x,
                                      point,
                                      boot_se,
                                      ci_type,
                                      alpha,
                                      crit){
  x <- x[!is.na(x)]
  if(length(x) == 0L){
    return(NA_real_)
  }

  if(identical(ci_type, "percentile")){
    return(unname(stats::quantile(x, probs = 1 - alpha / 2, na.rm = TRUE)))
  }

  point + crit * boot_se
}


#' @export
print.cbwsdid_bootstrap <- function(x, ...){
  cat("cbwsdid pipeline cluster bootstrap\n")
  cat("Cluster: ", x$cluster, "\n", sep = "")
  cat("Replications: ", x$B_effective, "/", x$B, sep = "")
  if(x$n_failed > 0L){
    cat(" (", x$n_failed, " failed)", sep = "")
  }
  cat("\n")
  if("n_eff" %in% names(x$summary) && nrow(x$summary) > 0L){
    min_eff <- suppressWarnings(min(x$summary$n_eff, na.rm = TRUE))
    if(is.finite(min_eff)){
      cat("Minimum term effective B: ", min_eff, "\n", sep = "")
    }
  }
  print(x$summary, ...)
  invisible(x)
}
