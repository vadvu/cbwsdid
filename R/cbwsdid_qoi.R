#' Quantities of Interest from `cbwsdid`
#'
#' Extracts aggregated quantities of interest from a fitted [cbwsdid()] model.
#'
#' @param model `fixest` model object from [cbwsdid()].
#' @param type Character. Quantity of interest to compute. `"dynamic"` returns
#'   an event-time table convenient for plotting. `"simple"` returns the
#'   equally weighted post-treatment average. `"pooled"` returns the
#'   placebo-corrected overall effect: the equally weighted post-treatment
#'   average minus the equally weighted pre-treatment placebo average over
#'   event times at or before `-2`. `"cohort"` returns cohort- or episode-level
#'   dynamic and simple effects.
#' @param model_type Character. Which second-stage model to use. `"overall"` uses the baseline unmoderated models. `"moderated"` uses the moderated second-stage models stored in the metadata.
#' @param post Optional. Integer array of non-negative event times to include in the simple ATT aggregation. By default, all available post-treatment event times are used.
#' @param newdata Optional data frame describing moderator values for moderated quantities of interest. If `NULL`, moderated effects are evaluated at the treated unit-by-stack mean moderator profile from the original fit.
#' @param level Confidence level for interval estimates.
#' @param se.method Standard-error method. `"conditional"` (the default) uses
#'   the cluster-robust covariance matrix from the fitted second-stage
#'   regression. `"analytic"` uses the unit-clustered influence function.
#'   `"bootstrap"` reruns the full estimation pipeline through
#'   [cbwsdid_bootstrap()].
#' @param data Optional raw panel used only for `se.method = "bootstrap"`. By
#'   default, it is taken from the fitted model when [cbwsdid()] was called with
#'   `keep_data = TRUE`.
#' @param B Number of pipeline-bootstrap replications.
#' @param seed Optional bootstrap random seed.
#' @param parallel Logical. If `TRUE`, run bootstrap replications in parallel.
#' @param workers Optional number of bootstrap workers.
#'
#' @details
#' Analytic inference clusters influence-function contributions by unit and
#' stores the full event-time covariance matrix in the `"vcov"` attribute of
#' the returned tibble. For smooth weighting without an outcome model, it
#' treats estimated balancing weights as fixed; the pipeline bootstrap also
#' captures weight re-estimation.
#'
#' For `model_type = "overall"`, the pooled estimate and all three uncertainty
#' methods use the same linear functional of the dynamic coefficients:
#' the equally weighted post-treatment average minus the equally weighted
#' pre-treatment placebo average over event times at or before `-2`. Thus the
#' overall pooled quantity does not require `pooled = TRUE`. Moderated pooled
#' effects continue to use the separately fitted pooled moderation model.
#'
#' Moderated QOIs are marginal effects
#' `theta_e(m) = beta_e + psi_e' m`, where `m` is an expanded moderator profile.
#' The moderated second-stage model includes both moderator-by-event-time main
#' effects and treatment-by-moderator-by-event-time interactions. The main
#' effects capture moderator-dependent counterfactual trends and are not part of
#' the treatment-effect profile combo. The `beta_e` coefficients are effects at
#' moderator value zero; use `newdata` profiles, the default treated-profile
#' mean, or [cbwsdid_qoi_moderation_slopes()] rather than reading them as
#' aggregate effects.
#'
#' For moderated bootstrap QOIs, moderator profiles are resolved once from the
#' original model and then frozen across bootstrap draws. With explicit
#' `newdata`, this is the requested fixed evaluation point. With the default
#' profile, the bootstrap targets the effect at the original treated mean
#' profile rather than re-estimating a new mean profile in each draw.
#'
#' For doubly robust nearest-neighbour matching, the analytic variance is
#' conditional on the realized matches because the Abadie--Imbens fixed-match
#' count term is not yet included. For pure matching without an outcome model,
#' `se.method = "analytic"` warns and falls back to the existing conditional
#' variance. Analytic inference is available for overall and cohort quantities.
#' For moderated QOIs, use conditional or bootstrap inference.
#' @return A tibble. For `type = "dynamic"`, one row per event time. For `type = "simple"` and `type = "pooled"`, one row per requested moderator profile. For `type = "cohort"`, a list with `dynamic` and `simple` tibbles.
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
#' cbwsdid_qoi(fit, type = "simple")
#' cbwsdid_qoi(fit, type = "dynamic")
#' cbwsdid_qoi(fit, type = "pooled")
#' cohort_qoi <- cbwsdid_qoi(fit, type = "cohort")
#' cohort_qoi$dynamic
#' cohort_qoi$simple
cbwsdid_qoi <- function(model,
                        type = c("dynamic", "simple", "pooled", "cohort"),
                        model_type = c("overall", "moderated"),
                        post = NULL,
                        newdata = NULL,
                        level = 0.95,
                        se.method = c("conditional", "analytic", "bootstrap"),
                        data = NULL,
                        B = 999L,
                        seed = NULL,
                        parallel = FALSE,
                        workers = NULL){
  type <- match.arg(type)
  model_type <- match.arg(model_type)
  se.method <- match.arg(se.method)

  if(!inherits(model, "fixest")){
    stop("`model` must be a fixest object returned by `cbwsdid()`.")
  }

  meta <- attr(model, "cbwsdid")
  alpha <- 1 - level
  crit <- stats::qnorm(1 - alpha / 2)

  if(identical(se.method, "bootstrap")){
    return(cbwsdid_qoi_bootstrap(
      model = model,
      type = type,
      model_type = model_type,
      post = post,
      newdata = newdata,
      level = level,
      data = data,
      B = B,
      seed = seed,
      parallel = parallel,
      workers = workers
    ))
  }

  if(identical(se.method, "analytic")){
    if(!identical(model_type, "overall")){
      stop(
        "`se.method = \"analytic\"` is not available for moderated QOIs; ",
        "use `se.method = \"conditional\"` or `se.method = \"bootstrap\"`."
      )
    }

    if(identical(se.method, "analytic") &&
       identical(meta$refinement.method, "matchit") &&
       !isTRUE(meta$dr)){
      warning(
        "Analytic unconditional variance is not valid for pure fixed-k ",
        "matching without bias correction; using conditional-on-matching ",
        "standard errors. Prefer doubly robust matching with bootstrap or ",
        "subsampling inference.",
        call. = FALSE
      )
      se.method <- "conditional"
    }
  }

  if(type == "cohort"){
    if(!identical(model_type, "overall")){
      stop("`type = \"cohort\"` is only available for `model_type = \"overall\"`.")
    }

    return(cbwsdid_qoi_cohort(
      model,
      post = post,
      level = level,
      se.method = se.method
    ))
  }

  if(identical(model_type, "overall")){
    if(type == "dynamic"){
      if(identical(se.method, "analytic")){
        return(cbwsdid_qoi_dynamic_analytic(model, level = level))
      }
      return(cbwsdid_qoi_dynamic_overall(model, level = level))
    }

    if(type == "simple"){
      if(identical(se.method, "analytic")){
        return(cbwsdid_qoi_simple_analytic(
          model,
          post = post,
          level = level
        ))
      }
      return(cbwsdid_qoi_simple_overall(model, post = post, level = level))
    }

    if(identical(se.method, "analytic")){
      return(cbwsdid_qoi_pooled_analytic(
        model = model,
        level = level
      ))
    }

    return(cbwsdid_qoi_pooled_overall(model, level = level))
  }

  if(is.null(meta$moderation.formula) || nrow(meta$moderation.spec) == 0){
    stop(
      "This `cbwsdid` object does not contain moderated second-stage models. ",
      "Re-estimate with `moderation.formula`."
    )
  }

  profiles <- cbwsdid_prepare_moderator_profiles(
    meta = meta,
    newdata = newdata
  )

  if(type == "dynamic"){
    moderation_model <- meta$moderation_model
    if(is.null(moderation_model)){
      stop("No moderated dynamic model is stored in this `cbwsdid` object.")
    }

    dynamic_terms <- stats::coef(moderation_model)
    dynamic_terms <- names(dynamic_terms)[grepl("^et::(-?[0-9]+):treated_sa$", names(dynamic_terms))]
    if(length(dynamic_terms) == 0){
      stop("No moderated dynamic treatment coefficients were found.")
    }

    dynamic_ets <- sort(as.numeric(stringr::str_match(dynamic_terms, "^et::(-?[0-9]+):treated_sa$")[, 2]))

    out <- purrr::map_dfr(seq_len(nrow(profiles$expanded)), function(i){
      profile_i <- profiles$expanded[i, , drop = FALSE]
      dynamic_i <- purrr::map_dfr(dynamic_ets, function(et_i){
        combo_i <- cbwsdid_qoi_moderated_dynamic_combo(
          model = moderation_model,
          et = et_i,
          moderator_values = profile_i
        )

        tibble::tibble(
          et = et_i,
          estimate = combo_i$estimate,
          std.error = combo_i$std.error,
          conf.low = combo_i$conf.low,
          conf.high = combo_i$conf.high
        )
      })

      dplyr::bind_cols(profiles$profile[i, , drop = FALSE], dynamic_i)
    }) %>%
      dplyr::mutate(
        model_type = "moderated",
        type = "dynamic",
        .before = 1
      )

    return(out)
  }

  if(type == "simple"){
    moderation_model <- meta$moderation_model
    if(is.null(moderation_model)){
      stop("No moderated dynamic model is stored in this `cbwsdid` object.")
    }

    overall_dynamic <- cbwsdid_qoi_dynamic_overall(model, level = level)

    if(is.null(post)){
      post <- overall_dynamic$et[overall_dynamic$et >= 0]
    }

    post <- sort(unique(post))
    if(length(post) == 0){
      stop("No post-treatment event times were selected.")
    }

    missing_post <- setdiff(post, overall_dynamic$et)
    if(length(missing_post) > 0){
      stop(
        "Some requested event times are not present in the model: ",
        paste(missing_post, collapse = ", ")
      )
    }

    out <- purrr::map_dfr(seq_len(nrow(profiles$expanded)), function(i){
      combo_i <- cbwsdid_qoi_moderated_simple_combo(
        model = moderation_model,
        post = post,
        moderator_values = profiles$expanded[i, , drop = FALSE],
        level = level
      )

      dplyr::bind_cols(
        profiles$profile[i, , drop = FALSE],
        tibble::tibble(
          model_type = "moderated",
          type = "simple",
          periods = paste0(min(post), " - ", max(post)),
          n_periods = length(post),
          estimate = combo_i$estimate,
          std.error = combo_i$std.error,
          conf.low = combo_i$conf.low,
          conf.high = combo_i$conf.high
        )
      )
    })

    return(out)
  }

  pooled_moderation_model <- meta$pooled_moderation_model
  if(is.null(pooled_moderation_model)){
    stop(
      "No moderated pooled model is stored in this `cbwsdid` object. ",
      "Re-estimate with `pooled = TRUE` and `moderation.formula`."
    )
  }

  out <- purrr::map_dfr(seq_len(nrow(profiles$expanded)), function(i){
    combo_i <- cbwsdid_qoi_moderated_pooled_combo(
      model = pooled_moderation_model,
      moderator_values = profiles$expanded[i, , drop = FALSE],
      level = level
    )

    dplyr::bind_cols(
      profiles$profile[i, , drop = FALSE],
      tibble::tibble(
        model_type = "moderated",
        type = "pooled",
        estimate = combo_i$estimate,
        std.error = combo_i$std.error,
        conf.low = combo_i$conf.low,
        conf.high = combo_i$conf.high
      )
    )
  })

  out
}


#' Moderation Slopes from a Moderated `cbwsdid` Model
#'
#' Extracts the event-time path of `psi_e` coefficients from the
#' moderated dynamic second-stage model.
#'
#' @param model `fixest` model object from [cbwsdid()] fitted with
#'   `moderation.formula`.
#' @param level Confidence level for interval estimates.
#'
#' @details
#' Each row is the coefficient on
#' `i(et, treated_sa * moderator, ref = -1)` for one expanded moderator basis
#' column. These slopes are identified as differences in moderator gradients of
#' outcome changes between treated units and weighted controls. Pre-treatment
#' rows are placebo checks for moderator-dependent untreated trends.
#'
#' @return A tibble with expanded moderator name, event time, estimate,
#'   cluster-robust standard error, and confidence interval.
#' @export
cbwsdid_qoi_moderation_slopes <- function(model, level = 0.95){
  if(!inherits(model, "fixest")){
    stop("`model` must be a fixest object returned by `cbwsdid()`.")
  }

  meta <- attr(model, "cbwsdid")
  if(is.null(meta) || is.null(meta$moderation.formula) ||
     nrow(meta$moderation.spec) == 0){
    stop(
      "This `cbwsdid` object does not contain moderated second-stage models. ",
      "Re-estimate with `moderation.formula`."
    )
  }

  moderation_model <- meta$moderation_model
  if(is.null(moderation_model)){
    stop("No moderated dynamic model is stored in this `cbwsdid` object.")
  }

  coef_all <- stats::coef(moderation_model)
  vcov_all <- stats::vcov(moderation_model)

  # Build slope term names exactly from the stored expanded moderator columns:
  # an open-ended regex on `treated_x_` would also match the main-effect terms
  # of a moderator whose own name starts with `treated_x_`.
  expanded_cols <- meta$moderation.expanded_cols
  et_pattern <- "^et::(-?[0-9]+):treated_sa$"
  et_terms <- names(coef_all)[grepl(et_pattern, names(coef_all))]
  ets <- sort(as.numeric(stringr::str_match(et_terms, et_pattern)[, 2]))

  slope_grid <- tidyr::expand_grid(moderator = expanded_cols, et = ets) %>%
    dplyr::mutate(term = paste0("et::", .data$et, ":treated_x_", .data$moderator)) %>%
    dplyr::filter(.data$term %in% names(coef_all))
  if(nrow(slope_grid) == 0L){
    stop("No moderated dynamic slope coefficients were found.")
  }

  slope_terms <- slope_grid$term
  alpha <- 1 - level
  crit <- stats::qnorm(1 - alpha / 2)

  tibble::tibble(
    model_type = "moderated",
    type = "moderation_slope",
    moderator = slope_grid$moderator,
    et = slope_grid$et,
    term = slope_terms,
    estimate = unname(coef_all[slope_terms]),
    std.error = sqrt(pmax(diag(vcov_all[slope_terms, slope_terms, drop = FALSE]), 0))
  ) %>%
    dplyr::arrange(.data$moderator, .data$et) %>%
    dplyr::mutate(
      conf.low = .data$estimate - crit * .data$std.error,
      conf.high = .data$estimate + crit * .data$std.error
    ) %>%
    dplyr::select(
      "model_type",
      "type",
      "moderator",
      "et",
      "estimate",
      "std.error",
      "conf.low",
      "conf.high"
    )
}


cbwsdid_qoi_bootstrap <- function(model,
                                  type,
                                  model_type,
                                  post,
                                  newdata,
                                  level,
                                  data,
                                  B,
                                  seed,
                                  parallel,
                                  workers){
  if(identical(type, "cohort")){
    stop(
      "`se.method = \"bootstrap\"` is not yet available for ",
      "`type = \"cohort\"`; use `cbwsdid_bootstrap()` for supported overall ",
      "quantities."
    )
  }
  if(identical(model_type, "overall") && !is.null(newdata)){
    stop("`newdata` is not used for overall bootstrap quantities.")
  }

  meta <- attr(model, "cbwsdid")
  if(is.null(data)){
    data <- meta$data
  }
  if(is.null(data)){
    stop(
      "No raw panel is stored in `model`. Refit with `keep_data = TRUE` ",
      "or pass the original panel through `data=`."
    )
  }

  bootstrap <- cbwsdid_bootstrap(
    model = model,
    data = data,
    B = B,
    qoi = type,
    post = post,
    model_type = model_type,
    newdata = newdata,
    level = level,
    parallel = parallel,
    workers = workers,
    seed = seed,
    progress = FALSE
  )
  out <- cbwsdid_qoi(
    model = model,
    type = type,
    model_type = model_type,
    post = post,
    newdata = newdata,
    level = level,
    se.method = "conditional"
  )
  boot.summary <- bootstrap$summary

  if(identical(model_type, "moderated")){
    out.term <- cbwsdid_qoi_moderated_bootstrap_terms(out, type = type)
    boot.row <- match(out.term, boot.summary$term)
  } else if(identical(type, "dynamic")){
    boot.row <- match(out$et, boot.summary$et)
  } else {
    boot.row <- rep(1L, nrow(out))
  }

  out$std.error <- boot.summary$boot_se[boot.row]
  out$conf.low <- boot.summary$conf.low[boot.row]
  out$conf.high <- boot.summary$conf.high[boot.row]
  attr(out, "bootstrap") <- bootstrap
  attr(out, "se.method") <- "bootstrap"
  out
}


cbwsdid_qoi_dynamic_overall <- function(model, level = 0.95){
  coef_all <- stats::coef(model)
  vcov_all <- stats::vcov(model)
  alpha <- 1 - level
  crit <- stats::qnorm(1 - alpha / 2)

  es_terms <- names(coef_all)[grepl("^et::(-?[0-9]+):treated_sa$", names(coef_all))]
  if(length(es_terms) == 0){
    stop("No event-study coefficients matching `treated_sa` were found.")
  }

  tibble::tibble(
    model_type = "overall",
    type = "dynamic",
    term = es_terms,
    et = as.numeric(stringr::str_match(es_terms, "^et::(-?[0-9]+):treated_sa$")[, 2]),
    estimate = unname(coef_all[es_terms]),
    std.error = sqrt(diag(vcov_all[es_terms, es_terms, drop = FALSE]))
  ) %>%
    dplyr::arrange(.data$et) %>%
    dplyr::mutate(
      conf.low = .data$estimate - crit * .data$std.error,
      conf.high = .data$estimate + crit * .data$std.error
    ) %>%
    dplyr::select("model_type",
                  "type",
                  "et",
                  "estimate",
                  "std.error",
                  "conf.low",
                  "conf.high")
}


cbwsdid_qoi_dynamic_analytic <- function(model, level = 0.95){
  out <- cbwsdid_qoi_dynamic_overall(model, level = level)
  analytic <- cbwsdid_vcov_analytic(model, event_times = out$et)
  crit <- stats::qnorm(1 - (1 - level) / 2)

  out$std.error <- unname(analytic$std.error[as.character(out$et)])
  out$conf.low <- out$estimate - crit * out$std.error
  out$conf.high <- out$estimate + crit * out$std.error
  attr(out, "vcov") <- analytic$vcov
  attr(out, "influence") <- analytic$influence
  attr(out, "se.method") <- "analytic"
  out
}


cbwsdid_qoi_simple_overall <- function(model, post = NULL, level = 0.95){
  dynamic_df <- cbwsdid_qoi_dynamic_overall(model, level = level)

  if(is.null(post)){
    post <- dynamic_df$et[dynamic_df$et >= 0]
  }

  post <- sort(unique(post))
  if(length(post) == 0){
    stop("No post-treatment event times were selected.")
  }

  simple_df <- dynamic_df %>%
    dplyr::filter(.data$et %in% post)

  if(nrow(simple_df) != length(post)){
    missing_post <- setdiff(post, simple_df$et)
    stop(
      "Some requested event times are not present in the model: ",
      paste(missing_post, collapse = ", ")
    )
  }

  weight_vec <- rep(1 / length(post), length(post))
  names(weight_vec) <- paste0("et::", simple_df$et, ":treated_sa")

  combo <- cbwsdid_qoi_linear_combo(
    model = model,
    weights = weight_vec,
    level = level
  )

  tibble::tibble(
    model_type = "overall",
    type = "simple",
    periods = paste0(min(post), " - ", max(post)),
    n_periods = length(post),
    estimate = combo$estimate,
    std.error = combo$std.error,
    conf.low = combo$conf.low,
    conf.high = combo$conf.high
  )
}


cbwsdid_qoi_simple_analytic <- function(model, post = NULL, level = 0.95){
  dynamic <- cbwsdid_qoi_dynamic_overall(model, level = level)
  if(is.null(post)){
    post <- dynamic$et[dynamic$et >= 0]
  }

  post <- sort(unique(post))
  if(length(post) == 0L){
    stop("No post-treatment event times were selected.")
  }

  missing.post <- setdiff(post, dynamic$et)
  if(length(missing.post) > 0L){
    stop(
      "Some requested event times are not present in the model: ",
      paste(missing.post, collapse = ", ")
    )
  }

  out <- cbwsdid_qoi_simple_overall(model, post = post, level = level)
  analytic <- cbwsdid_vcov_analytic(model, event_times = dynamic$et)
  post.columns <- match(as.character(post), colnames(analytic$influence))
  simple.influence <- rowMeans(
    analytic$influence[, post.columns, drop = FALSE]
  )
  simple.se <- sqrt(sum(simple.influence^2))
  crit <- stats::qnorm(1 - (1 - level) / 2)

  out$std.error <- simple.se
  out$conf.low <- out$estimate - crit * simple.se
  out$conf.high <- out$estimate + crit * simple.se
  attr(out, "vcov") <- analytic$vcov
  attr(out, "influence") <- simple.influence
  attr(out, "se.method") <- "analytic"
  out
}


cbwsdid_qoi_pooled_weights <- function(dynamic){
  pre <- dynamic$et <= -2
  post <- dynamic$et >= 0

  if(!any(pre)){
    stop(
      "`type = \"pooled\"` requires at least one pre-treatment period ",
      "(`kappa_pre <= -2`); use `type = \"simple\"` for the ",
      "post-treatment average."
    )
  }
  if(!any(post)){
    stop("`type = \"pooled\"` requires at least one post-treatment period.")
  }

  weights <- ifelse(
    post,
    1 / sum(post),
    -1 / sum(pre)
  )
  weights[!pre & !post] <- 0
  stats::setNames(
    weights,
    paste0("et::", dynamic$et, ":treated_sa")
  )
}


cbwsdid_qoi_pooled_overall <- function(model, level = 0.95){
  dynamic <- cbwsdid_qoi_dynamic_overall(model, level = level)
  pooled.weights <- cbwsdid_qoi_pooled_weights(dynamic)
  pooled.estimate <- cbwsdid_qoi_linear_combo(
    model,
    weights = pooled.weights,
    level = level
  )

  tibble::tibble(
    model_type = "overall",
    type = "pooled",
    estimate = pooled.estimate$estimate,
    std.error = pooled.estimate$std.error,
    conf.low = pooled.estimate$conf.low,
    conf.high = pooled.estimate$conf.high
  )
}


cbwsdid_qoi_pooled_analytic <- function(model, level = 0.95){
  out <- cbwsdid_qoi_pooled_overall(model, level = level)
  dynamic <- cbwsdid_qoi_dynamic_overall(model, level = level)
  pooled.weights <- unname(cbwsdid_qoi_pooled_weights(dynamic))

  analytic <- cbwsdid_vcov_analytic(model, event_times = dynamic$et)
  pooled.influence <- as.numeric(analytic$influence %*% pooled.weights)
  pooled.se <- sqrt(sum(pooled.influence^2))
  crit <- stats::qnorm(1 - (1 - level) / 2)

  out$std.error <- pooled.se
  out$conf.low <- out$estimate - crit * pooled.se
  out$conf.high <- out$estimate + crit * pooled.se
  attr(out, "vcov") <- analytic$vcov
  attr(out, "influence") <- pooled.influence
  attr(out, "se.method") <- "analytic"
  out
}


cbwsdid_qoi_cohort <- function(model,
                               post = NULL,
                               level = 0.95,
                               se.method = "conditional"){
  meta <- attr(model, "cbwsdid")
  if(is.null(meta)){
    stop("`model` does not contain `cbwsdid` metadata.")
  }

  stacked.data <- meta$stacked.data
  if(is.null(stacked.data)){
    stop(
      "`type = \"cohort\"` requires the stacked estimation sample. ",
      "Re-estimate with `keep_data = TRUE`."
    )
  }

  cohort_fit <- cbwsdid_fit_cohort_model(stacked.data)
  dynamic <- cbwsdid_qoi_cohort_dynamic(cohort_fit, meta = meta, level = level)
  simple <- cbwsdid_qoi_cohort_simple(
    cohort_fit = cohort_fit,
    dynamic = dynamic,
    meta = meta,
    post = post,
    level = level
  )

  if(identical(se.method, "analytic")){
    analytic <- cbwsdid_qoi_cohort_analytic(
      model = model,
      dynamic = dynamic,
      simple = simple,
      post = post,
      level = level
    )
    dynamic <- analytic$dynamic
    simple <- analytic$simple
  }

  out <- list(
    dynamic = dynamic,
    simple = simple
  )
  attr(out, "se.method") <- se.method
  out
}


cbwsdid_qoi_cohort_analytic <- function(model,
                                        dynamic,
                                        simple,
                                        post = NULL,
                                        level = 0.95){
  influence <- cbwsdid_influence_cohort(
    model,
    event_times = sort(unique(dynamic$et))
  ) %>%
    dplyr::mutate(subexp.key = as.character(.data$subexp.id))
  crit <- stats::qnorm(1 - (1 - level) / 2)

  cell.variance <- influence %>%
    dplyr::group_by(.data$subexp.key, .data$et) %>%
    dplyr::summarise(
      variance = sum(.data$contribution^2),
      .groups = "drop"
    )
  dynamic.key <- paste(as.character(dynamic$subexp.id), dynamic$et, sep = "\r")
  variance.key <- paste(
    cell.variance$subexp.key,
    cell.variance$et,
    sep = "\r"
  )
  dynamic$std.error <- sqrt(pmax(
    cell.variance$variance[match(dynamic.key, variance.key)],
    0
  ))
  dynamic$conf.low <- dynamic$estimate - crit * dynamic$std.error
  dynamic$conf.high <- dynamic$estimate + crit * dynamic$std.error

  if(is.null(post)){
    post <- sort(unique(dynamic$et[dynamic$et >= 0]))
  } else {
    post <- sort(unique(post))
  }
  # Psi^s_i = |H+|^{-1} * sum_{e in post} Psi_{i,e}, with Psi_{i,e} = 0 when
  # unit i does not contribute at e (zero-fill, matching the overall simple
  # path) — dividing the partial sum by length(post) implements the zero-fill
  # without materializing absent rows.
  simple.influence <- influence %>%
    dplyr::filter(.data$et %in% post) %>%
    dplyr::group_by(.data$id, .data$subexp.key) %>%
    dplyr::summarise(
      contribution = sum(.data$contribution) / length(post),
      .groups = "drop"
    ) %>%
    dplyr::group_by(.data$subexp.key) %>%
    dplyr::summarise(
      variance = sum(.data$contribution^2),
      .groups = "drop"
    )
  simple$std.error <- sqrt(pmax(
    simple.influence$variance[
      match(as.character(simple$subexp.id), simple.influence$subexp.key)
    ],
    0
  ))
  simple$conf.low <- simple$estimate - crit * simple$std.error
  simple$conf.high <- simple$estimate + crit * simple$std.error

  attr(dynamic, "influence") <- influence
  attr(dynamic, "se.method") <- "analytic"
  attr(simple, "se.method") <- "analytic"

  list(
    dynamic = dynamic,
    simple = simple
  )
}


cbwsdid_fit_cohort_model <- function(stacked.data){
  required.cols <- c("id", "time", "y", "subexp.id", "et", "treated_sa", "Qsa")
  missing.cols <- setdiff(required.cols, names(stacked.data))
  if(length(missing.cols) > 0){
    stop(
      "`stacked.data` is missing required column(s): ",
      paste(missing.cols, collapse = ", ")
    )
  }

  model.data <- stacked.data %>%
    dplyr::mutate(subexp.id = as.character(.data$subexp.id))

  # One coefficient per treated (subexp.id, et) cell, with et = -1 as the
  # within-stack reference. The cells are encoded in a single factor column so
  # that fixest builds the interaction sparsely through i(), instead of
  # materializing one dense dummy column per cell. The two are algebraically
  # identical, but the dense version blows up memory for switch designs with
  # many episode types.
  cell.sep <- "@@"
  model.data <- model.data %>%
    dplyr::mutate(
      ce_cell = dplyr::if_else(
        .data$treated_sa == 1 & .data$et != -1,
        paste0(.data$subexp.id, cell.sep, .data$et),
        "0"
      )
    )

  term.map <- model.data %>%
    dplyr::filter(.data$treated_sa == 1, .data$et != -1) %>%
    dplyr::group_by(.data$subexp.id, .data$et) %>%
    dplyr::summarise(n_treated = dplyr::n_distinct(.data$id), .groups = "drop") %>%
    dplyr::arrange(.data$subexp.id, .data$et) %>%
    dplyr::mutate(term = paste0("ce_cell::", .data$subexp.id, cell.sep, .data$et))

  if(nrow(term.map) == 0){
    stop("No cohort-event coefficients are available outside the reference period `et = -1`.")
  }

  cohort.model <- feols(
    y ~ i(ce_cell, ref = "0") | id^subexp.id + time^subexp.id,
    data = model.data,
    cluster = ~id,
    weights = ~Qsa
  )

  list(
    model = cohort.model,
    term.map = term.map
  )
}


cbwsdid_qoi_cohort_dynamic <- function(cohort_fit, meta, level = 0.95){
  coef_all <- stats::coef(cohort_fit$model)
  vcov_all <- stats::vcov(cohort_fit$model)
  term.map <- cohort_fit$term.map
  alpha <- 1 - level
  crit <- stats::qnorm(1 - alpha / 2)

  out <- term.map %>%
    dplyr::mutate(
      estimate = purrr::map_dbl(.data$term, function(term_i){
        if(!term_i %in% names(coef_all)){
          return(NA_real_)
        }
        unname(coef_all[term_i])
      }),
      std.error = purrr::map_dbl(.data$term, function(term_i){
        if(!term_i %in% rownames(vcov_all) || !term_i %in% colnames(vcov_all)){
          return(NA_real_)
        }
        sqrt(pmax(unname(vcov_all[term_i, term_i]), 0))
      })
    ) %>%
    dplyr::left_join(
      cbwsdid_cohort_info(meta) %>% dplyr::select(-dplyr::any_of("n_treated")),
      by = "subexp.id"
    ) %>%
    dplyr::mutate(
      conf.low = .data$estimate - crit * .data$std.error,
      conf.high = .data$estimate + crit * .data$std.error,
      contribution = .data$cohort_weight * .data$estimate
    )

  cbwsdid_order_cohort_output(out, meta = meta) %>%
    dplyr::select(
      dplyr::any_of(c("subexp.id", "cohort", "switch_time", "history_signature")),
      "et",
      "estimate",
      "std.error",
      "conf.low",
      "conf.high",
      "cohort_weight",
      "contribution",
      "n_treated",
      "n_control",
      "control_mass"
    )
}


cbwsdid_qoi_cohort_simple <- function(cohort_fit,
                                      dynamic,
                                      meta,
                                      post = NULL,
                                      level = 0.95){
  if(is.null(post)){
    post <- dynamic$et[dynamic$et >= 0]
  }

  post <- sort(unique(post))
  if(length(post) == 0){
    stop("No post-treatment event times were selected.")
  }

  missing_post <- setdiff(post, unique(dynamic$et))
  if(length(missing_post) > 0){
    stop(
      "Some requested event times are not present in the cohort model: ",
      paste(missing_post, collapse = ", ")
    )
  }

  coef_all <- stats::coef(cohort_fit$model)
  vcov_all <- stats::vcov(cohort_fit$model)

  simple <- cohort_fit$term.map %>%
    dplyr::filter(.data$et %in% post) %>%
    dplyr::group_by(.data$subexp.id) %>%
    dplyr::group_modify(function(.x, .y){
      if(nrow(.x) != length(post)){
        return(
          tibble::tibble(
            periods = paste0(min(post), " - ", max(post)),
            n_periods = length(post),
            estimate = NA_real_,
            std.error = NA_real_,
            conf.low = NA_real_,
            conf.high = NA_real_
          )
        )
      }

      combo <- cbwsdid_qoi_terms_combo(
        coef_all = coef_all,
        vcov_all = vcov_all,
        terms = .x$term,
        weights = rep(1 / length(post), length(post)),
        level = level
      )

      tibble::tibble(
        periods = paste0(min(post), " - ", max(post)),
        n_periods = length(post),
        estimate = combo$estimate,
        std.error = combo$std.error,
        conf.low = combo$conf.low,
        conf.high = combo$conf.high
      )
    }) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(cbwsdid_cohort_info(meta), by = "subexp.id") %>%
    dplyr::mutate(contribution = .data$cohort_weight * .data$estimate)

  cbwsdid_order_cohort_output(simple, meta = meta) %>%
    dplyr::select(
      dplyr::any_of(c("subexp.id", "cohort", "switch_time", "history_signature")),
      "periods",
      "n_periods",
      "estimate",
      "std.error",
      "conf.low",
      "conf.high",
      "cohort_weight",
      "contribution",
      "n_treated",
      "n_control",
      "control_mass"
    )
}


cbwsdid_qoi_terms_combo <- function(coef_all,
                                    vcov_all,
                                    terms,
                                    weights,
                                    level = 0.95){
  alpha <- 1 - level
  crit <- stats::qnorm(1 - alpha / 2)
  missing_terms <- setdiff(terms, names(coef_all))

  if(length(missing_terms) > 0 ||
     any(!terms %in% rownames(vcov_all)) ||
     any(!terms %in% colnames(vcov_all))){
    return(
      list(
        estimate = NA_real_,
        std.error = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_
      )
    )
  }

  estimate <- sum(weights * unname(coef_all[terms]))
  variance <- as.numeric(t(weights) %*% vcov_all[terms, terms, drop = FALSE] %*% weights)
  std.error <- sqrt(pmax(variance, 0))

  list(
    estimate = estimate,
    std.error = std.error,
    conf.low = estimate - crit * std.error,
    conf.high = estimate + crit * std.error
  )
}


cbwsdid_cohort_info <- function(meta){
  cohort.weights <- meta$cohort.weights
  if(is.null(cohort.weights) || nrow(cohort.weights) == 0){
    stop("No cohort weights are stored in this `cbwsdid` object.")
  }

  out <- cohort.weights %>%
    dplyr::mutate(
      subexp.id = as.character(.data$subexp.id),
      cohort_weight = .data$Nda / sum(.data$Nda, na.rm = TRUE),
      n_treated = .data$Nda,
      n_control = .data$Nca,
      control_mass = .data$Bca
    )

  if(identical(meta$design, "absorbing")){
    out <- out %>%
      dplyr::mutate(
        cohort = suppressWarnings(as.numeric(.data$subexp.id)),
        .before = "subexp.id"
      )
  }

  if(identical(meta$design, "switch01") || identical(meta$design, "switch10")){
    out <- out %>%
      dplyr::mutate(
        switch_time = suppressWarnings(as.numeric(stringr::str_extract(.data$subexp.id, "^[^_]+"))),
        history_signature = stringr::str_replace(.data$subexp.id, "^[^_]+__", ""),
        .before = "subexp.id"
      )
  }

  out %>%
    dplyr::select(
      dplyr::any_of(c("subexp.id", "cohort", "switch_time", "history_signature")),
      "cohort_weight",
      "n_treated",
      "n_control",
      "control_mass"
    )
}


cbwsdid_order_cohort_output <- function(x, meta){
  if(identical(meta$design, "absorbing") && "cohort" %in% names(x)){
    return(
      x %>%
        dplyr::arrange(.data$cohort, dplyr::across(dplyr::any_of("et")))
    )
  }

  if((identical(meta$design, "switch01") || identical(meta$design, "switch10")) &&
     "switch_time" %in% names(x)){
    return(
      x %>%
        dplyr::arrange(.data$switch_time, .data$history_signature, dplyr::across(dplyr::any_of("et")))
    )
  }

  x %>%
    dplyr::arrange(.data$subexp.id, dplyr::across(dplyr::any_of("et")))
}


cbwsdid_qoi_linear_combo <- function(model, weights, level = 0.95){
  coef_all <- stats::coef(model)
  vcov_all <- stats::vcov(model)
  alpha <- 1 - level
  crit <- stats::qnorm(1 - alpha / 2)

  l_vec <- rep(0, length(coef_all))
  names(l_vec) <- names(coef_all)

  common_terms <- intersect(names(weights), names(l_vec))
  if(length(common_terms) > 0){
    l_vec[common_terms] <- weights[common_terms]
  }

  estimate <- sum(l_vec * coef_all)
  variance <- as.numeric(t(l_vec) %*% vcov_all %*% l_vec)
  std.error <- sqrt(variance)

  list(
    estimate = estimate,
    std.error = std.error,
    conf.low = estimate - crit * std.error,
    conf.high = estimate + crit * std.error
  )
}


cbwsdid_qoi_linear_combo_strict <- function(model, weights, level = 0.95){
  if(anyNA(weights)){
    return(list(
      estimate = NA_real_,
      std.error = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_
    ))
  }

  coef_all <- stats::coef(model)
  vcov_all <- stats::vcov(model)
  required <- names(weights)[abs(weights) > sqrt(.Machine$double.eps)]
  missing_terms <- setdiff(required, names(coef_all))

  if(length(missing_terms) > 0L ||
     any(!required %in% rownames(vcov_all)) ||
     any(!required %in% colnames(vcov_all))){
    return(list(
      estimate = NA_real_,
      std.error = NA_real_,
      conf.low = NA_real_,
      conf.high = NA_real_
    ))
  }

  cbwsdid_qoi_linear_combo(model = model, weights = weights, level = level)
}


cbwsdid_prepare_moderator_profiles <- function(meta, newdata = NULL){
  moderator_cols <- meta$moderation.spec$feature_name
  reference <- meta$moderation.reference

  if(is.null(reference) || !nrow(reference)){
    stop("No moderator reference data are stored in this `cbwsdid` object.")
  }

  if(is.null(newdata)){
    expanded_ref <- cbwsdid_expand_moderator_values(
      data = reference,
      reference = reference,
      moderator_cols = moderator_cols
    )

    if(length(expanded_ref$cols) == 0){
      expanded_mean <- tibble::tibble()
    } else {
      expanded_mean <- expanded_ref$data[expanded_ref$cols] %>%
        lapply(function(x) mean(x, na.rm = TRUE)) %>%
        tibble::as_tibble()
    }

    profile <- tibble::tibble(profile = "average")
    return(list(
      profile = profile,
      expanded = expanded_mean
    ))
  }

  newdata <- as.data.frame(newdata)
  defaults <- cbwsdid_default_moderator_values(reference, moderator_cols)

  for(col in moderator_cols){
    if(!col %in% names(newdata)){
      newdata[[col]] <- defaults[[col]]
    }
  }

  expanded <- cbwsdid_expand_moderator_values(
    data = newdata,
    reference = reference,
    moderator_cols = moderator_cols
  )

  profile <- tibble::as_tibble(newdata[moderator_cols])
  profile <- dplyr::mutate(profile, profile = paste0("profile_", dplyr::row_number()), .before = 1)

  list(
    profile = profile,
    expanded = tibble::as_tibble(expanded$data[expanded$cols])
  )
}


cbwsdid_default_moderator_values <- function(reference, moderator_cols){
  out <- list()

  for(col in moderator_cols){
    x <- reference[[col]]

    if(is.character(x) || is.factor(x)){
      x_chr <- as.character(x)
      x_chr <- x_chr[!is.na(x_chr)]
      if(length(x_chr) == 0){
        out[[col]] <- NA_character_
      } else {
        tab <- sort(table(x_chr), decreasing = TRUE)
        out[[col]] <- names(tab)[1]
      }
    } else {
      out[[col]] <- mean(x, na.rm = TRUE)
    }
  }

  out
}


cbwsdid_expand_moderator_values <- function(data, reference, moderator_cols){
  data <- as.data.frame(data)
  expanded_cols <- character()

  for(col in moderator_cols){
    ref_x <- reference[[col]]

    if(is.character(ref_x) || is.factor(ref_x)){
      ref_levels <- levels(as.factor(ref_x))
      data_x <- data[[col]]
      data_x <- factor(as.character(data_x), levels = ref_levels)

      if(any(is.na(data_x) & !is.na(data[[col]]))){
        bad_levels <- unique(as.character(data[[col]][is.na(data_x) & !is.na(data[[col]])]))
        stop(
          "Unknown level(s) supplied for moderator `", col, "`: ",
          paste(bad_levels, collapse = ", ")
        )
      }

      if(length(ref_levels) <= 1){
        next
      }

      for(level_i in ref_levels[-1]){
        col_i <- paste0(col, "__", make.names(level_i))
        data[[col_i]] <- as.integer(data_x == level_i)
        data[[col_i]][is.na(data[[col_i]])] <- 0
        expanded_cols <- c(expanded_cols, col_i)
      }
    } else {
      data[[col]] <- as.numeric(data[[col]])
      expanded_cols <- c(expanded_cols, col)
    }
  }

  list(
    data = data,
    cols = expanded_cols
  )
}


cbwsdid_qoi_moderated_bootstrap_terms <- function(out, type){
  if(!"profile" %in% names(out)){
    stop("Moderated QOI output is missing the `profile` column.")
  }

  if(identical(type, "dynamic")){
    return(paste0("dynamic:", out$profile, ":et=", out$et))
  }

  paste0(type, ":", out$profile)
}


cbwsdid_qoi_moderated_dynamic_combo <- function(model,
                                                et,
                                                moderator_values,
                                                level = 0.95,
                                                strict = FALSE){
  weights <- stats::setNames(1, paste0("et::", et, ":treated_sa"))

  if(ncol(moderator_values) > 0){
    for(col in names(moderator_values)){
      weights[paste0("et::", et, ":treated_x_", col)] <- as.numeric(moderator_values[[col]][1])
    }
  }

  if(isTRUE(strict)){
    return(cbwsdid_qoi_linear_combo_strict(
      model = model,
      weights = weights,
      level = level
    ))
  }

  cbwsdid_qoi_linear_combo(model = model, weights = weights, level = level)
}


cbwsdid_qoi_moderated_simple_combo <- function(model,
                                               post,
                                               moderator_values,
                                               level = 0.95,
                                               strict = FALSE){
  post <- sort(unique(post))
  weights <- c()

  for(et_i in post){
    weights[paste0("et::", et_i, ":treated_sa")] <- 1 / length(post)

    if(ncol(moderator_values) > 0){
      for(col in names(moderator_values)){
        weights[paste0("et::", et_i, ":treated_x_", col)] <- as.numeric(moderator_values[[col]][1]) / length(post)
      }
    }
  }

  if(isTRUE(strict)){
    return(cbwsdid_qoi_linear_combo_strict(
      model = model,
      weights = weights,
      level = level
    ))
  }

  cbwsdid_qoi_linear_combo(model = model, weights = weights, level = level)
}


cbwsdid_qoi_moderated_pooled_combo <- function(model,
                                               moderator_values,
                                               level = 0.95,
                                               strict = FALSE){
  weights <- stats::setNames(1, "treated_post")

  if(ncol(moderator_values) > 0){
    for(col in names(moderator_values)){
      weights[paste0("treated_post_x_", col)] <- as.numeric(moderator_values[[col]][1])
    }
  }

  if(isTRUE(strict)){
    return(cbwsdid_qoi_linear_combo_strict(
      model = model,
      weights = weights,
      level = level
    ))
  }

  cbwsdid_qoi_linear_combo(model = model, weights = weights, level = level)
}
