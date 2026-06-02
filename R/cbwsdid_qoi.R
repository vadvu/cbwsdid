#' Quantities of Interest from `cbwsdid`
#'
#' Extracts aggregated quantities of interest from a fitted [cbwsdid()] model.
#'
#' @param model `fixest` model object from [cbwsdid()].
#' @param type Character. Quantity of interest to compute. `"dynamic"` returns an event-time table convenient for plotting. `"simple"` returns a single post-treatment ATT aggregated across selected event times. `"pooled"` returns the pooled ATT from the additional pooled regression. `"cohort"` returns cohort- or episode-level dynamic and simple effects. The difference between `"simple"` and `"pooled"` is in reference period: in the first it is `-1` average, while in the second it is an average of all pre-treatment periods.  
#' @param model_type Character. Which second-stage model to use. `"overall"` uses the baseline unmoderated models. `"moderated"` uses the moderated second-stage models stored in the metadata.
#' @param post Optional. Integer array of non-negative event times to include in the simple ATT aggregation. By default, all available post-treatment event times are used.
#' @param newdata Optional data frame describing moderator values for moderated quantities of interest. If `NULL`, moderated effects are averaged over the observed moderator distribution used in estimation.
#' @param level Confidence level for interval estimates.
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
                        level = 0.95){
  type <- match.arg(type)
  model_type <- match.arg(model_type)
  
  if(!inherits(model, "fixest")){
    stop("`model` must be a fixest object returned by `cbwsdid()`.")
  }
  
  meta <- attr(model, "cbwsdid")
  alpha <- 1 - level
  crit <- stats::qnorm(1 - alpha / 2)

  if(type == "cohort"){
    if(!identical(model_type, "overall")){
      stop("`type = \"cohort\"` is only available for `model_type = \"overall\"`.")
    }
    
    return(cbwsdid_qoi_cohort(model, post = post, level = level))
  }
  
  if(identical(model_type, "overall")){
    if(type == "dynamic"){
      return(cbwsdid_qoi_dynamic_overall(model, level = level))
    }
    
    if(type == "simple"){
      return(cbwsdid_qoi_simple_overall(model, post = post, level = level))
    }
    
    pooled_model <- meta$pooled_model
    if(is.null(pooled_model)){
      stop(
        "No pooled model is stored in this `cbwsdid` object. ",
        "Re-estimate with `pooled = TRUE`."
      )
    }
    
    pooled_est <- cbwsdid_qoi_linear_combo(
      pooled_model,
      weights = stats::setNames(1, "treated_post"),
      level = level
    )
    
    return(
      tibble::tibble(
        model_type = "overall",
        type = "pooled",
        estimate = pooled_est$estimate,
        std.error = pooled_est$std.error,
        conf.low = pooled_est$conf.low,
        conf.high = pooled_est$conf.high
      )
    )
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


cbwsdid_qoi_cohort <- function(model, post = NULL, level = 0.95){
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
  
  term.map <- model.data %>% 
    dplyr::distinct(.data$subexp.id, .data$et) %>% 
    dplyr::filter(.data$et != -1) %>% 
    dplyr::arrange(.data$subexp.id, .data$et) %>% 
    dplyr::mutate(term = paste0("cbwsdid_ce_", dplyr::row_number()))
  
  if(nrow(term.map) == 0){
    stop("No cohort-event coefficients are available outside the reference period `et = -1`.")
  }
  
  for(i in seq_len(nrow(term.map))){
    model.data[[term.map$term[i]]] <- as.integer(
      model.data$treated_sa == 1 &
        model.data$subexp.id == term.map$subexp.id[i] &
        model.data$et == term.map$et[i]
    )
  }
  
  rhs <- paste(term.map$term, collapse = " + ")
  cohort.formula <- stats::as.formula(
    paste("y ~", rhs, "| id^subexp.id + time^subexp.id")
  )
  
  cohort.model <- feols(
    data = model.data,
    cohort.formula,
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
    dplyr::left_join(cbwsdid_cohort_info(meta), by = "subexp.id") %>% 
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


cbwsdid_qoi_moderated_dynamic_combo <- function(model, et, moderator_values, level = 0.95){
  weights <- stats::setNames(1, paste0("et::", et, ":treated_sa"))
  
  if(ncol(moderator_values) > 0){
    for(col in names(moderator_values)){
      weights[paste0("et::", et, ":treated_x_", col)] <- as.numeric(moderator_values[[col]][1])
    }
  }
  
  cbwsdid_qoi_linear_combo(model = model, weights = weights, level = level)
}


cbwsdid_qoi_moderated_simple_combo <- function(model,
                                               post,
                                               moderator_values,
                                               level = 0.95){
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
  
  cbwsdid_qoi_linear_combo(model = model, weights = weights, level = level)
}


cbwsdid_qoi_moderated_pooled_combo <- function(model,
                                               moderator_values,
                                               level = 0.95){
  weights <- stats::setNames(1, "treated_post")
  
  if(ncol(moderator_values) > 0){
    for(col in names(moderator_values)){
      weights[paste0("treated_post_x_", col)] <- as.numeric(moderator_values[[col]][1])
    }
  }
  
  cbwsdid_qoi_linear_combo(model = model, weights = weights, level = level)
}
