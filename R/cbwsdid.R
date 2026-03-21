
#' Covariate-Balanced Weighted Stacked DID
#'
#' Fits covariate-balanced weighted stacked difference-in-differences on a panel dataset and returns a `fixest` model object.
#'
#' @param data Panel dataframe.
#' @param y Character. Name of the outcome variable.
#' @param d Character. Name of the binary treatment indicator.
#' @param id character array. A length-2 character vector giving the unit and time identifiers. For example, `c("id", "time")`
#' @param kappa Numeric array. A length-2 numeric vector with the event-time window `c(kappa_pre, kappa_post)`. For example, `c(-2, 2)`
#' @param design Design type. One of `"absorbing"` (staggered adoption), `"switch01"`, or `"switch10"`.
#' @param post_path Only used if `design != "absorbing"`. One of `"stable"` or
#'   `"relaxed"`. Under `"stable"`, treated episodes must remain treated
#'   throughout the post-treatment window. Under `"relaxed"`, treated episodes
#'   are defined only by treatment onset at event time zero and may subsequently
#'   reverse within the event window.
#' @param history.length Only used if design != `"absorbing"`. Length of treatment history used to define admissible switch-on or switch-off episode types. Required for switch designs.
#' @param first_switch_only Only used if design != `"absorbing"`. Logical. If `TRUE`, keep only the first admissible switch episode for each unit in switch designs.
#' @param refinement.method Refinement method. One of `"none"` (for weighted stacked did without refinement), `"matchit"`, or `"weightit"`.
#' @param covs.formula Optional formula describing refinement covariates. For example `~ lag(y, 1:4) + lag(x, 1:4)`.
#' @param exact.formula Optional formula describing exact matching or exact stratification variables. for example, `~ region`.
#' @param refinement.args Optional named list of arguments passed to the chosen refinement backend. See [MatchIt::matchit()] ot [WeightIt::weightit()]. 
#' @param allow_treated_drop Logical. If `TRUE`, allow the refinement stage to drop treated units or episodes, which changes the estimand to a matched-sample or overlap-trimmed ATT.
#'
#' @details
#' The `design` argument determines how primary subexperiments are constructed.
#' `"absorbing"` implements the standard staggered-adoption design. `"switch01"`
#' constructs subexperiments around repeated `0 -> 1` switch-on episodes.
#' `"switch10"` applies the same switch-on logic after internally recoding the
#' treatment as `1 - d`, which makes it convenient to study repeated `1 -> 0`
#' switch-off episodes with the same pipeline.
#'
#' For switch-based designs, `post_path` controls the post-treatment path
#' restriction on treated episodes. `"stable"` requires treated episodes to
#' remain treated throughout the event window and therefore targets effects of
#' sustained exposure within that window. `"relaxed"` requires only treatment
#' onset at event time zero and allows later treatment reversals within the
#' window, so post-treatment coefficients should be interpreted as effects of
#' treatment onset averaged over subsequent realized treatment paths.
#'
#' @return A `fixest` model object from [fixest::feols()] with additional `cbwsdid` metadata stored in `attr(x, "cbwsdid")`.
#' @export
#'
#' @examples
#' \dontrun{
#' fit <- cbwsdid(
#'   data = panel_df,
#'   y = "y",
#'   d = "d",
#'   id = c("unit", "time"),
#'   kappa = c(-4, 4),
#'   refinement.method = "none"
#' )
#' }
cbwsdid <- function(data = data, 
                    y = "y", 
                    d = "d", 
                    id = c("country", "year"),
                    kappa = c(-10, 10),
                    design = c("absorbing", "switch01", "switch10"),
                    post_path = c("stable", "relaxed"),
                    history.length = NULL,
                    first_switch_only = FALSE,
                    refinement.method = c("none", "matchit", "weightit"),
                    covs.formula = NULL,
                    exact.formula = NULL,
                    refinement.args = NULL,
                    allow_treated_drop = FALSE
                    ){
  design <- match.arg(design)
  post_path <- match.arg(post_path)
  refinement.method <- match.arg(refinement.method)
  refinement.vars <- unique(c(all.vars(covs.formula), all.vars(exact.formula)))
  refinement.source.vars <- setdiff(refinement.vars, c("id", "time", "y", "d"))
  keep.vars <- unique(c(id, y, d, refinement.source.vars))
  
  db <- data %>% 
    select(all_of(keep.vars)) %>% 
    rename(id = all_of(id[1]),
           time = all_of(id[2]),
           y = all_of(y),
           d = all_of(d)) %>% 
    arrange(id, time)
  
  if(design == "switch10"){
    db <- db %>% 
      mutate(d = 1 - d)
  }
  
  # Checks: td
  
  if(design == "absorbing"){
    db <- db %>% 
      group_by(id) %>% 
      mutate(treated = ifelse(sum(d) == 0, 0, 1),
             A = ifelse(treated == 0, Inf, min(ifelse(d == 1, time, Inf))),
             ) %>% 
      ungroup()
    
    candidates <- sort(unique(db$A[is.finite(db$A)]))
    
    subexp.results <- lapply(
      candidates,
      FUN = function(x){
        sub.i <- subexp.row.construct(db = db,
                                      a = x,
                                      kappa = kappa)
        subexp.refine(subexp = sub.i,
                      covs.formula = covs.formula,
                      refinement.method = refinement.method,
                      exact.formula = exact.formula,
                      refinement.args = refinement.args,
                      allow_treated_drop = allow_treated_drop)
      }
    )
  } else {
    if(is.null(history.length) || length(history.length) != 1 || history.length < 1){
      stop("For `design = \"switch01\"` or `design = \"switch10\"`, `history.length` must be a positive integer.")
    }
    
    history.length <- as.integer(history.length)
    
    candidates <- switch01_candidates(
      db = db,
      kappa = kappa,
      history.length = history.length,
      first_switch_only = first_switch_only,
      post_path = post_path
    )
    
    if(nrow(candidates) == 0){
      stop("No admissible 0->1 episode types were found for the requested design and window.")
    }
    
    subexp.results <- lapply(
      seq_len(nrow(candidates)),
      FUN = function(i){
        cand.i <- candidates[i, ]
        
        sub.i <- subexp.switch01.construct(
          db = db,
          tau = cand.i$tau,
          h = unlist(cand.i$h),
          kappa = kappa,
          history.length = history.length,
          first_switch_only = first_switch_only,
          post_path = post_path
        )
        
        subexp.refine(subexp = sub.i,
                      covs.formula = covs.formula,
                      refinement.method = refinement.method,
                      exact.formula = exact.formula,
                      refinement.args = refinement.args,
                      allow_treated_drop = allow_treated_drop)
      }
    )
  }
  
  subexp.all <- subexp.results %>% 
    purrr::map("subexp.refined") %>% 
    list_rbind()
  
  if(nrow(subexp.all) == 0){
    stop("All admissible sub-experiments were dropped before estimation.")
  }
  
  subexp.all.weighted <- compute_stack_weights(subexp.all)
  
  subexp.all.data <- subexp.all.weighted$subexp.weighted

  
  model <- feols(data = subexp.all.data, 
                 y ~ i(et, treated_sa, ref = -1) | id^subexp.id + time^subexp.id,
                 cluster = ~id, 
                 weights = ~Qsa
  )
  
  refinement.details <- subexp.results %>% 
    purrr::map(function(x){
      list(
        subexp.id = unique(x$design.sample$subexp.id),
        backend = x$backend,
        design.sample = x$design.sample,
        design.weights = x$design.weights,
        refinement.stat = x$refinement.stat,
        feature.spec = x$feature.spec
      )
    })
  
  attr(model, "cbwsdid") <- list(
    call = match.call(),
    kappa = kappa,
    design = design,
    treatment_recode = identical(design, "switch10"),
    post_path = post_path,
    history.length = history.length,
    first_switch_only = first_switch_only,
    candidates = candidates,
    refinement.method = refinement.method,
    covs.formula = covs.formula,
    exact.formula = exact.formula,
    refinement.args = refinement.args,
    allow_treated_drop = allow_treated_drop,
    design.weights = subexp.all.weighted$design.weights,
    cohort.weights = subexp.all.weighted$cohort.weights,
    weight.summary = subexp.all.weighted$weight.summary,
    refinement.details = refinement.details
  )
  
  return(model)
}




subexp.row.construct <- function(db, a, kappa){
  event.t.tot <- kappa[2] - kappa[1] + 1

  subexp <- db %>% 
    filter(A == a | A > (a + kappa[2])) %>% 
    mutate(et = time - a) %>% 
    filter(dplyr::between(et, kappa[1], kappa[2])) %>% 
    add_count(id, name = "n_obs_window") %>% 
    filter(n_obs_window == event.t.tot) %>% 
    mutate(treated_sa = as.integer(A == a)) %>% 
    mutate(subexp.id = a)
  
  cohort.size <- subexp %>% 
    summarise(Nda = n_distinct(id[treated_sa == 1]),
              Nca = n_distinct(id[treated_sa == 0]))
  
  subexp <- subexp %>% 
    mutate(Nda = cohort.size$Nda,
           Nca = cohort.size$Nca) %>% 
    select(-n_obs_window)
  
  return(subexp)
}


switch01_candidates <- function(db,
                                kappa = c(-4, 2),
                                history.length = 4,
                                first_switch_only = FALSE,
                                post_path = c("stable", "relaxed")){
  post_path <- match.arg(post_path)

  if(history.length < 1){
    stop("`history.length` must be a positive integer.")
  }
  
  history.cols <- paste0("history_", seq_len(history.length))
  post.cols <- paste0("post_", 0:kappa[2])
  history.lags <- history.length:1
  
  history.exprs <- purrr::map(history.lags, ~ rlang::expr(dplyr::lag(.data$d, n = !!.x)))
  names(history.exprs) <- history.cols
  
  post.exprs <- purrr::map(0:kappa[2], ~ rlang::expr(dplyr::lead(.data$d, n = !!.x)))
  names(post.exprs) <- post.cols
  
  anchor <- db %>% 
    group_by(id) %>% 
    arrange(time, .by_group = TRUE) %>% 
    mutate(
      !!!history.exprs,
      !!!post.exprs,
      switch01 = dplyr::lag(.data$d, n = 1) == 0 & .data$d == 1,
      switch01_number = cumsum(dplyr::coalesce(switch01, FALSE))
    ) %>% 
    ungroup() %>% 
    filter(
      .data$switch01,
      if_all(all_of(history.cols), ~ !is.na(.x)),
      if_all(all_of(post.cols), ~ !is.na(.x))
    )

  if(identical(post_path, "stable")){
    anchor <- anchor %>%
      filter(if_all(all_of(post.cols), ~ .x == 1))
  }
  
  if(first_switch_only){
    anchor <- anchor %>% 
      filter(.data$switch01_number == 1)
  }
  
  anchor %>% 
    rowwise() %>% 
    mutate(
      h = list(as.integer(c_across(all_of(history.cols)))),
      history_signature = paste0(unlist(.data$h), collapse = "")
    ) %>% 
    ungroup() %>% 
    distinct(time, history_signature, .keep_all = TRUE) %>% 
    transmute(
      tau = .data$time,
      history_signature = .data$history_signature,
      h = .data$h
    ) %>% 
    arrange(.data$tau, .data$history_signature)
}


subexp.switch01.construct <- function(db,
                                      tau,
                                      h,
                                      kappa,
                                      history.length = length(h),
                                      first_switch_only = FALSE,
                                      post_path = c("stable", "relaxed")){
  post_path <- match.arg(post_path)

  if(is.null(h) || length(h) == 0){
    stop("`h` must be a non-empty 0/1 vector describing recent treatment history.")
  }
  
  h <- as.integer(h)
  
  if(any(!h %in% c(0L, 1L))){
    stop("`h` must contain only 0/1 values.")
  }
  
  if(length(h) != history.length){
    stop("`history.length` must match `length(h)`.")
  }
  
  if(kappa[2] < 0){
    stop("For the switch-on design, `kappa[2]` must be non-negative.")
  }
  
  if(h[history.length] != 0L){
    stop(
      "For a 0->1 switch-on design, the most recent lag in `h` must be 0, ",
      "because treated episodes satisfy D_(tau-1) = 0."
    )
  }
  
  history.cols <- paste0("history_", seq_len(history.length))
  history.lags <- history.length:1
  post.cols <- paste0("post_", 0:kappa[2])
  support.start <- min(kappa[1], -history.length)
  support.width <- kappa[2] - support.start + 1
  history.signature <- paste(h, collapse = "")
  subexp.label <- paste0(tau, "__", history.signature)
  history.exprs <- purrr::map(history.lags, ~ rlang::expr(dplyr::lag(.data$d, n = !!.x)))
  names(history.exprs) <- history.cols
  post.exprs <- purrr::map(0:kappa[2], ~ rlang::expr(dplyr::lead(.data$d, n = !!.x)))
  names(post.exprs) <- post.cols
  
  db.switch <- db %>% 
    group_by(id) %>% 
    arrange(time, .by_group = TRUE) %>% 
    mutate(
      !!!history.exprs,
      !!!post.exprs,
      switch01 = dplyr::lag(.data$d, n = 1) == 0 & .data$d == 1,
      switch01_number = cumsum(dplyr::coalesce(switch01, FALSE))
    ) %>% 
    ungroup()
  
  support.ids <- db.switch %>% 
    filter(dplyr::between(time, tau + support.start, tau + kappa[2])) %>% 
    add_count(id, name = "n_support_obs") %>% 
    filter(n_support_obs == support.width) %>% 
    distinct(id)
  
  anchor.sample <- db.switch %>% 
    filter(time == tau) %>% 
    semi_join(support.ids, by = "id")
  
  if(nrow(anchor.sample) == 0){
    return(
      db %>% 
        slice(0) %>% 
        mutate(
          et = numeric(),
          treated_sa = integer(),
          subexp.id = character(),
          switch_time = numeric(),
          history_signature = character(),
          episode_role = character(),
          Nda = integer(),
          Nca = integer()
        )
    )
  }
  
  history.matrix <- anchor.sample %>% 
    select(all_of(history.cols)) %>% 
    as.matrix()
  
  history.match <- apply(history.matrix, 1, function(x){
    all(!is.na(x) & x == h)
  })
  
  post.matrix <- anchor.sample %>% 
    select(all_of(post.cols)) %>% 
    as.matrix()
  
  stable.on <- apply(post.matrix, 1, function(x){
    all(!is.na(x) & x == 1)
  })
  
  stable.off <- apply(post.matrix, 1, function(x){
    all(!is.na(x) & x == 0)
  })
  
  treated_episode <- if(identical(post_path, "stable")){
    history.match & anchor.sample$switch01 & stable.on
  } else {
    history.match & anchor.sample$switch01
  }

  anchor.sample <- anchor.sample %>% 
    mutate(
      history_match = history.match,
      stable_on = stable.on,
      stable_off = stable.off,
      treated_episode = treated_episode,
      control_episode = history_match & d == 0 & stable_off
    )
  
  if(first_switch_only){
    anchor.sample <- anchor.sample %>% 
      mutate(
        treated_episode = treated_episode & switch01_number == 1
      )
  }
  
  eligible.ids <- anchor.sample %>% 
    filter(treated_episode | control_episode) %>% 
    transmute(
      id,
      treated_sa = as.integer(treated_episode),
      episode_role = dplyr::if_else(treated_episode, "treated", "control")
    )
  
  subexp <- db.switch %>% 
    semi_join(eligible.ids, by = "id") %>% 
    mutate(et = time - tau) %>% 
    filter(dplyr::between(et, kappa[1], kappa[2])) %>% 
    inner_join(eligible.ids, by = "id") %>% 
    mutate(
      subexp.id = subexp.label,
      switch_time = tau,
      history_signature = history.signature
    )
  
  cohort.size <- subexp %>% 
    summarise(
      Nda = n_distinct(id[treated_sa == 1]),
      Nca = n_distinct(id[treated_sa == 0])
    )
  
  subexp <- subexp %>% 
    mutate(
      Nda = cohort.size$Nda,
      Nca = cohort.size$Nca
    ) %>% 
    select(-any_of(c("history_match",
                     "stable_on",
                     "stable_off",
                     "treated_episode",
                     "control_episode",
                     "switch01",
                     "switch01_number")))
  
  return(subexp)
}



subexp.refine <- function(subexp,
                          covs.formula = NULL,
                          refinement.method = c("none", "matchit", "weightit"),
                          exact.formula = NULL,
                          refinement.args = NULL,
                          allow_treated_drop = FALSE){
  
  refinement.method <- match.arg(refinement.method)
  if(is.null(refinement.args)){
    refinement.args <- switch(
      refinement.method,
      matchit = list(
        method = "nearest",
        distance = "mahalanobis",
        ratio = 1,
        replace = FALSE
      ),
      weightit = list(
        method = "glm"
      ),
      list()
    )
  }
  
  # ---------------------------------------------------------------------------
  # 1. Parse user formulas in a PanelMatch style but without I()
  # ---------------------------------------------------------------------------
  #
  # Supported syntax for now:
  #   ~ lag(var, 0:4) + lag(y, 1:4) + z
  #
  # Interpretation inside sub-experiment a:
  #   lag(var, 0) = value at et = 0
  #   lag(var, 1) = value at et = -1
  #   lag(var, 2) = value at et = -2
  #
  # A bare variable z is interpreted as its value at et = -1.
  # ---------------------------------------------------------------------------
  split_formula_terms <- function(expr){
    if(is.null(expr)){
      return(list())
    }
    if(rlang::is_call(expr, "+")){
      return(c(split_formula_terms(expr[[2]]), split_formula_terms(expr[[3]])))
    }
    if(rlang::is_call(expr, "I")){
      return(split_formula_terms(expr[[2]]))
    }
    list(expr)
  }
  
  parse_feature_formula <- function(formula.obj){
    if(is.null(formula.obj)){
      return(tibble(source_var = character(),
                    lag_k = integer(),
                    et_lookup = integer(),
                    feature_name = character(),
                    term_type = character()))
    }
    
    rhs <- rlang::f_rhs(formula.obj)
    raw.terms <- split_formula_terms(rhs)
    
    feature.spec <- purrr::map_dfr(raw.terms, function(term){
      if(is.numeric(term) && length(term) == 1 && term == 1){
        return(tibble(source_var = character(),
                      lag_k = integer(),
                      et_lookup = integer(),
                      feature_name = character(),
                      term_type = character()))
      }
      
      if(is.symbol(term)){
        var.name <- as.character(term)
        return(tibble(source_var = var.name,
                      lag_k = 0L,
                      et_lookup = -1L,
                      feature_name = var.name,
                      term_type = "bare"))
      }
      
      if(rlang::is_call(term, "lag")){
        var.name <- as.character(term[[2]])
        lag.values <- rlang::eval_bare(term[[3]], env = baseenv()) %>% 
          as.integer()
        
        return(tibble(source_var = var.name,
                      lag_k = lag.values,
                      et_lookup = -lag.values,
                      feature_name = paste0(var.name, "_l", lag.values),
                      term_type = "lag"))
      }
      
      stop("Unsupported term in formula. Use terms like `lag(y, 1:4)` or bare variables.")
    })
    
    feature.spec %>% 
      distinct(feature_name, .keep_all = TRUE)
  }
  
  cov.spec <- parse_feature_formula(covs.formula)
  exact.spec <- parse_feature_formula(exact.formula)
  all.spec <- bind_rows(cov.spec, exact.spec) %>% 
    distinct(feature_name, .keep_all = TRUE)
  
  if(any(all.spec$term_type == "lag" & all.spec$lag_k == 0L)){
    warning(
      "By including 0 in lag(), you asked refinement to use post-treatment information (et = 0). ",
      "This is a brutal thing. I do not appreciate it."
    )
  }
  
  
  # ---------------------------------------------------------------------------
  # 2. Build one-row-per-unit design sample
  # ---------------------------------------------------------------------------
  design.sample <- subexp %>% 
    filter(et == -1) %>% 
    select(id, subexp.id, treated_sa) %>% 
    distinct()
  
  if(nrow(all.spec) > 0){
    feature.tables <- purrr::map(seq_len(nrow(all.spec)), function(i){
      spec.i <- all.spec[i, ]
      
      tmp <- subexp %>% 
        filter(et == spec.i$et_lookup) %>% 
        transmute(id,
                  subexp.id,
                  feature_value = .data[[spec.i$source_var]]) %>% 
        distinct()
      
      names(tmp)[names(tmp) == "feature_value"] <- spec.i$feature_name
      tmp
    })
    
    design.sample <- purrr::reduce(
      feature.tables,
      .f = left_join,
      .init = design.sample,
      by = c("id", "subexp.id")
    )
  }
  
  if(nrow(all.spec) > 0){
    design.sample <- design.sample %>% 
      filter(stats::complete.cases(across(all_of(all.spec$feature_name))))
  }
  
  empty.design.weights <- design.sample %>% 
    slice(0) %>% 
    transmute(id, subexp.id, treated_sa, bsa = numeric())
  
  if(nrow(design.sample) == 0 ||
     !any(design.sample$treated_sa == 1) ||
     !any(design.sample$treated_sa == 0)){
    subexp.refined <- subexp %>% 
      inner_join(empty.design.weights, by = c("id", "subexp.id", "treated_sa"))
    
    return(list(
      backend = refinement.method,
      design.sample = design.sample,
      feature.spec = list(covariates = cov.spec,
                          exact = exact.spec),
      design.weights = empty.design.weights,
      refinement.stat = tibble(treated_sa = integer(),
                               n_ids = integer(),
                               weight_mass = numeric()),
      subexp.refined = subexp.refined
    ))
  }
  
  cov.cols <- cov.spec$feature_name
  exact.cols <- exact.spec$feature_name
  
  build_exact_or_unweighted <- function(design.sample,
                                        exact.cols,
                                        allow_treated_drop){
    if(length(exact.cols) == 0){
      return(
        design.sample %>% 
          transmute(id, subexp.id, treated_sa, bsa = 1)
      )
    }
    
    design.sample.ex <- design.sample %>% 
      tidyr::unite("stratum", all_of(exact.cols), sep = "___", remove = FALSE)
    
    strata.stat <- design.sample.ex %>% 
      count(stratum, treated_sa, name = "n") %>% 
      tidyr::pivot_wider(names_from = treated_sa,
                         values_from = n,
                         names_prefix = "n_",
                         values_fill = 0) %>% 
      rename(n_control = n_0,
             n_treated = n_1)
    
    unsupported.treated.strata <- strata.stat %>% 
      filter(n_treated > 0, n_control == 0) %>% 
      pull(stratum)
    
    if(length(unsupported.treated.strata) > 0 && !allow_treated_drop){
      stop(
        "Exact refinement would drop treated units because some treated strata have no controls: ",
        paste(utils::head(unsupported.treated.strata, 5), collapse = ", "),
        ". Set `allow_treated_drop = TRUE` if you really want that."
      )
    }
    
    if(length(unsupported.treated.strata) > 0 && allow_treated_drop){
      warning(
        "Exact refinement dropped treated units in unsupported strata. ",
        "The resulting estimand is no longer the full treated-cohort ATT."
      )
    }
    
    strata.stat <- strata.stat %>% 
      filter(n_control > 0, n_treated > 0)
    
    design.sample.ex %>% 
      inner_join(strata.stat, by = "stratum") %>% 
      transmute(id,
                subexp.id,
                treated_sa,
                bsa = ifelse(treated_sa == 1, 1, n_treated / n_control))
  }
  
  if(length(cov.cols) == 0){
    cov.cols.use <- cov.cols
  } else {
    varying.cols <- cov.cols[purrr::map_lgl(cov.cols, function(v){
      x <- design.sample[[v]]
      x <- x[!is.na(x)]
      dplyr::n_distinct(x) > 1
    })]
    
    cov.cols.use <- varying.cols
    
    if(length(cov.cols.use) > 1 &&
       all(purrr::map_lgl(design.sample[cov.cols.use], is.numeric))){
      X <- as.matrix(design.sample[cov.cols.use])
      qr.X <- qr(X)
      max.rank <- min(qr.X$rank, nrow(X) - 1)
      
      if(max.rank <= 0){
        cov.cols.use <- character()
      } else {
        cov.cols.use <- cov.cols.use[qr.X$pivot[seq_len(max.rank)]]
      }
    }
  }
  
  
  # ---------------------------------------------------------------------------
  # 3. Estimate raw design weights b_sa
  # ---------------------------------------------------------------------------

  if(refinement.method == "none"){
    design.weights <- build_exact_or_unweighted(
      design.sample = design.sample,
      exact.cols = exact.cols,
      allow_treated_drop = allow_treated_drop
    )
    
  } else if (refinement.method == "matchit"){
    matchit.distance <- refinement.args$distance %||% "mahalanobis"
    
    if(length(cov.cols.use) == 0){
      warning(
        "No varying covariates remained in this sub-experiment after screening. ",
        "Falling back to exact / unweighted refinement."
      )
      
      design.weights <- build_exact_or_unweighted(
        design.sample = design.sample,
        exact.cols = exact.cols,
        allow_treated_drop = allow_treated_drop
      )
      
    } else if (identical(matchit.distance, "mahalanobis") && nrow(design.sample) <= 2){
      warning(
        "Mahalanobis refinement is ill-defined in this sub-experiment because only ",
        nrow(design.sample),
        " units remain in the design sample. Falling back to exact / unweighted refinement."
      )
      
      design.weights <- build_exact_or_unweighted(
        design.sample = design.sample,
        exact.cols = exact.cols,
        allow_treated_drop = allow_treated_drop
      )
      
    } else {
      match.formula <- stats::as.formula(
        paste("treated_sa ~", paste(cov.cols.use, collapse = " + "))
      )
      
      exact.constraint <- NULL
      if(length(exact.cols) > 0){
        exact.constraint <- stats::as.formula(
          paste("~", paste(exact.cols, collapse = " + "))
        )
      }
      
      forbidden.args <- c("formula", "data", "estimand", "exact")
      if(any(names(refinement.args) %in% forbidden.args)){
        stop(
          "`refinement.args` cannot contain: ",
          paste(intersect(names(refinement.args), forbidden.args), collapse = ", "),
          ". These are controlled internally by `subexp.refine()`."
        )
      }
      
      match.call.args <- c(
        list(
          formula = match.formula,
          data = design.sample,
          estimand = "ATT",
          exact = exact.constraint
        ),
        refinement.args
      )
      
      match.out <- do.call(MatchIt::matchit, match.call.args)
      
      matched.data <- MatchIt::match.data(match.out) %>% 
        filter(weights > 0) %>% 
        transmute(id, subexp.id, treated_sa, bsa = weights)
      
      treated.original <- design.sample %>% 
        filter(treated_sa == 1) %>% 
        pull(id) %>% 
        unique()
      
      treated.matched <- matched.data %>% 
        filter(treated_sa == 1) %>% 
        pull(id) %>% 
        unique()
      
      treated.dropped <- setdiff(treated.original, treated.matched)
      
      if(length(treated.dropped) > 0 && !allow_treated_drop){
        stop(
          "MatchIt refinement dropped treated units. This changes the estimand away from the full treated-cohort ATT. ",
          "Set `allow_treated_drop = TRUE` if you want to allow a matched-sample ATT."
        )
      }
      
      if(length(treated.dropped) > 0 && allow_treated_drop){
        warning(
          "MatchIt refinement dropped treated units. The resulting estimand is a matched-sample / overlap-trimmed ATT, ",
          "not the full treated-cohort ATT."
        )
      }
      
      design.weights <- matched.data
    }
    
  } else if (refinement.method == "weightit"){
    
    if(length(cov.cols.use) == 0){
      warning(
        "No varying covariates remained in this sub-experiment after screening. ",
        "Falling back to exact / unweighted refinement."
      )
      
      design.weights <- build_exact_or_unweighted(
        design.sample = design.sample,
        exact.cols = exact.cols,
        allow_treated_drop = allow_treated_drop
      )
      
    } else {
      forbidden.args <- c("formula", "data", "estimand")
      if(any(names(refinement.args) %in% forbidden.args)){
        stop(
          "`refinement.args` cannot contain: ",
          paste(intersect(names(refinement.args), forbidden.args), collapse = ", "),
          ". These are controlled internally by `subexp.refine()`."
        )
      }
      
      weight.formula <- stats::as.formula(
        paste("treated_sa ~", paste(cov.cols.use, collapse = " + "))
      )
      
      if(length(exact.cols) == 0){
        design.sample.ex <- design.sample %>% 
          mutate(stratum = "__all__")
        
        strata.stat <- design.sample.ex %>% 
          summarise(stratum = "__all__",
                    n_control = sum(treated_sa == 0),
                    n_treated = sum(treated_sa == 1))
      } else {
        design.sample.ex <- design.sample %>% 
          tidyr::unite("stratum", all_of(exact.cols), sep = "___", remove = FALSE)
        
        strata.stat <- design.sample.ex %>% 
          count(stratum, treated_sa, name = "n") %>% 
          tidyr::pivot_wider(names_from = treated_sa,
                             values_from = n,
                             names_prefix = "n_",
                             values_fill = 0) %>% 
          rename(n_control = n_0,
                 n_treated = n_1)
      }
      
      unsupported.treated.strata <- strata.stat %>% 
        filter(n_treated > 0, n_control == 0) %>% 
        pull(stratum)
      
      if(length(unsupported.treated.strata) > 0 && !allow_treated_drop){
        stop(
          "WeightIt refinement would drop treated units because some treated strata have no controls: ",
          paste(utils::head(unsupported.treated.strata, 5), collapse = ", "),
          ". Set `allow_treated_drop = TRUE` if you really want that."
        )
      }
      
      if(length(unsupported.treated.strata) > 0 && allow_treated_drop){
        warning(
          "WeightIt refinement dropped treated units in unsupported strata. ",
          "The resulting estimand is no longer the full treated-cohort ATT."
        )
      }
      
      supported.strata <- strata.stat %>% 
        filter(n_control > 0, n_treated > 0) %>% 
        pull(stratum)
      
      if(length(supported.strata) == 0){
        design.weights <- empty.design.weights
      } else {
        design.weights <- purrr::map_dfr(supported.strata, function(stratum.i){
          data.i <- design.sample.ex %>% 
            filter(stratum == stratum.i)
          
          weight.call.args <- c(
            list(
              formula = weight.formula,
              data = data.i,
              estimand = "ATT"
            ),
            refinement.args
          )
          
          weight.out <- do.call(WeightIt::weightit, weight.call.args)
          
          tibble(id = data.i$id,
                 subexp.id = data.i$subexp.id,
                 treated_sa = data.i$treated_sa,
                 bsa = weight.out$weights)
        })
      }
    }
    
  }
  
  
  # ---------------------------------------------------------------------------
  # 4. Merge design weights back to the stacked sub-experiment
  # ---------------------------------------------------------------------------
  subexp.refined <- subexp %>% 
    inner_join(design.weights, by = c("id", "subexp.id", "treated_sa"))
  
  refinement.stat <- design.weights %>% 
    group_by(treated_sa) %>% 
    summarise(n_ids = n_distinct(id),
              weight_mass = sum(bsa),
              .groups = "drop")
  
  return(list(
    backend = refinement.method,
    design.sample = design.sample,
    feature.spec = list(covariates = cov.spec,
                        exact = exact.spec),
    design.weights = design.weights,
    refinement.stat = refinement.stat,
    subexp.refined = subexp.refined
  ))
}


compute_stack_weights <- function(subexp.all,
                                  b_var = "bsa"){
  
  # ---------------------------------------------------------------------------
  # 1. Collapse to one row per (id, subexp.id)
  # ---------------------------------------------------------------------------
  #
  # Final stacked weights are objects at the unit-within-subexperiment level.
  # We therefore first strip the event-time duplication out of the stacked data.
  # ---------------------------------------------------------------------------
  design.weights <- subexp.all %>% 
    select(any_of(c("id", "subexp.id", "treated_sa", b_var))) %>% 
    distinct()
  
  if(!(b_var %in% names(design.weights))){
    design.weights <- design.weights %>% 
      mutate(bsa = 1)
  } else {
    design.weights <- design.weights %>% 
      rename(bsa = all_of(b_var))
  }
  
  
  # ---------------------------------------------------------------------------
  # 2. Cohort-level totals
  # ---------------------------------------------------------------------------
  #
  # Nda = number of treated ids in subexperiment a
  # Nca = number of control ids in subexperiment a
  # Bca = sum of raw refinement weights over controls in subexperiment a
  # ---------------------------------------------------------------------------
  cohort.weights <- design.weights %>% 
    group_by(subexp.id) %>% 
    summarise(
      Nda = n_distinct(id[treated_sa == 1]),
      Nca = n_distinct(id[treated_sa == 0]),
      Bca = sum(bsa[treated_sa == 0], na.rm = TRUE),
      .groups = "drop"
    )
  
  Nd <- sum(cohort.weights$Nda)
  Nc <- sum(cohort.weights$Nca)
  Bc <- sum(cohort.weights$Bca)
  
  
  # ---------------------------------------------------------------------------
  # 3. Final stacked weight
  # ---------------------------------------------------------------------------
  #
  # Final weight always has the covariate-balanced weighted stacked DID form:
  #   treated units get weight 1
  #   controls get bsa times the corrective factor
  #
  # If refinement.method = "none", then bsa = 1 and the formula collapses to
  # the usual weighted stacked DID.
  # ---------------------------------------------------------------------------
  design.weights <- design.weights %>% 
    left_join(cohort.weights, by = "subexp.id") %>% 
    mutate(stack_weight = ifelse(treated_sa == 1,
                                 1,
                                 bsa * (Nda / Nd) / (Bca / Bc)))
  
  
  # ---------------------------------------------------------------------------
  # 4. Merge weights back to every stacked row
  # ---------------------------------------------------------------------------
  subexp.weighted <- subexp.all %>% 
    select(-any_of(c("stack_weight", "Qsa"))) %>% 
    left_join(
      design.weights %>% 
        select(id, subexp.id, treated_sa, bsa, stack_weight),
      by = c("id", "subexp.id", "treated_sa")
    ) %>% 
    mutate(Qsa = stack_weight)
  
  weight.summary <- design.weights %>% 
    group_by(subexp.id, treated_sa) %>% 
    summarise(
      n_ids = n_distinct(id),
      raw_mass = sum(bsa, na.rm = TRUE),
      final_mass = sum(stack_weight, na.rm = TRUE),
      .groups = "drop"
    )
  
  return(list(
    subexp.weighted = subexp.weighted,
    design.weights = design.weights,
    cohort.weights = cohort.weights,
    totals = tibble(Nd = Nd, Nc = Nc, Bc = Bc),
    weight.summary = weight.summary
  ))
}
