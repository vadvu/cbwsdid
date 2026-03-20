#' Balance Diagnostics for `cbwsdid`
#'
#' Computes standardized mean difference tables for the refinement covariates used in [cbwsdid()]. Diagnostics are computed within subexperiments using the first-stage refinement weights rather than the final stacked weights.
#'
#' @param model `fixest` model object returned by [cbwsdid()].
#' @param include_exact Logical. If `TRUE`, include variables from `exact.formula` in the balance output.
#' @param which Which standardized mean difference table to return: `"both"`, `"adj"`, or `"un"`.
#' @param return_cobalt Logical. If `TRUE`, also return the underlying [cobalt::bal.tab()] objects.
#'
#' @return A list with compact balance tables for the pooled sample and for each
#'   subexperiment.
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
#'
#' cbwsdid_balance(fit)
#' }
cbwsdid_balance <- function(model,
                            include_exact = FALSE,
                            which = c("both", "adj", "un"),
                            return_cobalt = FALSE){
  which <- match.arg(which)
  
  if(!inherits(model, "fixest")){
    stop("`model` must be a fixest object returned by `cbwsdid()`.")
  }
  
  meta <- attr(model, "cbwsdid")
  if(is.null(meta)){
    stop("`model` does not contain `cbwsdid` metadata.")
  }
  
  details <- meta$refinement.details
  if(length(details) == 0){
    stop("No refinement details were stored in the model object.")
  }
  
  make_empty_table <- function(){
    tibble::tibble(et = numeric())
  }
  
  build_smd_table <- function(df, value_col){
    if(nrow(df) == 0){
      return(make_empty_table())
    }
    
    df %>% 
      dplyr::transmute(
        et = .data$period,
        variable = .data$variable,
        value = .data[[value_col]]
      ) %>% 
      dplyr::group_by(.data$et, .data$variable) %>% 
      dplyr::summarise(value = dplyr::first(.data$value), .groups = "drop") %>% 
      dplyr::arrange(.data$et, .data$variable) %>% 
      tidyr::pivot_wider(names_from = "variable", values_from = "value")
  }
  
  format_output_tables <- function(df){
    if(which == "adj"){
      return(build_smd_table(df, "Diff.Adj"))
    }
    
    if(which == "un"){
      return(build_smd_table(df, "Diff.Un"))
    }
    
    list(
      smd.un = build_smd_table(df, "Diff.Un"),
      smd.adj = build_smd_table(df, "Diff.Adj")
    )
  }
  
  balance_one_subexp <- function(detail){
    subexp.id <- if(length(detail$subexp.id) == 0) NA_character_ else as.character(detail$subexp.id[1])
    
    feature.spec <- detail$feature.spec$covariates
    if(include_exact){
      feature.spec <- dplyr::bind_rows(feature.spec, detail$feature.spec$exact)
    }
    
    feature.spec <- feature.spec %>% 
      dplyr::distinct(.data$feature_name, .keep_all = TRUE) %>% 
      dplyr::transmute(
        covariate = .data$feature_name,
        variable = .data$source_var,
        period = .data$et_lookup
      )
    
    design.sample <- detail$design.sample
    design.weights <- detail$design.weights
    
    empty_out <- list(
      subexp.id = subexp.id,
      n_treated = if(nrow(design.sample) == 0) 0 else sum(design.sample$treated_sa == 1),
      long = tibble::tibble(
        subexp.id = character(),
        n_treated = integer(),
        covariate = character(),
        variable = character(),
        period = numeric(),
        Diff.Un = numeric(),
        Diff.Adj = numeric()
      ),
      cobalt = NULL
    )
    
    if(nrow(feature.spec) == 0 || nrow(design.sample) == 0 || nrow(design.weights) == 0){
      return(empty_out)
    }
    
    bal.data <- dplyr::left_join(
      design.sample,
      design.weights,
      by = c("id", "subexp.id", "treated_sa")
    ) %>% 
      dplyr::filter(!is.na(.data$bsa))
    
    if(nrow(bal.data) == 0 ||
       sum(bal.data$treated_sa == 1) == 0 ||
       sum(bal.data$treated_sa == 0) == 0){
      return(empty_out)
    }
    
    control.weights <- bal.data$bsa[bal.data$treated_sa == 0]
    tol <- sqrt(.Machine$double.eps)
    
    if(length(control.weights) == 0 ||
       all(is.na(control.weights)) ||
       sum(control.weights, na.rm = TRUE) <= tol ||
       all(control.weights <= tol, na.rm = TRUE)){
      return(empty_out)
    }
    
    bal.formula <- stats::reformulate(
      termlabels = feature.spec$covariate,
      response = "treated_sa"
    )
    
    bal.obj <- cobalt::bal.tab(
      bal.formula,
      data = bal.data,
      weights = bal.data$bsa,
      estimand = "ATT",
      un = TRUE,
      quick = FALSE
    )
    
    bal.df <- as.data.frame(bal.obj$Balance) %>% 
      tibble::rownames_to_column(var = "covariate") %>% 
      tibble::as_tibble() %>% 
      dplyr::inner_join(feature.spec, by = "covariate") %>% 
      dplyr::mutate(
        subexp.id = subexp.id,
        n_treated = sum(design.sample$treated_sa == 1),
        .before = 1
      ) %>% 
      dplyr::select("subexp.id",
                    "n_treated",
                    "covariate",
                    "variable",
                    "period",
                    dplyr::any_of(c("Diff.Un", "Diff.Adj")))
    
    return(list(
      subexp.id = subexp.id,
      n_treated = sum(design.sample$treated_sa == 1),
      long = bal.df,
      cobalt = if(return_cobalt) bal.obj else NULL
    ))
  }
  
  balance.results <- lapply(details, balance_one_subexp)
  
  balance.long <- dplyr::bind_rows(purrr::map(balance.results, "long"))
  
  if(nrow(balance.long) == 0){
    out <- list(
      call = match.call(),
      refinement.method = meta$refinement.method,
      covs.formula = meta$covs.formula,
      exact.formula = meta$exact.formula,
      overall = format_output_tables(balance.long),
      by_subexp = list()
    )
    
    if(return_cobalt){
      out$cobalt <- purrr::map(balance.results, "cobalt")
    }
    
    class(out) <- "cbwsdid_balance"
    return(out)
  }
  
  overall.long <- balance.long %>% 
    dplyr::group_by(.data$period, .data$variable) %>% 
    dplyr::summarise(
      Diff.Un = stats::weighted.mean(.data$Diff.Un,
                                     w = .data$n_treated,
                                     na.rm = TRUE),
      Diff.Adj = stats::weighted.mean(.data$Diff.Adj,
                                      w = .data$n_treated,
                                      na.rm = TRUE),
      .groups = "drop"
    )
  
  subexp.ids <- unique(balance.long$subexp.id)
  by_subexp <- stats::setNames(
    lapply(subexp.ids, function(a){
      df.a <- balance.long %>% 
        dplyr::filter(.data$subexp.id == a)
      format_output_tables(df.a)
    }),
    nm = as.character(subexp.ids)
  )
  
  out <- list(
    call = match.call(),
    refinement.method = meta$refinement.method,
    covs.formula = meta$covs.formula,
    exact.formula = meta$exact.formula,
    overall = format_output_tables(overall.long),
    by_subexp = by_subexp
  )
  
  if(return_cobalt){
    out$cobalt <- purrr::map(balance.results, "cobalt")
  }
  
  class(out) <- "cbwsdid_balance"
  return(out)
}
