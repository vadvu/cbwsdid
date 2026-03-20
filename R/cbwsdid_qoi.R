#' Quantities of Interest from `cbwsdid`
#'
#' Extracts aggregated quantities of interest from a fitted [cbwsdid()] model.
#'
#' @param model `fixest` model object from [cbwsdid()].
#' @param type Character. Quantity of interest to compute. `"simple"` returns a single post-treatment ATT aggregated across selected event times. `"dynamic"` returns an event-time table convenient for plotting.
#' @param post Optional. Integer array of non-negative event times to include in the simple ATT aggregation. By default, all available post-treatment event times are used.
#' @param level Confidence level for interval estimates.
#' @return A tibble. For `type = "simple"`, one row with the aggregated ATT and confidence interval. For `type = "dynamic"`, one row per event time.
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
#' cbwsdid_qoi(fit, type = "simple")
#' cbwsdid_qoi(fit, type = "dynamic")
#' }
cbwsdid_qoi <- function(model,
                        type = c("simple", "dynamic"),
                        post = NULL,
                        level = 0.95){
  type <- match.arg(type)
  
  if(!inherits(model, "fixest")){
    stop("`model` must be a fixest object returned by `cbwsdid()`.")
  }
  
  coef_all <- stats::coef(model)
  vcov_all <- stats::vcov(model)
  
  es_terms <- names(coef_all)[grepl("treated_sa", names(coef_all))]
  if(length(es_terms) == 0){
    stop("No event-study coefficients matching `treated_sa` were found.")
  }
  
  es_df <- tibble::tibble(
    term = es_terms,
    et = as.numeric(stringr::str_match(es_terms, "et::(-?[0-9]+):treated_sa")[, 2]),
    estimate = unname(coef_all[es_terms]),
    std.error = sqrt(diag(vcov_all[es_terms, es_terms, drop = FALSE]))
  ) %>% 
    dplyr::arrange(.data$et)
  
  alpha <- 1 - level
  crit <- stats::qnorm(1 - alpha / 2)
  
  if(type == "dynamic"){
    out <- es_df %>% 
      dplyr::mutate(
        conf.low = .data$estimate - crit * .data$std.error,
        conf.high = .data$estimate + crit * .data$std.error
      ) %>% 
      dplyr::select("et",
                    "estimate",
                    "std.error",
                    "conf.low",
                    "conf.high")
    
    return(out)
  }
  
  if(is.null(post)){
    post <- es_df$et[es_df$et >= 0]
  }
  
  post <- sort(unique(post))
  if(length(post) == 0){
    stop("No post-treatment event times were selected.")
  }
  
  simple_df <- es_df %>% 
    dplyr::filter(.data$et %in% post)
  
  if(nrow(simple_df) != length(post)){
    missing_post <- setdiff(post, simple_df$et)
    stop(
      "Some requested event times are not present in the model: ",
      paste(missing_post, collapse = ", ")
    )
  }
  
  weight_vec <- rep(1 / length(post), length(post))
  names(weight_vec) <- simple_df$term
  
  l_vec <- rep(0, length(coef_all))
  names(l_vec) <- names(coef_all)
  l_vec[names(weight_vec)] <- weight_vec
  
  estimate <- sum(l_vec * coef_all)
  variance <- as.numeric(t(l_vec) %*% vcov_all %*% l_vec)
  std.error <- sqrt(variance)
  
  tibble::tibble(
    type = "simple",
    periods = paste(post, collapse = ", "),
    n_periods = length(post),
    estimate = estimate,
    std.error = std.error,
    conf.low = estimate - crit * std.error,
    conf.high = estimate + crit * std.error
  )
}
