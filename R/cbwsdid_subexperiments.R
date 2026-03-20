#' Summarize Subexperiments Used by `cbwsdid`
#'
#' Produces compact summaries of the raw and refined donor pools across subexperiments used in a fitted [cbwsdid()] model.
#'
#' @param model `fixest` model object returned by [cbwsdid()].
#'
#' @return A list with:
#' - `overall` -- a summary-style tibble of design-level diagnostics;
#' - `by_subexp` -- a compact tibble with one row per subexperiment.
#'
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
#' cbwsdid_subexperiments(fit)
#' }
cbwsdid_subexperiments <- function(model){
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
  
  kappa <- meta$kappa
  if(is.null(kappa) || length(kappa) != 2){
    window_length <- NA_integer_
  } else {
    window_length <- as.integer(kappa[2] - kappa[1] + 1)
  }
  
  get_subexp_id <- function(detail){
    vals <- unique(c(detail$design.sample$subexp.id, detail$design.weights$subexp.id))
    vals <- vals[!is.na(vals)]
    
    if(length(vals) == 0){
      return(NA_character_)
    }
    
    as.character(vals[[1]])
  }
  
  by_subexp <- purrr::map_dfr(details, function(detail){
    design.sample <- detail$design.sample
    design.weights <- detail$design.weights
    
    tibble::tibble(
      subexp.id = get_subexp_id(detail),
      n_design_treated_ids = sum(design.sample$treated_sa == 1, na.rm = TRUE),
      n_design_control_ids = sum(design.sample$treated_sa == 0, na.rm = TRUE),
      n_refined_treated_ids = sum(design.weights$treated_sa == 1, na.rm = TRUE),
      n_refined_control_ids = sum(design.weights$treated_sa == 0, na.rm = TRUE)
    )
  })
  
  final_counts <- meta$weight.summary %>% 
    dplyr::mutate(
      subexp.id = as.character(.data$subexp.id),
      role = dplyr::if_else(.data$treated_sa == 1, "treated", "control")
    ) %>% 
    dplyr::select("subexp.id", "role", "n_ids", "raw_mass", "final_mass") %>% 
    tidyr::pivot_wider(
      names_from = "role",
      values_from = c("n_ids", "raw_mass", "final_mass"),
      names_sep = "_",
      values_fill = 0
    )
  
  by_subexp <- by_subexp %>% 
    dplyr::left_join(final_counts, by = "subexp.id") %>% 
    dplyr::mutate(
      dplyr::across(
        dplyr::starts_with("n_ids_"),
        ~ dplyr::coalesce(.x, 0L)
      ),
      dplyr::across(
        dplyr::starts_with("raw_mass_") | dplyr::starts_with("final_mass_"),
        ~ dplyr::coalesce(.x, 0)
      ),
      included_in_estimation = (.data$n_ids_treated > 0 & .data$n_ids_control > 0)
    )
  
  if(identical(meta$design, "switch01")){
    by_subexp <- by_subexp %>% 
      dplyr::mutate(
        switch_time = suppressWarnings(as.integer(stringr::str_extract(.data$subexp.id, "^[0-9]+"))),
        .before = 1
      )
  }
  
  if(identical(meta$design, "switch10")){
    by_subexp <- by_subexp %>% 
      dplyr::mutate(
        switch_time = suppressWarnings(as.integer(stringr::str_extract(.data$subexp.id, "^[0-9]+"))),
        .before = 1
      )
  }
  
  if(identical(meta$design, "absorbing")){
    by_subexp <- by_subexp %>% 
      dplyr::mutate(
        cohort = suppressWarnings(as.numeric(.data$subexp.id)),
        .before = 1
      )
  }
  
  by_subexp <- by_subexp %>% 
    dplyr::arrange(dplyr::desc(.data$included_in_estimation), .data$subexp.id)
  
  subexp_lookup <- by_subexp %>% 
    dplyr::transmute(subexp = dplyr::row_number(),
                     subexp.id = .data$subexp.id)
  
  by_subexp <- by_subexp %>% 
    dplyr::mutate(subexp = dplyr::row_number(), .before = 1)
  
  if(identical(meta$design, "switch01") || identical(meta$design, "switch10")){
    by_subexp <- by_subexp %>% 
      dplyr::transmute(
        subexp = .data$subexp,
        switch_time = .data$switch_time,
        included = .data$included_in_estimation,
        n_treated_raw = .data$n_design_treated_ids,
        n_control_raw = .data$n_design_control_ids,
        n_treated_refined = .data$n_ids_treated,
        n_control_refined = .data$n_ids_control,
        mass_control = .data$final_mass_control
      )
  }
  
  if(identical(meta$design, "absorbing")){
    by_subexp <- by_subexp %>% 
      dplyr::transmute(
        subexp = .data$subexp,
        cohort = .data$cohort,
        included = .data$included_in_estimation,
        n_treated_raw = .data$n_design_treated_ids,
        n_control_raw = .data$n_design_control_ids,
        n_treated_refined = .data$n_ids_treated,
        n_control_refined = .data$n_ids_control,
        mass_control = .data$final_mass_control
      )
  }
  
  n_candidates <- if(is.null(meta$candidates)) NA_integer_ else nrow(meta$candidates)
  n_included <- sum(by_subexp$included, na.rm = TRUE)
  n_nonempty_design <- sum((by_subexp$n_treated_raw + by_subexp$n_control_raw) > 0, na.rm = TRUE)
  n_tiny_design <- sum((by_subexp$n_treated_raw + by_subexp$n_control_raw) <= 2 &
                         (by_subexp$n_treated_raw + by_subexp$n_control_raw) > 0,
                       na.rm = TRUE)
  total_final_treated_ids <- sum(by_subexp$n_treated_refined, na.rm = TRUE)
  total_final_control_ids <- sum(by_subexp$n_control_refined, na.rm = TRUE)
  
  overall <- tibble::tibble(
    statistic = c(
      "design",
      "kappa_pre",
      "kappa_post",
      "window_length",
      "n_candidates",
      "n_included",
      "n_nonempty_design",
      "n_tiny_design",
      "total_treated_refined",
      "total_control_refined"
    ),
    value = c(
      as.character(meta$design),
      as.character(meta$kappa[1]),
      as.character(meta$kappa[2]),
      as.character(window_length),
      as.character(n_candidates),
      as.character(n_included),
      as.character(n_nonempty_design),
      as.character(n_tiny_design),
      as.character(total_final_treated_ids),
      as.character(total_final_control_ids)
    )
  )
  
  out <- list(
    call = match.call(),
    overall = overall,
    by_subexp = by_subexp
  )
  
  attr(out, "subexp_lookup") <- subexp_lookup
  
  class(out) <- "cbwsdid_subexperiments"
  return(out)
}
