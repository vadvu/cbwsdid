#' Simulated Panel Data for `cbwsdid` Examples
#'
#' A simulated panel dataset for illustrating covariate-balanced weighted
#' stacked difference-in-differences with staggered absorbing treatment.
#'
#' @format A data frame with 6,500 rows and 7 variables:
#' \describe{
#'   \item{id}{Unit identifier.}
#'   \item{year}{Time period.}
#'   \item{outcome}{Observed outcome.}
#'   \item{adopt_year}{First treatment year; `NA` for never-treated units.}
#'   \item{D}{Binary treatment indicator.}
#'   \item{x1}{Continuous baseline covariate.}
#'   \item{x2}{Binary baseline covariate.}
#' }
#'
#' @details
#' Treatment timing is staggered and related to observed covariates. Untreated
#' potential outcome trends also depend on those covariates, so the dataset is
#' useful for demonstrating weighted stacked DID, covariate refinement, balance
#' diagnostics, and post-estimation quantities of interest.
#'
#' @source Simulated for examples accompanying Ustyuzhanin, Vadim (2026),
#' "Covariate-Balanced Weighted Stacked Difference-in-Differences", arXiv
#' preprint. \url{https://doi.org/10.48550/arXiv.2604.02293}.
#'
#' @examples
#' data(cbwsdid_sim)
#' head(cbwsdid_sim)
#'
#' @name cbwsdid_sim
#' @docType data
#' @keywords datasets
"cbwsdid_sim"
