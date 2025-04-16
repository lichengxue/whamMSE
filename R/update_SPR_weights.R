#' Update SPR Weights for Spatially Explicit Assessment Models
#'
#' Updates the weights used to calculate spawning potential ratio (SPR)-based biological reference points
#' in spatially explicit assessment models. This function sets region- or stock-specific weights
#' for averaging metrics like weight-at-age and maturity-at-age in the reference point calculations.
#'
#' @param em_input List. An input list for the estimation model (typically created by \code{make_em_input}).
#' Must include \code{em_input$data$agg_catch}, \code{em_input$data$fleet_regions}, and other standard data structures.
#' @param method Character. Specifies how weights are assigned to regions or stocks. Options include:
#'   \describe{
#'     \item{\code{"equal"}}{Assigns equal weights across all stocks or regions (default).}
#'     \item{\code{"fleet_region"}}{Assigns weights based on average catch per region over the last \code{weight_years}. Requires complete fleet coverage across regions.}
#'     \item{\code{"index_region"}}{Assigns weights based on average index values over the last \code{weight_years}. Requires \code{index_pointer}.}
#'   }
#' @param weight_years Integer. Number of years at the end of the time series over which to average catch or index values (default = 1).
#' @param index_pointer Integer. Index number used when \code{method = "index_region"}.
#'
#' @return An updated version of \code{em_input}, with SPR weighting information added to \code{em_input$data$SPR_weights}
#' and \code{em_input$data$SPR_weight_type}.
#'
#' @details
#' This function is intended for use only in spatially explicit estimation models. It modifies the weighting applied
#' to regional data when computing spawning potential ratio (SPR) and associated reference points. The default method
#' applies equal weights, while alternative methods use historical catch or index data to determine weights dynamically.
#'
#' @export
#' 

update_SPR_weights <- function(em_input, method = "equal", weight_years = 1, index_pointer = NULL) {
  em_input$data$SPR_weight_type = 1
  n_stocks = em_input$data$n_stocks
  n_regions = em_input$data$n_regions
  final_year = length(em_input$years)
  avg_years <- (final_year - weight_years + 1):final_year
  if (method == "equal") {
    em_input$data$SPR_weights = rep(1/n_stocks, n_stocks)
  }
  if (method == "fleet_region") {
    fleet_regions = em_input$data$fleet_regions
    if (length(unique(fleet_regions)) != n_regions) stop("Current version does not allow missing fleet in a region!")
    catch <- colMeans(em_input$data$agg_catch[avg_years, , drop = FALSE])
    em_input$data$SPR_weights <- tapply(catch, fleet_regions, sum) / sum(catch)
  }
  if (method == "index_region") {
    if(is.null(index_pointer)) stop("index_pointer must be specified!")
    catch <- colMeans(om$input$data$agg_indices[avg_years, index_pointer, drop = FALSE])
    em_input$data$SPR_weights <- tapply(catch, 1:n_regions, sum) / sum(catch)
  }
  return(em_input)
}