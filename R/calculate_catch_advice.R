#' Calculate Catch Advice for Region- and Gear-Specific Fleets
#'
#' This function distributes total catch into region- and gear-specific fleets 
#' based on different weighting methods: `"fleet_region"`, `"fleet_gear"`, `"fleet_combined"`, `"fleet_catch"`, `"index_equal"`, `"index_gear"`, `"index_multiple"`, or user-specified weights.
#'
#' @param om List. The operating model containing input data, fleet, and index information.
#' @param advice Matrix or vector. Catch advice of dimension (`n_years, n_fleets`), 
#'   where each column represents a gear type. If `advice` is a vector, it represents total annual catch.
#' @param aggregate_catch_info List. Contains fleet allocation information:
#'   \itemize{
#'     \item `$fleet_pointer` - Vector indicating the gear type assignment for each fleet.
#'   }
#' @param aggregate_index_info List. Contains index region allocation information:
#'   \itemize{
#'     \item `$index_pointer` - Vector indicating the index region assignment.
#'   }
#' @param final_year Integer. The last year to consider when calculating historical catch weights.
#' @param catch_alloc List. Contains specifications for catch allocation:
#'   \itemize{
#'     \item `$weight_type` - Integer (1–4) indicating the type of weighting used.
#'       \itemize{
#'         \item `1` - Equal weighting across fleets.
#'         \item `2` - Weighting based on past catch data.
#'         \item `3` - Weighting based on past survey index data.
#'         \item `4` - User-specified weights.
#'       }
#'     \item `$method` - String. Specifies the catch allocation method:
#'       \itemize{
#'         \item `"fleet_region"` - Uses regional catch totals to compute weights.
#'         \item `"fleet_gear"` - Uses gear-type catch totals to compute weights.
#'         \item `"fleet_combined"` - First allocates to gear-specific total catch, then splits into regions.
#'         \item `"fleet_catch"` - Uses fleet-specific catch to compute weights.
#'         \item `"index_equal"` - Uses survey-based regional weighting, then assigns equally among fleets in the same region.
#'         \item `"index_gear"` - Uses survey-based regional weighting, then assigns based on gear-specific weights.
#'         \item `"index_multiple"` - Uses multiple survey-based regional weighting, then assigns equally among fleets in the same region.
#'         \item `"user_defined"` - Uses manually specified weights for each region or each fleet.
#'       }
#'     \item `$user_weights` - Numeric vector (`n_regions` or `n_fleets`). Optional. User-defined weights summing to 1.
#'     \item `$weight_years` - Integer. Number of years to average for calculating historical catch weights.
#'     \item `$survey_pointer` - Integer/Vector. Specifies which survey index type to use for weighting.
#'   }
#'
#' @return A matrix of dimension (`n_years, n_fleets`) representing region-specific 
#'   and gear-specific catch advice.
#'
#' @export
#'

calculate_catch_advice <- function(om, 
                                   advice, 
                                   aggregate_catch_info, 
                                   aggregate_index_info, 
                                   final_year,
                                   catch_alloc = list(weight_type = 1, method = "fleet_region", user_weights = NULL, weight_years = 1, survey_pointer = 1)) {
  
  if (is.null(catch_alloc)) catch_alloc <- list(weight_type = 1, method = "fleet_region", user_weights = NULL, weight_years = 1)
  
  weight_type <- catch_alloc$weight_type
  method <- catch_alloc$method
  user_weights <- catch_alloc$user_weights
  
  fleet_regions <- om$input$data$fleet_regions
  
  # if(any(aggregate_index_info$index_pointer == 0)) warnings("Note: Some indices are excluded!")
  
  index_regions <- om$input$data$index_regions
  n_fleets <- length(fleet_regions)
  n_indices <- length(index_regions)
  n_years <- ifelse(is.matrix(advice), nrow(advice), length(advice))
  
  fleet_pointer <- if (!is.null(aggregate_catch_info$fleet_pointer)) aggregate_catch_info$fleet_pointer else rep(1, n_fleets)
  index_pointer <- if (!is.null(aggregate_index_info$index_pointer)) aggregate_index_info$index_pointer else rep(1, n_indices)
  
  if (weight_type == 1) { 
    # Equal weighting across fleets
    weight_matrix <- matrix(1 / n_fleets, n_years, n_fleets, byrow = TRUE)
    final_advice <- weight_matrix * if (is.matrix(advice)) rowSums(advice) else advice
  }
  
  if (weight_type == 2) { 
    
    if(!method %in% c("fleet_region","fleet_gear","fleet_combined", "fleet_catch")) warning("method is not specified correctly!")
    
    # Weighting based on past catch data
    avg_years <- (final_year - catch_alloc$weight_years + 1):final_year
    catch <- colMeans(om$input$data$agg_catch[which(om$years %in% avg_years), , drop = FALSE])
    
    region_weights <- tapply(catch, fleet_regions, sum) / sum(catch)
    gear_weights <- tapply(catch, fleet_pointer, sum) / sum(catch)
    
    fleet_weights <- rep(0, n_fleets)
    if (method == "fleet_region") {
      for (region in unique(fleet_regions)) {
        fleets_in_region <- which(fleet_regions == region)
        fleet_weights[fleets_in_region] <- region_weights[region] / length(fleets_in_region)
      }
    } else if (method == "fleet_gear") {
      for (gear in unique(fleet_pointer)) {
        fleets_with_gear <- which(fleet_pointer == gear)
        fleet_weights[fleets_with_gear] <- gear_weights[gear] / length(fleets_with_gear)
      }
    } else if (method == "fleet_combined") {
      fleet_weights <- region_weights[fleet_regions] * gear_weights[fleet_pointer]
    } else if (method == "fleet_catch") {
      fleet_weights <- tapply(catch, 1:n_fleets, sum) / sum(catch)
    }
    
    weight_matrix <- matrix(fleet_weights, n_years, n_fleets, byrow = TRUE)
    final_advice <- weight_matrix * if (is.matrix(advice)) rowSums(advice) else advice
  }
  
  if (weight_type == 3) { 
    
    if(!method %in% c("index_equal","index_gear","index_multiple")) warning("method is not specified correctly!")
    
    # Weighting based on past survey index data
    avg_years <- (final_year - catch_alloc$weight_years + 1):final_year
    survey_catch <- colMeans(om$input$data$agg_indices[which(om$years %in% avg_years), index_pointer == catch_alloc$survey_pointer, drop = FALSE])
    
    region_weights <- survey_catch / sum(survey_catch)

    gear_catch <- colMeans(om$input$data$agg_catch[which(om$years %in% avg_years), , drop = FALSE])
    
    gear_weights <- tapply(gear_catch, fleet_pointer, sum) / sum(gear_catch)
    
    fleet_weights <- rep(0, n_fleets)
    if (method == "index_equal") {
      for (region in unique(fleet_regions)) {
        fleets_in_region <- which(fleet_regions == region)
        fleet_weights[fleets_in_region] <- region_weights[region] / length(fleets_in_region)
      }
    } else if (method == "index_gear") {
      fleet_weights <- region_weights[fleet_regions] * gear_weights[fleet_pointer]
    } else if (method == "index_multiple") {
      if (length(catch_alloc$survey_pointer) == 1) stop("survey_pointer must be > 1!")
      gear_weights_matrix <- matrix(0, nrow = length(unique(fleet_pointer)), ncol = length(catch_alloc$survey_pointer))
      for (i in 1:length(catch_alloc$survey_pointer)) {
        survey_catch <- colMeans(om$input$data$agg_indices[which(om$years %in% avg_years), index_pointer == catch_alloc$survey_pointer[i], drop = FALSE])
        gear_weights_matrix[,i] <- survey_catch / sum(survey_catch)
      }
      gear_weights <- rowMeans(gear_weights_matrix)
      for (region in unique(fleet_regions)) {
        fleets_in_region <- which(fleet_regions == region)
        fleet_weights[fleets_in_region] <- gear_weights[region] / length(fleets_in_region)
      }
    }
    
    weight_matrix <- matrix(fleet_weights, n_years, n_fleets, byrow = TRUE)
    final_advice <- weight_matrix * if (is.matrix(advice)) rowSums(advice) else advice
  }
  
  if (!is.null(user_weights)) {
    if (length(user_weights) == n_fleets) {
      weight_matrix <- matrix(user_weights, n_years, n_fleets, byrow = TRUE)
      final_advice <- weight_matrix * if (is.matrix(advice)) rowSums(advice) else advice
    }
    if (length(user_weights) == n_regions) {
      for (region in unique(fleet_regions)) {
        fleets_in_region <- which(fleet_regions == region)
        fleet_weights[fleets_in_region] <- user_weights[region] / length(fleets_in_region)
      }
      weight_matrix <- matrix(fleet_weights, n_years, n_fleets, byrow = TRUE)
      final_advice <- weight_matrix * if (is.matrix(advice)) rowSums(advice) else advice
    }
  } 
  
  return(final_advice)
}
