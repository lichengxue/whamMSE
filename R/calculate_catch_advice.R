#' Calculate Catch Advice for Region- and Gear-Specific Fleets
#'
#' This function distributes total catch into region- and gear-specific fleets 
#' based on different weighting methods: "fleet_region", "fleet_gear", "fleet_combined", "fleet_catch", "index_equal", "index_gear", "multiple_index_equal", "multiple_index_gear", or user-specified weights.
#'
#' @param om List. The operating model containing input data, fleet, and index information.
#' @param advice Matrix or vector. Catch advice of dimension (n_years, n_fleets), 
#'   where each column represents a gear type. If advice is a vector, it represents total annual catch.
#' @param aggregate_catch_info List. Contains fleet allocation information:
#'   \itemize{
#'     \item $fleet_pointer - Vector indicating the gear type assignment for each fleet.
#'   }
#' @param aggregate_index_info List. Contains index region allocation information:
#'   \itemize{
#'     \item $index_pointer - Vector indicating the index region assignment.
#'   }
#' @param final_year Integer. The last year to consider when calculating historical catch weights.
#' @param catch_alloc List. Contains specifications for catch allocation:
#'   \itemize{
#'     \item $weight_type - Integer (1–4) indicating the type of weighting used.
#'       \itemize{
#'         \item 1 - Equal weighting across fleets.
#'         \item 2 - Weighting based on past catch data.
#'         \item 3 - Weighting based on past survey index data.
#'         \item 4 - User-specified weights.
#'       }
#'     \item $method - String. Specifies the catch allocation method:
#'       \itemize{
#'         \item "equal" - Total catch is equally allocated to all fleets across regions.
#'         \item "gear_equal" - Gear-specific total catch is equally allocated to fleets using the same gear type.
#'         \item "fleet_region" - Total catch is allocated to each region based on historical regional catch. Then, the regional catch is equally distributed among all fleets within the region.
#'         \item "fleet_gear" - Total catch is allocated to each gear type based on historical gear-specific catch. Then, the gear-specific catch is equally distributed among regions that use this gear.
#'         \item "fleet_combined" - Total catch is allocated based on historical regional and gear-specific catch.
#'         \item "fleet_catch" - Total catch is allocated based on historical fleet-specific catch within each region.
#'         \item "index_equal" - Total catch is allocated to each region based on historical survey catch. Then, the regional catch is distributed among fleets based on historical gear-specific catch.
#'         \item "index_gear" - Total catch is allocated based on historical region-specific survey catch and then distributed among fleets based on historical gear-specific catch.
#'         \item "multiple_index_equal" - Total catch is allocated to each region based on historical survey catch from multiple surveys. Then, the regional catch is equally distributed among fleets within the region.
#'         \item "multiple_index_gear" - Total catch is allocated to each region based on multiple historical survey-based catch estimates, then distributed among fleets within the region based on historical gear-specific catch.
#'         \item "user_defined_fleets" - Catch allocation is based on user-defined weights for each fleet.
#'         \item "user_defined_regions" - Catch is allocated based on user-defined weights for regions, then equally distributed among fleets within the region.
#'       }
#'     \item $user_weights - Numeric vector (n_regions or n_fleets). Optional. User-defined weights summing to 1 (use when weight_type = 4).
#'     \item $weight_years - Integer. Number of years to average for calculating historical catch weights.
#'     \item $survey_pointer - Integer/Vector. Specifies which survey index type to use for weighting (use when weight_type = 3).
#'   }
#'
#' @return A matrix of dimension (n_years, n_fleets) representing region-specific 
#'   and gear-specific catch advice.
#'
#' @export

calculate_catch_advice <- function(om, 
                                   advice, 
                                   aggregate_catch_info, 
                                   aggregate_index_info, 
                                   final_year,
                                   catch_alloc = list(weight_type = 1, method = "equal", user_weights = NULL, weight_years = 1, survey_pointer = 1)) {
  
  if (is.null(catch_alloc)) catch_alloc <- list(weight_type = 1, method = "equal", user_weights = NULL, weight_years = 1)
  if (is.null(catch_alloc$survey_pointer) && catch_alloc$weight_type == 3) {
    cat("\nsurvey_pointer must be specified when weight_type = 3!\n")
    catch_alloc$survey_pointer = 1
    cat("\nsurvey_pointer is set as 1.\n")
  }
  
  if (is.null(catch_alloc$user_weights) && catch_alloc$weight_type == 4) {
    cat("\ncatch_alloc$user_weights must be specified when weight_type = 4!\n")
    catch_alloc$user_weights = rep(1/length(fleet_regions),length(fleet_regions))
    cat("\nsurvey_pointer is set as equal weights.\n")
  }
  
  weight_type <- catch_alloc$weight_type
  method <- catch_alloc$method
  user_weights <- catch_alloc$user_weights
  
  fleet_regions <- om$input$data$fleet_regions
  
  # if(any(aggregate_index_info$index_pointer == 0)) warnings("Note: Some indices are excluded!")
  
  index_regions <- om$input$data$index_regions
  n_fleets <- length(fleet_regions)
  n_indices <- length(index_regions)
  n_years <- ifelse(is.matrix(advice), nrow(advice), 1)
  
  fleet_pointer <- if (!is.null(aggregate_catch_info$fleet_pointer)) aggregate_catch_info$fleet_pointer else rep(1, n_fleets)
  index_pointer <- if (!is.null(aggregate_index_info$index_pointer)) aggregate_index_info$index_pointer else rep(1, n_indices)
  
  if (weight_type == 1) { 
    if(!method %in% c("equal","gear_equal"))
      warning("method is not specified correctly!")
    weight_matrix <- matrix(1 / n_fleets, n_years, n_fleets, byrow = TRUE)
    if (method == "equal") {
      final_advice <- weight_matrix * if (is.matrix(advice)) rowSums(advice) else sum(advice) # completely equal
    } else if (method == "gear_equal") {
      if (is.vector(advice)) {
        # Initialize advice matrix
        advice_matrix <- matrix(0, nrow = n_years, ncol = n_fleets)
        
        # Loop through each region and assign advice equally across fleets in that region
        for (f in unique(fleet_pointer)) {
          fleets_in_region <- which(fleet_pointer == f)
          if (length(fleets_in_region) > 0) {
            advice_matrix[, fleets_in_region] <- advice[f] / length(fleets_in_region)
          }
        }
      } else {
        # Convert region-level matrix advice into fleet-level matrix
        advice_matrix <- matrix(0, nrow = n_years, ncol = n_fleets)
        
        for (f in unique(fleet_pointer)) {
          fleets_in_region <- which(fleet_pointer == f)
          if (length(fleets_in_region) > 0) {
            advice_matrix[, fleets_in_region] <- advice[, f] / length(fleets_in_region)
          }
        }
      } 
      final_advice <- advice_matrix
    }
  }
  
  if (weight_type == 2) {
    
    if(!method %in% c("fleet_region","fleet_gear","fleet_combined", "fleet_catch"))
      warning("method is not specified correctly!")
    
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
    
    # Handle both matrix and vector advice cases
    if (is.matrix(advice)) {
      final_advice <- weight_matrix * rowSums(advice)
    } else if (is.vector(advice)) {
      final_advice <- weight_matrix * sum(advice)
    }
  }
  
  if (weight_type == 3) {
    
    if(!method %in% c("index_equal","index_gear","multiple_index_equal","multiple_index_gear"))
      warning("method is not specified correctly!")
    
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
    } else if (method == "multiple_index_equal") {
      if (length(catch_alloc$survey_pointer) == 1)
        stop("survey_pointer must be > 1!")
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
    } else if (method == "multiple_index_gear") {
      if (length(catch_alloc$survey_pointer) == 1)
        stop("survey_pointer must be > 1!")
      
      gear_weights_matrix <- matrix(0, nrow = length(unique(fleet_pointer)), ncol = length(catch_alloc$survey_pointer))
      for (i in 1:length(catch_alloc$survey_pointer)) {
        survey_catch <- colMeans(om$input$data$agg_indices[which(om$years %in% avg_years), index_pointer == catch_alloc$survey_pointer[i], drop = FALSE])
        gear_weights_matrix[,i] <- survey_catch / sum(survey_catch)
      }
      region_weights <- rowMeans(gear_weights_matrix)  # Final region weights
      
      gear_catch <- colMeans(om$input$data$agg_catch[which(om$years %in% avg_years), , drop = FALSE])
      gear_weights <- tapply(gear_catch, fleet_pointer, sum) / sum(gear_catch)
      
      fleet_weights <- rep(0, n_fleets)
      for (region in unique(fleet_regions)) {
        fleets_in_region <- which(fleet_regions == region)
        if (length(fleets_in_region) > 0) {
          fleet_weights[fleets_in_region] <- region_weights[region] * gear_weights[fleet_pointer[fleets_in_region]]
        }
      }
    }
   
    weight_matrix <- matrix(fleet_weights, n_years, n_fleets, byrow = TRUE)

    # Apply weight matrix to catch advice
    if (is.matrix(advice)) {
      final_advice <- weight_matrix * rowSums(advice)
    } else {
      final_advice <- weight_matrix * sum(advice)
    }
  }

  if (weight_type == 4 && method %in% c("user_defined_fleets","user_defined_regions") && !is.null(user_weights)) {
    n_regions = om$input$data$n_regions
    if (sum(user_weights) != 1) stop("Sum of weights is not 1!")
    if (length(user_weights) == n_fleets) {
      weight_matrix <- matrix(user_weights, n_years, n_fleets, byrow = TRUE)
    } else if (length(user_weights) == n_regions) {
      fleet_weights <- rep(0, n_fleets)
      for (region in unique(fleet_regions)) {
        fleets_in_region <- which(fleet_regions == region)
        fleet_weights[fleets_in_region] <- user_weights[region] / length(fleets_in_region)
      }
      weight_matrix <- matrix(fleet_weights, n_years, n_fleets, byrow = TRUE)
    }
    
    # Handle both matrix and vector advice cases
    if (is.matrix(advice)) {
      final_advice <- weight_matrix * rowSums(advice)
    } else {
      final_advice <- weight_matrix * sum(advice)
    } 
  }
  
  # if (weight_type == 2) { 
  #   
  #   if(!method %in% c("fleet_region","fleet_gear","fleet_combined", "fleet_catch")) warning("method is not specified correctly!")
  #   
  #   # Weighting based on past catch data
  #   avg_years <- (final_year - catch_alloc$weight_years + 1):final_year
  #   catch <- colMeans(om$input$data$agg_catch[which(om$years %in% avg_years), , drop = FALSE])
  #   
  #   region_weights <- tapply(catch, fleet_regions, sum) / sum(catch)
  #   gear_weights <- tapply(catch, fleet_pointer, sum) / sum(catch)
  #   
  #   fleet_weights <- rep(0, n_fleets)
  #   if (method == "fleet_region") {
  #     for (region in unique(fleet_regions)) {
  #       fleets_in_region <- which(fleet_regions == region)
  #       fleet_weights[fleets_in_region] <- region_weights[region] / length(fleets_in_region)
  #     }
  #   } else if (method == "fleet_gear") {
  #     for (gear in unique(fleet_pointer)) {
  #       fleets_with_gear <- which(fleet_pointer == gear)
  #       fleet_weights[fleets_with_gear] <- gear_weights[gear] / length(fleets_with_gear)
  #     }
  #   } else if (method == "fleet_combined") {
  #     fleet_weights <- region_weights[fleet_regions] * gear_weights[fleet_pointer]
  #   } else if (method == "fleet_catch") {
  #     fleet_weights <- tapply(catch, 1:n_fleets, sum) / sum(catch)
  #   }
  #   
  #   weight_matrix <- matrix(fleet_weights, n_years, n_fleets, byrow = TRUE)
  #   final_advice <- weight_matrix * if (is.matrix(advice)) rowSums(advice) else advice
  # }
  # 
  # if (weight_type == 3) { 
  #   
  #   if(!method %in% c("index_equal","index_gear","multiple_index_equal")) warning("method is not specified correctly!")
  #   
  #   # Weighting based on past survey index data
  #   avg_years <- (final_year - catch_alloc$weight_years + 1):final_year
  #   survey_catch <- colMeans(om$input$data$agg_indices[which(om$years %in% avg_years), index_pointer == catch_alloc$survey_pointer, drop = FALSE])
  #   
  #   region_weights <- survey_catch / sum(survey_catch)
  # 
  #   gear_catch <- colMeans(om$input$data$agg_catch[which(om$years %in% avg_years), , drop = FALSE])
  #   
  #   gear_weights <- tapply(gear_catch, fleet_pointer, sum) / sum(gear_catch)
  #   
  #   fleet_weights <- rep(0, n_fleets)
  #   if (method == "index_equal") {
  #     for (region in unique(fleet_regions)) {
  #       fleets_in_region <- which(fleet_regions == region)
  #       fleet_weights[fleets_in_region] <- region_weights[region] / length(fleets_in_region)
  #     }
  #   } else if (method == "index_gear") {
  #     fleet_weights <- region_weights[fleet_regions] * gear_weights[fleet_pointer]
  #   } else if (method == "multiple_index_equal") {
  #     if (length(catch_alloc$survey_pointer) == 1) stop("survey_pointer must be > 1!")
  #     gear_weights_matrix <- matrix(0, nrow = length(unique(fleet_pointer)), ncol = length(catch_alloc$survey_pointer))
  #     for (i in 1:length(catch_alloc$survey_pointer)) {
  #       survey_catch <- colMeans(om$input$data$agg_indices[which(om$years %in% avg_years), index_pointer == catch_alloc$survey_pointer[i], drop = FALSE])
  #       gear_weights_matrix[,i] <- survey_catch / sum(survey_catch)
  #     }
  #     gear_weights <- rowMeans(gear_weights_matrix)
  #     for (region in unique(fleet_regions)) {
  #       fleets_in_region <- which(fleet_regions == region)
  #       fleet_weights[fleets_in_region] <- gear_weights[region] / length(fleets_in_region)
  #     }
  #   }
  #   
  #   weight_matrix <- matrix(fleet_weights, n_years, n_fleets, byrow = TRUE)
  #   final_advice <- weight_matrix * if (is.matrix(advice)) rowSums(advice) else advice
  # }
  # 
  # if (!is.null(user_weights)) {
  #   if (length(user_weights) == n_fleets) {
  #     weight_matrix <- matrix(user_weights, n_years, n_fleets, byrow = TRUE)
  #     final_advice <- weight_matrix * if (is.matrix(advice)) rowSums(advice) else advice
  #   }
  #   if (length(user_weights) == n_regions) {
  #     for (region in unique(fleet_regions)) {
  #       fleets_in_region <- which(fleet_regions == region)
  #       fleet_weights[fleets_in_region] <- user_weights[region] / length(fleets_in_region)
  #     }
  #     weight_matrix <- matrix(fleet_weights, n_years, n_fleets, byrow = TRUE)
  #     final_advice <- weight_matrix * if (is.matrix(advice)) rowSums(advice) else advice
  #   }
  # } 
  
  return(final_advice)
}
