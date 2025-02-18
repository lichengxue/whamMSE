#' Calculate Catch Advice for Region- and Gear-Specific Fleets
#'
#' This function distributes total catch into region- and gear-specific fleets 
#' based on different weighting methods: `"fleet_region"`, `"fleet_gear"`, `"fleet_combined"`, or user-specified weights.
#'
#' @param om List. The operating model containing input data, fleet, and index information.
#' @param advice Matrix or vector. Catch advice of dimension (`n_years, n_gear_types`), 
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
#'         \item `"region"` - Uses regional catch totals to compute weights.
#'         \item `"gear"` - Uses gear-type catch totals to compute weights.
#'         \item `"combined"` - First allocates to gear-specific total catch, then splits into regions.
#'       }
#'     \item `$user_weights` - Numeric vector (`n_regions`). Optional. User-defined weights summing to 1.
#'     \item `$weight_years` - Integer. Number of years to average for calculating historical catch weights.
#'   }
#'
#' @return A matrix of dimension (`n_years, n_fleets`) representing region-specific 
#'   and gear-specific catch advice.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' om <- some_operating_model
#' advice <- matrix(rnorm(30), nrow = 10, ncol = 3) # Example catch advice for 10 years, 3 gear types
#' catch_alloc <- list(weight_type = 2, method = "region", weight_years = 5)
#' final_year <- 2025
#' 
#' catch_allocation <- calculate_catch_advice(om, advice, 
#'                                            aggregate_catch_info = list(fleet_pointer = c(1, 2, 3)), 
#'                                            aggregate_index_info = list(index_pointer = c(1, 2, 3)), 
#'                                            final_year = final_year, 
#'                                            catch_alloc = catch_alloc)
#' }

calculate_catch_advice <- function(om, 
                                   advice, 
                                   aggregate_catch_info, 
                                   aggregate_index_info, 
                                   final_year,
                                   catch_alloc = list(weight_type = 1, method = "region", user_weights = NULL, weight_years = 1, survey_pointer = 1)
                                   ) {
  
  if(is.null(catch_alloc)) catch_alloc = list(weight_type = 1, method = "region", user_weights = NULL, weight_years = 1)
  
  weight_type = catch_alloc$weight_type
  method = catch_alloc$method
  user_weights = catch_alloc$user_weights
  
  fleet_regions <- om$input$data$fleet_regions
  index_regions <- om$input$data$index_regions

  n_fleets <- length(fleet_regions) # Total number of fleets, but be cautious if one region is not fishing
  n_indices <- length(index_regions)
  n_years <- ifelse(is.matrix(advice),nrow(advice), length(advice))
  
  if (is.null(aggregate_catch_info$fleet_pointer)) {
    warnings("aggregate_catch_info$fleet_pointer is not specified!")
    fleet_pointer <- rep(1, n_fleets)
  } else {
    fleet_pointer <- aggregate_catch_info$fleet_pointer
  }
  
  if (is.null(aggregate_index_info$index_pointer)) {
    warnings("aggregate_index_info$index_pointer is not specified!")
    index_pointer <- rep(1, n_indices)
  } else {
    index_pointer <- aggregate_index_info$index_pointer
  }
  
  if (weight_type == 1) {
    
    weights <- rep(1 / n_fleets, n_fleets)
    weight_matrix <- matrix(weights, n_years, n_fleets, byrow = T)
    final_advice = matrix(0, nrow = n_years, ncol = n_fleets)
    
    if (is.matrix(advice)) {
      final_advice = weight_matrix*rowSums(advice)
    } else {
      final_advice = weight_matrix*advice
    }
  }
  
  if (weight_type == 2) {
    
    method = catch_alloc$method
    weight_years = catch_alloc$weight_years
    if (is.null(catch_alloc$weight_years)) weight_years = 1
    if (is.null(final_year)) final_year = tail(om$years,1)
    
    unique_regions <- unique(fleet_regions)   # Extract unique regions
    unique_gears <- unique(fleet_pointer)     # Extract unique gear types
    
    if (length(unique_regions) != om$input$data$n_regions) stop("Fleet is not operated in every region!")
    
    avg_years <- (final_year - weight_years + 1):final_year
    catch <- colMeans(om$input$data$agg_catch[which(om$years %in% avg_years), , drop = FALSE])
    
    # Sum total catch by region
    region_catch = NULL
    for (r in unique_regions) {
      region_catch[r] <- sum(catch[which(fleet_regions == unique_regions[r])])
    }
    
    # Sum total catch by gear type
    gear_catch = NULL
    for (g in unique_gears) {
      gear_catch[g] <- sum(catch[which(fleet_pointer == unique_gears[g])])
    }
    
    # regional weights
    region_weights <- region_catch / sum(region_catch)
    
    # gear-type weights
    gear_weights <- gear_catch / sum(gear_catch)  # Gear-based weight matrix (gear types)
    
    weight_matrix <- matrix(0, n_years, n_fleets)
    final_advice = matrix(0, nrow = n_years, ncol = n_fleets)
    
    if (!method %in% c("fleet_region","fleet_gear","fleet_combined")) {
      warnings("method should be one of the following region, fleet_gear, combined!\n")
      warnings("method is forced to be region")
      method = "fleet_region"
    }
    
    if (method == "fleet_region") {
      
      region_counts <- table(fleet_regions) 
      fleet_weights <- numeric(length(fleet_regions)) 
      
      for (region in names(region_counts)) {
        region_index <- as.numeric(region) 
        fleets_in_region <- which(fleet_regions == region_index) 
        num_fleets <- length(fleets_in_region) # Split the region weight equally among fleets in the region 
        fleet_weights[fleets_in_region] <- region_weights[region_index] / num_fleets 
      } 
      
      weight_matrix = matrix(fleet_weights, n_years, length(fleet_weights), byrow = TRUE)
      
      if (is.matrix(advice)) {
        final_advice = weight_matrix*rowSums(advice)
      } else {
        final_advice = weight_matrix*rowSums(as.matrix(advice))
      }
    }
    
    if (method == "fleet_gear") {
      
      gear_counts <- table(fleet_pointer) 
      fleet_weights <- numeric(length(gear_counts)) 
      
      for (gear in names(gear_counts)) {
        gear_index <- as.numeric(gear) 
        gears_in_region <- which(fleet_pointer == gear_index) 
        num_gears <- length(gears_in_region) # Split the region weight equally among fleets in the region 
        fleet_weights[gears_in_region] <- gear_weights[gear_index] / num_gears 
      } 
      
      weight_matrix = matrix(fleet_weights, n_years, length(fleet_weights), byrow = TRUE)
      
      if (is.matrix(advice)) {
        final_advice = weight_matrix*rowSums(advice)
      } else {
        final_advice = weight_matrix*rowSums(as.matrix(advice))
      }
    }
    
    if (method == "fleet_combined") {
      
      fleet_weights = NULL
      for (r in 1:unique(length(fleet_regions))){
        pointer2 = fleet_pointer[which(fleet_regions == r)]
        weights = region_weights[r]*prop.table(gear_weights[pointer2])
        fleet_weights = c(fleet_weights,weights)
      }
      
      weight_matrix = matrix(fleet_weights, n_years, length(fleet_weights), byrow = TRUE)
      
      if (is.matrix(advice)) {
        final_advice = weight_matrix*rowSums(advice)
      } else {
        final_advice = weight_matrix*rowSums(as.matrix(advice))
      }
    }
  }
  
  if (weight_type == 3) {
    
    method = catch_alloc$method
    weight_years = catch_alloc$weight_years
    if (is.null(catch_alloc$weight_years)) weight_years = 1
    if (is.null(final_year)) final_year = tail(om$years,1)
    
    unique_regions <- unique(fleet_regions)   # Extract unique regions
    unique_gears <- unique(fleet_pointer)     # Extract unique gear types
    
    if (length(unique_regions) != om$input$data$n_regions) stop("Fleet is not operated in every region!")
    
    avg_years <- (final_year - weight_years + 1):final_year
    
    if (is.null(catch_alloc$survey_pointer)) stop("Specify `survey_pointer` to determine which survey to use for weighting.")
    survey_pointer <- catch_alloc$survey_pointer # User-defined survey type to focus on
    
    index_subset <- which(index_pointer == survey_pointer) # Get index positions that match the survey type
    if (length(index_subset) == 0) stop("No matching survey indices found for the given `survey_pointer`.")
    
    survey_catch <- colMeans(om$input$data$agg_indices[which(om$years %in% avg_years), index_subset, drop = FALSE])
    
    if (length(unique(index_regions[index_subset])) != length(unique_regions)) stop("Survey type does not cover all regions.")
    
    region_weights <- tapply(survey_catch, index_regions[index_subset], sum) / sum(survey_catch)
    
    fleet_weights <- numeric(length(fleet_regions))
    
    if (method == "index_equal") {
      
      # Equal assignment: divide region weight by the number of fleets in that region
      region_counts <- table(fleet_regions)
      for (region in names(region_counts)) {
        region_index <- as.numeric(region)
        fleets_in_region <- which(fleet_regions == region_index)
        num_fleets <- length(fleets_in_region)
        fleet_weights[fleets_in_region] <- region_weights[region_index] / num_fleets
      }
      
    } else if (method == "index_gear") {
      
      catch <- colMeans(om$input$data$agg_catch[which(om$years %in% avg_years), , drop = FALSE])
      
      gear_catch = NULL
      for (g in unique_gears) {
        gear_catch[g] <- sum(catch[which(fleet_pointer == unique_gears[g])])
      }
      
      gear_weights <- gear_catch / sum(gear_catch) 
      
      region_counts <- table(fleet_regions)
      for (region in unique_regions) {
        fleets_in_region <- which(fleet_regions == region)
        num_fleets <- length(fleets_in_region)
        
        gear_types_in_region <- fleet_pointer[fleets_in_region] 
        
        # Compute fleet weights using regional survey weights and gear weights
        for (fleet in fleets_in_region) {
          gear_type <- fleet_pointer[fleet] # Get gear type
          fleet_weights[fleet] <- region_weights[region] * gear_weights[gear_type]
        }
      }
    }
    
    # Construct weight matrix
    weight_matrix <- matrix(fleet_weights, n_years, length(fleet_weights), byrow = TRUE)
    
    # Compute final catch advice
    if (is.matrix(advice)) {
      final_advice <- weight_matrix * rowSums(advice)
    } else {
      final_advice <- weight_matrix * sum(advice)
    }
  }

  if (weight_type == 4) {
    
    user_weights <- catch_alloc$user_weights
    if (!is.null(user_weights) & length(user_weights) == n_regions & sum(user_weights) == 1) {
      weights <- user_weights
    } else {
      warnings("user_weights is not specified correctly!\nEqual weights are assumed.")
      user_weights == 1/n_regions
    }
    
    weight_matrix <- matrix(weights, n_years, n_regions, byrow = T)
    
    if (is.matrix(advice)) {
      final_advice = weight_matrix*rowSums(advice)
    } else {
      final_advice = weight_matrix*sum(advice)
    }
    
  }
  
  return(final_advice)
}
