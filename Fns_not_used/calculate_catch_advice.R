#' Calculate Catch Advice for Region- and Gear-Specific Fleets
#'
#' This function distributes total catch into region- and gear-specific fleets
#' based on different weighting methods: "regional", "gear", or "combined".
#'
#' @param advice A matrix of dimension (n_years, n_gear_types), where each column represents a gear type.
#' @param fleet_regions A vector indicating the region assignment for each fleet.
#' @param fleet_pointer A vector indicating the gear type assignment for each fleet.
#' @param method A string specifying the weighting method:
#'   - `"regional"`: Uses regional catch totals to compute weights.
#'   - `"gear"`: Uses gear-type catch totals to compute weights.
#'   - `"combined"`: First allocates to gear-specific total catch, then splits into regions.
#'
#' @return A matrix of dimension (n_years, n_fleets) representing region-specific and gear-specific catch advice.
#' @export
calculate_catch_advice <- function(om, 
                                   advice, 
                                   aggregate_catch_info, 
                                   aggregate_index_info, 
                                   catch_alloc = list(weight_type = 1, method = "region", user_weights = NULL, weight_years = 1, final_year = NULL)
                                   ) {
  
  if(is.null(catch_alloc)) catch_alloc = list(weight_type = 1, method = "region", user_weights = NULL, weight_years = 1, final_year = NULL)
  
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
    final_year = catch_alloc$final_year
    if (is.null(catch_alloc$weight_years)) weight_years = 1
    if (is.null(catch_alloc$final_year)) final_year = tail(om$years,1)
    
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
    
    if (!method %in% c("region","gear","combined")) {
      warnings("method should be one of the following region, gear, combined!\n")
      warnings("method is forced to be region")
      method = "region"
    }
    
    if (method == "region") {
      
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
    
    if (method == "gear") {
      
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
    
    if (method == "combined") {
      
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
    
    # if (method == "region") {
    #   for (f in 1:n_fleets) {
    #     region_idx <- which(unique_regions == fleet_regions[f])  # Get region index
    #     gear_idx <- which(unique_gears == fleet_pointer[f])  # Get gear index
    #     idx <- which(fleet_regions == region_idx & fleet_pointer == gear_idx)
    #     weight_matrix[, idx] <- region_weights[gear_idx]/sum(fleet_pointer == gear_idx)
    #   }
    #   if (is.matrix(advice)) {
    #     final_advice = weight_matrix*rowSums(advice)
    #   } else {
    #     final_advice = weight_matrix*sum(advice)
    #   }
    # }
    
    # if (method == "gear") {
    #   for (f in 1:n_fleets) {
    #     region_idx <- which(unique_regions == fleet_regions[f])  # Get region index
    #     gear_idx <- which(unique_gears == fleet_pointer[f])  # Get gear index
    #     idx <- which(fleet_regions == region_idx & fleet_pointer == gear_idx)
    #     weight_matrix[, idx] <- gear_weights[gear_idx]/sum(fleet_pointer == gear_idx)
    #   }
    #   if (is.matrix(advice)) {
    #     final_advice = weight_matrix*rowSums(advice)
    #   } else {
    #     final_advice = weight_matrix*sum(advice)
    #   }
    # }
    
    # if (method == "combined") {
    #   region_advice = matrix(0, nrow = nrow(advice), ncol = n_regions)
    #   for (f in 1:n_fleets) {
    #     region_idx <- which(unique_regions == fleet_regions[f])  # Get region index
    #     if (is.matrix(advice)) {
    #       region_advice[,region_idx] = rowSums(advice) * region_weights[region_idx]
    #     } else {
    #       region_advice[,region_idx] = sum(advice) * region_weights[region_idx]
    #     }
    #     gear_idx <- which(unique_gears == fleet_pointer[f])  # Get gear index
    #     idx <- which(fleet_regions == region_idx & fleet_pointer == gear_idx)
    #     final_advice[, idx] <- region_advice[,region_idx]*gear_weights[gear_idx]
    #   }
    # }
  }
  
  if (weight_type == 3) {
    
    unique_regions <- unique(fleet_regions)   # Extract unique regions
    unique_gears <- unique(fleet_pointer)     # Extract unique gear types
    
    avg_years <- (final_year - weight_years + 1):final_year
    
    if (length(which(index_pointer == index_id)) != n_regions) stop("Survey is not conducted in every region!")
    catch <- colMeans(om$input$data$agg_indices[which(om$years %in% avg_years),  which(index_pointer == index_id), drop = FALSE])
    
    index_weights <- catch / sum(catch)
    
    count = NULL
    for (i in 1:n_regions) count[i] = sum(fleet_regions == i)
    
    weight_matrix <- matrix(0, n_years, n_fleets)
    for (f in 1:n_fleets) {
      region_idx <- which(unique_regions == fleet_regions[f])  # Get region index
      gear_idx <- which(unique_gears == fleet_pointer[f])  # Get gear index
      idx <- which(fleet_regions == region_idx & fleet_pointer == gear_idx)
      weight_matrix[, idx] <- index_weights[region_idx]/count[region_idx]
    }
    
    if (is.matrix(advice)) {
      final_advice = weight_matrix*rowSums(advice)
    } else {
      final_advice = weight_matrix*sum(advice)
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
