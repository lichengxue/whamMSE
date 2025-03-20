#' Filter Fleets and Indices for an Estimation Model
#'
#' This function filters fleets and indices based on region specifications
#' and removes indices marked for exclusion while ensuring correct renumbering.
#' Additionally, it extracts and correctly assembles the weight-at-age (WAA) matrix
#' for each region based on fleet, region, index, and stock assignments.
#'
#' @param em_info List. The original estimation model input data.
#' @param fleet_regions Integer vector. Specifies the region assignment of each fleet.
#' @param index_regions Integer vector. Specifies the region assignment of each index.
#' @param filter_indices Integer vector (optional). Indicates which indices to keep (1) or exclude (0).
#'
#' @return A list containing modified `em_info` lists for each region with filtered fleets, indices, renumbered index regions, and extracted weight-at-age (WAA) matrices.
#'
#' @export
#'
#'
filter_and_generate_em_info <- function(em_info, fleet_regions, index_regions, filter_indices = NULL, em.opt) {
  
  # Helper function to subset values if they are a vector of expected length, otherwise use first value
  subset_or_use_first <- function(param, idx, expected_length) {
    if (length(param) == expected_length) {
      return(param[idx])  # Subset correctly
    } else {
      return(rep(param[1], length(idx)))  # Use first value if length is mismatched
    }
  }
  
  # Helper function to filter and renumber index regions
  filter_and_renumber_indices <- function(index_regions, filter_indices) {
    if (!is.null(filter_indices)) {
      if (length(index_regions) != length(filter_indices))
        stop("index_regions and filter_indices must have the same length!")
      
      # Apply filtering (remove indices marked as 0)
      filtered_indices <- index_regions[filter_indices != 0]
    } else {
      # If filter_indices is NULL, use all indices as they are
      filtered_indices <- index_regions
    }
    
    # Renumber indices sequentially if any were removed
    if (length(filtered_indices) > 0) {
      unique_values <- sort(unique(filtered_indices))
      new_values <- seq_along(unique_values)
      mapping <- setNames(new_values, unique_values)
      renumbered_indices <- as.integer(mapping[as.character(filtered_indices)])
    } else {
      renumbered_indices <- integer(0)  # Return empty if all indices are removed
    }
    
    return(renumbered_indices)
  }
  
  if (em.opt$separate.em == FALSE) em.opt$separate.em.type = 0
  
  if (em.opt$separate.em.type == 2) {
    
    n_regions <- em_info$basic_info$n_regions
    
    # Create a new em_info structure for each region
    em_info_new <- em_info
    
    # Apply `filter_indices` (remove excluded indices) only if provided
    if (!is.null(filter_indices) && any(filter_indices != 0)) {
      relevant_indices <- which(filter_indices != 0)
      em_info_new$basic_info$n_indices = length(relevant_indices)
    } else {
      relevant_indices <- 1:em_info_new$par_inputs$n_indices
    }
    
    n_indices <- length(relevant_indices)
    
    if (n_indices > 0) {
      # Subset index-related parameters
      em_info_new$par_inputs$index_cv <- subset_or_use_first(em_info$par_inputs$index_cv, relevant_indices, length(index_regions))
      em_info_new$par_inputs$index_Neff <- subset_or_use_first(em_info$par_inputs$index_Neff, relevant_indices, length(index_regions))
      em_info_new$par_inputs$fracyr_indices <- subset_or_use_first(em_info$par_inputs$fracyr_indices, relevant_indices, length(index_regions))
      em_info_new$par_inputs$q <- subset_or_use_first(em_info$par_inputs$q, relevant_indices, length(index_regions))
      em_info_new$par_inputs$use_indices <- subset_or_use_first(em_info$par_inputs$use_indices, relevant_indices, length(index_regions))
      em_info_new$par_inputs$use_index_paa <- subset_or_use_first(em_info$par_inputs$use_index_paa, relevant_indices, length(index_regions))
      em_info_new$par_inputs$units_indices <- subset_or_use_first(em_info$par_inputs$units_indices, relevant_indices, length(index_regions))
      em_info_new$par_inputs$units_index_paa <- subset_or_use_first(em_info$par_inputs$units_index_paa, relevant_indices, length(index_regions))
      
      # Renumber index_regions for remaining indices
      em_info_new$par_inputs$index_regions <- filter_and_renumber_indices(index_regions[relevant_indices], filter_indices[relevant_indices])
    } else {
      em_info_new$par_inputs$index_regions <- NULL  # No indices left
    }
    
    aggregated_region_waa = list()
    for (f in 1) { # Assuming n_regions = 1
      ind = length(fleet_regions) + 1:em_info$basic_info$n_regions
      if(is.null(em_info$par_inputs$user_waa)) {
        aggregated_region_waa[[f]] <- em_info$basic_info$waa[ind,1, ]
        if (is.matrix(aggregated_region_waa[[f]])) aggregated_region_waa[[f]] <- colMeans(as.matrix(aggregated_region_waa[[f]]))
      } else {
        aggregated_region_waa[[f]] <- em_info$par_inputs$user_waa[, ]
        if (is.matrix(aggregated_region_waa[[f]])) aggregated_region_waa[[f]] <- colMeans(as.matrix(aggregated_region_waa[[f]]))
      }
    }
    
    aggregated_stock_waa = aggregated_region_waa
    
    # Extract WAA matrix for the region
    n_fleets = length(fleet_regions)
    waa_indices <- c(1:n_fleets,
                     n_fleets + 1,
                     n_fleets + n_regions + relevant_indices,
                     n_fleets + n_regions + length(index_regions) + 1)
    
    em_info_new$par_inputs$user_waa = em_info$basic_info$waa[waa_indices, 1, ] # only first year will be used for simplicity!
    em_info_new$par_inputs$user_waa[n_fleets+1, ] = aggregated_region_waa[[1]]
    em_info_new$par_inputs$user_waa[dim(em_info_new$par_inputs$user_waa)[1], ] = aggregated_region_waa[[1]]
    
    em_info_new$basic_info$waa = em_info$basic_info$waa[waa_indices, , ]
    em_info_new$basic_info$waa_pointer_fleets   <- 1:n_fleets
    em_info_new$basic_info$waa_pointer_totcatch <- n_fleets + 1 # n_region must be 1 for FAA!
    em_info_new$basic_info$waa_pointer_indices  <- (n_fleets + 1 + 1):(n_fleets + 1 + n_indices)
    em_info_new$basic_info$waa_pointer_ssb      <- (n_fleets + 1 + n_indices + 1):(n_fleets + 1 + n_indices + 1)
    em_info_new$basic_info$waa_pointer_M        <- em_info_new$basic_info$waa_pointer_ssb

    em_info_new$par_inputs$n_stocks = em_info_new$par_inputs$n_regions = 1
    em_info_new$par_inputs$n_indices = length(relevant_indices)
    em_info_new$par_inputs$index_regions = rep(1,length(relevant_indices))
    
    # Store in list
    em_info_list = em_info_new
  }
  
  # Important: I don't think this model can handle missing survey (using filter_indices) #
  if (em.opt$separate.em.type == 3) {
    # Initialize a list to store per-region em_info
    em_info_list <- list()
    n_regions <- length(unique(fleet_regions))
    
    for (r in 1:n_regions) {
      
      # Create a new em_info structure for each region
      em_info_new <- em_info
      
      # Extract relevant fleets for region `r`
      relevant_fleets <- which(fleet_regions == r)
      n_fleets <- length(relevant_fleets)
      
      # Subset fleet-related parameters
      em_info_new$par_inputs$catch_cv <- subset_or_use_first(em_info$par_inputs$catch_cv, relevant_fleets, length(fleet_regions))
      em_info_new$par_inputs$catch_Neff <- subset_or_use_first(em_info$par_inputs$catch_Neff, relevant_fleets, length(fleet_regions))
      em_info_new$par_inputs$use_agg_catch <- subset_or_use_first(em_info$par_inputs$use_agg_catch, relevant_fleets, length(fleet_regions))
      em_info_new$par_inputs$use_catch_paa <- subset_or_use_first(em_info$par_inputs$use_catch_paa, relevant_fleets, length(fleet_regions))
      
      # Extract relevant indices for region `r`
      relevant_indices <- which(index_regions == r)
      
      # Apply `filter_indices` (remove excluded indices) only if provided
      if (is.null(filter_indices)) filter_indices = 1:relevant_indices
      
      if (!is.null(filter_indices)) {
        relevant_indices <- relevant_indices[filter_indices[relevant_indices] != 0]
      } 
      
      n_indices <- length(relevant_indices)
      
      if (n_indices > 0) {
        # Subset index-related parameters
        em_info_new$par_inputs$index_cv <- subset_or_use_first(em_info$par_inputs$index_cv, relevant_indices, length(index_regions))
        em_info_new$par_inputs$index_Neff <- subset_or_use_first(em_info$par_inputs$index_Neff, relevant_indices, length(index_regions))
        em_info_new$par_inputs$fracyr_indices <- subset_or_use_first(em_info$par_inputs$fracyr_indices, relevant_indices, length(index_regions))
        em_info_new$par_inputs$q <- subset_or_use_first(em_info$par_inputs$q, relevant_indices, length(index_regions))
        em_info_new$par_inputs$use_indices <- subset_or_use_first(em_info$par_inputs$use_indices, relevant_indices, length(index_regions))
        em_info_new$par_inputs$use_index_paa <- subset_or_use_first(em_info$par_inputs$use_index_paa, relevant_indices, length(index_regions))
        em_info_new$par_inputs$units_indices <- subset_or_use_first(em_info$par_inputs$units_indices, relevant_indices, length(index_regions))
        em_info_new$par_inputs$units_index_paa <- subset_or_use_first(em_info$par_inputs$units_index_paa, relevant_indices, length(index_regions))
        
        # Renumber index_regions for remaining indices
        em_info_new$par_inputs$index_regions <- filter_and_renumber_indices(index_regions[relevant_indices], filter_indices[relevant_indices])
      } else {
        em_info_new$par_inputs$index_regions <- NULL  # No indices left
      }
      
      # Extract WAA matrix for the region
      waa_indices <- c(relevant_fleets,
                       length(fleet_regions) + r,
                       length(fleet_regions) + n_regions + relevant_indices,
                       length(fleet_regions) + n_regions + length(index_regions) + r)
      
      em_info_new$par_inputs$user_waa = em_info$basic_info$waa[waa_indices, 1, ] # only first year will be used for simplicity!
      
      em_info_new$basic_info$waa = em_info$basic_info$waa[waa_indices, , ]
      em_info_new$basic_info$waa_pointer_fleets   <- 1:n_fleets
      em_info_new$basic_info$waa_pointer_totcatch <- n_fleets + 1
      em_info_new$basic_info$waa_pointer_indices  <- (n_fleets + 1 + 1):(n_fleets + 1 + n_indices)
      em_info_new$basic_info$waa_pointer_ssb      <- (n_fleets + 1 + n_indices + 1):(n_fleets + 1 + n_indices + 1)
      em_info_new$basic_info$waa_pointer_M        <- em_info_new$basic_info$waa_pointer_ssb
      
      em_info_new$par_inputs$n_indices = length(relevant_indices)
      em_info_new$par_inputs$index_regions = rep(1,length(relevant_indices))
      
      em_info_new$par_inputs$fleet_regions = rep(1,length(relevant_fleets))
      em_info_new$par_inputs$n_fleets = length(relevant_fleets)
      
      em_info_new$par_inputs$n_stocks = em_info_new$par_inputs$n_regions = 1
      
      # Store in list
      em_info_list[[r]] <- em_info_new
    }
  }
  
  if (em.opt$separate.em == FALSE) {
    
    n_regions <- em_info$basic_info$n_regions
    
    # Create a new em_info structure for each region
    em_info_new <- em_info
    
    # Apply `filter_indices` (remove excluded indices) only if provided
    if (!is.null(filter_indices) && any(filter_indices != 0)) {
      relevant_indices <- which(filter_indices != 0)
      em_info_new$basic_info$n_indices = length(relevant_indices)
    } else {
      relevant_indices <- 1:em_info_new$par_inputs$n_indices
    }
    
    n_indices <- length(relevant_indices)
    
    if (n_indices > 0) {
      # Subset index-related parameters
      em_info_new$par_inputs$n_indices <- n_indices
      em_info_new$par_inputs$index_cv <- subset_or_use_first(em_info$par_inputs$index_cv, relevant_indices, length(index_regions))
      em_info_new$par_inputs$index_Neff <- subset_or_use_first(em_info$par_inputs$index_Neff, relevant_indices, length(index_regions))
      em_info_new$par_inputs$fracyr_indices <- subset_or_use_first(em_info$par_inputs$fracyr_indices, relevant_indices, length(index_regions))
      em_info_new$par_inputs$q <- subset_or_use_first(em_info$par_inputs$q, relevant_indices, length(index_regions))
      em_info_new$par_inputs$use_indices <- subset_or_use_first(em_info$par_inputs$use_indices, relevant_indices, length(index_regions))
      em_info_new$par_inputs$use_index_paa <- subset_or_use_first(em_info$par_inputs$use_index_paa, relevant_indices, length(index_regions))
      em_info_new$par_inputs$units_indices <- subset_or_use_first(em_info$par_inputs$units_indices, relevant_indices, length(index_regions))
      em_info_new$par_inputs$units_index_paa <- subset_or_use_first(em_info$par_inputs$units_index_paa, relevant_indices, length(index_regions))
      
      # Renumber index_regions for remaining indices
      em_info_new$par_inputs$index_regions <- filter_and_renumber_indices(index_regions[relevant_indices], filter_indices[relevant_indices])
    } else {
      em_info_new$par_inputs$index_regions <- NULL  # No indices left
    }
    
    # Extract WAA matrix for the region
    n_fleets = length(fleet_regions)
    
    if (!is.null(filter_indices) && any(filter_indices != 0)) {
      max.dim = n_fleets + n_regions + length(index_regions) + n_regions
      index_keep = n_fleets + n_regions + relevant_indices
      index_dim = (n_fleets + n_regions + 1):(n_fleets + n_regions + length(index_regions))
      index_rm = setdiff(index_dim,index_keep)
      em_info_new$par_inputs$user_waa = em_info$basic_info$waa[-index_rm, 1, ]
      em_info_new$basic_info$waa = em_info$basic_info$waa[-index_rm, , ]
      em_info_new$basic_info$waa_pointer_fleets   <- 1:n_fleets
      em_info_new$basic_info$waa_pointer_totcatch <- n_fleets + 1:n_regions
      em_info_new$basic_info$waa_pointer_indices  <- (n_fleets + n_regions + 1):(n_fleets + n_regions + n_indices)
      em_info_new$basic_info$waa_pointer_ssb      <- (n_fleets + n_regions + n_indices + 1):(n_fleets + n_regions + n_indices + n_regions)
      em_info_new$basic_info$waa_pointer_M        <- em_info_new$basic_info$waa_pointer_ssb
      
    } 
    
    # Store in list
    em_info_list = em_info_new
  }
  
  return(em_info_list)
}