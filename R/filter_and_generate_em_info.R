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
#' @param reduce_region_info List (optional). Specifies modifications for regions. If `NULL`, no modifications are applied.
#'   The expected components include:
#'   \itemize{
#'     \item `$remove_regions`
#'       Specifies which regions should be removed from the model.
#'
#'     \item `$reassign`
#'       Specifies reassignment of surveys from removed regions to non-removed regions.
#'
#'     \item `$NAA_where`
#'       Specifies recruitment assignments after region reduction.
#'
#'     \item `$sel_em`, `$M_em`, `$NAA_re_em`, `$move_em` 
#'       Model settings that may change based on region reduction.
#'
#'     \item `$onto_move_list`
#'       Contains movement-related parameters with the following elements:
#'       \itemize{
#'         \item `$onto_move` (array, dimension: `n_stocks × n_regions × (n_regions-1)`)  
#'           Specifies movement rules between stocks and regions. Default = NULL.
#'
#'         \item `$onto_move_pars` (array, dimension: `n_stocks × n_regions × (n_regions-1) × 4`)  
#'           Specifies movement parameters. Default = NULL.
#'
#'         \item `$age_mu_devs` (array, dimension: `n_stocks × n_regions × (n_regions-1) × n_ages`)  
#'           Stores age-based movement deviations. If `onto_move == 5`, values are extracted from `basic_info`.
#'       }
#'   }
#'   
#' @return A list containing modified `em_info` lists for each region with filtered fleets, indices, renumbered index regions, and extracted weight-at-age (WAA) matrices.
#' 
#' @export
#'
#'
filter_and_generate_em_info <- function(em_info, fleet_regions, index_regions, filter_indices = NULL, 
                                        reduce_region_info = NULL, em.opt) {
  
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
    
    if (!is.null(filter_indices) && any(filter_indices != 0) && !all(filter_indices == 1)) {
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
    
    if (!is.null(reduce_region_info$remove_regions)) {
      
      remove_regions <- reduce_region_info$remove_regions
      prev_n_fleets <- n_fleets  # Keep original fleet count
      prev_n_regions <- n_regions
      prev_n_indices <- n_indices
      
      # Identify fleets and indices associated with the regions to remove
      fleets_to_remove <- which(fleet_regions %in% which(remove_regions == 0))  # Get fleet indices
      indices_to_remove <- which(index_regions %in% which(remove_regions == 0))  # Get index indices
      
      # Remove fleets from fleet_regions if any
      if (length(fleets_to_remove) > 0) {
        fleet_regions <- fleet_regions[-fleets_to_remove]
      } else {
        fleet_regions <- fleet_regions
      }
      
      # Remove indices from index_regions if any
      if (length(indices_to_remove) > 0) {
        index_regions <- index_regions[-indices_to_remove]
      } else {
        index_regions = index_regions
      }
      
      # If just reassign index to one remaining region...
      if (!is.null(reduce_region_info$reassign)) {
        index_regions <- c(index_regions,reduce_region_info$reassign)
        indices_to_remove <- 0
      }
      
      waa_indices_to_remove <- c()
      
      if (length(fleets_to_remove) > 0) {
        waa_indices_to_remove <- c(waa_indices_to_remove, fleets_to_remove)  # Remove fleet-related WAA
      }
      if (length(remove_regions) > 0) {
        waa_indices_to_remove <- c(waa_indices_to_remove, prev_n_fleets + which(remove_regions == 0))  # Region-specific total catch
      }
      if (length(indices_to_remove) > 0) {
        waa_indices_to_remove <- c(waa_indices_to_remove, prev_n_fleets + prev_n_regions + indices_to_remove)  # Indices
      }
      if (length(remove_regions) > 0) {
        waa_indices_to_remove <- c(waa_indices_to_remove, 
                                   prev_n_fleets + prev_n_regions + prev_n_indices + which(remove_regions == 0))  # SSB and M
      }
      
      # Ensure valid range before subsetting
      max_dim <- dim(em_info_new$basic_info$waa)[1]
      waa_indices_to_remove <- waa_indices_to_remove[waa_indices_to_remove <= max_dim]
      
      # Remove WAA information for the removed regions (only if there's data to remove)
      if (length(waa_indices_to_remove) > 0) {
        em_info_new$par_inputs$user_waa <- em_info$basic_info$waa[-waa_indices_to_remove, 1, ]
        em_info_new$basic_info$waa <- em_info$basic_info$waa[-waa_indices_to_remove, , ]
      }
      
      # Update number of remaining regions
      n_regions <- length(unique(fleet_regions))  # Get unique remaining regions
      
      # Identify WAA indices to remove
      n_fleets <- length(fleet_regions)  # Updated fleet count
      n_indices <- length(index_regions) # Updated index count
      
      # Update pointers for remaining data
      em_info_new$basic_info$waa_pointer_fleets   <- 1:n_fleets
      em_info_new$basic_info$waa_pointer_totcatch <- n_fleets + 1:n_regions
      em_info_new$basic_info$waa_pointer_indices <- (n_fleets + n_regions + 1):(n_fleets + n_regions + n_indices)
      em_info_new$basic_info$waa_pointer_ssb <- (n_fleets + n_regions + n_indices + 1):(n_fleets + n_regions + n_indices + n_regions)
      em_info_new$basic_info$waa_pointer_M <- em_info_new$basic_info$waa_pointer_ssb
      
      # Also remove the region from the metadata
      em_info_new$par_inputs$n_regions <- n_regions
      em_info_new$par_inputs$n_stocks <- n_regions
      em_info_new$par_inputs$n_fleets <- n_fleets
      em_info_new$par_inputs$n_indices <- n_indices
      
      # Subset index-related parameters
      if (length(indices_to_remove) > 0) {
        if (is.null(reduce_region_info$reassign)) {
          em_info_new$par_inputs$index_cv <- em_info_new$par_inputs$index_cv[-indices_to_remove]
          em_info_new$par_inputs$index_Neff <- em_info_new$par_inputs$index_Neff[-indices_to_remove]
          em_info_new$par_inputs$fracyr_indices <- em_info_new$par_inputs$fracyr_indices[-indices_to_remove]
          em_info_new$par_inputs$q <- em_info_new$par_inputs$q[-indices_to_remove]
          em_info_new$par_inputs$use_indices <- em_info_new$par_inputs$use_indices[-indices_to_remove]
          em_info_new$par_inputs$use_index_paa <- em_info_new$par_inputs$use_index_paa[-indices_to_remove]
          em_info_new$par_inputs$units_indices <- em_info_new$par_inputs$units_indices[-indices_to_remove]
          em_info_new$par_inputs$units_index_paa <- em_info_new$par_inputs$units_index_paa[-indices_to_remove]
        } 
      }
      
      if(is.matrix(em_info_new$par_inputs$user_maturity)) em_info_new$par_inputs$user_maturity = em_info_new$par_inputs$user_maturity[which(remove_regions == 1),]
      
      if (length(fleets_to_remove) > 0) {
        em_info_new$par_inputs$catch_Neff <- em_info_new$par_inputs$catch_Neff[-fleets_to_remove]
        em_info_new$par_inputs$catch_cv <- em_info_new$par_inputs$catch_cv[-fleets_to_remove]
        em_info_new$par_inputs$use_agg_catch <- em_info_new$par_inputs$use_agg_catch[-fleets_to_remove]
        em_info_new$par_inputs$use_catch_paa <- em_info_new$par_inputs$use_catch_paa[-fleets_to_remove]
      }
      
      em_info_new$par_inputs$fleet_regions = fleet_regions
      em_info_new$par_inputs$index_regions = index_regions
      
      cat("If ontogenetic movement is assumed, make sure it is re-defined after regions are dropped from the model!")
      if(!is.null(reduce_region_info$onto_move_list)){
        em_info_new$par_inputs$onto_move <- reduce_region_info$onto_move_list$onto_move
        em_info_new$par_inputs$onto_move_pars <- reduce_region_info$onto_move_list$onto_move_pars
        em_info_new$par_inputs$age_mu_devs <- reduce_region_info$onto_move_list$age_mu_devs
      }
      if(all(em_info_new$par_inputs$onto_move == 0))em_info_new$par_inputs$onto_move = array(0, dim=c(n_stocks,n_regions,n_regions-1))
      if(all(em_info_new$par_inputs$onto_move_pars == 1))em_info_new$par_inputs$onto_move_pars = array(0, dim=c(n_stocks,n_regions,n_regions-1,4))
      if(all(em_info_new$par_inputs$age_mu_devs == 0))em_info_new$par_inputs$age_mu_devs = array(0, dim=c(n_stocks,n_regions,n_regions-1,dim(em_info_new$par_inputs$age_mu_devs)[4]))
      
      em_info_new$par_inputs$fleets_to_remove = fleets_to_remove
      em_info_new$par_inputs$indices_to_remove = indices_to_remove
    }
    
    # Store in list
    em_info_list = em_info_new
  }
  
  return(em_info_list)
}