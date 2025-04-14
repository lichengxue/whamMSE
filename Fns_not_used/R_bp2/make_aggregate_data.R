#' Aggregate Data for Estimation Model
#'
#' Aggregates catch and index data from the operating model (`om`) and updates
#' the estimation model (`em_info`) with the aggregated data. It also adjusts
#' `par_inputs` and weight-at-age (`user_waa`) to reflect the modified fleet,
#' region, index, and stock structure.
#'
#' @param om List. The operating model containing observed data, including:
#' @param em_info List. The estimation model information, including:
#' @param aggregate_catch_info List (optional). User-specified catch aggregation settings for panmictic models using aggregate catch.
#'   \itemize{
#'     \item `$n_fleets` Integer. Number of fleets.
#'     \item `$catch_cv` Numeric vector (`n_fleets`). CVs for annual aggregate catches by fleet.
#'     \item `$catch_Neff` Numeric vector (`n_fleets`). Effective sample sizes for fleet catches.
#'     \item `$use_agg_catch` Integer vector (`n_fleets`). 0/1 values flagging whether to use aggregate catches.
#'     \item `$use_catch_paa` Integer vector (`n_fleets`). 0/1 values flagging whether to use proportions at age observations.
#'     \item `$fleet_pointer` Integer vector (`n_fleets`). Defines fleet grouping (0 = exclude).
#'     \item `$use_catch_weighted_waa` Logical. Whether to use weighted weight-at-age based on fleet catches.
#'   }
#' @param aggregate_index_info List (optional). User-specified index aggregation settings for panmictic models using aggregate indices.
#'   \itemize{
#'     \item `$n_indices` Integer. Number of indices.
#'     \item `$index_cv` Numeric vector (`n_indices`). CVs for annual aggregate index catches.
#'     \item `$index_Neff` Numeric vector (`n_indices`). Effective sample sizes for survey indices.
#'     \item `$fracyr_indices` Numeric vector (`n_indices`). Fraction of the year for each survey index.
#'     \item `$q` Numeric vector (`n_indices`). Initial survey catchabilities.
#'     \item `$use_indices` Integer vector (`n_indices`). 0/1 values flagging whether to use survey indices.
#'     \item `$use_index_paa` Integer vector (`n_indices`). 0/1 values flagging whether to use proportions at age observations.
#'     \item `$units_indices` Integer vector (`n_indices`). 1/2 values flagging whether aggregate observations are biomass (1) or numbers (2).
#'     \item `$units_index_paa` Integer vector (`n_indices`). 1/2 values flagging whether composition observations are biomass (1) or numbers (2).
#'     \item `$index_pointer` Integer vector (`n_indices`). Defines index grouping (0 = exclude).
#'     \item `$use_index_weighted_waa` Logical. Whether to use weighted weight-at-age based on survey catches.
#'      }
#' @param aggregate_weights_info List (optional). Specifies how to compute weighted averages for
#' weight-at-age (`waa`) and maturity-at-age during data aggregation (For panmictic and fleets-as-areas models only).
#' Used for aggregating across fleets or indices to form total/stock/regional summaries.
#' \itemize{
#'   \item `$ssb_waa_weights` List. Settings for weighting weight-at-age used for spawning stock biomass.
#'     \itemize{
#'       \item `$fleet` Logical. Whether to use fleet-specific weights.
#'       \item `$index` Logical. Whether to use index-specific weights.
#'       \item `$pointer` Integer. Index pointing to the selected fleet or index group (e.g., 1 means the first valid fleet group).
#'     }
#'   \item `$maturity_weights` List. Settings for weighting maturity-at-age used for spawning stock biomass.
#'     \itemize{
#'       \item `$fleet` Logical. Whether to use fleet-specific weights.
#'       \item `$index` Logical. Whether to use index-specific weights.
#'       \item `$pointer` Integer. Index pointing to the selected fleet or index group.
#'     }
#' }
#' @param ind_em Vector. Indices specifying the years for which the estimation model should use data.
#'
#' @return List. An updated `em_info` object with:
#'   \itemize{
#'     \item Aggregated `catch_info` (`agg_catch`, `catch_paa`).
#'     \item Aggregated `index_info` (`agg_indices`, `index_paa`).
#'     \item Updated `par_inputs` with fleet/index aggregations.
#'     \item Updated weight-at-age (`user_waa`) and `basic_info$waa`.
#'   }
#'
#' @details
#' - Fleets and indices with `fleet_pointer = 0` or `index_pointer = 0` are ignored during aggregation.
#' - **Aggregated parameters** (`catch_cv`, `index_cv`, `q`, etc.) will return **vectors of length equal to the number of aggregated fleets or indices**.
#' - If `aggregate_catch_info` or `aggregate_index_info` are missing values, they are replaced with the first element of `em_info$par_inputs`.
#' - Weight-at-age (`waa`) is aggregated based on fleet or index weights, using either **catch-weighted averages or equal proportions**.
#'
#' @export
#' 
make_aggregate_data <- function(om, em_info, ind_em, aggregate_catch_info, aggregate_index_info, aggregate_weights_info) {
  
  data <- om$input$data
  
  aggregate_parameters <- function(values, pointers, default_values) {
    # Ignore 0s in pointers
    valid_pointers <- pointers[pointers > 0]
    unique_groups <- unique(valid_pointers)
    
    # If `values` is NULL, replace with the first element of `default_values`
    
    if (is.null(values)) {
      values <- rep(default_values[1], length(pointers))
    }
    
    # If values already match unique pointers, return directly
    
    if (length(values) == length(unique_groups)) {
      return(values)
    }
    
    # Ensure `values` and `pointers` have the same length before proceeding
    
    if (length(values) != length(pointers)) {
      stop("Length of `values` and `pointers` must be the same")
    }
    
    # Otherwise, compute aggregated values
    
    aggregated_values <- numeric(length(unique_groups))
    
    for (i in seq_along(unique_groups)) {
      group_indices <- which(pointers == unique_groups[i])  # Find all matching pointers
      aggregated_values[i] <- mean(values[group_indices], na.rm = TRUE)  # Compute mean for each group
    }
    
    return(aggregated_values)
  }
  
  # ---- Define fleet_pointer & index_pointer BEFORE Aggregation ---- #
  
  fleet_pointer <- aggregate_catch_info$fleet_pointer
  valid_fleets <- unique(fleet_pointer[fleet_pointer > 0])
  
  index_pointer <- aggregate_index_info$index_pointer
  valid_indices <- unique(index_pointer[index_pointer > 0])
  
  # Aggregate catch-related parameters
  
  em_info$par_inputs$catch_cv <- aggregate_parameters(aggregate_catch_info$catch_cv, fleet_pointer, em_info$par_inputs$catch_cv)
  em_info$par_inputs$catch_Neff <- aggregate_parameters(aggregate_catch_info$catch_Neff, fleet_pointer, em_info$par_inputs$catch_Neff)
  em_info$par_inputs$use_agg_catch <- aggregate_parameters(aggregate_catch_info$use_agg_catch, fleet_pointer, em_info$par_inputs$use_agg_catch)
  em_info$par_inputs$use_catch_paa <- aggregate_parameters(aggregate_catch_info$use_catch_paa, fleet_pointer, em_info$par_inputs$use_catch_paa)
  
  # Aggregate index-related parameters
  em_info$par_inputs$index_cv <- aggregate_parameters(aggregate_index_info$index_cv, index_pointer, em_info$par_inputs$index_cv)
  em_info$par_inputs$index_Neff <- aggregate_parameters(aggregate_index_info$index_Neff, index_pointer, em_info$par_inputs$index_Neff)
  em_info$par_inputs$q <- aggregate_parameters(aggregate_index_info$q, index_pointer, em_info$par_inputs$q)
  em_info$par_inputs$use_indices <- aggregate_parameters(aggregate_index_info$use_indices, index_pointer, em_info$par_inputs$use_indices)
  em_info$par_inputs$use_index_paa <- aggregate_parameters(aggregate_index_info$use_index_paa, index_pointer, em_info$par_inputs$use_index_paa)
  em_info$par_inputs$units_indices <- aggregate_parameters(aggregate_index_info$units_indices, index_pointer, em_info$par_inputs$units_indices)
  em_info$par_inputs$fracyr_indices <- aggregate_parameters(aggregate_index_info$fracyr_indices, index_pointer, em_info$par_inputs$fracyr_indices)
  em_info$par_inputs$units_index_paa <- aggregate_parameters(aggregate_index_info$units_index_paa, index_pointer, em_info$par_inputs$units_index_paa)
  
  # ---- Set Defaults for Missing Parameters (Using `em_info$par_inputs`) ---- #
  
  if (is.null(aggregate_catch_info$n_fleets)) aggregate_catch_info$n_fleets <- em_info$par_inputs$n_fleets
  if (is.null(aggregate_index_info$n_indices)) aggregate_index_info$n_indices <- em_info$par_inputs$n_indices
  
  n_fleets <- aggregate_catch_info$n_fleets
  n_indices <- aggregate_index_info$n_indices
  n_ages <- em_info$par_inputs$n_ages
  
  # ---- Aggregate Catch Data ---- #
  
  fleet_pointer <- aggregate_catch_info$fleet_pointer
  valid_fleets <- unique(fleet_pointer[fleet_pointer > 0])
  
  agg_catch <- matrix(NA, length(ind_em), length(valid_fleets))
  agg_catch_paa <- array(0, dim = c(length(valid_fleets), length(ind_em), n_ages))  # Initialize with zeros
  use_catch_paa <- em_info$par_inputs$use_catch_paa[valid_fleets]  
  
  for (f in seq_along(valid_fleets)) {
    fleet_group <- which(fleet_pointer == valid_fleets[f])  # Get fleets belonging to group f
    fleet_catch <- data$agg_catch[ind_em, fleet_group, drop = FALSE]  # Extract fleet catch
    total_catch <- rowSums(fleet_catch, na.rm = TRUE)  # Sum total catch for fleet group
    agg_catch[, f] <- total_catch  # Store aggregated catch
    
    if (use_catch_paa[f] == 1) {
      catch_paa <- data$catch_paa[fleet_group, ind_em, , drop = FALSE]  # Extract PAA for relevant fleets
      
      # Initialize an empty matrix for weighted PAA aggregation
      weighted_paa <- matrix(0, nrow = length(ind_em), ncol = n_ages)
      
      # Compute weighted PAA
      for (i in seq_along(fleet_group)) {
        if (sum(fleet_catch[, i], na.rm = TRUE) > 0) {  # Avoid division by zero
          weighted_paa <- weighted_paa + (catch_paa[i, , ] * fleet_catch[, i])
        }
      }
      
      # Normalize to ensure proportions sum to 1
      for (j in 1:length(ind_em)) {
        total_paa <- sum(weighted_paa[j, ], na.rm = TRUE)
        if (total_paa > 0) {
          agg_catch_paa[f, j, ] <- weighted_paa[j, ] / total_paa  # Normalize proportions
        } else {
          agg_catch_paa[f, j, ] <- rep(0, n_ages)  # Set all to zero if no valid data
        }
      }
    }
  }
  
  # ---- Aggregate Index Data ---- #
  
  index_pointer <- aggregate_index_info$index_pointer
  valid_indices <- unique(index_pointer[index_pointer > 0])
  
  agg_indices <- matrix(NA, length(ind_em), length(valid_indices))
  agg_index_paa <- array(0, dim = c(length(valid_indices), length(ind_em), n_ages))  # Initialize with zeros
  use_index_paa <- em_info$par_inputs$use_index_paa[valid_indices]  
  
  for (f in seq_along(valid_indices)) {
    index_group <- which(index_pointer == valid_indices[f])  # Get indices belonging to group f
    index_catch <- data$agg_indices[ind_em, index_group, drop = FALSE]  # Extract index catch
    total_index_catch <- rowSums(index_catch, na.rm = TRUE)  # Sum total index catch for index group
    agg_indices[, f] <- total_index_catch  # Store aggregated index catch
    
    if (use_index_paa[f] == 1) {
      index_paa <- data$index_paa[index_group, ind_em, , drop = FALSE]  # Extract PAA for relevant indices
      
      # Initialize an empty matrix for weighted index PAA aggregation
      weighted_index_paa <- matrix(0, nrow = length(ind_em), ncol = n_ages)
      
      # Compute weighted PAA
      for (i in seq_along(index_group)) {
        if (sum(index_catch[, i], na.rm = TRUE) > 0) {  # Avoid division by zero
          weighted_index_paa <- weighted_index_paa + (index_paa[i, , ] * index_catch[, i])
        }
      }
      
      # Normalize to ensure proportions sum to 1
      for (j in 1:length(ind_em)) {
        total_paa <- sum(weighted_index_paa[j, ], na.rm = TRUE)
        if (total_paa > 0) {
          agg_index_paa[f, j, ] <- weighted_index_paa[j, ] / total_paa  # Normalize proportions
        } else {
          agg_index_paa[f, j, ] <- rep(0, n_ages)  # Set all to zero if no valid data
        }
      }
    }
  }
  
  em_info$par_inputs$n_fleets <- n_fleets
  em_info$par_inputs$fleet_regions <- rep(1, n_fleets)
  
  em_info$par_inputs$n_indices <- n_indices
  em_info$par_inputs$index_regions <- rep(1, n_indices)
  
  # ---- Exclude fleets and indices marked with `0` in pointers ---- #
  if (is.null(aggregate_catch_info$fleet_pointer)) warnings("aggregate_catch_info$fleet_pointer is not specified!")
  if (is.null(aggregate_catch_info$fleet_pointer)) aggregate_catch_info$fleet_pointer = rep(1,aggregate_catch_info$n_fleets)
  if (length(aggregate_catch_info$fleet_pointer) != om$input$data$n_fleets) stop("Length of aggregate_catch_info$fleet_pointer is wrong!")
  fleet_pointer <- aggregate_catch_info$fleet_pointer
  valid_fleets <- unique(fleet_pointer[fleet_pointer > 0])
  
  if (is.null(aggregate_index_info$index_pointer)) stop("aggregate_index_info$index_pointer is not specified!")
  if (length(aggregate_index_info$index_pointer) != om$input$data$n_indices) stop("Length of aggregate_index_info$index_pointer is wrong!")
  index_pointer <- aggregate_index_info$index_pointer
  valid_indices <- unique(index_pointer[index_pointer > 0])
  
  # ---- Compute Fleet WAA Weights ---- #
  
  fleet_weights <- list()
  
  for (f in valid_fleets) {
    
    if (aggregate_catch_info$use_catch_weighted_waa) {
      
      agg_catch_filtered <- data$agg_catch[ind_em, which(fleet_pointer == f), drop = FALSE]
      
      if (ncol(agg_catch_filtered) > 1) {
        fleet_weights[[f]] <- agg_catch_filtered/rowSums(agg_catch_filtered)
      } else {
        fleet_weights[[f]] <- matrix(1, nrow = nrow(agg_catch_filtered), ncol = 1)
      }
    } else {
      
      if (ncol(agg_catch_filtered) > 1) {
        fleet_weights[[f]][] <- 1/ncol(agg_catch_filtered)
      } else {
        fleet_weights[[f]][] <- 1
      }
    }
  }
  
  # ---- Compute Index WAA Weights ---- #
  
  index_weights <- list()
  
  for (f in valid_indices) {
    
    if (aggregate_index_info$use_catch_weighted_waa) {
      
      agg_indices_filtered <- data$agg_indices[ind_em, which(index_pointer == f), drop = FALSE]
      
      if (ncol(agg_indices_filtered) > 1) {
        index_weights[[f]] <- agg_indices_filtered/rowSums(agg_indices_filtered)
      } else {
        index_weights[[f]] <- matrix(1, nrow = nrow(agg_indices_filtered), ncol = 1)
      }
    } else {
      
      if (ncol(agg_indices_filtered) > 1) {
        index_weights[[f]][] <- 1/ncol(agg_indices_filtered)
      } else {
        index_weights[[f]][] <- 1
      }
    }
  }
  
  
  # ---- Aggregate Fleet WAA (Weighted Average) ---- #
  
  aggregated_fleet_waa = list()
  
  for (f in unique(valid_fleets)) {
    
    waa_pointer_fleets <- em_info$par_inputs$user_waa$waa_pointer_fleets
    waa_mat <- em_info$par_inputs$user_waa$waa[waa_pointer_fleets,ind_em,]
    waa_mat <- waa_mat[fleet_pointer == f,,]
    
    weights <- fleet_weights[[f]]
    
    if (length(dim(waa_mat)) == 2) {
      
      waa_mat <- array(waa_mat, dim = c(1, dim(waa_mat)))
      
    }
    
    weighted_waa <- waa_mat * array(weights, dim = dim(waa_mat))
    
    aggregated <- apply(weighted_waa, c(2,3), sum)
    
    aggregated_fleet_waa[[f]] <- array(aggregated, dim = c(1, dim(aggregated)))
    
  }
  
  if (f == 1) {
    aggregated_fleet_waa = aggregated_fleet_waa[[1]]
  } else {
    aggregated_fleet_waa = abind::abind(aggregated_fleet_waa[1:f], along = 1)
  }
  
  # ---- Aggregate Index WAA (Weighted Average) ---- #
  
  aggregated_index_waa = list()
  
  for (f in unique(valid_indices)) {
    
    waa_pointer_indices <- em_info$par_inputs$user_waa$waa_pointer_indices
    waa_mat <- em_info$par_inputs$user_waa$waa[waa_pointer_indices,ind_em,]
    waa_mat <- waa_mat[index_pointer == f,,]

    weights <- index_weights[[f]]
    
    if (length(dim(waa_mat)) == 2) {
      
      waa_mat <- array(waa_mat, dim = c(1, dim(waa_mat)))
      
    }
    
    weighted_waa <- waa_mat * array(weights, dim = dim(waa_mat))
    
    aggregated <- apply(weighted_waa, c(2,3), sum)
    
    aggregated_index_waa[[f]] <- array(aggregated, dim = c(1, dim(aggregated)))
    
  }
  
  if (f == 1) {
    aggregated_index_waa = aggregated_index_waa[[1]]
  } else {
    aggregated_index_waa = abind::abind(aggregated_index_waa[1:f], along = 1)
  }
    
  # ---- Aggregate Regional WAA (Weighted Average) ---- #
  
  waa_pointer_totcatch <- em_info$par_inputs$user_waa$waa_pointer_totcatch
  waa_mat <- em_info$par_inputs$user_waa$waa[waa_pointer_totcatch,ind_em,]
  
  if (length(dim(waa_mat)) == 2) {
    
    waa_mat <- array(waa_mat, dim = c(1, dim(waa_mat)))
    
  }
  
  ssb_waa_weights = aggregate_weights_info$ssb_waa_weights
  
  if(is.null(ssb_waa_weights)) {
    weights <- matrix(1/dim(waa_mat)[1], nrow = dim(waa_mat)[1], ncol = dim(waa_mat)[2])
  } else {
    if(ssb_waa_weights$fleet) {
      pointer = ssb_waa_weights$pointer
      weights = t(fleet_weights[[pointer]])
    }
    if(ssb_waa_weights$index) {
      pointer = ssb_waa_weights$pointer
      weights = t(index_weights[[pointer]])
    }
  }
  
  weighted_waa <- waa_mat * array(weights, dim = dim(waa_mat))
  
  aggregated <- apply(weighted_waa, c(2,3), sum)
  
  aggregated_region_waa <- array(aggregated, dim = c(1, dim(aggregated)))
  
  # ---- Aggregate Stock WAA (Weighted Average) ---- #
  
  waa_pointer_ssb <- em_info$par_inputs$user_waa$waa_pointer_ssb
  waa_mat <- em_info$par_inputs$user_waa$waa[waa_pointer_ssb,ind_em,]
  
  if (length(dim(waa_mat)) == 2) {
    
    waa_mat <- array(waa_mat, dim = c(1, dim(waa_mat)))
    
  }
  
  ssb_waa_weights = aggregate_weights_info$ssb_waa_weights
  
  if(is.null(ssb_waa_weights)) {
    weights <- matrix(1/dim(waa_mat)[1], nrow = dim(waa_mat)[1], ncol = dim(waa_mat)[2])
  } else {
    if(ssb_waa_weights$fleet) {
      pointer = ssb_waa_weights$pointer
      weights = t(fleet_weights[[pointer]])
    }
    if(ssb_waa_weights$index) {
      pointer = ssb_waa_weights$pointer
      weights = t(index_weights[[pointer]])
    }
  }

  weighted_waa <- waa_mat * array(weights, dim = dim(waa_mat))
  
  aggregated <- apply(weighted_waa, c(2,3), sum)
  
  aggregated_stock_waa <- array(aggregated, dim = c(1, dim(aggregated)))
  
  # Store the aggregated WAA back into `em_info`
  waa_list <- list(aggregated_fleet_waa, aggregated_region_waa,
                   aggregated_index_waa, aggregated_stock_waa)
  
  em_info$par_inputs$user_waa$waa <- abind::abind(waa_list, along = 1)
  em_info$basic_info$waa = em_info$par_inputs$user_waa$waa
  
  # Average MAA to get global MAA
  
  maturity <- em_info$par_inputs$user_maturity[,ind_em,]
  
  if (length(dim(maturity)) == 2) {
    
    maturity <- array(maturity, dim = c(1, dim(maturity)))
    
  }
  
  maturity_weights = aggregate_weights_info$maturity_weights
    
  if(is.null(maturity_weights)) {
    weights <- matrix(1/dim(maturity)[1], nrow = dim(maturity)[1], ncol = dim(maturity)[2])
  } else {
    if(maturity_weights$fleet) {
      pointer = maturity_weights$pointer
      weights = t(fleet_weights[[pointer]])
    }
    if(maturity_weights$index) {
      pointer = maturity_weights$pointer
      weights = t(index_weights[[pointer]])
    }
  }
  
  weighted_maturity <- maturity * array(weights, dim = dim(maturity))
  
  aggregated <- apply(weighted_maturity, c(2,3), sum)
  
  aggregated_maturity <- array(aggregated, dim = c(1, dim(aggregated)))
  
  em_info$par_inputs$user_maturity <- aggregated_maturity
  em_info$basic_info$maturity <- aggregated_maturity
  
  em_info$par_inputs$n_regions = n_regions = 1
  em_info$par_inputs$n_stocks = n_stocks = 1 
  em_info$par_inputs$n_fleets = n_fleets
  em_info$par_inputs$n_indices = n_indices
  
  em_info$par_inputs$agg_catch <- agg_catch
  em_info$par_inputs$catch_paa <- agg_catch_paa  
  
  em_info$par_inputs$agg_indices <- agg_indices
  em_info$par_inputs$index_paa <- agg_index_paa  
  
  em_info$basic_info$waa_pointer_fleets   <- 1:n_fleets
  em_info$basic_info$waa_pointer_totcatch <- n_fleets + n_regions
  em_info$basic_info$waa_pointer_indices  <- (n_fleets + n_regions + 1):(n_fleets + n_regions + n_indices)
  em_info$basic_info$waa_pointer_ssb      <- (n_fleets + n_regions + n_indices + 1):(n_fleets + n_regions + n_indices + 1)
  em_info$basic_info$waa_pointer_M        <- em_info$basic_info$waa_pointer_ssb
  
  em_info$par_inputs$user_waa$waa_pointer_fleets   <- 1:n_fleets
  em_info$par_inputs$user_waa$waa_pointer_totcatch <- n_fleets + n_regions
  em_info$par_inputs$user_waa$waa_pointer_indices  <- (n_fleets + n_regions + 1):(n_fleets + n_regions + n_indices)
  em_info$par_inputs$user_waa$waa_pointer_ssb      <- (n_fleets + n_regions + n_indices + 1):(n_fleets + n_regions + n_indices + 1)
  em_info$par_inputs$user_waa$waa_pointer_M        <- em_info$par_inputs$user_waa$waa_pointer_ssb
  
  # return(list(em_info = em_info, agg_catch = agg_catch, catch_paa = catch_paa, agg_indices = agg_indices, index_paa = index_paa))
  return(em_info)
}