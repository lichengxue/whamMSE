#' Aggregate Data for Estimation Model
#'
#' Aggregates catch and index data from the operating model (`om`) and updates
#' the estimation model (`em_info`) with the aggregated data. It also adjusts
#' `par_inputs` and weight-at-age (`user_waa`) to reflect the modified fleet,
#' region, index, and stock structure.
#'
#' @param om List. The operating model containing observed data, including:
#' @param em_info List. The estimation model information, including:
#' @param aggregate_catch_info List. Contains details for aggregating fleet-specific catch data:
#'   \itemize{
#'     \item `fleet_pointer` - Integer vector of length `n_fleets`, specifying how fleets should be aggregated.
#'     \item `catch_cv`, `catch_Neff` - Numeric vectors of length `n_fleets`, specifying catch CV and sample size.
#'     \item `use_agg_catch`, `use_catch_paa` - Boolean vectors of length `n_fleets`, indicating whether to aggregate.
#'   }
#' @param aggregate_index_info List. Contains details for aggregating index-specific data:
#'   \itemize{
#'     \item `index_pointer` - Integer vector of length `n_indices`, specifying how indices should be aggregated.
#'     \item `index_cv`, `index_Neff`, `q` - Numeric vectors of length `n_indices`, specifying index CV, effective sample size, and catchability coefficient.
#'     \item `use_indices`, `use_index_paa` - Boolean vectors of length `n_indices`, indicating whether to aggregate.
#'   }
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
#' @examples
#' \dontrun{
#' # Example usage
#' new_em_info <- make_aggregate_data(om, em_info, aggregate_catch_info, aggregate_index_info, ind_em)
#' }
#'
#' @export
#' 
make_aggregate_data <- function(om, em_info, aggregate_catch_info, aggregate_index_info, ind_em) {
  
  data <- om$input$data
  
  # Function to Aggregate Parameters (Handles Single Values & Group Averaging)
  # aggregate_parameters <- function(values, pointers, default_values) {
  #   # Ignore 0s in pointers
  #   valid_pointers <- pointers[pointers > 0]
  #   unique_groups <- unique(valid_pointers)
  #   
  #   # If `values` is NULL, replace with the first element of `default_values`
  #   if (is.null(values)) {
  #     values <- rep(default_values[1], length(pointers))
  #   }
  #   
  #   # Ensure `values` and `pointers` have the same length
  #   if (length(values) != length(pointers)) {
  #     stop("Length of `values` and `pointers` must be the same")
  #   }
  #   
  #   aggregated_values <- numeric(length(unique_groups))
  #   
  #   for (i in seq_along(unique_groups)) {
  #     group_indices <- which(pointers == unique_groups[i])  # Find all non-zero fleet groups
  #     aggregated_values[i] <- mean(values[group_indices], na.rm = TRUE)  # Compute mean for each group
  #   }
  #   
  #   return(aggregated_values)
  # }
  
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
  
  # for (f in seq_along(valid_fleets)) {
  #   fleet_group <- which(fleet_pointer == valid_fleets[f])
  #   fleet_catch <- data$agg_catch[ind_em, fleet_group, drop = FALSE]
  #   total_catch <- rowSums(fleet_catch)  
  #   agg_catch[, f] <- total_catch  
  #   
  #   if (use_catch_paa[f] == 1) {
  #     catch_paa <- data$catch_paa[fleet_group, ind_em, , drop = FALSE]  
  #     weighted_paa <- array(0, dim = c(length(ind_em), n_ages))  
  #     
  #     for (i in seq_along(fleet_group)) {
  #       weighted_paa <- weighted_paa + (catch_paa[i, , ] * fleet_catch[, i])  
  #     }
  #     
  #     agg_catch_paa[f, , ] <- t(apply(weighted_paa, 1, function(row) ifelse(sum(row) > 0, row / sum(row), 0)))  
  #   }
  # }
  
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
  
  # for (f in seq_along(valid_indices)) {
  #   index_group <- which(index_pointer == valid_indices[f])
  #   index_catch <- data$agg_indices[ind_em, index_group, drop = FALSE]
  #   total_index_catch <- rowSums(index_catch)  
  #   agg_indices[, f] <- total_index_catch  
  #   
  #   if (use_index_paa[f] == 1) {
  #     index_paa <- data$index_paa[index_group, ind_em, , drop = FALSE]  
  #     weighted_index_paa <- array(0, dim = c(length(ind_em), n_ages))  
  #     
  #     for (i in seq_along(index_group)) {
  #       weighted_index_paa <- weighted_index_paa + (index_paa[i, , ] * index_catch[, i])  
  #     }
  #     
  #     agg_index_paa[f, , ] <- t(apply(weighted_index_paa, 1, function(row) ifelse(sum(row) > 0, row / sum(row), 0)))  
  #   }
  # }
  
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
  
  # ---- Store in `em_info` ---- #
  # em_info$catch_info$agg_catch <- agg_catch
  # em_info$catch_info$catch_paa <- agg_catch_paa  
  # 
  # em_info$index_info$agg_indices <- agg_indices
  # em_info$index_info$index_paa <- agg_index_paa  
  
  # ---- Aggregate Parameters ---- #
  # em_info$par_inputs$catch_cv <- aggregate_catch_info$catch_cv[valid_fleets]
  # em_info$par_inputs$catch_Neff <- aggregate_catch_info$catch_Neff[valid_fleets]
  # em_info$par_inputs$use_agg_catch <- aggregate_catch_info$use_agg_catch[valid_fleets]
  # em_info$par_inputs$use_catch_paa <- use_catch_paa
  # 
  # em_info$par_inputs$index_cv <- aggregate_index_info$index_cv[valid_indices]
  # em_info$par_inputs$index_Neff <- aggregate_index_info$index_Neff[valid_indices]
  # em_info$par_inputs$use_indices <- aggregate_index_info$use_indices[valid_indices]
  # em_info$par_inputs$use_index_paa <- use_index_paa
  
  em_info$par_inputs$n_fleets <- n_fleets
  em_info$par_inputs$fleet_regions <- rep(1, n_fleets)
  
  em_info$par_inputs$n_indices <- n_indices
  em_info$par_inputs$index_regions <- rep(1, n_indices)
  
  # ---- Final Updates ---- #
  
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
  
  # ---- Compute WAA Weights ---- #
  fleet_weights <- list()
  for (f in valid_fleets) {
    if (aggregate_catch_info$use_catch_weighted_waa) {
      agg_catch_filtered <- data$agg_catch[ind_em, which(fleet_pointer == f), drop = FALSE]
      
      # Ensure weight vector has correct sum
      mean_catch <- if (is.matrix(agg_catch_filtered)) colMeans(agg_catch_filtered, na.rm = TRUE) else mean(agg_catch_filtered, na.rm = TRUE)
      fleet_weights[[f]] <- mean_catch / sum(mean_catch, na.rm = TRUE)  # Normalize weights so they sum to 1
    } else {
      fleet_weights[[f]] <- rep(1 / length(valid_fleets), length(valid_fleets))  # Equal weights
    }
  }
  
  index_weights <- list()
  for (f in valid_indices) {
    if (aggregate_index_info$use_catch_weighted_waa) {
      agg_indices_filtered <- data$agg_indices[ind_em, which(index_pointer == f), drop = FALSE]
      
      # Ensure weight vector has correct sum
      mean_index <- if (is.matrix(agg_indices_filtered)) colMeans(agg_indices_filtered, na.rm = TRUE) else mean(agg_indices_filtered, na.rm = TRUE)
      index_weights[[f]] <- mean_index / sum(mean_index, na.rm = TRUE)  # Normalize weights
    } else {
      index_weights[[f]] <- rep(1 / length(valid_indices), length(valid_indices))  # Equal weights
    }
  }
  
  # ---- Aggregate WAA ---- #
  aggregated_fleet_waa = list()
  for (f in unique(valid_fleets)) {
    if (is.null(em_info$par_inputs$user_waa)) {
      aggregated_fleet_waa[[f]] <- as.matrix(em_info$basic_info$waa[which(fleet_pointer == f), 1, ])
      if (ncol(aggregated_fleet_waa[[f]]) > 1) {
        aggregated_fleet_waa[[f]] <- colSums(aggregated_fleet_waa[[f]] * fleet_weights[[f]])
      } else {
        aggregated_fleet_waa[[f]] <- sum(aggregated_fleet_waa[[f]] * fleet_weights[[f]])
      }
    } else {
      aggregated_fleet_waa[[f]] <- as.matrix(em_info$par_inputs$user_waa[which(fleet_pointer == f), ])
      if (ncol(aggregated_fleet_waa[[f]]) > 1) {
        aggregated_fleet_waa[[f]] <- colSums(aggregated_fleet_waa[[f]] * fleet_weights[[f]])
      } else {
        aggregated_fleet_waa[[f]] <- sum(aggregated_fleet_waa[[f]] * fleet_weights[[f]])
      }
    }
  }
  
  aggregated_region_waa = list()
  for (f in 1) {  # Assuming n_regions = 1
    ind = length(fleet_pointer) + 1:em_info$basic_info$n_regions
    if (is.null(em_info$par_inputs$user_waa)) {
      aggregated_region_waa[[f]] <- as.matrix(em_info$basic_info$waa[ind, 1, ])
      if (ncol(aggregated_region_waa[[f]]) > 1) {
        aggregated_region_waa[[f]] <- colMeans(aggregated_region_waa[[f]])
      } else {
        aggregated_region_waa[[f]] <- mean(aggregated_region_waa[[f]])
      }
    } else {
      aggregated_region_waa[[f]] <- as.matrix(em_info$par_inputs$user_waa[ind, ])
      if (ncol(aggregated_region_waa[[f]]) > 1) {
        aggregated_region_waa[[f]] <- colMeans(aggregated_region_waa[[f]])
      } else {
        aggregated_region_waa[[f]] <- mean(aggregated_region_waa[[f]])
      }
    }
  }
  
  aggregated_index_waa = list()
  for (f in unique(valid_indices)) {
    ind = length(fleet_pointer) + 1:em_info$basic_info$n_regions + which(index_pointer == f)
    if (is.null(em_info$par_inputs$user_waa)) {
      aggregated_index_waa[[f]] <- as.matrix(em_info$basic_info$waa[ind, 1, ])
      if (ncol(aggregated_index_waa[[f]]) > 1) {
        aggregated_index_waa[[f]] <- colMeans(aggregated_index_waa[[f]])
      } else {
        aggregated_index_waa[[f]] <- mean(aggregated_index_waa[[f]])
      }
    } else {
      aggregated_index_waa[[f]] <- as.matrix(em_info$par_inputs$user_waa[ind, ])
      if (ncol(aggregated_index_waa[[f]]) > 1) {
        aggregated_index_waa[[f]] <- colSums(aggregated_index_waa[[f]] * index_weights[[f]])
      } else {
        aggregated_index_waa[[f]] <- sum(aggregated_index_waa[[f]] * index_weights[[f]])
      }
    }
  }
  
  aggregated_stock_waa = list()
  for (f in 1) {
    ind = length(fleet_pointer) + em_info$basic_info$n_regions + length(index_pointer) + 1:em_info$basic_info$n_stocks
    if (is.null(em_info$par_inputs$user_waa)) {
      aggregated_stock_waa[[f]] <- as.matrix(em_info$basic_info$waa[ind, 1, ])
      if (ncol(aggregated_stock_waa[[f]]) > 1) {
        aggregated_stock_waa[[f]] <- colMeans(aggregated_stock_waa[[f]])
      } else {
        aggregated_stock_waa[[f]] <- mean(aggregated_stock_waa[[f]])
      }
    } else {
      aggregated_stock_waa[[f]] <- as.matrix(em_info$par_inputs$user_waa[ind, ])
      if (ncol(aggregated_stock_waa[[f]]) > 1) {
        aggregated_stock_waa[[f]] <- colMeans(aggregated_stock_waa[[f]])
      } else {
        aggregated_stock_waa[[f]] <- mean(aggregated_stock_waa[[f]])
      }
    }
  }
  
  
  
  # ---- Compute WAA Weights ---- #
  # fleet_weights <- list()
  # for (f in valid_fleets) {
  #   if (aggregate_catch_info$use_catch_weighted_waa) {
  #     agg_catch_filtered <- data$agg_catch[ind_em, which(fleet_pointer == f), drop = FALSE]
  #     
  #     # Ensure weight vector has correct sum
  #     mean_catch <- if (is.matrix(agg_catch_filtered)) colMeans(agg_catch_filtered, na.rm = TRUE) else mean(agg_catch_filtered, na.rm = TRUE)
  #     fleet_weights[[f]] <- mean_catch / sum(mean_catch, na.rm = TRUE)  # Normalize weights so they sum to 1
  #   } else {
  #     fleet_weights[[f]] <- rep(1 / length(valid_fleets), length(valid_fleets))  # Equal weights
  #   }
  # }
  # 
  # index_weights <- list()
  # for (f in valid_indices) {
  #   if (aggregate_index_info$use_catch_weighted_waa) {
  #     agg_indices_filtered <- data$agg_indices[ind_em, which(index_pointer == f), drop = FALSE]
  #     
  #     # Ensure weight vector has correct sum
  #     mean_index <- if (is.matrix(agg_indices_filtered)) colMeans(agg_indices_filtered, na.rm = TRUE) else mean(agg_indices_filtered, na.rm = TRUE)
  #     index_weights[[f]] <- mean_index / sum(mean_index, na.rm = TRUE)  # Normalize weights
  #   } else {
  #     index_weights[[f]] <- rep(1 / length(valid_indices), length(valid_indices))  # Equal weights
  #   }
  # }
  # 
  # # ---- Aggregate WAA ---- #
  # aggregated_fleet_waa = list()
  # for (f in unique(valid_fleets)) {
  #   if(is.null(em_info$par_inputs$user_waa)) {
  #     aggregated_fleet_waa[[f]] <- em_info$basic_info$waa[which(fleet_pointer == f),1,]
  #     if (is.matrix(aggregated_fleet_waa[[f]])) aggregated_fleet_waa[[f]] <- colSums(em_info$basic_info$waa[which(fleet_pointer == f),1, ]*fleet_weights[[f]])
  #   } else {
  #     aggregated_fleet_waa[[f]] <- em_info$par_inputs$user_waa[which(fleet_pointer == f), ]
  #     if (is.matrix(aggregated_fleet_waa[[f]])) aggregated_fleet_waa[[f]] <- colSums(em_info$par_inputs$user_waa[which(fleet_pointer == f), ]*fleet_weights[[f]])
  #   }
  # }
  # 
  # aggregated_region_waa = list()
  # for (f in 1) { # Assuming n_regions = 1
  #   ind = length(fleet_pointer) + 1:em_info$basic_info$n_regions
  #   if(is.null(em_info$par_inputs$user_waa)) {
  #     aggregated_region_waa[[f]] <- em_info$basic_info$waa[ind,1, ]
  #     if (is.matrix(aggregated_region_waa[[f]])) aggregated_region_waa[[f]] <- colMeans(as.matrix(aggregated_region_waa[[f]]))
  #   } else {
  #     aggregated_region_waa[[f]] <- em_info$par_inputs$user_waa[, ]
  #     if (is.matrix(aggregated_region_waa[[f]])) aggregated_region_waa[[f]] <- colMeans(as.matrix(aggregated_region_waa[[f]]))
  #   }
  # }
  # 
  # aggregated_index_waa = list()
  # for (f in unique(valid_indices)) {
  #   ind = length(fleet_pointer) + 1:em_info$basic_info$n_regions + which(index_pointer == f)
  #   if(is.null(em_info$par_inputs$user_waa)) {
  #     aggregated_index_waa[[f]] <- em_info$basic_info$waa[ind,1, ]
  #     if (is.matrix(aggregated_index_waa[[f]])) aggregated_index_waa[[f]] <- colMeans(as.matrix(aggregated_index_waa[[f]]))
  #   } else {
  #     aggregated_index_waa[[f]] <- em_info$par_inputs$user_waa[ind, ]
  #     if (is.matrix(aggregated_index_waa[[f]])) aggregated_index_waa[[f]] <- colSums(em_info$par_inputs$user_waa[which(index_pointer == f), ]*index_weights[[f]])
  #   }
  # }
  # 
  # aggregated_stock_waa = list()
  # for (f in 1) {
  #   ind = length(fleet_pointer) + em_info$basic_info$n_regions + length(index_pointer) + 1:em_info$basic_info$n_stocks
  #   if(is.null(em_info$par_inputs$user_waa)) {
  #     aggregated_stock_waa[[f]] <- em_info$basic_info$waa[ind,1, ]
  #     if (is.matrix(aggregated_stock_waa[[f]])) aggregated_stock_waa[[f]] <- colMeans(as.matrix(aggregated_stock_waa[[f]]))
  #   } else {
  #     aggregated_stock_waa[[f]] <- em_info$par_inputs$user_waa[ind, ]
  #     if (is.matrix(aggregated_stock_waa[[f]])) aggregated_stock_waa[[f]] <- colMeans(as.matrix(aggregated_stock_waa[[f]]))
  #   }
  # }
  
  # Store the aggregated WAA back into `em_info`
  em_info$par_inputs$user_waa <- do.call(rbind,c(aggregated_fleet_waa,aggregated_region_waa,aggregated_index_waa,aggregated_stock_waa))
  em_info$basic_info$waa = em_info$par_inputs$user_waa
  
  # Average MAA to get global MAA
  maturity <- em_info$par_inputs$user_maturity
  if(is.vector(maturity)) em_info$par_inputs$user_maturity = maturity
  if(is.matrix(maturity)) em_info$par_inputs$user_maturity = colMeans(maturity)
  
  em_info$par_inputs$n_regions = n_regions = 1
  em_info$par_inputs$n_stocks = n_stocks = 1 
  em_info$par_inputs$n_fleets = n_fleets
  em_info$par_inputs$n_indices = n_indices
  
  em_info$par_inputs$agg_catch <- agg_catch
  em_info$par_inputs$catch_paa <- agg_catch_paa  
  
  em_info$par_inputs$agg_indices <- agg_indices
  em_info$par_inputs$index_paa <- agg_index_paa  
  
  em_info$basic_info$waa_pointer_fleets   <- n_fleets
  em_info$basic_info$waa_pointer_totcatch <- n_fleets + n_regions
  em_info$basic_info$waa_pointer_indices  <- (n_fleets + n_regions + 1):(n_fleets + n_regions + n_indices)
  em_info$basic_info$waa_pointer_ssb      <- (n_fleets + n_regions + n_indices + 1):(n_fleets + n_regions + n_indices + 1)
  em_info$basic_info$waa_pointer_M        <- em_info$basic_info$waa_pointer_ssb
  
  # return(list(em_info = em_info, agg_catch = agg_catch, catch_paa = catch_paa, agg_indices = agg_indices, index_paa = index_paa))
  return(em_info)
}