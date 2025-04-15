#' Aggregate Catch and Index Data for Estimation Model
#'
#' This function processes and aggregates catch and index data from the operating model
#' for use in the estimation model (`em_info`). It accounts for multiple fleets and indices
#' and adjusts aggregation based on fleet and index pointers.
#'
#' @param data List. The operating model output, including:
#'   \itemize{
#'     \item `$agg_catch` - Aggregated catch data.
#'     \item `$catch_paa` - Proportions-at-age for catch.
#'     \item `$agg_indices` - Aggregated index data.
#'     \item `$index_paa` - Proportions-at-age for indices.
#'   }
#' @param em_info List. The estimation model inputs, including `par_inputs` and fleet/index settings.
#' @param aggregate_catch_info List. Contains fleet aggregation details, including:
#'   \itemize{
#'     \item `$fleet_pointer` - Vector mapping original fleets to aggregated fleets (0 = exclude).
#'   }
#' @param aggregate_index_info List. Contains index aggregation details, including:
#'   \itemize{
#'     \item `$index_pointer` - Vector mapping original indices to aggregated indices (0 = exclude).
#'   }
#' @param ind_em Vector. Indexing the estimation model years within the dataset.
#' @param n_fleets Integer. The number of fleets in the estimation model.
#' @param n_indices Integer. The number of indices in the estimation model.
#'
#' @return List. A modified `info` list containing aggregated `catch_info` and `index_info`:
#'   \itemize{
#'     \item `$catch_info$agg_catch` - Aggregated catch data for the estimation model.
#'     \item `$catch_info$catch_paa` - Aggregated catch-at-age data.
#'     \item `$index_info$agg_indices` - Aggregated index data.
#'     \item `$index_info$index_paa` - Aggregated index-at-age data.
#'   }
#'
#' @export

aggregate_em_data <- function(data, info, aggregate_catch_info, aggregate_index_info, ind_em, n_fleets, n_indices) {
  
  # ---- Aggregate Catch Data ---- #
  if (n_fleets == 1) {
    info$catch_info$agg_catch <- data$agg_catch[ind_em, , drop = FALSE]
    n = ncol(info$catch_info$agg_catch)
    info$catch_info$agg_catch <- matrix(rowSums(info$catch_info$agg_catch), ncol = 1)
    info$catch_info$catch_paa <- data$catch_paa[, ind_em, , drop = FALSE]
    catch <- data$agg_catch[ind_em, , drop = FALSE]
    ratio <- data$catch_paa[, ind_em, ]
    
    if (n > 1) {
      result <- 0
      for (i in 1:dim(ratio)[1]) {
        tmp <- ratio[i, , ] * catch[, i]
        result <- result + tmp
      }
      result <- t(apply(result, 1, function(row) row / sum(row)))
      info$catch_info$catch_paa <- array(result, dim = c(1, nrow(result), ncol(result)))
    } 
    
  } else if (n_fleets > 1) {
    
    if (is.null(aggregate_catch_info$fleet_pointer)) {
      warnings("aggregate_catch_info$fleet_pointer is not specified!")
      fleet_pointer <- rep(1, em_info$par_inputs$n_fleets)
    } else {
      fleet_pointer <- aggregate_catch_info$fleet_pointer
    }
    
    valid_pointer <- unique(fleet_pointer[fleet_pointer > 0])
    
    info$catch_info$agg_catch <- matrix(NA, length(ind_em), length(valid_pointer))
    info$catch_info$catch_paa <- array(NA, dim = c(length(valid_pointer), length(ind_em), info$basic_info$n_ages))
    
    for (f in valid_pointer) {
      agg_catch.tmp <- data$agg_catch[ind_em, which(fleet_pointer == f), drop = FALSE]
      n = ncol(agg_catch.tmp)
      agg_catch.tmp <- matrix(rowSums(agg_catch.tmp), ncol = 1)
      catch_paa.tmp <- data$catch_paa[which(fleet_pointer == f), ind_em, , drop = FALSE]
      catch <- data$agg_catch[ind_em, which(fleet_pointer == f), drop = FALSE]
      ratio <- data$catch_paa[which(fleet_pointer == f), ind_em, , drop = FALSE]
      
      if (n > 1) {
        result <- 0
        for (i in 1:dim(ratio)[1]) {
          tmp <- ratio[i, , ] * catch[, i]
          result <- result + tmp
        }
        result <- t(apply(result, 1, function(row) row / sum(row)))
        catch_paa.tmp <- array(result, dim = c(1, nrow(result), ncol(result)))
      } 

      info$catch_info$agg_catch[, f] <- agg_catch.tmp
      info$catch_info$catch_paa[f, , ] <- catch_paa.tmp[1, , ]
    }
  }
  
  # ---- Aggregate Index Data ---- #
  if (n_indices == 1) {
    info$index_info$agg_indices <- data$agg_indices[ind_em, , drop = FALSE]
    n = ncol(info$index_info$agg_indices)
    info$index_info$agg_indices <- matrix(rowSums(info$index_info$agg_indices), ncol = 1)
    info$index_info$index_paa <- data$index_paa[, ind_em, , drop = FALSE]
    catch <- data$agg_indices[ind_em, , drop = FALSE]
    ratio <- data$index_paa[, ind_em, ]
    
    if (n > 1) {
      result <- 0
      for (i in 1:dim(ratio)[1]) {
        tmp <- ratio[i, , ] * catch[, i]
        result <- result + tmp
      }
      result <- t(apply(result, 1, function(row) row / sum(row)))
      info$index_info$index_paa <- array(result, dim = c(1, nrow(result), ncol(result)))
    }
    
  } else if (n_indices > 1) {
    
    if (is.null(aggregate_index_info$index_pointer)) {
      warnings("aggregate_index_info$index_pointer is not specified!")
      index_pointer <- rep(1, em_info$par_inputs$n_indices)
    } else {
      index_pointer <- aggregate_index_info$index_pointer
    }
    
    valid_pointer <- unique(index_pointer[index_pointer > 0])
    
    info$index_info$agg_indices <- matrix(NA, length(ind_em), length(valid_pointer))
    info$index_info$index_paa <- array(NA, dim = c(length(valid_pointer), length(ind_em), info$basic_info$n_ages))
    
    for (f in valid_pointer) {
      agg_indices.tmp <- data$agg_indices[ind_em, which(index_pointer == f), drop = FALSE]
      n <- ncol(agg_indices.tmp)
      agg_indices.tmp <- matrix(rowSums(agg_indices.tmp), ncol = 1)
      index_paa.tmp <- data$index_paa[which(index_pointer == f), ind_em, , drop = FALSE]
      catch <- data$agg_indices[ind_em, which(index_pointer == f), drop = FALSE]
      ratio <- data$index_paa[which(index_pointer == f), ind_em, , drop = FALSE]
      
      if (n > 1) {
        result <- 0
        for (i in 1:dim(ratio)[1]) {
          tmp <- ratio[i, , ] * catch[, i]
          result <- result + tmp
        }
        result <- t(apply(result, 1, function(row) row / sum(row)))
        index_paa.tmp <- array(result, dim = c(1, nrow(result), ncol(result)))
      }
      info$index_info$agg_indices[, f] <- agg_indices.tmp
      info$index_info$index_paa[f, , ] <- index_paa.tmp[1, , ]
    }
  }
  
  return(info)
}
