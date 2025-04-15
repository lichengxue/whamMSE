#' Update Index Information in Input List (with Optional Row Subsetting)
#'
#' Updates \code{input$data$index_Neff} and \code{input$data$agg_index_sigma} using provided index CV and Neff matrices.
#' Optionally removes index data (aggregated or PAA) for specified pointers and years, and updates internal flags.
#'
#' @param input List. A WHAM input list that must contain \code{input$data}.
#' @param agg_index_sigma Numeric matrix. Aggregated index standard deviations (CVs) to update. 
#' Can be full size (\code{n_years × n_indices}) or subset rows (\code{length(ind_em) × n_indices}).
#' @param index_Neff Numeric matrix. Effective sample size for index data. Must match dimensions of \code{agg_index_sigma}.
#' @param remove_agg Logical. Whether to remove aggregated index observations for specific pointers/years.
#' @param remove_agg_pointer Integer vector. Index pointers (columns) for which aggregated index data should be removed.
#' @param remove_agg_years Integer vector or matrix. Years (rows) to remove aggregated index data. If matrix, should be [years × pointers].
#' @param remove_paa Logical. Whether to remove index PAA observations for specific pointers/years.
#' @param remove_paa_pointer Integer vector. Index pointers (columns) for which index PAA data should be removed.
#' @param remove_paa_years Integer vector or matrix. Years (rows) to remove index PAA data. If matrix, should be [years × pointers].
#' @param ind_em Integer vector or NULL. If provided, only updates rows corresponding to these indices. If NULL, updates all rows.
#'
#' @return The modified \code{input} list with updated \code{agg_index_sigma}, \code{index_Neff}, and associated index flags.
#'
#' @details
#' If \code{ind_em} is specified, only those rows will be updated in the input matrices. 
#' Removal of index or PAA observations is handled via internal flags and array zeroing.
#' The function automatically calls \code{set_osa_obs()} at the end to update OSA flags.
#'
#' @export
#' 
update_input_index_info <- function(input, agg_index_sigma, index_Neff, 
                                    remove_agg = FALSE, remove_agg_pointer = NULL, remove_agg_years = NULL,
                                    remove_paa = FALSE, remove_paa_pointer = NULL, remove_paa_years = NULL,
                                    ind_em = NULL) {
  # Validate dimensions
  if (!all(dim(agg_index_sigma) == dim(index_Neff))) {
    stop("agg_index_sigma and index_Neff must have the same dimensions.")
  }
  
  if (is.null(ind_em)) {
    # Update full matrices
    input$data$agg_index_sigma <- agg_index_sigma
    input$data$index_Neff <- index_Neff
    
    if (remove_agg) {
      if (is.null(remove_agg_pointer)) stop("remove_paa_pointer must be specified!")
      if (is.null(remove_agg_years)) stop("remove_agg_years must be specified!")
      for (i in remove_agg_pointer) {
        if(is.vector(remove_agg_years)) {
          y = remove_agg_years
          input$data$use_indices[y,i] = 0
          input$data$agg_indices[y,i] = 0
        }
        
        if(is.matrix(remove_agg_years)) {
          y = remove_agg_years[,i]
          input$data$use_indices[y,i] = 0
          input$data$agg_indices[y,i] = 0
        }
      }
      
      # input$data$agg_indices[y,i] = -1
      # y = ind_em %in% remove_agg_years 
      # input$data$agg_indices[y,i] = -1
      # input = wham:::set_osa_obs(input)
      # input = fit_wham(input)
    }
    
    if (remove_paa) {
      if (is.null(remove_paa_pointer)) stop("remove_paa_pointer must be specified!")
      if (is.null(remove_paa_years)) stop("remove_paa_years must be specified!")
      
      for (j in remove_paa_pointer) {
        if(is.vector(remove_paa_years)) {
          y = remove_paa_years
          input$data$use_index_paa[y,j] = 0
          input$data$index_paa[j,y,] = 0
        }
        
        if(is.matrix(remove_paa_years)) {
          y = remove_paa_years[,j]
          input$data$use_index_paa[y,j] = 0
          input$data$index_paa[j,y,] = 0
        }
      }
    }
    
    
    # input <- wham::set_indices(input)
    input <- set_osa_obs(input)
    
    
  } else {
    
    # Update only selected rows
    if (ncol(agg_index_sigma) == 1) {
      input$data$agg_index_sigma = matrix(agg_index_sigma[ind_em, ], byrow = T)
    } else {
      input$data$agg_index_sigma <- agg_index_sigma[ind_em, ]
    }
    
    if (ncol(index_Neff) == 1) {
      input$data$index_Neff = matrix(index_Neff[ind_em, ], byrow = T)
    } else {
      input$data$index_Neff <- index_Neff[ind_em, ]
    }
    
    all_zero <- apply(input$data$index_paa, c(1, 2), function(x) all(x == 0))
    
    if(sum(all_zero > 0)) {
      id = which(all_zero, arr.ind = TRUE)
      for (i in 1:nrow(id)) {
        input$data$use_index_paa[id[i,2],id[i,1]] = 0
        # input$data$use_indices[id[i,2],id[i,1]] = 0
      }
    }
    
    all_zero <- apply(input$data$agg_indices, c(1, 2), function(x) all(x == 0))
    
    if(sum(all_zero > 0)) {
      id = which(all_zero, arr.ind = TRUE)
      for (i in 1:nrow(id)) {
        input$data$use_index_paa[id[i,1],id[i,2]] = 0
        input$data$use_indices[id[i,1],id[i,2]] = 0
      }
    }
    
    input <- set_osa_obs(input)
  }
  
  return(input)
}