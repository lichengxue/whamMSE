#' Update Index CV and Neff in Input List (with Optional Row Subsetting)
#'
#' Updates `input$data$index_Neff` and `input$data$agg_index_sigma` using provided index CV and Neff matrices.
#'
#' @param input List. Must contain `input$data`.
#' @param agg_index_sigma Matrix. Either full (n_years x n_indices) or subset (length(ind_em) x n_indices).
#' @param index_Neff Matrix. Same dimension as `agg_index_sigma`.
#' @param ind_em Integer vector or NULL. Row indices of years to update. If NULL, updates all rows.
#'
#' @return The updated `input` list.
#' @export
update_input_index_info <- function(input, agg_index_sigma, index_Neff, 
                                    remove.agg = FALSE, remove_agg_pointer = NULL, remove_agg_years = NULL,
                                    remove.paa = FALSE, remove_paa_pointer = NULL, remove_paa_years = NULL,
                                    ind_em = NULL) {
  # Validate dimensions
  if (!all(dim(agg_index_sigma) == dim(index_Neff))) {
    stop("agg_index_sigma and index_Neff must have the same dimensions.")
  }
  
  if (is.null(ind_em)) {
    # Update full matrices
    input$data$agg_index_sigma <- agg_index_sigma
    input$data$index_Neff <- index_Neff
  } else {
    # Check index validity
    # if (max(ind_em) > nrow(input$data$agg_index_sigma)) {
    #   stop("ind_em contains invalid row indices.")
    # }
    # if (length(ind_em) != nrow(agg_index_sigma)) {
    #   stop("Number of rows in agg_index_sigma and index_Neff must match length of ind_em.")
    # }
    
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
    
    if (remove_agg) {
      if (is.null(remove_agg_pointer)) stop("remove_paa_pointer must be specified!")
      if (is.null(remove_agg_years)) stop("remove_agg_years must be specified!")
      for (i in remove_agg_pointer) {
        input$data$agg_indices[remove_paa_years,i] = -1
      }
    }
    
    if (remove_paa) {
      if (is.null(remove_paa_pointer)) stop("remove_paa_pointer must be specified!")
      if (is.null(remove_paa_years)) stop("remove_paa_years must be specified!")
      for (i in remove_paa_pointer) {
        input$data$index_paa[i,remove_paa_years,] = 0
      }
    }

  }
  
  return(input)
}