#' Update Catch CV and Neff in Input List (with Optional Row Subsetting)
#'
#' Updates `input$data$catch_Neff` and `input$data$agg_catch_sigma` using provided catch CV and Neff matrices.
#'
#' @param input List. Must contain `input$data`.
#' @param agg_catch_sigma Matrix. Either full (n_years x n_fleets) or subset (length(ind_em) x n_fleets).
#' @param catch_Neff Matrix. Same dimension as `agg_catch_sigma`.
#' @param ind_em Integer vector or NULL. Row indices of years to update. If NULL, updates all rows.
#'
#' @return The updated `input` list.
#' @export
#' 
update_input_catch_info <- function(input, agg_catch_sigma, catch_Neff, ind_em = NULL) {
  # Validate dimensions
  if (!all(dim(agg_catch_sigma) == dim(catch_Neff))) {
    stop("agg_catch_sigma and catch_Neff must have the same dimensions.")
  }
  
  if (is.null(ind_em)) {
    # Update full matrices
    input$data$agg_catch_sigma <- agg_catch_sigma
    input$data$catch_Neff <- catch_Neff
  } else {
    # Check index validity
    # if (max(ind_em) > nrow(input$data$agg_catch_sigma)) {
    #   stop("ind_em contains invalid row indices.")
    # }
    # if (length(ind_em) != nrow(agg_catch_sigma)) {
    #   stop("Number of rows in agg_catch_sigma and catch_Neff must match length of ind_em.")
    # }
    
    # Update only selected rows

    if (ncol(agg_catch_sigma) == 1) {
      input$data$agg_catch_sigma = matrix(agg_catch_sigma[ind_em, ], byrow = T)
    } else {
      input$data$agg_catch_sigma <- agg_catch_sigma[ind_em, ]
    }
    
    if (ncol(catch_Neff) == 1) {
      input$data$catch_Neff = matrix(catch_Neff[ind_em, ], byrow = T)
    } else {
      input$data$catch_Neff <- catch_Neff[ind_em, ]
    }
    
  }
  
  return(input)
}