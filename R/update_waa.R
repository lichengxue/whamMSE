#' Update Weight-at-Age (WAA) Data in Input List
#'
#' This function updates the weight-at-age (WAA) data and associated pointers in
#' the input list based on provided `waa_info`. It replaces `input$data$waa`
#' and ensures that relevant pointers (`waa_pointer_ssb`, `waa_pointer_M`, etc.)
#' are correctly assigned.
#'
#' @param input A list containing `data` and optional `log` elements. `data` must include
#'   stock, fleet, and index information.
#' @param waa_info A list containing the weight-at-age (WAA) array and associated pointers.
#'   The list may include:
#'   \itemize{
#'     \item `waa`: A 3D array `(n_sources, n_years, n_ages)` representing weight-at-age data.
#'     \item `waa_pointer_ssb`: A numeric vector specifying WAA matrix indices for SSB.
#'     \item `waa_pointer_M`: A numeric vector specifying WAA matrix indices for natural mortality.
#'     \item `waa_pointer_fleets`: A numeric vector specifying WAA matrix indices for fleets.
#'     \item `waa_pointer_indices`: A numeric vector specifying WAA matrix indices for indices.
#'   }
#'
#' @return The modified `input` list with updated `data$waa` and associated pointers.
#' @export
#' 
update_waa <- function(input, waa_info = NULL) {
  data <- input$data  # Extract existing data
  
  # Ensure waa_info is provided
  if (!is.null(waa_info$waa)) {
    data$waa <- waa_info$waa  # Assign WAA
    
    # Validate dimensions
    dim_waa <- dim(data$waa)
    if (length(dim_waa) != 3)
      stop("waa_info$waa must be a 3D array: (n_sources, n_years, n_ages)")
    
    # Assign waa pointers if available
    if (!is.null(waa_info$waa_pointer_ssb)) {
      if (length(waa_info$waa_pointer_ssb) != data$n_stocks) {
        stop("Length of waa_pointer_ssb must match number of stocks.")
      }
      data$waa_pointer_ssb <- waa_info$waa_pointer_ssb
    } else {
      data$waa_pointer_ssb <- rep(1, data$n_stocks)  # Default if missing
    }
    
    if (!is.null(waa_info$waa_pointer_M)) {
      if (length(waa_info$waa_pointer_M) != data$n_stocks) {
        stop("Length of waa_pointer_M must match number of stocks.")
      }
      data$waa_pointer_M <- waa_info$waa_pointer_M
    } else {
      data$waa_pointer_M <- data$waa_pointer_ssb  # Default to SSB pointer
    }
    
    if (!is.null(waa_info$waa_pointer_totcatch)) {
      if (length(waa_info$waa_pointer_totcatch) != data$n_stocks) {
        stop("Length of waa_pointer_totcatch must match number of stocks.")
      }
      data$waa_pointer_totcatch <- waa_info$waa_pointer_totcatch
    } else {
      data$waa_pointer_totcatch <- data$waa_pointer_ssb  # Default to SSB pointer
    }
    
    if (!is.null(waa_info$waa_pointer_fleets)) {
      data$waa_pointer_fleets <- waa_info$waa_pointer_fleets
    }
    
    if (!is.null(waa_info$waa_pointer_indices)) {
      data$waa_pointer_indices <- waa_info$waa_pointer_indices
    }
    
    print("Updated input$data$waa with waa_info")
  } else {
    stop("waa_info$waa is missing, cannot update input$data$waa.")
  }
  
  # Update input$data with modified data
  input$data <- data
  input$options$waa <- waa_info
  
  return(input)
}