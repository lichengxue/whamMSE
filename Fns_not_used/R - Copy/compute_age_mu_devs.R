#' Compute Age-Specific Movement Deviations and Baseline Movement Rates
#'
#' This function calculates the deviation of movement rates at age from the minimum movement rate
#' for each region-to-region pair, ensuring that all deviations are greater than or equal to zero.
#' The function also saves the baseline movement rates, which can be used to reconstruct the
#' original movement matrix.
#'
#' @param move_matrix A 3D array of movement rates with dimensions `(n_regions, n_regions-1, n_ages)`.
#'        Each entry represents the movement rate from one region to another across different ages.
#' @param n_stocks An integer specifying the number of stocks for which deviations should be computed.
#'
#' @return A list containing:
#'   \describe{
#'     \item{`age_mu_devs`}{A 4D array of movement deviations with dimensions `(n_stocks, n_regions, n_regions-1, n_ages)`,
#'                           ensuring all values are non-negative.}
#'     \item{`baseline_move`}{A 3D array `(n_regions, n_regions-1, n_ages)` representing the minimum movement
#'                           rate for each region-pair across ages, used as a baseline.}
#'     \item{`mean_vals`}{A 4D array `(n_stocks, n_seasons, n_regions, n_regions-1)` representing the minimum movement
#'                           rate for each region-pair across stocks and season, used as a baseline mean movement rate.}
#'   }
#' @export
#' 
compute_age_mu_devs <- function(move_matrix, n_stocks, n_seasons) {
  # Validate input dimensions
  dims <- dim(move_matrix)
  if (length(dims) != 3) stop("move_matrix must be a 3D array with dimensions (n_regions, n_regions-1, n_ages)")
  
  n_regions <- dims[1]
  n_regions_minus1 <- dims[2]
  n_ages <- dims[3]
  
  # Find the minimum movement rate for each region-pair across ages
  min_move_matrix <- apply(move_matrix, c(1, 2), min) 
  
  move_matrix_expand <- array(rep(move_matrix, each = n_stocks),
                       dim = c(n_stocks, n_regions, n_regions_minus1, n_ages))
  
  mean_vals <- array(rep(min_move_matrix, each = n_stocks * n_seasons),
                     dim = c(n_stocks, n_seasons, n_regions, n_regions-1))
  
  baseline_move <- array(0, dim = c(n_stocks, n_seasons, n_regions, n_regions-1, n_ages))
  for (i in 1:n_ages) {
    baseline_move[,,,,i] <- mean_vals
  }

  baseline_move <- baseline_move[,1,,,]
  
  # Compute age_mu_devs (ensuring all values are ≥ 0)
  
  age_mu_devs <- move_matrix_expand - baseline_move
  
  # Return both age_mu_devs and baseline_move
  return(list(
    age_mu_devs = age_mu_devs,
    mean_vals = mean_vals
  ))
}