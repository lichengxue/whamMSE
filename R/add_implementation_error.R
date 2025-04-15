#' Apply Implementation Error to Catch Advice
#'
#' Adds implementation error to catch advice in a Management Strategy Evaluation (MSE) framework.
#' Supports log-normal, normal, uniform, or constant multipliers to simulate realistic deviations
#' between recommended and realized catch.
#'
#' @param catch_advice A numeric matrix of catch advice with dimensions \code{[n_years, n_fleets]}.
#'   Can also be a vector or a matrix with dimensions \code{(1,1)} (i.e., a single value).
#' @param method Character. Type of implementation error distribution to use. Options are:
#'   \itemize{
#'     \item \code{"lognormal"}: Multiplies catch advice by \code{exp(N(mean, sd))}, where \code{sd} is derived from \code{cv} if not provided.
#'     \item \code{"normal"}: Adds error to catch via \code{1 + N(mean, sd)}.
#'     \item \code{"uniform"}: Multiplies catch advice by a uniform multiplier between \code{min} and \code{max}.
#'     \item \code{"constant"}: Multiplies catch advice by a fixed value (\code{constant_value}).
#'   }
#' @param mean Numeric. The mean of the error distribution. For log-normal, this is on the log scale.
#' @param cv Numeric. Coefficient of variation. Used only for log-normal error if \code{sd} is not provided.
#' @param sd Numeric. Standard deviation of the error distribution. Required for \code{"normal"}, optional for \code{"lognormal"}.
#' @param min Numeric. Minimum multiplier. Used only for \code{"uniform"}.
#' @param max Numeric. Maximum multiplier. Used only for \code{"uniform"}.
#' @param constant_value Numeric. A fixed multiplier applied to all entries if \code{method = "constant"}.
#'   If \code{constant_value} is not provided and \code{cv} is specified, the multiplier will default
#'   to the mean of a log-normal distribution with the given CV. If neither is provided, it defaults to 1.
#' @param seed Integer (optional). Random seed for reproducibility (used for stochastic methods).
#'
#' @return A numeric matrix of the same dimensions as \code{catch_advice}, representing the realized catch after applying implementation error.
#'
#' @examples
#' # Matrix: 5 years, 2 fleets
#' advice <- matrix(1000, nrow = 5, ncol = 2)
#' add_implementation_error(advice, method = "lognormal", cv = 0.2, seed = 123)
#'
#' # Single value (1x1 matrix)
#' add_implementation_error(matrix(1000, 1, 1), method = "constant", constant_value = 0.9)
#'
#' @export
add_implementation_error <- function(
    catch_advice,
    method = c("lognormal", "normal", "uniform", "constant"),
    mean = 0,
    cv = NULL,
    sd = NULL,
    min = NULL,
    max = NULL,
    constant_value = NULL,
    seed = NULL
) {
  method <- match.arg(method)
  
  # Convert catch_advice to matrix (if not already)
  catch_advice <- as.matrix(catch_advice)
  
  if (!is.null(seed)) set.seed(seed)
  
  # Matrix dimensions (will also work for 1x1)
  n_years <- nrow(catch_advice)
  n_fleets <- ncol(catch_advice)
  n_total <- n_years * n_fleets
  
  # Generate the error matrix based on method
  error_matrix <- switch(method,
                         "lognormal" = {
                           if (is.null(sd)) {
                             if (is.null(cv)) stop("Must provide either 'cv' or 'sd' for lognormal method.")
                             sd <- sqrt(log(1 + cv^2))
                           }
                           matrix(exp(rnorm(n_total, mean = mean, sd = sd)), nrow = n_years, ncol = n_fleets)
                         },
                         "normal" = {
                           if (is.null(sd)) stop("Must provide 'sd' for normal method.")
                           matrix(1 + rnorm(n_total, mean = mean, sd = sd), nrow = n_years, ncol = n_fleets)
                         },
                         "uniform" = {
                           if (is.null(min) || is.null(max)) stop("Must provide 'min' and 'max' for uniform method.")
                           matrix(runif(n_total, min = min, max = max), nrow = n_years, ncol = n_fleets)
                         },
                         "constant" = {
                           if (is.null(constant_value)) {
                             if (!is.null(cv)) {
                               constant_value <- exp(0.5 * log(1 + cv^2))  # Mean of log-normal with CV
                             } else {
                               constant_value <- 1  # Default: no error
                             }
                           }
                           matrix(constant_value, nrow = n_years, ncol = n_fleets)
                         }
  )
  
  realized_catch <- catch_advice * error_matrix
  return(realized_catch)
}
