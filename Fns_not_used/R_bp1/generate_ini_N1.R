#' Generate initial numbers-at-age
#' 
#' @param log_N1 a vector (length = n_stocks) of log numbers-at-age1
#' @param basic_info (optional) list specifying options for numbers-at-age random effects, initial parameter values, and recruitment model (see details)
#' @param ini.opt N1_model
#'       \describe{
#'          \item{"age-specific-fe"}{(default) age- and region-specific fixed effects parameters}
#'          \item{"equilibrium"}{2 fixed effects parameters: an initial recruitment and an instantaneous fishing mortality rate to generate an equilibruim abundance at age.}
#'          \item{"iid-re"}{(default) age- and region-specific iid random effects parameters. 2 parameters: mean and sd for log NAA}
#'          \item{"ar1-re"}{(default) age- and region-specific random effects parameters. 3 parameters: mean and sd, and cor for log NAA}
#'       }
#' @return an array of numbers-at-age (n_stocks x n_regions x n_ages)
#' 
#' @export
#' 
#' @examples
#' \dontrun{
#' basic_info <- generate_basic_info()
#' N1_pars <- generate_ini_N1(log_N1 = c(10.6,10), basic_info, ini.opt = "equilibrium") 
#' }
#' 
generate_ini_N1 <- function(log_N1 = c(10.6, 10), 
                            basic_info,
                            ini.opt = "equilibrium") {
  
  # Validate inputs
  if (!is.list(basic_info)) stop("basic_info must be a list!")
  
  # Extract parameters from basic_info
  n_stocks  <- as.integer(basic_info$n_stocks)
  n_regions <- as.integer(basic_info$n_regions)
  n_ages    <- as.integer(basic_info$n_ages)
  NAA_where <- basic_info$NAA_where
  
  # Validate lengths
  if (length(log_N1) != n_stocks) stop("Length of log_N1 must be equal to n_stocks!")
  
  # Initialize the array for numbers-at-age
  N1_pars <- array(1, dim = c(n_stocks, n_regions, n_ages))
  
  for (s in 1:n_stocks) {
    if (ini.opt == "equilibrium") { # Equilibrium assumption, 2 pars per stock
      N1_pars[s, s, 1:2] <- c(exp(log_N1[s]), exp(log(0.1)))
    } else if (ini.opt == "age-specific-fe") { # Age-specific fixed effects
      if (is.null(NAA_where)) stop("NAA_where must be specified for age-specific-fe option!")
      
      init_NAA <- exp(log_N1[s]) * exp(-(0:(n_ages-1)) * 0.2)
      
      for (r in 1:n_regions) {
        for (a in 1:n_ages) {
          if (NAA_where[s, r, a] == 1) {
            N1_pars[s, r, a] <- init_NAA[a]
          }
        }
      }
    } else {
      stop("Unknown ini.opt value: ", ini.opt)
    }
  }
  
  return(N1_pars)
}
