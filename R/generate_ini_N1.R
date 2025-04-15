#' Generate initial numbers-at-age
#' 
#' @param basic_info list specifying options for numbers-at-age random effects, initial parameter values, and recruitment model (see details)
#' @param ini.opt N1_model
#' @param log_N1 a vector (length = n_stocks) of log number-at-age 1 at equilibrium
#' @param log_N1_F a vector (length = n_stocks) of fishing mortality at equilibrium
#'       \describe{
#'          \item{"age-specific-fe"}{(default) age- and region-specific fixed effects parameters}
#'          \item{"equilibrium"}{2 fixed effects parameters: an initial recruitment and an instantaneous fishing mortality rate to generate an equilibruim abundance at age.}
#'          \item{"iid-re"}{(default) age- and region-specific iid random effects parameters. 2 parameters: mean and sd for log NAA}
#'          \item{"ar1-re"}{(default) age- and region-specific random effects parameters. 3 parameters: mean and sd, and cor for log NAA}
#'       }
#' @param init_NAA a matrix (n_ages x n_stocks) of initial numbers-at-age
#' @return an array of numbers-at-age (n_stocks x n_regions x n_ages)
#' 
#' @export

generate_ini_N1 <- function(basic_info, ini.opt = "equilibrium",
                            log_N1 = NULL, log_N1_F = NULL,
                            init_NAA = NULL
                            ) {
  
  # Validate inputs
  if (!is.list(basic_info)) stop("basic_info must be a list!")
  
  # Extract parameters from basic_info
  n_stocks  <- as.integer(basic_info$n_stocks)
  n_regions <- as.integer(basic_info$n_regions)
  n_ages    <- as.integer(basic_info$n_ages)
  NAA_where <- basic_info$NAA_where
  
  # Initialize the array for numbers-at-age
  N1_pars <- array(1, dim = c(n_stocks, n_regions, n_ages))
  
  for (s in 1:n_stocks) {
    if (ini.opt == "equilibrium") { # Equilibrium assumption, 2 pars per stock
      if (is.null(log_N1) && is.null(log_N1_F)) {
        N1_pars[s, s, 1:2] <- c(exp(10), exp(log(0.1)))
      } else if (!is.null(log_N1) && is.null(log_N1_F)) {
        N1_pars[s, s, 1:2] <- c(exp(log_N1[s]), exp(log(0.1)))
      } else if (is.null(log_N1) && !is.null(log_N1_F)) {
        N1_pars[s, s, 1:2] <- c(exp(log_N1[s]), exp(log(log_N1_F[s])))
      } else {
        N1_pars[s, s, 1:2] <- c(exp(log_N1[s]), exp(log(log_N1_F[s])))
      }
    } else if (ini.opt == "age-specific-fe") { # Age-specific fixed effects
      if (is.null(NAA_where)) warnings("NAA_where must be specified for age-specific-fe option!")
      if (is.null(NAA_where)) basic_info = generate_NAA_where(basic_info = basic_info)
      
      if(is.null(init_NAA)) {
        init_NAA_tmp <- exp(log_N1[s]) * exp(-(0:(n_ages-1)) * 0.2)
      } else {
        init_NAA_tmp <- init_NAA[,s]
        if(length(init_NAA_tmp) != n_ages) stop("Length of user-specified NAA must be equal to the number of ages!")
      }
      
      for (r in 1:n_regions) {
        for (a in 1:n_ages) {
          if (NAA_where[s, r, a] == 1) {
            N1_pars[s, r, a] <- init_NAA_tmp[a]
          }
        }
      }
    } else {
      stop("Unknown ini.opt value: ", ini.opt)
    }
  }
  
  return(N1_pars)
}
