#' Update the Operating Model and Generate Simulated Data
#'
#' This function updates fishing mortality (F) in the operating model, projects the population forward,
#' and generates new simulated data based on the updated F time series.
#'
#' @param om A fitted operating model that includes both burn-in and feedback years.
#' @param interval.info A list containing projected catch advice for future years from the estimation model. It includes:
#'   \describe{
#'     \item{\code{years}}{Vector of projection years.}
#'     \item{\code{catch}}{Matrix (n_regions x n_years) of projected catch.}
#'   }
#' @param seed Integer. Seed used for generating simulated data to ensure reproducibility (default = 123).
#' @param random Character vector of process names that are treated as random effects in the operating model (default = "log_NAA").
#' @param method Character. Optimization method used for solving F. Options include:
#'   \itemize{
#'     \item \code{"nlminb"} (default): Uses `nlminb` optimization.
#'     \item \code{"BFGS"}: Uses `optim` with the BFGS method.
#'   }
#' @param by_fleet Logical. If TRUE, estimates F separately for each fleet. If FALSE, estimates a single global F (default = FALSE).
#' @param do.brps Logical. If TRUE, calculates reference points in the operating model (default = FALSE).
#' 
#' @return An updated operating model (`om`) with:
#'   \itemize{
#'     \item Updated F time series.
#'     \item Simulated observation data.
#'     \item Updated random effect parameters.
#'   }
#'
#' @export
#'
#' @seealso \code{\link{get_F_from_Catch}}, \code{\link{update_om_F}}
#'
update_om_fn <- function(om, interval.info = NULL, seed = 123, random = "log_NAA", method = "nlminb", by_fleet = TRUE, do.brps = FALSE) {
  
  if(!is.null(interval.info)){
    # Iterative update F in the OM using get_F_from_Catch_region function
    t <- 0
    for (y in interval.info$years) {
      
      year <- which(om$years == y)
      t <- t + 1
      Catch = interval.info$catch[t,]
      
      rep <- om$rep #generate the reported values given the parameters
      
      cat("\nNow calculating fleet-specific F in year ", y, "\n")
      
      Fsolve <- get_F_from_Catch(om, Catch, year, method = "nlminb", by_fleet = by_fleet)
      
      cat("\nFishing mortality is ", Fsolve, "\n")
      cat("\nNow updating F in OM for year ", y, "\n")
      
      om <- update_om_F(om, year, Fsolve)
      
      # Names of observation data to update
      obs_names <- c("agg_indices", "agg_catch", "catch_paa", "index_paa", "Ecov_obs", "obsvec")
      
      # Set seed for reproducibility
      set.seed(seed)
      
      cat("\nNow simulating data for year ", y, "\n")
      
      # Simulate the population and observations
      om_sim <- om$simulate(complete = TRUE)
      
      # Update the simulated data in the operating model
      om$input$data[obs_names] <- om_sim[obs_names]
      
      # Update the parameters in the operating model
      om$input$par[random] <- om_sim[random]
      
      cat("\nNow projecting data for years after ", y, "\n")
      
      # Fit the WHAM model without actually performing the fit (do.fit = FALSE)
      om <- fit_wham(om$input, do.fit = FALSE, do.brps = do.brps, MakeADFun.silent = TRUE)
      
    }
    
  } else {
    
    # Names of observation data to update
    obs_names <- c("agg_indices", "agg_catch", "catch_paa", "index_paa", "Ecov_obs", "obsvec")
    
    set.seed(seed)
    
    om_sim = om$simulate(complete=TRUE) #resimulate the population and observations
    
    om$input$data[obs_names] = om_sim[obs_names] #update any simulated data
    
    om$input$par[random] = om_sim[random] #update any simulated random effects
    
    cat("\nNow simulating data")
    # reset the om
    om <- fit_wham(om$input, do.fit = FALSE, do.brps = do.brps, MakeADFun.silent = TRUE)
  }
  
  return(om)
}