#' Update the operating model and generate data
#' 
#' Function to update F in the operating model (see \code{\link{update_om_F}}) and generate data given the updated F.
#' 
#' @param om Operating model with years including burn-in + feedback years.
#' @param interval.info Catch advice for a number of years projected from the estimation model.
#'   \itemize{
#'     \item \code{$years} Projection years.
#'     \item \code{$catch} Matrix (n_region x n_years) projected catch.
#'   }
#' @param seed Seed used to generate data.
#' @param random A vector of processes that are treated as random effects in the operating model
#' 
#' @return An operating model with simulated data and updated F time series.
#'   
#' @export
#'
#' @seealso \code{\link{get_F_from_Catch_region}}
#' 
update_om_fn <- function(om, interval.info = NULL, seed = 123, random = "log_NAA") {
  if(!is.null(interval.info)){
    # Iterative update F in the OM using get_F_from_Catch_region function
    t <- 0
    for (y in interval.info$years) {
      year <- which(om$years == y)
      t <- t + 1
      rep <- om$rep
      
      # Solve for fishing mortality rate (F)
      Fsolve <- get_F_from_Catch_region(
        Catch = interval.info$catch[t,], 
        NAA = rep$NAA[,,year,], 
        log_M = log(rep$MAA[,,year,]), 
        mu = rep$mu[,,,year,,], 
        L = rep$L[year,], 
        # sel = (rep$FAA / max(rep$FAA))[,year,]/max(max(rep$FAA)), 
        sel = rep$FAA[,year,]/max(rep$FAA[,year,]),
        fracyr_season = om$input$data$fracyr_seasons, 
        fleet_regions = om$input$data$fleet_regions, 
        fleet_seasons = om$input$data$fleet_seasons, 
        can_move = om$input$data$can_move, 
        mig_type = om$input$data$mig_type, 
        waacatch = om$input$data$waa[(om$input$data$n_fleets + 1):(om$input$data$n_fleets + om$input$data$n_regions), year,], 
        trace = TRUE, 
        F_init = rep(0.5, om$input$data$n_fleets)
      )
      
      cat(paste0("\nFsolve: ", round(Fsolve, 3), " in Year ", y, "\n"))
      
      # Update F parameters in the operating model
      om$input$par$F_pars[year,] <- log(Fsolve)
      
      # Names of observation data to update
      obs_names <- c("agg_indices", "agg_catch", "catch_paa", "index_paa", "Ecov_obs", "obsvec")
      
      # Set seed for reproducibility
      set.seed(seed)
      
      # Simulate the population and observations
      om_sim <- om$simulate(complete = TRUE)
      
      # Update the simulated data in the operating model
      om$input$data[obs_names] <- om_sim[obs_names]
      
      # Update the parameters in the operating model
      om$input$par[random] <- om_sim[random]
      
      # Fit the WHAM model without actually performing the fit (do.fit = FALSE)
      om <- fit_wham(om$input, do.fit = FALSE, MakeADFun.silent = TRUE)
      
      }
    
    } else {
      
      # Names of observation data to update
      obs_names <- c("agg_indices", "agg_catch", "catch_paa", "index_paa", "Ecov_obs", "obsvec")
      
      set.seed(seed)
      
      om_sim = om$simulate(complete=TRUE) #resimulate the population and observations
      
      om$input$data[obs_names] = om_sim[obs_names] #update any simulated data
      
      om$input$par[random] = om_sim[random] #update any simulated random effects
      
      # reset the om
      om <- fit_wham(om$input, do.fit = FALSE, MakeADFun.silent = TRUE)
    }
  
  return(om)
}
