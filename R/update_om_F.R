#' Update Fishing Mortality (F) in the Operating Model
#'
#' This function updates fishing mortality (F) in the operating model (`om`) for a specified year
#' and ensures consistency in future F projections.
#'
#' @param om A fitted operating model that includes burn-in and feedback years.
#' @param year Integer. The year for which fishing mortality is being updated.
#' @param Fsolve Numeric vector. Estimated fishing mortality values for each fleet.
#' 
#' @return An updated operating model (`om`) with:
#'   \itemize{
#'     \item Updated F time series.
#'     \item Adjusted F deviation parameters (\code{F_pars}) for model consistency.
#'   }
#'
#' @export
#'
#' @seealso \code{\link{get_F_from_Catch}}, \code{\link{update_om_fn}}
#'
update_om_F <- function(om, year, Fsolve) {
  
  if(length(Fsolve) > 1) {
    if(om$input$data$F_config==1) {
      if(year>1) {
        FAA_ym1 <- rbind(rep$FAA[,year-1,]) #n_fleets x n_ages
        F_fleet_ym1 <- apply(rbind(FAA_ym1),1,max) #Full F for each fleet in previous year
        om$input$par$F_pars[year,] <- log(F_fleet_y) - log(F_fleet_ym1) #change the F_dev to produce the right full F
        if(year< NROW(om$input$par$F_pars)){ #change F devs in later years to retain F in those years. Not really necessary for closed loop sims
          FAA_yp1 <- rbind(rep$FAA[,year+1,]) #n_fleets x n_ages
          F_fleet_yp1 <- apply(rbind(FAA_yp1),1,max) #Full F for each fleet in previous year
          om$input$par$F_pars[year+1,] <- log(F_fleet_yp1) - log(F_fleet_y) #change the F_dev to produce the right full F
        }
      } else om$input$par$F_pars[year,] <- log(F_fleet_y) #if year is the first year of the model, change F in year 1
    } else{ #alternative configuration of F_pars
      # cat("\nFleet-specific fishing mortality is ", Fsolve, "\n")
      om$input$par$F_pars[year,] <- log(Fsolve)
    }
  } else {
    rep = om$rep
    FAA <- rbind(rep$FAA[,year,]) #n_fleets x n_ages
    age_ind <- om$env$data$which_F_age[year] #which age is used by wham to define total "full F"
    old_max_F <- apply(FAA,2,sum)[age_ind] # n_ages
    selAA <- FAA/old_max_F #sum(selAA[i,]) = 1
    new_FAA <- Fsolve * selAA #updated FAA
    F_fleet_y <- apply(new_FAA, 1, max) #full F for each fleet
    # cat("\nFleet-specific fishing mortality is ", F_fleet_y, "\n")
    if(om$input$data$F_config==1) {
      if(year>1) {
        FAA_ym1 <- rbind(rep$FAA[,year-1,]) #n_fleets x n_ages
        F_fleet_ym1 <- apply(rbind(FAA_ym1),1,max) #Full F for each fleet in previous year
        om$input$par$F_pars[year,] <- log(F_fleet_y) - log(F_fleet_ym1) #change the F_dev to produce the right full F
        if(year< NROW(om$input$par$F_pars)){ #change F devs in later years to retain F in those years. Not really necessary for closed loop sims
          FAA_yp1 <- rbind(rep$FAA[,year+1,]) #n_fleets x n_ages
          F_fleet_yp1 <- apply(rbind(FAA_yp1),1,max) #Full F for each fleet in previous year
          om$input$par$F_pars[year+1,] <- log(F_fleet_yp1) - log(F_fleet_y) #change the F_dev to produce the right full F
        }
      } else om$input$par$F_pars[year,] <- log(F_fleet_y) #if year is the first year of the model, change F in year 1
    } else{ #alternative configuration of F_pars
      om$input$par$F_pars[year,] <- log(F_fleet_y)
    }
  }
  
  om <- fit_wham(om$input, do.fit = FALSE, do.brps = FALSE, MakeADFun.silent = TRUE)
  return(om)
}
  