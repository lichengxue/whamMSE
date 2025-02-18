log_catch_fleets_F_multi <- function(log_F, NAA, log_M, mu, L, sel, fracyr_season, fleet_regions, fleet_seasons, can_move, mig_type, waacatch, annual_Ps) {
  
  n_stocks <- dim(log_M)[1]
  n_regions <- dim(log_M)[2]
  n_seasons <- length(fracyr_season)
  n_ages <- dim(log_M)[3]
  n_fleets <- length(fleet_regions)

  # Initialize FAA_T
  FAA_T <- matrix(0, nrow = n_fleets, ncol = n_ages)
  
  for (f in 1:n_fleets) {
    for (a in 1:n_ages) {
      FAA_T[f, a] <- exp(log_F) * sel[f, a]  # Use global F
    }
  }
  
  # Initialize Catch Matrix
  catch_stock_fleet <- matrix(0, nrow = n_stocks, ncol = n_fleets)
  
  for (s in 1:n_stocks) {
    for (a in 1:n_ages) {
      P_ya <- annual_Ps[s,a,,]
      for (f in 1:n_fleets) {
        catch_stock_fleet[s, f] <- catch_stock_fleet[s, f] +
          NAA[s, fleet_regions[f], a] * P_ya[fleet_regions[f],n_regions + f] * waacatch[f, a]
      }
    }
  }
  
  # Compute Fleet-Specific or Global Catch
  if (n_fleets > 1) {
    Catch <- colSums(catch_stock_fleet)
  } else {
    Catch <- sum(catch_stock_fleet) 
  }
  
  return(log(Catch))
}


year = 30
Catch = om$input$data$agg_catch[year,]

get_F_from_Catch(Catch = Catch, year = 30, method = "nlminb", om = om)

get_F_from_Catch <- function(Catch, year, method = "nlminb", om) {
  
  rep = om$rep
  NAA = rep$NAA[,,year,]
  log_M = log(rep$MAA[,,year,])
  mu = rep$mu[,,,year,,]
  L = rep$L[year,]
  sel = rbind(rep$FAA[,year,]/max(exp(rep$log_FAA_tot[year,])))
  fracyr_season = om$input$data$fracyr_seasons
  fleet_regions = om$input$data$fleet_regions
  fleet_seasons = om$input$data$fleet_seasons
  can_move = om$input$data$can_move
  mig_type = om$input$data$mig_type
  waacatch = om$input$data$waa[(1:om$input$data$n_fleets + 1), year,]
  annual_Ps <- rep$annual_Ps[,year,,,]
      
  n_fleets <- length(Catch)
  
  F_init <- 0.1
  
  log_F_init <- log(F_init)
  log_F_fleet_nlminb <- log_F_init
  log_F_fleet_optim <- log_F_init
  
  obj_function <- function(log_F) {
    log_pred_Catch <- log_catch_fleets_F_multi(log_F, NAA, log_M, mu, L, sel,
                                               fracyr_season, fleet_regions, fleet_seasons,
                                               can_move, mig_type, waacatch, annual_Ps)
    res <- sum((log_pred_Catch - log(Catch))^2)
    if (!is.finite(res)) res <- 1e6  # Penalize NaN results
    return(res)
  }
  
  if (method == "nlminb") {
    ## ---- Method 2: nlminb Optimization ---- ##
    opt_nlminb <- nlminb(log_F_fleet_nlminb, objective = obj_function)
    log_F_fleet_nlminb <- opt_nlminb$par
    Fsolve <- exp(log_F_fleet_nlminb)
  }
  
  if (method == "BFGS") {
    ## ---- Method 3: BFGS Optimization ---- ##
    opt_optim <- optim(log_F_fleet_optim, obj_function, method = "BFGS")
    log_F_fleet_optim <- opt_optim$par
    Fsolve <- exp(log_F_fleet_optim)
  }
  
  # Apply upper bound to prevent extreme values
  Fsolve[Fsolve > 10] <- 10
  
  cat("\nEstimated F is \n", Fsolve, "\n")
  
  # Return estimated F
  return(Fsolve)
  
}


