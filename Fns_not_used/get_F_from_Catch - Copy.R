library(numDeriv)
library(expm)

#' Calculate fleet-specific F from fleet-specific catch
#'
#' @param Catch Vector of catch values for each fleet.
#' @param NAA Array (n_stocks x n_regions x ny x n_ages); numbers-at-age.
#' @param log_M Array (n_stocks x n_regions x ny x n_ages); natural mortality.
#' @param mu Array (n_stocks x n_ages x n_seasons x n_years x n_regions x n_regions); movement rates.
#' @param L Matrix (n_years x n_regions); "extra" mortality rate.
#' @param sel Matrix of selectivity values (n_fleets x n_ages).
#' @param fracyr_season Vector of length n_seasons; length of intervals for each season.
#' @param fleet_regions Vector of length n_fleets; indicates which region each fleet operates in.
#' @param fleet_seasons Matrix (n_fleets x n_seasons); 0/1 indicating fleet operation per season.
#' @param can_move Array (n_stocks x n_seasons x n_regions x n_regions); 0/1 for movement between regions.
#' @param mig_type Vector of length n_stocks; 0 = migration after survival, 1 = movement and mortality simultaneous.
#' @param waacatch Matrix (n_fleets x n_ages); weight-at-age.
#' @param F_init Vector of initial F values for each fleet.
#'
#' @return A vector of fleet-specific F values estimated using Newton's method.
#'
#' @export
#' 
#' 

year = 20
Catch = om$input$data$agg_catch[year,]

get_F_from_Catch(om = om, Catch = Catch, year = 20, method = "nlminb", by_fleet = 0)

get_F_from_Catch <- function(om, Catch, year, method = "nlminb", by_fleet = 0) {
  
  rep = om$rep
  NAA = rep$NAA[,,year,]
  log_M = log(rep$MAA[,,year,])
  mu = rep$mu[,,,year,,]
  L = rep$L[year,]
  sel = rbind(rep$FAA[,year,]/max(exp(rep$log_FAA_tot[year,])))
  # sel = get_avg_fleet_sel_as_array(rep$FAA, avg_years_ind, which_F_age, by_fleet = 1) 
  # sel = rep$sel_static
  fracyr_season = om$input$data$fracyr_seasons
  fleet_regions = om$input$data$fleet_regions
  fleet_seasons = om$input$data$fleet_seasons
  can_move = om$input$data$can_move
  mig_type = om$input$data$mig_type
  waacatch = om$input$data$waa[(1:om$input$data$n_fleets + 1), year,]
  
  n_iter <- 10
  n_fleets <- length(Catch)
  
  F_init <- if(by_fleet == 1) rep(0.1,n_fleets) else 0.1
  
  if (by_fleet == 1) {
    log_F_init <- log(F_init)
    log_F_fleet_newton <- matrix(log_F_init, nrow = n_fleets, ncol = n_iter)
    log_F_fleet_nlminb <- log_F_init
    log_F_fleet_optim <- log_F_init
  }
  
  if (by_fleet == 0) {
    log_F_init <- log(F_init)
    log_F_fleet_nlminb <- log_F_init
    log_F_fleet_optim <- log_F_init
  }
  
  obj_function <- function(log_F) {
    log_pred_Catch <- log_catch_fleets_F_multi(om, year, log_F, NAA, log_M, mu, L, sel,
                                               fracyr_season, fleet_regions, fleet_seasons,
                                               can_move, mig_type, waacatch, by_fleet)
    res <- sum((log_pred_Catch - log(Catch))^2)
    if (!is.finite(res)) res <- 1e6  # Penalize NaN results
    return(res)
  }
  
  if(is.null(method)) method = "nlminb"
  
  # ---- Method 1: Newton's Method ---- ##
  
  if (method == "Newton") {
    for (i in 1:(n_iter - 1)) {
      log_F_i <- log_F_fleet_newton[, i]
      
      # Compute Jacobian matrix
      jac_log_catch_at_F <- numDeriv::jacobian(
        func = function(log_F) log_catch_fleets_F_multi(om, year, log_F, NAA, log_M, mu, L, sel,
                                                        fracyr_season, fleet_regions, fleet_seasons,
                                                        can_move, mig_type, waacatch, by_fleet),
        x = log_F_i
      )
      
      if (qr(jac_log_catch_at_F)$rank < ncol(jac_log_catch_at_F)) {
        warning("Singular Jacobian matrix, stopping Newton's method.")
        break
      }
      
      inv_jac <- solve(jac_log_catch_at_F)
      diff <- log_catch_fleets_F_multi(om, year, log_F_i, NAA, log_M, mu, L, sel, fracyr_season,
                                       fleet_regions, fleet_seasons, can_move, mig_type,
                                       waacatch, by_fleet) - log(Catch)
      
      log_F_fleet_newton[, i + 1] <- log_F_fleet_newton[, i] - inv_jac %*% diff
    }
    F_solve <- exp(log_F_fleet_newton[, n_iter]) # F_newton
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

# Helper function to compute transition probabilities
get_P_t_base <- function(fleet_regions, can_move, mig_type, time, F, M, mu, L) {

  n_regions <- length(L)
  n_fleets <- length(fleet_regions)
  dim <- n_regions + n_fleets + 1
  P <- matrix(0, nrow = dim, ncol = dim)

  Z <- M + L
  for (f in 1:n_fleets) {
    Z[fleet_regions[f]] <- Z[fleet_regions[f]] + F[f]
  }

  if (n_regions == 1) {
    if (time < 1e-15) {
      P[1, 1] <- 1.0
    } else {
      P[1, 1] <- exp(-Z[1] * time)
      for (f in 1:n_fleets) {
        P[1, 1 + f] <- F[f] * (1 - exp(-Z[1] * time)) / Z[1]
      }
      P[1, 1 + n_fleets + 1] <- M[1] * (1 - exp(-Z[1] * time)) / Z[1]
    }
    diag(P)[-1] <- 1.0
  } else {
    A <- matrix(0, nrow = dim, ncol = dim)
    for (i in 1:n_regions) {
      A[i, dim] <- M[i] + L[i]
    }
    for (f in 1:n_fleets) {
      A[fleet_regions[f], n_regions + f] <- F[f]
    }
    for (i in 1:n_regions) {
      for (j in 1:n_regions) {
        if (i != j) A[i, j] <- mu[i, j]
      }
      A[i, i] <- -sum(A[i, ])
    }
    A <- A * time
    P <- expm(A)  # Matrix exponential for migration
  }

  return(P)
}

log_catch_fleets_F_multi <- function(om, year, log_F, NAA, log_M, mu, L, sel, fracyr_season, fleet_regions, fleet_seasons, can_move, mig_type, waacatch, by_fleet = 0) {
  
  n_stocks <- dim(log_M)[1]
  n_regions <- dim(log_M)[2]
  n_seasons <- length(fracyr_season)
  n_ages <- dim(log_M)[3]
  n_fleets <- length(fleet_regions)
  Pdim <- n_regions + n_fleets + 1
  
  # Initialize FAA_T
  FAA_T <- matrix(0, nrow = n_fleets, ncol = n_ages)
  
  if (length(log_F) == n_fleets) {
    for (f in 1:n_fleets) {
      for (a in 1:n_ages) {
        FAA_T[f, a] <- exp(log_F[f]) * sel[f, a] 
      }
    }
  } else {
    for (f in 1:n_fleets) {
      for (a in 1:n_ages) {
        FAA_T[f, a] <- exp(log_F[1]) * sel[f, a]  # Use global F
      }
    }
  }
  
  # Initialize Catch Matrix
  catch_stock_fleet <- matrix(0, nrow = n_stocks, ncol = n_fleets)
  annual_Ps <- om$rep$annual_Ps
  # I <- diag(Pdim)

  # for (s in 1:n_stocks) {
  #   for (a in 1:n_ages) {
  #     P_ya <- I
  #     for (t in 1:n_seasons) {
  #       P_ya <- P_ya %*% get_P_t_base(fleet_regions, can_move[s, t, , ], mig_type[s], fracyr_season[t], FAA_T, exp(log_M[s, , a]), mu[s, a, t, , ], L)
  #     }
  # 
  #     for (f in 1:n_fleets) {
  #       catch_stock_fleet[s, f] <- catch_stock_fleet[s, f] +
  #         NAA[s, fleet_regions[f], a] * P_ya[fleet_regions[f], n_regions + f] * waacatch[f, a]
  #     }
  #   }
  # }
  
  for (s in 1:n_stocks) {
    for (a in 1:n_ages) {
      P_ya <- annual_Ps[s,year,a,,]
      for (f in 1:n_fleets) {
        catch_stock_fleet[s, f] <- catch_stock_fleet[s, f] +
          NAA[s, fleet_regions[f], a] * P_ya[fleet_regions[f],n_regions + f] * waacatch[f, a]
      }
    }
  }
  
  # Compute Fleet-Specific or Global Catch
  if (length(log_F) == n_fleets) {
    Catch <- colSums(catch_stock_fleet)
  } else {
    if (n_fleets > 1 && by_fleet == 0) {
      Catch <- colSums(catch_stock_fleet)
    } else {
      Catch <- sum(catch_stock_fleet) 
    }
  }
  
  return(log(Catch))
}


get_avg_fleet_sel_as_array <- function(FAA, avg_years_ind, which_F_age, by_fleet = 0) {
  # Convert input dimensions
  n_fleets <- dim(FAA)[1]
  n_ages <- dim(FAA)[3]
  n_toavg <- length(avg_years_ind)
  
  # Initialize FAA_avg matrix
  FAA_avg <- matrix(0, nrow = n_fleets, ncol = n_ages)
  
  # Compute average FAA over selected years
  for (f in 1:n_fleets) {
    for (a in 1:n_ages) {
      FAA_avg[f, a] <- mean(FAA[f, avg_years_ind, a])  # Averaging across selected years
    }
  }
  
  # Compute total fully selected F
  FAA_avg_tot <- colSums(FAA_avg)
  F_full <- FAA_avg_tot[which_F_age]
  
  # Initialize selectivity array
  sel <- matrix(0, nrow = n_fleets, ncol = n_ages)
  
  # Compute selectivity based on by_fleet flag
  for (f in 1:n_fleets) {
    for (a in 1:n_ages) {
      if (by_fleet == 0) {
        sel[f, a] <- FAA_avg[f, a] / F_full
      } else {
        sel[f, a] <- FAA_avg[f, a] / FAA_avg[f, which_F_age]
      }
    }
  }
  
  return(sel)
}
