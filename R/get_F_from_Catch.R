#' Estimate Fleet-Specific or Global Fishing Mortality (F) from Catch
#'
#' This function estimates fleet-specific or global fishing mortality (F) using Newton's method.
#'
#' @param om A fitted operating model containing stock and fleet information.
#' @param Catch Numeric vector of catch values for each fleet.
#' @param year Integer. The year for which F is estimated.
#' @param method Character. The optimization method used for solving F. Options include:
#'   \itemize{
#'     \item \code{"nlminb"} (default): Uses `nlminb` optimization.
#'     \item \code{"BFGS"}: Uses `optim` with the BFGS method.
#'   }
#' @param by_fleet Logical. If TRUE, estimates F separately for each fleet. If FALSE, estimates a single global F (default = FALSE).
#'
#' @return A numeric vector of estimated fleet-specific F values.
#'
#' @export
#'
#'
get_F_from_Catch <- function(om, Catch, year, method = "nlminb", by_fleet = FALSE) {
  
  require(numDeriv)
  require(expm)
  
  rep = om$rep
  NAA = rep$NAA[,,year,]
  log_M = log(rep$MAA[,,year,])
  mu = rep$mu[,,,year,,]
  sel = rbind(rep$FAA[,year,]/max(exp(rep$log_FAA_tot[year,])))
  fracyr_seasons = om$input$data$fracyr_seasons
  if(length(fracyr_seasons) == 1) mu = array(mu, dim = c(dim(mu)[1:2],1,dim(mu)[3:4]))
  fleet_regions = om$input$data$fleet_regions
  fleet_seasons = om$input$data$fleet_seasons
  can_move = om$input$data$can_move
  mig_type = om$input$data$mig_type
  waacatch = om$input$data$waa[(1:om$input$data$n_fleets), year,]
  
  n_fleets <- length(Catch)
  
  F_init <- 0.1
  
  if (by_fleet) {
    log_F_init <- log(rep(F_init,n_fleets))
    log_F_fleet_nlminb <- log_F_init
    log_F_fleet_optim <- log_F_init
    for(i in 1:n_fleets){
      sel[i,] = rep$FAA[i,year,]/max(rep$FAA[i,year,])
    }
  } else {
    log_F_init <- log(F_init)
    log_F_fleet_nlminb <- log_F_init
    log_F_fleet_optim <- log_F_init
  }
  
  obj_function <- function(log_F) {
    log_pred_Catch <- log_catch_fleets_F_multi(log_F, NAA, log_M, mu, sel,
                                               fracyr_seasons, fleet_regions, fleet_seasons,
                                               can_move, mig_type, waacatch)
    # cat("\nPredicted catch is ", exp(log_pred_Catch), "\n")
    # cat("\nActual catch is ", Catch, "\n")
    res <- sum((log_pred_Catch - log(Catch))^2)
    if (!is.finite(res)) res <- 1e15 # Penalize NaN results
    return(res)
  }
  
  if (method == "nlminb") {
    ## ---- Method 2: nlminb Optimization ---- ##
    opt_nlminb <- nlminb(log_F_fleet_nlminb, objective = obj_function)
    log_F_fleet_nlminb <- opt_nlminb$par
    Fsolve <- exp(log_F_fleet_nlminb)
    # cat("\nEstimated global F is ", Fsolve, "\n")
  }
  
  if (method == "BFGS") {
    ## ---- Method 3: BFGS Optimization ---- ##
    opt_optim <- optim(log_F_fleet_optim, obj_function, method = "BFGS")
    log_F_fleet_optim <- opt_optim$par
    Fsolve <- exp(log_F_fleet_optim)
    # cat("\nEstimated global F is ", Fsolve, "\n")
  }
  
  # Apply upper bound to prevent extreme values
  Fsolve[Fsolve > 10] <- 10
  
  # Return estimated F
  return(Fsolve)
  
}

get_P_t_base <- function(fleet_regions, can_move, mig_type, time, F, M, mu) {
  
  n_regions <- length(M)
  n_fleets <- length(fleet_regions)
  dim <- n_regions + n_fleets + 1
  P <- matrix(0, nrow = dim, ncol = dim)
  
  Z <- M
  for (f in seq_along(F)) {
    Z[fleet_regions[f]] <- Z[fleet_regions[f]] + F[f]
  }
  
  if (n_regions == 1) {  # Single region case (Baranov's equation)
    if (time < 1e-15) {
      P[1, 1] <- 1.0
    } else {
      P[1, 1] <- exp(-Z[1] * time)
      for (f in seq_along(F)) {
        P[1, 1 + f] <- F[f] * (1 - exp(-Z[1] * time)) / Z[1]
      }
      P[1, dim] <- 1 - sum(P[1, 1:(dim-1)])
    }
    diag(P)[-1] <- 1.0
  } else {  # Multiple regions case
    if (sum(can_move) > 0) {  # Migration scenario
      if (mig_type == 0) {  # Migration after survival
        if (time < 1e-15) {
          P[1:n_regions, 1:n_regions] <- mu
        } else {
          for (f in seq_along(F)) {
            P[fleet_regions[f], n_regions + f] <- F[f] * (1.0 - exp(-Z[fleet_regions[f]] * time)) / Z[fleet_regions[f]]
          }
          for (i in seq_len(n_regions)) {
            for (j in seq_len(n_regions)) {
              P[i, j] <- exp(-Z[i] * time) * mu[i, j]
            }
            P[i, dim] <- 1 - sum(P[i, 1:(dim-1)])
          }
        }
        diag(P[(n_regions + 1):dim, (n_regions + 1):dim]) <- 1.0
      } else if (mig_type == 1) {  # Continuous migration
        if (time < 1e-15) {
          diag(P[1:n_regions, 1:n_regions]) <- 1.0
        } else {
          A <- matrix(0, nrow = dim, ncol = dim)
          for (i in seq_len(n_regions)) A[i, dim] <- M[i] + Z[i]
          for (f in seq_along(F)) A[fleet_regions[f], n_regions + f] <- F[f]
          for (i in seq_len(n_regions)) {
            for (j in seq_len(n_regions)) {
              if (i != j) A[i, j] <- mu[i, j]
            }
          }
          diag(A[1:n_regions, 1:n_regions]) <- -rowSums(A[1:n_regions, ])
          A <- A * time
          P <- expm(A)
        }
      }
    } else {  # No migration case
      if (time < 1e-15) {
        diag(P) <- 1.0
      } else {
        for (i in seq_len(n_regions)) {
          P[i, i] <- exp(-Z[i] * time)
          P[i, dim] <- M[i]*(1-exp(-Z[i]*time))/Z[i]
        }
        for (f in seq_along(F)) {
          P[fleet_regions[f], n_regions + f] <- F[f] * (1 - exp(-Z[fleet_regions[f]] * time)) / Z[fleet_regions[f]]
        }
        diag(P[(n_regions + 1):dim, (n_regions + 1):dim]) <- 1.0
      }
    }
  }
  
  return(P)
}

get_P_t <- function(age, stock, season, fleet_regions, fleet_seasons,
                    can_move, mig_type, time, FAA, log_M, mu) {
  
  n_regions <- dim(log_M)[2]  # Extract number of regions
  
  get_F_t <- function(fleet_season, age, FAA) {
    
    # Initialize F_t as a vector of zeros with the same length as the first dimension of FAA
    F_t <- rep(0, dim(FAA)[1])
    
    # Assign values where fleet is active
    for (f in seq_along(F_t)) {
      if (fleet_season[f] == 1) {  # Only assign if fleet is operating
        F_t[f] <- FAA[f, age]
      }
    }
    
    return(F_t)
  }
  
  # Extract fishing mortality for the given season and age
  F <- get_F_t(fleet_seasons[, season], age, FAA)
  
  # Initialize mortality and movement matrices
  M <- rep(0, n_regions)
  mu_stya <- matrix(0, n_regions, n_regions)
  
  # Compute natural mortality M
  for (r in 1:n_regions) {
    M[r] <- exp(log_M[stock, r, age])
  }
  
  # Initialize can_move and movement rate matrices
  can_move_sta <- matrix(0, n_regions, n_regions)
  for (r in 1:n_regions) {
    for (rr in 1:n_regions) {
      can_move_sta[r, rr] <- can_move[stock, season, r, rr]
      mu_stya[r, rr] <- mu[stock, age, season, r, rr]
    }
  }
  
  # Compute transition probability matrix
  P <- get_P_t_base(fleet_regions, can_move_sta, mig_type[stock], time, F, M, mu_stya)
  
  return(P)
}

get_annual_Ps <- function(fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons,
                          FAA, log_M, mu) {
  
  # Extract dimensions
  n_fleets <- dim(FAA)[1]
  n_seasons <- ncol(fleet_seasons)
  n_stocks <- dim(log_M)[1]
  n_regions <- dim(log_M)[2]
  n_ages <- dim(log_M)[3]
  P_dim <- n_regions + n_fleets + 1  # Probability transition matrix dimensions
  
  # Initialize annual probability transition array
  annual_Ps <- array(0, dim = c(n_stocks, n_ages, P_dim, P_dim))
  
  # Identity matrix for initialization
  I_mat <- diag(1, P_dim, P_dim)
  
  for (s in 1:n_stocks) {
    for (a in 1:n_ages) {
      
      P_y <- I_mat  # Reset for each year, age, stock
      
      for (t in 1:n_seasons) {
        # Compute probability transition matrix for season t
        P_t <- get_P_t(a, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons[t], FAA, log_M, mu)
        
        # Multiply season matrices together to get the annual transition matrix
        P_y <- P_y %*% P_t
        
      }
      
      # Store the annual transition matrix
      annual_Ps[s, a, , ] <- P_y
    }
  }
  
  return(annual_Ps)
}

log_catch_fleets_F_multi <- function(log_F, NAA, log_M, mu, sel, fracyr_seasons, fleet_regions, 
                                     fleet_seasons, can_move, mig_type, waacatch) {
  
  n_stocks <- dim(log_M)[1]
  n_regions <- dim(log_M)[2]
  n_seasons <- length(fracyr_seasons)
  n_ages <- dim(log_M)[3]
  n_fleets <- length(fleet_regions)
  
  # Initialize FAA_T
  FAA_T <- matrix(0, nrow = n_fleets, ncol = n_ages)
  
  if (length(log_F) == n_fleets) {
    for (f in 1:n_fleets) {
      for (a in 1:n_ages) {
        FAA_T[f, a] <- exp(log_F[f]) * sel[f, a] # use fleet-specific F
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
  
  Pdim <- n_regions + n_fleets + 1
  I <- diag(Pdim)
  
  for (s in 1:n_stocks) {
    for (a in 1:n_ages) {
      P_ya <- I
      for (t in 1:n_seasons) {
        P_ya <- P_ya %*% get_P_t(a, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_seasons[t], FAA_T, log_M, mu)
      }
      for (r in 1:n_regions) {
        for (f in 1:n_fleets) {
          catch_stock_fleet[s, f] <- catch_stock_fleet[s, f] + NAA[s, r, a] * P_ya[r, n_regions + f] * waacatch[f, a]
        }
      }
    }
  }
  
  # Should be same if directly using annual_Ps saved in the OM
  # for (s in 1:n_stocks) {
  #   for (a in 1:n_ages) {
  #     for (r in 1:n_regions) {
  #       for (f in 1:n_fleets) {
  #         catch_stock_fleet[s, f] <- catch_stock_fleet[s, f] + NAA[s, r, a] * annual_Ps[s, a, r, n_regions + f] * waacatch[f, a]
  #       }
  #     }
  #   }
  # }
  
  # Compute Fleet-Specific or Global Catch
  if (n_fleets > 1) {
    Catch <- colSums(catch_stock_fleet)
  } else {
    Catch <- sum(catch_stock_fleet) 
  }
  return(log(Catch))
}
