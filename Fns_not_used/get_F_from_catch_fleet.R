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
#' @param trace Logical; whether to print debug info.
#' @param F_init Vector of initial F values for each fleet.
#'
#' @return A vector of fleet-specific F values estimated using Newton's method.
#'
#' @export
# get_F_from_Catch_fleet <- function(Catch, NAA, log_M, mu, L, sel, fracyr_season, 
#                                    fleet_regions, fleet_seasons, can_move, mig_type, 
#                                    waacatch, trace = FALSE, F_init) {
#   n_iter <- 10
#   n_fleets <- length(Catch)
#   
#   log_F_fleet <- matrix(log(F_init), nrow = n_fleets, ncol = n_iter)
#   
#   for (i in 1:(n_iter - 1)) {
#     log_F_i <- log_F_fleet[, i]
#     
#     # Compute Jacobian matrix
#     jac_log_catch_at_F <- numDeriv::jacobian(
#       func = function(log_F) log_catch_fleets_F_multi(log_F, NAA, log_M, mu, L, sel, 
#                                                       fracyr_season, fleet_regions, fleet_seasons, 
#                                                       can_move, mig_type, waacatch, trace), 
#       x = log_F_i
#     )
#     
#     # Check for singular Jacobian
#     if (det(jac_log_catch_at_F) == 0) {
#       warning("Singular Jacobian matrix, returning last valid estimate.")
#       break
#     }
#     
#     # Compute Newton update
#     inv_jac <- solve(jac_log_catch_at_F)
#     diff <- log_catch_fleets_F_multi(log_F_i, NAA, log_M, mu, L, sel, fracyr_season, 
#                                      fleet_regions, fleet_seasons, can_move, mig_type, 
#                                      waacatch, trace) - log(Catch)
#     
#     log_F_fleet[, i + 1] <- log_F_fleet[, i] - inv_jac %*% diff
#   }
#   
#   F_iter <- exp(log_F_fleet[, n_iter])
#   F_iter[F_iter > 10] <- 10  # Cap maximum F
#   
#   return(F_iter)  # Fleet-specific F
# }

get_F_from_catch_fleet <- function(Catch, NAA, log_M, mu, L, sel, fracyr_season, fleet_regions, fleet_seasons, can_move, mig_type, waacatch, trace, F_init) {
  # n_fleets <- length(Catch)
  # log_F_init <- log(F_init)
  # 
  # # Define the objective function to minimize
  # objective_function <- function(log_F) {
  #   log_pred_Catch <- log_catch_fleets_F_multi(log_F, NAA, log_M, mu, L, sel, fracyr_season, fleet_regions, fleet_seasons, can_move, mig_type, waacatch, trace)
  #   return(sum((log_pred_Catch - log(Catch))^2))  # Minimize squared log-difference
  # }
  # 
  # # Run nlminb optimization for F in each fleet
  # opt_results <- nlminb(start = log_F_init, objective_function)
  
  n <- 10
  n_fleets <- length(Catch)
  log_F_init <- log(F_init)
  log_F_fleet <- matrix(NA, nrow = n_fleets, ncol = n)
  log_F_fleet[,1] <- log_F_init
  
  obj <- function(log_F, log_catch = log(Catch), trace = FALSE){
    if(trace) print(log_catch)
    log_pred_catch <- log_catch_fleets_F_multi(log_F, NAA, log_M, mu, L, sel, fracyr_season, fleet_regions, fleet_seasons, can_move, mig_type, waacatch, trace)
    if(trace) print(log_pred_catch)
    return(sum((log_pred_catch - log_catch)^2))
  }
  
  for (i in 2:n) {
    opt <- nlminb(log_F_fleet[,i-1], objective = function(log_F) obj(log_F, trace = FALSE))
    log_F_fleet[,i] <- opt$par
  }

  # for (i in 2:n) {
  #   opt <- optim(log_F_fleet[,i-1], obj, method = "BFGS")
  #   log_F_fleet[,i] <- opt$par
  # }
  
  F_iter <- exp(log_F_fleet)
  
  # # Extract optimized F values
  # F_opt <- exp(opt_results$par)
  
  # Apply upper bound to prevent extreme values
  F_iter[F_iter > 10] <- 10
  
  return(F_iter)
}


get_F_from_Catch_region <- function(Catch, NAA, log_M, mu, L, sel, fracyr_season, fleet_regions, fleet_seasons, can_move, mig_type, waacatch, trace, F_init) {
  n_iter <- 10
  n_fleets <- length(Catch)
  log_F_region <- matrix(log(F_init), nrow = n_fleets, ncol = n_iter)
  
  for (i in 1:(n_iter - 1)) {
    log_F_i <- log_F_region[, i]
    
    grad_log_catch_at_F <- numDeriv::grad(
      func = function(log_F) log_catch_fleets_F_multi(log_F, NAA, log_M, mu, L, sel, fracyr_season, fleet_regions, fleet_seasons, can_move, mig_type, waacatch, trace), 
      x = log_F_i
    )
    log_F_region[, i + 1] <- log_F_region[, i] - (log_catch_fleets_F_multi(log_F_i, NAA, log_M, mu, L, sel, fracyr_season, fleet_regions, fleet_seasons, can_move, mig_type, waacatch, trace) - log(Catch)) / grad_log_catch_at_F
  }
  
  F_iter <- exp(log_F_region[, n_iter])
  F_iter[F_iter > 10] <- 10
  
  return(F_iter)
}

# Helper function to compute transition probabilities
get_P_t_base <- function(fleet_regions, can_move, mig_type, time, F, M, mu, L, trace = FALSE) {
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

# log_catch_fleets_F_multi <- function(log_F, NAA, log_M, mu, L, sel, fracyr_season, fleet_regions, fleet_seasons, can_move, mig_type, waacatch, trace = FALSE) {
#   n_stocks <- dim(log_M)[1]
#   n_regions <- dim(log_M)[2]
#   n_seasons <- length(fracyr_season)
#   n_ages <- dim(log_M)[3]
#   n_fleets <- length(fleet_regions)
#   Pdim <- n_regions + n_fleets + 1
#   
#   FAA_T <- matrix(0, nrow = n_fleets, ncol = n_ages)
#   for (f in 1:n_fleets) {
#     for (a in 1:n_ages) {
#       FAA_T[f, a] <- exp(log_F[f]) * sel[f, a]
#     }
#   }
#   
#   logM_T <- log_M
#   mu_T <- mu
#   
#   catch_stock_fleet <- matrix(0, nrow = n_stocks, ncol = n_fleets)
#   I <- diag(Pdim)
#   
#   for (s in 1:n_stocks) {
#     for (a in 1:n_ages) {
#       P_ya <- I
#       for (t in 1:n_seasons) {
#         P_ya <- P_ya %*% get_P_t_base(fleet_regions, can_move[s, t, , ], mig_type[s], fracyr_season[t], FAA_T, exp(logM_T[s, , a]), mu_T[s, a, t, , ], L, trace)
#       }
#       
#       for (r in 1:n_regions) {
#         for (f in 1:n_fleets) {
#           catch_stock_fleet[s, f] <- catch_stock_fleet[s, f] + NAA[s, r, a] * P_ya[r, n_regions + f] * waacatch[f, a]
#         }
#       }
#     }
#   }
#   
#   Catch <- colSums(catch_stock_fleet)  # Ensure fleet-specific catch
#   return(log(Catch))
# }

log_catch_fleets_F_multi <- function(log_F, NAA, log_M, mu, L, sel, fracyr_season, fleet_regions, fleet_seasons, can_move, mig_type, waacatch, trace = FALSE) {
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
  I <- diag(Pdim)
  
  for (s in 1:n_stocks) {
    for (a in 1:n_ages) {
      P_ya <- I
      for (t in 1:n_seasons) {
        P_ya <- P_ya %*% get_P_t_base(fleet_regions, can_move[s, t, , ], mig_type[s], fracyr_season[t], FAA_T, exp(log_M[s, , a]), mu[s, a, t, , ], L, trace)
      }
      
      for (f in 1:n_fleets) {
        catch_stock_fleet[s, f] <- catch_stock_fleet[s, f] +
          NAA[s, fleet_regions[f], a] * P_ya[fleet_regions[f], n_regions + f] * waacatch[f, a]
      }
    }
  }
  
  # Compute Fleet-Specific or Global Catch
  if (length(log_F) == n_fleets) {
    Catch <- colSums(catch_stock_fleet)
  } else {
    Catch <- sum(catch_stock_fleet) 
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
