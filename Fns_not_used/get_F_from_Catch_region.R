#' Calculate region-specific F from region-specific catch
#' 
#' The get_F_from_Catch_region() function calculates region-specific F given the aggregate region-specific catch.
#' Note that some helper functions that support the get_F_from_Catch_region() function are included below.
#' 
#' @param Catch Vector of catch values for each region.
#' @param NAA Array (n_stocks x n_regions x ny x n_ages); numbers-at-age.
#' @param log_M Array (n_stocks x n_regions x ny x n_ages); natural mortality.
#' @param mu Array (n_stocks x n_ages x n_seasons x n_years x n_regions x n_regions); movement rates.
#' @param L Matrix (n_years x n_regions); "extra" mortality rate.
#' @param sel Matrix of selectivity values.
#' @param fracyr_season Vector of length n_seasons; length of intervals for each season.
#' @param fleet_regions Vector of length n_fleets; indicates which region each fleet is operating in.
#' @param fleet_seasons Matrix (n_fleets x n_seasons); 0/1 indicating whether fleet is operating in the season.
#' @param can_move Array (n_stocks x n_seasons x n_regions x n_regions); 0/1 indicating whether movement can occur from one region to another.
#' @param mig_type Vector of length n_stocks; 0 = migration after survival, 1 = movement and mortality simultaneous.
#' @param waacatch Matrix (n_fleets x n_ages); weight-at-age.
#' @param trace Logical; whether to print information to the screen.
#' @param F_init Vector of initial F values for each region.
#' 
#' @return A vector of region-specific F values estimated using Newton's method.
#'
#' @export
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

# Helper Functions
get_P_t_base <- function(fleet_regions, can_move, mig_type, time, F, M, mu, L, trace = 0) {
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
    if (sum(can_move) > 0) {
      if (mig_type == 0) {
        if (time < 1e-15) {
          P[1:n_regions, 1:n_regions] <- mu
        } else {
          for (f in 1:n_fleets) {
            P[fleet_regions[f], n_regions + f] <- F[f] * (1.0 - exp(-Z[fleet_regions[f]] * time)) / Z[fleet_regions[f]]
          }
          for (i in 1:n_regions) {
            for (j in 1:n_regions) {
              P[i, j] <- exp(-Z[i] * time) * mu[i, j]
            }
            P[i, dim] <- (M[i] + L[i]) * (1.0 - exp(-Z[i] * time)) / Z[i]
          }
        }
        diag(P)[(n_regions + 1):dim] <- 1.0
      } else {
        if (time < 1e-15) {
          diag(P)[1:n_regions] <- 1.0
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
          P <- expm::expm(A)
        }
      }
    } else {
      if (time < 1e-15) {
        diag(P) <- 1.0
      } else {
        for (i in 1:n_regions) {
          P[i, i] <- exp(-Z[i] * time)
          P[i, dim] <- (M[i] + L[i]) * (1.0 - exp(-Z[i] * time)) / Z[i]
        }
        for (f in 1:n_fleets) {
          P[fleet_regions[f], n_regions + f] <- F[f] * (1.0 - exp(-Z[fleet_regions[f]] * time)) / Z[fleet_regions[f]]
        }
        diag(P)[(n_regions + 1):dim] <- 1.0
      }
    }
  }
  return(P)
}

get_P_t <- function(age, stock, season, fleet_regions, fleet_seasons, can_move, mig_type, time, FAA, log_M, mu, L, trace = 0) {
  n_regions <- dim(log_M)[2]
  M <- exp(log_M[stock, , age])
  F <- FAA[, age] * as.vector(fleet_seasons[, season])
  mu_stya <- mu[stock, age, season, , ]
  can_move_sta <- can_move[stock, season, , ]
  P <- get_P_t_base(fleet_regions, can_move_sta, mig_type[stock], time, F, M, mu_stya, L, trace)
  return(P)
}

get_S <- function(P, n_regions) {
  S <- P[1:n_regions, 1:n_regions]
  return(S)
}

get_D <- function(P, n_regions, n_fleets) {
  D <- P[1:n_regions, (n_regions + 1):(n_regions + n_fleets)]
  return(D)
}

log_catch_fleets_F_multi <- function(log_F, NAA, log_M, mu, L, sel, fracyr_season, fleet_regions, fleet_seasons, can_move, mig_type, waacatch, trace) {
  n_stocks <- dim(log_M)[1]
  n_regions <- dim(log_M)[2]
  n_seasons <- length(fracyr_season)
  n_ages <- dim(log_M)[3]
  n_fleets <- dim(waacatch)[1]
  Pdim <- n_regions + n_fleets + 1
  
  FAA_T <- matrix(0, nrow = n_fleets, ncol = n_ages)
  for (f in 1:n_fleets) {
    for (a in 1:n_ages) {
      FAA_T[f, a] <- exp(log_F[fleet_regions[f]]) * sel[f, a]
    }
  }
  
  logM_T <- log_M
  mu_T <- mu
  
  catch_stock_fleet <- matrix(0, nrow = n_stocks, ncol = n_fleets)
  I <- diag(Pdim)
  
  for (s in 1:n_stocks) {
    for (a in 1:n_ages) {
      P_ya <- I
      for (t in 1:n_seasons) {
        P_ya <- P_ya %*% get_P_t(a, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_season[t], FAA_T, logM_T, mu_T, L, trace)
      }
      for (r in 1:n_regions) {
        for (f in 1:n_fleets) {
          catch_stock_fleet[s, f] <- catch_stock_fleet[s, f] + NAA[s, r, a] * P_ya[r, n_regions + f] * waacatch[f, a]
        }
      }
    }
  }
  
  Catch <- colSums(catch_stock_fleet)
  return(log(Catch))
}

get_seasonal_Ps <- function(log_F, NAA, log_M, mu, L, sel, fracyr_season, fleet_regions, fleet_seasons, can_move, mig_type, waacatch, trace) {
  n_stocks <- dim(log_M)[1]
  n_regions <- dim(log_M)[2]
  n_seasons <- length(fracyr_season)
  n_ages <- dim(log_M)[3]
  n_fleets <- dim(waacatch)[1]
  Pdim <- n_regions + n_fleets + 1
  
  FAA_T <- matrix(0, nrow = n_fleets, ncol = n_ages)
  for (f in 1:n_fleets) {
    for (a in 1:n_ages) {
      FAA_T[f, a] <- exp(log_F[fleet_regions[f]]) * sel[f, a]
    }
  }
  
  logM_T <- log_M
  mu_T <- mu
  
  seasonal_Ps <- array(NA, dim = c(n_stocks, n_seasons, n_ages, Pdim, Pdim))
  
  for (s in 1:n_stocks) {
    for (a in 1:n_ages) {
      for (t in 1:n_seasons) {
        seasonal_Ps[s, t, a, , ] <- get_P_t(a, s, t, fleet_regions, fleet_seasons, can_move, mig_type, fracyr_season[t], FAA_T, logM_T, mu_T, L, trace)
      }
    }
  }
  return(seasonal_Ps)
}
