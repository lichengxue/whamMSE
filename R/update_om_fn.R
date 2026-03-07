#' Update the Operating Model and Generate Simulated Data
#'
#' This function updates fishing mortality (F) in the operating model, projects the
#' population forward, and generates new simulated data based on the updated F time
#' series. It can either:
#'   1. behave like the original function and replace all simulated random effects, or
#'   2. keep the historical process fixed and only overwrite the future portion of
#'      `log_NAA` when `process_fix = 1`.
#'
#' @param om A fitted or unfitted operating model that includes both burn-in and feedback years.
#' @param interval.info A list containing projected catch advice for future years from the
#'   estimation model. It includes:
#'   \describe{
#'     \item{\code{years}}{Vector of projection years.}
#'     \item{\code{catch}}{Matrix (n_years x n_fleets or compatible) of projected catch.}
#'   }
#' @param seed Integer. Seed used for generating simulated data to ensure reproducibility
#'   (default = 123).
#' @param random Character vector of process names treated as random effects in the
#'   operating model (default = `"log_NAA"`).
#' @param method Character. Optimization method used for solving F.
#'   Options include `"nlminb"` (default) or `"BFGS"`.
#' @param by_fleet Logical. If TRUE, estimates F separately for each fleet.
#'   If FALSE, estimates a single global F (default = TRUE).
#' @param do.brps Logical. If TRUE, calculates reference points in the operating model
#'   (default = FALSE).
#' @param process_fix Integer or logical. If 1/TRUE, keep historical process fixed and only
#'   overwrite future `log_NAA`. If 0/FALSE, overwrite all simulated random effects as in
#'   the original function. Default = 0.
#' @param first_free_year Integer. First model year that is free to change when
#'   `process_fix = 1`. For example, `first_free_year = 51` means years 1:50 are fixed and
#'   only years 51:end are overwritten in `log_NAA`. Default = 51L.
#'
#' @return Updated operating model object.
#' @export
#'
#' @seealso \code{\link{get_F_from_Catch}}, \code{\link{update_om_F}}
update_om_fn <- function(om,
                         interval.info = NULL,
                         seed = 123,
                         random = "log_NAA",
                         method = "nlminb",
                         by_fleet = TRUE,
                         do.brps = FALSE,
                         process_fix = 0,
                         first_free_year = 1L) {
  
  obs_names <- c("agg_indices", "agg_catch", "catch_paa", "index_paa", "Ecov_obs", "obsvec")
  idx_free_start <- first_free_year - 1L
  
  update_future_only_4d <- function(target, sim, idx_free_start) {
    d <- dim(target)
    if (length(d) != 4) stop("Expected a 4D array (e.g., log_NAA).")
    if (!all(dim(sim) == d)) {
      stop(
        "Dimension mismatch between target and simulated array.\n",
        "target dim = ", paste(d, collapse = " x "),
        " ; sim dim = ", paste(dim(sim), collapse = " x ")
      )
    }
    idx_future <- idx_free_start:d[3]
    target[, , idx_future, ] <- sim[, , idx_future, ]
    target
  }
  
  update_random_effects <- function(om, om_sim, random, process_fix, idx_free_start) {
    if (!isTRUE(process_fix == 1)) {
      # original behavior: overwrite all requested simulated random effects
      if (length(random) == 1L) {
        om$input$par[[random]] <- om_sim[[random]]
        if (!is.null(om$parList) && !is.null(om$parList[[random]])) {
          om$parList[[random]] <- om_sim[[random]]
        }
      } else {
        for (nm in random) {
          om$input$par[[nm]] <- om_sim[[nm]]
          if (!is.null(om$parList) && !is.null(om$parList[[nm]])) {
            om$parList[[nm]] <- om_sim[[nm]]
          }
        }
      }
    } else {
      # process_fix = 1: only overwrite future log_NAA
      if (!identical(random, "log_NAA")) {
        stop("process_fix = 1 is currently implemented only for random = 'log_NAA'.")
      }
      if (is.null(om$input$par$log_NAA)) stop("om$input$par$log_NAA is NULL.")
      if (is.null(om_sim$log_NAA)) stop("om_sim$log_NAA is NULL from simulate().")
      
      om$input$par$log_NAA <- update_future_only_4d(
        target = om$input$par$log_NAA,
        sim    = om_sim$log_NAA,
        idx_free_start = idx_free_start
      )
      
      if (!is.null(om$parList) && !is.null(om$parList$log_NAA)) {
        om$parList$log_NAA <- update_future_only_4d(
          target = om$parList$log_NAA,
          sim    = om_sim$log_NAA,
          idx_free_start = idx_free_start
        )
      }
    }
    
    om
  }
  
  simulate_and_rebuild <- function(om, seed, random, process_fix, idx_free_start, do.brps) {
    set.seed(seed)
    
    om_sim <- om$simulate(complete = TRUE)
    
    # update simulated observations
    keep_obs <- intersect(obs_names, names(om_sim))
    om$input$data[keep_obs] <- om_sim[keep_obs]
    
    # update simulated random effects / process states
    om <- update_random_effects(
      om = om,
      om_sim = om_sim,
      random = random,
      process_fix = process_fix,
      idx_free_start = idx_free_start
    )
    
    # rebuild the OM reports with the updated state
    om <- fit_wham(
      om$input,
      do.fit = FALSE,
      do.brps = do.brps,
      MakeADFun.silent = TRUE
    )
    
    om
  }
  
  if (!is.null(interval.info)) {
    t <- 0L
    
    for (y in interval.info$years) {
      year <- which(om$years == y)
      t <- t + 1L
      Catch <- interval.info$catch[t, ]
      
      cat("\nNow calculating fleet-specific F in year ", y, "\n", sep = "")
      Fsolve <- get_F_from_Catch(om, Catch, year, method = method, by_fleet = by_fleet)
      
      cat("\nFishing mortality is ", Fsolve, "\n", sep = "")
      cat("\nNow updating F in OM for year ", y, "\n", sep = "")
      
      om <- update_om_F(om, year, Fsolve)
      
      cat("\nNow simulating (process_fix = ", process_fix,
          ", first_free_year = ", first_free_year, ") for year ", y, "\n", sep = "")
      
      om <- simulate_and_rebuild(
        om = om,
        seed = seed,
        random = random,
        process_fix = process_fix,
        idx_free_start = idx_free_start,
        do.brps = do.brps
      )
    }
  } else {
    cat("\nNow simulating (process_fix = ", process_fix,
        ", first_free_year = ", first_free_year, ")\n", sep = "")
    
    om <- simulate_and_rebuild(
      om = om,
      seed = seed,
      random = random,
      process_fix = process_fix,
      idx_free_start = idx_free_start,
      do.brps = do.brps
    )
  }
  
  return(om)
}