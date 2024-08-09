#' Generate basic information for multi-wham input
#' 
#' The generate_basic_info() function is used to generate biological and fishery information for 
#' simulation-estimation and management strategy evaluation. 
#' Note that not all information is included here. More details can be found in \code{\link{prepare_wham_input}}
#' 
#' @param n_stocks Number of stocks
#' @param n_regions Number of regions 
#' @param n_indices Number of indices
#' @param n_fleets Number of fleets
#' @param n_seasons Number of seasons
#' @param base.years Model years (or "burn-in" period in MSE)
#' @param n_feedback_years (optional) Number of years in the MSE feedback loop
#' @param life_history Fish life history (parameters are obtained from \href{https://doi.org/10.1139/cjfas-2016-0381}{Wiedenmann et al. 2017}):
#'   \itemize{
#'     \item \code{"median"} median-lived species (default)
#'     \item \code{"short"} short-lived species
#'     \item \code{"long"} long-lived species
#'     }
#' @param n_ages Number of ages 
#' @param Fbar_ages Ages to use to average total F at age defining fully selected F for reference points.
#' @param recruit_model Option to specify stock-recruit model
#'   \itemize{
#'     \item \code{1} {SCAA (without NAA_re option specified) or Random walk (if NAA_re$sigma specified), i.e. predicted recruitment in year i = recruitment in year i-1}
#'     \item \code{2} {(default) Random about mean, i.e. steepness = 1}
#'     \item \code{3} {Beverton-Holt}
#'     \item \code{4} {Ricker}
#'   }
#' @param q Survey catchability
#' @param F_info Historical fishing pressure
#'   \itemize{
#'     \item \code{$F.year1} fishing mortality in the first year
#'     \item \code{$Fhist}: \cr
#'     {"constant"} = constant across years \cr
#'     {"updown"}   = increase from F.year1 to Fmax until a change point and then decrease to Fmin \cr
#'     {"downup"}   = decrease from F.year1 to Fmin until a change point and then increase to Fmax \cr
#'     {"F-H-L"} = constant Fmax(a multiplier)xF.year1 until the change point and then constant Fmin(a multiplier)xF.year1
#'     \item \code{$Fmax} maximum F (or a multiplier when Fhist = "F-H-L")
#'     \item \code{$Fmin} minimum F (or a multiplier when Fhist = "F-H-L")
#'     \item \code{$change_time} a ratio that determines the change point
#'     }
#' @param catch_info Fleet information
#'   \itemize{
#'     \item \code{$catch_cv} cv for the catch
#'     \item \code{$catch_Neff} effective sample size for the catch
#'     }
#' @param index_info Survey information (see details)
#'   \itemize{
#'     \item \code{$index_cv} cv for the indices
#'     \item \code{$index_Neff} effective sample size for the indices
#'     \item \code{$fracyr_indices} fraction of the year when survey is conducted
#'     }
#' @param fracyr_spawn Fraction of the year when spawning occurs
#' @param bias.correct.process T/F process error is bias corrected
#' @param bias.correct.observation T/F observation error is bias corrected
#' @param bias.correct.BRPs T/F biological reference point is bias corrected
#' @param mig_type 0 = migration after survival, 1 = movement and mortality simultaneous
#' @param XSPR_R_opt
#'   \itemize{
#'     \item \code{1} use annual R estimates for annual SSB_XSPR
#'     \item \code{2} use average R estimates for annual SSB_XSPR (default)
#'     \item \code{3} use annual expected R for annual SSB_XSPR
#'     \item \code{4} use average expected R estimates for annual SSB_XSPR
#'     \item \code{5} use bias-corrected expected recruitment
#'     }
#' @return A list of information that will be used to prepare wham input
#'
#' @export
#'
#' @seealso \code{\link{prepare_wham_input}}
#'
#' @examples
#' \dontrun{
#' data <- generate_basic_info()
#' }
generate_basic_info <- function(n_stocks = 2, 
                                n_regions = 2, 
                                n_indices = 2, 
                                n_fleets = 2, 
                                n_seasons = 4, 
                                base.years = 1993:2022,
                                n_feedback_years = 0,
                                life_history = "medium",
                                n_ages = 12, 
                                Fbar_ages = 12, 
                                recruit_model = 2, 
                                F_info = list(F.year1 = 0.2, Fhist = "constant", Fmax = 0.2, Fmin = 0.2, change_time = 0.5),
                                catch_info = list(catch_cv = 0.1,catch_Neff = 100),
                                index_info = list(index_cv = 0.1,index_Neff = 100, fracyr_indices = 0.625, q = 0.2),
                                fracyr_spawn = 0.625,
                                bias.correct.process = FALSE,
                                bias.correct.observation = FALSE,
                                bias.correct.BRPs = FALSE,
                                mig_type = 0,
                                XSPR_R_opt = 2) {
  
  check_dimensions <- function(...){
    length(unique(c(...))) == 1
  }
  
  if (!check_dimensions(n_stocks, n_regions, n_indices, n_fleets)) cat("\nn_stocks, n_regions, n_fleets, n_indices are not the same!")
  
  basic_info = list()
  basic_info$bias_correct_process = bias.correct.process
  basic_info$bias_correct_observation = bias.correct.observation
  basic_info$bias_correct_BRPs = bias.correct.BRPs
  basic_info$mig_type = mig_type # 0: mortality and movement separate,  1: mortality and movement simultaneous
  basic_info$XSPR_R_opt = XSPR_R_opt 
  
  basic_info$years = as.integer(base.years[1] - 1 + 1:(length(base.years) + n_feedback_years))
  basic_info$ages = as.integer(1:n_ages)
  basic_info$n_ages = as.integer(length(basic_info$ages))
  
  na = n_ages
  ny = length(basic_info$years)
  
  basic_info$recruit_model = recruit_model
  
  basic_info$n_stocks  = as.integer(n_stocks)
  basic_info$n_regions = as.integer(n_regions)
  basic_info$n_indices = as.integer(n_indices)
  basic_info$n_fleets  = as.integer(n_fleets)
  basic_info$n_seasons = as.integer(n_seasons)
  
  basic_info$fracyr_seasons = rep(1/n_seasons, n_seasons) # Can be others/ User can define this
  
  basic_info$fracyr_SSB <- matrix(0, ny, n_stocks) # Assume fish is recruited at the beginning of the year
  basic_info$fracyr_spawn = rep(fracyr_spawn, n_stocks)
  # basic_info$spawn_regions = 1:n_stocks
  
  F_info_list <- list2env(F_info, envir = .GlobalEnv)
  
  # Fishing mortality 
  nby = length(base.years)
  if (is.null(F_info)) {
    F = matrix(0.2, nby, n_fleets)
  } else {
    if(is.null(F_info$F.year1)) stop("Users must specify initial F!")
    F.year1 = F_info$F.year1
    if(F_info$Fhist == "constant") {
      F = matrix(F.year1, nby, n_fleets)
    } else {
      Fmax = F_info$Fmax
      if(is.null(Fmax)) stop("Fmax must be specified!")
      Fmin = F_info$Fmin
      if(is.null(Fmin)) stop("Fmin must be specified!")
      change_time = F_info$change_time
      if(is.null(change_time)) stop("change_time must be specified!")
      mid <- ceiling(nby * change_time)
      if(F_info$Fhist == "updown") {
        F = matrix(c(seq(F.year1, Fmax, length.out = mid), seq(Fmax, Fmin, length.out = nby - mid + 1)[-1]), nby, n_fleets)
      } 
      if(F_info$Fhist == "downup") {
        F = matrix(c(seq(F.year1, Fmin, length.out = mid), seq(Fmin, Fmax, length.out = nby - mid + 1)[-1]), nby, n_fleets)
      } 
      if(F_info$Fhist == "F-H-L") {
        F <- matrix(F.year1 * Fmin, nby, n_fleets)
        F[1:mid,] = F.year1 * Fmax
      }
    }
  }
  
  if(n_feedback_years > 0) {
    if(n_fleets > 1) F <- rbind(F, F[rep(nby, n_feedback_years), , drop = FALSE])
    else F <- rbind(F, matrix(F[rep(nby, n_feedback_years)], ncol = 1))
  }
  
  F_info = list(F = F, F_config = 2)
  
  if(is.null(Fbar_ages)) {
    basic_info$Fbar_ages = as.integer(na)
  } else {
    basic_info$Fbar_ages = as.integer(Fbar_ages)
  }
  
  # Maturity at age
  maturity <- Generate_Maturity(life_history, na)
  maturity <- t(matrix(maturity, na, ny))
  basic_info$maturity <- array(NA, dim = c(n_stocks, ny, na))
  for (i in 1:n_stocks) basic_info$maturity[i, , ] <- maturity
  
  # Weight at age
  W <- Generate_WAA(life_history, na)
  nwaa <- n_fleets + n_regions + n_indices + n_stocks
  basic_info$waa <- array(NA, dim = c(nwaa, ny, na))
  for(i in 1:nwaa) basic_info$waa[i, , ] <- t(matrix(W, na, ny))
  
  basic_info$waa_pointer_fleets   <- 1:n_fleets
  basic_info$waa_pointer_totcatch <- (n_fleets + 1):(n_fleets + n_regions)
  basic_info$waa_pointer_indices  <- (n_fleets + n_regions + 1):(n_fleets + n_regions + n_indices)
  basic_info$waa_pointer_ssb      <- (n_fleets + n_regions + n_indices + 1):(n_fleets + n_regions + n_indices + n_stocks)
  basic_info$waa_pointer_M        <- basic_info$waa_pointer_ssb
  
  # Catch information
  # catch_info_list <- catch_info
  catch_info_list <- list2env(catch_info, envir = .GlobalEnv)
  
  if (is.null(catch_info)) {
    catch_cv.input = 0.1
    catch_Neff.input = 100
  } else {
    catch_cv.input = catch_info$catch_cv
    catch_Neff.input = catch_info$catch_Neff
  }
  
  catch_info <- list()
  catch_info$n_fleets = n_fleets
  catch_info$agg_catch = catch_info$agg_catch_sigma = catch_info$catch_Neff = matrix(NA, ny, n_fleets)
  catch_info$use_catch_paa = catch_info$use_agg_catch = matrix(NA, ny, n_fleets)
  catch_info$selblock_pointer_fleets = matrix(NA, ny, n_fleets)
  catch_info$catch_paa = array(NA, dim = c(n_fleets, ny, na))
  catch_info$agg_catch[] = 1000
  catch_info$catch_paa[] = 1 / na
  catch_info$catch_cv = catch_info_list$catch_cv
  catch_info$agg_catch_sigma[] = sqrt(log(catch_info_list$catch_cv^2 + 1))
  catch_info$catch_Neff[] = catch_info_list$catch_Neff
  catch_info$use_catch_paa[] = 1
  catch_info$use_agg_catch[] = 1 
  catch_info$selblock_pointer_fleets = t(matrix(1:n_fleets, n_fleets, ny))
  
  if(n_regions == 1) catch_info$fleet_regions = rep(1, n_fleets)
  if(n_regions > 1) {
    catch_info$fleet_regions = unlist(lapply(1:n_regions, function(r) rep(r, n_fleets / n_regions)))
  }
  
  # Index information
  # index_info_list <- index_info
  index_info_list <- list2env(index_info, envir = .GlobalEnv)
  
  index_info <- list()
  index_info$n_indices = n_indices
  index_info$agg_indices = index_info$agg_index_sigma = index_info$index_Neff = matrix(NA, ny, n_indices)
  index_info$index_paa = array(NA, dim = c(n_indices, ny, na))
  index_info$use_indices = index_info$use_index_paa = matrix(NA, ny, n_indices)
  index_info$units_indices = index_info$units_index_paa = rep(NA, n_indices)
  index_info$index_seasons = rep(NA, n_indices)
  index_info$q <- rep(index_info_list$q, n_indices)
  
  if(n_regions == 1) index_info$index_regions = rep(1, n_indices)
  if(n_regions > 1) {
    index_info$index_regions = unlist(lapply(1:n_regions, function(r) rep(r, n_indices / n_regions)))
  }
  
  index_info$index_cv = index_info_list$index_cv
  index_info$agg_indices[] = 1000
  index_info$agg_index_sigma[] = sqrt(log(index_info_list$index_cv^2 + 1))
  index_info$index_Neff[] = index_info_list$index_Neff
  index_info$index_paa[] = 1 / na
  index_info$use_indices[] = 1
  index_info$use_index_paa[] = 1
  index_info$selblock_pointer_indices = t(matrix(n_fleets + 1:n_indices, n_indices, ny))
  index_info$units_indices = rep(2, n_indices)
  index_info$units_index_paa = rep(2, n_indices)
  
  index_info$fracyr_indices = matrix(NA, ny, n_indices)
  for (s in 1:n_stocks) index_info$fracyr_indices[, s] = index_info_list$fracyr_indices
  
  for(s in 1:n_indices) {
    int_starts = cumsum(c(0, basic_info$fracyr_seasons))
    ind <- max(which(int_starts <= index_info_list$fracyr_indices))
    index_info$index_seasons[s] <- ind
    index_info$fracyr_indices[, s] = index_info_list$fracyr_indices - int_starts[ind]
  }
  
  for(s in 1:n_stocks) {
    int_starts = cumsum(c(0, basic_info$fracyr_seasons))
    ind <- max(which(int_starts <= fracyr_spawn))
    basic_info$spawn_seasons[s] <- ind
    basic_info$fracyr_SSB[, s] = fracyr_spawn - int_starts[ind] 
  }
  
  par_inputs = list(n_stocks = n_stocks, n_regions = n_regions, 
                    n_indices = n_indices, n_fleets = n_fleets, 
                    n_seasons = n_seasons, life_history = life_history,
                    n_ages = n_ages, Fbar_ages = Fbar_ages, 
                    recruit_model = recruit_model, 
                    F.year1 = F_info_list$F.year1, Fhist = F_info_list$Fhist, Fmax = F_info_list$Fmax, 
                    Fmin = F_info_list$Fmin, change_time = F_info_list$change_time,
                    catch_cv = catch_info_list$catch_cv, catch_Neff = catch_Neff.input,
                    index_cv = index_info_list$index_cv, index_Neff = index_info_list$index_Neff, fracyr_indices = index_info_list$fracyr_indices, q = index_info_list$q,
                    fracyr_spawn = fracyr_spawn,
                    bias.correct.process = bias.correct.process,
                    bias.correct.observation = bias.correct.observation,
                    bias.correct.BRPs = bias.correct.BRPs,
                    mig_type = mig_type, XSPR_R_opt = XSPR_R_opt)
  
  return(list(basic_info = basic_info, catch_info = catch_info, index_info = index_info, F = F_info, par_inputs = par_inputs))
}

Generate_Maturity <- function(life_history = NULL, na) {
  if (is.null(life_history)){
    warning("Life history is not specified and default is used!")
    maturity <- t(matrix(1 / (1 + exp(-1 * (1:na - na / 2))), na))
  } else if (life_history == "short") {
    m50 = 1.75; mslope = 1
    maturity <- t(matrix(1 / (1 + exp(-(1:na - m50) / mslope)), na))
  } else if (life_history == "medium") {
    m50 = 3.5; mslope = 1
    maturity <- t(matrix(1 / (1 + exp(-(1:na - m50) / mslope)), na))
  } else if (life_history == "long") {
    m50 = 7; mslope = 1
    maturity <- t(matrix(1 / (1 + exp(-(1:na - m50) / mslope)), na))
  }
  return(maturity)
}

Generate_Len <- function(Linf, k, n_ages) {
  Len <- Linf * (1 - exp(-k * 1:n_ages))
  return(Len)
}

Generate_WAA <- function(life_history = NULL, na) {
  if (is.null(life_history)){
    warning("Life history is not specified and default is used!")
    Len <- 100 * (1 - exp(-0.3 * (1:na - 0)))
    W <- 3e-6 * Len^3
  } else if (life_history == "short") {
    k = 0.27; Linf = 90
    Len <- Generate_Len(Linf, k, na)
    LWexp = 3; LWscaler = 3e-6
    W <- LWscaler * Len^LWexp
  } else if (life_history == "medium") {
    k = 0.13; Linf = 90
    Len <- Generate_Len(Linf, k, na)
    LWexp = 3; LWscaler = 3e-6
    W <- LWscaler * Len^LWexp
  } else if (life_history == "long") {
    k = 0.07; Linf = 90
    Len <- Generate_Len(Linf, k, na)
    LWexp = 3; LWscaler = 3e-6
    W <- LWscaler * Len^LWexp
  }
  return(W)
}

set_sel <- function(a50, k, ages){
  selAA <- 1 / (1 + exp(-(ages - a50) / k))
  return(selAA)
}
