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
#' @param F_info Historical fishing pressure or user-specified F values
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
#'     \item \code{$user_F} (optional) {matrix (n_years x n_fleets) of user-specified fishing mortality.}
#'     }
#' @param catch_info Fleet information
#'   \itemize{
#'     \item \code{$catch_cv} {vector (n_fleets) of CVs for fleet catch.}
#'     \item \code{$catch_Neff} {vector (n_fleets) of effective sample sizes for fleet catch.}
#'     \item \code{$use_agg_catch} {vector (n_fleets) of 0/1 values flagging whether whether to use aggregate catch.}
#'     \item \code{$use_catch_paa} {vector (n_fleets) of 0/1 values flagging whether to use proportions at age observations.}
#'     }
#' @param index_info Survey information (see details)
#'   \itemize{
#'     \item \code{$index_cv} {vector (n_indices) of CVs for survey indices.}
#'     \item \code{$index_Neff} {vector (n_indices) of effective sample sizes for survey indices.}
#'     \item \code{$fracyr_indices} {vector (n_indices) of fractions of the year when each survey is conducted.}
#'     \item \code{$q} {vector (n_indices) of survey catchabilities.}
#'     \item \code{$use_indices} {vector (n_indices) of 0/1 values flagging whether to use aggregate observations.}
#'     \item \code{$use_index_paa} {vector (n_indices) of 0/1 values flagging whether to use proportions at age observations.}
#'     \item \code{$units_indices} {vector (n_indices) of 1/2 values flagging whether aggregate observations are biomass (1) or numbers (2).}
#'     \item \code{$units_index_paa} {vector (n_indices) of 1/2 values flagging whether composition observations are biomass (1) or numbers (2).}
#'     }
#' @param fracyr_spawn Fraction of the year when spawning occurs
#' @param fracyr_seasons (optional) User-defined fraction of the year for each season. Should sum to 1.
#' @param fleet_pointer (optional) User-defined vector specifying the region allocation for each fleet. Default is every region has 1 fleet.
#' @param index_pointer (optional) User-defined vector specifying the region allocation for each index. Default is every region has 1 index.
#' @param user_waa (optional) User-defined weight-at-age vector (length = n_ages) or matrix (c(n_fleets + n_regions + n_indices + n_stocks) x n_ages) to override default WAA values. 
#' @param user_maturity (optional) User-defined maturity-at-age vector to override default maturity values.
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
#' @param move_dyn Movement dynamics. 0 = natal homing, 1 = meta-population
#' @param onto_move Array of age-specific movement type with dims = n_stocks x n_regions (from region x) x n_regions-1 (to region y)
#'   \itemize{
#'     \item \code{0} same mean movement rate across ages (default)
#'     \item \code{1} increasing logistic (2 parameters)
#'     \item \code{2} decreasing logistic (2 parameters)
#'     \item \code{3} double-logistic (4 parameters)
#'     \item \code{4} double-normal (4 parameters)
#'     \item \code{5} user specified age-specific movement rate
#'     }
#' 
#' @param onto_move_pars Parameters for age-specific movement, array of dim = c(n_stocks, n_regions, 4) 
#' @param apply_re_trend Indicator to apply recruitment trend (default is 0)
#' @param trend_re_rate Rate for recruitment trend (required if apply_re_trend = 1)
#' @param apply_mu_trend Indicator to apply movement trend (default is 0)
#' @param trend_mu_rate Rate for movement trend (required if apply_mu_trend = 1)
#' @param age_mu_devs Movement deviations used when onto_move = 5 (user specified age-specific movement rate), array of dim = c(n_stocks, n_regions, n_regions-1, n_ages)
#' 
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
                                n_seasons = 5,
                                base.years = 1:20,
                                n_feedback_years = 0,
                                life_history = "medium",
                                n_ages = 12,
                                Fbar_ages = 12,
                                recruit_model = 2,
                                F_info = list(F.year1 = 0.2, Fhist = "constant", Fmax = 0.2, Fmin = 0.2, change_time = 0.5, user_F = NULL),
                                catch_info = list(catch_cv = 0.1, catch_Neff = 100, use_agg_catch = 1, use_catch_paa = 1),
                                index_info = list(index_cv = 0.1, index_Neff = 100, fracyr_indices = 0.5, q = 0.2, use_indices = 1, use_index_paa = 1, units_indices = 2, units_index_paa = 2),
                                fracyr_spawn = 0.5,
                                fracyr_seasons = NULL,
                                fleet_pointer = NULL,
                                index_pointer = NULL,
                                user_waa = NULL,
                                user_maturity = NULL,
                                bias.correct.process = FALSE,
                                bias.correct.observation = FALSE,
                                bias.correct.BRPs = FALSE,
                                mig_type = 0,
                                XSPR_R_opt = 2,
                                move_dyn = 0,
                                onto_move = 0,
                                onto_move_pars = NULL,
                                apply_re_trend = 0,
                                trend_re_rate = NULL,
                                apply_mu_trend = 0,
                                trend_mu_rate = NULL,
                                age_mu_devs = NULL) {
  
  check_dimensions <- function(...) {
    length(unique(c(...))) == 1
  }
  
  if (!check_dimensions(n_stocks, n_regions, n_indices, n_fleets)) cat("\nn_stocks, n_regions, n_fleets, n_indices are not the same!")
  
  basic_info = list()
  basic_info$bias_correct_process = bias.correct.process
  basic_info$bias_correct_observation = bias.correct.observation
  basic_info$bias_correct_BRPs = bias.correct.BRPs
  basic_info$mig_type = mig_type # 0: mortality and movement separate,  1: mortality and movement simultaneous
  basic_info$XSPR_R_opt = XSPR_R_opt
  basic_info$move_dyn = move_dyn
  basic_info$onto_move = onto_move
  basic_info$onto_move_pars = onto_move_pars
  basic_info$age_mu_devs = age_mu_devs
  basic_info$apply_re_trend = apply_re_trend
  basic_info$trend_re_rate = ifelse(apply_re_trend == 1, trend_re_rate, 0)
  basic_info$apply_mu_trend = apply_mu_trend
  basic_info$trend_mu_rate = ifelse(apply_mu_trend == 1, trend_mu_rate, 0)
  
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
  
  if (is.null(fracyr_seasons)) {
    basic_info$fracyr_seasons = rep(1 / n_seasons, n_seasons) # Default: even split across seasons
  } else {
    if (length(fracyr_seasons) != n_seasons || sum(fracyr_seasons) != 1) {
      stop("fracyr_seasons must have length equal to n_seasons and sum to 1.")
    }
    basic_info$fracyr_seasons = fracyr_seasons
  }
  
  basic_info$fracyr_SSB <- matrix(0, ny, n_stocks) # Assume fish is recruited at the beginning of the year
  basic_info$fracyr_spawn = fracyr_spawn # rep(fracyr_spawn, n_stocks) different spawning season???
  
  catch_info_list <- catch_info
  index_info_list <- index_info
  
  # Fishing mortality
  nby <- length(base.years)
  
  # Extract relevant values directly from F_info
  F.year1 <- F_info$F.year1
  Fhist <- F_info$Fhist
  Fmax <- F_info$Fmax
  Fmin <- F_info$Fmin
  change_time <- F_info$change_time
  user_F <- F_info$user_F
  
  if (!is.null(user_F)) {
    # Check user-specified F matrix
    if (!all(dim(user_F) == c(nby, n_fleets))) {
      stop("user_F must be a matrix with dimensions n_years x n_fleets.")
    }
    F <- user_F
  } else {
    if (is.null(F.year1)) stop("Users must specify initial F!")
    
    if (Fhist == "constant") {
      F <- matrix(F.year1, nby, n_fleets)
    } else {
      if (is.null(Fmax)) stop("Fmax must be specified!")
      if (is.null(Fmin)) stop("Fmin must be specified!")
      if (is.null(change_time)) stop("change_time must be specified!")
      
      mid <- ceiling(nby * change_time)
      
      if (Fhist == "updown") {
        F <- matrix(c(seq(F.year1, Fmax, length.out = mid),
                      seq(Fmax, Fmin, length.out = nby - mid + 1)[-1]), nby, n_fleets)
      } else if (Fhist == "downup") {
        F <- matrix(c(seq(F.year1, Fmin, length.out = mid),
                      seq(Fmin, Fmax, length.out = nby - mid + 1)[-1]), nby, n_fleets)
      } else if (Fhist == "F-H-L") {
        F <- matrix(F.year1 * Fmin, nby, n_fleets)
        F[1:mid, ] <- F.year1 * Fmax
      }
    }
  }
  
  if (n_feedback_years > 0) {
    if (n_fleets > 1) F <- rbind(F, F[rep(nby, n_feedback_years), , drop = FALSE])
    else F <- rbind(F, matrix(F[rep(nby, n_feedback_years)], ncol = 1))
  }
  
  F_info = list(F = F, F_config = 2)
  
  if (is.null(Fbar_ages)) {
    basic_info$Fbar_ages = as.integer(na)
  } else {
    basic_info$Fbar_ages = as.integer(Fbar_ages)
  }
  
  # Maturity at age
  if (is.null(user_maturity)) {
    maturity <- Generate_Maturity(life_history, na)
  } else {
    if (length(user_maturity) != na) {
      stop("user_maturity must have length equal to n_ages.")
    }
    maturity <- user_maturity
  }
  maturity <- t(matrix(maturity, na, ny))
  basic_info$maturity <- array(NA, dim = c(n_stocks, ny, na))
  for (i in 1:n_stocks) basic_info$maturity[i, , ] <- maturity
  
  # Weight at age
  nwaa <- n_fleets + n_regions + n_indices + n_stocks
  basic_info$waa <- array(NA, dim = c(nwaa, ny, na))
  
  if (is.null(user_waa)) {
    W <- Generate_WAA(life_history, na)
    for (i in 1:nwaa) basic_info$waa[i, , ] <- t(matrix(W, na, ny))
  } else {
    if (length(user_waa) == na) {
      W <- user_waa
      for (i in 1:nwaa) basic_info$waa[i, , ] <- do.call(rbind, replicate(ny, W, simplify = FALSE))
      } else if (is.matrix(user_waa) && dim(user_waa)[1] == nwaa && dim(user_waa)[2] == na) {
        W <- user_waa
        for (i in 1:nwaa) basic_info$waa[i, , ] <- do.call(rbind, replicate(ny, W[i,], simplify = FALSE))
      } else {
        warnings("Dimension of W should be either a vector of n_ages or a matrix with nrow = c(n_fleets + n_regions + n_indices + n_stocks) and ncol = n_ages!")
      }
  }

  basic_info$waa_pointer_fleets   <- 1:n_fleets
  basic_info$waa_pointer_totcatch <- (n_fleets + 1):(n_fleets + n_regions)
  basic_info$waa_pointer_indices  <- (n_fleets + n_regions + 1):(n_fleets + n_regions + n_indices)
  basic_info$waa_pointer_ssb      <- (n_fleets + n_regions + n_indices + 1):(n_fleets + n_regions + n_indices + n_stocks)
  basic_info$waa_pointer_M        <- basic_info$waa_pointer_ssb
  
  # Catch information
  if (is.null(catch_info)) {
    # If catch_info is not provided, set default values
    catch_cv.input <- rep(0.1, n_fleets)
    catch_Neff.input <- rep(100, n_fleets)
  } else {
    # Handle catch_cv
    if (length(catch_info$catch_cv) == 1) {
      catch_cv.input <- rep(catch_info$catch_cv, n_fleets)
    } else if (length(catch_info$catch_cv) == n_fleets) {
      catch_cv.input <- catch_info$catch_cv
    } else {
      stop("catch_cv must be either a single value or a vector of length n_fleets.")
    }
    
    # Handle catch_Neff
    if (length(catch_info$catch_Neff) == 1) {
      catch_Neff.input <- rep(catch_info$catch_Neff, n_fleets)
    } else if (length(catch_info$catch_Neff) == n_fleets) {
      catch_Neff.input <- catch_info$catch_Neff
    } else {
      stop("catch_Neff must be either a single value or a vector of length n_fleets.")
    }

  }

  if (is.null(catch_info$use_agg_catch)) {
    use_agg_catch <- matrix(1, ny, n_fleets)
  } else {
    if(length(catch_info$use_agg_catch) == n_fleets) stop("Length of use_agg_catch should be n_fleets!")
    use_agg_catch <- matrix(catch_info$use_agg_catch, ny, n_fleets, byrow = TRUE)
  }
  
  if (is.null(catch_info$use_catch_paa)) {
    use_catch_paa <- matrix(1, ny, n_fleets)
  } else {
    if(length(catch_info$use_catch_paa) == n_fleets) stop("Length of use_catch_paa should be n_fleets!")
    use_catch_paa <- matrix(catch_info$use_catch_paa, ny, n_fleets, byrow = TRUE)
  }
  
  # Create catch_info list with all necessary attributes
  if(is.null(catch_info)) catch_info <- list()
  
  catch_info$n_fleets <- n_fleets
  catch_info$use_agg_catch <- use_agg_catch
  catch_info$use_catch_paa <- use_catch_paa
  catch_info$agg_catch <- matrix(1000, ny, n_fleets)
  catch_info$agg_catch_cv <- matrix(sqrt(log(catch_cv.input^2 + 1)), ny, n_fleets, byrow = TRUE)
  catch_info$catch_Neff <- matrix(catch_Neff.input, ny, n_fleets, byrow = TRUE)
  catch_info$selblock_pointer_fleets <- t(matrix(1:n_fleets, n_fleets, ny))
  catch_info$catch_paa <- array(1 / na, dim = c(n_fleets, ny, na))
  catch_info$catch_cv <- catch_cv.input
  
  # Handle fleet_regions
  if (is.null(fleet_pointer)) {
    if (n_regions == 1) {
      catch_info$fleet_regions <- rep(1, n_fleets)
    } else if (n_regions > 1) {
      # Check if n_fleets is evenly divisible by n_regions
      if (n_fleets %% n_regions != 0) {
        warning("The number of fleets is not evenly divisible by the number of regions. Please specify 'fleet_pointer' manually.")
      } else {
        # Assign fleet_regions sequentially as 1, 1, 2, 2
        catch_info$fleet_regions <- rep(1:n_regions, each = n_fleets / n_regions)
      }
    }
  } else {
    if (length(fleet_pointer) != n_fleets) {
      stop("fleet_pointer must have length equal to n_fleets.")
    }
    catch_info$fleet_regions <- fleet_pointer
  }
  
  # Index information
  if (is.null(index_info)) {
    # Set default values if index_info is not provided
    index_cv.input <- rep(0.1, n_indices)
    index_Neff.input <- rep(100, n_indices)
    fracyr_indices.input <- rep(0.2, n_indices)  # Default value for fracyr_indices
    q.input <- rep(0.1, n_indices)               # Default value for q
  } else {
    # Handle index_cv
    if (length(index_info$index_cv) == 1) {
      index_cv.input <- rep(index_info$index_cv, n_indices)
    } else if (length(index_info$index_cv) == n_indices) {
      index_cv.input <- index_info$index_cv
    } else {
      stop("index_cv must be either a single value or a vector of length n_indices.")
    }
    
    # Handle index_Neff
    if (length(index_info$index_Neff) == 1) {
      index_Neff.input <- rep(index_info$index_Neff, n_indices)
    } else if (length(index_info$index_Neff) == n_indices) {
      index_Neff.input <- index_info$index_Neff
    } else {
      stop("index_Neff must be either a single value or a vector of length n_indices.")
    }
    
    # Handle fracyr_indices
    if (is.null(index_info$fracyr_indices)) {
      fracyr_indices.input <- rep(0.2, n_indices)  # Default value if not provided
    } else if (length(index_info$fracyr_indices) == 1) {
      fracyr_indices.input <- rep(index_info$fracyr_indices, n_indices)
    } else if (length(index_info$fracyr_indices) == n_indices) {
      fracyr_indices.input <- index_info$fracyr_indices
    } else {
      stop("fracyr_indices must be either a single value or a vector of length n_indices.")
    }
    
    # Handle q
    if (is.null(index_info$q)) {
      q.input <- rep(0.1, n_indices)  # Default value if not provided
    } else if (length(index_info$q) == 1) {
      q.input <- rep(index_info$q, n_indices)
    } else if (length(index_info$q) == n_indices) {
      q.input <- index_info$q
    } else {
      stop("q must be either a single value or a vector of length n_indices.")
    }
  }
  
  # Create the index_info list with necessary values
  if(is.null(index_info)) index_info <- list()
  index_info$n_indices <- n_indices
  index_info$agg_indices <- matrix(1000, ny, n_indices)
  index_info$index_paa <- array(1 / na, dim = c(n_indices, ny, na))
  index_info$agg_index_cv <- matrix(sqrt(log(index_cv.input^2 + 1)), ny, n_indices, byrow = TRUE)
  index_info$index_Neff <- matrix(index_Neff.input, ny, n_indices, byrow = TRUE)
  index_info$q <- q.input
  
  if (is.null(index_info$use_indices)) {
    use_indices <- matrix(1, ny, n_indices)
  } else {
    if(length(index_info$use_indices) == 1) {
      index_info$use_indices <- rep(index_info$use_indices, n_indices)
      use_indices <- matrix(index_info$use_indices, ny, n_indices, byrow = TRUE)
    } else if (length(index_info$use_indices) == n_indices) {
      use_indices <- matrix(index_info$use_indices, ny, n_indices, byrow = TRUE)
    } else {
      stop("Length of use_indices should be n_indices!")
    }
  }
  
  if (is.null(index_info$use_index_paa)) {
    use_index_paa <- matrix(1, ny, n_indices)
  } else {
    if(length(index_info$use_index_paa) == 1) {
      index_info$use_index_paa <- rep(index_info$use_index_paa, n_indices)
      use_index_paa <- matrix(index_info$use_index_paa, ny, n_indices, byrow = TRUE)
    } else if (length(index_info$use_index_paa) == n_indices) {
      use_index_paa <- matrix(index_info$use_index_paa, ny, n_indices, byrow = TRUE)
    } else {
      stop("Length of use_index_paa should be n_indices!")
    }
    use_index_paa[,which(index_info$use_indices==0)] = 0
  }
  
  
  if (is.null(index_info$units_indices)) {
    units_indices <- rep(2, n_indices) # biomass (1) or numbers (2)
  } else {
    if(length(index_info$units_indices) == 1) {
      units_indices <- rep(index_info$units_indices, n_indices)
    } else if (length(index_info$units_indices) == n_indices) {
      units_indices <- index_info$units_indices
    } else {
      stop("Length of units_indices should be n_indices!")
    }
  }
  
  if (is.null(index_info$units_index_paa)) {
    units_index_paa <- rep(2, n_indices) # biomass (1) or numbers (2)
  } else {
    if(length(index_info$units_index_paa) == 1) {
      units_index_paa <- rep(index_info$units_index_paa, n_indices)
    } else if (length(index_info$units_index_paa) == n_indices) {
      units_index_paa <- index_info$units_index_paa
    } else {
      stop("Length of units_index_paa should be n_indices!")
    }
  }
  
  index_info$use_indices <- use_indices
  index_info$use_index_paa <- use_index_paa
  index_info$units_indices <- units_indices
  index_info$units_index_paa <- units_index_paa
  
  # Handle index_regions
  if (is.null(index_pointer)) {
    if (n_regions == 1) {
      index_info$index_regions <- rep(1, n_indices)
    } else if (n_regions > 1) {
      # Check if n_indices is evenly divisible by n_regions
      if (n_indices %% n_regions != 0) {
        warning("The number of indices is not evenly divisible by the number of regions. Please specify 'index_pointer' manually.")
      } else {
        # Assign index_regions sequentially as 1, 1, 2, 2
        index_info$index_regions <- rep(1:n_regions, each = n_indices / n_regions)
      }
    }
  } else {
    if (length(index_pointer) != n_indices) {
      stop("index_pointer must have length equal to n_indices.")
    }
    index_info$index_regions <- index_pointer
  }
  
  # Set remaining index information
  index_info$index_cv <- index_cv.input
  
  # Calculate index seasons and fracyr_indices using fracyr_indices.input, similar to the fracyr_SSB logic
  index_info$index_seasons <- rep(NA, n_indices)
  index_info$fracyr_indices <- matrix(NA, ny, n_indices)
  
  for (s in 1:n_indices) {
    # Find the cumulative starting times for seasons
    int_starts <- cumsum(c(0, basic_info$fracyr_seasons))
    
    # Determine the season corresponding to the survey fraction for each index
    valid_indices <- which(int_starts <= fracyr_indices.input[s])
    
    if (length(valid_indices) > 0) {
      ind <- max(valid_indices)
    } else {
      stop(paste("No valid season found for index", s, ". Check fracyr_indices.input and fracyr_seasons."))
    }
    
    # Assign season and fractional year value to index_info
    index_info$index_seasons[s] <- ind
    index_info$fracyr_indices[, s] <- fracyr_indices.input[s] - int_starts[ind]
    
    # Check if any value in index_info$fracyr_indices[, s] is equal to 0
    if (any(index_info$fracyr_indices[, s] == 0)) {
      warning(paste("Survey fraction for index", s, "is equal to 0, indicating the survey is happening at the edge of a season. This can cause issues in calculations."))
    }
  }
  
  for (s in 1:n_stocks) {
    # Find the cumulative starting times for seasons
    int_starts <- cumsum(c(0, basic_info$fracyr_seasons))
    
    # Determine the season corresponding to the spawning fraction for each stock
    valid_indices <- which(int_starts <= fracyr_spawn)
    
    if (length(valid_indices) > 0) {
      ind <- max(valid_indices)
    } else {
      stop(paste("No valid season found for spawning for stock", s, ". Check fracyr_spawn and fracyr_seasons."))
    }
    
    # Assign season and fractional year value to basic_info
    basic_info$spawn_seasons[s] <- ind
    basic_info$fracyr_SSB[, s] <- fracyr_spawn - int_starts[ind]
    
    # Check if any value in basic_info$fracyr_SSB[, s] is equal to 0
    if (any(basic_info$fracyr_SSB[, s] == 0)) {
      warning(paste("Spawning fraction for stock", s, "is equal to 0, indicating spawning is happening at the edge of a season. This can cause issues in calculations."))
    }
  }
  
  # Handling onto_move
  if (is.null(basic_info$onto_move)) {
    onto_move = array(0, dim = c(n_stocks,n_regions,n_regions-1))
  } else {
    if (sum(basic_info$onto_move) == 0) {
      onto_move = array(0, dim = c(n_stocks,n_regions,n_regions-1))
    } else {
      if (all(dim(basic_info$onto_move) == c(n_stocks, n_regions, n_regions-1)) && is.array(basic_info$onto_move)) {
        if (all(basic_info$onto_move %in% 0:5)) {
          onto_move = basic_info$onto_move
        } else {
          stop("onto_move must only contain integers between 1 and 5.")
        }
      } else{
        if (basic_info$onto_move %in% 0:5) {
          onto_move = array(basic_info$onto_move, dim = c(n_stocks,n_regions,n_regions-1))
        } else {
          stop("\n Dimention is not correct!")
        }
      }
    }
  }
  
  # Handling onto_move_pars
  if(is.null(onto_move_pars)) onto_move_pars = array(1, dim = c(n_stocks, n_regions, n_regions-1, 4))
  
  if (any(onto_move != 0)) {
    if (is.null(basic_info$onto_move_pars)) {
      # Default: set to 1 if not provided
      onto_move_pars = array(1, dim = c(n_stocks, n_regions, n_regions-1, 4))
    } else if (is.array(basic_info$onto_move_pars)) {
      # Check if the provided array matches the expected dimensions
      if (dim(basic_info$onto_move_pars) != c(n_stocks, n_regions, n_regions-1, 4)) {
        stop("onto_move_pars array must have dimensions (n_stocks, n_regions, n_regions-1, 4).")
      } else {
        onto_move_pars = basic_info$onto_move_pars
      }
    } else if (length(basic_info$onto_move_pars) == 1) {
      onto_move_pars = array(rep(basic_info$onto_move_pars, each = n_stocks * n_regions), dim = c(n_stocks, n_regions, n_regions-1, 4))
    } else if (length(basic_info$onto_move_pars) == 2) {
      onto_move_pars = array(rep(basic_info$onto_move_pars, each = n_stocks * n_regions), dim = c(n_stocks, n_regions, n_regions-1, 4))
    } else if (length(basic_info$onto_move_pars) == 4) {
      onto_move_pars = array(rep(basic_info$onto_move_pars, each = n_stocks * n_regions), dim = c(n_stocks, n_regions, n_regions-1, 4))
    } else {
      stop("onto_move_pars must be a vector of length = 1,2,4 or an array with dimensions (n_stocks, n_regions, n_regions-1, 4).")
    }
  }
  
  age_mu_devs <- array(0, dim = c(n_stocks, n_regions, n_regions-1, n_ages))
  
  if (any(onto_move == 5)) {
    if (!is.null(basic_info$age_mu_devs)) {
      if (all(dim(basic_info$age_mu_devs) == c(n_stocks, n_regions, n_regions-1, n_ages))) {
        if (is.array(basic_info$age_mu_devs)) {
          for (s in 1:n_stocks) {
            for (r in 1:n_regions) {
              for (rr in 1:(n_regions - 1)) {
                if (onto_move[s, r, rr] == 5) {
                  age_mu_devs[s, r, rr, ] <- basic_info$age_mu_devs[s, r, rr, ]
                }
              }
            }
          }
        } else if (length(basic_info$age_mu_devs) == n_ages) {
          for (s in 1:n_stocks) {
            for (r in 1:n_regions) {
              for (rr in 1:(n_regions - 1)) {
                if (onto_move[s, r, rr] == 5) {
                  age_mu_devs[s, r, rr, ] <- basic_info$age_mu_devs
                }
              }
            }
          }
        } else if (length(basic_info$age_mu_devs) == 1) {
          for (s in 1:n_stocks) {
            for (r in 1:n_regions) {
              for (rr in 1:(n_regions - 1)) {
                if (onto_move[s, r, rr] == 5) {
                  age_mu_devs[s, r, rr, ] <- rep(basic_info$age_mu_devs, n_ages)
                }
              }
            }
          }
        } else {
          stop("age_mu_devs should be one of the followings:\n 1) an array with dim = n_stocks, n_regions, n_regions-1, n_ages
               \n 2) a vector with length = n_ages which will apply to all stocks and regions that have onto_move = 5
               \n 3) a single value which will apply to all ages, stocks, and regions that have onto_move = 5")
        }
      }
    } else {
      warning("onto_move = 5 but age_mu_devs is not specified!")
    }
  }
  
  basic_info$onto_move = onto_move
  basic_info$onto_move_pars = onto_move_pars
  basic_info$age_mu_devs = age_mu_devs
  
  par_inputs <- list(
    n_stocks = n_stocks,
    n_regions = n_regions,
    n_indices = n_indices,
    n_fleets = n_fleets,
    n_seasons = n_seasons,
    life_history = life_history,
    n_ages = n_ages,
    Fbar_ages = Fbar_ages,
    recruit_model = recruit_model,
    
    # Fishing Mortality Information (from F_info)
    F.year1 = F.year1,
    Fhist = Fhist,
    Fmax = Fmax,
    Fmin = Fmin,
    change_time = change_time,
    user_F = user_F,
    
    # Catch Information 
    
    catch_cv = catch_info_list$catch_cv,
    catch_Neff = catch_info_list$catch_Neff,
    use_agg_catch = catch_info_list$use_agg_catch,
    use_catch_paa = catch_info_list$use_catch_paa,
    
    # Index Information 
    index_cv = index_info_list$index_cv,
    index_Neff = index_info_list$index_Neff,
    fracyr_indices = index_info_list$fracyr_indices,
    q = index_info_list$q,
    use_indices = index_info_list$use_indices,
    use_index_paa = index_info_list$use_index_paa,
    units_indices = index_info_list$units_indices,
    units_index_paa = index_info_list$units_index_paa,

    # Spawning and Seasons Information
    fracyr_spawn = fracyr_spawn,
    fracyr_seasons = fracyr_seasons,
    
    # Pointers and User-Defined Inputs
    fleet_pointer = fleet_pointer,
    index_pointer = index_pointer,
    user_waa = user_waa,
    user_maturity = user_maturity,
    
    # Bias Correction Flags
    bias_correct_process = bias.correct.process,
    bias_correct_observation = bias.correct.observation,
    bias_correct_BRPs = bias.correct.BRPs,
    
    # Migration and Meta-population
    mig_type = mig_type,
    XSPR_R_opt = XSPR_R_opt,
    move_dyn = move_dyn,
    onto_move = onto_move,
    onto_move_pars = onto_move_pars,
    
    # Recruitment and Movement Trends
    apply_re_trend = apply_re_trend,
    trend_re_rate = trend_re_rate,
    apply_mu_trend = apply_mu_trend,
    trend_mu_rate = trend_mu_rate,
    age_mu_devs = age_mu_devs
  )
  
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
