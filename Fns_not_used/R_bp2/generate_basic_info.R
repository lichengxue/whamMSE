#' Generate Basic Information for Multi-WHAM Input
#'
#' The `generate_basic_info()` function creates biological and fishery information necessary 
#' for simulation-estimation and management strategy evaluation (MSE). This function defines 
#' stock structure, fleet operations, life history parameters, fishing mortality, and 
#' movement dynamics, among other essential components.
#'
#' @param n_stocks Integer. Number of stocks.
#' @param n_regions Integer. Number of regions.
#' @param n_indices Integer. Number of indices.
#' @param n_fleets Integer. Number of fleets.
#' @param n_seasons Integer. Number of seasons.
#' @param base.years Integer vector. Model years (or "burn-in" period in MSE).
#' @param n_feedback_years Integer. Optional. Number of years in the MSE feedback loop.
#' @param life_history Character. Fish life history, parameters obtained from 
#'   \href{https://doi.org/10.1139/cjfas-2016-0381}{Wiedenmann et al. 2017}. Options:
#'   \itemize{
#'     \item `"medium"` - Median-lived species (default).
#'     \item `"short"` - Short-lived species.
#'     \item `"long"` - Long-lived species.
#'   }
#' @param n_ages Integer. Number of ages.
#' @param F_info List. Historical fishing pressure or user-specified fishing mortality values:
#'   \itemize{
#'     \item `$F.year1` - Fishing mortality in the first year.
#'     \item `$Fhist` - Pattern of historical fishing mortality. Options include:
#'       \itemize{
#'         \item `"constant"` - Constant across years.
#'         \item `"updown"` - Increase to `Fmax` until change point, then decrease to `Fmin`.
#'         \item `"downup"` - Decrease to `Fmin` until change point, then increase to `Fmax`.
#'         \item `"F-H-L"` - Constant `F.year1*Fmax` until change point, followed by constant `F.year1*Fmin`.
#'         \item `"F-L-H"` - Constant `F.year1*Fmin` until change point, followed by constant `F.year1*Fmax`.
#'       }
#'     \item `$Fmax` - Maximum F (or multiplier when `Fhist = "F-H-L"`).
#'     \item `$Fmin` - Minimum F (or multiplier when `Fhist = "F-H-L"`).
#'     \item `$change_time` - Proportion of the time series (0–1) indicating when the change in fishing mortality occurs.
#'     \item `$user_F` - Optional matrix (`n_years x n_fleets`) of user-specified fishing mortality values.
#'   }
#' @param catch_info List. Fleet catch information:
#'   \itemize{
#'     \item `"$catch_cv"` – Vector of length `n_fleets`, or a single value applied to all fleets, specifying CV for each fleet catch.
#'     \item `"$catch_Neff"` – Vector of length `n_fleets`, or a single value applied to all fleets, specifying effective sample sizes for fleet catch age composition.
#'     \item `"$use_agg_catch"` – Vector of length `n_fleets`, or a single value applied to all fleets, with 0/1 flags indicating whether aggregate catch data are used.
#'     \item `"$use_catch_paa"` – Vector of length `n_fleets`, or a single value applied to all fleets, with 0/1 flags indicating whether age composition data for catch are used.
#'   }
#' @param index_info List. Survey index information:
#'   \itemize{
#'     \item `"$index_cv"` – Vector of length `n_indices`, or a single value applied to all indices, specifying CV for each survey index.
#'     \item `"$index_Neff"` – Vector of length `n_indices`, or a single value applied to all indices, specifying effective sample sizes for survey age composition.
#'     \item `"$fracyr_indices"` – Vector of length `n_indices`, or a single value applied to all indices, specifying the fraction of the year when surveys occur.
#'     \item `"$q"` – Vector of length `n_indices`, or a single value applied to all indices, specifying survey catchabilities.
#'     \item `"$use_indices"` – Vector of length `n_indices`, or a single value applied to all indices, with 0/1 flags indicating whether aggregate index data are used.
#'     \item `"$use_index_paa"` – Vector of length `n_indices`, or a single value applied to all indices, with 0/1 flags indicating whether age composition data for indices are used.
#'     \item `"$units_indices"` – Vector of length `n_indices`, or a single value applied to all indices, where 1 = biomass and 2 = numbers for aggregate index data.
#'     \item `"$units_index_paa"` – Vector of length `n_indices`, or a single value applied to all indices, where 1 = biomass and 2 = numbers for index age composition data.
#'   }
#' @param fracyr_spawn Numeric. Fraction of the year when spawning occurs.
#' @param fracyr_seasons Numeric vector. Optional. User-defined seasonal fractions summing to 1.
#' @param fleet_regions Integer vector. Optional. User-defined region allocation for each fleet.
#' @param index_regions Integer vector. Optional. User-defined region allocation for each index.
#' @param user_waa Numeric vector/matrix. Optional. User-defined weight-at-age vector/matrix.
#' @param user_maturity Numeric vector. Optional. User-defined maturity-at-age vector/matrix.
#' @param bias.correct.process Logical. Whether to apply process error bias correction.
#' @param bias.correct.observation Logical. Whether to apply observation error bias correction.
#' @param bias.correct.BRPs Logical. Whether to apply biological reference point bias correction.
#' @param mig_type Integer. Migration type:
#'   \itemize{
#'     \item `0` - Movement occurs after survival.
#'     \item `1` - Movement and mortality occur simultaneously.
#'   }
#' @param XSPR_R_opt Integer. Recruitment that will be used to calculate SPR-based reference points:
#'   \itemize{
#'     \item `1` - Annual recruitment estimates.
#'     \item `2` - Average recruitment estimates (Default).
#'     \item `3` - Annual expected recruitment.
#'     \item `4` - Average expected recruitment.
#'     \item `5` - Bias-corrected expected recruitment.
#'   }
#' @param move_dyn Integer. Movement dynamics:
#'   \itemize{
#'     \item `0` - Natal homing.
#'     \item `1` - Meta-population.
#'   }
#' @param onto_move Integer array (`n_stocks x n_regions x (n_regions - 1)`). Age-specific movement type:
#'   \itemize{
#'     \item `0` - Same movement rate across ages (default).
#'     \item `1` - Increasing logistic (2 parameters).
#'     \item `2` - Decreasing logistic (2 parameters).
#'     \item `3` - Double-logistic (4 parameters).
#'     \item `4` - Double-normal (4 parameters).
#'     \item `5` - User-specified age-specific movement rate.
#'   }
#' @param onto_move_pars Numeric array (`n_stocks x n_regions x (n_regions - 1) x 4`). Parameters controlling age-specific movement patterns. 
#' @param apply_re_trend Logical. Whether to apply a temporal trend to the movement random effects (i.e., year-to-year deviations).
#' @param trend_re_rate Numeric. Annual rate of change applied to the random effects when `apply_re_trend = TRUE`.
#' @param apply_mu_trend Logical. Whether to apply a temporal trend to the movement mean (`mu`), representing the expected movement rate over time.
#' @param trend_mu_rate Numeric. Annual rate of change applied to the movement mean when `apply_mu_trend = TRUE`.
#' @param age_mu_devs Numeric array (`n_stocks x n_regions x (n_regions - 1) x n_ages`). User-specified deviations in mean movement rates across ages.
#' @return A list containing:
#'   \itemize{
#'     \item `$basic_info` - General biological and fishery details.
#'     \item `$catch_info` - Fleet-specific catch data.
#'     \item `$index_info` - Survey-related data.
#'     \item `$F` - Fishing mortality settings.
#'     \item `$par_inputs` - A list of input parameters for model setup.
#'   }
#'
#' @export
#'
#' @seealso \code{\link{prepare_wham_input}}
#'
#' @examples
#' \dontrun{
#' data <- generate_basic_info(n_stocks = 2, n_regions = 2, n_fleets = 3, base.years = 1:30)
#' }

generate_basic_info <- function(n_stocks = 2,
                                n_regions = 2,
                                n_fleets = 2,
                                n_indices = 2,
                                n_seasons = 5,
                                base.years = 1:20,
                                n_feedback_years = 0,
                                life_history = "medium",
                                n_ages = 12,
                                F_info = list(F.year1 = 0.2, Fhist = "constant", Fmax = 0.2, Fmin = 0.2, change_time = 0.5, user_F = NULL, F_feedback = NULL),
                                catch_info = list(catch_cv = 0.1, catch_Neff = 100, use_agg_catch = 1, use_catch_paa = 1),
                                index_info = list(index_cv = 0.1, index_Neff = 100, fracyr_indices = 0.5, q = 0.2,
                                                  use_indices = 1, use_index_paa = 1, units_indices = 2, units_index_paa = 2),
                                user_waa = NULL, 
                                fracyr_spawn = 0.5,
                                fracyr_seasons = NULL,
                                fleet_regions = NULL,
                                index_regions = NULL,
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
  
  if (!check_dimensions(n_stocks, n_regions, n_fleets, n_indices)) cat("\n Note: n_stocks, n_regions, n_fleets, n_indices are not the same.\n")
  
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
  
  # basic_info$recruit_model = recruit_model
  
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
  
  fracyr_seasons = basic_info$fracyr_seasons
  basic_info$fracyr_SSB <- matrix(0, ny, n_stocks) # Assume fish is recruited at the beginning of the year
  basic_info$fracyr_spawn = fracyr_spawn # rep(fracyr_spawn, n_stocks) different spawning season???

  if(is.null(catch_info$catch_cv)) catch_info$catch_cv = 0.1
  if(is.null(catch_info$catch_Neff)) catch_info$catch_Neff = 100
  if(is.null(catch_info$use_agg_catch)) catch_info$use_agg_catch = 1
  if(is.null(catch_info$use_catch_paa)) catch_info$use_catch_paa = 1
  
  if(is.null(index_info$index_cv)) index_info$index_cv = 0.1
  if(is.null(index_info$index_Neff)) index_info$index_Neff = 100
  if(is.null(index_info$fracyr_indices)) index_info$fracyr_indices = 0.5
  if(is.null(index_info$q)) index_info$q = 0.2
  if(is.null(index_info$use_indices)) index_info$use_indices = 1
  if(is.null(index_info$use_index_paa)) index_info$use_index_paa = 1
  if(is.null(index_info$units_indices)) index_info$units_indices = 2
  if(is.null(index_info$units_index_paa)) index_info$units_index_paa = 2
  
  catch_info_list <- catch_info
  index_info_list <- index_info
  
  # Fishing mortality
  nby <- length(base.years)
  
  # Extract relevant values directly from F_info
  F.year1     <- F_info$F.year1
  Fhist       <- F_info$Fhist
  Fmax        <- F_info$Fmax
  Fmin        <- F_info$Fmin
  change_time <- F_info$change_time
  user_F      <- F_info$user_F
  F_feedback  <- F_info$F_feedback
  
  if (!is.null(user_F)) {
    # Check user-specified F matrix
    if (!all(dim(user_F) == c(nby, n_fleets))) {
      stop("user_F must be a matrix with dimensions n_years x n_fleets.")
    }
    F_vals <- user_F
  } else {
    if (is.null(F.year1)) stop("Users must specify initial F!")
    
    if (Fhist == "constant") {
      F_vals <- matrix(F.year1, nby, n_fleets)
    } else {
      if (is.null(Fmax)) stop("Fmax must be specified!")
      if (is.null(Fmin)) stop("Fmin must be specified!")
      if (is.null(change_time)) stop("change_time must be specified!")
      
      mid <- ceiling(nby * change_time)
      
      if (Fhist == "updown") {
        F_vals <- matrix(c(seq(F.year1, Fmax, length.out = mid),
                           seq(Fmax, Fmin, length.out = nby - mid + 1)[-1]), nby, n_fleets)
      } else if (Fhist == "downup") {
        F_vals <- matrix(c(seq(F.year1, Fmin, length.out = mid),
                           seq(Fmin, Fmax, length.out = nby - mid + 1)[-1]), nby, n_fleets)
      } else if (Fhist == "F-H-L") {
        F_vals <- matrix(F.year1 * Fmin, nby, n_fleets)
        F_vals[1:mid, ] <- F.year1 * Fmax
      }
      else if (Fhist == "F-L-H") {
        F_vals <- matrix(F.year1 * Fmin, nby, n_fleets)
        F_vals[-(1:mid), ] <- F.year1 * Fmax
      }
    }
  }
  
  if (n_feedback_years > 0) {
    if (n_fleets > 1) F_vals <- rbind(F_vals, F_vals[rep(nby, n_feedback_years), , drop = FALSE])
    else F_vals <- rbind(F_vals, matrix(F_vals[rep(nby, n_feedback_years)], ncol = 1))
    if (is.null(F_feedback)) {
      F_vals[(nby+1):(nby+n_feedback_years),] = 0.1 # assume 0.1
    } else {
      F_vals[(nby+1):(nby+n_feedback_years),] = F_feedback
    }
  }
  
  F_info = list(F = F_vals, F_config = 2)
  
  # if (is.null(Fbar_ages)) {
  #   basic_info$Fbar_ages = as.integer(na)
  # } else {
  #   basic_info$Fbar_ages = as.integer(Fbar_ages)
  # }
  basic_info$Fbar_ages = as.integer(na)
  
  # Maturity at age
  if (is.null(user_maturity)) {
    maturity <- Generate_Maturity(life_history, na)
    basic_info$maturity <- array(NA, dim = c(n_stocks, ny, na))
    for (i in 1:n_stocks) basic_info$maturity[i, , ] <- t(matrix(maturity, na, ny))
  } else {
    if (is.matrix(user_maturity) && ncol(user_maturity) == na && nrow(user_maturity) == n_stocks) {
      maturity <- user_maturity
      basic_info$maturity <- array(NA, dim = c(n_stocks, ny, na))
      for (i in 1:n_stocks) basic_info$maturity[i, , ] <- t(matrix(maturity[i,], na, ny))
    }
    else if (is.vector(user_maturity) && length(user_maturity) == na) {
      maturity <- user_maturity
      basic_info$maturity <- array(NA, dim = c(n_stocks, ny, na))
      for (i in 1:n_stocks) basic_info$maturity[i, , ] <- t(matrix(maturity, na, ny))
    } 
    else if (is.array(user_maturity) && all(dim(user_maturity) == c(n_stocks, ny, na))) {
      basic_info$maturity = user_maturity
    } else {
      user_maturity = array(NA, dim = c(n_stocks, ny, na))
      warnings("user_maturity must have correct dimensions!")
    }
  } 
  
  user_maturity = basic_info$maturity
  
  # Weight at age
  nwaa <- n_fleets + n_regions + n_indices + n_stocks
  basic_info$waa <- array(NA, dim = c(nwaa, ny, na))
  
  if (is.null(user_waa)) {
    
    W <- Generate_WAA(life_history, na)
    for (i in 1:nwaa) basic_info$waa[i, , ] <- t(matrix(W, na, ny))
    basic_info$waa_pointer_fleets   <- 1:n_fleets
    basic_info$waa_pointer_totcatch <- (n_fleets + 1):(n_fleets + n_regions)
    basic_info$waa_pointer_indices  <- (n_fleets + n_regions + 1):(n_fleets + n_regions + n_indices)
    basic_info$waa_pointer_ssb      <- (n_fleets + n_regions + n_indices + 1):(n_fleets + n_regions + n_indices + n_stocks)
    basic_info$waa_pointer_M        <- basic_info$waa_pointer_ssb
    
    user_waa$waa                  <- basic_info$waa
    user_waa$waa_pointer_fleets   <- basic_info$waa_pointer_fleets
    user_waa$waa_pointer_totcatch <- basic_info$waa_pointer_totcatch
    user_waa$waa_pointer_indices  <- basic_info$waa_pointer_indices
    user_waa$waa_pointer_ssb      <- basic_info$waa_pointer_ssb 
    user_waa$waa_pointer_M        <- basic_info$waa_pointer_M 

  } else {
    
    if(is.null(user_waa$waa)) stop("waa must be specified! (dim = c(n_fleets + n_regions + n_indices + n_stocks, n_years_model, n_ages))")
    if(is.null(user_waa$waa_pointer_fleets)) stop("waa_pointer_fleets must be specified! (dim = n_fleets)")
    if(is.null(user_waa$waa_pointer_totcatch)) stop("waa_pointer_totcatch must be specified! (dim = n_regions)")
    if(is.null(user_waa$waa_pointer_indices)) stop("waa_pointer_indices must be specified! (dim = n_indices)")
    if(is.null(user_waa$waa_pointer_ssb)) stop("waa_pointer_ssb must be specified! (dim = n_regions)")
    if(is.null(user_waa$waa_pointer_M)) stop("waa_pointer_M must be specified! (dim = n_regions)")
    
    basic_info$waa                  <- user_waa$waa
    basic_info$waa_pointer_fleets   <- user_waa$waa_pointer_fleets
    basic_info$waa_pointer_totcatch <- user_waa$waa_pointer_totcatch
    basic_info$waa_pointer_indices  <- user_waa$waa_pointer_indices
    basic_info$waa_pointer_ssb      <- user_waa$waa_pointer_ssb
    basic_info$waa_pointer_M        <- user_waa$waa_pointer_M
  }
  
  # else {
  #   if (length(user_waa) == na) {
  #     W <- user_waa
  #     for (i in 1:nwaa) basic_info$waa[i, , ] <- do.call(rbind, replicate(ny, W, simplify = FALSE))
  #   } else if (is.matrix(user_waa) && dim(user_waa)[1] == nwaa && dim(user_waa)[2] == na) {
  #     W <- user_waa
  #     for (i in 1:nwaa) basic_info$waa[i, , ] <- do.call(rbind, replicate(ny, W[i,], simplify = FALSE))
  #   } else {
  #     warnings("Dimension of W should be either a vector of n_ages or a matrix with nrow = c(n_fleets + n_regions + n_indices + n_stocks) and ncol = n_ages!")
  #   }
  # }
  
  #waa_info = list(waa = NULL, waa_pointer_fleets = NULL, waa_pointer_totcatch = NULL, waa_pointer_indices = NULL, waa_pointer_ssb = NULL, waa_pointer_M = NULL)

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
    # if(length(catch_info$use_agg_catch) != n_fleets) stop("Length of use_agg_catch should be n_fleets!")
    use_agg_catch <- matrix(catch_info$use_agg_catch, ny, n_fleets, byrow = TRUE)
  }
  
  if (is.null(catch_info$use_catch_paa)) {
    use_catch_paa <- matrix(1, ny, n_fleets)
  } else {
    # if(length(catch_info$use_catch_paa) != n_fleets) stop("Length of use_catch_paa should be n_fleets!")
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
  if (is.null(fleet_regions)) {
    if (n_regions == 1) {
      catch_info$fleet_regions <- rep(1, n_fleets)
    } else if (n_regions > 1) {
      # Check if n_fleets is evenly divisible by n_regions
      if (n_fleets %% n_regions != 0) {
        warning("The number of fleets is not evenly divisible by the number of regions. Please specify 'fleet_regions' manually.")
      } else {
        # Assign fleet_regions sequentially as 1, 1, 2, 2
        catch_info$fleet_regions <- rep(1:n_regions, each = n_fleets / n_regions)
      }
    }
  } else {
    if (length(fleet_regions) != n_fleets) {
      stop("fleet_regions must have length equal to n_fleets.")
    }
    catch_info$fleet_regions <- fleet_regions
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
  index_info$index_cv <- index_cv.input
  
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
  if (is.null(index_regions)) {
    if (n_regions == 1) {
      index_info$index_regions <- rep(1, n_indices)
    } else if (n_regions > 1) {
      # Check if n_indices is evenly divisible by n_regions
      if (n_indices %% n_regions != 0) {
        warning("The number of indices is not evenly divisible by the number of regions. Please specify 'index_regions' manually.")
      } else {
        # Assign index_regions sequentially as 1, 1, 2, 2
        index_info$index_regions <- rep(1:n_regions, each = n_indices / n_regions)
      }
    }
  } else {
    if (length(index_regions) != n_indices) {
      stop("index_regions must have length equal to n_indices.")
    }
    index_info$index_regions <- index_regions
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
      cat(paste("\nSurvey fraction for index", s, "is equal to 0, indicating the survey is happening at the edge of a season.\n"))
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
      cat("Spawning fraction for stock", s, "is equal to 0. This may cause issues when having multiple seasons.\n")
    }
  }
  
  # Handling onto_move
  if (is.null(basic_info$onto_move)) {
    onto_move = array(0, dim = c(n_stocks,n_regions,n_regions-1))
    # onto_move = NULL
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
    fleet_regions = catch_info$fleet_regions,
    index_regions = index_info$index_regions,
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