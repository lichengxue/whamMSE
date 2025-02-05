#' Generate input data for the assessment model
#'
#' A function to generate the input for the estimation model for management strategy evaluation.
#'
#' @param om Operating model
#' @param em_info A list of information used to generate the operating model
#' @param M_em Natural mortality random effects
#' @param sel_em Selectivity random effects
#' @param NAA_re_em Numbers-at-age random effects
#' @param move_em Movement random effects
#' @param em.opt Movement random effects options
#'   \itemize{
#'     \item \code{$separate.em} TRUE = No Global SPR, FALSE = Global SPR
#'     \item \code{$separate.em.type} only if separate.em = TRUE \cr
#'     {=1} panmictic (spatially-aggregated) \cr
#'     {=2} fleets-as-areas \cr
#'     {=3} n single assessment models (n = n_regions) \cr
#'     \item \code{$do.move} T/F movement is included (use if separate.em = FALSE)
#'     \item \code{$est.move} T/F movement rate is estimated (use if separate.em = FALSE)
#'   }
#' @param em_years Years used in the assessment model
#' @param year.use Number of years used in the assessment model
#' @param age_comp_em Likelihood distribution of age composition data in the assessment model
#'   \itemize{
#'     \item \code{"multinomial"} Default
#'     \item \code{"dir-mult"}
#'     \item \code{"dirichlet-miss0"}
#'     \item \code{"dirichlet-pool0"}
#'     \item \code{"logistic-normal-miss0"}
#'     \item \code{"logistic-normal-ar1-miss0"}
#'     \item \code{"logistic-normal-pool0"}
#'     \item \code{"logistic-normal-01-infl"}
#'     \item \code{"logistic-normal-01-infl-2par"}
#'     \item \code{"mvtweedie"}
#'     \item \code{"dir-mult-linear"}
#'   }
#' @param aggregate_catch_info (optional) User specified list of information for aggregate catch (default: use first fleet). Only when panmictic model with aggregate catch is used.
#'   \itemize{
#'     \item \code{$catch_cv} {vector (n_fleets) of CVs for fleet catch.}
#'     \item \code{$catch_Neff} {vector (n_fleets) of effective sample sizes for fleet catch.}
#'     \item \code{$use_agg_catch} {vector (n_fleets) of 0/1 values flagging whether whether to use aggregate catch.}
#'     \item \code{$use_catch_paa} {vector (n_fleets) of 0/1 values flagging whether to use proportions at age observations.}
#'     \item \code{$fleet_pointer} {vector (n_fleets) of which fleets should be combined. Use 0 to not include that fleet.}
#'     }
#' @param aggregate_index_info (optional) User specified list of information for aggregate index (default: use first index). Only when panmictic model with aggregate index is used.
#'   \itemize{
#'     \item \code{$index_cv} {vector (n_indices) of CVs for survey indices.}
#'     \item \code{$index_Neff} {vector (n_indices) of effective sample sizes for survey indices.}
#'     \item \code{$fracyr_indices} {vector (n_indices) of fractions of the year when each survey is conducted.}
#'     \item \code{$q} {vector (n_indices) of survey catchabilities.}
#'     \item \code{$use_indices} {vector (n_indices) of 0/1 values flagging whether to use aggregate observations.}
#'     \item \code{$use_index_paa} {vector (n_indices) of 0/1 values flagging whether to use proportions at age observations.}
#'     \item \code{$units_indices} {vector (n_indices) of 1/2 values flagging whether aggregate observations are biomass (1) or numbers (2).}
#'     \item \code{$units_index_paa} {vector (n_indices) of 1/2 values flagging whether composition observations are biomass (1) or numbers (2).}
#'     \item \code{$index_pointer} {vector (n_indices) of which indices should be combined. Use 0 to not include that index.}
#'     }
#'     
#' @return a wham input
#'
#' @export
#'
#' @seealso \code{\link{loop_through_fn}}
#' @examples
#' \dontrun{
#' data <- generate_basic_info()
#' input <- prepare_wham_input(basic_info = data)
#' mod <- fit_wham(input, do.fit = FALSE)
#' data <- make_em_input(mod, em.opt = list(separate.em = TRUE, separate.em.type = 1, do.move = FALSE, est.move = FALSE), year.use = 10, em_years = 1973:2022)
#' }

make_em_input <- function(om, em_info, M_em, sel_em, NAA_re_em, move_em,
                          em.opt, em_years, year.use,
                          age_comp_em = "multinomial",
                          aggregate_catch_info = NULL,
                          aggregate_index_info = NULL) {
  
  if (is.null(em.opt)) stop("em.opt must be specified!")
  
  # Determine movement type based on options
  if (em.opt$separate.em) {
    em.opt$do.move <- FALSE
    move.type <- NULL
  } else if (!em.opt$do.move) {
    move.type <- 3 # no movement
  } else if (all(move_em$stock_move)) {
    move.type <- 2 # bidirectional
  } else {
    move.type <- 1 # unidirectional
  }
  
  data <- om$input$data
  
  # Determine which years to use in the estimation model
  if (!is.null(year.use)) {
    if (year.use > length(em_years)) {
      warning("year.use must be <= em_years! Setting year.use to the length of em_years.")
      year.use <- length(em_years)
    }
    ind_em <- (length(em_years) - year.use + 1):length(em_years)
    em_years <- tail(em_years, year.use)
  } else {
    year.use <- length(em_years)
    ind_em <- (length(em_years) - year.use + 1):length(em_years)
  }
  
  if (em.opt$separate.em) {
    if (em.opt$separate.em.type == 1) {
      
      # Placeholder for handling one-region model
      if (is.null(aggregate_catch_info)) {
        em_info$par_inputs$catch_cv <- em_info$par_inputs$catch_cv[1]
        em_info$par_inputs$catch_Neff <- em_info$par_inputs$catch_Neff[1]
        em_info$par_inputs$use_agg_catch <- em_info$par_inputs$use_agg_catch[1]
        em_info$par_inputs$use_catch_paa <- em_info$par_inputs$use_catch_paa[1]
      } else {
        em_info$par_inputs$catch_cv <- aggregate_catch_info$catch_cv
        em_info$par_inputs$catch_Neff <- aggregate_catch_info$catch_Neff
        em_info$par_inputs$use_agg_catch <- em_info$par_inputs$use_agg_catch
        em_info$par_inputs$use_catch_paa <- em_info$par_inputs$use_catch_paa
      }
      if (is.null(aggregate_index_info)) {
        em_info$par_inputs$index_cv <- em_info$par_inputs$index_cv[1]
        em_info$par_inputs$index_Neff <- em_info$par_inputs$index_Neff[1]
        em_info$par_inputs$fracyr_indices <- em_info$par_inputs$fracyr_indices[1]
        em_info$par_inputs$q <- em_info$par_inputs$q[1]
        em_info$par_inputs$use_indices <- em_info$par_inputs$use_indices[1]
        em_info$par_inputs$use_index_paa <- em_info$par_inputs$use_index_paa[1]
        em_info$par_inputs$units_indices <- em_info$par_inputs$units_indices[1]
        em_info$par_inputs$units_index_paa <- em_info$par_inputs$units_index_paa[1]
      } else {
        em_info$par_inputs$index_cv <- aggregate_index_info$index_cv
        em_info$par_inputs$index_Neff <- aggregate_index_info$index_Neff
        em_info$par_inputs$fracyr_indices <- aggregate_index_info$fracyr_indices
        em_info$par_inputs$q <- aggregate_index_info$q
        em_info$par_inputs$use_indices <- aggregate_index_info$use_indices
        em_info$par_inputs$use_index_paa <- aggregate_index_info$use_index_paa
        em_info$par_inputs$units_indices <- aggregate_index_info$units_indices
        em_info$par_inputs$units_index_paa <- aggregate_index_info$units_index_paa
      }
      
      n_fleets <- 1
      n_indices <- 1
      if(!is.null(aggregate_index_info$n_fleets)) n_fleets <- aggregate_index_info$n_fleets 
      if(!is.null(aggregate_index_info$n_indices)) n_indices <- aggregate_index_info$n_indices  
      
      info <- generate_basic_info_em(em_info, em_years, n_stocks = 1, n_regions = 1, n_fleets = n_fleets, n_indices = n_indices)
      basic_info <- info$basic_info
      
      # Override any movement or trend information
      basic_info$move_dyn <- 0
      basic_info$onto_move <- matrix(0)
      basic_info$apply_re_trend <- 0
      basic_info$apply_mu_trend <- 0
      
      # Fill in the data from the operating model simulation
      if (n_fleets == 1) {
        info$catch_info$agg_catch <- data$agg_catch[ind_em, , drop = FALSE]
        info$catch_info$agg_catch <- matrix(rowSums(info$catch_info$agg_catch), ncol = 1)
        info$catch_info$catch_paa <- data$catch_paa[, ind_em, , drop = FALSE]
        catch <- data$agg_catch[ind_em, , drop = FALSE]
        ratio <- data$catch_paa[, ind_em, ]
        result <- 0
        for (i in 1:dim(ratio)[1]) {
          tmp <- ratio[i, , ] * catch[, i]
          result <- result + tmp
        }
        result <- t(apply(result, 1, function(row) ifelse(sum(row) == 0, 0, row / sum(row))))
        info$catch_info$catch_paa <- array(result, dim = c(1, nrow(result), ncol(result)))
      } else if (n_fleets > 1) {
        fleet_pointer <- aggregate_catch_info$fleet_pointer
        valid_pointer = unique(fleet_pointer[fleet_pointer > 0])
        
        info$catch_info$agg_catch = matrix(NA, length(ind_em), length(unique(valid_pointer)))
        info$catch_info$catch_paa = array(NA, dim = c(length(unique(valid_pointer)), length(ind_em), basic_info$n_ages))
        
        for (f in unique(valid_pointer)) {
          agg_catch.tmp <- data$agg_catch[ind_em, which(fleet_pointer == f), drop = FALSE]
          agg_catch.tmp <- matrix(rowSums(agg_catch.tmp), ncol = 1)
          catch_paa.tmp <- data$catch_paa[which(fleet_pointer == f), ind_em, , drop = FALSE]
          catch <- data$agg_catch[ind_em, which(fleet_pointer == f), drop = FALSE]
          ratio <- data$catch_paa[which(fleet_pointer == f), ind_em, ,drop = FALSE]
          result <- 0
          for (i in 1:dim(ratio)[1]) {
            tmp <- ratio[i, , ] * catch[, i]
            result <- result + tmp
          }
          result <- t(apply(result, 1, function(row) ifelse(sum(row) == 0, 0, row / sum(row))))
          catch_paa.tmp <- array(result, dim = c(1, nrow(result), ncol(result)))
          
          info$catch_info$agg_catch[,f] = agg_catch.tmp
          info$catch_info$catch_paa[f,,] = catch_paa.tmp[1,,]
        }
      }
        
      if (n_indices == 1) {
        info$index_info$agg_indices <- data$agg_indices[ind_em, , drop = FALSE]
        info$index_info$agg_indices <- matrix(rowSums(info$index_info$agg_indices), ncol = 1)
        info$index_info$index_paa <- data$index_paa[, ind_em, , drop = FALSE]
        catch <- data$agg_indices[ind_em, , drop = FALSE]
        ratio <- data$index_paa[, ind_em, ]
        result <- 0
        for (i in 1:dim(ratio)[1]) {
          tmp <- ratio[i, , ] * catch[, i]
          result <- result + tmp
        }
        result <- t(apply(result, 1, function(row) ifelse(sum(row) == 0, 0, row / sum(row))))
        info$index_info$index_paa <- array(result, dim = c(1, nrow(result), ncol(result)))
      } else if (n_indices > 1) {
        
        index_pointer <- aggregate_index_info$index_pointer
        valid_pointer = unique(index_pointer[index_pointer > 0])
        
        info$index_info$agg_indices = matrix(NA, length(ind_em), length(unique(index_pointer)))
        info$index_info$index_paa = array(NA, dim = c(length(unique(index_pointer)), length(ind_em), basic_info$n_ages))
        
          for (f in unique(valid_pointer)) {
            agg_indices.tmp <- data$agg_indices[ind_em, which(index_pointer == f), drop = FALSE]
            agg_indices.tmp <- matrix(rowSums(agg_indices.tmp), ncol = 1)
            index_paa.tmp <- data$index_paa[which(index_pointer == f), ind_em, , drop = FALSE]
            catch <- data$agg_indices[ind_em, which(index_pointer == f), drop = FALSE]
            ratio <- data$index_paa[which(index_pointer == f), ind_em, , drop = FALSE]
            result <- 0
            for (i in 1:dim(ratio)[1]) {
              tmp <- ratio[i, , ] * catch[, i]
              result <- result + tmp
            }
            result <- t(apply(result, 1, function(row) ifelse(sum(row) == 0, 0, row / sum(row))))
            index_paa.tmp <- array(result, dim = c(1, nrow(result), ncol(result)))
            
            info$index_info$agg_indices[,f] = agg_indices.tmp
            info$index_info$index_paa[f,,] = index_paa.tmp[1,,]
          }
      }
      
      em_input <- prepare_wham_input(
        basic_info = basic_info,
        selectivity = sel_em,
        M = M_em,
        NAA_re = NAA_re_em,
        move = NULL,
        age_comp = age_comp_em,
        catch_info = info$catch_info,
        index_info = info$index_info,
        F = info$F
      )
    }
    
    if (em.opt$separate.em.type == 2) {
      # Fleets-as-areas
      n_fleets <- data$n_fleets
      n_indices <- data$n_indices
      
      info <- generate_basic_info_em(em_info, em_years, n_stocks = 1, n_regions = 1, n_fleets = n_fleets, n_indices = n_indices)
      basic_info <- info$basic_info
      
      # Override any movement or trend information
      basic_info$move_dyn <- 0
      basic_info$onto_move <- matrix(0)
      basic_info$apply_re_trend <- 0
      basic_info$apply_mu_trend <- 0
      
      # Fill in the data from the operating model simulation
      info$catch_info$agg_catch <- data$agg_catch[ind_em, , drop = FALSE]
      info$index_info$agg_indices <- data$agg_indices[ind_em, , drop = FALSE]
      info$catch_info$catch_paa <- data$catch_paa[, ind_em, , drop = FALSE]
      info$index_info$index_paa <- data$index_paa[, ind_em, , drop = FALSE]
      
      em_input <- prepare_wham_input(
        basic_info = basic_info,
        selectivity = sel_em,
        M = M_em,
        NAA_re = NAA_re_em,
        move = NULL,
        age_comp = age_comp_em,
        catch_info = info$catch_info,
        index_info = info$index_info,
        F = info$F
      )
    }
    
    if (em.opt$separate.em.type == 3) {
      # Multiple regions and stocks scenario
      n_stocks <- data$n_stocks
      fleet_regions <- em_info$catch_info$fleet_regions
      index_regions <- em_info$index_info$index_regions
      
      em_input <- list()
      
      for (s in 1:n_stocks) {
        em_info_new <- em_info
        # Extract relevant fleets and indices
        relevant_fleets <- which(fleet_regions == s)
        n_fleets <- length(relevant_fleets)
        em_info_new$par_inputs$catch_cv <- em_info$par_inputs$catch_cv[relevant_fleets]
        em_info_new$par_inputs$catch_Neff <- em_info$par_inputs$catch_Neff[relevant_fleets]
        em_info_new$par_inputs$use_agg_catch <- em_info$par_inputs$use_agg_catch[relevant_fleets]
        em_info_new$par_inputs$use_catch_paa <- em_info$par_inputs$use_catch_paa[relevant_fleets]
        
        relevant_indices <- which(index_regions == s)
        n_indices <- length(relevant_indices)
        em_info_new$par_inputs$index_cv <- em_info$par_inputs$index_cv[relevant_indices]
        em_info_new$par_inputs$index_Neff <- em_info$par_inputs$index_Neff[relevant_indices]
        em_info_new$par_inputs$fracyr_indices <- em_info$par_inputs$fracyr_indices[relevant_indices]
        em_info_new$par_inputs$q <- em_info$par_inputs$q[relevant_indices]
        em_info_new$par_inputs$use_indices <- em_info$par_inputs$use_indices[relevant_indices]
        em_info_new$par_inputs$use_index_paa <- em_info$par_inputs$use_index_paa[relevant_indices]
        em_info_new$par_inputs$units_indices <- em_info$par_inputs$units_indices[relevant_indices]
        em_info_new$par_inputs$units_index_paa <- em_info$par_inputs$units_index_paa[relevant_indices]
        
        # Generate basic info for current stock
        info <- generate_basic_info_em(em_info, em_years, n_stocks = 1, n_regions = 1, n_fleets = n_fleets, n_indices = n_indices)
        basic_info <- info$basic_info
        
        # Override any movement or trend information
        basic_info$move_dyn <- 0
        basic_info$onto_move <- matrix(0)
        basic_info$apply_re_trend <- 0
        basic_info$apply_mu_trend <- 0
        
        # Fill in the data from operating model simulation
        info$catch_info$agg_catch <- data$agg_catch[ind_em, relevant_fleets, drop = FALSE]
        info$index_info$agg_indices <- data$agg_indices[ind_em, relevant_indices, drop = FALSE]
        info$catch_info$catch_paa <- data$catch_paa[relevant_fleets, ind_em, , drop = FALSE]
        info$index_info$index_paa <- data$index_paa[relevant_indices, ind_em, , drop = FALSE]
        
        em_input[[s]] <- prepare_wham_input(
          basic_info = basic_info,
          selectivity = sel_em,
          M = M_em,
          NAA_re = NAA_re_em,
          move = NULL,
          age_comp = age_comp_em,
          catch_info = info$catch_info,
          index_info = info$index_info,
          F = info$F
        )
      }
    }
  } else {
    # Generic case for estimation model generation
    info <- generate_basic_info_em(em_info, em_years)
    basic_info <- info$basic_info
    
    # Fill in the data from operating model simulation
    info$catch_info$agg_catch <- data$agg_catch[ind_em, , drop = FALSE]
    info$index_info$agg_indices <- data$agg_indices[ind_em, , drop = FALSE]
    info$catch_info$catch_paa <- data$catch_paa[, ind_em, , drop = FALSE]
    info$index_info$index_paa <- data$index_paa[, ind_em, , drop = FALSE]
    
    if (em.opt$do.move) {
      
      em_input <- prepare_wham_input(
        basic_info = basic_info,
        selectivity = sel_em,
        M = M_em,
        NAA_re = NAA_re_em,
        move = move_em,
        age_comp = age_comp_em,
        catch_info = info$catch_info,
        index_info = info$index_info,
        F = info$F
      )
    } else {
      
      # Override any movement or trend information
      basic_info$move_dyn <- 0
      basic_info$onto_move <- matrix(0)
      basic_info$apply_re_trend <- 0
      basic_info$apply_mu_trend <- 0
      
      # No movement
      em_input <- prepare_wham_input(
        basic_info = basic_info,
        selectivity = sel_em,
        M = M_em,
        NAA_re = NAA_re_em,
        move = NULL,
        age_comp = age_comp_em,
        catch_info = info$catch_info,
        index_info = info$index_info,
        F = info$F
      )
    }
  }
  
  return(em_input)
}

# Helper function to generate basic info with multiple parameters
generate_basic_info_em <- function(em_info, em_years, n_stocks = NULL, n_regions = NULL, n_fleets = NULL, n_indices = NULL, n_seasons = NULL) {
  # Default to values from em_info if not provided
  if (is.null(n_stocks)) n_stocks <- em_info$par_inputs$n_stocks
  if (is.null(n_regions)) n_regions <- em_info$par_inputs$n_regions
  if (is.null(n_fleets)) n_fleets <- em_info$par_inputs$n_fleets
  if (is.null(n_indices)) n_indices <- em_info$par_inputs$n_indices
  if (is.null(n_seasons)) n_seasons <- em_info$par_inputs$n_seasons
  
  # Generate basic info using provided or default values from em_info
  basic_info <- generate_basic_info(
    n_stocks = n_stocks,
    n_regions = n_regions,
    n_indices = n_indices,
    n_fleets = n_fleets,
    n_seasons = n_seasons,
    base.years = em_years,
    life_history = em_info$par_inputs$life_history,
    n_ages = em_info$par_inputs$n_ages,
    Fbar_ages = em_info$par_inputs$Fbar_ages,
    recruit_model = em_info$par_inputs$recruit_model,
    
    # Using F_info from em_info to construct F_info list
    F_info = list(
      F.year1 = em_info$par_inputs$F.year1,
      Fhist = em_info$par_inputs$Fhist,
      Fmax = em_info$par_inputs$Fmax,
      Fmin = em_info$par_inputs$Fmin,
      change_time = em_info$par_inputs$change_time,
      user_F = em_info$par_inputs$user_F
    ),
    
    # Construct catch_info list
    catch_info = list(
      catch_cv = em_info$par_inputs$catch_cv,
      catch_Neff = em_info$par_inputs$catch_Neff,
      use_agg_catch = em_info$par_inputs$use_agg_catch, 
      use_catch_paa = em_info$par_inputs$use_catch_paa
    ),
    
    # Construct index_info list
    index_info = list(
      index_cv = em_info$par_inputs$index_cv,
      index_Neff = em_info$par_inputs$index_Neff,
      fracyr_indices = em_info$par_inputs$fracyr_indices,
      q = em_info$par_inputs$q,
      use_indices = em_info$par_inputs$use_indices, 
      use_index_paa = em_info$par_inputs$use_index_paa, 
      units_indices = em_info$par_inputs$units_indices, 
      units_index_paa = em_info$par_inputs$units_index_paa
    ),
    
    fracyr_spawn = em_info$par_inputs$fracyr_spawn,
    fracyr_seasons = em_info$par_inputs$fracyr_seasons, # Pass user-defined season fractions if available
    fleet_pointer = em_info$par_inputs$fleet_pointer, # Fleet region allocation
    index_pointer = em_info$par_inputs$index_pointer, # Index region allocation
    user_waa = em_info$par_inputs$user_waa, # Use user-defined weight-at-age if provided
    user_maturity = em_info$par_inputs$user_maturity, # Use user-defined maturity-at-age if provided
    bias.correct.process = em_info$par_inputs$bias_correct_process,
    bias.correct.observation = em_info$par_inputs$bias_correct_observation,
    bias.correct.BRPs = em_info$par_inputs$bias_correct_BRPs,
    mig_type = em_info$par_inputs$mig_type,
    XSPR_R_opt = em_info$par_inputs$XSPR_R_opt,
    
    # Add meta-population and movement parameters if specified
    move_dyn = em_info$par_inputs$move_dyn,
    onto_move = em_info$par_inputs$onto_move,
    onto_move_pars = em_info$par_inputs$onto_move_pars,
    
    # Add trend options if specified
    apply_re_trend = em_info$par_inputs$apply_re_trend,
    trend_re_rate = em_info$par_inputs$trend_re_rate,
    apply_mu_trend = em_info$par_inputs$apply_mu_trend,
    trend_mu_rate = em_info$par_inputs$trend_mu_rate,
    age_mu_devs = em_info$par_inputs$age_mu_devs
  )
  
  return(basic_info)
}
