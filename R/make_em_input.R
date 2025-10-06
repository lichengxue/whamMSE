#' Generate Input Data for the Estimation Model
#'
#' This function generates input data for the estimation model used in 
#' management strategy evaluation. It constructs the necessary parameters 
#' and structures for the assessment model.
#'
#' @param om List. The operating model containing observed data and simulation outputs.
#' @param em_info List. A set of estimation model parameters used to define the model structure.
#' @param M_em List. Natural mortality random effects.
#' @param sel_em List. Selectivity random effects.
#' @param NAA_re_em List. Numbers-at-age random effects.
#' @param move_em List. Movement random effects.
#' @param em.opt List. Options for movement in the estimation model.
#'   \itemize{
#'     \item `$separate.em` Logical. `TRUE` indicates no global SPR, `FALSE` allows global SPR.
#'     \item `$separate.em.type` Integer (1–3). Specifies assessment model type if `separate.em = TRUE`:
#'       \itemize{
#'         \item `1` - Panmictic (spatially aggregated).
#'         \item `2` - Fleets-as-areas.
#'         \item `3` - Separate assessment models for each region (`n_regions` models).
#'       }
#'     \item `$do.move` Logical. Whether movement is included (if `separate.em = FALSE`).
#'     \item `$est.move` Logical. Whether movement rates are estimated (if `separate.em = FALSE`).
#'   }
#' @param em_years Vector. Years used in the assessment model.
#' @param year.use Integer. Number of years used in the assessment model.
#' @param age_comp_em Character. Likelihood distribution for age composition data.
#'   \itemize{
#'     \item `"multinomial"` (default)
#'     \item `"dir-mult"`
#'     \item `"dirichlet-miss0"`
#'     \item `"dirichlet-pool0"`
#'     \item `"logistic-normal-miss0"`
#'     \item `"logistic-normal-ar1-miss0"`
#'     \item `"logistic-normal-pool0"`
#'     \item `"logistic-normal-01-infl"`
#'     \item `"logistic-normal-01-infl-2par"`
#'     \item `"mvtweedie"`
#'     \item `"dir-mult-linear"`
#'   }
#' @param aggregate_catch_info List (optional). User-specified catch aggregation settings for panmictic models using aggregate catch.
#'   \itemize{
#'     \item `$n_fleets` Integer. Number of fleets.
#'     \item `$catch_cv` Numeric vector (`n_fleets`). CVs for annual aggregate catches by fleet.
#'     \item `$catch_Neff` Numeric vector (`n_fleets`). Effective sample sizes for fleet catches.
#'     \item `$use_agg_catch` Integer vector (`n_fleets`). 0/1 values flagging whether to use aggregate catches.
#'     \item `$use_catch_paa` Integer vector (`n_fleets`). 0/1 values flagging whether to use proportions at age observations.
#'     \item `$fleet_pointer` Integer vector (`n_fleets`). Defines fleet grouping (0 = exclude).
#'     \item `$use_catch_weighted_waa` Logical. Whether to use weighted weight-at-age based on fleet catches.
#'   }
#' @param aggregate_index_info List (optional). User-specified index aggregation settings for panmictic models using aggregate indices.
#'   \itemize{
#'     \item `$n_indices` Integer. Number of indices.
#'     \item `$index_cv` Numeric vector (`n_indices`). CVs for annual aggregate index catches.
#'     \item `$index_Neff` Numeric vector (`n_indices`). Effective sample sizes for survey indices.
#'     \item `$fracyr_indices` Numeric vector (`n_indices`). Fraction of the year for each survey index.
#'     \item `$q` Numeric vector (`n_indices`). Initial survey catchabilities.
#'     \item `$use_indices` Integer vector (`n_indices`). 0/1 values flagging whether to use survey indices.
#'     \item `$use_index_paa` Integer vector (`n_indices`). 0/1 values flagging whether to use proportions at age observations.
#'     \item `$units_indices` Integer vector (`n_indices`). 1/2 values flagging whether aggregate observations are biomass (1) or numbers (2).
#'     \item `$units_index_paa` Integer vector (`n_indices`). 1/2 values flagging whether composition observations are biomass (1) or numbers (2).
#'     \item `$index_pointer` Integer vector (`n_indices`). Defines index grouping (0 = exclude).
#'     \item `$use_index_weighted_waa` Logical. Whether to use weighted weight-at-age based on survey catches.
#'      }
#' @param aggregate_weights_info List (optional). Specifies how to compute weighted averages for
#' weight-at-age (`waa`) and maturity-at-age during data aggregation (For panmictic and fleets-as-areas models only).
#' Used for aggregating across fleets or indices to form total/stock/regional summaries.
#' \itemize{
#'   \item `$ssb_waa_weights` List. Settings for weighting weight-at-age used for spawning stock biomass.
#'     \itemize{
#'       \item `$fleet` Logical. Whether to use fleet-specific weights.
#'       \item `$index` Logical. Whether to use index-specific weights.
#'       \item `$pointer` Integer. Index pointing to the selected fleet or index group (e.g., 1 means the first valid fleet group).
#'     }
#'   \item `$maturity_weights` List. Settings for weighting maturity-at-age used for spawning stock biomass.
#'     \itemize{
#'       \item `$fleet` Logical. Whether to use fleet-specific weights.
#'       \item `$index` Logical. Whether to use index-specific weights.
#'       \item `$pointer` Integer. Index pointing to the selected fleet or index group.
#'     }
#' }
#' @param filter_indices Integer (0/1) vector (optional). User-specified which indices are excluded from the assessment model 
#' @param reduce_region_info List (optional). Remove specific regions from the assessment model. If `NULL`, no modifications are applied.
#'   The expected components include:
#'   \itemize{
#'     \item `$remove_regions` Integer vector (`n_regions`). 0/1 values flagging whether to include regions.
#'     \item `$reassign` Numeric. Specify reassignment for surveys from the "removed" regions to "remaining" regions.
#'     \item `$NAA_where_em`Integer array (`n_stocks × n_regions × n_ages`) 0/1 values flagging whether the stock from a certain age can be present in a specific region.
#'     \item `$sel_em` Selectivity configuration for the "reduced" assessment model.
#'     \item `$M_em` Natural mortality configuration for the "reduced" assessment model.
#'     \item `$NAA_re_em` Numbers-at-age configuration for the "reduced" assessment model.
#'     \item `$catchability_em` Cathability configuration for the "reduced" assessment model.
#'     \item `$move_em` Movement configuration for the "reduced" assessment model.
#'     \item `$ecov_em` Environmental covariate configuration for the "reduced" assessment model.
#'     \item `$onto_move_list` List of ontogenetic movement information for the "reduced" assessment model:
#'       \itemize{
#'         \item `$onto_move` (array, dimension: `n_stocks × n_regions × (n_regions-1)`)  
#'           Age-specific movement type.
#'         \item `$onto_move_pars` (array, dimension: `n_stocks × n_regions × (n_regions-1) × 4`)  
#'           Parameters controlling age-specific movement patterns.
#'         \item `$age_mu_devs` (array, dimension: `n_stocks × n_regions × (n_regions-1) × n_ages`)  
#'           User-specified deviations in mean movement rates across ages.
#'       }
#'   }
#' @param update_catch_info List (optional). Update catch CV and Neff in the EM.
#'   The expected components include:
#'   \itemize{
#'     \item `$agg_catch_sigma` Matrix. Either full (n_years x n_fleets) or subset (length(ind_em) x n_fleets).
#'     \item `$catch_Neff` Matrix. Same dimension as `agg_catch_sigma`.
#'   }
#' @param update_index_info List (optional). Update index CV and Neff in the EM.
#'   The expected components include:
#'   \itemize{
#'     \item `$agg_index_sigma` Matrix. Either full (n_years x n_indices) or subset (length(ind_em) x n_indices).
#'     \item `$index_Neff` Matrix. Same dimension as `agg_index_sigma`.
#'   }
#' @return List. A `wham` input object prepared for stock assessment.
#' 
#' @export
#'
#' @seealso \code{\link{loop_through_fn}}
#' 

make_em_input <- function(om, 
                          em_info, 
                          M_em, 
                          sel_em, 
                          NAA_re_em, 
                          move_em, 
                          catchability_em,
                          ecov_em,
                          em.opt, 
                          em_years, 
                          year.use,
                          age_comp_em,
                          aggregate_catch_info = NULL,
                          aggregate_index_info = NULL,
                          aggregate_weights_info = NULL,
                          filter_indices = NULL,
                          reduce_region_info = NULL,
                          update_catch_info = NULL,
                          update_index_info = NULL,
                          ecov_em_opts = NULL) {
  
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
  
  if (em.opt$separate.em) {    # Non-spatial or Spatially implicit  

    if (em.opt$separate.em.type == 1) {
      
      n_fleets <- ifelse(is.null(aggregate_catch_info$n_fleets), 1, aggregate_catch_info$n_fleets)
      n_indices <- ifelse(is.null(aggregate_index_info$n_indices), 1, aggregate_index_info$n_indices)
      n_stocks = n_regions = 1
      
      em_info = make_aggregate_data(om, em_info, ind_em, aggregate_catch_info, aggregate_index_info, aggregate_weights_info)
      
      agg_catch = em_info$par_inputs$agg_catch
      catch_paa = em_info$par_inputs$catch_paa
      agg_indices = em_info$par_inputs$agg_indices
      index_paa = em_info$par_inputs$index_paa
      
      # Override any movement or trend information
      em_info$par_inputs$move_dyn <- 0
      em_info$par_inputs$onto_move <-  matrix(0)
      em_info$par_inputs$apply_re_trend <- 0
      em_info$par_inputs$apply_mu_trend <- 0
      
      info <- generate_basic_info_em(em_info, em_years, n_stocks = 1, n_regions = 1, n_fleets = n_fleets, n_indices = n_indices)
      
      basic_info <- info$basic_info
      catch_info <- info$catch_info
      index_info <- info$index_info
      F_info     <- info$F
      
      # Fill in the data from the operating model simulation
      catch_info$agg_catch <- agg_catch
      catch_info$catch_paa <- catch_paa
      index_info$agg_indices <- agg_indices
      index_info$index_paa <- index_paa
      
      # Ecov Effect 
      if(!is.null(ecov_em)) {
        ecov_em_new <- ecov_em
        ecov_em_new$year <- ecov_em_new$year
        ecov_mean <- om$input$data$Ecov_obs
        ecov_em_new$mean <- ecov_mean
        if (any(ecov_em_new$logsigma %in% c("est_1", "est_re"))) {
          ecov_em_new$logsigma = ecov_em_new$logsigma
        } else {
          ecov_em_new$logsigma <- ecov_em_new$logsigma
        }
        ecov_em_new$use_obs <- ecov_em_new$use_obs
      } else {
        ecov_em_new <- NULL
      }
      
      em_input <- prepare_wham_input(
        basic_info = basic_info,
        selectivity = sel_em,
        M = M_em,
        NAA_re = NAA_re_em,
        move = NULL,
        catchability = catchability_em,
        ecov = ecov_em_new,
        age_comp = age_comp_em,
        catch_info = catch_info,
        index_info = index_info,
        F = F_info
      )
      
      waa_info <- info$par_inputs$user_waa
      em_input <- update_waa(em_input, waa_info = waa_info)
      
    }
    
    if (em.opt$separate.em.type == 2) {
      
      # Fleets-as-areas
      n_fleets <- data$n_fleets
      n_indices <- data$n_indices
      if(!is.null(filter_indices) && length(filter_indices) != n_indices) stop("Length of filter_indices must = n_indices!")
      fleet_regions <- em_info$catch_info$fleet_regions
      index_regions <- em_info$index_info$index_regions

      em_info <- filter_and_generate_em_info(em_info, 
                                             em.opt = em.opt, 
                                             ind_em, 
                                             fleet_regions, 
                                             index_regions, 
                                             filter_indices = filter_indices, 
                                             aggregate_weights_info = aggregate_weights_info)
      
      if (!is.null(filter_indices) & any(filter_indices == 0)) {
        n_indices = sum(filter_indices != 0)
        idx = which(filter_indices != 0)
      } 
      
      info <- generate_basic_info_em(em_info, em_years, n_stocks = 1, n_regions = 1, n_fleets = n_fleets, n_indices = n_indices, filter_indices = filter_indices)
      
      basic_info <- info$basic_info

      info$catch_info$fleet_regions[] = 1
      info$index_info$index_regions[] = 1
      
      # Fill in the data from the operating model simulation
      info$catch_info$agg_catch <- data$agg_catch[ind_em, , drop = FALSE]
      info$catch_info$catch_paa <- data$catch_paa[, ind_em, , drop = FALSE]
      info$catch_info$use_agg_catch <- data$use_agg_catch[ind_em, , drop = FALSE]
      info$catch_info$use_catch_paa <- data$use_catch_paa[ind_em, , drop = FALSE]
      
      # if (!is.null(filter_indices) & any(filter_indices == 0)) {
      #   info$index_info$agg_indices <- data$agg_indices[ind_em, idx, drop = FALSE]
      #   info$index_info$index_paa <- data$index_paa[idx, ind_em, , drop = FALSE]
      # } else {
      #   info$index_info$agg_indices <- data$agg_indices[ind_em, , drop = FALSE]
      #   info$index_info$index_paa <- data$index_paa[, ind_em, , drop = FALSE]
      # }
      
      if (!is.null(filter_indices) & any(filter_indices == 0)) {
        info$index_info$agg_indices <- data$agg_indices[ind_em, idx, drop = FALSE]
        info$index_info$index_paa <- data$index_paa[idx, ind_em, , drop = FALSE]
        info$index_info$use_indices <- data$use_indices[ind_em, idx, drop = FALSE]
        info$index_info$use_index_paa <- data$use_index_paa[ind_em, idx, drop = FALSE]
      } else {
        info$index_info$agg_indices <- data$agg_indices[ind_em, , drop = FALSE]
        info$index_info$index_paa <- data$index_paa[, ind_em, , drop = FALSE]
        info$index_info$use_indices <- data$use_indices[ind_em, , drop = FALSE]
        info$index_info$use_index_paa <- data$use_index_paa[ind_em, ,  drop = FALSE]
      }
      
      # Override any movement or trend information
      basic_info$move_dyn <- 0
      basic_info$onto_move <- matrix(0)
      basic_info$apply_re_trend <- 0
      basic_info$apply_mu_trend <- 0
      
      if(is.null(age_comp_em)) age_comp_em = "multinomial"
      
      # Ecov Effect 
      if(!is.null(ecov_em)) {
        ecov_em_new <- ecov_em
        ecov_em_new$year <- ecov_em_new$year
        ecov_mean <- om$input$data$Ecov_obs
        ecov_em_new$mean <- ecov_mean
        if (any(ecov_em_new$logsigma %in% c("est_1", "est_re"))) {
          ecov_em_new$logsigma = ecov_em_new$logsigma
        } else {
          ecov_em_new$logsigma <- ecov_em_new$logsigma
        }
        ecov_em_new$use_obs <- ecov_em_new$use_obs
      } else {
        ecov_em_new <- NULL
      }
      
      em_input <- prepare_wham_input(
        basic_info = basic_info,
        selectivity = sel_em,
        M = M_em,
        NAA_re = NAA_re_em,
        move = NULL,
        catchability = catchability_em,
        ecov = ecov_em_new,
        age_comp = age_comp_em,
        catch_info = info$catch_info,
        index_info = info$index_info,
        F = info$F
      )
      
      waa_info <- info$par_inputs$user_waa
      em_input <- update_waa(em_input, waa_info = waa_info)
      
    }
    
    if (em.opt$separate.em.type == 3) {
      # Multiple regions and stocks scenario

      # fleet_regions <- em_info$catch_info$fleet_regions
      # index_regions <- em_info$index_info$index_regions
      fleet_regions <- data$fleet_regions
      index_regions <- data$index_regions
      
      em_input <- list()
      em_info_new <- filter_and_generate_em_info(em_info, 
                                                 em.opt = em.opt, 
                                                 ind_em, 
                                                 fleet_regions, 
                                                 index_regions, 
                                                 filter_indices)
      
      for (r in 1:data$n_regions) {
        
        # Generate basic info for current stock
        info <- generate_basic_info_em(em_info_new[[r]], em_years, n_stocks = 1, n_regions = 1, 
                                       n_fleets = em_info_new$par_inputs$n_fleets, 
                                       n_indices = em_info_new$par_inputs$n_indices, filter_indices)
        basic_info <- info$basic_info
        
        # Override any movement or trend information
        basic_info$move_dyn <- 0
        basic_info$onto_move <- matrix(0)
        basic_info$apply_re_trend <- 0
        basic_info$apply_mu_trend <- 0
        
        relevant_fleets <- which(fleet_regions == r)
        
        # Fill in the data from operating model simulation
        info$catch_info$agg_catch <- data$agg_catch[ind_em, relevant_fleets, drop = FALSE]
        info$catch_info$catch_paa <- data$catch_paa[relevant_fleets, ind_em, , drop = FALSE]
        info$catch_info$use_agg_catch <- data$use_agg_catch[ind_em, relevant_fleets, drop = FALSE]
        info$catch_info$use_catch_paa <- data$use_catch_paa[ind_em, relevant_fleets, drop = FALSE]
        
        relevant_indices <- which(index_regions == r)
        
        # Apply `filter_indices` (remove excluded indices) only if provided
        if (!is.null(filter_indices)) {
          relevant_indices <- relevant_indices[filter_indices[relevant_indices] != 0]
        }
        
        info$index_info$agg_indices <- data$agg_indices[ind_em, relevant_indices, drop = FALSE]
        info$index_info$index_paa <- data$index_paa[relevant_indices, ind_em, , drop = FALSE]
        info$index_info$use_indices <- data$use_indices[ind_em, relevant_indices, drop = FALSE]
        info$index_info$use_index_paa <- data$use_index_paa[ind_em, relevant_indices, drop = FALSE]
        
        # Ecov Effect 
        if(!is.null(ecov_em)) {
          ecov_em_new <- ecov_em
          ecov_em_new$year <- ecov_em_new$year
          ecov_mean <- om$input$data$Ecov_obs
          ecov_em_new$mean <- ecov_mean
          if (any(ecov_em_new$logsigma %in% c("est_1", "est_re"))) {
            ecov_em_new$logsigma = ecov_em_new$logsigma
          } else {
            ecov_em_new$logsigma <- ecov_em_new$logsigma
          }
          ecov_em_new$use_obs <- ecov_em_new$use_obs
        } else {
          ecov_em_new <- NULL
        }
        
        em_input[[r]] <- prepare_wham_input(
          basic_info = basic_info,
          selectivity = sel_em,
          M = M_em,
          NAA_re = NAA_re_em,
          move = NULL,
          catchability = catchability_em,
          ecov = ecov_em_new,
          age_comp = age_comp_em,
          catch_info = info$catch_info,
          index_info = info$index_info,
          F = info$F
        )
        
        waa_info <- basic_info[grepl("waa", names(basic_info))]
        em_input[[r]] <- update_waa(em_input[[r]], waa_info = waa_info)
      }
    }
  } 
  
  ##############################################
  ############# Spatially explicit #############
  ##############################################
  
  if (!em.opt$separate.em) { 
    
    n_fleets <- data$n_fleets
    n_indices <- data$n_indices
    fleet_regions <- em_info$catch_info$fleet_regions
    index_regions <- em_info$index_info$index_regions
    
    remove_regions <- reduce_region_info$remove_regions
    
    em_info <- filter_and_generate_em_info(em_info, 
                                           em.opt = em.opt, 
                                           ind_em, 
                                           fleet_regions, 
                                           index_regions, 
                                           filter_indices = filter_indices,
                                           reduce_region_info = reduce_region_info)
    
    if (!is.null(filter_indices) & any(filter_indices == 0)) {
      n_indices = sum(filter_indices != 0)
      idx = which(filter_indices != 0)
    } 
    
    # If remove_regions = NULL
    if (is.null(remove_regions)) {
      
      n_stocks = om$input$data$n_stocks
      n_regions = om$input$data$n_regions
      
      info <- generate_basic_info_em(em_info, 
                                     em_years,
                                     n_stocks = n_stocks, 
                                     n_regions = n_regions, 
                                     n_fleets = n_fleets, 
                                     n_indices = n_indices, 
                                     filter_indices = filter_indices)
      
      basic_info <- info$basic_info
      
      info$catch_info$agg_catch <- data$agg_catch[ind_em, , drop = FALSE]
      info$catch_info$catch_paa <- data$catch_paa[, ind_em, , drop = FALSE]
      info$catch_info$use_agg_catch <- data$use_agg_catch[ind_em, , drop = FALSE]
      info$catch_info$use_catch_paa <- data$use_catch_paa[ind_em, , drop = FALSE]
      
      # if (!is.null(filter_indices) & any(filter_indices == 0)) {
      #   info$index_info$agg_indices <- data$agg_indices[ind_em, idx, drop = FALSE]
      #   info$index_info$index_paa <- data$index_paa[idx, ind_em, , drop = FALSE]
      # } else {
      #   info$index_info$agg_indices <- data$agg_indices[ind_em, , drop = FALSE]
      #   info$index_info$index_paa <- data$index_paa[, ind_em, , drop = FALSE]
      # }
      
      if (!is.null(filter_indices) & any(filter_indices == 0)) {
        info$index_info$agg_indices <- data$agg_indices[ind_em, idx, drop = FALSE]
        info$index_info$index_paa <- data$index_paa[idx, ind_em, , drop = FALSE]
        info$index_info$use_indices <- data$use_indices[ind_em, idx, drop = FALSE]
        info$index_info$use_index_paa <- data$use_index_paa[ind_em, idx, drop = FALSE]
      } else {
        info$index_info$agg_indices <- data$agg_indices[ind_em, , drop = FALSE]
        info$index_info$index_paa <- data$index_paa[, ind_em, , drop = FALSE]
        info$index_info$use_indices <- data$use_indices[ind_em, , drop = FALSE]
        info$index_info$use_index_paa <- data$use_index_paa[ind_em, ,  drop = FALSE]
      }
      
      # Ecov Effect 
      if(!is.null(ecov_em)) {
        ecov_em_new <- ecov_em
        ecov_em_new$year <- ecov_em_new$year
        ecov_mean <- om$input$data$Ecov_obs
        ecov_em_new$mean <- ecov_mean
        if(!is.null(ecov_em_opts) && ecov_em_opts$use_ecov_em) ecov_em_new$mean[ecov_em_opts$period,] <- ecov_em[ecov_em_opts$period,]
        if (any(ecov_em_new$logsigma %in% c("est_1", "est_re"))) {
          ecov_em_new$logsigma = ecov_em_new$logsigma
        } else {
          ecov_em_new$logsigma <- ecov_em_new$logsigma
        }
        ecov_em_new$use_obs <- ecov_em_new$use_obs
      } else {
        ecov_em_new <- NULL
      }
      
      if (em.opt$do.move) {
        basic_info$NAA_where = om$input$data$NAA_where
        em_input <- prepare_wham_input(
          basic_info = basic_info,
          selectivity = sel_em,
          M = M_em,
          NAA_re = NAA_re_em,
          move = move_em,
          catchability = catchability_em,
          ecov = ecov_em_new,
          age_comp = age_comp_em,
          catch_info = info$catch_info,
          index_info = info$index_info,
          F = info$F
        )
        
        waa_info <- info$par_inputs$user_waa
        em_input <- update_waa(em_input, waa_info = waa_info)
        
        if(!is.null(update_catch_info)) {
          agg_catch_sigma = update_catch_info$agg_catch_sigma
          catch_Neff = update_catch_info$catch_Neff
          remove_agg = update_catch_info$remove_agg
          remove_agg_pointer = update_catch_info$remove_agg_pointer 
          remove_agg_years = update_catch_info$remove_agg_years
          remove_paa = update_catch_info$remove_paa
          remove_paa_pointer = update_catch_info$remove_paa_pointer 
          remove_paa_years = update_catch_info$remove_paa_years
          em_input = update_input_catch_info(input = em_input, 
                                             agg_catch_sigma = agg_catch_sigma, 
                                             catch_Neff = catch_Neff,
                                             remove_agg = remove_agg,
                                             remove_agg_pointer = remove_agg_pointer,
                                             remove_agg_years = remove_agg_years,
                                             remove_paa = remove_paa,
                                             remove_paa_pointer = remove_paa_pointer, 
                                             remove_paa_years = remove_paa_years,
                                             ind_em = ind_em)
        }
        
        if(!is.null(update_index_info)) {
          agg_index_sigma = update_index_info$agg_index_sigma
          index_Neff = update_index_info$index_Neff
          remove_agg = update_index_info$remove_agg
          remove_agg_pointer = update_index_info$remove_agg_pointer 
          remove_agg_years = update_index_info$remove_agg_years
          remove_paa = update_index_info$remove_paa
          remove_paa_pointer = update_index_info$remove_paa_pointer 
          remove_paa_years = update_index_info$remove_paa_years
          em_input = update_input_index_info(input = em_input, 
                                             agg_index_sigma = agg_index_sigma, 
                                             index_Neff = index_Neff, 
                                             remove_agg = remove_agg,
                                             remove_agg_pointer = remove_agg_pointer,
                                             remove_agg_years = remove_agg_years,
                                             remove_paa = remove_paa,
                                             remove_paa_pointer = remove_paa_pointer, 
                                             remove_paa_years = remove_paa_years,
                                             ind_em = ind_em)
        }
        
      } else {
        
        # No movement
        # Override any movement or trend information
        basic_info$move_dyn <- 0
        basic_info$onto_move <- matrix(0)
        basic_info$apply_re_trend <- 0
        basic_info$apply_mu_trend <- 0
        
        basic_info$NAA_where = NULL
        
        # No movement
        em_input <- prepare_wham_input(
          basic_info = basic_info,
          selectivity = sel_em,
          M = M_em,
          NAA_re = NAA_re_em,
          move = NULL,
          catchability = catchability_em,
          ecov = ecov_em_new,
          age_comp = age_comp_em,
          catch_info = info$catch_info,
          index_info = info$index_info,
          F = info$F
        )
        
        waa_info <- info$par_inputs$user_waa
        em_input <- update_waa(em_input, waa_info = waa_info)
        
        if(!is.null(update_catch_info)) {
          agg_catch_sigma = update_catch_info$agg_catch_sigma
          catch_Neff = update_catch_info$catch_Neff
          remove_agg = update_catch_info$remove_agg
          remove_agg_pointer = update_catch_info$remove_agg_pointer 
          remove_agg_years = update_catch_info$remove_agg_years
          remove_paa = update_catch_info$remove_paa
          remove_paa_pointer = update_catch_info$remove_paa_pointer 
          remove_paa_years = update_catch_info$remove_paa_years
          em_input = update_input_catch_info(input = em_input, 
                                             agg_catch_sigma = agg_catch_sigma, 
                                             catch_Neff = catch_Neff,
                                             remove_agg = remove_agg,
                                             remove_agg_pointer = remove_agg_pointer,
                                             remove_agg_years = remove_agg_years,
                                             remove_paa = remove_paa,
                                             remove_paa_pointer = remove_paa_pointer, 
                                             remove_paa_years = remove_paa_years,
                                             ind_em = ind_em)
        }
        
        if(!is.null(update_index_info)) {
          agg_index_sigma = update_index_info$agg_index_sigma
          index_Neff = update_index_info$index_Neff
          remove_agg = update_index_info$remove_agg
          remove_agg_pointer = update_index_info$remove_agg_pointer 
          remove_agg_years = update_index_info$remove_agg_years
          remove_paa = update_index_info$remove_paa
          remove_paa_pointer = update_index_info$remove_paa_pointer 
          remove_paa_years = update_index_info$remove_paa_years
          em_input = update_input_index_info(input = em_input, 
                                             agg_index_sigma = agg_index_sigma, 
                                             index_Neff = index_Neff, 
                                             remove_agg = remove_agg,
                                             remove_agg_pointer = remove_agg_pointer,
                                             remove_agg_years = remove_agg_years,
                                             remove_paa = remove_paa,
                                             remove_paa_pointer = remove_paa_pointer, 
                                             remove_paa_years = remove_paa_years,
                                             ind_em = ind_em)
        }
        
      }
    } 
    
    # If there is "remove_regions" 
    if (!is.null(remove_regions)) {
      n_stocks = if(!is.null(em_info$par_inputs$n_stocks))em_info$par_inputs$n_stocks
      n_regions = if(!is.null(em_info$par_inputs$n_regions))em_info$par_inputs$n_regions
      n_fleets = if(!is.null(em_info$par_inputs$n_fleets))em_info$par_inputs$n_fleets
      n_indices = if(!is.null(em_info$par_inputs$n_indices))em_info$par_inputs$n_indices
      
      info <- generate_basic_info_em(em_info = em_info, 
                                     em_years = em_years,
                                     n_stocks = n_stocks, 
                                     n_regions = n_regions, 
                                     n_fleets = n_fleets, 
                                     n_indices = n_indices, 
                                     filter_indices = filter_indices)
                                            
      basic_info <- info$basic_info
      id_fleets = info$fleets_to_remove
      id_indices = info$indices_to_remove
      if(id_indices == 0 || is.null(id_indices)) id_indices = numeric(0)
      
      if(length(id_fleets) > 0) {
        info$catch_info$agg_catch <- data$agg_catch[ind_em, -id_fleets, drop = FALSE]
        info$catch_info$catch_paa <- data$catch_paa[-id_fleets, ind_em, , drop = FALSE]
        info$catch_info$use_agg_catch <- info$catch_info$use_agg_catch[ind_em, -id_fleets, drop = FALSE]
        info$catch_info$use_catch_paa <- info$catch_info$use_catch_paa[ind_em, -id_fleets, drop = FALSE]
      } else {
        info$catch_info$agg_catch <- data$agg_catch[ind_em, , drop = FALSE]
        info$catch_info$catch_paa <- data$catch_paa[, ind_em, , drop = FALSE]
        info$catch_info$use_agg_catch <- info$catch_info$use_agg_catch[ind_em, , drop = FALSE]
        info$catch_info$use_catch_paa <- info$catch_info$use_catch_paa[ind_em, , drop = FALSE]
      }

      if (!is.null(filter_indices) & any(filter_indices == 0)) {
        info$index_info$agg_indices <- data$agg_indices[ind_em, idx, drop = FALSE]
        info$index_info$index_paa <- data$index_paa[idx, ind_em, , drop = FALSE]
        info$index_info$use_indices <- data$use_indices[ind_em, idx, drop = FALSE]
        info$index_info$use_index_paa <- data$use_index_paa[ind_em, idx, drop = FALSE]
      } else {
        info$index_info$agg_indices <- data$agg_indices[ind_em, , drop = FALSE]
        info$index_info$index_paa <- data$index_paa[, ind_em, , drop = FALSE]
        info$index_info$use_indices <- data$use_indices[ind_em, , drop = FALSE]
        info$index_info$use_index_paa <- data$use_index_paa[ind_em, , drop = FALSE]
      }
      
      if(length(id_indices) > 0) {
        info$index_info$agg_indices <- data$agg_indices[ind_em, -id_indices, drop = FALSE]
        info$index_info$index_paa <- data$index_paa[-id_indices, ind_em, , drop = FALSE]
        info$index_info$use_indices <- data$use_indices[ind_em, -id_indices, drop = FALSE]
        info$index_info$use_index_paa <- data$use_index_paa[ind_em, -id_indices, drop = FALSE]
      } else {
        info$index_info$agg_indices <- data$agg_indices[ind_em, , drop = FALSE]
        info$index_info$index_paa <- data$index_paa[, ind_em, , drop = FALSE]
        info$index_info$use_indices <- data$use_indices[ind_em, , drop = FALSE]
        info$index_info$use_index_paa <- data$use_index_paa[ind_em, , drop = FALSE]
      }
      
      if(is.null(reduce_region_info)) Stop("Users must prepare a list of new model configuration (NAA_where, sel_em, M_em, NAA_re_em, move_em, onto_move_list) if some areas are dropped from the model!")
      
      basic_info$NAA_where = reduce_region_info$NAA_where_em
      sel_em               = reduce_region_info$sel_em
      M_em                 = reduce_region_info$M_em
      NAA_re_em            = reduce_region_info$NAA_re_em
      move_em              = reduce_region_info$move_em
      catchability_em      = reduce_region_info$catchability_em
      ecov_em              = reduce_region_info$ecov_em
      onto_move_list       = reduce_region_info$onto_move_list
      
      if(n_regions == 1) move_em = NULL # This must be NULL when n_regions = 1
      if(n_regions == 1) basic_info$NAA_where = NULL # This must be NULL when n_regions = 1
      if(n_regions == 1) {
        basic_info$move_dyn <- 0
        basic_info$onto_move <- matrix(0)
        basic_info$apply_re_trend <- 0
        basic_info$apply_mu_trend <- 0
      }  
      
      if (em.opt$do.move) {
        
        # Ecov Effect 
        if(!is.null(ecov_em)) {
          ecov_em_new <- ecov_em
          ecov_em_new$year <- ecov_em_new$year
          ecov_mean <- om$input$data$Ecov_obs
          ecov_em_new$mean <- ecov_mean
          if(!is.null(ecov_em_opts) && ecov_em_opts$use_ecov_em) ecov_em_new$mean[ecov_em_opts$period,] <- ecov_em[ecov_em_opts$period,]
          if (any(ecov_em_new$logsigma %in% c("est_1", "est_re"))) {
            ecov_em_new$logsigma = ecov_em_new$logsigma
          } else {
            ecov_em_new$logsigma <- ecov_em_new$logsigma
          }
          ecov_em_new$use_obs <- ecov_em_new$use_obs
        } else {
          ecov_em_new <- NULL
        }
        
        em_input <- prepare_wham_input(
          basic_info = basic_info,
          selectivity = sel_em,
          M = M_em,
          NAA_re = NAA_re_em,
          move = move_em,
          catchability = catchability_em,
          ecov = ecov_em_new,
          age_comp = age_comp_em,
          catch_info = info$catch_info,
          index_info = info$index_info,
          F = info$F
        )
        
        em_input <- update_waa(em_input, waa_info = em_info$par_inputs$user_waa)
        
        if(!is.null(update_catch_info)) {
          agg_catch_sigma = update_catch_info$agg_catch_sigma
          catch_Neff = update_catch_info$catch_Neff
          remove_agg = update_catch_info$remove_agg
          remove_agg_pointer = update_catch_info$remove_agg_pointer 
          remove_agg_years = update_catch_info$remove_agg_years
          remove_paa = update_catch_info$remove_paa
          remove_paa_pointer = update_catch_info$remove_paa_pointer 
          remove_paa_years = update_catch_info$remove_paa_years
          em_input = update_input_catch_info(input = em_input, 
                                             agg_catch_sigma = agg_catch_sigma, 
                                             catch_Neff = catch_Neff,
                                             remove_agg = remove_agg,
                                             remove_agg_pointer = remove_agg_pointer,
                                             remove_agg_years = remove_agg_years,
                                             remove_paa = remove_paa,
                                             remove_paa_pointer = remove_paa_pointer, 
                                             remove_paa_years = remove_paa_years,
                                             ind_em = ind_em)
        }
        
        if(!is.null(update_index_info)) {
          agg_index_sigma = update_index_info$agg_index_sigma
          index_Neff = update_index_info$index_Neff
          remove_agg = update_index_info$remove_agg
          remove_agg_pointer = update_index_info$remove_agg_pointer 
          remove_agg_years = update_index_info$remove_agg_years
          remove_paa = update_index_info$remove_paa
          remove_paa_pointer = update_index_info$remove_paa_pointer 
          remove_paa_years = update_index_info$remove_paa_years
          em_input = update_input_index_info(input = em_input, 
                                             agg_index_sigma = agg_index_sigma, 
                                             index_Neff = index_Neff, 
                                             remove_agg = remove_agg,
                                             remove_agg_pointer = remove_agg_pointer,
                                             remove_agg_years = remove_agg_years,
                                             remove_paa = remove_paa,
                                             remove_paa_pointer = remove_paa_pointer, 
                                             remove_paa_years = remove_paa_years,
                                             ind_em = ind_em)
        }
        
      } else {
        
        # No movement
        # Override any movement or trend information
        basic_info$move_dyn <- 0
        basic_info$onto_move <- matrix(0)
        basic_info$apply_re_trend <- 0
        basic_info$apply_mu_trend <- 0
        
        basic_info$NAA_where = NULL
        
        if(!is.null(ecov_em)) {
          ecov_em_new <- ecov_em
          ecov_em_new$year <- ecov_em_new$year
          ecov_mean <- om$input$data$Ecov_obs
          ecov_em_new$mean <- ecov_mean
          if(!is.null(ecov_em_opts) && ecov_em_opts$use_ecov_em) ecov_em_new$mean[ecov_em_opts$period,] <- ecov_em[ecov_em_opts$period,]
          if (any(ecov_em_new$logsigma %in% c("est_1", "est_re"))) {
            ecov_em_new$logsigma = ecov_em_new$logsigma
          } else {
            ecov_em_new$logsigma <- ecov_em_new$logsigma
          }
          ecov_em_new$use_obs <- ecov_em_new$use_obs
        } else {
          ecov_em_new <- NULL
        }
        
        # No movement
        em_input <- prepare_wham_input(
          basic_info = basic_info,
          selectivity = sel_em,
          M = M_em,
          NAA_re = NAA_re_em,
          move = NULL,
          catchability = catchability_em,
          ecov = ecov_em_new,
          age_comp = age_comp_em,
          catch_info = info$catch_info,
          index_info = info$index_info,
          F = info$F
        )
        
        em_input <- update_waa(em_input, waa_info = em_info$par_inputs$user_waa)
        
        if(!is.null(update_catch_info)) {
          agg_catch_sigma = update_catch_info$agg_catch_sigma
          catch_Neff = update_catch_info$catch_Neff
          remove_agg = update_catch_info$remove_agg
          remove_agg_pointer = update_catch_info$remove_agg_pointer 
          remove_agg_years = update_catch_info$remove_agg_years
          remove_paa = update_catch_info$remove_paa
          remove_paa_pointer = update_catch_info$remove_paa_pointer 
          remove_paa_years = update_catch_info$remove_paa_years
          em_input = update_input_catch_info(input = em_input, 
                                             agg_catch_sigma = agg_catch_sigma, 
                                             catch_Neff = catch_Neff,
                                             remove_agg = remove_agg,
                                             remove_agg_pointer = remove_agg_pointer,
                                             remove_agg_years = remove_agg_years,
                                             remove_paa = remove_paa,
                                             remove_paa_pointer = remove_paa_pointer, 
                                             remove_paa_years = remove_paa_years,
                                             ind_em = ind_em)
        }
        
        if(!is.null(update_index_info)) {
          agg_index_sigma = update_index_info$agg_index_sigma
          index_Neff = update_index_info$index_Neff
          remove_agg = update_index_info$remove_agg
          remove_agg_pointer = update_index_info$remove_agg_pointer 
          remove_agg_years = update_index_info$remove_agg_years
          remove_paa = update_index_info$remove_paa
          remove_paa_pointer = update_index_info$remove_paa_pointer 
          remove_paa_years = update_index_info$remove_paa_years
          em_input = update_input_index_info(input = em_input, 
                                             agg_index_sigma = agg_index_sigma, 
                                             index_Neff = index_Neff, 
                                             remove_agg = remove_agg,
                                             remove_agg_pointer = remove_agg_pointer,
                                             remove_agg_years = remove_agg_years,
                                             remove_paa = remove_paa,
                                             remove_paa_pointer = remove_paa_pointer, 
                                             remove_paa_years = remove_paa_years,
                                             ind_em = ind_em)
        }
        
      }
    }
    
  }
  
  return(em_input)
}
