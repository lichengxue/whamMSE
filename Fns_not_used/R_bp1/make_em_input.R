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
#'   
#' @param aggregate_catch_info List (optional). User-specified catch aggregation settings 
#'   for panmictic models using aggregate catch.
#'   \itemize{
#'     \item `$n_fleets` Integer. Number of fleets.
#'     \item `$catch_cv` Numeric vector (`n_fleets`). Catch CVs for each fleet.
#'     \item `$catch_Neff` Numeric vector (`n_fleets`). Effective sample sizes for catch.
#'     \item `$use_agg_catch` Integer vector (`n_fleets`). 0/1 flag for using aggregate catch.
#'     \item `$use_catch_paa` Integer vector (`n_fleets`). 0/1 flag for using age composition.
#'     \item `$fleet_pointer` Integer vector (`n_fleets`). Defines fleet grouping (0 = exclude).
#'     \item `$use_catch_weighted_waa` Logical. Whether to weight-at-age using catch.
#'   }
#' @param aggregate_index_info List (optional). User-specified index aggregation settings for panmictic models using aggregate indices.
#'   \itemize{
#'     \item `$n_indices` Integer. Number of indices.
#'     \item `$index_cv` Numeric vector (`n_indices`). CVs for survey indices.
#'     \item `$index_Neff` Numeric vector (`n_indices`). Effective sample sizes for indices.
#'     \item `$fracyr_indices` Numeric vector (`n_indices`). Fraction of the year for each survey.
#'     \item `$q` Numeric vector (`n_indices`). Survey catchability coefficients.
#'     \item `$use_indices` Integer vector (`n_indices`). 0/1 flag for using survey indices.
#'     \item `$use_index_paa` Integer vector (`n_indices`). 0/1 flag for using age composition.
#'     \item `$units_indices` Integer vector (`n_indices`). 1 = Biomass, 2 = Numbers.
#'     \item `$units_index_paa` Integer vector (`n_indices`). 1 = Biomass, 2 = Numbers.
#'     \item `$index_pointer` Integer vector (`n_indices`). Defines index grouping (0 = exclude).
#'     \item `$use_index_weighted_waa` Logical. Whether to weight-at-age using indices.
#'      }
#' @param filter_indices Integer (0/1) vector (optional). User-specified which indices are excluded from the assessment model 
#' @return List. A `wham` input object prepared for stock assessment.
#'
#' @export
#'
#' @seealso \code{\link{loop_through_fn}}
#' 

make_em_input <- function(om, em_info, M_em, sel_em, NAA_re_em, move_em,
                          em.opt, em_years, year.use,
                          age_comp_em,
                          aggregate_catch_info = NULL,
                          aggregate_index_info = NULL,
                          filter_indices = NULL,
                          global_waa = TRUE
                          ) {
  
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
      
      n_fleets <- ifelse(is.null(aggregate_catch_info$n_fleets), 1, aggregate_catch_info$n_fleets)
      n_indices <- ifelse(is.null(aggregate_index_info$n_indices), 1, aggregate_index_info$n_indices)
      n_stocks = n_regions = 1
      
      em_info = make_aggregate_data(om, em_info, aggregate_catch_info, aggregate_index_info, ind_em)
      
      agg_catch = em_info$par_inputs$agg_catch
      catch_paa = em_info$par_inputs$catch_paa
      agg_indices = em_info$par_inputs$agg_indices
      index_paa = em_info$par_inputs$index_paa
      
      # if (is.null(aggregate_catch_info$n_fleets)) aggregate_catch_info$n_fleets = n_fleets
      # if (is.null(aggregate_index_info$n_indices)) aggregate_index_info$n_indices = n_indices
      # if (is.null(aggregate_catch_info$fleet_pointer)) aggregate_catch_info$fleet_pointer = rep(1,om$input$data$n_fleets)
      # if (is.null(aggregate_index_info$index_pointer)) aggregate_index_info$index_pointer = rep(1,om$input$data$n_indices)
      # if (is.null(aggregate_catch_info$use_catch_weighted_waa)) aggregate_catch_info$use_catch_weighted_waa = TRUE 
      # if (is.null(aggregate_index_info$use_catch_weighted_waa)) aggregate_index_info$use_catch_weighted_waa = TRUE 
      
      # em_info = info
      # em_info <- update_em_with_basic_info(om, em_info, aggregate_catch_info, aggregate_index_info, ind_em)
      #   
      # em_info <- aggregate_em_data(data, em_info, aggregate_catch_info, aggregate_index_info, ind_em, n_fleets, n_indices)
      
      # Override any movement or trend information
      em_info$par_inputs$move_dyn <- 0
      em_info$par_inputs$onto_move <- matrix(0)
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
      
      em_input <- prepare_wham_input(
        basic_info = basic_info,
        selectivity = sel_em,
        M = M_em,
        NAA_re = NAA_re_em,
        move = NULL,
        age_comp = age_comp_em,
        catch_info = catch_info,
        index_info = index_info,
        F = F_info
      )
      
      if (!global_waa) {
        waa_info <- basic_info[grepl("waa", names(basic_info))]
        em_input <- update_waa(em_input, waa_info = waa_info)
      }
      
      # if (global_maa) {
      #   maturity <- basic_info$maturity
      #   em_input$data$mature[1,,] = colMeans(maturity[,1,])
      # } else {
      #   stop("For a panmictic assessment model, global_maa is forced to be TRUE")
      # }
    }
    
    if (em.opt$separate.em.type == 2) {
      
      # Fleets-as-areas
      n_fleets <- data$n_fleets
      n_indices <- data$n_indices
      fleet_regions <- em_info$catch_info$fleet_regions
      index_regions <- em_info$index_info$index_regions

      em_info <- filter_and_generate_em_info(em_info, fleet_regions, index_regions, filter_indices = filter_indices, em.opt)
      
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
      
      if (!is.null(filter_indices) & any(filter_indices == 0)) {
        info$index_info$agg_indices <- data$agg_indices[ind_em, idx, drop = FALSE]
        info$index_info$index_paa <- data$index_paa[idx, ind_em, , drop = FALSE]
      } else {
        info$index_info$agg_indices <- data$agg_indices[ind_em, , drop = FALSE]
        info$index_info$index_paa <- data$index_paa[, ind_em, , drop = FALSE]
      }
      
      # Override any movement or trend information
      basic_info$move_dyn <- 0
      basic_info$onto_move <- matrix(0)
      basic_info$apply_re_trend <- 0
      basic_info$apply_mu_trend <- 0
      
      if(is.null(age_comp_em)) age_comp_em = "multinomial"
        
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
      
      if (!global_waa) {
        waa_info <- basic_info[grepl("waa", names(basic_info))]
        em_input <- update_waa(em_input, waa_info = waa_info)
      }
      
      # if (global_maa) {
      #   maturity <- basic_info$maturity
      #   em_input$data$mature[1,,] = colMeans(maturity[,1,])
      # } else {
      #   stop("For a FAA assessment model, global_maa is forced to be TRUE")
      # }
      
    }
    
    if (em.opt$separate.em.type == 3) {
      # Multiple regions and stocks scenario

      # fleet_regions <- em_info$catch_info$fleet_regions
      # index_regions <- em_info$index_info$index_regions
      fleet_regions <- data$fleet_regions
      index_regions <- data$index_regions
      
      em_input <- list()
      em_info_new <- filter_and_generate_em_info(em_info, fleet_regions, index_regions, filter_indices, em.opt)
      
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
        
        relevant_indices <- which(index_regions == r)
        
        # Apply `filter_indices` (remove excluded indices) only if provided
        if (!is.null(filter_indices)) {
          relevant_indices <- relevant_indices[filter_indices[relevant_indices] != 0]
        }
        
        info$index_info$agg_indices <- data$agg_indices[ind_em, relevant_indices, drop = FALSE]
        info$index_info$index_paa <- data$index_paa[relevant_indices, ind_em, , drop = FALSE]
        
        em_input[[r]] <- prepare_wham_input(
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
        
        if (!global_waa) {
          waa_info <- basic_info[grepl("waa", names(basic_info))]
          em_input[[r]] <- update_waa(em_input[[r]], waa_info = waa_info)
        }
      }
    }
  } else {
    # Generic case for estimation model generation
    # info <- generate_basic_info_em(em_info, em_years)

    n_fleets <- data$n_fleets
    n_indices <- data$n_indices
    fleet_regions <- em_info$catch_info$fleet_regions
    index_regions <- em_info$index_info$index_regions
    
    em_info <- filter_and_generate_em_info(em_info, fleet_regions, index_regions, filter_indices = filter_indices, em.opt)
    
    if (!is.null(filter_indices) & any(filter_indices == 0)) {
      n_indices = sum(filter_indices != 0)
      idx = which(filter_indices != 0)
    } 
    
    info <- generate_basic_info_em(em_info, em_years,
                                   n_stocks = om$input$data$n_stocks, 
                                   n_regions = om$input$data$n_regions, 
                                   n_fleets = n_fleets, 
                                   n_indices = n_indices, 
                                   filter_indices = filter_indices)
    basic_info <- info$basic_info
    
    # Fill in the data from operating model simulation
    info$catch_info$agg_catch <- data$agg_catch[ind_em, , drop = FALSE]
    info$catch_info$catch_paa <- data$catch_paa[, ind_em, , drop = FALSE]
    
    if (!is.null(filter_indices) & any(filter_indices == 0)) {
      info$index_info$agg_indices <- data$agg_indices[ind_em, idx, drop = FALSE]
      info$index_info$index_paa <- data$index_paa[idx, ind_em, , drop = FALSE]
    } else {
      info$index_info$agg_indices <- data$agg_indices[ind_em, , drop = FALSE]
      info$index_info$index_paa <- data$index_paa[, ind_em, , drop = FALSE]
    }
    
    if (em.opt$do.move) {
      basic_info$NAA_where = om$input$data$NAA_where
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
      
      if (!global_waa) {
        waa_info <- basic_info[grepl("waa", names(basic_info))]
        em_input <- update_waa(em_input, waa_info = waa_info)
      }
      
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
      
      if (!global_waa) {
        waa_info <- basic_info[grepl("waa", names(basic_info))]
        em_input <- update_waa(em_input, waa_info = waa_info)
      }
    }
  }
  
  return(em_input)
}
