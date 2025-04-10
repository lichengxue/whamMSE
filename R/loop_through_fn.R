#' Perform Management Strategy Evaluation (MSE)
#' 
#' A wrapper function to iterate through the feedback period in a management strategy evaluation (MSE), 
#' updating the operating model (OM), refitting the assessment model (AM), and generating catch advice.
#' 
#' @param om An operating model that contains pseudo data.
#' @param em_info A list that contains information used to generate the assessment model (default = NULL).
#' @param random Character vector of processes treated as random effects in the operating model (default = "log_NAA").
#' @param M_em Configuration of natural mortality in the assessment model.
#' @param sel_em Configuration of selectivity in the assessment model.
#' @param NAA_re_em Configuration of numbers-at-age (NAA) in the assessment model.
#' @param move_em Configuration of movement in the assessment model.
#' @param age_comp_em Character. The likelihood distribution of age composition data in the assessment model. Options include:
#'   \itemize{
#'     \item \code{"multinomial"} (default)
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
#' @param em.opt List of options for the assessment model, including:
#'   \describe{
#'     \item{\code{separate.em}}{Logical. Whether to separate assessment models (default = FALSE).}
#'     \item{\code{separate.em.type}}{Integer. Used only if \code{separate.em = TRUE}:}
#'       \itemize{
#'         \item 1 = Panmictic (spatially-aggregated model).
#'         \item 2 = Spatially implicit assessment model (fleets-as-areas).
#'         \item 3 = Separate region-specific assessment models (n = number of regions).
#'       }
#'     \item{\code{do.move}}{Logical. Whether to include movement in the assessment model (default = FALSE).}
#'     \item{\code{est.move}}{Logical. Whether to estimate movement rates in the assessment model (default = FALSE).}
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
#' @param reduce_region_info List (optional). Remove specific regions from the assessment model. If `NULL`, no modifications are applied.
#'   The expected components include:
#'   \itemize{
#'     \item `$remove_regions` Integer vector (`n_regions`). 0/1 values flagging whether to include regions.
#'     \item `$reassign` Numeric. Specify reassignment for surveys from the "removed" regions to "remaining" regions.
#'     \item `$NAA_where`Integer array (`n_stocks × n_regions × n_ages`) 0/1 values flagging whether the stock from a certain age can be present in a specific region.
#'     \item `$sel_em` Selectivity configuration for the "reduced" assessment model.
#'     \item `$M_em` Natural mortality configuration for the "reduced" assessment model.
#'     \item `$NAA_re_em` Numbers-at-age configuration for the "reduced" assessment model.
#'     \item `$move_em` Movement configuration for the "reduced" assessment model.
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
#' @param filter_indices Integer (0/1) vector (optional). User-specified which indices are excluded from the assessment model. For example, c(1,0,1,1) indicates Index 1 (include) and 2 (exclude) in region 1, Index 3 (include) and 4 (include) in region 2
#' @param agg_index_sigma Matrix. Either full (n_years x n_indices) or subset (length(ind_em) x n_indices).
#' @param index_Neff Matrix. Same dimension as `agg_index_sigma`.
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
#' @param assess_years Vector of years when assessments are conducted.
#' @param assess_interval Integer. The interval between assessments in the MSE feedback loop.
#' @param base_years Vector of years used in the burn-in period.
#' @param year.use Integer. Number of years included in the assessment model (default = 20).
#' @param add.years Logical. Whether to use entire time series of data in the assessment model (default = FALSE).
#' @param by_fleet Logical. Whether to calculate fleet-specific fishing mortality (default = TRUE).
#' @param FXSPR_init Numeric. Change initial fishing mortality for estimating reference point in the assessment model (default = NULL).
#' @param hcr A list specifying the harvest control rule (HCR) type and options. It should contain:
#'   \itemize{
#'     \item{\code{hcr.type}}{Integer. The type of harvest control rule to use. Options include:}
#'     \itemize{
#'       \item \code{1}: Annual projected catch based on 75\% of F40\% (default).
#'       \item \code{2}: Constant catch based on 75\% of F40\%.
#'       \item \code{3}: "Hockey stick" approach based on stock status.
#'     }
#'     \item{\code{hcr.opts}}{A list of additional options for the harvest control rule, including:}
#'     \itemize{
#'       \item{\code{use.FXSPR}}{Logical. Calculate and use F at X\% SPR for projections (default = TRUE).}
#'       \item{\code{percentSPR}}{Numeric. The X\% SPR used to calculate reference points (default = 40).}
#'       \item{\code{percentFXSPR}}{Numeric. The percentage of F\_XSPR to use for calculating catch in projections (default = 75).}
#'       \item{\code{use.FMSY}}{Logical. Calculate and use FMSY for projections (default = FALSE).}
#'       \item{\code{percentFMSY}}{Numeric. Percentage of F_MSY to use for catch projections. Used only if use.FMSY = TRUE (default = 100).}
#'       \item{\code{avg.yrs}}{Integer. The number of years over which mean values of MAA, WAA, M, and F will be averaged in projections (default = 5).}
#'       \item{\code{max_percent}}{Numeric. Maximum percentage of F\_XSPR to use in projections (default = 75).}
#'       \item{\code{min_percent}}{Numeric. Minimum percentage of F\_XSPR to use in projections (default = 0.01).}
#'       \item{\code{BThresh_up}}{Numeric. Upper bound of overfished biomass threshold (default = 0.5).}
#'       \item{\code{BThresh_low}}{Numeric. Lower bound of overfished biomass threshold (default = 0.1).}
#'       \item{\code{cont.M.re}}{Logical. Whether to continue natural mortality random effects (e.g., AR1_y or 2D AR1) during projections (default = FALSE).}
#'       \item{\code{cont.move.re}}{Logical. Whether to continue movement random effects during projections (default = FALSE).}
#'     }
#'   }
#'  
#' @param catch_alloc List that contains specifications for catch allocation:
#'   \itemize{
#'     \item `$weight_type` - Integer (1–4) indicating the type of weighting used.
#'       \itemize{
#'         \item `1` - Equal weighting.
#'         \item `2` - Weighting based on past catch data.
#'         \item `3` - Weighting based on past survey index data.
#'         \item `4` - User-specified weights.
#'       }
#'     \item `$method` - String. Specify the catch allocation method:
#'       \itemize{
#'         \item `"equal"` - Use equal weights, split total catch into fleet-specific catches (use when weight_type = 1).
#'         \item `"fleet_equal"` - Use equal weights, split total fleet catch into fleet-specific catches (use when weight_type = 1).
#'         \item `"fleet_region"` - Use regional catch totals to compute weights (use when weight_type = 2).
#'         \item `"fleet_gear"` - Use gear-type catch totals to compute weights (use when weight_type = 2).
#'         \item `"fleet_combined"` - First allocates to gear-specific total catch, then splits into regions (use when weight_type = 2).
#'         \item `"fleet_catch"` - Use fleet-specific catch to compute weights (use when weight_type = 2).
#'         \item `"index_equal"` - Use survey-based regional weighting, then assigns equally among fleets in the same region (use when weight_type = 3).
#'         \item `"index_gear"` - Use survey-based regional weighting, then assigns based on gear-specific weights (use when weight_type = 3).
#'         \item `"index_multiple"` - Use multiple survey-based regional weighting, then assigns equally among fleets in the same region (use when weight_type = 3).
#'         \item `"user_defined_fleets"` - Catch allocation is based on user-defined weights for each fleet.
#'         \item `"user_defined_regions"` - Catch is allocated based on user-defined weights for regions, then equally distributed among fleets within the region.
##'       }
#'     \item `$user_weights` - Numeric vector (`n_regions` or `n_fleets`). Optional. User-defined weights summing to 1 (use when weight_type = 4).
#'     \item `$weight_years` - Integer. Number of years to average for calculating historical catch weights.
#'     \item `$survey_pointer` - Integer vector. Specify which survey index type to use for weighting (use when weight_type = 3).
#'   }
#' @param do.retro Logical. Whether to perform retrospective analysis for the assessment model (default = FALSE).
#' @param do.osa Logical. Whether to calculate one-step-ahead (OSA) residuals for the assessment model (default = FALSE).
#' @param do.brps Logical. Whether to calculate reference points in the operating model (default = FALSE).
#' @param seed Integer. The random seed used for stochastic processes.
#' @param save.sdrep Logical. Whether to save the results of every assessment model (memory-intensive).
#' @param save.last.em Logical. Whether to save the last estimation model (default = FALSE).
#' 
#' @return A list containing:
#'   \describe{
#'     \item{\code{om}}{Updated operating model.}
#'     \item{\code{em_list}}{List of estimation model results across assessment years.}
#'     \item{\code{par.est}}{Estimated parameter values.}
#'     \item{\code{par.se}}{Standard errors of parameter estimates.}
#'     \item{\code{adrep.est}}{Estimated values from ADMB report.}
#'     \item{\code{adrep.se}}{Standard errors from ADMB report.}
#'     \item{\code{opt_list}}{Optimization results from WHAM.}
#'     \item{\code{converge_list}}{Convergence diagnostics.}
#'     \item{\code{catch_advice}}{Projected catch advice over assessment years.}
#'     \item{\code{em_full}}{List of full estimation model outputs.}
#'     \item{\code{runtime}}{Elapsed time for function execution.}
#'     \item{\code{seed.save}}{Random seed.}
#'     \item{\code{em_input}}{List of input information for the last assessment.}
#'   }
#'   
#' @export
#' 
#' @seealso \code{\link{make_em_input}}, \code{\link{update_om_fn}}, \code{\link{advice_fn}}


loop_through_fn <- function(om, 
                            em_info = NULL, 
                            random = NULL, 
                            M_em = NULL, 
                            sel_em = NULL, 
                            NAA_re_em = NULL, 
                            move_em = NULL, 
                            age_comp_em = "multinomial", 
                            em.opt = list(separate.em = TRUE, separate.em.type = 1,
                                          do.move = FALSE, est.move = FALSE), 
                            aggregate_catch_info = NULL,
                            aggregate_index_info = NULL,
                            aggregate_weights_info = NULL,
                            reduce_region_info = NULL,
                            filter_indices = NULL,
                            update_catch_info = NULL,
                            update_index_info = NULL,
                            assess_years = NULL, 
                            assess_interval = NULL, 
                            base_years = NULL, 
                            year.use = 20, 
                            add.years = FALSE,
                            by_fleet = TRUE,
                            FXSPR_init = NULL,
                            hcr = list(hcr.type = 1, hcr.opts = NULL),
                            catch_alloc = list(weight_type = 1, method = "equal", user_weights = NULL, weight_years = 1),
                            do.retro = FALSE, 
                            do.osa = FALSE, 
                            do.brps = FALSE,
                            seed = 123, 
                            save.sdrep = FALSE, 
                            save.last.em = FALSE
                            ) {
  
  start.time <- Sys.time()
  
  # Helper function to check convergence
  check_conv <- function(em) {
    conv <- as.logical(1 - em$opt$convergence)
    pdHess <- as.logical(if (em$na_sdrep == FALSE & !is.na(em$na_sdrep)) 1 else 0)
    if (!conv | !pdHess) warnings("Assessment model is not converged!")
    list(conv = conv, pdHess = pdHess)
  }
  
  if (is.null(em.opt)) stop("em.opt has to be specified!")
  if (!is.null(move_em) & em.opt$separate.em) stop("move_em must be NULL if em.opt$separate.em = TRUE!")
  if (em.opt$separate.em) move_em <- NULL
  
  em_list <- list()
  par.est <- list()
  par.se <- list()
  adrep.est <- list()
  adrep.se <- list()
  opt_list <- list()
  converge_list <- list()
  catch_advice <- list()
  em_full <- list()
  
  if(is.null(age_comp_em)) age_comp_em = "multinomial"
  if(is.null(em_info)) stop("em_info must be specified!")
  
  if (em.opt$separate.em) {
    
    for (y in assess_years) {
      
      cat(paste0("\nNow conducting stock assessment for year ", y, "\n"))
      i <- which(assess_years == y)
      em.years <- base_years[1]:y
      
      if (add.years && i != 1) year.use = year.use + assess_interval
        
      em_input <- make_em_input(om = om, em_info = em_info, M_em = M_em, sel_em = sel_em,
                                NAA_re_em = NAA_re_em, move_em = move_em, em.opt = em.opt,
                                em_years = em.years, year.use = year.use, age_comp_em = age_comp_em,
                                aggregate_catch_info = aggregate_catch_info,
                                aggregate_index_info = aggregate_index_info,
                                aggregate_weights_info = aggregate_weights_info,
                                filter_indices = filter_indices)
      
      #cat("\nNow agg Catch is..", em_input$data$agg_catch[1,])
      
      if(!is.null(FXSPR_init)) em_input$data$FXSPR_init[] = FXSPR_init
      
      n_stocks <- om$input$data$n_stocks
      
      if (em.opt$separate.em.type == 1) {
        
        cat("\nNow fitting assessment model...\n")
        em <- fit_wham(em_input, do.retro = do.retro, do.osa = do.osa, do.brps = TRUE, MakeADFun.silent = TRUE)
     
        cat("\nNow checking convergence of assessment model...\n")
        conv <- check_conv(em)$conv
        pdHess <- check_conv(em)$pdHess
        if (conv & pdHess) cat("\nAssessment model is converged.\n") else warnings("\nAssessment model is not converged!\n")
        
        cat("\nNow using the EM to project catch...\n")
        
        em.advice <- advice_fn(em, pro.yr = assess_interval, hcr)
        
        if(is.vector(em.advice)) em.advice = matrix(em.advice, byrow = TRUE)
        
        cat("\nProject catch from assessment model is ", em.advice, "\n")
        
        cat("\nNow allocating catch...\n")
        
        advice <- calculate_catch_advice(om, em.advice, 
                                         aggregate_catch_info, 
                                         aggregate_index_info, 
                                         final_year = y,
                                         catch_alloc)
        
        colnames(advice) <- paste0("Fleet_", 1:om$input$data$n_fleets)
        rownames(advice) <- paste0("Year_", y + 1:assess_interval)
        
        cat("\nNow generating catch advice...\n")
        print(advice)
        
        interval.info <- list(catch = advice, years = y + 1:assess_interval)
        
        cat("\nNow calculating F at age in the OM given the catch advice...\n")
        
        om <- update_om_fn(om, interval.info, seed = seed, random = random, method = "nlminb", by_fleet = by_fleet, do.brps = do.brps)
        
        em_list[[i]] <- em$rep
        par.est[[i]] <- as.list(em$sdrep, "Estimate")
        par.se[[i]] <- as.list(em$sdrep, "Std. Error")
        adrep.est[[i]] <- as.list(em$sdrep, "Estimate", report = TRUE)
        adrep.se[[i]] <- as.list(em$sdrep, "Std. Error", report = TRUE)
        opt_list[[i]] <- em$opt
        converge_list[[i]] <- conv + pdHess
        catch_advice[[i]] <- advice
        
        if (save.sdrep) {
          em_full[[i]] <- em
        } else {
          if (y == assess_years[length(assess_years)]) {
            em_full[[1]] <- em
          }
          if (!save.last.em) em_full[[1]] <- list()
        }
      } else if (em.opt$separate.em.type == 2) {
        
        cat("\nNow fitting assessment model...\n")
        em <- fit_wham(em_input, do.retro = do.retro, do.osa = do.osa, do.brps = TRUE, MakeADFun.silent = TRUE)

        cat("\nNow checking convergence of assessment model...\n")
        conv <- check_conv(em)$conv
        pdHess <- check_conv(em)$pdHess
        if (conv & pdHess) cat("\nAssessment model is converged.\n") else warnings("\nAssessment model is not converged!\n")
        
        cat("\nNow using the EM to project catch...\n")
        # advice <- advice_fn(em, pro.yr = assess_interval, hcr.type = hcr.type, hcr.opts = hcr.opts)
        advice <- advice_fn(em, pro.yr = assess_interval, hcr)
        
        if(is.vector(advice)) advice <- as.matrix(t(advice))
        colnames(advice) <- paste0("Fleet_", 1:om$input$data$n_fleets)
        rownames(advice) <- paste0("Year_", y + 1:assess_interval)
        
        cat("\nNow generating catch advice...\n")
        print(advice)
        
        interval.info <- list(catch = advice, years = y + 1:assess_interval)
        
        cat("\nNow calculating F at age in the OM given the catch advice...\n")
        
        om <- update_om_fn(om, interval.info, seed = seed, random = random, method = "nlminb", by_fleet = by_fleet, do.brps = do.brps)
        
        em_list[[i]] <- em$rep
        par.est[[i]] <- as.list(em$sdrep, "Estimate")
        par.se[[i]] <- as.list(em$sdrep, "Std. Error")
        adrep.est[[i]] <- as.list(em$sdrep, "Estimate", report = TRUE)
        adrep.se[[i]] <- as.list(em$sdrep, "Std. Error", report = TRUE)
        opt_list[[i]] <- em$opt
        converge_list[[i]] <- conv + pdHess
        catch_advice[[i]] <- advice
        
        if (save.sdrep) {
          em_full[[i]] <- em
        } else {
          if (y == assess_years[length(assess_years)]) {
            em_full[[1]] <- em
          }
          if (!save.last.em) em_full[[1]] <- list()
        }
      } else if (em.opt$separate.em.type == 3) {
        
        em_list[[i]] <- list()
        par.est[[i]] <- list()
        par.se[[i]] <- list()
        adrep.est[[i]] <- list()
        adrep.se[[i]] <- list()
        opt_list[[i]] <- list()
        converge_list[[i]] <- list()
        em_full[[i]] <- list()
        
        advice <- NULL
        em <- list()
        conv <- rep(0, n_stocks)
        pdHess <- rep(0, n_stocks)
        
        cat("\nNow generating catch advice...\n")
        for (s in 1:n_stocks) {
          
          cat("\nNow fitting assessment model...\n")
          em[[s]] <- fit_wham(em_input[[s]], do.retro = do.retro, do.osa = do.osa, do.brps = TRUE, MakeADFun.silent = TRUE)
          
          cat("\nNow checking convergence of assessment model...\n")
          conv <- check_conv(em[[s]])$conv
          pdHess <- check_conv(em[[s]])$pdHess
          if (conv & pdHess) cat("\nAssessment model is converged.\n") else warnings("\nAssessment model is not converged!\n")
          
          tmp <- advice_fn(em[[s]], pro.yr = assess_interval, hcr)
          advice <- cbind(advice, tmp)
        }
        
        if(is.vector(advice)) advice <- as.matrix(t(advice))
        colnames(advice) <- paste0("Fleet_", 1:om$input$data$n_fleets)
        rownames(advice) <- paste0("Year_", assess_years[i] + 1:assess_interval)
        
        print(advice)
        
        # set the catch for the next assess_interval years
        interval.info <- list(catch = advice, years = assess_years[i] + 1:assess_interval)
        
        cat("\nNow calculating F at age in the OM given the catch advice...\n")
        om <- update_om_fn(om, interval.info, seed = seed, random = random, method = "nlminb", by_fleet = by_fleet, do.brps = do.brps)
        
        for (s in 1:n_stocks) {
          em_list[[i]][[s]] <- em[[s]]$rep
          par.est[[i]][[s]] <- as.list(em[[s]]$sdrep, "Estimate")
          par.se[[i]][[s]] <- as.list(em[[s]]$sdrep, "Std. Error")
          adrep.est[[i]][[s]] <- as.list(em[[s]]$sdrep, "Estimate", report = TRUE)
          adrep.se[[i]][[s]] <- as.list(em[[s]]$sdrep, "Std. Error", report = TRUE)
          opt_list[[i]][[s]] <- em[[s]]$opt
          converge_list[[i]] <- sum(conv, pdHess)
          catch_advice[[i]] <- advice
          
          if (save.sdrep) {
            em_full[[i]][[s]] <- em[[s]]
          } else {
            if (y == assess_years[length(assess_years)]) {
              em_full[[1]][[s]] <- em[[s]]
            }
            if (!save.last.em) em_full[[1]][[s]] <- list()
          }
        }
      }
    }
  } else {
    for (y in assess_years) {
     
      cat(paste0("\nNow conducting stock assessment for year ", y, "\n"))
      
      i <- which(assess_years == y)
      em.years <- base_years[1]:y
      
      if (add.years && i != 1) year.use = year.use + assess_interval
      
      em_input <- make_em_input(om = om, em_info = em_info, M_em = M_em, sel_em = sel_em,
                                NAA_re_em = NAA_re_em, move_em = move_em, em.opt = em.opt,
                                em_years = em.years, year.use = year.use, age_comp_em = age_comp_em,
                                aggregate_catch_info = aggregate_catch_info,
                                aggregate_index_info = aggregate_index_info,
                                filter_indices = filter_indices,
                                reduce_region_info = reduce_region_info,
                                update_catch_info = update_catch_info,
                                update_index_info = update_index_info) 
      
      if(!is.null(FXSPR_init)) em_input$data$FXSPR_init[] = FXSPR_init
      
      cat("\nNow fitting assessment model...\n")
      
      if (em.opt$do.move) {
        if (em.opt$est.move) {
          em <- fit_wham(em_input, do.retro = do.retro, do.osa = do.osa, do.brps = TRUE, MakeADFun.silent = TRUE)
        } else {
          em_input <- fix_move(em_input)
          em <- fit_wham(em_input, do.retro = do.retro, do.osa = do.osa, do.brps = TRUE, MakeADFun.silent = TRUE)
        }
      } else {
        em <- fit_wham(em_input, do.retro = do.retro, do.osa = do.osa, do.brps = TRUE, MakeADFun.silent = TRUE)
      }
      
      cat("\nNow checking convergence of assessment model...\n")
      conv <- check_conv(em)$conv
      pdHess <- check_conv(em)$pdHess
      if (conv & pdHess) cat("\nAssessment model is converged.\n") else warnings("\nAssessment model is not converged!\n")
      
      cat("\nNow generating catch advice...\n")
      advice <- advice_fn(em, pro.yr = assess_interval, hcr)
      if(!is.null(reduce_region_info$remove_regions)) {
        remove_regions = reduce_region_info$remove_regions
        fleets_to_remove <- which(om$input$data$fleet_regions %in% which(remove_regions == 0))  # Get fleet indices
        fleets_to_keep <- which(!om$input$data$fleet_regions %in% which(remove_regions == 0))  # Get fleet indices
        advice.tmp <- matrix(0, nrow = assess_interval, ncol = length(om$input$data$fleet_regions))
        advice.tmp[,fleets_to_keep] = advice
        if (!is.null(reduce_region_info$fleet_catch)){
          advice.tmp[,fleets_to_remove] = reduce_region_info$fleet_catch # this can be a vector then the value will be filled by vertical and then by horizontal, OR can be a matrix (nrow = assess_interval x ncol = n_fleets_to_keep).
        }
        advice <- advice.tmp
      }
      
      if(is.vector(advice)) {
        if(assess_interval == 1) {
          advice <- as.matrix(t(advice))
        } else {
          advice <- matrix(advice, byrow = T)
        }
      }
        
      colnames(advice) <- paste0("Fleet_", 1:om$input$data$n_fleets)
      rownames(advice) <- paste0("Year_", y + 1:assess_interval)
      
      cat("Catch Advice \n",advice,"\n")
      
      interval.info <- list(catch = advice, years = y + 1:assess_interval)
      
      cat("\nNow calculating F at age in the OM given the catch advice...\n")
      om <- update_om_fn(om, interval.info, seed = seed, random = random, method = "nlminb", by_fleet = by_fleet, do.brps = do.brps)
      
      em_list[[i]] <- em$rep
      par.est[[i]] <- as.list(em$sdrep, "Estimate")
      par.se[[i]] <- as.list(em$sdrep, "Std. Error")
      adrep.est[[i]] <- as.list(em$sdrep, "Estimate", report = TRUE)
      adrep.se[[i]] <- as.list(em$sdrep, "Std. Error", report = TRUE)
      opt_list[[i]] <- em$opt
      converge_list[[i]] <- conv + pdHess
      catch_advice[[i]] <- advice
      
      if (save.sdrep) {
        em_full[[i]] <- em
      } else {
        if (y == assess_years[length(assess_years)]) {
          em_full[[1]] <- em
        }
        if (!save.last.em) em_full[[1]] <- list()
      }
    }
  }
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  cat("Please ignore Warning in check_projF(proj_mod).")
  cat("\nTotal Runtime = ", time.taken,"\n")
  
  return(list(om = om, em_list = em_list, par.est = par.est, par.se = par.se, 
              adrep.est = adrep.est, adrep.se = adrep.se, opt_list = opt_list, 
              converge_list = converge_list, catch_advice = catch_advice, em_full = em_full,
              runtime = time.taken, seed.save = seed, em_input = em_input))
}
