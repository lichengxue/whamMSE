#' Perform Management Strategy Evaluation (MSE)
#'
#' A wrapper function to iterate through the feedback period in a management strategy evaluation (MSE),
#' updating the operating model (OM), refitting the assessment model (EM), and generating catch advice.
#'
#' @param om List. The operating model containing pseudo data and simulation configurations.
#' @param em_info List. Information used to generate the assessment model (default = NULL).
#' @param random Character vector. Processes treated as random effects in the operating model (default = "log_NAA").
#' @param M_em Configuration for natural mortality in the assessment model.
#' @param sel_em Configuration for selectivity in the assessment model.
#' @param NAA_re_em Configuration for NAA in the assessment model.
#' @param move_em Configuration for movement in the assessment model.
#' @param catchability_em Configuration for survey catchability in the assessment model.
#' @param ecov_em Configuration for environmental covariates in the assessment model.
#' @param age_comp_em Character. Likelihood distribution for age composition data in the assessment model.
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
#'
#' @param em.opt List. EM model options.
#'   \describe{
#'     \item{\code{separate.em}}{Logical. Use region-separated EMs (default = FALSE).}
#'     \item{\code{separate.em.type}}{Integer. Type of separation:
#'       \itemize{
#'         \item 1: Panmictic
#'         \item 2: Fleets-as-areas
#'         \item 3: Separate regional EMs
#'       }
#'     }
#'     \item{\code{do.move}, \code{est.move}}{Logical. Include/estimate movement in EM (default = FALSE).}
#'   }
#'
#' @param aggregate_catch_info List (optional). Aggregation settings for catch data.
#'   \describe{
#'     \item{\code{n_fleets}}{Integer. Number of fleets.}
#'     \item{\code{catch_cv}}{Numeric vector. CV for each fleet.}
#'     \item{\code{catch_Neff}}{Numeric vector. Effective sample sizes.}
#'     \item{\code{use_agg_catch}}{Integer vector. 0/1 for using aggregated catch.}
#'     \item{\code{use_catch_paa}}{Integer vector. 0/1 for including catch PAA.}
#'     \item{\code{fleet_pointer}}{Integer vector. Grouping fleets for aggregation (0 = exclude).}
#'     \item{\code{use_catch_weighted_waa}}{Logical. Use fleet-catch-weighted WAA.}
#'   }
#'
#' @param aggregate_index_info List (optional). Aggregation settings for index data.
#'   \describe{
#'     \item{\code{n_indices}}{Integer. Number of indices.}
#'     \item{\code{index_cv}}{Numeric vector. CV for each index.}
#'     \item{\code{index_Neff}}{Numeric vector. Effective sample sizes.}
#'     \item{\code{fracyr_indices}}{Numeric vector. Fraction of the year for each index.}
#'     \item{\code{q}}{Numeric vector. Initial catchability values.}
#'     \item{\code{use_indices}}{Integer vector. 0/1 for using the index.}
#'     \item{\code{use_index_paa}}{Integer vector. 0/1 for including index PAA.}
#'     \item{\code{units_indices}}{Integer vector. 1 = biomass, 2 = numbers.}
#'     \item{\code{units_index_paa}}{Integer vector. 1 = biomass, 2 = numbers.}
#'     \item{\code{index_pointer}}{Integer vector. Grouping indices (0 = exclude).}
#'     \item{\code{use_index_weighted_waa}}{Logical. Use survey-index-weighted WAA.}
#'   }
#'
#' @param aggregate_weights_info List (optional). Used to compute weighted average weight-at-age and maturity-at-age.
#'   \describe{
#'     \item{\code{ssb_waa_weights}}{List.
#'       \describe{
#'         \item{\code{fleet}}{Logical. Use fleet-specific weights.}
#'         \item{\code{index}}{Logical. Use index-specific weights.}
#'         \item{\code{pointer}}{Integer. Points to fleet/index group to use.}
#'       }
#'     }
#'     \item{\code{maturity_weights}}{List.
#'       \describe{
#'         \item{\code{fleet}}{Logical.}
#'         \item{\code{index}}{Logical.}
#'         \item{\code{pointer}}{Integer.}
#'       }
#'     }
#'   }
#'
#' @param reduce_region_info List (optional). Reduce regions and reassign data.
#'   \describe{
#'     \item{\code{remove_regions}}{Integer vector. 0/1 flag to remove regions.}
#'     \item{\code{reassign}}{Numeric. Region reassignment code.}
#'     \item{\code{NAA_where}}{3D array (stock × region × age) indicating stock presence.}
#'     \item{\code{sel_em}, \code{M_em}, \code{NAA_re_em}, \code{move_em}}{EM configs for reduced regions.}
#'     \item{\code{onto_move_list}}{List.
#'       \describe{
#'         \item{\code{onto_move}}{3D array (stock × region × region).}
#'         \item{\code{onto_move_pars}}{4D array (stock × region × region × 4).}
#'         \item{\code{age_mu_devs}}{4D array (stock × region × region × age).}
#'       }
#'     }
#'   }
#'
#' @param catch_alloc List. Catch allocation specifications.
#'   \describe{
#'     \item{\code{weight_type}}{Integer (1–4). Weighting method.}
#'     \item{\code{method}}{Character. Allocation method:
#'       \itemize{
#'         \item \code{"equal"}
#'         \item \code{"fleet_equal"}
#'         \item \code{"fleet_region"}
#'         \item \code{"fleet_gear"}
#'         \item \code{"fleet_combined"}
#'         \item \code{"fleet_catch"}
#'         \item \code{"index_equal"}
#'         \item \code{"index_gear"}
#'         \item \code{"multiple_index_equal"}
#'         \item \code{"multiple_index_gear"}
#'         \item \code{"user_defined_fleets"}
#'         \item \code{"user_defined_regions"}
#'       }
#'     }
#'     \item{\code{user_weights}}{Numeric vector. User-defined weights (must sum to 1).}
#'     \item{\code{weight_years}}{Integer. Number of past years used to average weights.}
#'     \item{\code{survey_pointer}}{Integer vector. Indices used for weighting (when \code{weight_type = 3}).}
#'   }
#'
#' @param add_implementation_error List (optional). Adds variability to catch realization.
#'   \describe{
#'     \item{\code{method}}{Character. Distribution type: \code{"lognormal"}, \code{"normal"}, \code{"uniform"}, or \code{"constant"}.}
#'     \item{\code{mean}}{Numeric. Mean of error distribution (log scale if lognormal).}
#'     \item{\code{cv}}{Numeric. Coefficient of variation (required for lognormal).}
#'     \item{\code{sd}}{Numeric. Standard deviation (required for normal).}
#'     \item{\code{min}}{Numeric. Lower bound for uniform distribution.}
#'     \item{\code{max}}{Numeric. Upper bound for uniform distribution.}
#'     \item{\code{constant_value}}{Numeric. Fixed multiplier (used if \code{method = "constant"}).}
#'   }
#'
#' @param filter_indices Integer vector (optional). 0/1 vector to exclude survey indices by region.
#'
#' @param update_catch_info List (optional). Update catch CV and Neff values in the assessment model.
#'   \describe{
#'     \item{\code{agg_catch_sigma}}{Matrix. Catch CVs (standard deviation of log catch).}
#'     \item{\code{catch_Neff}}{Matrix. Effective sample sizes for catch. Must match dimensions of \code{agg_catch_sigma}.}
#'   }
#'
#' @param update_index_info List (optional). Update index CV and Neff values in the assessment model.
#'   \describe{
#'     \item{\code{agg_index_sigma}}{Matrix. Aggregated survey index CVs.}
#'     \item{\code{index_Neff}}{Matrix. Effective sample sizes for index data.}
#'     \item{\code{remove_agg}}{Logical. Whether to remove aggregated index observations for specific pointers/years.}
#'     \item{\code{remove_agg_pointer}}{Integer vector. Index pointers to remove.}
#'     \item{\code{remove_agg_years}}{Integer vector or matrix. Years to remove.}
#'     \item{\code{remove_paa}}{Logical. Whether to remove index PAA observations.}
#'     \item{\code{remove_paa_pointer}}{Integer vector. Index pointers to remove PAA.}
#'     \item{\code{remove_paa_years}}{Integer vector or matrix. Years to remove index PAA data.}
#'   }
#'
#' @param hcr List. Harvest control rule.
#'   \describe{
#'     \item{\code{hcr.type}}{Integer (1: F\_XSPR, 2: constant catch, 3: hockey stick).}
#'     \item{\code{hcr.opts}}{List of options:
#'       \itemize{
#'         \item \code{use_FXSPR}
#'         \item \code{percentSPR}
#'         \item \code{percentFXSPR}
#'         \item \code{use_FMSY}
#'         \item \code{percentFMSY}
#'         \item \code{avg_yrs}
#'         \item \code{max_percent}
#'         \item \code{min_percent}
#'         \item \code{BThresh_up}
#'         \item \code{BThresh_low}
#'         \item \code{cont.M.re}
#'         \item \code{cont.move.re}
#'       }
#'     }
#'   }
#'
#' @param user_SPR_weights_info List. Update SPR weights for calculating biological reference points.
#'   \describe{
#'     \item{\code{method}}{Character. Specifies how weights are assigned to regions or stocks.}
#'     \item{\code{weight_years}}{Integer. Number of years to average catch or index values.}
#'     \item{\code{index_pointer}}{Integer. Index number used when \code{method = "index_region"}.}
#'   }
#'
#' @param assess_years Integer vector. Assessment years.
#' @param assess_interval Integer. Interval (years) between assessments.
#' @param base_years Integer vector. Burn-in period years.
#' @param year.use Integer. Years used in each EM (default = 20).
#' @param add.years Logical. Use entire time series in EM (default = FALSE).
#' @param by_fleet Logical. Calculate fleet-specific F (default = TRUE).
#' @param FXSPR_init Numeric. Initial F for reference point calculation (optional).
#' @param do.retro Logical. Run retrospective analysis (default = FALSE).
#' @param do.osa Logical. Calculate OSA residuals (default = FALSE).
#' @param do.brps Logical. Calculate reference points in OM (default = FALSE).
#' @param seed Integer. Random seed.
#' @param save.sdrep Logical. Save sdrep objects (default = FALSE).
#' @param save.last.em Logical. Save final EM (default = FALSE).
#'
#' @return List with elements.
#'   \describe{
#'     \item{\code{om}}{Updated OM.}
#'     \item{\code{em_list}}{List of EM results.}
#'     \item{\code{par.est}, \code{par.se}}{Parameter estimates and SEs.}
#'     \item{\code{adrep.est}, \code{adrep.se}}{ADMB report values and SEs.}
#'     \item{\code{opt_list}}{Optimization results.}
#'     \item{\code{converge_list}}{Convergence checks.}
#'     \item{\code{catch_advice}, \code{catch_realized}}{Advised and realized catch.}
#'     \item{\code{em_full}}{Full EM outputs.}
#'     \item{\code{em_input}}{Final EM input.}
#'     \item{\code{runtime}}{Elapsed time.}
#'     \item{\code{seed.save}}{Final random seed.}
#'   }
#'
#' @seealso \code{\link{make_em_input}}, \code{\link{update_om_fn}}, \code{\link{advice_fn}}
#' @export


loop_through_fn <- function(om, 
                            em_info = NULL, 
                            random = NULL, 
                            M_em = NULL, 
                            sel_em = NULL, 
                            NAA_re_em = NULL, 
                            move_em = NULL, 
                            catchability_em = NULL,
                            ecov_em = NULL,
                            ecov_em_opts = NULL,
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
                            user_SPR_weights_info = NULL,
                            assess_years = NULL, 
                            assess_interval = NULL, 
                            base_years = NULL, 
                            year.use = 20, 
                            add.years = FALSE,
                            by_fleet = TRUE,
                            FXSPR_init = NULL,
                            hcr = list(hcr.type = 1, hcr.opts = NULL),
                            proj.opts = list(),   # <-- NEW
                            catch_alloc = list(weight_type = 1, method = "equal", user_weights = NULL, weight_years = 1),
                            implementation_error = NULL,
                            do.retro = FALSE, 
                            do.osa = FALSE, 
                            do.brps = FALSE,
                            seed = 123, 
                            save.sdrep = FALSE, 
                            save.last.em = FALSE) {
  
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
  catch_realized <- list()
  em_full <- list()
  em_input_list <- list()
  
  if(is.null(age_comp_em)) age_comp_em = "multinomial"
  if(is.null(em_info)) stop("em_info must be specified!")
  
  if (em.opt$separate.em) {
    
    for (y in assess_years) {
      
      cat(paste0("\nNow conducting stock assessment for year ", y, "\n"))
      i <- which(assess_years == y)
      em.years <- base_years[1]:y
      
      if (add.years && i != 1) year.use = year.use + assess_interval
        
      em_input <- make_em_input(om = om, em_info = em_info, 
                                M_em = M_em, sel_em = sel_em,
                                NAA_re_em = NAA_re_em, move_em = move_em, 
                                catchability_em = catchability_em, ecov_em = ecov_em,
                                em.opt = em.opt, em_years = em.years, year.use = year.use, 
                                age_comp_em = age_comp_em,
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
        
        advice <- advice_fn(em = em,
                            pro.yr = assess_interval,
                            hcr = hcr,
                            proj.opts = proj.opts)
        
        if(is.vector(em.advice)) em.advice = matrix(em.advice, byrow = TRUE)
        
        cat("\nProject catch from assessment model is\n")
        print(em.advice)
        
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
        
        if(!is.null(implementation_error)) {
          cat("\nNow generating implementation error on catch advice...\n")
          method = implementation_error$method
          catch_mean = implementation_error$mean
          catch_cv = implementation_error$cv
          catch_sd = implementation_error$sd
          catch_min = implementation_error$min
          catch_max = implementation_error$max
          constant_value = implementation_error$constant_value
          real_catch = add_implementation_error(catch_advice = advice,
                                                method = method,
                                                mean = catch_mean,
                                                cv = catch_cv,
                                                sd = catch_sd,
                                                min = catch_min,
                                                max = catch_max,
                                                constant_value = constant_value,
                                                seed = seed)
          cat("\nRealized catch is...\n")
          print(real_catch)
          interval.info <- list(catch = real_catch, years = y + 1:assess_interval)
        } else {
          interval.info <- list(catch = advice, years = y + 1:assess_interval)
          real_catch = advice
        }

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
        catch_realized[[i]] <- real_catch
        
        if (save.sdrep) {
          em_full[[i]] <- em
        } else {
          if (y == assess_years[length(assess_years)]) {
            em_full[[1]] <- em
          }
          if (!save.last.em) em_full[[1]] <- list()
        }
        
        em_input_list[[i]] <- em_input
        
      } else if (em.opt$separate.em.type == 2) {
        
        cat("\nNow fitting assessment model...\n")
        em <- fit_wham(em_input, do.retro = do.retro, do.osa = do.osa, do.brps = TRUE, MakeADFun.silent = TRUE)

        cat("\nNow checking convergence of assessment model...\n")
        conv <- check_conv(em)$conv
        pdHess <- check_conv(em)$pdHess
        if (conv & pdHess) cat("\nAssessment model is converged.\n") else warnings("\nAssessment model is not converged!\n")
        
        cat("\nNow using the EM to project catch...\n")
        
        advice <- advice_fn(em = em,
                            pro.yr = assess_interval,
                            hcr = hcr,
                            proj.opts = proj.opts)
        
        
        if(is.vector(advice)) advice <- as.matrix(t(advice))
        colnames(advice) <- paste0("Fleet_", 1:om$input$data$n_fleets)
        rownames(advice) <- paste0("Year_", y + 1:assess_interval)
        
        cat("\nNow generating catch advice...\n")
        print(advice)
        
        if(!is.null(implementation_error)) {
          cat("\nNow generating implementation error on catch advice...\n")
          method = implementation_error$method
          catch_mean = implementation_error$mean
          catch_cv = implementation_error$cv
          catch_sd = implementation_error$sd
          catch_min = implementation_error$min
          catch_max = implementation_error$max
          constant_value = implementation_error$constant_value
          real_catch = add_implementation_error(catch_advice = advice,
                                                method = method,
                                                mean = catch_mean,
                                                cv = catch_cv,
                                                sd = catch_sd,
                                                min = catch_min,
                                                max = catch_max,
                                                constant_value = constant_value,
                                                seed = seed)
          cat("\nRealized catch is...\n")
          print(real_catch)
          interval.info <- list(catch = real_catch, years = y + 1:assess_interval)
        } else {
          interval.info <- list(catch = advice, years = y + 1:assess_interval)
          real_catch = advice
        }
        
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
        catch_realized[[i]] <- real_catch
        
        if (save.sdrep) {
          em_full[[i]] <- em
        } else {
          if (y == assess_years[length(assess_years)]) {
            em_full[[1]] <- em
          }
          if (!save.last.em) em_full[[1]] <- list()
        }
        
        em_input_list[[i]] <- em_input
        
      } else if (em.opt$separate.em.type == 3) {
        
        em_list[[i]] <- list()
        par.est[[i]] <- list()
        par.se[[i]] <- list()
        adrep.est[[i]] <- list()
        adrep.se[[i]] <- list()
        opt_list[[i]] <- list()
        converge_list[[i]] <- list()
        em_full[[i]] <- list()
        em_input_list[[i]] <- list()
        
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
          
          tmp <- advice_fn(em = em[[s]],
                           pro.yr = assess_interval,
                           hcr = hcr,
                           proj.opts = proj.opts)
          
          advice <- cbind(advice, tmp)
        }
        
        if(is.vector(advice)) advice <- as.matrix(t(advice))
        colnames(advice) <- paste0("Fleet_", 1:om$input$data$n_fleets)
        rownames(advice) <- paste0("Year_", assess_years[i] + 1:assess_interval)
        
        print(advice)
        
        if(!is.null(implementation_error)) {
          cat("\nNow generating implementation error on catch advice...\n")
          method = implementation_error$method
          catch_mean = implementation_error$mean
          catch_cv = implementation_error$cv
          catch_sd = implementation_error$sd
          catch_min = implementation_error$min
          catch_max = implementation_error$max
          constant_value = implementation_error$constant_value
          real_catch = add_implementation_error(catch_advice = advice,
                                                method = method,
                                                mean = catch_mean,
                                                cv = catch_cv,
                                                sd = catch_sd,
                                                min = catch_min,
                                                max = catch_max,
                                                constant_value = constant_value,
                                                seed = seed)
          cat("\nRealized catch is...\n")
          print(real_catch)
          interval.info <- list(catch = real_catch, years = assess_years[i] + 1:assess_interval)
        } else {
          interval.info <- list(catch = advice, years = assess_years[i] + 1:assess_interval)
          real_catch = advice
        }
        
        # set the catch for the next assess_interval years
        # interval.info <- list(catch = advice, years = assess_years[i] + 1:assess_interval)
        
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
          catch_realized[[i]] <- real_catch
          
          if (save.sdrep) {
            em_full[[i]][[s]] <- em[[s]]
          } else {
            if (y == assess_years[length(assess_years)]) {
              em_full[[1]][[s]] <- em[[s]]
            }
            if (!save.last.em) em_full[[1]][[s]] <- list()
          }
          em_input_list[[i]][[s]] <- em_input[[s]]
        }
      }
    }
  } else {
    for (y in assess_years) {
     
      cat(paste0("\nNow conducting stock assessment for year ", y, "\n"))
      
      i <- which(assess_years == y)
      em.years <- base_years[1]:y
      
      if (add.years && i != 1) year.use = year.use + assess_interval
      
      em_input <- make_em_input(om = om, em_info = em_info, 
                                M_em = M_em, sel_em = sel_em,
                                NAA_re_em = NAA_re_em, move_em = move_em, 
                                catchability_em = catchability_em, ecov_em = ecov_em,
                                em.opt = em.opt, em_years = em.years, year.use = year.use, 
                                age_comp_em = age_comp_em,
                                aggregate_catch_info = aggregate_catch_info,
                                aggregate_index_info = aggregate_index_info,
                                filter_indices = filter_indices,
                                reduce_region_info = reduce_region_info,
                                update_catch_info = update_catch_info,
                                update_index_info = update_index_info,
                                ecov_em_opts = ecov_em_opts) 
      
      if(!is.null(user_SPR_weights_info)) {
        if(is.null(user_SPR_weights_info$method)) user_SPR_weights_info$method = "equal"
        if(is.null(user_SPR_weights_info$weight_years)) user_SPR_weights_info$weight_years = 1
        if(is.null(user_SPR_weights_info$index_pointer)) user_SPR_weights_info$index_pointer = NULL
        em_input <- update_SPR_weights(em_input, 
                                       method = user_SPR_weights$weights_method,
                                       weight_years = user_SPR_weights_info$weight_years, 
                                       index_pointer = user_SPR_weights_info$index_pointer
                                       )
      }
      
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
      
      if (assess_interval != 0) {
        
        cat("\nNow checking convergence of assessment model...\n")
        
        conv <- check_conv(em)$conv
        pdHess <- check_conv(em)$pdHess
        if (conv & pdHess) cat("\nAssessment model is converged.\n") else warnings("\nAssessment model is not converged!\n")
        
        cat("\nNow generating catch advice...\n")
        
        advice <- advice_fn(em = em,
                            pro.yr = assess_interval,
                            hcr = hcr,
                            proj.opts = proj.opts)
        
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
        
        if(!is.null(implementation_error)) {
          cat("\nNow generating implementation error on catch advice...\n")
          method = implementation_error$method
          catch_mean = implementation_error$mean
          catch_cv = implementation_error$cv
          catch_sd = implementation_error$sd
          catch_min = implementation_error$min
          catch_max = implementation_error$max
          constant_value = implementation_error$constant_value
          real_catch = add_implementation_error(catch_advice = advice,
                                                method = method,
                                                mean = catch_mean,
                                                cv = catch_cv,
                                                sd = catch_sd,
                                                min = catch_min,
                                                max = catch_max,
                                                constant_value = constant_value,
                                                seed = seed)
          cat("\nRealized catch is...\n")
          print(real_catch)
          interval.info <- list(catch = real_catch, years = y + 1:assess_interval)
        } else {
          interval.info <- list(catch = advice, years = y + 1:assess_interval)
          real_catch = advice
        }
        
        cat("\nNow calculating F at age in the OM given the catch advice...\n")
        om <- update_om_fn(om, interval.info, seed = seed, random = random, method = "nlminb", by_fleet = by_fleet, do.brps = do.brps)
      } else {
        cat("\nNow performing simulation-estimation experiments...\n")
        conv = NULL
        pdHess = NULL
        advice = NULL
      }

      em_list[[i]] <- em$rep
      par.est[[i]] <- as.list(em$sdrep, "Estimate")
      par.se[[i]] <- as.list(em$sdrep, "Std. Error")
      adrep.est[[i]] <- as.list(em$sdrep, "Estimate", report = TRUE)
      adrep.se[[i]] <- as.list(em$sdrep, "Std. Error", report = TRUE)
      opt_list[[i]] <- em$opt
      converge_list[[i]] <- conv + pdHess
      catch_advice[[i]] <- advice
      catch_realized[[i]] <- real_catch
      
      if (save.sdrep) {
        em_full[[i]] <- em
      } else {
        if (y == assess_years[length(assess_years)]) {
          em_full[[1]] <- em
        }
        if (!save.last.em) em_full[[1]] <- list()
      }
      
      em_input_list[[i]] <- em_input
    }
  }
  
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  cat("Please ignore Warning in check_projF(proj_mod).")
  cat("\nTotal Runtime = ", time.taken,"\n")
  
  return(list(om = om, em_list = em_list, par.est = par.est, par.se = par.se, 
              adrep.est = adrep.est, adrep.se = adrep.se, opt_list = opt_list, 
              converge_list = converge_list, catch_advice = catch_advice, catch_realized = catch_realized, 
              em_full = em_full, em_input = em_input_list, runtime = time.taken, seed.save = seed))
}
