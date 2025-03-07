#' Perform Management Strategy Evaluation (MSE)
#' 
#' A wrapper function to iterate through the feedback period in a management strategy evaluation (MSE), 
#' updating the operating model (OM), refitting the estimation model (EM), and generating catch advice.
#' 
#' @param om A fitted operating model containing pseudo data.
#' @param em_info A list containing information used to generate the estimation model (default = NULL).
#' @param random Character vector of processes treated as random effects in the operating model (default = "log_NAA").
#' @param M_em Configuration of natural mortality in the estimation model.
#' @param sel_em Configuration of selectivity in the estimation model.
#' @param NAA_re_em Configuration of numbers-at-age (NAA) in the estimation model.
#' @param move_em Configuration of movement in the estimation model.
#' @param age_comp_em Character. The likelihood distribution of age composition data in the estimation model. Options include:
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
#'     \item{\code{separate.em}}{Logical. If TRUE, separate estimation models are used for different stocks (default = FALSE).}
#'     \item{\code{separate.em.type}}{Integer. Used only if \code{separate.em = TRUE}:}
#'       \itemize{
#'         \item 1 = Panmictic (spatially-aggregated model).
#'         \item 2 = Fleets-as-areas model.
#'         \item 3 = Separate assessment models for each region (n = number of regions).
#'       }
#'     \item{\code{do.move}}{Logical. If TRUE, movement is included (default = FALSE).}
#'     \item{\code{est.move}}{Logical. If TRUE, movement rates are estimated (default = FALSE).}
#'   }
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
#' @param ind_em Vector. Indices specifying the years for which the estimation model should use data.
#' @param filter_indices Integer (0/1) vector (optional). User-specified which indices are excluded from the assessment model. For example, c(1,0,1,1) indicates Index 1 (include) and 2 (exclude) in region 1, Index 3 (include) and 4 (include) in region 2
#' @param assess_years Vector of years when assessments are conducted.
#' @param assess_interval Integer. The interval between stock assessments in the MSE feedback loop.
#' @param base_years Vector of years used in the burn-in period.
#' @param year.use Integer. Number of years included in the estimation model (default = 20).
#' @param add.years Logical. Whether or not using entire time series of data in the assessment model (default = FALSE).
#' @param by_fleet Logical. Whether calculate fleet-specific fishing mortality or global fishing mortality (default = TRUE).
#' @param do.brps Logical. Whether or not update reference point in the OM (default = FALSE).
#' @param hcr List containing harvest control rule (HCR) settings:
#'   \describe{
#'     \item{\code{hcr.type}}{Integer. The type of harvest control rule:}
#'       \itemize{
#'         \item \code{1} Annual projected catch based on 75\% of F40\% (default).
#'         \item \code{2} Constant catch based on 75\% of F40\%.
#'         \item \code{3} "Hockey stick" catch based on stock status.
#'       }
#'     \item{\code{hcr.opts}}{(Only for \code{hcr.type = 3}) A list containing additional HCR options:}
#'       \describe{
#'         \item{\code{max_percent}}{Maximum percentage of F\_XSPR used for catch projections (default = 75).}
#'         \item{\code{min_percent}}{Minimum percentage of F\_XSPR used for catch projections (default = 0.01).}
#'         \item{\code{BThresh_up}}{Upper bound of overfished biomass threshold (default = 0.5).}
#'         \item{\code{BThresh_low}}{Lower bound of overfished biomass threshold (default = 0.1).}
#'       }
#'   }
#' @param catch_alloc List. Contains specifications for catch allocation:
#'   \itemize{
#'     \item `$weight_type` - Integer (1–4) indicating the type of weighting used.
#'       \itemize{
#'         \item `1` - Equal weighting across fleets.
#'         \item `2` - Weighting based on past catch data.
#'         \item `3` - Weighting based on past survey index data.
#'         \item `4` - User-specified weights.
#'       }
#'     \item `$method` - String. Specifies the catch allocation method:
#'       \itemize{
#'         \item `"equal"` - Use equal weights, split total catch into fleet-specific catches (use when weight_type = 1).
#'         \item `"fleet_equal"` - Use equal weights, split total fleet catch into fleet-specific catches (use when weight_type = 1).
#'         \item `"fleet_region"` - Uses regional catch totals to compute weights (use when weight_type = 2).
#'         \item `"fleet_gear"` - Uses gear-type catch totals to compute weights (use when weight_type = 2).
#'         \item `"fleet_combined"` - First allocates to gear-specific total catch, then splits into regions (use when weight_type = 2).
#'         \item `"fleet_catch"` - Uses fleet-specific catch to compute weights (use when weight_type = 2).
#'         \item `"index_equal"` - Uses survey-based regional weighting, then assigns equally among fleets in the same region (use when weight_type = 3).
#'         \item `"index_gear"` - Uses survey-based regional weighting, then assigns based on gear-specific weights (use when weight_type = 3).
#'         \item `"index_multiple"` - Uses multiple survey-based regional weighting, then assigns equally among fleets in the same region (use when weight_type = 3).
#'         \item `"user_defined"` - Uses manually specified weights for each region or each fleet (use when weight_type = 4).
#'       }
#'     \item `$user_weights` - Numeric vector (`n_regions` or `n_fleets`). Optional. User-defined weights summing to 1 (use when weight_type = 4).
#'     \item `$weight_years` - Integer. Number of years to average for calculating historical catch weights.
#'     \item `$survey_pointer` - Integer/Vector. Specifies which survey index type to use for weighting (use when weight_type = 3).
#'   }
#' @param do.retro Logical. If TRUE, performs retrospective analysis for the assessment model (default = TRUE).
#' @param do.osa Logical. If TRUE, calculates one-step-ahead (OSA) residuals for the assessment model (default = TRUE).
#' @param do.brps Logical. If TRUE, calculates reference points in the operating model (default = FALSE).
#' @param seed Integer. The random seed used for stochastic processes.
#' @param save.sdrep Logical. If TRUE, saves the results of every assessment model (memory-intensive).
#' @param save.last.em Logical. If TRUE, saves only the last estimation model (default = FALSE).
#' @param by_fleet Logical. If TRUE, estimates F separately for each fleet. If FALSE, estimates a single global F (default = TRUE).
#' @param FXSPR_init Numeric. change initial F for estimating reference point in the assessment model (default = NULL).
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
                            age_comp_em = NULL, 
                            em.opt = list(separate.em = TRUE, separate.em.type = 1,
                                          do.move = FALSE, est.move = FALSE), 
                            global_waa = TRUE,
                            aggregate_catch_info = NULL,
                            aggregate_index_info = NULL,
                            filter_indices = NULL,
                            assess_years = NULL, 
                            assess_interval = NULL, 
                            base_years = NULL, 
                            year.use = 20, 
                            add.years = FALSE,
                            hcr = list(hcr.type = 1, hcr.opts = NULL),
                            catch_alloc = list(weight_type = 1, method = "equal", user_weights = NULL, weight_years = 1),
                            do.retro = FALSE, 
                            do.osa = FALSE, 
                            do.brps = FALSE,
                            seed = 123, 
                            save.sdrep = FALSE, 
                            save.last.em = FALSE,
                            by_fleet = TRUE,
                            FXSPR_init = NULL) {
  
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
        
      # Note em_info$par_inputs$user_waa for fleet 1 is wrong!
      em_input <- make_em_input(om = om, em_info = em_info, M_em = M_em, sel_em = sel_em,
                                NAA_re_em = NAA_re_em, move_em = move_em, em.opt = em.opt,
                                em_years = em.years, year.use = year.use, age_comp_em = age_comp_em,
                                aggregate_catch_info = aggregate_catch_info,
                                aggregate_index_info = aggregate_index_info,
                                filter_indices = filter_indices,
                                global_waa) # turn off for panmictic model for now
      
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
                                global_waa) # turn off for panmictic model for now
      
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
      
      if(is.vector(advice)) advice <- as.matrix(t(advice))
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
