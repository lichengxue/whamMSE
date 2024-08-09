#' Function to perform management strategy evaluation
#'
#' A wrapper function to step through the feedback period, update the operating model,
#' refit the estimation model, and generate catch advice.
#'
#' @param om Operating model including pseudo data
#' @param em_info A list of information used to generate the estimation model (default = NULL)
#' @param random A vector of processes that are treated as random effects in the operating model (default = "log_NAA")
#' @param M_em Natural mortality configuration in the assessment model
#' @param sel_em Selectivity configuration in the assessment model
#' @param NAA_re_em Numbers-at-age (NAA) configuration in the assessment model
#' @param move_em Movement configuration in the assessment model
#' @param age_comp_em Likelihood distribution of age composition data in the assessment model
#'   \itemize{
#'     \item \code{"multinomial"} Default
#'     \item \code{"dir-mult"}
#'     \item \code{"dirichlet-miss0"}
#'     \item \code{"dirichlet-pool0"}
#'     \item \code{"dir-mult"}
#'     \item \code{"logistic-normal-miss0"}    
#'     \item \code{"logistic-normal-ar1-miss0"}
#'     \item \code{"logistic-normal-pool0"}
#'     \item \code{"logistic-normal-01-infl"}
#'     \item \code{"logistic-normal-pool0"}
#'     \item \code{"logistic-normal-01-infl-2par"}
#'     \item \code{"mvtweedie"}
#'     \item \code{"dir-mult-linear"}
#'     }
#' @param em.opt List of options for the assessment model
#'   \itemize{
#'     \item \code{$separate.em} TRUE = No Global SPR, FALSE = Global SPR
#'     \item \code{$separate.em.type} only if separate.em = TRUE \cr
#'     {=1} panmictic (spatially-aggregated) \cr
#'     {=2} fleets-as-areas \cr
#'     {=3} n single assessment models (n = n_regions) \cr
#'     \item \code{$do.move} T/F movement is included (use if separate.em = FALSE)
#'     \item \code{$est.move} T/F movement rate is estimated (use if separate.em = FALSE)
#'     }
#' @param assess_years Year in which the assessment is conducted
#' @param assess_interval Assessment interval used in the MSE feedback loop
#' @param base_years Years used in the burn-in period
#' @param year.use Number of years used in the assessment model
#' @param hcr.type Type of harvest control rule
#'   \itemize{
#'     \item \code{"1"} Annual projected catch based on 75% of F40% (default)
#'     \item \code{"2"} Constant catch based on 75% of F40%
#'     \item \code{"3"} "Hockey stick" catch based on stock status
#'     }
#' @param hcr.opts only used if hcr.type = 3
#'   \itemize{
#'     \item \code{"max_percent"} maximum percent of F_XSPR to use for calculating catch in projections (default = 75)
#'     \item \code{"min_percent"} minimum percent of F_XSPR to use for calculating catch in projections (default = 0.01)
#'     \item \code{"BThresh_up"} Upper bound of overfished level (default = 0.5) 
#'     \item \code{"BThresh_low"} Lower bound of overfished level (default = 0.1)
#'     }
#' @param do.retro T/F Do retrospective analysis? Default = TRUE
#' @param do.osa T/F Calculate one-step-ahead (OSA) residuals? Default = TRUE
#' @param seed Seed used for generating data
#' @param save.sdrep T/F Save every assessment model (memory intensive)
#'   \itemize{
#'     \item \code{"TRUE"} Save every assessment model 
#'     \item \code{"FALSE"} Save the last assessment model (default)
#'     }
#' @param save.last.em T/F Save the last assessment model (Default = FALSE) (use if separate.em = FALSE)
#' @return a list of model output
#' @export
#' @seealso \code{\link{make_em_input}}, \code{\link{update_om_fn}}, \code{\link{advice_fn}}

loop_through_fn <- function(om, em_info = NULL, random = "log_NAA", M_em, sel_em, NAA_re_em, move_em, 
                            age_comp_em = "multinomial",em.opt = NULL, 
                            assess_years = NULL, assess_interval = NULL, base_years = NULL,
                            year.use = 30, hcr.type = 1, hcr.opts = NULL, do.retro = FALSE,
                            do.osa = FALSE, seed = 123, save.sdrep = FALSE, save.last.em = FALSE) {
  
  # Helper function to check convergence
  check_conv <- function(em) {
    conv <- as.logical(1 - em$opt$convergence)
    pdHess <- as.logical(if (em$na_sdrep == FALSE & !is.na(em$na_sdrep)) 1 else 0)
    if (!conv | !pdHess) warning("Assessment model is not converged!")
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
  
  if (em.opt$separate.em) {
    
    for (y in assess_years) {
      
      cat(paste0("\n-----\nStock Assessment in Year ", y, "\n"))
      
      i <- which(assess_years == y)
      
      em.years <- base_years[1]:y
      
      em_input <- make_em_input(om = om, em_info = em_info, M_em = M_em, sel_em = sel_em, NAA_re_em = NAA_re_em,
                                move_em = move_em, em.opt = em.opt, em_years = em.years,
                                year.use = year.use, age_comp = age_comp_em)
      
      n_stocks <- om$input$data$n_stocks
      
      if (em.opt$separate.em.type == 1) {
        
        em <- fit_wham(em_input, do.retro = FALSE, do.osa = FALSE, do.brps = TRUE, MakeADFun.silent = TRUE)
        
        conv <- check_conv(em)$conv
        pdHess <- check_conv(em)$pdHess
        
        advice <- advice_fn(em, pro.yr = assess_interval, hcr.type = hcr.type, hcr.opts = hcr.opts)
        
        # Weights are calculated based on the survey catch from the previous year
        weights <- om$input$data$agg_indices[which(om$years == y),]/sum(om$input$data$agg_indices[which(om$years == y),])
        
        catch_matrix <- matrix(advice, ncol = 1)
        
        apply_weights <- function(advice) {
          advice*weights
        }
        
        advice <- t(apply(catch_matrix, 1, apply_weights))
        
        colnames(advice) <- paste0("Region_", 1:om$input$data$n_fleets)
        rownames(advice) <- paste0("Year_", y + 1:assess_interval)
        
        cat("\n---------------------------\n")
        print(advice)
        cat("\n---------------------------\n")
        
        interval.info <- list(catch = advice, years = y + 1:assess_interval)
        om <- update_om_fn(om, interval.info, seed = seed, random = random)
        
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
          if(!save.last.em) em_full[[1]] <- list()
        }
      }
      
      if (em.opt$separate.em.type == 2) {
        
        em <- fit_wham(em_input, do.retro = FALSE, do.osa = FALSE, do.brps = TRUE, MakeADFun.silent = TRUE)
        
        conv <- check_conv(em)$conv
        pdHess <- check_conv(em)$pdHess
        
        advice <- advice_fn(em, pro.yr = assess_interval, hcr.type = hcr.type, hcr.opts = hcr.opts)
        
        colnames(advice) <- paste0("Region_", 1:om$input$data$n_fleets)
        rownames(advice) <- paste0("Year_", y + 1:assess_interval)
        
        cat("\n---------------------------\n")
        print(advice)
        cat("\n---------------------------\n")
        
        interval.info <- list(catch = advice, years = y + 1:assess_interval)
        om <- update_om_fn(om, interval.info, seed = seed, random = random)
        
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
          if(!save.last.em) em_full[[1]] <- list()
        }
        
      } 
      
      if (em.opt$separate.em.type == 3) {
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
        
        for (s in 1:n_stocks) {
          em[[s]] <- fit_wham(em_input[[s]], do.retro = FALSE, do.osa = FALSE, do.brps = TRUE, MakeADFun.silent = TRUE)
          convergence <- check_conv(em[[s]])
          conv[s] <- convergence$conv
          pdHess[s] <- convergence$pdHess
          tmp <- advice_fn(em[[s]], pro.yr = assess_interval, hcr.type = hcr.type, hcr.opts = hcr.opts)
          advice <- cbind(advice, tmp)
        }
        
        colnames(advice) <- paste0("Region_", 1:om$input$data$n_fleets)
        rownames(advice) <- paste0("Year_", assess_years[i] + 1:assess_interval)
        
        cat("\n---------------------------\n")
        print(advice)
        cat("\n---------------------------\n")
        
        # set the catch for the next assess_interval years
        interval.info <- list(catch = advice, years = assess_years[i] + 1:assess_interval)
        om <- update_om_fn(om, interval.info, seed = seed, random = random)
        
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
            if(!save.last.em) em_full[[1]][[s]] <- list()
          }
        }
      } 
      
    }
  } else {
    
    for (y in assess_years) {
      
      cat(paste0("\n-----\nStock Assessment in Year ", y, "\n"))
      
      i <- which(assess_years == y)
      
      em.years <- base_years[1]:y
      
      em_input <- make_em_input(om = om, em_info = em_info, M_em = M_em, sel_em = sel_em, NAA_re_em = NAA_re_em,
                                move_em = move_em, em.opt = em.opt, em_years = em.years,
                                year.use = year.use, age_comp = age_comp_em)
      
      if (em.opt$do.move) {
        
        if (em.opt$est.move) {
          
          em <- fit_wham(em_input, do.retro = FALSE, do.osa = FALSE, do.brps = TRUE, MakeADFun.silent = TRUE)
          
        } else {
          
          em_input <- fix_move(em_input)
          
          em <- fit_wham(em_input, do.retro = FALSE, do.osa = FALSE, do.brps = TRUE, MakeADFun.silent = TRUE)
          
        }
        
      } else {
        
        em <- fit_wham(em_input, do.retro = FALSE, do.osa = FALSE, do.brps = TRUE, MakeADFun.silent = TRUE)
        
      }
      
      conv <- check_conv(em)$conv
      pdHess <- check_conv(em)$pdHess
      
      advice <- advice_fn(em, pro.yr = assess_interval, hcr.type = hcr.type, hcr.opts = hcr.opts)
      
      colnames(advice) <- paste0("Region_", 1:om$input$data$n_fleets)
      rownames(advice) <- paste0("Year_", y + 1:assess_interval)
      
      cat("\n---------------------------\n")
      print(advice)
      cat("\n---------------------------\n")
      
      interval.info <- list(catch = advice, years = y + 1:assess_interval)
      om <- update_om_fn(om, interval.info, seed = seed, random = random)
      
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
        if(!save.last.em) em_full[[1]] <- list()
      }
    }
  }
  
  return(list(om = om, em_list = em_list, par.est = par.est, par.se = par.se,
              adrep.est = adrep.est, adrep.se = adrep.se, opt_list = opt_list,
              converge_list = converge_list, catch_advice = catch_advice, em_full = em_full))
}
