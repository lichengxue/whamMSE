#' Function to perform management strategy evaluation
#'
#' A wrapper function to step through the feedback period, update the operating model,
#' refit the estimation model, and generate catch advice.
#'
#' @param om Operating model including pseudo data
#' @param random A vector of processes that are treated as random effects in the operating model
#' @param M_om Natural mortality configuration in the operating model
#' @param sel_om Selectivity configuration in the operating model
#' @param NAA_re_om Numbers-at-age (NAA) configuration in the operating model
#' @param move_om Movement configuration in the operating model
#' @param age_comp_om Likelihood distribution of age composition data in the operating model
#' @param M_em Natural mortality configuration in the assessment model
#' @param sel_em Selectivity configuration in the assessment model
#' @param NAA_re_em Numbers-at-age (NAA) configuration in the assessment model
#' @param move_em Movement configuration in the assessment model
#' @param em.opt List of options for the assessment model
#' @param age_comp_em Likelihood distribution of age composition data in the assessment model
#' @param assess_years Year in which the assessment is conducted
#' @param assess_interval Assessment interval used in the MSE feedback loop
#' @param base_years Years used in the burn-in period
#' @param year.use Number of years used in the assessment model
#' @param hcr.type Type of harvest control rule
#' @param hcr.opts A list of HCR information, only used if hcr.type = 3
#' @param do.retro T/F Do retrospective analysis? Default = TRUE
#' @param do.osa T/F Calculate one-step-ahead (OSA) residuals? Default = TRUE
#' @param seed Seed used for generating data
#' @param save.sdrep T/F Save every assessment model (memory intensive)
#'   \itemize{
#'     \item \code{"TRUE"} Save every assessment model 
#'     \item \code{"FALSE"} Save the last assessment model (default)
#'     }
#'        
#' @return a list of model output
#' @export
#' @seealso \code{\link{make_em_input}}, \code{\link{update_om_fn}}, \code{\link{advice_fn}}

loop_through_fn <- function(om, random, M_om, sel_om, NAA_re_om, 
                            move_om, age_comp_om = "multinomial", M_em, sel_em,
                            NAA_re_em, move_em, em.opt = NULL, age_comp_em = "multinomial",
                            assess_years = NULL, assess_interval = NULL, base_years = NULL,
                            year.use = 30, hcr.type = 1, hcr.opts = NULL, do.retro = FALSE,
                            do.osa = FALSE, seed = 123, save.sdrep = FALSE) {
  
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
      
      em_input <- make_em_input(om = om, M_em = M_em, sel_em = sel_em, NAA_re_em = NAA_re_em,
                                move_em = move_em, em.opt = em.opt, em_years = em.years,
                                year.use = year.use, age_comp = age_comp_em)
      
      n_stocks <- om$input$data$n_stocks
      
      if (em.opt$separate.em.type == 1) {
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
          converge_list[[i]][[s]] <- sum(conv, pdHess)
          catch_advice[[i]] <- advice
          
          if (save.sdrep) {
            em_full[[i]][[s]] <- em[[s]]
          } else {
            if (y == assess_years[length(assess_years)]) {
              em_full[[1]][[s]] <- em[[s]]
            }
          }
        }
      } else if (em.opt$separate.em.type == 2) {
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
        }
        
      } else if (em.opt$separate.em.type == 3) {
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
        }
      }
    }
  } else {
    for (y in assess_years) {
      cat(paste0("\n-----\nStock Assessment in Year ", y, "\n"))
      i <- which(assess_years == y)
      em.years <- base_years[1]:y
      
      em_input <- make_em_input(om = om, M_em = M_em, sel_em = sel_em, NAA_re_em = NAA_re_em,
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
      }
    }
  }
  
  return(list(om = om, em_list = em_list, par.est = par.est, par.se = par.se,
              adrep.est = adrep.est, adrep.se = adrep.se, opt_list = opt_list,
              converge_list = converge_list, catch_advice = catch_advice, em_full = em_full))
}
