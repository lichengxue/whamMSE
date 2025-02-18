#' Generate Catch Advice for Management Strategy Evaluation
#'
#' A function to generate catch advice based on different harvest control rules (HCRs) for use in management strategy evaluation.
#'
#' @param em An estimation model fitted using WHAM.
#' @param pro.yr Integer. The number of years for projection.
#' @param hcr A list specifying the harvest control rule (HCR) type and options. It should contain:
#'   \describe{
#'     \item{\code{hcr.type}}{Integer. The type of harvest control rule to use. Options include:}
#'     \itemize{
#'       \item \code{1}: Annual projected catch based on 75\% of F40\% (default).
#'       \item \code{2}: Constant catch based on 75\% of F40\%.
#'       \item \code{3}: "Hockey stick" approach based on stock status.
#'     }
#'     \item{\code{hcr.opts}}{A list of additional options for the harvest control rule, including:}
#'     \describe{
#'       \item{\code{percentFXSPR}}{Numeric. The percentage of F\_XSPR to use for calculating catch in projections (default = 75).}
#'       \item{\code{percentSPR}}{Numeric. The X\% SPR used to calculate reference points (default = 40).}
#'       \item{\code{avg.yrs}}{Integer. The number of years over which mean values of MAA, WAA, M, and F will be averaged in projections (default = 5).}
#'       \item{\code{max_percent}}{Numeric. Maximum percentage of F\_XSPR to use in projections (default = 75).}
#'       \item{\code{min_percent}}{Numeric. Minimum percentage of F\_XSPR to use in projections (default = 0.01).}
#'       \item{\code{BThresh_up}}{Numeric. Upper bound of overfished biomass threshold (default = 0.5).}
#'       \item{\code{BThresh_low}}{Numeric. Lower bound of overfished biomass threshold (default = 0.1).}
#'     }
#'   }
#'
#' @return A matrix containing the projected catch advice for the specified number of projection years.
#'
#' @export
#'
#' @seealso \code{\link{project_wham}}
#'
#' @examples
#' \dontrun{
#' data <- generate_basic_info(n_stocks = 2, n_regions = 2, n_indices = 2, n_fleets = 2, base.years = 2003:2022)
#' NAA_re <- list(N1_model = c("equilibrium", "equilibrium"), sigma = c("rec", "rec"), cor = c("iid", "iid"), recruit_model = 2)
#' input <- prepare_wham_input(basic_info = data, NAA_re = NAA_re)
#' mod <- fit_wham(input, do.fit = FALSE)
#' input$data <- mod$simulate(complete = TRUE)
#' mod <- fit_wham(input, do.osa = FALSE, do.retro = FALSE, MakeADFun.silent = FALSE)
#' catch <- advice_fn(mod, pro.yr = 3, hcr = list(hcr.type = 1))
#' }

advice_fn <- function(em, pro.yr = assess.interval, hcr) {
  
  hcr.type = ifelse(is.null(hcr$hcr.type), 1, hcr$hcr.type) 
  hcr.opts = hcr$hcr.opts
  
  cat(paste0("\nHarvest Control Rule type ", hcr.type, "\n"))
  
  if(is.null(hcr.opts$percentFXSPR)) {
    percentFXSPR = 75
  } else {
    percentFXSPR = hcr.opts$percentFXSPR
  }
  
  if(is.null(hcr.opts$percentSPR)) {
    percentSPR = 40
  } else {
    percentSPR = hcr.opts$percentSPR
  }
  
  if(!is.null(hcr.opts$avg.yrs) & length(hcr.opts$avg.yrs) > length(em$years)) avg.yrs = length(em$years)
  
  if(is.null(hcr.opts$avg.yrs)) avg.yrs = 5

  proj_opts <- list(
    n.yrs = pro.yr,
    use.FXSPR = TRUE,
    avg.yrs = tail(em$years, avg.yrs),
    percentFXSPR = percentFXSPR,
    percentSPR = percentSPR
  )
  
  if (hcr.type %in% 1:2) em_proj <- project_wham(em, proj.opts = proj_opts, MakeADFun.silent = TRUE) # Projected version of the em
  
  if (hcr.type == 1) {
    advice <- em_proj$rep$pred_catch[length(em_proj$years) + 1:pro.yr,]
  } 
  if (hcr.type == 2) {
    advice <- em_proj$rep$pred_catch[length(em_proj$years) + 1:pro.yr,]
    if (nrow(advice) != 1) {
      advice <- colMeans(advice) # Mean of the projected catch over the next 5 years fishing at F40
    }
    advice <- matrix(rep(advice, pro.yr), ncol = length(advice), byrow = TRUE)
  } 
  if (hcr.type == 3) {
    
    if (is.null(hcr.opts$max_percent)) {
      max_percent = 75
    } else {
      max_percent = hcr.opts$max_percent
    }
    
    if (is.null(hcr.opts$min_percent)) {
      min_percent = 0.01
    } else {
      min_percent = hcr.opts$min_percent
    }
    
    if (is.null(hcr.opts$BThresh_up)) {
      BThresh_up = 0.5
    } else {
      BThresh_up = hcr.opts$BThresh_up
    }
    
    if (is.null(hcr.opts$BThresh_low)) {
      BThresh_low = 0.1
    } else {
      BThresh_low = hcr.opts$BThresh_low
    }
    
    if (ncol(em$rep$log_SSB_FXSPR) == ncol(em$rep$SSB) + 1) {
      print("Global SPR is calculated in the model")
      SSB40 <- exp(tail(em$rep$log_SSB_FXSPR, 1))[, ncol(em$rep$log_SSB_FXSPR)]
      SSB_t <- sum(tail(em$rep$SSB, 1))
      
      ratio <- SSB_t / SSB40
      
      if (ratio >= BThresh_up) {
        cat(paste0("SSB_t/SSB40 = ", round(ratio, 3), "\n"))
        proj_opts$percentFXSPR <- max_percent
      } else if (ratio < BThresh_up & ratio > BThresh_low) {
        cat(paste0("SSB_t/SSB40 = ", round(ratio, 3), "\n"))
        slope <- (max_percent - min_percent) / (BThresh_up - BThresh_low)
        percent <- slope * (ratio - BThresh_low) + min_percent
        proj_opts$percentFXSPR <- percent
      } else if (ratio <= BThresh_low) {
        cat(paste0("SSB_t/SSB40 = ", round(ratio, 3), "\n"))
        proj_opts$percentFXSPR <- min_percent
      }
      
      em_proj <- project_wham(em, proj.opts = proj_opts, MakeADFun.silent = TRUE) # Projected version of the em
      advice <- em_proj$rep$pred_catch[length(em_proj$years) + 1:pro.yr,]
    }
  }
  
  return(advice)
}
