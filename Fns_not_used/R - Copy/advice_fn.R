#' Generate Catch Advice for Management Strategy Evaluation
#'
#' A function to generate catch advice based on different harvest control rules (HCRs) for use in management strategy evaluation.
#'
#' @param em An estimation model fitted using WHAM.
#' @param pro.yr Integer. The number of years for projection.
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
#' @return A matrix containing the projected catch advice for the specified number of projection years.
#'
#' @export
#'
#' @seealso \code{\link{project_wham}}
#' 
advice_fn <- function(em, pro.yr = assess.interval, hcr) {
  
  hcr.type = ifelse(is.null(hcr$hcr.type), 1, hcr$hcr.type) 
  hcr.opts = hcr$hcr.opts
  
  cat(paste0("\nHarvest Control Rule type ", hcr.type, "\n"))
  
  if(is.null(hcr.opts$use.FXSPR)) {
    use.FXSPR = TRUE
    if(is.null(hcr.opts$percentFXSPR)) {
      percentFXSPR = 75
    } else {
      percentFXSPR = hcr.opts$percentFXSPR
    }
  }
  
  if(is.null(hcr.opts$use.FMSY)) {
    use.FMSY = FALSE
  } else {
    use.FMSY = hcr.opts$use.FMSY
  }
  
  if(is.null(hcr.opts$percentFMSY)) {
    percentFMSY = 75
  } else {
    percentFMSY = hcr.opts$percentFMSY
  }
  
  if(!is.null(hcr.opts$avg.yrs) & length(hcr.opts$avg.yrs) > length(em$years)) avg.yrs = length(em$years)
  
  if(is.null(hcr.opts$avg.yrs)) avg.yrs = 5
  
  if(is.null(hcr.opts$cont.M.re)) {
    cont.M.re = FALSE
  } else {
    cont.M.re = hcr.opts$cont.M.re
  }
  
  if(is.null(hcr.opts$cont.move.re)) {
    cont.move.re = NULL
  } else {
    if(em$input$data$n_regions == 1) {
      cont.move.re = NULL
    } else {
      cont.move.re = hcr.opts$cont.move.re
    }
  }
      
  proj_opts <- list(
    n.yrs = pro.yr,
    use.FXSPR = use.FXSPR,
    percentFXSPR = percentFXSPR,
    use.FMSY = use.FMSY,
    percentFMSY = percentFMSY,
    avg.yrs = tail(em$years, avg.yrs),
    cont.M.re = cont.M.re,
    cont.move.re = cont.move.re
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
    
    if (use.FXSPR) {
      if (ncol(em$rep$log_SSB_FXSPR) == ncol(em$rep$SSB) + 1) {
        print("Global SPR is calculated in the model")
        SSB_x <- exp(tail(em$rep$log_SSB_FXSPR, 1))[, ncol(em$rep$log_SSB_FXSPR)]
        SSB_t <- sum(tail(em$rep$SSB, 1))
        
        ratio <- SSB_t / SSB_x
        
        if (ratio >= BThresh_up) {
          cat(paste0("SSB_t/SSB_x = ", round(ratio, 3), "\n"))
          proj_opts$percentFXSPR <- max_percent
        } else if (ratio < BThresh_up & ratio > BThresh_low) {
          cat(paste0("SSB_t/SSB_x = ", round(ratio, 3), "\n"))
          slope <- (max_percent - min_percent) / (BThresh_up - BThresh_low)
          percent <- slope * (ratio - BThresh_low) + min_percent
          proj_opts$percentFXSPR <- percent
        } else if (ratio <= BThresh_low) {
          cat(paste0("SSB_t/SSB_x = ", round(ratio, 3), "\n"))
          proj_opts$percentFXSPR <- min_percent
        }
        
        em_proj <- project_wham(em, proj.opts = proj_opts, MakeADFun.silent = TRUE) # Projected version of the em
        advice <- em_proj$rep$pred_catch[length(em_proj$years) + 1:pro.yr,]
      }
    }
    
    if (use.FMSY) {
      if (ncol(em$rep$log_SSB_MSY) == ncol(em$rep$SSB) + 1) {
        print("Global SPR is calculated in the model")
        SSB_x <- exp(tail(em$rep$log_SSB_MSY, 1))[, ncol(em$rep$log_SSB_MSY)]
        SSB_t <- sum(tail(em$rep$SSB, 1))
        
        ratio <- SSB_t / SSB_x
        
        if (ratio >= BThresh_up) {
          cat(paste0("SSB_t/SSB_x = ", round(ratio, 3), "\n"))
          proj_opts$percentFMSY <- max_percent
        } else if (ratio < BThresh_up & ratio > BThresh_low) {
          cat(paste0("SSB_t/SSB_x = ", round(ratio, 3), "\n"))
          slope <- (max_percent - min_percent) / (BThresh_up - BThresh_low)
          percent <- slope * (ratio - BThresh_low) + min_percent
          proj_opts$percentFMSY <- percent
        } else if (ratio <= BThresh_low) {
          cat(paste0("SSB_t/SSB_x = ", round(ratio, 3), "\n"))
          proj_opts$percentFMSY <- min_percent
        }
        
        em_proj <- project_wham(em, proj.opts = proj_opts, MakeADFun.silent = TRUE) # Projected version of the em
        advice <- em_proj$rep$pred_catch[length(em_proj$years) + 1:pro.yr,]
      }
    }
  }
  
  return(advice)
}
