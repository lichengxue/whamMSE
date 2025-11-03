#' Generate Catch Advice for Management Strategy Evaluation
#'
#' A function to generate catch advice based on different harvest control rules (HCRs)
#' for use in management strategy evaluation (MSE).
#'
#' @param em A fitted WHAM estimation model.
#' @param pro.yr Integer. Number of years for projection.
#' @param hcr List. Harvest control rule:
#'   \describe{
#'     \item{\code{hcr.type}}{Integer (1: F\_XSPR, 2: constant catch, 3: hockey stick).}
#'     \item{\code{hcr.opts}}{List of options, including:
#'       \itemize{
#'         \item \code{use_FXSPR}: Logical. Use F\%SPR-based projections (default: TRUE).
#'         \item \code{percentFXSPR}: Numeric. % of F\%SPR to apply (default: 75).
#'         \item \code{use_FMSY}: Logical. Use FMSY-based projections (default: FALSE).
#'         \item \code{percentFMSY}: Numeric. % of FMSY to apply (default: 75).
#'         \item \code{avg_yrs}: Integer. Number of years to average over (default: 5).
#'         \item \code{cont.M.re}: Logical. Continue M random effects (default: FALSE).
#'         \item \code{cont.move.re}: Logical. Continue movement random effects (default: NULL if 1 region).
#'         \item \code{max_percent}: Numeric. Max % of F\%SPR or FMSY (default: 75).
#'         \item \code{min_percent}: Numeric. Min % of F\%SPR or FMSY (default: 0.01).
#'         \item \code{BThresh_up}: Numeric. Upper biomass threshold (default: 0.5).
#'         \item \code{BThresh_low}: Numeric. Lower biomass threshold (default: 0.1).
#'       }
#'     }
#'   }
#' @param ecov_em_opts List (optional). Options for projecting environmental
#'   covariates in the estimation model during catch advice. Allows explicit
#'   control of Ecov inputs in the projection period. The expected components include:
#'   \itemize{
#'     \item `$use_ecov_em` Logical. If `TRUE`, the EM uses user-specified or
#'       projected Ecov values instead of directly continuing the OM process.
#'     \item `$lag` Integer. Specifies the lag (in years) to align Ecov values when used in
#'       recruitment or other population processes. Must be provided if `$use_ecov_em = TRUE`.
#'     \item `$period` Integer vector (optional). If provided, specifies the projection
#'       years (indices relative to the Ecov time series) to override with user-defined
#'       values. If `NULL`, then the entire projection period is replaced.
#'   }
#'   If `NULL`, Ecov values for projections are inherited directly from the OM/EM state without override.
#' @return A matrix containing the projected catch advice for \code{pro.yr} years.
#'
#' @examples
#' # Example: Define a harvest control rule using F_40%SPR with 100% target
#' hcr <- list(
#'   hcr.type = 1,
#'   hcr.opts = list(
#'     use_FXSPR = TRUE,
#'     percentFXSPR = 100,  # Use 100% of F_40%SPR
#'     avg_yrs = 5,         # Average over the last 5 years
#'     cont.M.re = FALSE,   # Do not continue M random effects
#'     cont.move.re = FALSE # Do not continue movement random effects
#'   )
#' )
#'
#' # Run advice_fn with a WHAM model `em`
#' # advice <- advice_fn(em, pro.yr = 5, hcr = hcr)
#'
#' @seealso \code{\link{project_wham}}
#' @export
advice_fn <- function(em, pro.yr = assess.interval, hcr = NULL, ecov_em_opts = NULL) {
  
  # Ensure hcr and hcr.opts are always lists
  if (is.null(hcr)) {
    hcr <- list(hcr.type = 1, hcr.opts = list())
  }
  
  hcr.type <- ifelse(is.null(hcr$hcr.type), 1, hcr$hcr.type)
  hcr.opts <- if (is.null(hcr$hcr.opts)) list() else hcr$hcr.opts
  
  cat(paste0("\nHarvest Control Rule type ", hcr.type, "\n"))
  
  # Set defaults safely
  use_FXSPR <- ifelse(is.null(hcr.opts$use_FXSPR), TRUE, hcr.opts$use_FXSPR)
  percentFXSPR <- ifelse(is.null(hcr.opts$percentFXSPR), 75, hcr.opts$percentFXSPR)
  use_FMSY <- ifelse(is.null(hcr.opts$use_FMSY), FALSE, hcr.opts$use_FMSY)
  percentFMSY <- ifelse(is.null(hcr.opts$percentFMSY), 75, hcr.opts$percentFMSY)
  avg_yrs <- ifelse(is.null(hcr.opts$avg_yrs), 5, hcr.opts$avg_yrs)
  
  # Adjust avg_yrs if needed
  if (!is.null(hcr.opts$avg_yrs) && hcr.opts$avg_yrs > length(em$years)) {
    avg_yrs <- length(em$years)
  }
  
  cont.M.re <- ifelse(is.null(hcr.opts$cont.M.re), FALSE, hcr.opts$cont.M.re)
  
  if (is.null(hcr.opts$cont.move.re) || em$input$data$n_regions == 1) {
    cont.move.re <- NULL
  } else {
    cont.move.re <- hcr.opts$cont.move.re
  }

  proj_ecov = NULL
  if(!is.null(ecov_em_opts) && ecov_em_opts$use_ecov_em) {
    if(!is.null(ecov_em_opts$lag)) {
      id = length(em$input$years)-ecov_em_opts$lag
      proj_ecov = em$rep$Ecov_x[id:(id+pro.yr-ecov_em_opts$lag),, drop = FALSE]
    } else {
      warnings("Must specify lag when use_ecov_em = TRUE")
    }
  }
  
  # Set up projection options
  proj_opts <- list(
    n.yrs = pro.yr,
    use.FXSPR = use_FXSPR,
    percentFXSPR = percentFXSPR,
    use.FMSY = use_FMSY,
    percentFMSY = percentFMSY,
    avg.yrs = tail(em$years, avg_yrs),
    cont.M.re = cont.M.re,
    cont.move.re = cont.move.re,
    proj.ecov = proj_ecov
  )
  
  # --- HCR type 1 & 2 ---
  if (hcr.type %in% 1:2) {
    em_proj <- project_wham(em, proj.opts = proj_opts, MakeADFun.silent = TRUE)
  }
  
  if (hcr.type == 1) {
    advice <- em_proj$rep$pred_catch[
      (length(em_proj$years) + 1):(length(em_proj$years) + pro.yr), 
    ]
  } 
  
  if (hcr.type == 2) {
    advice <- em_proj$rep$pred_catch[
      (length(em_proj$years) + 1):(length(em_proj$years) + pro.yr), 
    ]
    if (nrow(advice) != 1) {
      advice <- colMeans(advice)
    }
    advice <- matrix(rep(advice, pro.yr), ncol = length(advice), byrow = TRUE)
  }
  
  # --- HCR type 3 ---
  if (hcr.type == 3) {
    
    max_percent <- ifelse(is.null(hcr.opts$max_percent), 75, hcr.opts$max_percent)
    min_percent <- ifelse(is.null(hcr.opts$min_percent), 0.01, hcr.opts$min_percent)
    BThresh_up <- ifelse(is.null(hcr.opts$BThresh_up), 0.5, hcr.opts$BThresh_up)
    BThresh_low <- ifelse(is.null(hcr.opts$BThresh_low), 0.1, hcr.opts$BThresh_low)
    
    # --- FXSPR path ---
    if (use_FXSPR) {
      if (ncol(em$rep$log_SSB_FXSPR) == ncol(em$rep$SSB) + 1) {
        cat("Global SPR is calculated in the model\n")
        SSB_x <- exp(tail(em$rep$log_SSB_FXSPR, 1))[, ncol(em$rep$log_SSB_FXSPR)]
        SSB_t <- sum(tail(em$rep$SSB, 1))
        
        ratio <- SSB_t / SSB_x
        cat(paste0("SSB_t/SSB_x = ", round(ratio, 3), "\n"))
        
        if (ratio >= BThresh_up) {
          proj_opts$percentFXSPR <- max_percent
        } else if (ratio < BThresh_up & ratio > BThresh_low) {
          slope <- (max_percent - min_percent) / (BThresh_up - BThresh_low)
          percent <- slope * (ratio - BThresh_low) + min_percent
          proj_opts$percentFXSPR <- percent
        } else if (ratio <= BThresh_low) {
          proj_opts$percentFXSPR <- min_percent
        }
        
        em_proj <- project_wham(em, proj.opts = proj_opts, MakeADFun.silent = TRUE)
        advice <- em_proj$rep$pred_catch[
          (length(em_proj$years) + 1):(length(em_proj$years) + pro.yr), 
        ]
      }
    }
    
    # --- FMSY path ---
    if (use_FMSY) {
      if (ncol(em$rep$log_SSB_MSY) == ncol(em$rep$SSB) + 1) {
        cat("Global SPR is calculated in the model\n")
        SSB_x <- exp(tail(em$rep$log_SSB_MSY, 1))[, ncol(em$rep$log_SSB_MSY)]
        SSB_t <- sum(tail(em$rep$SSB, 1))
        
        ratio <- SSB_t / SSB_x
        cat(paste0("SSB_t/SSB_x = ", round(ratio, 3), "\n"))
        
        if (ratio >= BThresh_up) {
          proj_opts$percentFMSY <- max_percent
        } else if (ratio < BThresh_up & ratio > BThresh_low) {
          slope <- (max_percent - min_percent) / (BThresh_up - BThresh_low)
          percent <- slope * (ratio - BThresh_low) + min_percent
          proj_opts$percentFMSY <- percent
        } else if (ratio <= BThresh_low) {
          proj_opts$percentFMSY <- min_percent
        }
        
        em_proj <- project_wham(em, proj.opts = proj_opts, MakeADFun.silent = TRUE)
        advice <- em_proj$rep$pred_catch[
          (length(em_proj$years) + 1):(length(em_proj$years) + pro.yr), 
        ]
      }
    }
  }
  
  # --- Print the projection options nicely ---
  cat("Projection Options:\n")
  for (proj_name in names(proj_opts)) {
    cat(sprintf(" %s: %s\n", proj_name, toString(proj_opts[[proj_name]])))
  }
  
  return(advice)
}
