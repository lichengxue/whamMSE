#' Generate catch advice
#' 
#' A function to generate catch advice for management strategy evaluation. 
#' 
#' @param em Estimation model
#' @param pro.yr Number of years in projection
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
#'     
#' @return a list of catch advice
#'   
#' @export
#' 
#' @seealso \code{\link{project_wham}}
#' 
#' @examples
#' \dontrun{
#' data <- generate_basic_info(n_stocks = 2, n_regions = 2, n_indices = 2, n_fleets = 2, base.years = 2003:2022)
#' NAA_re <- list(N1_model=c("equilibrium","equilibrium"),sigma=c("rec","rec"),cor=c("iid","iid"),recruit_model = 2)
#' input <- prepare_wham_input(basic_info = data, NAA_re = NAA_re)
#' mod <- fit_wham(input, do.fit = FALSE)
#' input$data = mod$simulate(complete=TRUE)
#' mod <- fit_wham(input, do.osa = FALSE, do.retro = FALSE, MakeADFun.silent = FALSE)
#' catch <- advice_fn(mod, pro.yr = 3, type = 1)
#' }

advice_fn <- function(em, pro.yr = assess.interval, hcr.type = 1, hcr.opts = NULL) {
  cat(paste0("----------------\nHarvest Control Rule type ", hcr.type, " is in use!\n----------------\n"))
  
  proj_opts <- list(
    n.yrs = pro.yr,
    use.FXSPR = TRUE,
    avg.yrs = tail(em$years, 5),
    percentFXSPR = 75,
    percentSPR = 40
  )
  
  em_proj <- project_wham(em, proj.opts = proj_opts, MakeADFun.silent = TRUE) # Projected version of the em
  
  if (hcr.type == 1) {
    advice <- em_proj$rep$pred_catch[length(em_proj$years) + 1:pro.yr,]
  } else if (hcr.type == 2) {
    advice <- em_proj$rep$pred_catch[length(em_proj$years) + 1:pro.yr,]
    if (nrow(advice) != 1) {
      advice <- colMeans(advice) # Mean of the projected catch over the next 5 years fishing at F40
    }
    advice <- matrix(rep(advice, pro.yr), ncol = length(advice), byrow = TRUE)
  } else if (hcr.type == 3) {
    if (is.null(hcr.opts)) {
      max_percent <- 75
      min_percent <- 0.01
      BThresh_up <- 0.5
      BThresh_low <- 0.1
    } else {
      max_percent <- hcr.opts$max_percent
      min_percent <- hcr.opts$min_percent
      BThresh_up <- hcr.opts$BThresh_up
      BThresh_low <- hcr.opts$BThresh_low
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
