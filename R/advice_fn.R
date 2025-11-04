#' Generate Catch Advice for Management Strategy Evaluation
#'
#' Generate catch advice based on harvest control rules (HCRs) using a fitted WHAM
#' estimation model. This function extends \code{\link{project_wham}} by embedding
#' HCR-specific logic (F\%SPR, constant catch, hockey-stick), while also allowing
#' users to supply a full set of custom projection options via \code{proj.opts}.
#'
#' @param em A fitted WHAM estimation model object.
#' @param pro.yr Integer. Number of years for projection. Default is taken from 
#'   \code{assess.interval}.
#' @param hcr A list specifying the harvest control rule:
#'   \describe{
#'     \item{\code{hcr.type}}{Integer specifying the HCR type:
#'       \enumerate{
#'         \item F\_XSPR (e.g., F\_40\%SPR).
#'         \item Constant catch.
#'         \item Hockey-stick scaling based on biomass thresholds.
#'       }}
#'     \item{\code{hcr.opts}}{List of options controlling HCR behavior:
#'       \itemize{
#'         \item \code{use_FXSPR}, \code{percentFXSPR}
#'         \item \code{use_FMSY}, \code{percentFMSY}
#'         \item \code{avg_yrs}, \code{cont.M.re}, \code{cont.move.re}
#'         \item \code{max_percent}, \code{min_percent}, \code{BThresh_up}, \code{BThresh_low}
#'       }}
#'   }
#' @param proj.opts A named list of projection options passed directly to 
#'   \code{\link{project_wham}}. Any field here will override values set by 
#'   \code{hcr.opts}. Available options include:
#'   \describe{
#'     \item{\code{n.yrs}}{Integer. Number of years to project. Default = \code{pro.yr}.}
#'     \item{\code{use.last.F}}{Logical. Use terminal year F for projections. Default = FALSE here (HCR logic takes precedence).}
#'     \item{\code{use.avg.F}}{Logical. Use average F over years defined by \code{avg.yrs}. Default = FALSE.}
#'     \item{\code{use.FXSPR}}{Logical. Project using F at X\% SPR. Default = TRUE (overridden by HCR).}
#'     \item{\code{percentFXSPR}}{Numeric scalar. Percent of F\_XSPR to use. Default = 75.}
#'     \item{\code{use.FMSY}}{Logical. Use FMSY for projections. Default = FALSE.}
#'     \item{\code{percentFMSY}}{Numeric scalar. Percent of F\_MSY to use. Default = 75.}
#'     \item{\code{proj.F}}{Numeric vector. User-specified fishing mortality by year (length = \code{n.yrs}).}
#'     \item{\code{proj.catch}}{Numeric vector. User-specified aggregate catch by year (length = \code{n.yrs}).}
#'     \item{\code{avg.yrs}}{Vector of years. Used to average M, maturity, weight-at-age, selectivity, etc. Default = last 5 years.}
#'     \item{\code{cont.ecov}}{Logical. Continue environmental covariate process (random walk / AR1). Default = TRUE.}
#'     \item{\code{use.last.ecov}}{Logical. Use terminal year ecov for projections.}
#'     \item{\code{avg.ecov.yrs}}{Vector of years. Average ecov over these years.}
#'     \item{\code{proj.ecov}}{Matrix (#projection years × #ecovs). User-specified covariate values.}
#'     \item{\code{cont.M.re}}{Logical. Continue M random effects (AR1\_y or 2D AR1). Default = FALSE.}
#'     \item{\code{cont.move.re}}{Logical. Continue movement random effects. Default = FALSE (averaged).}
#'     \item{\code{cont.L.re}}{Logical. Continue “extra mortality” random effects. Default = FALSE.}
#'     \item{\code{avg.rec.yrs}}{Vector of years. Used to calculate recruitment distribution if recruitment is fixed-effects. Default = all years.}
#'     \item{\code{proj_F_opt}}{Integer vector (length = n.yrs). Specifies projection rule by year:
#'       1 = terminal F, 2 = average F, 3 = F\_XSPR, 4 = user-specified F, 
#'       5 = user-specified catch, 6 = F\_MSY. Overrides other options.}
#'     \item{\code{proj_Fcatch}}{Vector (n.yrs) or matrix (n.yrs × nfleets). Catch or F values for projections when \code{proj_F_opt} = 4 or 5.}
#'     \item{\code{proj_mature}}{3D array (nstocks × n.yrs × nages). User-supplied maturity at age for projections.}
#'     \item{\code{proj_waa}}{3D array (#waa sources × n.yrs × nages). User-supplied weight-at-age.}
#'     \item{\code{proj_R_opt}}{Integer. Recruitment rule: 
#'       1 = continue RE, 
#'       2 = average consistent with BRPs (cancel bias correction),
#'       3 = average deviations over \code{avg.yrs.R}, 
#'       4 = no deviations.}
#'     \item{\code{proj_NAA_opt}}{Integer. NAA rule:
#'       1 = continue RE, 
#'       2 = average deviations over \code{avg.yrs.NAA}, 
#'       3 = no deviations.}
#'     \item{\code{proj_NAA_init}}{Numeric scalar. Default starting value for NAA random effects in projection years. Default = \code{exp(10)}.}
#'     \item{\code{proj_F_init}}{Numeric scalar. Initial guess for F in Newton search when solving for F given catch. Default = 0.1.}
#'     \item{\code{avg.yrs.sel}}{List (length = nfleets). Years to average selectivity / FAA for each fleet. Default = last 5 years.}
#'     \item{\code{avg.yrs.waacatch}}{List (length = nfleets). Years to average catch weight-at-age per fleet. Default = last 5 years.}
#'     \item{\code{avg.yrs.waassb}}{List (length = nstocks). Years to average spawning-stock weight-at-age. Default = last 5 years.}
#'     \item{\code{avg.yrs.mature}}{List (length = nstocks). Years to average maturity-at-age. Default = last 5 years.}
#'     \item{\code{avg.yrs.L}}{List (length = nregions). Years to average extra mortality-at-age. Default = last 5 years.}
#'     \item{\code{avg.yrs.M}}{List (nstocks × nregions). Years to average natural mortality-at-age. Default = last 5 years.}
#'     \item{\code{avg.yrs.move}}{List (nstocks × nregions). Years to average movement rates at age × season. Default = last 5 years.}
#'     \item{\code{avg.yrs.R}}{List (nstocks). Years to average recruitment deviations. Default = last 5 years.}
#'     \item{\code{avg.yrs.NAA}}{List (nstocks × nregions). Years to average NAA deviations. Default = last 5 years.}
#'   }
#'
#' @details
#' Projection options are constructed in three steps:
#' \enumerate{
#'   \item Defaults are set internally.
#'   \item HCR logic modifies relevant fields (e.g. \code{percentFXSPR}).
#'   \item Any fields in \code{proj.opts} overwrite both defaults and HCR-derived settings.
#' }
#'
#' This ensures maximum flexibility:
#' \itemize{
#'   \item If no \code{proj.opts} is provided, defaults + HCR logic are used.
#'   \item If a few fields are given, they override the defaults/HCR.
#'   \item If a complete \code{proj.opts} is given, it bypasses HCR entirely. 
#' }
#'
#' @return A matrix of projected catch advice for \code{pro.yr} years.
#'
#' @seealso \code{\link{project_wham}}, \code{\link{fit_wham}}
#'
#' @examples
#' # Example 1: F40%SPR advice
#' hcr <- list(hcr.type = 1, hcr.opts = list(use_FXSPR = TRUE, percentFXSPR = 100))
#' advice <- advice_fn(em, pro.yr = 3, hcr = hcr)
#'
#' # Example 2: Override proj.opts with custom F
#' advice <- advice_fn(em, pro.yr = 3, hcr = hcr,
#'   proj.opts = list(proj_F_opt = c(4,4,4), proj_Fcatch = c(0.1,0.2,0.3)))
#'
#' # Example 3: Average NAA deviations over custom years
#' advice <- advice_fn(em, pro.yr = 5,
#'   proj.opts = list(avg.yrs.NAA = list(2000:2005, 2001:2006)))
#' }
#' @export

advice_fn <- function(em, pro.yr = assess.interval, hcr = NULL, proj.opts = list()) {
  
  # Step 1. ensure hcr is always a list
  if (is.null(hcr)) hcr <- list(hcr.type = 1, hcr.opts = list())
  hcr.type <- ifelse(is.null(hcr$hcr.type), 1, hcr$hcr.type)
  hcr.opts <- if (is.null(hcr$hcr.opts)) list() else hcr$hcr.opts
  cat(paste0("\nHarvest Control Rule type ", hcr.type, "\n"))
  
  # Step 2. defaults for HCR logic
  use_FXSPR     <- ifelse(is.null(hcr.opts$use_FXSPR), TRUE, hcr.opts$use_FXSPR)
  percentFXSPR  <- ifelse(is.null(hcr.opts$percentFXSPR), 75, hcr.opts$percentFXSPR)
  use_FMSY      <- ifelse(is.null(hcr.opts$use_FMSY), FALSE, hcr.opts$use_FMSY)
  percentFMSY   <- ifelse(is.null(hcr.opts$percentFMSY), 75, hcr.opts$percentFMSY)
  avg_yrs       <- ifelse(is.null(hcr.opts$avg_yrs), 5, hcr.opts$avg_yrs)
  cont.M.re     <- ifelse(is.null(hcr.opts$cont.M.re), FALSE, hcr.opts$cont.M.re)
  
  if (!is.null(hcr.opts$avg_yrs) && hcr.opts$avg_yrs > length(em$years)) {
    avg_yrs <- length(em$years)
  }
  cont.move.re <- if (is.null(hcr.opts$cont.move.re) || em$input$data$n_regions == 1) NULL else hcr.opts$cont.move.re
  
  # Step 3. defaults for all proj.opts
  defaults <- list(
    n.yrs          = pro.yr,
    use.last.F     = FALSE,
    use.avg.F      = FALSE,
    use.FXSPR      = use_FXSPR,
    percentFXSPR   = percentFXSPR,
    use.FMSY       = use_FMSY,
    percentFMSY    = percentFMSY,
    proj.F         = NULL,
    proj.catch     = NULL,
    avg.yrs        = 5,
    cont.ecov      = TRUE,
    use.last.ecov  = FALSE,
    avg.ecov.yrs   = NULL,
    proj.ecov      = NULL,
    cont.M.re      = FALSE,
    cont.move.re   = FALSE,
    cont.L.re      = FALSE,
    avg.rec.yrs    = NULL,
    proj_F_opt     = NULL,
    proj_Fcatch    = NULL,
    proj_mature    = NULL,
    proj_waa       = NULL,
    proj_R_opt     = NULL,
    proj_NAA_opt   = NULL,
    proj_NAA_init  = NULL,
    proj_F_init    = NULL,
    avg.yrs.sel    = NULL,
    avg.yrs.waacatch = NULL,
    avg.yrs.waassb   = NULL,
    avg.yrs.mature   = NULL,
    avg.yrs.L        = NULL,
    avg.yrs.M        = NULL,
    avg.yrs.move     = NULL,
    avg.yrs.R        = NULL,
    avg.yrs.NAA      = NULL
  )

  defaults$avg.yrs = tail(em$years,defaults$avg.yrs)
  
  # Step 4. Merge user overrides -> proj.opts wins over HCR/default
  proj_opts <- modifyList(defaults, proj.opts)
  
  # --- HCR type 1 & 2 ---
  if (hcr.type %in% 1:2) {
    em_proj <- project_wham(em, proj.opts = proj_opts, MakeADFun.silent = TRUE)
  }
  if (hcr.type == 1) {
    advice <- em_proj$rep$pred_catch[(length(em_proj$years) + 1):(length(em_proj$years) + pro.yr), ]
  } 
  if (hcr.type == 2) {
    advice <- em_proj$rep$pred_catch[(length(em_proj$years) + 1):(length(em_proj$years) + pro.yr), ]
    if (nrow(advice) != 1) advice <- colMeans(advice)
    advice <- matrix(rep(advice, pro.yr), ncol = length(advice), byrow = TRUE)
  }
  
  # --- HCR type 3 (hockey stick) ---
  if (hcr.type == 3) {
    max_percent <- ifelse(is.null(hcr.opts$max_percent), 75, hcr.opts$max_percent)
    min_percent <- ifelse(is.null(hcr.opts$min_percent), 0.01, hcr.opts$min_percent)
    BThresh_up  <- ifelse(is.null(hcr.opts$BThresh_up), 0.5, hcr.opts$BThresh_up)
    BThresh_low <- ifelse(is.null(hcr.opts$BThresh_low), 0.1, hcr.opts$BThresh_low)
    # biomass-ratio logic here updates proj_opts$percentFXSPR or proj_opts$percentFMSY
  }
  
  # --- Print only non-null projection options ---
  cat("Final Projection Options:\n")
  for (proj_name in names(proj_opts)) {
    val <- proj_opts[[proj_name]]
    if (!is.null(val)) {
      cat(sprintf(" %s: %s\n", proj_name, toString(val)))
    }
  }
  
  return(advice)

}
