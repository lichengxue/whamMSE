#' Generate Catch Advice for Management Strategy Evaluation
#'
#' Generate catch advice based on harvest control rules (HCRs) using a fitted WHAM
#' estimation model. This function extends \code{\link{project_wham}} by embedding
#' HCR-specific logic (F%SPR, constant catch, hockey-stick), while also allowing
#' users to supply a full set of custom projection options via \code{proj.opts}.
#'
#' @param em A fitted WHAM estimation model object.
#' @param pro.yr Integer. Number of years for projection. Default is taken from
#'   \code{assess.interval}.
#' @param hcr A list specifying the harvest control rule:
#'   \describe{
#'     \item{\code{hcr.type}}{Integer specifying the HCR type:
#'       \enumerate{
#'         \item F\_XSPR (e.g., F\_40\%SPR) or F\_MSY if \code{use.FMSY=TRUE}.
#'         \item Constant catch (replicate projected catch; see Details).
#'         \item Hockey-stick scaling based on biomass thresholds (sets percent of F\_XSPR).
#'       }}
#'     \item{\code{hcr.opts}}{List of options controlling HCR behavior:
#'       \itemize{
#'         \item \code{use_FXSPR}, \code{percentFXSPR}
#'         \item \code{use_FMSY}, \code{percentFMSY}
#'         \item \code{avg_yrs}, \code{cont.M.re}, \code{cont.move.re}
#'         \item \code{max_percent}, \code{min_percent}, \code{BThresh_up}, \code{BThresh_low} (for HCR 3)
#'       }}
#'   }
#' @param proj.opts A named list of projection options passed directly to
#'   \code{\link{project_wham}}. Any field here will override values set by
#'   \code{hcr.opts}. Common options include:
#'   \describe{
#'     \item{\code{n.yrs}}{Integer. Number of years to project. Default = \code{pro.yr}.}
#'     \item{\code{use.last.F}}{Logical. Use terminal year F. Default = FALSE (HCR logic takes precedence).}
#'     \item{\code{use.avg.F}}{Logical. Use average F over years defined by \code{avg.yrs}.}
#'     \item{\code{use.FXSPR}}{Logical. Project using F at X\% SPR.}
#'     \item{\code{percentFXSPR}}{Numeric. Percent of F\_XSPR to use.}
#'     \item{\code{use.FMSY}}{Logical. Use F\_MSY for projections.}
#'     \item{\code{percentFMSY}}{Numeric. Percent of F\_MSY to use.}
#'     \item{\code{proj.F}}{Numeric vector length \code{n.yrs}. User F by year.}
#'     \item{\code{proj.catch}}{Numeric vector length \code{n.yrs}. User aggregate catch by year.}
#'     \item{\code{avg.yrs}}{Vector of years used to average inputs (M, maturity, WAA, selectivity, etc.). Default = last 5 years.}
#'     \item{\code{cont.ecov}}{Logical. Continue environmental covariate process. Default = TRUE.}
#'     \item{\code{use.last.ecov}}{Logical. Hold ecov at terminal year.}
#'     \item{\code{avg.ecov.yrs}}{Vector. Average ecov over these years.}
#'     \item{\code{proj.ecov}}{Matrix (#projection years × #ecovs). User covariates.}
#'     \item{\code{cont.M.re}}{Logical. Continue M REs.}
#'     \item{\code{cont.move.re}}{Logical. Continue movement REs.}
#'     \item{\code{cont.L.re}}{Logical. Continue “extra mortality” REs.}
#'     \item{\code{avg.rec.yrs}}{Vector. Years for recruitment distribution if fixed-effects recruitment.}
#'     \item{\code{proj_F_opt}}{Integer vector (length = n.yrs): 1=terminal F, 2=avg F, 3=F\_XSPR, 4=user F, 5=user catch, 6=F\_MSY.}
#'     \item{\code{proj_Fcatch}}{Vector (n.yrs) or matrix (n.yrs × nfleets) when \code{proj_F_opt} = 4 or 5.}
#'     \item{\code{proj_mature}}{3D array (nstocks × n.yrs × nages) of maturity-at-age.}
#'     \item{\code{proj_waa}}{3D array (#waa sources × n.yrs × nages) of weight-at-age.}
#'     \item{\code{proj_R_opt}}{Integer: 1=continue RE, 2=avg consistent w/ BRPs (cancel BC), 3=avg devs over \code{avg.yrs.R}, 4=no devs.}
#'     \item{\code{proj_NAA_opt}}{Integer: 1=continue RE, 2=avg devs over \code{avg.yrs.NAA}, 3=no devs.}
#'     \item{\code{proj_NAA_init}}{Numeric. Default starting value for NAA REs in projection years.}
#'     \item{\code{proj_F_init}}{Numeric. Initial guess for F when solving F given catch.}
#'     \item{\code{avg.yrs.sel}}{List length nfleets. Years to average selectivity/FAA.}
#'     \item{\code{avg.yrs.waacatch}}{List length nfleets. Years to average catch WAA.}
#'     \item{\code{avg.yrs.waassb}}{List length nstocks. Years to average SSB WAA.}
#'     \item{\code{avg.yrs.mature}}{List length nstocks. Years to average maturity.}
#'     \item{\code{avg.yrs.L}}{List length nregions. Years to average extra mortality-at-age.}
#'     \item{\code{avg.yrs.M}}{List (nstocks × nregions). Years to average natural mortality-at-age.}
#'     \item{\code{avg.yrs.move}}{List (nstocks × nregions). Years to average movement.}
#'     \item{\code{avg.yrs.R}}{List (nstocks). Years to average recruitment deviations.}
#'     \item{\code{avg.yrs.NAA}}{List (nstocks × nregions). Years to average NAA deviations.}
#'   }
#'
#' @details
#' Projection options are constructed in three steps:
#' \enumerate{
#'   \item Defaults are set internally.
#'   \item HCR logic modifies relevant fields (e.g. \code{percentFXSPR} for HCR 3).
#'   \item Any fields in \code{proj.opts} overwrite both defaults and HCR-derived settings.
#' }
#'
#' \strong{HCR 1}: Uses F\_XSPR (or F\_MSY if \code{use_FMSY=TRUE}) and returns projected catches.
#'
#' \strong{HCR 2}: “Constant catch” advisory derived by projecting, then repeating the
#' projected catch vector (if multiple projection years, the mean of those projected
#' catches is repeated across \code{pro.yr} years).
#'
#' \strong{HCR 3}: Hockey-stick: compute ratio \eqn{SSB_t / SSB40} (global) and map it
#' linearly to \code{percentFXSPR} between \code{BThresh_low} and \code{BThresh_up},
#' clamped by \code{min_percent} and \code{max_percent}; then project at that F\_XSPR%.
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
#' @export
advice_fn <- function(em, pro.yr = assess.interval, hcr = NULL, proj.opts = list()) {
  
  ## ---- Step 1. normalize HCR input
  if (is.null(hcr)) hcr <- list(hcr.type = 1, hcr.opts = list())
  hcr.type <- if (is.null(hcr$hcr.type)) 1 else hcr$hcr.type
  hcr.opts <- if (is.null(hcr$hcr.opts)) list() else hcr$hcr.opts
  cat(paste0("\nHarvest Control Rule type ", hcr.type, "\n"))
  
  ## ---- Step 2. HCR-level defaults
  use_FXSPR    <- if (is.null(hcr.opts$use_FXSPR)) TRUE  else hcr.opts$use_FXSPR
  percentFXSPR <- if (is.null(hcr.opts$percentFXSPR)) 75 else hcr.opts$percentFXSPR
  use_FMSY     <- if (is.null(hcr.opts$use_FMSY)) FALSE else hcr.opts$use_FMSY
  percentFMSY  <- if (is.null(hcr.opts$percentFMSY)) 75 else hcr.opts$percentFMSY
  avg_yrs_n    <- if (is.null(hcr.opts$avg_yrs)) 5 else hcr.opts$avg_yrs
  if (!is.null(hcr.opts$avg_yrs) && hcr.opts$avg_yrs > length(em$years)) {
    avg_yrs_n <- length(em$years)
  }
  cont.M.re_h  <- if (is.null(hcr.opts$cont.M.re)) FALSE else hcr.opts$cont.M.re
  cont.move.re_h <- if (is.null(hcr.opts$cont.move.re) || em$input$data$n_regions == 1) {
    NULL
  } else {
    hcr.opts$cont.move.re
  }
  
  ## ---- Step 3. defaults for all proj.opts
  defaults <- list(
    n.yrs          = pro.yr,
    use.last.F     = FALSE,   # let HCR set the F-mode unless user overrides
    use.avg.F      = FALSE,
    use.FXSPR      = use_FXSPR,
    percentFXSPR   = percentFXSPR,
    use.FMSY       = use_FMSY,
    percentFMSY    = percentFMSY,
    proj.F         = NULL,
    proj.catch     = NULL,
    avg.yrs        = NULL,    # filled below as tail(em$years, avg_yrs_n)
    cont.ecov      = TRUE,
    use.last.ecov  = FALSE,
    avg.ecov.yrs   = NULL,
    proj.ecov      = NULL,
    cont.M.re      = cont.M.re_h,
    cont.move.re   = cont.move.re_h,
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
    avg.yrs.sel      = NULL,
    avg.yrs.waacatch = NULL,
    avg.yrs.waassb   = NULL,
    avg.yrs.mature   = NULL,
    avg.yrs.L        = NULL,
    avg.yrs.M        = NULL,
    avg.yrs.move     = NULL,
    avg.yrs.R        = NULL,
    avg.yrs.NAA      = NULL
  )
  
  ## default averaging window for inputs
  defaults$avg.yrs <- tail(em$years, avg_yrs_n)
  
  ## ---- Step 4. Merge user overrides -> user proj.opts wins
  proj_opts <- modifyList(defaults, proj.opts)
  
  ## If user specified avg.yrs as an integer count, convert to tail years
  if (!is.null(proj.opts$avg.yrs) && is.numeric(proj.opts$avg.yrs) && length(proj.opts$avg.yrs) == 1) {
    n_tail <- min(as.integer(proj.opts$avg.yrs), length(em$years))
    proj_opts$avg.yrs <- tail(em$years, n_tail)
  }
  
  ## Drop NULLs so absent options don't appear
  proj_opts <- Filter(Negate(is.null), proj_opts)
  
  ## Helper: enforce exactly one F-spec method is active
  flags <- c(
    use.last.F = isTRUE(proj_opts$use.last.F),
    use.avg.F  = isTRUE(proj_opts$use.avg.F),
    use.FXSPR  = isTRUE(proj_opts$use.FXSPR),
    use.FMSY   = isTRUE(proj_opts$use.FMSY),
    proj.F     = !is.null(proj_opts$proj.F),
    proj.catch = !is.null(proj_opts$proj.catch)
  )
  if (sum(flags) != 1 && hcr.type %in% c(1,2)) {
    stop(sprintf("Exactly one F-spec must be set. You set: %s",
                 paste(names(flags)[flags], collapse = ", ")), call. = FALSE)
  }
  
  ## ---- HCR type 1 & 2: project first
  if (hcr.type %in% 1:2) {
    em_proj <- project_wham(em, proj.opts = proj_opts, MakeADFun.silent = TRUE)
    
    pc <- em_proj$rep$pred_catch
    nr <- nrow(pc)
    if (pro.yr > nr) stop("pro.yr exceeds available rows in pred_catch.", call. = FALSE)
    proj_rows <- (nr - pro.yr + 1):nr
    advice_mat <- pc[proj_rows, , drop = FALSE]
    
    if (hcr.type == 1) {
      advice <- advice_mat
    } else {
      ## HCR 2: constant catch — replicate the (columnwise) mean of projected catches
      cc <- if (nrow(advice_mat) == 1) as.numeric(advice_mat) else colMeans(advice_mat)
      advice <- matrix(rep(cc, each = pro.yr), nrow = pro.yr, byrow = TRUE)
      colnames(advice) <- colnames(advice_mat)
      rownames(advice) <- NULL
    }
  }
  
  ## ---- HCR type 3: hockey-stick on SSB / SSB40 -> percentFXSPR, then project
  if (hcr.type == 3) {
    max_percent <- if (is.null(hcr.opts$max_percent)) 75    else hcr.opts$max_percent
    min_percent <- if (is.null(hcr.opts$min_percent)) 0.01  else hcr.opts$min_percent
    BThresh_up  <- if (is.null(hcr.opts$BThresh_up))  0.5   else hcr.opts$BThresh_up
    BThresh_low <- if (is.null(hcr.opts$BThresh_low)) 0.1   else hcr.opts$BThresh_low
    
    if (is.null(em$rep$log_SSB_FXSPR))
      stop("HCR 3 requires em$rep$log_SSB_FXSPR.", call. = FALSE)
    
    lspr <- em$rep$log_SSB_FXSPR
    ## Typical WHAM: last column is GLOBAL when ncol == ncol(SSB)+1
    if (ncol(lspr) == ncol(em$rep$SSB) + 1) {
      SSB40_global <- exp(lspr[nrow(lspr), ncol(lspr)])
    } else {
      SSB40_global <- exp(lspr[nrow(lspr), ncol(lspr)])
    }
    SSB_t <- sum(em$rep$SSB[nrow(em$rep$SSB), ])
    ratio <- SSB_t / SSB40_global
    
    if (ratio >= BThresh_up) {
      percent <- max_percent
    } else if (ratio > BThresh_low) {
      slope   <- (max_percent - min_percent) / (BThresh_up - BThresh_low)
      percent <- slope * (ratio - BThresh_low) + min_percent
    } else {
      percent <- min_percent
    }
    cat(sprintf("SSB_t/SSB40 = %.3f -> percentFXSPR = %.2f\n", ratio, percent))
    
    ## Enforce FXSPR only
    proj_opts$use.last.F   <- FALSE
    proj_opts$use.avg.F    <- FALSE
    proj_opts$use.FMSY     <- FALSE
    proj_opts$proj.F       <- NULL
    proj_opts$proj.catch   <- NULL
    proj_opts$use.FXSPR    <- TRUE
    proj_opts$percentFXSPR <- as.numeric(percent)
    
    flags <- c(
      use.last.F = isTRUE(proj_opts$use.last.F),
      use.avg.F  = isTRUE(proj_opts$use.avg.F),
      use.FXSPR  = isTRUE(proj_opts$use.FXSPR),
      use.FMSY   = isTRUE(proj_opts$use.FMSY),
      proj.F     = !is.null(proj_opts$proj.F),
      proj.catch = !is.null(proj_opts$proj.catch)
    )
    if (sum(flags) != 1) {
      stop(sprintf("Exactly one F-spec must be set (HCR 3). You set: %s",
                   paste(names(flags)[flags], collapse = ", ")), call. = FALSE)
    }
    
    em_proj <- project_wham(em, proj.opts = proj_opts, MakeADFun.silent = TRUE)
    
    pc <- em_proj$rep$pred_catch
    nr <- nrow(pc)
    if (pro.yr > nr) stop("pro.yr exceeds available rows in pred_catch.", call. = FALSE)
    proj_rows <- (nr - pro.yr + 1):nr
    advice <- pc[proj_rows, , drop = FALSE]
  }
  
  ## ---- Print only non-null projection options
  cat("Final Projection Options:\n")
  for (proj_name in names(proj_opts)) {
    val <- proj_opts[[proj_name]]
    if (!is.null(val)) {
      cat(sprintf(" %s: %s\n", proj_name, toString(val)))
    }
  }
  
  return(advice)
}
