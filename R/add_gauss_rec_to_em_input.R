#' Add Gaussian environmental recruitment settings to EM input
#'
#' Adds the data, parameter initial values, and parameter mapping needed to fit
#' a Gaussian environmental effect on recruitment in the assessment model (EM).
#' This function is intended for WHAM/TMB model versions that include the
#' corresponding Gaussian recruitment objects in the compiled model.
#'
#' The Gaussian recruitment effect uses an environmental covariate column
#' specified by \code{gauss_rec_em$Ecov_rec_T_col}, along with the parameters
#' \code{Topt_rec}, \code{width_rec}, and \code{beta_T_rec}.
#'
#' When \code{gauss_rec_em$use = FALSE}, this function inserts harmless
#' placeholder values and maps the Gaussian recruitment parameters out of the
#' model so that they remain fixed.
#'
#' \strong{Important:} this function should only be used with WHAM/TMB branches
#' that include the Gaussian recruitment parameters and data objects in the
#' compiled model. For branches without Gaussian recruitment functionality,
#' users should skip calling this function entirely.
#'
#' @param em_input A WHAM EM input list, typically produced by
#'   \code{\link{make_em_input}}. Must contain at least
#'   \code{em_input$data$n_stocks} and \code{em_input$data$n_Ecov}, and should
#'   contain \code{data}, \code{par}, and optionally \code{map} components.
#'
#' @param gauss_rec_em A list specifying whether and how to include a Gaussian
#'   environmental effect on recruitment. Expected elements include:
#'   \describe{
#'     \item{\code{use}}{Logical. Whether to activate the Gaussian recruitment
#'     effect. Default is \code{FALSE} if not provided.}
#'     \item{\code{Ecov_rec_T_col}}{Integer. The environmental covariate column
#'     to use for recruitment, supplied as a 1-based R index.}
#'     \item{\code{Topt_rec}}{Numeric. The optimal environmental value at which
#'     recruitment is maximized.}
#'     \item{\code{width_rec}}{Numeric. Width parameter of the Gaussian effect
#'     on the original environmental scale. Must be positive.}
#'     \item{\code{beta_T_rec}}{Numeric scalar or vector. Magnitude of the
#'     Gaussian environmental effect on recruitment. If a scalar is supplied, it
#'     is recycled across stocks.}
#'     \item{\code{estimate}}{Logical. If \code{TRUE}, the Gaussian recruitment
#'     parameters are estimated. If \code{FALSE}, they are fixed at the supplied
#'     values using \code{map}.}
#'   }
#'
#' @return The modified \code{em_input} list with Gaussian recruitment data,
#'   parameter initial values, and map entries added.
#'
#' @details
#' This function performs the following steps:
#' \enumerate{
#'   \item If \code{gauss_rec_em$use = FALSE}, it adds placeholder values for
#'   the Gaussian recruitment objects and maps the parameters out.
#'   \item If \code{gauss_rec_em$use = TRUE}, it validates the supplied Gaussian
#'   recruitment settings.
#'   \item Converts the environmental covariate column from 1-based R indexing
#'   to 0-based TMB indexing.
#'   \item Adds the required data inputs:
#'   \code{use_gauss_T_rec} and \code{Ecov_rec_T_col}.
#'   \item Adds the required parameter initial values:
#'   \code{Topt_rec}, \code{log_width_rec}, and \code{beta_T_rec}.
#'   \item Updates \code{em_input$map} so the Gaussian recruitment parameters
#'   are either estimated or fixed.
#' }
#'
#' @seealso \code{\link{make_em_input}}
#'
#' @export
add_gauss_rec_to_em_input <- function(em_input, gauss_rec_em) {

  if (is.null(gauss_rec_em$use)) gauss_rec_em$use <- FALSE

  n_stocks <- em_input$data$n_stocks
  n_Ecov   <- em_input$data$n_Ecov

  # Always pass required data objects because TMB expects them
  if (!isTRUE(gauss_rec_em$use)) {
    em_input$data$use_gauss_T_rec <- 0L
    em_input$data$Ecov_rec_T_col  <- 0L

    # harmless placeholders
    em_input$par$Topt_rec      <- 0
    em_input$par$log_width_rec <- 0
    em_input$par$beta_T_rec    <- rep(0, n_stocks)

    # optional: map out so they are fixed
    if (is.null(em_input$map)) em_input$map <- list()
    em_input$map$Topt_rec      <- factor(NA)
    em_input$map$log_width_rec <- factor(NA)
    em_input$map$beta_T_rec    <- factor(rep(NA, n_stocks))

    return(em_input)
  }

  # ---------- checks ----------
  if (is.null(gauss_rec_em$Ecov_rec_T_col)) {
    stop("gauss_rec_em$Ecov_rec_T_col must be provided when gauss_rec_em$use = TRUE.")
  }

  # allow user to pass either 1-based R index or 0-based TMB index
  # here I recommend user supplies 1-based R index
  Ecov_col_R <- gauss_rec_em$Ecov_rec_T_col

  if (Ecov_col_R < 1 || Ecov_col_R > n_Ecov) {
    stop("gauss_rec_em$Ecov_rec_T_col is out of range for em_input$data$n_Ecov.")
  }

  Ecov_col_TMB <- as.integer(Ecov_col_R - 1L)

  if (is.null(gauss_rec_em$Topt_rec)) {
    stop("gauss_rec_em$Topt_rec must be provided when gauss_rec_em$use = TRUE.")
  }

  if (is.null(gauss_rec_em$width_rec)) {
    stop("gauss_rec_em$width_rec must be provided when gauss_rec_em$use = TRUE.")
  }

  if (gauss_rec_em$width_rec <= 0) {
    stop("gauss_rec_em$width_rec must be > 0.")
  }

  if (is.null(gauss_rec_em$beta_T_rec)) {
    gauss_rec_em$beta_T_rec <- rep(0, n_stocks)
  }

  if (length(gauss_rec_em$beta_T_rec) == 1) {
    gauss_rec_em$beta_T_rec <- rep(gauss_rec_em$beta_T_rec, n_stocks)
  }

  if (length(gauss_rec_em$beta_T_rec) != n_stocks) {
    stop("gauss_rec_em$beta_T_rec must have length 1 or n_stocks.")
  }

  # ---------- add data ----------
  em_input$data$use_gauss_T_rec <- 1L
  em_input$data$Ecov_rec_T_col  <- Ecov_col_TMB

  # ---------- add parameter initials ----------
  em_input$par$Topt_rec      <- gauss_rec_em$Topt_rec
  em_input$par$log_width_rec <- log(gauss_rec_em$width_rec)
  em_input$par$beta_T_rec    <- as.numeric(gauss_rec_em$beta_T_rec)

  # ---------- map logic ----------
  if (is.null(em_input$map)) em_input$map <- list()

  if (isTRUE(gauss_rec_em$estimate)) {
    # leave them free unless user already provided maps
    if (is.null(em_input$map$Topt_rec)) {
      em_input$map$Topt_rec <- factor(1)
    }
    if (is.null(em_input$map$log_width_rec)) {
      em_input$map$log_width_rec <- factor(1)
    }
    if (is.null(em_input$map$beta_T_rec)) {
      em_input$map$beta_T_rec <- factor(seq_len(n_stocks))
    }
  } else {
    # fix them at supplied values
    em_input$map$Topt_rec      <- factor(NA)
    em_input$map$log_width_rec <- factor(NA)
    em_input$map$beta_T_rec    <- factor(rep(NA, n_stocks))
  }

  return(em_input)
}
