#' Plot terminal-year SSB and F bias by assessment
#'
#' For each assessment in the feedback loop (each element of \code{em_list}),
#' compute the terminal-year global SSB and F bias for every estimation model
#' and realization, relative to the corresponding operating model (OM) values
#' for that year. SSB is row-summed across regions for both OM and EM; F uses
#' the global F column (last column of \code{Fbar}) in both OM and EM.
#'
#' For assessment \eqn{k}, the function:
#' \itemize{
#'   \item Takes the EM's terminal-year global SSB (sum over regions) and aligns
#'         it with the OM global SSB in the same year.
#'   \item Takes the EM's terminal-year global F (last column of \code{Fbar})
#'         and aligns it with the OM global F in the same year.
#'   \item Computes relative bias:
#'         \deqn{\text{Bias} = \frac{\text{EM}}{\text{OM}} - 1}
#'         separately for SSB and F.
#'   \item Collects these across estimation models and realizations and
#'         produces a boxplot- or median-IQR-style summary, faceted by
#'         \code{"SSB"} and \code{"F"}.
#' }
#'
#' One PNG file is saved per assessment in:
#' \code{file.path(main.dir, sub.dir, "Diagnostics")}, with names like
#' \code{"Terminal_bias_assessment_01.png"}.
#'
#' @param mods A list of model outputs. For single-simulation runs
#'   (\code{is.nsim = FALSE}) this is typically
#'   \code{list(Mod1, Mod2, ...)} where each \code{ModX} has components
#'   \code{$om} and \code{$em_list}. For multi-simulation runs
#'   (\code{is.nsim = TRUE}) this is typically
#'   \code{list(sim1 = list(Mod1, Mod2, ...), sim2 = list(...), ...)}.
#' @param is.nsim Logical; \code{FALSE} if \code{mods} is a simple list of
#'   models, \code{TRUE} if \code{mods} is a list of simulations each
#'   containing a list of models.
#' @param main.dir Main output directory.
#' @param sub.dir Subdirectory (relative to \code{main.dir}) where the
#'   \code{"Diagnostics"} folder will be created and PNGs saved.
#' @param width,height Figure size passed to \code{ggsave()} (in inches).
#' @param dpi Resolution (dots per inch) for the saved PNG.
#' @param col.opt Viridis color option (passed to
#'   \code{scale_color_viridis_d(option = ...)}).
#' @param plot.style Either \code{"boxplot"} or \code{"median_iqr"} to control
#'   the style of the summary graphic.
#' @param outlier.opt Passed to \code{geom_boxplot(outlier.shape = ...)} when
#'   \code{plot.style = "boxplot"}. Use \code{NA} to hide outliers.
#' @param show.whisker Logical; if \code{TRUE} and
#'   \code{plot.style = "median_iqr"}, draw whiskers at 1.5×IQR (clipped to
#'   observed min/max).
#' @param new_model_names Optional character vector of new model names to use
#'   on the x-axis and in the legend. Length must match the number of models.
#'
#' @return A named list of \code{ggplot} objects (one per assessment) is
#'   returned (invisibly). PNG files are written to
#'   \code{file.path(main.dir, sub.dir, "Diagnostics")}.
#'
#' @details
#' The function assumes:
#' \itemize{
#'   \item OM SSB is in \code{om$rep$SSB}, with columns = regions.
#'   \item OM F is in \code{om$rep$Fbar}, with columns = fleets, regions,
#'         and a final global column; the last column is taken as global F.
#'   \item Each EM object in \code{em_list[[k]]} has:
#'         \code{SSB} (or \code{rep$SSB}) and \code{Fbar} (or
#'         \code{rep$Fbar}) with the same structure as the OM.
#'   \item The EM terminal-year index (number of rows in EM SSB) is interpreted
#'         as the year index in the OM; the function aligns OM and EM by this
#'         index for each assessment.
#' }
#'
#' @examples
#' \dontrun{
#' # Suppose `mods` is a multi-simulation object from whamMSE with feedback EMs:
#' plot_terminal_ssb_f_bias_by_assessment(
#'   mods           = mods,
#'   is.nsim        = TRUE,
#'   main.dir       = "Results",
#'   sub.dir        = "Report",
#'   width          = 10,
#'   height         = 7,
#'   dpi            = 300,
#'   col.opt        = "D",
#'   plot.style     = "median_iqr",
#'   outlier.opt    = NA,
#'   show.whisker   = TRUE,
#'   new_model_names = c("Equal", "FleetCatch", "IndexEqual")
#' )
#' }
#'
#' @export
plot_terminal_ssb_f_bias_by_assessment <- function(
    mods,
    is.nsim,
    main.dir,
    sub.dir,
    width           = 10,
    height          = 7,
    dpi             = 300,
    col.opt         = "D",
    plot.style      = "median_iqr",
    outlier.opt     = NA,
    show.whisker    = TRUE,
    new_model_names = NULL
) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(rlang)
  
  ## ------------------------------------------------------------------
  ## 1. Helper to pull EM SSB/Fbar regardless of where they are stored
  ## ------------------------------------------------------------------
  get_em_ssb <- function(em) {
    if (!is.null(em$SSB)) {
      return(as.matrix(em$SSB))
    } else if (!is.null(em$rep$SSB)) {
      return(as.matrix(em$rep$SSB))
    } else {
      stop("EM object does not contain SSB or rep$SSB.")
    }
  }
  
  get_em_fbar <- function(em) {
    if (!is.null(em$Fbar)) {
      return(as.matrix(em$Fbar))
    } else if (!is.null(em$rep$Fbar)) {
      return(as.matrix(em$rep$Fbar))
    } else {
      stop("EM object does not contain Fbar or rep$Fbar.")
    }
  }
  
  ## ------------------------------------------------------------------
  ## 2. Determine number of assessments
  ## ------------------------------------------------------------------
  if (!is.nsim) {
    if (is.null(mods[[1]]$em_list))
      stop("mods[[1]]$em_list is NULL; cannot find assessments.")
    n_assess <- length(mods[[1]]$em_list)
  } else {
    if (is.null(mods[[1]][[1]]$em_list))
      stop("mods[[1]][[1]]$em_list is NULL; cannot find assessments.")
    n_assess <- length(mods[[1]][[1]]$em_list)
  }
  
  if (n_assess == 0) {
    stop("No assessments found (em_list is empty).")
  }
  
  ## ------------------------------------------------------------------
  ## 3. Build plots for each assessment
  ## ------------------------------------------------------------------
  plot_list  <- vector("list", length = n_assess)
  names(plot_list) <- sprintf("assessment_%02d", seq_len(n_assess))
  
  # Output folder (Diagnostics)
  diag_dir <- file.path(main.dir, sub.dir, "Diagnostics")
  if (!dir.exists(diag_dir)) dir.create(diag_dir, recursive = TRUE)
  
  # Helper for ordinal label ("1st", "2nd", "3rd", ...)
  ordinal <- function(k) {
    if (k %% 100 %in% c(11, 12, 13)) return(paste0(k, "th"))
    last <- k %% 10
    if (last == 1) return(paste0(k, "st"))
    if (last == 2) return(paste0(k, "nd"))
    if (last == 3) return(paste0(k, "rd"))
    paste0(k, "th")
  }
  
  for (k in seq_len(n_assess)) {
    df_k <- NULL
    
    if (!is.nsim) {
      ## ---------------- Single-sim case ----------------
      for (m in seq_along(mods)) {
        om_obj   <- mods[[m]]$om
        em_list  <- mods[[m]]$em_list
        if (is.null(em_list) || length(em_list) < k || is.null(em_list[[k]])) next
        
        em <- em_list[[k]]
        
        # OM global SSB (rowsum over regions)
        om_ssb_mat    <- as.matrix(om_obj$rep$SSB)
        om_ssb_global <- rowSums(om_ssb_mat, na.rm = TRUE)
        
        # EM global SSB
        em_ssb_mat    <- get_em_ssb(em)
        em_ssb_global <- rowSums(em_ssb_mat, na.rm = TRUE)
        
        # Use terminal year from EM (kth assessment)
        nT <- length(em_ssb_global)
        if (nT > length(om_ssb_global)) {
          nT <- length(om_ssb_global)  # safety clip
        }
        
        ssb_em_T <- em_ssb_global[nT]
        ssb_om_T <- om_ssb_global[nT]
        ssb_bias <- ssb_em_T / ssb_om_T - 1
        
        # OM global F (last column)
        om_f_mat    <- as.matrix(om_obj$rep$Fbar)
        om_f_global <- om_f_mat[, ncol(om_f_mat)]
        
        # EM global F
        em_f_mat    <- get_em_fbar(em)
        em_f_global <- em_f_mat[, ncol(em_f_mat)]
        
        if (length(em_f_global) < nT) {
          nT_f <- length(em_f_global)
        } else {
          nT_f <- nT
        }
        
        f_em_T <- em_f_global[nT_f]
        f_om_T <- om_f_global[nT_f]
        f_bias <- f_em_T / f_om_T - 1
        
        df_k <- bind_rows(
          df_k,
          data.frame(
            Assessment  = k,
            Model       = paste0("Model", m),
            Realization = 1L,
            SSB_bias    = ssb_bias,
            F_bias      = f_bias
          )
        )
      }
      
    } else {
      ## ---------------- Multi-sim case ----------------
      for (r in seq_along(mods)) {
        for (m in seq_along(mods[[r]])) {
          om_obj   <- mods[[r]][[m]]$om
          em_list  <- mods[[r]][[m]]$em_list
          if (is.null(em_list) || length(em_list) < k || is.null(em_list[[k]])) next
          
          em <- em_list[[k]]
          
          # OM global SSB (rowsum over regions)
          om_ssb_mat    <- as.matrix(om_obj$rep$SSB)
          om_ssb_global <- rowSums(om_ssb_mat, na.rm = TRUE)
          
          # EM global SSB
          em_ssb_mat    <- get_em_ssb(em)
          em_ssb_global <- rowSums(em_ssb_mat, na.rm = TRUE)
          
          nT <- length(em_ssb_global)
          if (nT > length(om_ssb_global)) {
            nT <- length(om_ssb_global)  # safety clip
          }
          
          ssb_em_T <- em_ssb_global[nT]
          ssb_om_T <- om_ssb_global[nT]
          ssb_bias <- ssb_em_T / ssb_om_T - 1
          
          # OM global F (last column)
          om_f_mat    <- as.matrix(om_obj$rep$Fbar)
          om_f_global <- om_f_mat[, ncol(om_f_mat)]
          
          # EM global F
          em_f_mat    <- get_em_fbar(em)
          em_f_global <- em_f_mat[, ncol(em_f_mat)]
          
          if (length(em_f_global) < nT) {
            nT_f <- length(em_f_global)
          } else {
            nT_f <- nT
          }
          
          f_em_T <- em_f_global[nT_f]
          f_om_T <- om_f_global[nT_f]
          f_bias <- f_em_T / f_om_T - 1
          
          df_k <- bind_rows(
            df_k,
            data.frame(
              Assessment  = k,
              Model       = paste0("Model", m),
              Realization = r,
              SSB_bias    = ssb_bias,
              F_bias      = f_bias
            )
          )
        }
      }
    }
    
    # If nothing for this assessment, skip
    if (is.null(df_k) || nrow(df_k) == 0) {
      warning("No data for assessment ", k, "; skipping plot.")
      next
    }
    
    ## ----------------------------------------------------------------
    ## Optional relabeling of models
    ## ----------------------------------------------------------------
    if (!is.null(new_model_names)) {
      old_levels <- sort(unique(as.character(df_k$Model)))
      if (length(new_model_names) != length(old_levels)) {
        stop("Length of new_model_names must match the number of models.")
      }
      df_k$Model <- factor(
        df_k$Model,
        levels = old_levels,
        labels = new_model_names
      )
    } else {
      df_k$Model <- factor(df_k$Model)
    }
    
    ## ----------------------------------------------------------------
    ## Long format: one column 'Metric' (SSB/F) and 'Bias'
    ## ----------------------------------------------------------------
    df_long <- df_k %>%
      pivot_longer(
        cols      = c("SSB_bias", "F_bias"),
        names_to  = "Metric",
        values_to = "Bias"
      ) %>%
      mutate(
        Metric = dplyr::recode(Metric,
                               SSB_bias = "SSB",
                               F_bias   = "F")
      )
    
    ## ----------------------------------------------------------------
    ## Build plot (same styles as other performance functions)
    ## ----------------------------------------------------------------
    title_txt <- paste0(
      "Terminal-year SSB and F bias (",
      ordinal(k),
      " assessment)"
    )
    
    if (plot.style == "boxplot") {
      p <- ggplot(df_long, aes(x = Model, y = Bias, color = Model)) +
        geom_boxplot(lwd = 0.8, outlier.shape = outlier.opt) +
        facet_grid(Metric ~ ., scales = "free") +
        scale_color_viridis_d(option = col.opt) +
        ggtitle(title_txt) +
        ylab("Relative bias (EM / OM - 1)") +
        xlab("Estimation model") +
        theme_bw()
      
    } else if (plot.style == "median_iqr") {
      # Summary stats with 1.5x IQR whiskers, per Model x Metric
      res_summary <- df_long %>%
        group_by(Model, Metric) %>%
        summarise(
          q1  = quantile(Bias, 0.25, na.rm = TRUE),
          med = median(Bias, na.rm = TRUE),
          q3  = quantile(Bias, 0.75, na.rm = TRUE),
          iqr = q3 - q1,
          .groups = "drop"
        ) %>%
        mutate(
          x    = as.numeric(factor(Model)),
          ymin = if (show.whisker) q1 - 1.5 * iqr else NA_real_,
          ymax = if (show.whisker) q3 + 1.5 * iqr else NA_real_
        )
      
      # Clip whiskers to observed range
      res_limits <- df_long %>%
        group_by(Model, Metric) %>%
        summarise(
          min_val = min(Bias, na.rm = TRUE),
          max_val = max(Bias, na.rm = TRUE),
          .groups = "drop"
        )
      
      res_summary <- left_join(res_summary, res_limits,
                               by = c("Model", "Metric")) %>%
        mutate(
          ymin = pmax(ymin, min_val),
          ymax = pmin(ymax, max_val)
        )
      
      p <- ggplot(res_summary, aes(x = x, color = Model)) +
        { if (show.whisker)
          geom_segment(aes(x = x, xend = x, y = ymin, yend = q1))
        } +
        { if (show.whisker)
          geom_segment(aes(x = x, xend = x, y = q3, yend = ymax))
        } +
        geom_rect(
          aes(xmin = x - 0.3, xmax = x + 0.3,
              ymin = q1, ymax = q3, col = Model),
          fill = NA, linewidth = 0.8
        ) +
        geom_segment(
          aes(x = x - 0.3, xend = x + 0.3,
              y = med, yend = med, color = Model),
          linewidth = 0.8
        ) +
        scale_x_continuous(
          breaks = res_summary$x,
          labels = res_summary$Model
        ) +
        facet_grid(Metric ~ ., scales = "free") +
        scale_color_viridis_d(option = col.opt) +
        ggtitle(title_txt) +
        ylab("Relative bias (EM / OM - 1)") +
        xlab("Estimation model") +
        theme_bw()
      
    } else {
      stop("Unknown plot.style. Choose 'boxplot' or 'median_iqr'.")
    }
    
    ## ----------------------------------------------------------------
    ## Save figure for this assessment
    ## ----------------------------------------------------------------
    plot_name <- sprintf("Terminal_bias_assessment_%02d.png", k)
    ggsave(
      filename = file.path(diag_dir, plot_name),
      plot     = p,
      width    = width,
      height   = height,
      dpi      = dpi
    )
    
    plot_list[[k]] <- p
  }
  
  invisible(plot_list)
}
