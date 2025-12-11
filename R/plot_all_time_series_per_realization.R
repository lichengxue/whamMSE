#' Plot SSB, Fbar, and Catch time series by realization
#'
#' Create time-series line plots of operating model SSB, global Fbar, and Catch
#' for each realization in an MSE run, overlaid across estimation models.  
#' Plots are written to a single multi-page PDF, with three pages per realization:
#' \enumerate{
#'   \item SSB time series,
#'   \item global Fbar time series, and
#'   \item Catch time series.
#' }
#'
#' The function works for both single- and multi-realization objects:
#' \itemize{
#'   \item If \code{mods} is a simple list of models (single realization),
#'         it is internally wrapped as a one-element list.
#'   \item If \code{mods} is a list of realizations, each realization is
#'         assumed to be a list of models with an \code{$om} component.
#' }
#'
#' @param mods A list of MSE results, either:
#'   \itemize{
#'     \item a simple list of models for a single realization
#'           (\code{list(Mod1, Mod2, ...)}), or
#'     \item a list of realizations, where each realization is a list of models
#'           (\code{list(real1 = list(Mod1, Mod2, ...), real2 = list(...), ...)}).
#'   }
#'   Each model is expected to contain an \code{$om} component with:
#'   \code{$om$years}, \code{$om$rep$SSB}, \code{$om$rep$Fbar},
#'   and \code{$om$rep$pred_catch}, as well as
#'   \code{$om$input$data$n_fleets} and \code{$om$input$data$n_regions}.
#' @param main.dir Character path to the main output directory.
#' @param sub.dir Character folder name inside \code{main.dir} where the PDF
#'   will be written. The directory is created if it does not exist.
#' @param width Numeric, width of each PDF page in inches (default: \code{10}).
#' @param height Numeric, height of each PDF page in inches (default: \code{7}).
#' @param dpi Resolution hint (dots per inch). Included for interface
#'   consistency with other plotting functions; not directly used by the
#'   \code{\link[grDevices]{pdf}} device.
#' @param col.opt Character string passed to \code{scale_color_viridis_d(option = ...)}
#'   to control the viridis palette option (default: \code{"D"}).
#' @param new_model_names Optional character vector of custom model names.  
#'   If provided, must have length equal to the number of models in each
#'   realization. These labels replace \code{"Model1"}, \code{"Model2"}, etc.
#' @param pdf_name Name of the output PDF file (default:
#'   \code{"All_TS_by_realization.pdf"}).
#' @param realizations_to_plot Optional integer vector of realization indices
#'   to include. If \code{NULL}, all realizations are plotted
#'   (\code{seq_len(n_realizations)}).
#'
#' @details
#' The function first detects whether \code{mods} represents a single or
#' multiple realizations using:
#' \preformatted{
#' is.nsim <- is.list(mods[[1]]) && !is.null(mods[[1]][[1]]$om)
#' }
#' If \code{is.nsim} is \code{FALSE}, a single realization is assumed and
#' \code{mods} is wrapped as \code{list(mods)} for consistent indexing.
#'
#' For each realization \code{r} in \code{realizations_to_plot}:
#' \itemize{
#'   \item SSB: \code{mods[[r]][[m]]$om$rep$SSB} is extracted for each model
#'         \code{m}, then plotted as lines over \code{Year}.
#'   \item Fbar: the global Fbar column is assumed to be the last column of
#'         \code{om$rep$Fbar}, i.e.
#'         \code{n_fleets + n_regions + 1}, where
#'         \code{n_fleets <- om$input$data$n_fleets[1]} and
#'         \code{n_regions <- om$input$data$n_regions[1]}.
#'   \item Catch: \code{mods[[r]][[m]]$om$rep$pred_catch} (typically multiple
#'         columns for fleets/regions) is pivoted to long format and plotted
#'         by \code{Label} (column name) and \code{Model}.
#' }
#'
#' Each of the three plots (SSB, Fbar, Catch) is printed as a separate page in
#' the PDF for that realization.
#'
#' @return
#' Invisibly returns \code{NULL}. The main side effect is the creation of a
#' multi-page PDF file at:
#' \code{file.path(main.dir, sub.dir, pdf_name)}.
#'
#' @seealso
#' \code{\link{plot_NAA_tile_by_realization}} for NAA deviation tiles by
#' realization, and other MSE summary plotting functions in this package.
#'
#' @examples
#' \dontrun{
#' # Suppose `mods` is a list of realizations produced by whamMSE:
#'
#' # 1) Plot all realizations with default model names:
#' plot_all_time_series_per_realization(
#'   mods      = mods,
#'   main.dir  = "Results",
#'   sub.dir   = "Diagnostics"
#' )
#'
#' # 2) Plot only realizations 1 and 3 with custom model names:
#' plot_all_time_series_per_realization(
#'   mods                 = mods,
#'   main.dir             = "Results",
#'   sub.dir              = "Diagnostics",
#'   new_model_names      = c("Equal", "FleetCatch", "IndexEqual"),
#'   realizations_to_plot = c(1, 3)
#' )
#' }
plot_all_time_series_per_realization <- function(mods, main.dir, sub.dir,
                                                 width = 10, height = 7, dpi = 300,
                                                 col.opt = "D", new_model_names = NULL,
                                                 pdf_name = "All_TS_by_realization.pdf",
                                                 realizations_to_plot = NULL) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(grDevices)
  
  # --- Detect single vs multiple realizations ---
  is.nsim <- is.list(mods[[1]]) && !is.null(mods[[1]][[1]]$om)
  
  if (is.nsim) {
    n_realizations <- length(mods)
    Years <- mods[[1]][[1]]$om$years
    n_models <- length(mods[[1]])
  } else {
    n_realizations <- 1
    Years <- mods[[1]]$om$years
    n_models <- length(mods)
    mods <- list(mods)  # Wrap single realization as list for consistency
  }
  
  # Default to all realizations
  if (is.null(realizations_to_plot)) {
    realizations_to_plot <- seq_len(n_realizations)
  }
  
  # Safety check
  if (any(realizations_to_plot > n_realizations | realizations_to_plot < 1)) {
    stop("One or more values in realizations_to_plot are out of bounds.")
  }
  
  # Create output directory if needed
  dir.create(file.path(main.dir, sub.dir), recursive = TRUE, showWarnings = FALSE)
  
  # Open multi-page PDF
  pdf(file = file.path(main.dir, sub.dir, pdf_name), width = width, height = height)
  
  for (r in realizations_to_plot) {
    # === SSB ===
    df_ssb <- lapply(seq_len(n_models), function(m) {
      data.frame(
        SSB   = mods[[r]][[m]]$om$rep$SSB,
        Model = paste0("Model", m),
        Year  = Years
      )
    }) %>% bind_rows()
    
    if (!is.null(new_model_names)) {
      df_ssb$Model <- factor(df_ssb$Model,
                             levels = paste0("Model", seq_len(n_models)),
                             labels = new_model_names)
    }
    
    df_ssb <- pivot_longer(df_ssb, cols = starts_with("SSB"),
                           names_to = "Label", values_to = "SSB")
    
    p_ssb <- ggplot(df_ssb, aes(x = Year, y = SSB, color = Model, group = Model)) +
      geom_line(size = 0.8) +
      scale_color_viridis_d(option = col.opt) +
      facet_grid(Label ~ ., scales = "free") +
      theme_bw() +
      ggtitle(paste0("SSB from realization ", r)) +
      ylab("SSB")
    
    print(p_ssb)
    
    # === Fbar (Global) ===
    n_fleets  <- mods[[r]][[1]]$om$input$data$n_fleets[1]
    n_regions <- mods[[r]][[1]]$om$input$data$n_regions[1]
    index_range <- n_fleets + n_regions + 1  # Global Fbar
    
    df_fbar <- lapply(seq_len(n_models), function(m) {
      tmp <- mods[[r]][[m]]$om$rep$Fbar[, index_range, drop = FALSE]
      tmp <- as.data.frame(tmp)
      names(tmp) <- "Global"
      tmp$Model <- paste0("Model", m)
      tmp$Year  <- Years
      tmp
    }) %>% bind_rows()
    
    if (!is.null(new_model_names)) {
      df_fbar$Model <- factor(df_fbar$Model,
                              levels = paste0("Model", seq_len(n_models)),
                              labels = new_model_names)
    }
    
    df_fbar <- pivot_longer(df_fbar, cols = starts_with("Global"),
                            names_to = "Label", values_to = "Fbar")
    
    p_fbar <- ggplot(df_fbar, aes(x = Year, y = Fbar, color = Model, group = Model)) +
      geom_line(size = 0.8) +
      scale_color_viridis_d(option = col.opt) +
      facet_grid(Label ~ ., scales = "free") +
      theme_bw() +
      ggtitle(paste0("Fbar from realization ", r)) +
      ylab("Fbar")
    
    print(p_fbar)
    
    # === Catch ===
    df_catch <- lapply(seq_len(n_models), function(m) {
      data.frame(
        Catch = mods[[r]][[m]]$om$rep$pred_catch,
        Model = paste0("Model", m),
        Year  = Years
      )
    }) %>% bind_rows()
    
    if (!is.null(new_model_names)) {
      df_catch$Model <- factor(df_catch$Model,
                               levels = paste0("Model", seq_len(n_models)),
                               labels = new_model_names)
    }
    
    df_catch <- pivot_longer(df_catch, cols = starts_with("Catch"),
                             names_to = "Label", values_to = "Catch")
    
    p_catch <- ggplot(df_catch, aes(x = Year, y = Catch, color = Model, group = Model)) +
      geom_line(size = 0.8) +
      scale_color_viridis_d(option = col.opt) +
      facet_grid(Label ~ ., scales = "free") +
      theme_bw() +
      ggtitle(paste0("Catch from realization ", r)) +
      ylab("Catch")
    
    print(p_catch)
  }
  
  dev.off()
}
