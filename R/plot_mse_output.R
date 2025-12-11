#' Generate MSE Output Plots in PNG, HTML, or PDF Format
#'
#' This function creates a comprehensive set of plots to visualize MSE (Management Strategy Evaluation)
#' results from a list of model simulations. Outputs can be saved as individual PNG files or as a single
#' interactive HTML or PDF report.
#'
#' @param mods A list of model outputs from multiple simulations and estimation models.
#' @param main_dir The main directory path where results should be saved.
#' @param output_dir Folder name inside main_dir where plots/reports will be saved (default: "Report").
#' @param output_format Output type: "png", "html", or "pdf".
#' @param width Width of each plot in inches.
#' @param height Height of each plot in inches.
#' @param dpi Resolution of saved plots in dots per inch.
#' @param col.opt Color palette option for `viridis`-based scales.
#' @param method method = NULL: global median, method = "mean": median of means, method = "median": median of medians
#' @param outlier.opt outlier.opt = NA (default): no outliers will be shown, outlier.opt = 0 - 25 (26 different shapes to use)
#' @param plot.style plot.style = "median_iqr" (Default): only median, IQR will be displayed in the boxplot (no outliers will be displayed), method = "boxplot": boxplot including median, IQR, and outliers
#' @param show.whisker Logical, Whether to show 1.5IQR as a whisker
#' @param f.ymin Numeric value specifying the lower limit of the y-axis for F-based performance plots (Default = NULL)
#' @param f.ymax Numeric value specifying the upper limit of the y-axis for F-based performance plots. (Default = NULL)
#' @param start.years The first year used in short-term performance summaries.
#' @param use.n.years.first Number of years to summarize at the beginning of the time series.
#' @param use.n.years.last Number of years to summarize at the end of the time series.
#' @param new_model_names Optional character vector of custom model names.
#' @param base.model Model name used as a baseline for relative difference plots.

#'
#' @return Saves output plots/reports to disk; returns nothing explicitly.
#' @export
#' 

plot_mse_output <- function(mods,
                            main_dir,
                            output_dir = "Report",
                            output_format = c("png", "html", "pdf"),
                            width = 10, 
                            height = 7, 
                            dpi = 300,
                            col.opt = "D",
                            method = "mean",
                            outlier.opt = NA,
                            plot.style = "median_iqr",
                            show.whisker = TRUE,
                            f.ymin = NULL, 
                            f.ymax = NULL, 
                            start.years = 1,
                            use.n.years.first = 5,
                            use.n.years.last = 5,
                            new_model_names = NULL,
                            base.model = "Model1") {
  
  # Match output format
  output_format <- match.arg(output_format)
  
  # Build full path
  full_output_dir <- file.path(main_dir, output_dir)
  if (!dir.exists(full_output_dir)) dir.create(full_output_dir, recursive = TRUE)
  
  # Check if it's multi-simulation
  is.nsim <- if (!is.list(mods[[1]][[1]][[1]])) FALSE else TRUE
  
  if (is.null(f.ymin)) f.ymin = NA
  if (is.null(f.ymax)) f.ymax = NA
  
  if (output_format == "png") {
    cat("\nGenerating PNG plots...\n")
    
    # List of plot functions
    plot_calls <- list(
      function() plot_ssb_time_series(mods, is.nsim, main_dir, output_dir, "SSB", width, height, dpi, col.opt, new_model_names),
      function() plot_fbar_time_series(mods, is.nsim, main_dir, output_dir, "Fbar", width, height, dpi, col.opt, new_model_names),
      function() plot_catch_time_series(mods, is.nsim, main_dir, output_dir, "Catch", width, height, dpi, col.opt, new_model_names),
      function() plot_relative_trajectories(mods, is.nsim, main_dir, output_dir, base.model, new_model_names, width, height, dpi, col.opt),
      function() plot_ssb_performance(mods, is.nsim, main_dir, output_dir, "SSB", width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, use.n.years.last, new_model_names),
      function() plot_ssb_performance(mods, is.nsim, main_dir, output_dir, "SSB", width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, use.n.years.last, new_model_names, base.model),
      function() plot_fbar_performance(mods, is.nsim, main_dir, output_dir, "Fbar", width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, f.ymin, f.ymax, use.n.years.last, new_model_names),
      function() plot_fbar_performance(mods, is.nsim, main_dir, output_dir, "Fbar", width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, f.ymin, f.ymax, use.n.years.last, new_model_names, base.model),
      function() plot_catch_performance(mods, is.nsim, main_dir, output_dir, "Catch", width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, use.n.years.last, new_model_names),
      function() plot_catch_performance(mods, is.nsim, main_dir, output_dir, "Catch", width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, use.n.years.last, new_model_names, base.model),
      function() plot_total_catch_performance(mods, is.nsim, main_dir, output_dir, "Catch", width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, use.n.years.last, new_model_names),
      function() plot_total_catch_performance(mods, is.nsim, main_dir, output_dir, "Catch", width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, use.n.years.last, new_model_names, base.model),
      function() plot_ssb_performance2(mods, is.nsim, main_dir, output_dir, "SSB", width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, use.n.years.first, start.years, new_model_names),
      function() plot_ssb_performance2(mods, is.nsim, main_dir, output_dir, "SSB", width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, use.n.years.first, start.years, new_model_names, base.model),
      function() plot_fbar_performance2(mods, is.nsim, main_dir, output_dir, "Fbar", width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, f.ymin, f.ymax, use.n.years.first, start.years, new_model_names),
      function() plot_fbar_performance2(mods, is.nsim, main_dir, output_dir, "Fbar", width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, f.ymin, f.ymax, use.n.years.first, start.years, new_model_names, base.model),
      function() plot_catch_performance2(mods, is.nsim, main_dir, output_dir, "Catch", width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, use.n.years.first, start.years, new_model_names),
      function() plot_catch_performance2(mods, is.nsim, main_dir, output_dir, "Catch", width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, use.n.years.first, start.years, new_model_names, base.model),
      function() plot_total_catch_performance2(mods, is.nsim, main_dir, output_dir, "Catch", width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, use.n.years.first, start.years, new_model_names),
      function() plot_total_catch_performance2(mods, is.nsim, main_dir, output_dir, "Catch", width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, use.n.years.first, start.years, new_model_names, base.model),
      function() plot_catch_variation(mods, is.nsim, main_dir, output_dir, "Catch", width, height, dpi, col.opt, outlier.opt, new_model_names),
      function() plot_catch_variation(mods, is.nsim, main_dir, output_dir, "Catch", width, height, dpi, col.opt, outlier.opt, new_model_names, base.model),
      function() plot_ssb_variation(mods, is.nsim, main_dir, output_dir, "SSB", width, height, dpi, col.opt, outlier.opt, new_model_names),
      function() plot_ssb_variation(mods, is.nsim, main_dir, output_dir, "SSB", width, height, dpi, col.opt, outlier.opt, new_model_names, base.model),
      function() plot_fbar_variation(mods, is.nsim, main_dir, output_dir, "Fbar", width, height, dpi, col.opt, outlier.opt, new_model_names),
      function() plot_fbar_variation(mods, is.nsim, main_dir, output_dir, "Fbar", width, height, dpi, col.opt, outlier.opt, new_model_names, base.model),
      function() plot_ssb_status(mods, is.nsim, main_dir, output_dir, "SSB_status", width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, use.n.years.last, new_model_names),
      function() plot_ssb_status(mods, is.nsim, main_dir, output_dir, "SSB_status", width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, use.n.years.last, new_model_names, base.model),
      function() plot_ssb_status2(mods, is.nsim, main_dir, output_dir, "SSB_status", width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, use.n.years.first, start.years, new_model_names),
      function() plot_ssb_status2(mods, is.nsim, main_dir, output_dir, "SSB_status", width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, use.n.years.first, start.years, new_model_names, base.model),
      function() plot_fbar_status(mods, is.nsim, main_dir, output_dir, "Fbar_status", width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, f.ymin, f.ymax, use.n.years.last, new_model_names),
      function() plot_fbar_status(mods, is.nsim, main_dir, output_dir, "Fbar_status", width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, f.ymin, f.ymax, use.n.years.last, new_model_names, base.model),
      function() plot_fbar_status2(mods, is.nsim, main_dir, output_dir, "Fbar_status", width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, f.ymin, f.ymax, use.n.years.first, start.years, new_model_names),
      function() plot_fbar_status2(mods, is.nsim, main_dir, output_dir, "Fbar_status", width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, f.ymin, f.ymax, use.n.years.first, start.years, new_model_names, base.model),
      function() plot_kobe_status(mods, is.nsim, main_dir, output_dir, width, height, dpi, col.opt, new_model_names, use.n.years.last, show_density = FALSE),
      function() plot_kobe_status(mods, is.nsim, main_dir, output_dir, width, height, dpi, col.opt, new_model_names, use.n.years.last, show_density = TRUE),
      function() plot_model_performance_radar(mods, is.nsim, main_dir, output_dir, width, height, dpi, col.opt, method, use.n.years.first, use.n.years.last, start.years, new_model_names),
      function() plot_model_performance_radar2(mods, is.nsim, main_dir, output_dir, width, height, dpi, col.opt, method, use.n.years.first, use.n.years.last, start.years, new_model_names),
      function() plot_model_performance_radar3(mods, is.nsim, main_dir, output_dir, width, height, dpi, col.opt, method, use.n.years.first, use.n.years.last, start.years, new_model_names),
      function() plot_model_performance_triangle(mods, is.nsim, main_dir, output_dir, width, height, dpi, col.opt, method, new_model_names, use.n.years.first, use.n.years.last, start.years),
      function() plot_model_performance_triangle2(mods, is.nsim, main_dir, output_dir, width, height, dpi, col.opt, method, new_model_names, use.n.years.first, use.n.years.last, start.years, density_ridges = TRUE),
      function() plot_model_performance_triangle2(mods, is.nsim, main_dir, output_dir, width, height, dpi, col.opt, method, new_model_names, use.n.years.first, use.n.years.last, start.years, density_ridges = FALSE),
      function() plot_status_triangle(mods, is.nsim, main_dir, output_dir, width, height, dpi, col.opt, method, new_model_names, use.n.years.first, use.n.years.last, start.years),
      function() plot_model_performance_bar(mods, is.nsim, main_dir, output_dir, new_model_names, width, height, dpi, col.opt, method, use.n.years.first, use.n.years.last, start.years),
      function() plot_model_performance_bar2(mods, is.nsim, main_dir, output_dir, new_model_names, width, height, dpi, col.opt, method, use.n.years.first, use.n.years.last, start.years),
      function() plot_model_performance_bar3(mods, is.nsim, main_dir, output_dir, new_model_names, width, height, dpi, col.opt, method, use.n.years.first, use.n.years.last, start.years),
      function() plot_model_performance_bar4(mods, is.nsim, main_dir, output_dir, new_model_names, width, height, dpi, col.opt, method, use.n.years.first, use.n.years.last, start.years),
      function() plot_AAV_performance(mods, is.nsim, main_dir, output_dir, width, height, dpi, new_model_names, col.opt),
      function() plot_mean_rec_par(mods, is.nsim, main_dir, output_dir, width, height, dpi, col.opt, outlier.opt, new_model_names),
      function() plot_NAA_sigma_par(mods, is.nsim, main_dir, output_dir, width, height, dpi, col.opt, outlier.opt, new_model_names),
      function() plot_terminal_ssb_f_bias_by_assessment(mods, is.nsim, main_dir, output_dir, width, height, dpi, col.opt, plot.style, outlier.opt, show.whisker, new_model_names)
    )
    
    # Run all plot functions
    for (call in plot_calls) {
      tryCatch(call(), error = function(e) message("Error in plot: ", e$message))
    }
    
  } else if (output_format %in% c("html", "pdf")) {
    cat(paste0("\nGenerating ", toupper(output_format), " report...\n"))
    message("\nPlease ensure the biological reference point in the OM is calculated.\n")
    message("Otherwise, the PDF/HTML file will not be available.\n")
    message("But figures will still be saved as PNGs in the output folder.\n")
    
    rmd_path <- file.path(full_output_dir, "mse_report.Rmd")
    report_file <- paste0("mse_report.", output_format)
    report_path <- file.path(full_output_dir, report_file)
    
    # Rmd content
    rmd_lines <- c(
      "---",
      paste0("title: \"MSE Output Report\""),
      paste0("output: ", ifelse(output_format == "html", "html_document", "pdf_document")),
      "params:",
      "  mods: NULL",
      "---",
      "",
      "```{r setup, include=FALSE}",
      "knitr::opts_chunk$set(echo = FALSE, warning = FALSE, message = FALSE)",
      "library(ggplot2); library(dplyr); library(tidyr); library(fmsb)",
      "library(ggpubr); library(gridExtra); library(ggtern)",
      "mods <- params$mods",
      "is.nsim <- if (!is.list(mods[[1]][[1]][[1]])) FALSE else TRUE",
      paste0("output_dir <- '.'"),
      paste0("width <- ", width),
      paste0("height <- ", height),
      paste0("dpi <- ", dpi),
      paste0("col.opt <- '", col.opt, "'"),
      paste0("method <- '", method, "'"),
      paste0("outlier.opt <- ", outlier.opt), 
      paste0("plot.style <- '", plot.style, "'"),
      paste0("show.whisker <- ", deparse(show.whisker)), 
      paste0("f.ymin <- ", f.ymin),
      paste0("f.ymax <- ", f.ymax),
      paste0("start.years <- ", start.years),
      paste0("use.n.years.first <- ", use.n.years.first),
      paste0("use.n.years.last <- ", use.n.years.last),
      paste0("new_model_names <- ", if (is.null(new_model_names)) "NULL" else deparse(new_model_names)),
      paste0("base.model <- '", base.model, "'"),
      "```",
      "",
      "```{=latex}",
      "\\newpage",
      "```",
      "",
      "## Relative Trajectories of SSB and Catch",
      "",
      "```{=latex}",
      "\\newpage",
      "```",
      "```{r trajectories}",
      "plot_relative_trajectories1(mods, is.nsim, '.', output_dir, base.model, new_model_names, width, height, dpi, col.opt)",
      "plot_relative_trajectories2(mods, is.nsim, '.', output_dir, base.model, new_model_names, width, height, dpi, col.opt)",
      "```",

      "```{=latex}",
      "\\newpage",
      "```",
      "",
      "# Realized Time Series of SSB, F, and Catch", 
      "",
      "```{=latex}",
      "\\newpage",
      "```",
      "```{r time_series}",
      "plot_ssb_time_series(mods, is.nsim, '.', output_dir, 'SSB', width, height, dpi, col.opt, new_model_names)",
      "plot_fbar_time_series(mods, is.nsim, '.', output_dir, 'Fbar', width, height, dpi, col.opt, new_model_names)",
      "plot_catch_time_series(mods, is.nsim, '.', output_dir, 'Catch', width, height, dpi, col.opt, new_model_names)",
      "```",
      
      "```{=latex}",
      "\\newpage",
      "```",
      "",
      "# Short-term Performance of SSB, F, and Catch", 
      "",
      "```{=latex}",
      "\\newpage",
      "```",
      "```{r Short-term Performance}",
      "plot_ssb_performance2(mods, is.nsim, '.', output_dir, 'SSB', width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, use.n.years.first, start.years, new_model_names)",
      "plot_fbar_performance2(mods, is.nsim, '.', output_dir, 'Fbar', width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, f.ymin, f.ymax, use.n.years.first, start.years, new_model_names)",
      "plot_catch_performance2(mods, is.nsim, '.', output_dir, 'Catch', width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, use.n.years.first, start.years, new_model_names)",
      "plot_total_catch_performance2(mods, is.nsim, '.', output_dir, 'Catch', width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, use.n.years.first, start.years, new_model_names)",
      "```",

      "```{=latex}",
      "\\newpage",
      "```",
      "",
      "# Short-term Performance of SSB, F, and Catch (Relative)",  
      "",
      "```{=latex}",
      "\\newpage",
      "```",
      "```{r Short-term Relative Performance}",
      "plot_ssb_performance2(mods, is.nsim, '.', output_dir, 'SSB', width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, use.n.years.first, start.years, new_model_names, base.model)",
      "plot_fbar_performance2(mods, is.nsim, '.', output_dir, 'Fbar', width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, f.ymin, f.ymax, use.n.years.first, start.years, new_model_names, base.model)",
      "plot_catch_performance2(mods, is.nsim, '.', output_dir, 'Catch', width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, use.n.years.first, start.years, new_model_names, base.model)",
      "plot_total_catch_performance2(mods, is.nsim, '.', output_dir, 'Catch', width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, use.n.years.first, start.years, new_model_names, base.model)",
      "```",

      "```{=latex}",
      "\\newpage",
      "```",
      "",
      "# Long-term Performance of SSB, F, and Catch", 
      "",
      "```{=latex}",
      "\\newpage",
      "```",
      "```{r Long-term Performance}",
      "plot_ssb_performance(mods, is.nsim, '.', output_dir, 'SSB', width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, use.n.years.last, new_model_names)",
      "plot_fbar_performance(mods, is.nsim, '.', output_dir, 'Fbar', width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, f.ymin, f.ymax, use.n.years.last, new_model_names)",
      "plot_catch_performance(mods, is.nsim, '.', output_dir, 'Catch', width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, use.n.years.last, new_model_names)",
      "plot_total_catch_performance(mods, is.nsim, '.', output_dir, 'Catch', width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, use.n.years.last, new_model_names)",
      "```",
      
      "```{=latex}",
      "\\newpage",
      "```",
      "",
      "# Long-term Performance of SSB, F, and Catch (Relative)", 
      "",
      "```{=latex}",
      "\\newpage",
      "```",
      "```{r Relative Long-term Performance}",
      "plot_ssb_performance(mods, is.nsim, '.', output_dir, 'SSB', width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, use.n.years.last, new_model_names, base.model)",
      "plot_fbar_performance(mods, is.nsim, '.', output_dir, 'Fbar', width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, f.ymin, f.ymax, use.n.years.last, new_model_names, base.model)",
      "plot_catch_performance(mods, is.nsim, '.', output_dir, 'Catch', width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, use.n.years.last, new_model_names, base.model)",
      "plot_total_catch_performance(mods, is.nsim, '.', output_dir, 'Catch', width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, use.n.years.last, new_model_names, base.model)",
      "```",

      "```{=latex}",
      "\\newpage",
      "```",
      "",
      "# Annual Variation in SSB, F, and Catch",  
      "",
      "```{=latex}",
      "\\newpage",
      "```",
      "```{r variation}",
      "plot_ssb_variation(mods, is.nsim, '.', output_dir, 'SSB', width, height, dpi, col.opt, outlier.opt, new_model_names)",
      "plot_fbar_variation(mods, is.nsim, '.', output_dir, 'Fbar', width, height, dpi, col.opt, outlier.opt, new_model_names)",
      "plot_catch_variation(mods, is.nsim, '.', output_dir, 'Catch', width, height, dpi, col.opt, outlier.opt, new_model_names)",
      "```",

      "```{=latex}",
      "\\newpage",
      "```",
      "",
      "# Annual Variation in SSB, F, and Catch (Relative)",  
      "",
      "```{=latex}",
      "\\newpage",
      "```",
      "```{r Relative variation}",
      "plot_ssb_variation(mods, is.nsim, '.', output_dir, 'SSB', width, height, dpi, col.opt, outlier.opt, new_model_names, base.model)",
      "plot_fbar_variation(mods, is.nsim, '.', output_dir, 'Fbar', width, height, dpi, col.opt, outlier.opt, new_model_names, base.model)",
      "plot_catch_variation(mods, is.nsim, '.', output_dir, 'Catch', width, height, dpi, col.opt, outlier.opt, new_model_names, base.model)",
      "```",
      
      "```{=latex}",
      "\\newpage",
      "```",
      "",
      "# Stock Status and Probability of Overfishing and Overfished", 
      "",
      "```{=latex}",
      "\\newpage",
      "```",
      "```{r status}",
      "plot_ssb_status(mods, is.nsim, '.', output_dir, 'SSB_status', width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, use.n.years.last, new_model_names, base.model = NULL)",
      "plot_ssb_status2(mods, is.nsim, '.', output_dir, 'SSB_status', width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, use.n.years.first, start.years, new_model_names, base.model = NULL)",
      "plot_fbar_status(mods, is.nsim, '.', output_dir, 'Fbar_status', width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, f.ymin, f.ymax, use.n.years.last, new_model_names, base.model = NULL)",
      "plot_fbar_status2(mods, is.nsim, '.', output_dir, 'Fbar_status', width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, f.ymin, f.ymax, use.n.years.first, start.years, new_model_names, base.model = NULL)",
      "```",
      
      "```{=latex}",
      "\\newpage",
      "```",
      "",
      "# Stock Status in the Final Years (Relative)",  
      "",
      "```{=latex}",
      "\\newpage",
      "```",
      "```{r Relative status}",
      "plot_ssb_status(mods, is.nsim, '.', output_dir, 'SSB_status', width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, use.n.years.last, new_model_names, base.model, plot_prob = FALSE)",
      "plot_ssb_status2(mods, is.nsim, '.', output_dir, 'SSB_status', width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, use.n.years.first, start.years, new_model_names, base.model, plot_prob = FALSE)",
      "plot_fbar_status(mods, is.nsim, '.', output_dir, 'Fbar_status', width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, f.ymin, f.ymax, use.n.years.last, new_model_names, base.model, plot_prob = FALSE)",
      "plot_fbar_status2(mods, is.nsim, '.', output_dir, 'Fbar_status', width, height, dpi, col.opt, method, outlier.opt, plot.style, show.whisker, f.ymin, f.ymax, use.n.years.first, start.years, new_model_names, base.model, plot_prob = FALSE)",
      "```",
      
      "```{=latex}",
      "\\newpage",
      "```",
      "",
      "# Stock Status in the Final Years (KOBE Plot)", 
      "",
      "```{=latex}",
      "\\newpage",
      "```",
      "```{r KOBE}",
      "plot_kobe_status(mods, is.nsim, '.', output_dir, width, height, dpi, col.opt, new_model_names, use.n.years.last, show_density = FALSE)",
      "plot_kobe_status(mods, is.nsim, '.', output_dir, width, height, dpi, col.opt, new_model_names, use.n.years.last, show_density = TRUE)",
      "```",

      "```{=latex}",
      "\\newpage",
      "```",
      "",
      "# Holistic Performance", 
      "",
      "```{=latex}",
      "\\newpage",
      "```",
      "```{r Holistic Performance}",
      "plot_model_performance_radar(mods, is.nsim, '.', output_dir, width, height, dpi, col.opt, method, use.n.years.first, use.n.years.last, start.years, new_model_names)",
      "plot_model_performance_radar2(mods, is.nsim, '.', output_dir, width, height, dpi, col.opt, method, use.n.years.first, use.n.years.last, start.years, new_model_names)",
      "plot_model_performance_radar3(mods, is.nsim, '.', output_dir, width, height, dpi, col.opt, method, use.n.years.first, use.n.years.last, start.years, new_model_names)",
      "plot_model_performance_bar(mods, is.nsim, '.', output_dir, new_model_names, width, height, dpi, col.opt, method, use.n.years.first, use.n.years.last, start.years)",
      "plot_model_performance_bar2(mods, is.nsim, '.', output_dir, new_model_names, width, height, dpi, col.opt, method, use.n.years.first, use.n.years.last, start.years)",
      "plot_model_performance_bar3(mods, is.nsim, '.', output_dir, new_model_names, width, height, dpi, col.opt, method, use.n.years.first, use.n.years.last, start.years)",
      "plot_model_performance_bar4(mods, is.nsim, '.', output_dir, new_model_names, width, height, dpi, col.opt, method, use.n.years.first, use.n.years.last, start.years)",
      "plot_AAV_performance(mods, is.nsim,'.', output_dir, width, height, dpi, new_model_names, col.opt)",
      "```",
      
      "```{=latex}",
      "\\newpage",
      "```",
      "",
      "# Model Diagnostics", 
      "",
      "```{=latex}",
      "\\newpage",
      "```",
      "```{r Diagnostics}",
      "plot_mean_rec_par(mods, is.nsim, '.', output_dir, width, height, dpi, col.opt, outlier.opt, new_model_names)",
      "plot_NAA_sigma_par(mods, is.nsim, '.', output_dir, width, height, dpi, col.opt, outlier.opt, new_model_names)",
      "```"
    )
    
    writeLines(rmd_lines, con = rmd_path)
    
    # Render report
    rmarkdown::render(
      rmd_path,
      output_file = report_file,
      output_dir = full_output_dir,
      params = list(mods = mods),
      envir = new.env(),
      quiet = TRUE
    )
    cat(paste0("\nSaved ", toupper(output_format), " report to: ", report_path, "\n"))
    if (output_format == "html") browseURL(report_path)
  }
}
