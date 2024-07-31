#' Generate MSE Output Plots and Report
#'
#' This function generates plots for MSE output and creates a report in png/pdf/HTML format.
#'
#' @param mods List of models to compare.
#' @param main.dir Main directory for output files. Default is the current working directory.
#' @param use.n.years Number of years to use for comparison. Default is 10.
#' @param base.mod Base model for comparison. Default is 1.
#' @param short.term List with first.n.years, year_start, and year_end for short-term analysis. if = NULL, results will be summarized for the first 2 years in the feedback loop 
#' @param dpi Resolution for plots. Default is 150.
#' @param out.type Output type for plots.
#'   \itemize{
#'     \item \code{"png"} median-lived species (default)
#'     \item \code{"pdf"} short-lived species
#'     \item \code{"html"} long-lived species
#'     }
#' @return Generates plots in a report.
#' @export
#'
#' @examples
#' plot_mse_output(mods, main.dir = "your/directory/path", out.type = 'pdf')
#' 

plot_mse_output <- function(mods, main.dir = getwd(), 
                            use.n.years = 10, 
                            base.mod = 1, 
                            short.term = NULL,
                            dpi = 150,
                            out.type = 'png') {
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(fmsb)
  library(ggpubr)
  library(gridExtra)
  library(rmarkdown)
  
  sub.dir <- 'Report'
  output_dir <- file.path(main.dir,sub.dir)
  
  if(!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  extract_var <- function(mod, var){
    if (var != "Catch_r") {
      data <- data.frame(mod$om$rep[var])
      data[paste0(var,"_total")] <- rowSums(data)
      if (var == "Catch_s"){
        var = "pred_catch"
        data <- data.frame(mod$om$rep[var])
        data[paste0(var,"_total")] <- rowSums(data)
        names <- gsub("pred_catch","Catch_s",colnames(data))
        colnames(data) <- names
      }
    } else {
      var = "pred_stock_catch"
      data <- mod$om$rep$pred_stock_catch
      Catch_list = list()
      for (i in 1:nrow(data)) {
        Catch_list[[paste0("Catch_r.",i)]] <- apply(data, MARGIN = 3, FUN = colSums)[i,]
      }
      Catch_r.total = apply(data, MARGIN = 3, FUN = sum)
      data <- data.frame(Catch_list,Catch_r.total)
    }
    return(data)
  }
  
  generate_plots <- function(mods, Years, use.n.years) {
    vars <- list(
      list(var = "SSB", title = "SSB", ylab = "SSB"),
      list(var = "Fbar", title = "Fleet-specific F", ylab = "F"),
      list(var = "Catch_s", title = "Stock-specific Catch", ylab = "Catch"),
      list(var = "Catch_r", title = "Region-specific Catch", ylab = "Catch")
    )
    
    All_plots <- lapply(vars, function(v) {
      res = NULL
      for (i in 1:length(mods)){
        for (j in 1:length(mods[[i]])) {
          tmp <- extract_var(mods[[i]][[j]], v$var)
          tmp$nsim <- i
          tmp$Model <- paste("Model", j)
          tmp$Year <- Years
          tmp <- tail(tmp, use.n.years)
          res <- rbind(res, tmp)
        }
      }
      res <- pivot_longer(res, cols = starts_with(v$var), names_to = "Label", values_to = v$var)
      
      p <- ggplot(res, aes(x = Model, y = !! rlang::sym(v$var), col = Model)) +
        geom_boxplot(outlier.shape = NA) +
        facet_grid(Label ~., scales = "free") + 
        scale_colour_brewer(palette = "Set2") + 
        ggtitle(paste0(v$title, " (over last ", use.n.years, " years)")) +
        ylab(v$ylab) +
        theme_bw()
      return(p)
    })
    
    return(All_plots)
  }
  
  save_plots_to_pdf <- function(plot_list, file_path) {
    pdf(file_path, width = 10, height = 7)
    for (plot in plot_list) {
      print(plot)
    }
    dev.off()
  }
  
  generate_html_report <- function(main.dir, sub.dir, All_plots) {
    report_file <- file.path(main.dir, sub.dir, "report.Rmd")
    html_file <- file.path(main.dir, sub.dir, "report.html")
    
    rmd_content <- "
---
title: 'Performance Metrics Report'
output: html_document
params:
  All_plots: !r NULL
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = FALSE)
```

```{r plots, results='asis'}
for (p in params$All_plots) {
  print(p)
}
```"
  
  writeLines(rmd_content, report_file)
  render(report_file, output_file = html_file, params = list(All_plots = All_plots))
  browseURL(html_file)
  }

if (is.null(use.n.years)) use.n.years = 10

if (!is.list(mods[[1]][[1]][[1]])) {
  
  Years = mods[[1]]$om$years
  cat("Comparing models for one realization...\n")
  
  extract_var <- function(mod, var){
    if (var != "Catch_r") {
      data <- data.frame(mod$om$rep[var])
      data[paste0(var,"_total")] <- rowSums(data)
      if (var == "Catch_s"){
        var = "pred_catch"
        data <- data.frame(mod$om$rep[var])
        data[paste0(var,"_total")] <- rowSums(data)
        names <- gsub("pred_catch","Catch_s",colnames(data))
        colnames(data) <- names
      }
    } else {
      var = "pred_stock_catch"
      data <- mod$om$rep$pred_stock_catch
      Catch_list = list()
      for (i in 1:nrow(data)) {
        Catch_list[[paste0("Catch_r.",i)]] <- apply(data, MARGIN = 3, FUN = colSums)[i,]
      }
      Catch_r.total = apply(data, MARGIN = 3, FUN = sum)
      data <- data.frame(Catch_list,Catch_r.total)
    }
    return(data)
  }
  
  var = "SSB"
  res = NULL
  for (i in 1:length(mods)){
    tmp <- extract_var(mods[[i]],var)
    tmp$Model <- paste("Model",i)
    tmp$Year <- Years
    res <- rbind(res, tmp)
  }
  res <- pivot_longer(res,cols = starts_with(var),names_to = "Label", values_to = var)
  
  p1 <- ggplot(res, aes(x = Year, y=!! rlang::sym(var),col = Model)) +
    geom_line(size = 0.8, alpha = 0.8) +
    facet_grid(Label ~., scales = "free") + 
    scale_colour_brewer(palette = "Set2") + 
    ggtitle("SSB") +
    ylab("SSB") +
    theme_bw()
  ggsave(file.path(main.dir,sub.dir,paste0(var,".PNG")), p1, width = 10, height = 7, dpi = dpi)
  
  var = "Fbar"
  res = NULL
  for (i in 1:length(mods)){
    tmp <- extract_var(mods[[i]],var)
    tmp$Model <- paste("Model",i)
    tmp$Year <- Years
    res <- rbind(res, tmp)
  }
  res <- pivot_longer(res,cols = starts_with(var),names_to = "Label", values_to = var)
  
  p2 <- ggplot(res, aes(x = Year, y=!! rlang::sym(var),col = Model)) +
    geom_line(size = 0.8, alpha = 0.8) +
    facet_grid(Label ~., scales = "free") + 
    scale_colour_brewer(palette = "Set2") + 
    ggtitle("Fleet-specific F") +
    ylab("F") +
    theme_bw()
  ggsave(file.path(main.dir,sub.dir,paste0(var,".PNG")), p2, width = 10, height = 7, dpi = dpi)
  
  var = "Catch_s"
  res = NULL
  for (i in 1:length(mods)){
    tmp <- extract_var(mods[[i]],var)
    tmp$Model <- paste("Model",i)
    tmp$Year <- Years
    res <- rbind(res, tmp)
  }
  res <- pivot_longer(res,cols = starts_with(var),names_to = "Label", values_to = var)
  
  p3 <- ggplot(res, aes(x = Year, y=!! rlang::sym(var),col = Model)) +
    geom_line(size = 0.8, alpha = 0.8) +
    facet_grid(Label ~., scales = "free") + 
    scale_colour_brewer(palette = "Set2") + 
    ggtitle("Stock-specific Catch") +
    ylab("Catch") +
    theme_bw()
  ggsave(file.path(main.dir,sub.dir,paste0(var,".PNG")), p3, width = 10, height = 7, dpi = dpi)
  
  var = "Catch_r"
  res = NULL
  for (i in 1:length(mods)){
    tmp <- extract_var(mods[[i]],var)
    tmp$Model <- paste("Model",i)
    tmp$Year <- Years
    res <- rbind(res, tmp)
  }
  res <- pivot_longer(res,cols = starts_with(var),names_to = "Label", values_to = var)
  
  p4 <- ggplot(res, aes(x = Year, y=!! rlang::sym(var),col = Model)) +
    geom_line(size = 0.8, alpha = 0.8) +
    facet_grid(Label ~., scales = "free") + 
    scale_colour_brewer(palette = "Set2") + 
    ggtitle("Region-specific Catch") +
    ylab("Catch") +
    theme_bw()
  ggsave(file.path(main.dir,sub.dir,paste0(var,".PNG")), p4, width = 10, height = 7, dpi = dpi)
  
  # Plot them together
  p <- ggarrange(p1,p2,p3,p4,common.legend = TRUE,legend = "right")
  ggsave(file.path(main.dir,sub.dir,"Performance_Metrics.PNG"), p, width = 15, height = 10, dpi = dpi)
  
  # ------------------------------------------------------------
  # ------------------------------------------------------------
  # ------------------------- Box Plot ------------------------
  # ----------------- SSB, F, Catch_r, Catch_s -----------------
  # ------------------------------------------------------------
  # ------------------------------------------------------------
  
  # Performance metrics summarized over the last n years
  var = "SSB"
  res = NULL
  for (i in 1:length(mods)){
    tmp <- extract_var(mods[[i]],var)
    if (is.null(use.n.years)) {
      tmp$Model <- paste("Model",i)
      tmp$Year <- Years
      tmp <- tail(tmp,5)
      res <- rbind(res, tmp) 
    } else {
      tmp$Model <- paste("Model",i)
      tmp$Year <- Years
      tmp <- tail(tmp,use.n.years)
      res <- rbind(res, tmp) 
    }
  }
  
  res <- pivot_longer(res,cols = starts_with(var),names_to = "Label", values_to = var)
  
  p1 <- ggplot(res, aes(x = Model, y=!! rlang::sym(var),col = Model)) +
    geom_boxplot(outlier.shape = NA) +
    facet_grid(Label ~., scales = "free") + 
    scale_colour_brewer(palette = "Set2") + 
    ggtitle(paste0("SSB"," (over last ",use.n.years," years)")) +
    ylab("SSB") +
    theme_bw()
  ggsave(file.path(main.dir,sub.dir,paste0(var,"_last_",use.n.years,"_years.PNG")), p1, width = 10, height = 7, dpi = dpi)
  
  var = "Fbar"
  res = NULL
  for (i in 1:length(mods)){
    tmp <- extract_var(mods[[i]],var)
    if (is.null(use.n.years)) {
      tmp$Model <- paste("Model",i)
      tmp$Year <- Years
      tmp <- tail(tmp,5)
      res <- rbind(res, tmp) 
    } else {
      tmp$Model <- paste("Model",i)
      tmp$Year <- Years
      tmp <- tail(tmp,use.n.years)
      res <- rbind(res, tmp) 
    }
  }
  
  res <- pivot_longer(res,cols = starts_with(var),names_to = "Label", values_to = var)
  
  p2 <- ggplot(res, aes(x = Model, y=!! rlang::sym(var),col = Model)) +
    geom_boxplot(outlier.shape = NA) +
    facet_grid(Label ~., scales = "free") + 
    scale_colour_brewer(palette = "Set2") + 
    ggtitle(paste0("Fleet-specific F"," (over last ",use.n.years," years)")) +
    ylab("F") +
    theme_bw()
  ggsave(file.path(main.dir,sub.dir,paste0(var,"_last_",use.n.years,"_years.PNG")), p2, width = 10, height = 7, dpi = dpi)
  
  var = "Catch_s"
  res = NULL
  for (i in 1:length(mods)){
    tmp <- extract_var(mods[[i]],var)
    if (is.null(use.n.years)) {
      tmp$Model <- paste("Model",i)
      tmp$Year <- Years
      tmp <- tail(tmp,5)
      res <- rbind(res, tmp) 
    } else {
      tmp$Model <- paste("Model",i)
      tmp$Year <- Years
      tmp <- tail(tmp,use.n.years)
      res <- rbind(res, tmp) 
    }
  }
  
  res <- pivot_longer(res,cols = starts_with(var),names_to = "Label", values_to = var)
  
  p3 <- ggplot(res, aes(x = Model, y=!! rlang::sym(var),col = Model)) +
    geom_boxplot(outlier.shape = NA) +
    facet_grid(Label ~., scales = "free") + 
    scale_colour_brewer(palette = "Set2") + 
    ggtitle(paste0("Stock-specific Catch"," (over last ",use.n.years," years)")) +
    ylab("Catch") +
    theme_bw()
  ggsave(file.path(main.dir,sub.dir,paste0(var,"_last_",use.n.years,"_years.PNG")), p3, width = 10, height = 7, dpi = dpi)
  
  var = "Catch_r"
  res = NULL
  for (i in 1:length(mods)){
    tmp <- extract_var(mods[[i]],var)
    if (is.null(use.n.years)) {
      tmp$Model <- paste("Model",i)
      tmp$Year <- Years
      tmp <- tail(tmp,5)
      res <- rbind(res, tmp) 
    } else {
      tmp$Model <- paste("Model",i)
      tmp$Year <- Years
      tmp <- tail(tmp,use.n.years)
      res <- rbind(res, tmp) 
    }
  }
  res <- pivot_longer(res,cols = starts_with(var),names_to = "Label", values_to = var)
  
  p4 <- ggplot(res, aes(x = Model, y=!! rlang::sym(var),col = Model)) +
    geom_boxplot(outlier.shape = NA) +
    facet_grid(Label ~., scales = "free") + 
    scale_colour_brewer(palette = "Set2") + 
    ggtitle(paste0("Region-specific Catch"," (over last ",use.n.years," years)")) +
    ylab("Catch") +
    coord_cartesian(ylim = quantile(res[[var]],c(0.025,0.975))) + 
    theme_bw()
  
  ggsave(file.path(main.dir,sub.dir,paste0(var,"_last_",use.n.years,"_years.PNG")), p4, width = 10, height = 7, dpi = dpi)
  
  # Plot them together
  p <- ggarrange(p1,p2,p3,p4,common.legend = TRUE,legend = "right")
  ggsave(file.path(main.dir,sub.dir,paste0("Performance_last_",use.n.years,"_years.PNG")), p, width = 15, height = 10, dpi = dpi)
  
  cat("\n----------------------------------\n| Congrats! Your report is done! |\n----------------------------------\n")
  cat(paste0("\nReport has been saved in ",file.path(main.dir,sub.dir)))
  
} else {
  Years = mods[[1]][[1]]$om$years
  cat("Comparing models for multiple realizations...\n")
  tmp1 = mods[[1]][[1]][[1]]
  tmp2 = mods[[1]][[2]][[1]]
  if (sum(tmp1$input$data$agg_catch[1:3] - tmp2$input$data$agg_catch[1:3]) != 0) {
    warning("Results are not comparable because the realization is different.\nPlease check the seed!")
  } else {
    cat("Results are model comparisons for n realizations!\n")
    
    generate_box_plot <- function(var, title, ylab, ylim = NULL) {
      res = NULL
      for (i in 1:length(mods)){
        for (j in 1:length(mods[[1]])){
          tmp <- extract_var(mods[[i]][[j]], var)
          tmp$nsim <- i
          tmp$Model <- paste("Model",j)
          tmp$Year <- Years
          tmp <- tail(tmp, use.n.years)
          res <- rbind(res, tmp)
        }
      }
      res <- pivot_longer(res, cols = starts_with(var), names_to = "Label", values_to = var)
      
      p <- ggplot(res, aes(x = Model, y = !! rlang::sym(var), col = Model)) +
        geom_boxplot(outlier.shape = NA) +
        facet_grid(Label ~., scales = "free") + 
        scale_colour_brewer(palette = "Set2") + 
        ggtitle(paste0(title," (over last ",use.n.years," years)")) +
        ylab(ylab) +
        coord_cartesian(ylim = ylim) + 
        theme_bw()
      return(p)
    }
    
    p1 <- generate_box_plot("SSB", "SSB", "SSB")
    p2 <- generate_box_plot("Fbar", "Fleet-specific F", "F")
    p3 <- generate_box_plot("Catch_s", "Stock-specific Catch", "Catch")
    p4 <- generate_box_plot("Catch_r", "Region-specific Catch", "Catch")
    
    performance_plot <- list(p1, p2, p3, p4)
    
    All_plots <- list()
    
    if (out.type == 'png') {
      filenames <- c("SSB", "Fbar", "Catch_s", "Catch_r")
      mapply(function(p, fname) ggsave(file.path(main.dir, sub.dir, paste0(fname, "_last_", use.n.years, "_years.png")), p, width = 10, height = 7, dpi = dpi), performance_plot, filenames)
      combined_performance_plot <- ggarrange(plotlist = performance_plot, common.legend = TRUE, legend = "right")
      ggsave(file.path(main.dir, sub.dir, paste0("Performance_last_",use.n.years,"_years.png")), combined_performance_plot, width = 15, height = 10, dpi = dpi)
    } else if (out.type == 'pdf') {
      pdf(file.path(main.dir, sub.dir, paste0("Performance_last_",use.n.years,"_years.pdf")), width = 10, height = 7)
      lapply(performance_plot, print)
      dev.off()
    } else if (out.type == 'html') {
      generate_html_report(main.dir, sub.dir, performance_plot)
    }
    
    All_plots <- c(All_plots, performance_plot)
    
    # Additional plots
    cat("\nGenerating additional plots for multiple realizations...\n")
    
    # Box Plot for SSB Stock Status
    res = NULL
    for (i in 1:length(mods)){
      for (j in 1:length(mods[[1]])){
        tmp <- mods[[i]][[j]]$om$rep$SSB
        tmp <- cbind(tmp, rowSums(tmp))
        tmp <- tmp / exp(mods[[i]][[j]]$om$rep$log_SSB_FXSPR)
        tmp <- as.data.frame(tmp)
        names(tmp) <- paste0("SSB/SSB40%.s",1:ncol(tmp))
        names(tmp)[ncol(tmp)] <- "SSB/SSB40%_Global"
        tmp$nsim <- i
        tmp$Model <- paste("Model",j)
        tmp$Year <- Years
        tmp <- tail(tmp, use.n.years)
        res <- rbind(res, tmp)
      }
    }
    var = "SSB/SSB40%"
    res <- pivot_longer(res, cols = starts_with(var), names_to = "Label", values_to = var)
    p1 <- ggplot(res, aes(x = Model, y = !! rlang::sym(var), col = Model)) +
      geom_boxplot(outlier.shape = NA) +
      facet_grid(Label ~., scales = "free") + 
      scale_colour_brewer(palette = "Set2") + 
      ggtitle(paste0("SSB/SSB40%"," (over last ",use.n.years," years)")) +
      ylab("SSB/SSB40%") +
      theme_bw()
    
    ggsave(file.path(main.dir, sub.dir, paste0("SSB_status_last_",use.n.years,"_years.PNG")), p1, width = 10, height = 7, dpi = dpi)
    
    # All_plots <- c(All_plots, list(p1))
    
    # Box Plot for SSB Stock Status (Terminal)
    res = NULL
    for (i in 1:length(mods)){
      for (j in 1:length(mods[[1]])){
        tmp <- mods[[i]][[j]]$om$rep$SSB
        tmp <- cbind(tmp, rowSums(tmp))
        tmp <- tmp / exp(mods[[i]][[j]]$om$rep$log_SSB_FXSPR)
        tmp <- as.data.frame(tmp)
        names(tmp) <- paste0("SSB/SSB40%.s",1:ncol(tmp))
        names(tmp)[ncol(tmp)] <- "SSB/SSB40%_Global"
        tmp$nsim <- i
        tmp$Model <- paste("Model",j)
        tmp$Year <- Years
        tmp <- tail(tmp, 1)
        res <- rbind(res, tmp)
      }
    }
    var = "SSB/SSB40%"
    res <- pivot_longer(res, cols = starts_with(var), names_to = "Label", values_to = var)
    p2 <- ggplot(res, aes(x = Model, y = !! rlang::sym(var), col = Model)) +
      geom_boxplot(outlier.shape = NA) +
      facet_grid(Label ~., scales = "free") + 
      scale_colour_brewer(palette = "Set2") + 
      ggtitle("SSB/SSB40% in the terminal year") +
      ylab("SSB/SSB40%") +
      theme_bw()
    
    ggsave(file.path(main.dir, sub.dir, paste0("SSB_status_in_the_terminal_year.PNG")), p1, width = 10, height = 7, dpi = dpi)
    
    # All_plots <- c(All_plots, list(p2))
    
    # Box Plot for F Stock Status 
    res = NULL
    for (i in 1:length(mods)){
      for (j in 1:length(mods[[1]])){
        Fbar <- mods[[i]][[j]]$om$rep$Fbar
        Fbar <- cbind(Fbar, rowSums(Fbar))
        n_fleets <- mods[[i]][[j]]$om$input$data$n_fleets
        n_regions <- mods[[i]][[j]]$om$input$data$n_regions
        Fbar_age <- max(mods[[i]][[j]]$om$input$data$Fbar_ages)
        Fbar_XSPR <- exp(mods[[i]][[j]]$om$rep$log_FAA_XSPR)[,,Fbar_age]
        Fbar_XSPR <- t(rbind(Fbar_XSPR[1:n_fleets,], Fbar_XSPR[n_fleets + n_regions + 1,]))
        tmp <- as.data.frame(Fbar / Fbar_XSPR)
        names(tmp) <- paste0("F/F40%.r",1:ncol(tmp))
        names(tmp)[ncol(tmp)] <- "F/F40%_Global"
        tmp$nsim <- i
        tmp$Model <- paste("Model",j)
        tmp$Year <- Years
        tmp <- tail(tmp, use.n.years)
        res <- rbind(res, tmp)
      }
    }
    
    var = "F/F40%"
    res <- pivot_longer(res, cols = starts_with(var), names_to = "Label", values_to = var)
    p3 <- ggplot(res, aes(x = Model, y = !! rlang::sym(var), col = Model)) +
      geom_boxplot(outlier.shape = NA) +
      facet_grid(Label ~., scales = "free") + 
      scale_colour_brewer(palette = "Set2") + 
      ggtitle(paste0("F/F40%"," (over last ",use.n.years," years)")) +
      ylab("F/F40%") +
      theme_bw()
    
    ggsave(file.path(main.dir, sub.dir, paste0("F_status_last_",use.n.years,"_years.PNG")), p1, width = 10, height = 7, dpi = dpi)
    
    # All_plots <- c(All_plots, list(p3))
    
    # Box Plot for F Stock Status (Terminal)
    res = NULL
    for (i in 1:length(mods)){
      for (j in 1:length(mods[[1]])){
        Fbar <- mods[[i]][[j]]$om$rep$Fbar
        Fbar <- cbind(Fbar, rowSums(Fbar))
        n_fleets <- mods[[i]][[j]]$om$input$data$n_fleets
        n_regions <- mods[[i]][[j]]$om$input$data$n_regions
        Fbar_age <- max(mods[[i]][[j]]$om$input$data$Fbar_ages)
        Fbar_XSPR <- exp(mods[[i]][[j]]$om$rep$log_FAA_XSPR)[,,Fbar_age]
        Fbar_XSPR <- t(rbind(Fbar_XSPR[1:n_fleets,], Fbar_XSPR[n_fleets + n_regions + 1,]))
        tmp <- as.data.frame(Fbar / Fbar_XSPR)
        names(tmp) <- paste0("F/F40%.r",1:ncol(tmp))
        names(tmp)[ncol(tmp)] <- "F/F40%_Global"
        tmp$nsim <- i
        tmp$Model <- paste("Model",j)
        tmp$Year <- Years
        tmp <- tail(tmp, 1)
        res <- rbind(res, tmp)
      }
    }
    
    var = "F/F40%"
    res <- pivot_longer(res, cols = starts_with(var), names_to = "Label", values_to = var)
    p4 <- ggplot(res, aes(x = Model, y = !! rlang::sym(var), col = Model)) +
      geom_boxplot(outlier.shape = NA) +
      facet_grid(Label ~., scales = "free") + 
      scale_colour_brewer(palette = "Set2") + 
      ggtitle("F/F40% in the terminal year") +
      ylab("F/F40%") +
      theme_bw()
    
    ggsave(file.path(main.dir, sub.dir, paste0("F_status_in_the_terminal_year.PNG")), p1, width = 10, height = 7, dpi = dpi)
    
    stock_status_plots <- list(p1,p2,p3,p4)
    
    if (out.type == 'png') {
      filenames <- c(paste0("SSB_status_last_",use.n.years,"_years"), "SSB_status_in_the_terminal_year", 
                     paste0("F_status_last_",use.n.years,"_years"), "F_status_in_the_terminal_year")
      mapply(function(p, fname) ggsave(file.path(main.dir, sub.dir, paste0(fname, ".png")), p, width = 10, height = 7, dpi = dpi), stock_status_plots, filenames)
      combined_relative_diff_plot <- ggarrange(plotlist = stock_status_plots, common.legend = TRUE, legend = "right")
      ggsave(file.path(main.dir, sub.dir, "Stock_Status.png"), combined_relative_diff_plot, width = 15, height = 10, dpi = dpi)
    } else if (out.type == 'pdf') {
      pdf(file.path(main.dir, sub.dir, "Stock_Status.pdf"), width = 10, height = 7)
      lapply(stock_status_plots, print)
      dev.off()
    } else if (out.type == 'html') {
      generate_html_report(main.dir, sub.dir, stock_status_plots)
    }
    
    All_plots <- c(All_plots,stock_status_plots)
    
    # Radar Chart
    res1 <- NULL
    res2 <- NULL
    res3 <- NULL
    res4 <- NULL
    
    for (i in 1:length(mods)) {
      for (j in 1:length(mods[[1]])) {
        tmp <- extract_var(mods[[i]][[j]], "SSB")
        tmp <- data.frame(t(colSums(tail(tmp, use.n.years))))
        tmp$nsim <- i
        tmp$Model <- paste("Model", j)
        res1 <- rbind(res1, tmp)
        
        tmp <- extract_var(mods[[i]][[j]], "Catch_r")
        tmp <- data.frame(t(colSums(tail(tmp, use.n.years))))
        tmp$nsim <- i
        tmp$Model <- paste("Model", j)
        res2 <- rbind(res2, tmp)
        
        tmp <- mods[[i]][[j]]$om$rep$SSB
        tmp <- cbind(tmp, rowSums(tmp))
        tmp <- tmp / exp(mods[[i]][[j]]$om$rep$log_SSB_FXSPR)
        tmp <- tail(tmp, use.n.years)
        tmp <- apply(tmp, 2, function(x) sum(x > 0.5) / use.n.years)
        tmp <- data.frame(t(tmp))
        names(tmp) <- paste0("PNO_OFED_S", 1:ncol(tmp))
        names(tmp)[ncol(tmp)] <- "PNO_OFED_G"
        tmp$nsim <- i
        tmp$Model <- paste("Model", j)
        res3 <- rbind(res3, tmp)
        
        Fbar <- mods[[i]][[j]]$om$rep$Fbar
        Fbar <- cbind(Fbar, rowSums(Fbar))
        n_fleets <- mods[[i]][[j]]$om$input$data$n_fleets
        n_regions <- mods[[i]][[j]]$om$input$data$n_regions
        Fbar_age <- max(mods[[i]][[j]]$om$input$data$Fbar_ages)
        Fbar_XSPR <- exp(mods[[i]][[j]]$om$rep$log_FAA_XSPR)[,,Fbar_age]
        Fbar_XSPR <- t(rbind(Fbar_XSPR[1:n_fleets,], Fbar_XSPR[n_fleets + n_regions + 1,]))
        tmp <- as.data.frame(Fbar / Fbar_XSPR)
        names(tmp) <- paste0("PNO_OFING_S", 1:ncol(tmp))
        names(tmp)[ncol(tmp)] <- "PNO_OFING_G"
        tmp <- tail(tmp, use.n.years)
        tmp <- apply(tmp, 2, function(x) sum(x < 1) / use.n.years)
        tmp <- data.frame(t(tmp))
        names(tmp) <- paste0("PNO_OFING_S", 1:ncol(tmp))
        names(tmp)[ncol(tmp)] <- "PNO_OFING_G"
        tmp$nsim <- i
        tmp$Model <- paste("Model", j)
        res4 <- rbind(res4, tmp)
      }
    }
    
    res <- rbind(
      pivot_longer(res1, cols = -c(nsim, Model), names_to = "Index", values_to = "Value"),
      pivot_longer(res2, cols = -c(nsim, Model), names_to = "Index", values_to = "Value"),
      pivot_longer(res3, cols = -c(nsim, Model), names_to = "Index", values_to = "Value"),
      pivot_longer(res4, cols = -c(nsim, Model), names_to = "Index", values_to = "Value")
    )
    
    df_agg <- res %>% group_by(Model, Index) %>% summarise(Value = mean(Value), .groups = 'drop')
    df <- data.frame(spread(df_agg, key = Index, value = "Value"))
    rownames(df) <- df$Model
    
    normalize <- function(x) {
      return((x - min(x)) / (max(x) - min(x)))
    }
    
    df.tmp <- as.data.frame(lapply(df[, 2:ncol(df)], normalize))
    df.tmp[is.na(df.tmp)] <- 1
    
    max_values <- rep(1, ncol(df))
    min_values <- rep(0, ncol(df))
    
    legend.names <- rownames(df)
    df <- data.frame(Model = rownames(df), df.tmp)
    df <- rbind(max_values, min_values, df)
    row.names(df) <- c("Maximum", "Minimum", legend.names)
    
    my_colors <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")[1:length(legend.names)]
    colors_border <- c(rgb(1, 1, 1, 0), rgb(1, 1, 1, 0), my_colors)
    
    # Create a radar chart plot and store it
    radar_chart_plot <- function() {
      radarchart(
        df[, -1], axistype = 0, maxmin = F,
        plwd = 4, plty = 1, pcol = colors_border,
        cglcol = "grey", cglty = 1, axislabcol = "black", cglwd = 0.8,
        vlcex = 0.8
      )
      legend(x = -1.4, y = 1.25, legend = legend.names, col = my_colors,
             pch = 20, text.col = "grey", cex = 0.8, pt.cex = 2)
    }
    
    # Save radar chart as a plot object
    radar_plot <- function() {
      png(file.path(main.dir, sub.dir, "Radar_Chart.png"), width = 7, height = 7, units = "in", res = dpi)
      radar_chart_plot()
      dev.off()
      
      # Capture the radar chart plot
      radar_chart_plot()
    }
    
    # Create radar plot
    radar_plot()
    
    # Include radar chart in the list of plots for PDF
    All_plots <- c(All_plots, list(radar_chart_plot))
    
    # Save all plots to a single PDF
    save_plots_to_pdf <- function(plots, pdf_file) {
      pdf(pdf_file, width = 10, height = 7)
      for (plot in plots) {
        if (is.function(plot)) {
          plot()
        } else {
          print(plot)
        }
      }
      dev.off()
    }
    
    if (out.type == 'pdf') {
      pdf_file <- file.path(main.dir, sub.dir, "Performance_Metrics.pdf")
      save_plots_to_pdf(All_plots, pdf_file)
      cat(paste0("\nReport has been saved in ", pdf_file, "\n"))
    }
    
    # Calculate Relative Difference
    generate_relative_diff_plot <- function(var, title) {
      res = NULL
      for (i in 1:length(mods)){
        for (j in 1:length(mods[[1]])){
          tmp  <- extract_var(mods[[i]][[j]],var)
          true <- extract_var(mods[[i]][[base.mod]],var)
          tmp  <- tmp/true - 1
          tmp$nsim <- i
          tmp$Model <- paste("Model", j)
          tmp$Year <- Years
          tmp <- tail(tmp, use.n.years)
          res <- rbind(res, tmp)
        }
      }
      res <- pivot_longer(res, cols = starts_with(var), names_to = "Label", values_to = var)
      
      p <- ggplot(res, aes(x = Model, y = !! rlang::sym(var), col = Model)) +
        geom_boxplot(outlier.shape = NA) +
        facet_grid(Label ~., scales = "free") + 
        scale_colour_brewer(palette = "Set2") + 
        ggtitle(paste0(title, " (over last ", use.n.years, " years)")) +
        ylab("Relative Difference") +
        theme_bw()
      return(p)
    }
    
    p1 <- generate_relative_diff_plot("SSB", "SSB")
    p2 <- generate_relative_diff_plot("Fbar", "Fleet-specific F")
    p3 <- generate_relative_diff_plot("Catch_s", "Stock-specific Catch")
    p4 <- generate_relative_diff_plot("Catch_r", "Region-specific Catch")
    
    relative_diff_plots <- list(p1, p2, p3, p4)
    
    if (out.type == 'png') {
      filenames <- c("SSB_diff", "Fbar_diff", "Catch_s_diff", "Catch_r_diff")
      mapply(function(p, fname) ggsave(file.path(main.dir, sub.dir, paste0(fname, "_last_", use.n.years, "years.png")), p, width = 10, height = 7, dpi = dpi), relative_diff_plots, filenames)
      combined_relative_diff_plot <- ggarrange(plotlist = relative_diff_plots, common.legend = TRUE, legend = "right")
      ggsave(file.path(main.dir, sub.dir, "Performance_Relative_Diff.png"), combined_relative_diff_plot, width = 15, height = 10, dpi = dpi)
    } else if (out.type == 'pdf') {
      pdf(file.path(main.dir, sub.dir, "Performance_Relative_Diff.pdf"), width = 10, height = 7)
      lapply(relative_diff_plots, print)
      dev.off()
    } else if (out.type == 'html') {
      generate_html_report(main.dir, sub.dir, relative_diff_plots)
    }
    
    All_plots <- c(All_plots, relative_diff_plots)
    
    # KOBE Plot
    res1 <- NULL
    res2 <- NULL
    
    for (i in 1:length(mods)){
      for (j in 1:length(mods[[1]])){
        tmp <- mods[[i]][[j]]$om$rep$SSB
        tmp <- cbind(tmp,rowSums(tmp))
        tmp <- data.frame(tmp/exp(mods[[i]][[j]]$om$rep$log_SSB_FXSPR))
        tmp <- tail(tmp, 1)
        names(tmp) <- paste0("S", 1:ncol(tmp))
        tmp$nsim <- i
        tmp$Model <- paste("Model", j)
        res1 <- rbind(res1, tmp)
        
        Fbar <- mods[[i]][[j]]$om$rep$Fbar
        Fbar_XSPR = exp(mods[[i]][[j]]$om$rep$log_FAA_XSPR)[,,max(mods[[i]][[j]]$om$input$data$Fbar_ages)]
        Fbar_XSPR = t(rbind(Fbar_XSPR[1:ncol(Fbar),],Fbar_XSPR[nrow(Fbar_XSPR),]))
        Fbar <- cbind(Fbar,rowSums(Fbar))
        tmp <- data.frame(Fbar/Fbar_XSPR)
        names(tmp) <- paste0("S",1:ncol(tmp))
        tmp$nsim <- i
        tmp$Model <- paste("Model", j)
        tmp <- tail(tmp, 1)
        res2 <- rbind(res2, tmp)
      }
    }
    
    temp1 <- pivot_longer(res1, cols = starts_with("S"), names_to = "Index", values_to = "Overfished")
    temp2 <- pivot_longer(res2, cols = starts_with("S"), names_to = "Index", values_to = "Overfishing")
    
    temp <- cbind(temp1, temp2)
    temp <- temp %>% select(unique(names(.))) %>% select(-nsim)
    
    df2 <- temp[temp$Model == paste0("Model ", base.mod),]
    
    n.col = length(unique(temp$Model))
    my_colors <- c("#66C2A5", "#FC8D62", "#8DA0CB", "#E78AC3", "#A6D854", "#FFD92F", "#E5C494", "#B3B3B3")[1:n.col]
    
    p <- ggplot(data = temp, aes(x = Overfished ,y = Overfishing)) + 
      facet_wrap(~ Index) +
      geom_point(aes(color = Model), size = 3, alpha = 0.8) +
      scale_color_manual(values=my_colors) +
      annotate('rect', xmin = 0.5, xmax = 100, ymin = -100, ymax = 1, alpha = 0.2, fill = "yellow") +
      theme_bw() +
      xlab(bquote(paste("SSB/", SSB[paste(.(40),"%")]))) +
      ylab(bquote(paste("F/", F[paste(.(40),"%")]))) +
      ggtitle("Stock Status in the Terminal Year") +
      theme(axis.text = element_text(size = 12),
            axis.title = element_text(size = 20),
            plot.title = element_text(size = 12),
            strip.text = element_text(size = 12),
            legend.text = element_text(size = 12),
            legend.title = element_text(size = 12),
            aspect.ratio = 1) +
      theme(strip.text.x = element_text(size = 12, color = "black")) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) +
      geom_vline(xintercept = 0.5, linetype = "dashed", color = "red", size = 1) + 
      geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 1)
    
    ylim = c(0,quantile(temp[['Overfishing']], 0.95))
    xlim = c(0,quantile(temp[['Overfished']], 0.95))
    p <- p + coord_cartesian(ylim = ylim, xlim = xlim)
    
    ggsave(file.path(main.dir, sub.dir, "Population_Status_all.png"), p, width = 15, height = 10, dpi = dpi)
    
    All_plots <- c(All_plots, list(p))
    
    # Calculate First N years
    if(!is.null(short.term)) {
      first.n.years = short.term$first.n.years
      year_start = short.term$year_start
      year_end = short.term$year_end
      cat("\nShort-term performance is now calculated! Double check if year_start and year_end are specified correctly!\n")
      cat(paste0("\nFirst year of burn-in period is ",year_start,"\n"))
      cat(paste0("\nLast year of burn-in period is ",year_end,"\n"))
      cat(paste0("\nFirst ",first.n.years," years in the feedback loop are used here\n"))
      base.years = year_end - year_start + 1
    } else {
      tmp = which(mods[[1]][[1]]$om$rep$SSB[,1] - mods[[1]][[2]]$om$rep$SSB[,1] != 0)[1]
      year_start = Years[1]
      first.n.years = 2
      year_end = Years[tmp]
      base.years = year_end - year_start + 1
    }
    
    generate_first_n_years_plot <- function(var, title) {
      res = NULL
      for (i in 1:length(mods)){
        for (j in 1:length(mods[[1]])){
          tmp  <- extract_var(mods[[i]][[j]], var)
          true <- extract_var(mods[[i]][[base.mod]], var)
          tmp  <- tmp/true - 1
          tmp$nsim <- i
          tmp$Model <- paste("Model", j)
          tmp$Year <- Years
          tmp <- tmp[base.years + 1:first.n.years,]
          res <- rbind(res, tmp)
        }
      }
      res <- pivot_longer(res, cols = starts_with(var), names_to = "Label", values_to = var)
      
      p <- ggplot(res, aes(x = Model, y = !! rlang::sym(var), col = Model)) +
        geom_boxplot(outlier.shape = NA) +
        facet_grid(Label ~., scales = "free") + 
        scale_colour_brewer(palette = "Set2") + 
        ggtitle(paste0(title, " (over first ", first.n.years, " years)")) +
        ylab("Relative Difference") +
        theme_bw()
      return(p)
    }
    
    p1 <- generate_first_n_years_plot("SSB", "SSB")
    p2 <- generate_first_n_years_plot("Fbar", "Fleet-specific F")
    p3 <- generate_first_n_years_plot("Catch_s", "Stock-specific Catch")
    p4 <- generate_first_n_years_plot("Catch_r", "Region-specific Catch")
    
    first_n_years_plots <- list(p1, p2, p3, p4)
    
    if (out.type == 'png') {
      filenames <- c("SSB_first_n_years", "Fbar_first_n_years", "Catch_s_first_n_years", "Catch_r_first_n_years")
      mapply(function(p, fname) ggsave(file.path(main.dir, sub.dir, paste0(fname, "_first_", first.n.years, "years.png")), p, width = 10, height = 7, dpi = dpi), first_n_years_plots, filenames)
      combined_first_n_years_plot <- ggarrange(plotlist = first_n_years_plots, common.legend = TRUE, legend = "right")
      ggsave(file.path(main.dir, sub.dir, "Performance_First_N_Years.png"), combined_first_n_years_plot, width = 15, height = 10, dpi = dpi)
    } else if (out.type == 'pdf') {
      pdf(file.path(main.dir, sub.dir, "Performance_First_N_Years.pdf"), width = 10, height = 7)
      lapply(first_n_years_plots, print)
      dev.off()
    } else if (out.type == 'html') {
      generate_html_report(main.dir, sub.dir, first_n_years_plots)
    }
    
    All_plots <- c(All_plots, first_n_years_plots)
    
    # Simulation-Estimation: Mean_rec, Sigma
    res <- NULL
    for (i in 1:length(mods)) {
      for (j in 1:length(mods[[1]])) {
        tmp <- mods[[i]][[j]]
        k <- length(tmp$par.est)
        if (any(names(tmp$par.est[[k]]) == "mean_rec_pars")) {
          temp <- exp(tmp$par.est[[k]]$mean_rec_pars[, 1])
        } else {
          m <- length(tmp$par.est[[k]])
          temp <- NULL
          for (n in 1:m) {
            temp1 <- exp(tmp$par.est[[k]][[n]]$mean_rec_pars[, 1])
            temp <- c(temp, temp1)
          }
        }
        nsim <- i
        Model <- j
        res1 <- data.frame(Model, nsim, Value = temp)
        res1$Var <- paste0("Mean_Rec_", 1:length(temp))
        if (length(temp) == 1) res1$Var <- paste0("Mean_Rec")
        res <- rbind(res, res1)
      }
    }
    res$Model <- as.factor(res$Model)
    
    p1 <- ggplot(res, aes(x = Model, y = Value, col = Model)) +
      geom_boxplot(outlier.shape = NA) +
      facet_grid(Var ~., scales = "free") + 
      scale_colour_brewer(palette = "Set2") + 
      ggtitle(paste0("Mean Recruitment from the Last EM")) +
      ylab("") +
      theme_bw() +
      theme(axis.text = element_text(size = 10),
            axis.title = element_text(size = 15),
            plot.title = element_text(size = 12),
            strip.text = element_text(size = 8),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 10)) +
      theme(strip.text.x = element_text(size = 10, color = "black"))
    ggsave(file.path(main.dir, sub.dir, "Mean_rec_par.png"), p1, width = 10, height = 10, dpi = dpi)
    
    res <- NULL
    for (i in 1:length(mods)) {
      for (j in 1:length(mods[[1]])) {
        tmp <- mods[[i]][[j]]
        k <- length(tmp$par.est)
        if (any(names(tmp$par.est[[k]]) == "log_NAA_sigma")) {
          rec_sig <- exp(tmp$par.est[[k]]$log_NAA_sigma[, 1, 1])
          naa_sig <- exp(tmp$par.est[[k]]$log_NAA_sigma[, 1, 2])
          temp <- c(rec_sig, naa_sig)
        } else {
          m <- length(tmp$par.est[[k]])
          temp <- NULL
          for (n in 1:m) {
            rec_sig <- exp(tmp$par.est[[k]][[n]]$log_NAA_sigma[, 1, 1])
            naa_sig <- exp(tmp$par.est[[k]][[n]]$log_NAA_sigma[, 1, 2])
            temp <- c(rec_sig, naa_sig)
          }
        }
        nsim <- i
        Model <- j
        res1 <- data.frame(Model, nsim, Value = temp)
        if (length(temp) == 2) {
          res1$Var <- c("Rec_sigma", "NAA_sigma")
        } else {
          res1$Var <- c(paste0("Rec_sigma", 1:length(rec_sig)), paste0("NAA_sigma", 1:length(naa_sig)))
        }
        res <- rbind(res, res1)
      }
    }
    res$Model <- as.factor(res$Model)
    
    p2 <- ggplot(res, aes(x = Model, y = Value, col = Model)) +
      geom_boxplot(outlier.shape = NA) +
      facet_grid(Var ~., scales = "free") + 
      scale_colour_brewer(palette = "Set2") + 
      ggtitle(paste0("Variance of NAA from the Last EM")) +
      ylab("") +
      theme_bw() +
      theme(axis.text = element_text(size = 10),
            axis.title = element_text(size = 15),
            plot.title = element_text(size = 12),
            strip.text = element_text(size = 8),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 10)) +
      theme(strip.text.x = element_text(size = 10, color = "black"))
    ggsave(file.path(main.dir, sub.dir, "Variance_Para_NAA.png"), p2, width = 10, height = 15, dpi = dpi)
    
    # AIC, Pearson Residuals
    res <- NULL
    for (i in 1:length(mods)){
      for (j in 1:length(mods[[1]])){
        tmp <- mods[[i]][[j]]
        k <- length(tmp$par.est)
        if(any(names(tmp$par.est[[k]]) == "mean_rec_pars")){
          temp <- -2*(tmp$opt_list[[k]]$obj + length(tmp$opt_list[[k]]$par))
        } else {
          m <- length(tmp$par.est[[k]])
          temp <- NULL
          for (n in 1:m) {
            temp1 <- -2*(tmp$opt_list[[k]][[n]]$obj + length(tmp$opt_list[[k]][[n]]$par))
            temp <- sum(temp,temp1)
          }
        }
        nsim <- i
        Model <- j
        res1 <- data.frame(Model, nsim, Value = temp)
        res1$Var <- paste0("AIC_", 1:length(temp))
        if(length(temp) == 1) res1$Var <- paste0("AIC")
        res <- rbind(res, res1)
      }
    }
    res$Model <- as.factor(res$Model)
    
    p3 <- ggplot(res, aes(x = Model, y = Value, col = Model)) +
      geom_boxplot(outlier.shape = NA) +
      facet_grid(Var ~., scales = "free") + 
      scale_colour_brewer(palette = "Set2") + 
      ggtitle(paste0("AIC")) +
      ylab("") +
      theme_bw() +
      theme(axis.text = element_text(size = 10),
            axis.title = element_text(size = 15),
            plot.title = element_text(size = 12),
            strip.text = element_text(size = 8),
            legend.text = element_text(size = 10),
            legend.title = element_text(size = 10)) +
      theme(strip.text.x = element_text(size = 10, color = "black"))
    ggsave(file.path(main.dir, sub.dir, "AIC.png"), p3, width = 5, height = 5, dpi = dpi)
    
    sim_estimation_plots <- list(p1, p2, p3)
    
    if (out.type == 'png') {
      filenames <- c("Mean_rec_par","Variance_Para_NAA","AIC")
      mapply(function(p, fname) ggsave(file.path(main.dir, sub.dir, paste0(fname, ".png")), p, width = 10, height = 10, dpi = dpi), sim_estimation_plots, filenames)
      combined_sim_estimation_plot <- ggarrange(plotlist = sim_estimation_plots, common.legend = TRUE, legend = "right")
      ggsave(file.path(main.dir, sub.dir, "Diagnostic_Results.png"), combined_sim_estimation_plot, width = 15, height = 10, dpi = dpi)
    } else if (out.type == 'pdf') {
      pdf(file.path(main.dir, sub.dir, "Diagnostic_Results.pdf"), width = 10, height = 7)
      lapply(sim_estimation_plots, print)
      dev.off()
    } else if (out.type == 'html') {
      generate_html_report(main.dir, sub.dir, sim_estimation_plots)
    }
    
    All_plots <- c(All_plots, sim_estimation_plots)
  }
  
  if (out.type == 'pdf') {
    pdf_file <- file.path(main.dir, sub.dir, "Performance_Metrics.pdf")
    save_plots_to_pdf(All_plots, pdf_file)
    cat(paste0("\nReport has been saved in ", pdf_file, "\n"))
  } else if (out.type == 'html') {
    All_plots = All_plots[-9] # remove radar chart here!
    generate_html_report(main.dir, sub.dir, All_plots)
    cat(paste0("\nHTML report has been generated at ", file.path(main.dir, sub.dir, "report.html"), "\n"))
  } else {
    cat("Individual png figure is produced and the full report will not be provided when out.type = png is used!")
    # filenames <- sapply(All_plots, function(p, var) file.path(main.dir, sub.dir, paste0(var, ".png")), vars)
    # mapply(ggsave, filenames, All_plots, MoreArgs = list(width = 10, height = 7, dpi = dpi))
    # combined_box_plot <- ggarrange(plotlist = All_plots, common.legend = TRUE, legend = "right")
    # ggsave(file.path(main.dir, sub.dir, "Performance_Full_Report.png"), combined_box_plot, width = 15, height = 10, dpi = dpi)
  }
  cat("\n----------------------------------\n| Congrats! Your report is done! |\n----------------------------------\n")
  cat(paste0("\nReport has been saved in ",file.path(main.dir,sub.dir)))
}
}

# Use the function like this:
# plot_mse_output(mods, main.dir = "your/directory/path", out.type = 'pdf')
