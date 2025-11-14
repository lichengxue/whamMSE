plot_ssb_time_series <- function(mods, is.nsim, main.dir, sub.dir, var = "SSB",
                                 width = 10, height = 7, dpi = 300, col.opt = "D",
                                 new_model_names = NULL) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)

  if (!is.nsim) {

    Years = mods[[1]]$om$years

    # Only one realization
    res <- lapply(seq_along(mods), function(i) {
      data.frame(
        SSB = mods[[i]]$om$rep$SSB,
        Model = paste0("Model", i),
        Year = Years,
        Realization = 1  # default realization
      )
    }) %>% bind_rows()
  } else {

    Years = mods[[1]][[1]]$om$years

    # Multiple realizations
    res <- lapply(seq_along(mods), function(r) {
      lapply(seq_along(mods[[r]]), function(m) {
        data.frame(
          SSB = mods[[r]][[m]]$om$rep$SSB,
          Model = paste0("Model", m),
          Year = Years,
          Realization = rep(r, length(Years))
        )
      }) %>% bind_rows()
    }) %>% bind_rows()

  }

  # Allow renaming Model names
  if (!is.null(new_model_names)) {
    if (length(new_model_names) != length(unique(res$Model))) {
      stop("Length of new_model_names must match the number of models.")
    }
    res$Model <- factor(res$Model,
                        levels = paste0("Model", seq_along(new_model_names)),
                        labels = new_model_names)
  }

  # Pivot longer if needed (for multiple SSB types)
  res <- pivot_longer(res, cols = starts_with("SSB"), names_to = "Label", values_to = "SSB")

  # Plot
  p <- ggplot(res, aes(x = Year, y = SSB, color = Model, group = interaction(Model, Realization))) +
    geom_line(size = 0.3, alpha = 0.5) +
    facet_grid(Label ~ ., scales = "free") +
    scale_color_viridis_d(option = col.opt) +
    ggtitle(var) +
    ylab("SSB") +
    theme_bw()

  # Create the new subfolder if it doesn't exist
  new_sub_dir <- file.path(main.dir, sub.dir, "Time_Series")
  
  if (!file.exists(new_sub_dir)){
    dir.create(new_sub_dir)
  }
  
  # Save the figure inside the new subfolder
  ggsave(file.path(new_sub_dir, paste0(var, ".PNG")), plot = p, width = width, height = height, dpi = dpi)
  
  return(p)  # Return the plot if you want to print or modify later
}

plot_fbar_time_series <- function(mods, is.nsim, main.dir, sub.dir, var = "Fbar",
                                  width = 10, height = 7, dpi = 300, col.opt = "D",
                                  new_model_names = NULL) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)

  make_plot_data <- function(index_range, label_prefix) {

    if (!is.nsim) {

      Years = mods[[1]]$om$years

      # Only one realization
      lapply(seq_along(mods), function(i) {
        tmp <- mods[[i]]$om$rep$Fbar[, index_range, drop = FALSE]
        tmp <- as.data.frame(tmp)
        names(tmp) <- paste0(label_prefix, seq_along(index_range))
        tmp$Model <- paste0("Model", i)
        tmp$Year <- Years
        tmp$Realization <- 1
        tmp
      }) %>% bind_rows()
    } else {

      Years = mods[[1]][[1]]$om$years

      # Multiple realizations
      lapply(seq_along(mods), function(r) {
        lapply(seq_along(mods[[r]]), function(m) {
          tmp <- mods[[r]][[m]]$om$rep$Fbar[, index_range, drop = FALSE]
          tmp <- as.data.frame(tmp)
          names(tmp) <- paste0(label_prefix, seq_along(index_range))
          tmp$Model <- paste0("Model", m)
          tmp$Year <- Years
          tmp$Realization <- r
          tmp
        }) %>% bind_rows()
      }) %>% bind_rows()
    }
  }

  if (!is.nsim) {
    # Get number of fleets and regions
    n_fleets <- mods[[1]]$om$input$data$n_fleets[1]
    n_regions <- mods[[1]]$om$input$data$n_regions[1]
  } else {
    # Get number of fleets and regions
    n_fleets <- mods[[1]][[1]]$om$input$data$n_fleets[1]
    n_regions <- mods[[1]][[1]]$om$input$data$n_regions[1]
  }

  # Fleet-level Fbar
  res_fleet <- make_plot_data(1:n_fleets, "Fleet_")

  # Region-level Fbar
  res_region <- make_plot_data((n_fleets + 1):(n_fleets + n_regions), "Region_")

  # Global Fbar (only one column)
  res_global <- make_plot_data(n_fleets + n_regions + 1, "Global")

  # Function to handle renaming and plotting
  plot_and_save <- function(res, title, ylab_text, filename) {
    if (!is.null(new_model_names)) {
      if (length(new_model_names) != length(unique(res$Model))) {
        stop("Length of new_model_names must match the number of models.")
      }
      res$Model <- factor(res$Model,
                          levels = paste0("Model", seq_along(new_model_names)),
                          labels = new_model_names)
    }

    res_long <- pivot_longer(res, cols = starts_with(c("Fleet_", "Region_", "Global")),
                             names_to = "Label", values_to = "Fbar")

    p <- ggplot(res_long, aes(x = Year, y = Fbar, color = Model, group = interaction(Model, Realization))) +
      geom_line(size = 0.3, alpha = 0.5) +
      facet_grid(Label ~ ., scales = "free") +
      scale_color_viridis_d(option = col.opt) +
      ggtitle(title) +
      ylab(ylab_text) +
      theme_bw()
    
    # Create the new subfolder if it doesn't exist
    new_sub_dir <- file.path(main.dir, sub.dir, "Time_Series")
    
    if (!file.exists(new_sub_dir)){
      dir.create(new_sub_dir)
    }
    
    # Save the figure inside the new subfolder
    ggsave(file.path(new_sub_dir, paste0(filename, ".PNG")), plot = p, width = width, height = height, dpi = dpi)

    return(p)
  }

  # Create and save each plot
  p_fleet <- plot_and_save(res_fleet, "Fbar by Fleet", "Fbar", paste0(var, "_fleet"))
  p_region <- plot_and_save(res_region, "Fbar by Region", "Fbar", paste0(var, "_region"))
  p_global <- plot_and_save(res_global, "Global Fbar", "Fbar", paste0(var, "_global"))

  return(list(fleet = p_fleet, region = p_region, global = p_global))
}

plot_catch_time_series <- function(mods, is.nsim, main.dir, sub.dir, var = "Catch",
                                   width = 10, height = 7, dpi = 300, col.opt = "D",
                                   new_model_names = NULL) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  # Prepare data
  if (!is.nsim) {
    
    Years = mods[[1]]$om$years
    
    # Only one realization
    res <- lapply(seq_along(mods), function(i) {
      data.frame(
        Catch = mods[[i]]$om$rep$pred_catch,
        Model = paste0("Model", i),
        Year = Years,
        Realization = 1
      )
    }) %>% bind_rows()
  } else {
    
    Years = mods[[1]][[1]]$om$years
    
    # Multiple realizations
    res <- lapply(seq_along(mods), function(r) {
      lapply(seq_along(mods[[r]]), function(m) {
        data.frame(
          Catch = mods[[r]][[m]]$om$rep$pred_catch,
          Model = paste0("Model", m),
          Year = Years,
          Realization = r
        )
      }) %>% bind_rows()
    }) %>% bind_rows()
  }
  
  # Allow renaming Model names
  if (!is.null(new_model_names)) {
    if (length(new_model_names) != length(unique(res$Model))) {
      stop("Length of new_model_names must match the number of models.")
    }
    res$Model <- factor(res$Model,
                        levels = paste0("Model", seq_along(new_model_names)),
                        labels = new_model_names)
  }
  
  # Pivot longer if needed (for multiple Catch types)
  res <- pivot_longer(res, cols = starts_with("Catch"), names_to = "Label", values_to = "Catch")
  
  # Plot
  p <- ggplot(res, aes(x = Year, y = Catch, color = Model, group = interaction(Model, Realization))) +
    geom_line(size = 0.3, alpha = 0.5) +
    facet_grid(Label ~ ., scales = "free") +
    scale_color_viridis_d(option = col.opt) +
    ggtitle(var) +
    ylab("Catch") +
    theme_bw()
  
  # Create the new subfolder if it doesn't exist
  new_sub_dir <- file.path(main.dir, sub.dir, "Time_Series")
  
  if (!file.exists(new_sub_dir)){
    dir.create(new_sub_dir)
  }
  
  # Save the figure inside the new subfolder
  ggsave(file.path(new_sub_dir, paste0(var, ".PNG")), plot = p, width = width, height = height, dpi = dpi)
  
  return(p)  # Return the plot if you want to print or modify later
}

plot_ssb_performance <- function(mods, is.nsim, main.dir, sub.dir, var = "SSB",
                                 width = 10, height = 7, dpi = 300, col.opt = "D",
                                 method = NULL,
                                 outlier.opt = NA,
                                 plot.style = "median_iqr",
                                 show.whisker = TRUE,
                                 use.n.years = NULL,
                                 new_model_names = NULL,
                                 base.model = NULL) {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(rlang)
  
  if (is.null(use.n.years)) {
    cat("\nuse.n.years is not specified, so default (5 years) is used here!\n")
    use.n.years <- 5
  }
  
  Years <- if (!is.nsim) mods[[1]]$om$years else mods[[1]][[1]]$om$years
  
  # 1) Build wide-format res for LAST N years (same structure as your original code)
  res <- if (!is.nsim) {
    lapply(seq_along(mods), function(i) {
      data.frame(
        SSB = mods[[i]]$om$rep$SSB,   # <-- IMPORTANT: keep this exactly like your working code
        Model = paste0("Model", i),
        Year = Years,
        Realization = 1
      ) %>% 
        tail(use.n.years)
    }) %>% bind_rows()
  } else {
    lapply(seq_along(mods), function(r) {
      lapply(seq_along(mods[[r]]), function(m) {
        data.frame(
          SSB = mods[[r]][[m]]$om$rep$SSB,  # same here
          Model = paste0("Model", m),
          Year = Years,
          Realization = r
        ) %>% 
          tail(use.n.years)
      }) %>% bind_rows()
    }) %>% bind_rows()
  }
  
  # 2) Add total/global SSB in WIDE format (exactly like performance2)
  #    At this point res may have SSB, SSB.1, SSB.2, ...; all of these start with "SSB"
  res <- res %>%
    rowwise() %>%
    mutate(SSB_Global = sum(c_across(starts_with("SSB")), na.rm = TRUE)) %>%
    ungroup()
  
  # 3) Handle model renaming (still wide)
  if (!is.null(new_model_names)) {
    if (length(new_model_names) != length(unique(res$Model))) {
      stop("Length of new_model_names must match the number of models.")
    }
    res$Model <- factor(
      res$Model,
      levels = paste0("Model", seq_along(new_model_names)),
      labels = new_model_names
    )
    
    if (!is.null(base.model)) {
      if (!(base.model %in% new_model_names)) {
        warning("base.model does not match any of the new_model_names. Make sure it's the renamed version (e.g., 'M1').")
      }
    }
  }
  
  # 4) Pivot longer so regions + SSB_Global all become Labels
  res <- pivot_longer(
    res,
    cols = starts_with(var),  # var = "SSB": includes SSB, SSB.1, ..., SSB_Global
    names_to  = "Label",
    values_to = var
  )
  
  # 5) Compute relative difference from base.model if specified
  if (!is.null(base.model)) {
    base_df <- res %>% 
      filter(Model == base.model) %>%
      rename(base_val = !!sym(var)) %>%
      select(Realization, Year, Label, base_val)
    
    res <- res %>%
      left_join(base_df, by = c("Realization", "Year", "Label")) %>%
      mutate(!!var := (!!sym(var)) / base_val - 1)
  }
  
  # 6) Apply mean or median summarization across the last N years if method is specified
  if (!is.null(method)) {
    res <- res %>%
      group_by(Model, Realization, Label) %>%
      summarise(
        !!var := if (method == "mean") mean(!!sym(var)) else median(!!sym(var)),
        .groups = "drop"
      )
  }
  
  # 7) Plot
  if (plot.style == "boxplot") {
    
    p <- ggplot(res, aes(x = Model, y = !!sym(var), color = Model)) +
      geom_boxplot(lwd = 0.8, outlier.shape = outlier.opt) +
      facet_grid(Label ~ ., scales = "free") +
      scale_color_viridis_d(option = col.opt) +
      ggtitle(paste0(
        ifelse(is.null(base.model), var, paste0("Relative ", var, " vs ", base.model)),
        ": Last ", use.n.years, " Years"
      )) +
      ylab(ifelse(is.null(base.model), var, paste0("Relative ", var, " Difference"))) +
      xlab("Model") +
      theme_bw()
    
  } else if (plot.style == "median_iqr") {
    
    # Compute summary statistics with 1.5x IQR whiskers
    res_summary <- res %>%
      group_by(Model, Label) %>%
      summarise(
        q1  = quantile(!!sym(var), 0.25),
        med = median(!!sym(var)),
        q3  = quantile(!!sym(var), 0.75),
        iqr = q3 - q1,
        .groups = "drop"
      ) %>%
      mutate(
        x    = as.numeric(factor(Model)),
        ymin = if (show.whisker) q1 - 1.5 * iqr else NA_real_,
        ymax = if (show.whisker) q3 + 1.5 * iqr else NA_real_
      )
    
    # Clip whiskers to observed range
    res_limits <- res %>%
      group_by(Model, Label) %>%
      summarise(
        min_val = min(!!sym(var), na.rm = TRUE),
        max_val = max(!!sym(var), na.rm = TRUE),
        .groups = "drop"
      )
    
    res_summary <- res_summary %>%
      left_join(res_limits, by = c("Model", "Label")) %>%
      mutate(
        ymin = pmax(ymin, min_val),
        ymax = pmin(ymax, max_val)
      )
    
    p <- ggplot(res_summary, aes(x = x, color = Model)) +
      { if (show.whisker) geom_segment(aes(x = x, xend = x, y = ymin, yend = q1)) } +
      { if (show.whisker) geom_segment(aes(x = x, xend = x, y = q3, yend = ymax)) } +
      
      # IQR box (no fill)
      geom_rect(
        aes(xmin = x - 0.3, xmax = x + 0.3, ymin = q1, ymax = q3, col = Model),
        fill = NA, linewidth = 0.8
      ) +
      
      # Median line
      geom_segment(
        aes(x = x - 0.3, xend = x + 0.3, y = med, yend = med, color = Model),
        linewidth = 0.8
      ) +
      
      scale_x_continuous(breaks = res_summary$x, labels = res_summary$Model) +
      facet_grid(Label ~ ., scales = "free") +
      scale_color_viridis_d(option = col.opt) +
      ggtitle(paste0(
        ifelse(is.null(base.model), var, paste0("Relative ", var, " vs ", base.model)),
        ": Last ", use.n.years, " Years"
      )) +
      ylab(ifelse(is.null(base.model), var, paste0("Relative ", var, " Difference"))) +
      xlab("Model") +
      theme_bw()
    
  } else {
    stop("Unknown plot.style. Choose 'boxplot' or 'median_iqr'.")
  }
  
  # 8) Save
  plot_name <- paste0(
    var,
    ifelse(is.null(base.model), "", "_Relative"),
    "_last_", use.n.years, "_years.PNG"
  )
  
  new_sub_dir <- file.path(main.dir, sub.dir, "Performance_Boxplot")
  if (!file.exists(new_sub_dir)) dir.create(new_sub_dir, recursive = TRUE)
  
  ggsave(
    file.path(new_sub_dir, plot_name),
    plot   = p,
    width  = width,
    height = height,
    dpi    = dpi
  )
  
  return(p)
}


plot_ssb_performance2 <- function(mods, is.nsim, main.dir, sub.dir, var = "SSB",
                                  width = 10, height = 7, dpi = 300, col.opt = "D",
                                  method = NULL,
                                  outlier.opt = NA,
                                  plot.style = "median_iqr",
                                  show.whisker = TRUE,
                                  use.n.years = NULL,
                                  start.years = NULL,
                                  new_model_names = NULL,
                                  base.model = NULL
                                  ) {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(rlang)
  
  if (is.null(use.n.years)) {
    cat("\nuse.n.years is not specified, so default (5 years) is used here!\n")
    use.n.years <- 5
  }
  
  if (is.null(start.years)) {
    cat("\nstart.years is not specified, so default (1st year in historical period) is used here!\n")
    start.years <- 1
  }
  
  res <- NULL
  
  Years <- if (!is.nsim) mods[[1]]$om$years else mods[[1]][[1]]$om$years
  
  if (!is.nsim) {
    res <- lapply(seq_along(mods), function(i) {
      tmp <- data.frame(
        SSB = mods[[i]]$om$rep$SSB,
        Model = paste0("Model", i),
        Year = Years,
        Realization = 1
      )
      start_idx <- start.years
      end_idx <- min(start.years + use.n.years - 1, nrow(tmp))
      tmp[start_idx:end_idx, ]
    }) %>% bind_rows()
  } else {
    res <- lapply(seq_along(mods), function(r) {
      lapply(seq_along(mods[[r]]), function(m) {
        tmp <- data.frame(
          SSB = mods[[r]][[m]]$om$rep$SSB,
          Model = paste0("Model", m),
          Year = Years,
          Realization = r
        )
        start_idx <- start.years
        end_idx <- min(start.years + use.n.years - 1, nrow(tmp))
        tmp[start_idx:end_idx, ]
      }) %>% bind_rows()
    }) %>% bind_rows()
  }
  
  # Add the last column: total SSB by Year and Realization
  res <- res %>%
    rowwise() %>%
    mutate(SSB_Global = sum(c_across(starts_with("SSB")), na.rm = TRUE)) %>%
    ungroup()
  
  # Rename Model names if needed
  if (!is.null(new_model_names)) {
    if (length(new_model_names) != length(unique(res$Model))) {
      stop("Length of new_model_names must match the number of models.")
    }
    res$Model <- factor(res$Model,
                        levels = paste0("Model", seq_along(new_model_names)),
                        labels = new_model_names)
    # if (!is.null(base.model)) base.model <- new_model_names[base.model]
    if (!is.null(base.model)) {
      if (!(base.model %in% new_model_names)) {
        warning("base.model does not match any of the new_model_names.")
      }
    }
  }
  
  # Pivot longer if multiple SSB types
  res <- pivot_longer(res, cols = starts_with(var), names_to = "Label", values_to = var)
  
  # Apply relative comparison if base.model is specified
  if (!is.null(base.model)) {
    base_df <- res %>% filter(Model == base.model) %>%
      rename(base_val = !!sym(var)) %>%
      select(Realization, Year, Label, base_val)
    
    res <- left_join(res, base_df, by = c("Realization", "Year", "Label")) %>%
      mutate(!!var := (!!sym(var)) / base_val - 1)
  }
  
  # Apply mean or median summarization if method is specified
  if (!is.null(method)) {
    res <- res %>%
      group_by(Model, Realization, Label) %>%
      summarise(!!var := if (method == "mean") mean(!!sym(var)) else median(!!sym(var)),
                .groups = "drop")
  }
  
  # Plot
  if (plot.style == "boxplot") {
    p <- ggplot(res, aes(x = Model, y = SSB, color = Model)) +
      geom_boxplot(lwd = 0.8, outlier.shape = outlier.opt) +
      facet_grid(Label ~ ., scales = "free") +
      scale_color_viridis_d(option = col.opt) +
      ggtitle(paste0(ifelse(is.null(base.model), var, paste0("Relative ", var, " vs ", base.model)),
                     ": Years ", start.years, " to ", start.years + use.n.years - 1)) +
      ylab(ifelse(is.null(base.model), var, paste0("Relative ", var, " Difference"))) +
      xlab("Model") + 
      theme_bw()
  } else if (plot.style == "median_iqr") {
    # Compute summary statistics with 1.5x IQR whiskers
    res_summary <- res %>%
      group_by(Model, Label) %>%
      summarise(
        q1 = quantile(!!sym(var), 0.25),
        med = median(!!sym(var)),
        q3 = quantile(!!sym(var), 0.75),
        iqr = q3 - q1,
        .groups = "drop"
      ) %>%
      mutate(
        x = as.numeric(factor(Model)),
        ymin = if (show.whisker) q1 - 1.5 * iqr else NA_real_,
        ymax = if (show.whisker) q3 + 1.5 * iqr else NA_real_
      )
    
    # Clip whiskers to the observed range
    res_limits <- res %>%
      group_by(Model, Label) %>%
      summarise(
        min_val = min(!!sym(var), na.rm = TRUE),
        max_val = max(!!sym(var), na.rm = TRUE),
        .groups = "drop"
      )
    
    res_summary <- left_join(res_summary, res_limits, by = c("Model", "Label")) %>%
      mutate(
        ymin = pmax(ymin, min_val),
        ymax = pmin(ymax, max_val)
      )
    
    p <- ggplot(res_summary, aes(x = x, color = Model)) +
      # Whiskers (1.5 x IQR, clipped to observed range)
      {if (show.whisker) geom_segment(aes(x = x, xend = x, y = ymin, yend = q1)) } +
      {if (show.whisker) geom_segment(aes(x = x, xend = x, y = q3, yend = ymax)) } +
      
      # IQR box (no fill)
      geom_rect(aes(xmin = x - 0.3, xmax = x + 0.3, ymin = q1, ymax = q3, col = Model),
                fill = NA, linewidth = 0.8) +
      
      # Median line
      geom_segment(aes(x = x - 0.3, xend = x + 0.3, y = med, yend = med, color = Model),
                   linewidth = 0.8) +
      
      scale_x_continuous(breaks = res_summary$x, labels = res_summary$Model) +
      facet_grid(Label ~ ., scales = "free") +
      scale_color_viridis_d(option = col.opt) +
      ggtitle(paste0(ifelse(is.null(base.model), var, paste0("Relative ", var, " vs ", base.model)),
                     ": Years ", start.years, " to ", start.years + use.n.years - 1)) +
      ylab(ifelse(is.null(base.model), var, paste0("Relative ", var, " Difference"))) +
      xlab("Model") + 
      theme_bw()
  } else {
    stop("Unknown plot.style. Choose 'boxplot' or 'median_iqr'.")
  }
  
  # Save the plot
  plot_name <- paste0(var, ifelse(is.null(base.model), "", "_Relative"),
                      "_first_", use.n.years, "_years.PNG")
  
  # Create the new subfolder if it doesn't exist
  new_sub_dir <- file.path(main.dir, sub.dir, "Performance_Boxplot")
  
  if (!file.exists(new_sub_dir)){
    dir.create(new_sub_dir)
  }
  
  # Save the figure inside the new subfolder
  ggsave(file.path(new_sub_dir, plot_name), plot = p, width = width, height = height, dpi = dpi)
  
  return(p)
}

plot_fbar_performance <- function(mods, is.nsim, main.dir, sub.dir, var = "Fbar",
                                  width = 10, height = 7, dpi = 300, col.opt = "D",
                                  method = NULL,
                                  outlier.opt = NA,
                                  plot.style = "median_iqr", 
                                  show.whisker = TRUE,
                                  f.ymin = NULL, 
                                  f.ymax = NULL,
                                  use.n.years = NULL,
                                  new_model_names = NULL,
                                  base.model = NULL
                                  ) {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(rlang)
  
  if (is.null(use.n.years)) {
    cat("\nuse.n.years is not specified, so default (5 years) is used here!\n")
    use.n.years <- 5
  }
  
  make_plot_data <- function(index_range, label_prefix) {
    if (!is.nsim) {
      Years <- mods[[1]]$om$years
      n_fleets <- mods[[1]]$om$input$data$n_fleets[1]
      n_regions <- mods[[1]]$om$input$data$n_regions[1]
      lapply(seq_along(mods), function(i) {
        tmp <- mods[[i]]$om$rep$Fbar[, index_range, drop = FALSE]
        tmp <- as.data.frame(tmp)
        names(tmp) <- paste0(label_prefix, seq_along(index_range))
        tmp$Model <- paste0("Model", i)
        tmp$Year <- Years
        tmp$Realization <- 1
        tmp <- tail(tmp, use.n.years)
      }) %>% bind_rows()
    } else {
      Years <- mods[[1]][[1]]$om$years
      lapply(seq_along(mods), function(r) {
        lapply(seq_along(mods[[r]]), function(m) {
          tmp <- mods[[r]][[m]]$om$rep$Fbar[, index_range, drop = FALSE]
          tmp <- as.data.frame(tmp)
          names(tmp) <- paste0(label_prefix, seq_along(index_range))
          tmp$Model <- paste0("Model", m)
          tmp$Year <- Years
          tmp$Realization <- r
          tmp <- tail(tmp, use.n.years)
        }) %>% bind_rows()
      }) %>% bind_rows()
    }
  }
  
  if (!is.nsim) {
    # Get number of fleets and regions
    n_fleets <- mods[[1]]$om$input$data$n_fleets[1]
    n_regions <- mods[[1]]$om$input$data$n_regions[1]
  } else {
    n_fleets <- mods[[1]][[1]]$om$input$data$n_fleets[1]
    n_regions <- mods[[1]][[1]]$om$input$data$n_regions[1]
  }
  
  # Fleet-level Fbar
  res_fleet <- make_plot_data(1:n_fleets, "Fleet_")
  
  # Region-level Fbar
  res_region <- make_plot_data((n_fleets + 1):(n_fleets + n_regions), "Region_")
  
  # Global Fbar (only one column)
  res_global <- make_plot_data(n_fleets + n_regions + 1, "Global")
  
  plot_and_save <- function(res, title, ylab_text, filename) {
    if (!is.null(new_model_names)) {
      if (length(new_model_names) != length(unique(res$Model))) {
        stop("Length of new_model_names must match the number of models.")
      }
      res$Model <- factor(res$Model,
                          levels = paste0("Model", seq_along(new_model_names)),
                          labels = new_model_names)
      if (!is.null(base.model)) {
        if (!(base.model %in% new_model_names)) {
          warning("base.model does not match any of the new_model_names.")
        }
      }
    }
    
    res_long <- pivot_longer(res, cols = starts_with(c("Fleet_", "Region_", "Global")),
                             names_to = "Label", values_to = "Fbar")
    
    # Apply relative comparison if requested
    if (!is.null(base.model)) {
      base_df <- res_long %>% filter(Model == base.model) %>%
        rename(base_val = Fbar) %>%
        select(Realization, Year, Label, base_val)
      
      res_long <- left_join(res_long, base_df, by = c("Realization", "Year", "Label")) %>%
        mutate(Fbar = Fbar / base_val - 1)
    }
    
    if (!is.null(base.model)) {
      if(!is.null(f.ymin)) y1 = f.ymin else y1 = -1
      if(!is.null(f.ymax)) y2 = f.ymax else y2 = 2
    } else {
      if(!is.null(f.ymin)) y1 = f.ymin else y1 = 0
      if(!is.null(f.ymax)) y2 = f.ymax else y2 = 2
    }
    
    # Apply mean or median summarization if method is specified
    if (!is.null(method)) {
      res_long <- res_long %>%
        group_by(Model, Realization, Label) %>%
        summarise(!!var := if (method == "mean") mean(!!sym(var)) else median(!!sym(var)),
                  .groups = "drop")
    }
    
    # Plot
    if (plot.style == "boxplot") {
      p <- ggplot(res_long, aes(x = Model, y = Fbar, color = Model)) +
        geom_boxplot(lwd = 0.8, outlier.shape = outlier.opt) +
        coord_cartesian(ylim = c(y1, y2)) + 
        facet_grid(Label ~ ., scales = "free") +
        scale_color_viridis_d(option = col.opt) +
        ggtitle(paste0(ifelse(is.null(base.model), title, paste0("Relative ", title, " vs ", base.model)),
                       ": Last ", use.n.years, " Years")) +
        ylab(ifelse(is.null(base.model), ylab_text, "Relative Difference in Fbar")) +
        xlab("Model") + 
        theme_bw()
    } else if (plot.style == "median_iqr") {
      # Compute summary statistics with 1.5x IQR whiskers
      res_summary <- res_long %>%
        group_by(Model, Label) %>%
        summarise(
          q1 = quantile(!!sym(var), 0.25),
          med = median(!!sym(var)),
          q3 = quantile(!!sym(var), 0.75),
          iqr = q3 - q1,
          .groups = "drop"
        ) %>%
        mutate(
          x = as.numeric(factor(Model)),
          ymin = if (show.whisker) q1 - 1.5 * iqr else NA_real_,
          ymax = if (show.whisker) q3 + 1.5 * iqr else NA_real_
        )
      
      # Clip whiskers to the observed range
      res_limits <- res_long %>%
        group_by(Model, Label) %>%
        summarise(
          min_val = min(!!sym(var), na.rm = TRUE),
          max_val = max(!!sym(var), na.rm = TRUE),
          .groups = "drop"
        )
      
      res_summary <- left_join(res_summary, res_limits, by = c("Model", "Label")) %>%
        mutate(
          ymin = pmax(ymin, min_val),
          ymax = pmin(ymax, max_val)
        )
      
      p <- ggplot(res_summary, aes(x = x, color = Model)) +
        # Whiskers (1.5 x IQR, clipped to observed range)
        {if (show.whisker) geom_segment(aes(x = x, xend = x, y = ymin, yend = q1)) } +
        {if (show.whisker) geom_segment(aes(x = x, xend = x, y = q3, yend = ymax)) } +
        
        # IQR box (no fill)
        geom_rect(aes(xmin = x - 0.3, xmax = x + 0.3, ymin = q1, ymax = q3, col = Model),
                  fill = NA, linewidth = 0.8) +
        
        # Median line
        geom_segment(aes(x = x - 0.3, xend = x + 0.3, y = med, yend = med, color = Model),
                     linewidth = 0.8) +
        
        scale_x_continuous(breaks = res_summary$x, labels = res_summary$Model) +
        facet_grid(Label ~ ., scales = "free") +
        scale_color_viridis_d(option = col.opt) +
        ggtitle(paste0(ifelse(is.null(base.model), title, paste0("Relative ", title, " vs ", base.model)),
                       ": Last ", use.n.years, " Years")) +
        ylab(ifelse(is.null(base.model), ylab_text, "Relative Difference in Fbar")) +
        xlab("Model") + 
        theme_bw()
    } else {
      stop("Unknown plot.style. Choose 'boxplot' or 'median_iqr'.")
    }
    
    plot_name <- paste0(filename, ifelse(is.null(base.model), "", "_Relative"), ".PNG")
    
    # Create the new subfolder if it doesn't exist
    new_sub_dir <- file.path(main.dir, sub.dir, "Performance_Boxplot")
    
    if (!file.exists(new_sub_dir)){
      dir.create(new_sub_dir)
    }
    
    # Save the figure inside the new subfolder
    ggsave(file.path(new_sub_dir, plot_name), plot = p, width = width, height = height, dpi = dpi)
    
    return(p)
  }
  
  p_fleet <- plot_and_save(res_fleet, "Fbar by Fleet", "Fbar", paste0(var, "_fleet_last_", use.n.years, "_years"))
  p_region <- plot_and_save(res_region, "Fbar by Region", "Fbar", paste0(var, "_region_last_", use.n.years, "_years"))
  p_global <- plot_and_save(res_global, "Global Fbar", "Fbar", paste0(var, "_global_last_", use.n.years, "_years"))
  
  return(list(fleet = p_fleet, region = p_region, global = p_global))
}


plot_fbar_performance2 <- function(mods, is.nsim, main.dir, sub.dir, var = "Fbar",
                                   width = 10, height = 7, dpi = 300, col.opt = "D",
                                   method = NULL,
                                   outlier.opt = NA,
                                   plot.style = "median_iqr", 
                                   show.whisker = TRUE,
                                   f.ymin = NULL, 
                                   f.ymax = NULL, 
                                   use.n.years = NULL,
                                   start.years = NULL,
                                   new_model_names = NULL,
                                   base.model = NULL
                                   ) {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(rlang)
  
  if (is.null(use.n.years)) {
    cat("\nuse.n.years is not specified, so default (5 years) is used here!\n")
    use.n.years <- 5
  }
  
  if (is.null(start.years)) {
    cat("\nstart.years is not specified, so default (1st year in historical period) is used here!\n")
    start.years <- 1
  }
  
  make_plot_data <- function(index_range, label_prefix) {
    
    if (!is.nsim) {
      Years <- mods[[1]]$om$years
      lapply(seq_along(mods), function(i) {
        tmp <- mods[[i]]$om$rep$Fbar[, index_range, drop = FALSE]
        tmp <- as.data.frame(tmp)
        names(tmp) <- paste0(label_prefix, seq_along(index_range))
        tmp$Model <- paste0("Model", i)
        tmp$Year <- Years
        tmp$Realization <- 1
        
        start_idx <- start.years
        end_idx <- min(start.years + use.n.years - 1, nrow(tmp))
        tmp[start_idx:end_idx, ]
      }) %>% bind_rows()
    } else {
      Years <- mods[[1]][[1]]$om$years
      lapply(seq_along(mods), function(r) {
        lapply(seq_along(mods[[r]]), function(m) {
          tmp <- mods[[r]][[m]]$om$rep$Fbar[, index_range, drop = FALSE]
          tmp <- as.data.frame(tmp)
          names(tmp) <- paste0(label_prefix, seq_along(index_range))
          tmp$Model <- paste0("Model", m)
          tmp$Year <- Years
          tmp$Realization <- r
          
          start_idx <- start.years
          end_idx <- min(start.years + use.n.years - 1, nrow(tmp))
          tmp[start_idx:end_idx, ]
        }) %>% bind_rows()
      }) %>% bind_rows()
    }
  }
  
  # Get fleet and region info
  if (!is.nsim) {
    n_fleets <- mods[[1]]$om$input$data$n_fleets[1]
    n_regions <- mods[[1]]$om$input$data$n_regions[1]
  } else {
    n_fleets <- mods[[1]][[1]]$om$input$data$n_fleets[1]
    n_regions <- mods[[1]][[1]]$om$input$data$n_regions[1]
  }
  
  res_fleet <- make_plot_data(1:n_fleets, "Fleet_")
  res_region <- make_plot_data((n_fleets + 1):(n_fleets + n_regions), "Region_")
  res_global <- make_plot_data(n_fleets + n_regions + 1, "Global")
  
  plot_and_save <- function(res, title, ylab_text, filename) {
    if (!is.null(new_model_names)) {
      if (length(new_model_names) != length(unique(res$Model))) {
        stop("Length of new_model_names must match the number of models.")
      }
      res$Model <- factor(res$Model,
                          levels = paste0("Model", seq_along(new_model_names)),
                          labels = new_model_names)
      # if (!is.null(base.model)) base.model <- new_model_names[base.model]
      if (!is.null(base.model)) {
        if (!(base.model %in% new_model_names)) {
          warning("base.model does not match any of the new_model_names.")
        }
      }
    }
    
    res_long <- pivot_longer(res, cols = starts_with(c("Fleet_", "Region_", "Global")),
                             names_to = "Label", values_to = "Fbar")
    
    # Apply relative comparison
    if (!is.null(base.model)) {
      base_df <- res_long %>%
        filter(Model == base.model) %>%
        rename(base_val = Fbar) %>%
        select(Realization, Year, Label, base_val)
      
      res_long <- left_join(res_long, base_df, by = c("Realization", "Year", "Label")) %>%
        mutate(Fbar = Fbar / base_val - 1)
    }
    
    if (!is.null(base.model)) {
      if(!is.null(f.ymin)) y1 = f.ymin else y1 = -1
      if(!is.null(f.ymax)) y2 = f.ymax else y2 = 2
    } else {
      if(!is.null(f.ymin)) y1 = f.ymin else y1 = 0
      if(!is.null(f.ymax)) y2 = f.ymax else y2 = 2
    }
    
    # Apply mean or median summarization if method is specified
    if (!is.null(method)) {
      res_long <- res_long %>%
        group_by(Model, Realization, Label) %>%
        summarise(!!var := if (method == "mean") mean(!!sym(var)) else median(!!sym(var)),
                  .groups = "drop")
    }
    
    # Plot
    if (plot.style == "boxplot") {
      p <- ggplot(res_long, aes(x = Model, y = Fbar, color = Model)) +
        geom_boxplot(lwd = 0.8, outlier.shape = outlier.opt) +
        coord_cartesian(ylim = c(y1, y2)) + 
        facet_grid(Label ~ ., scales = "free") +
        scale_color_viridis_d(option = col.opt) +
        ggtitle(paste0(ifelse(is.null(base.model), title, paste0("Relative ", title, " vs ", base.model)),
                       ": Years ", start.years, " to ", start.years + use.n.years - 1)) +
        ylab(ifelse(is.null(base.model), ylab_text, "Relative Difference in Fbar")) +
        xlab("Model") + 
        theme_bw()
    } else if (plot.style == "median_iqr") {
      # Compute summary statistics with 1.5x IQR whiskers
      res_summary <- res_long %>%
        group_by(Model, Label) %>%
        summarise(
          q1 = quantile(!!sym(var), 0.25),
          med = median(!!sym(var)),
          q3 = quantile(!!sym(var), 0.75),
          iqr = q3 - q1,
          .groups = "drop"
        ) %>%
        mutate(
          x = as.numeric(factor(Model)),
          ymin = if (show.whisker) q1 - 1.5 * iqr else NA_real_,
          ymax = if (show.whisker) q3 + 1.5 * iqr else NA_real_
        )
      
      # Clip whiskers to the observed range
      res_limits <- res_long %>%
        group_by(Model, Label) %>%
        summarise(
          min_val = min(!!sym(var), na.rm = TRUE),
          max_val = max(!!sym(var), na.rm = TRUE),
          .groups = "drop"
        )
      
      res_summary <- left_join(res_summary, res_limits, by = c("Model", "Label")) %>%
        mutate(
          ymin = pmax(ymin, min_val),
          ymax = pmin(ymax, max_val)
        )
      
      p <- ggplot(res_summary, aes(x = x, color = Model)) +
        # Whiskers (1.5 x IQR, clipped to observed range)
        {if (show.whisker) geom_segment(aes(x = x, xend = x, y = ymin, yend = q1)) } +
        {if (show.whisker) geom_segment(aes(x = x, xend = x, y = q3, yend = ymax)) } +
        
        # IQR box (no fill)
        geom_rect(aes(xmin = x - 0.3, xmax = x + 0.3, ymin = q1, ymax = q3, col = Model),
                  fill = NA, linewidth = 0.8) +
        
        # Median line
        geom_segment(aes(x = x - 0.3, xend = x + 0.3, y = med, yend = med, color = Model),
                     linewidth = 0.8) +
        
        scale_x_continuous(breaks = res_summary$x, labels = res_summary$Model) +
        facet_grid(Label ~ ., scales = "free") +
        scale_color_viridis_d(option = col.opt) +
        ggtitle(paste0(ifelse(is.null(base.model), title, paste0("Relative ", title, " vs ", base.model)),
                       ": Years ", start.years, " to ", start.years + use.n.years - 1)) +
        ylab(ifelse(is.null(base.model), ylab_text, "Relative Difference in Fbar")) +
        xlab("Model") + 
        theme_bw()
    } else {
      stop("Unknown plot.style. Choose 'boxplot' or 'median_iqr'.")
    }
    
    plot_name <- paste0(filename, ifelse(is.null(base.model), "", "_Relative"), ".PNG")
    
    # Create the new subfolder if it doesn't exist
    new_sub_dir <- file.path(main.dir, sub.dir, "Performance_Boxplot")
    
    if (!file.exists(new_sub_dir)){
      dir.create(new_sub_dir)
    }
    
    # Save the figure inside the new subfolder
    ggsave(file.path(new_sub_dir, plot_name), plot = p, width = width, height = height, dpi = dpi)
    
    return(p)
  }
  
  p_fleet <- plot_and_save(res_fleet, "Fbar by Fleet", "Fbar", paste0(var, "_fleet_first_", use.n.years, "_years"))
  p_region <- plot_and_save(res_region, "Fbar by Region", "Fbar", paste0(var, "_region_first_", use.n.years, "_years"))
  p_global <- plot_and_save(res_global, "Global Fbar", "Fbar", paste0(var, "_global_first_", use.n.years, "_years"))
  
  return(list(fleet = p_fleet, region = p_region, global = p_global))
}


plot_catch_performance <- function(mods, is.nsim, main.dir, sub.dir, var = "Catch",
                                   width = 10, height = 7, dpi = 300, col.opt = "D",
                                   method = NULL,
                                   outlier.opt = NA,
                                   plot.style = "median_iqr", 
                                   show.whisker = TRUE,
                                   use.n.years = NULL,
                                   new_model_names = NULL,
                                   base.model = NULL
                                   ) {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(rlang)
  
  if (is.null(use.n.years)) {
    cat("\nuse.n.years is not specified, so default (5 years) is used here!\n")
    use.n.years <- 5
  }
  
  res <- NULL
  
  # Prepare data
  if (!is.nsim) {
    Years <- mods[[1]]$om$years
    res <- lapply(seq_along(mods), function(i) {
      data.frame(
        Catch = mods[[i]]$om$rep$pred_catch,
        Model = paste0("Model", i),
        Year = Years,
        Realization = 1
      ) %>% tail(use.n.years)
    }) %>% bind_rows()
  } else {
    Years <- mods[[1]][[1]]$om$years
    res <- lapply(seq_along(mods), function(r) {
      lapply(seq_along(mods[[r]]), function(m) {
        data.frame(
          Catch = mods[[r]][[m]]$om$rep$pred_catch,
          Model = paste0("Model", m),
          Year = Years,
          Realization = r
        ) %>% tail(use.n.years)
      }) %>% bind_rows()
    }) %>% bind_rows()
  }
  
  # Add the last column: total Catch by Year and Realization
  res <- res %>%
    rowwise() %>%
    mutate(Catch_Global = sum(c_across(starts_with("Catch")), na.rm = TRUE)) %>%
    ungroup()
  
  # Rename models if specified
  if (!is.null(new_model_names)) {
    if (length(new_model_names) != length(unique(res$Model))) {
      stop("Length of new_model_names must match the number of models.")
    }
    res$Model <- factor(res$Model,
                        levels = paste0("Model", seq_along(new_model_names)),
                        labels = new_model_names)
    # if (!is.null(base.model)) base.model <- new_model_names[base.model]
    if (!is.null(base.model)) {
      if (!(base.model %in% new_model_names)) {
        warning("base.model does not match any of the new_model_names.")
      }
    }
  }
  
  # Pivot longer if needed (for multiple Catch columns)
  res <- pivot_longer(res, cols = starts_with(var), names_to = "Label", values_to = var)
  
  # Relative difference from base.model
  if (!is.null(base.model)) {
    base_df <- res %>% filter(Model == base.model) %>%
      rename(base_val = !!sym(var)) %>%
      select(Realization, Year, Label, base_val)
    
    res <- left_join(res, base_df, by = c("Realization", "Year", "Label")) %>%
      mutate(!!var := (!!sym(var)) / base_val - 1) 
  }
  
  # Apply mean or median summarization if method is specified
  if (!is.null(method)) {
    res <- res %>%
      group_by(Model, Realization, Label) %>%
      summarise(!!var := if (method == "mean") mean(!!sym(var)) else median(!!sym(var)),
                .groups = "drop")
  }
  
  # Plot
  if (plot.style == "boxplot") {
    p <- ggplot(res, aes(x = Model, y = Catch, color = Model)) +
      geom_boxplot(lwd = 0.8, outlier.shape = outlier.opt) +
      facet_grid(Label ~ ., scales = "free") +
      scale_color_viridis_d(option = col.opt) +
      ggtitle(paste0(ifelse(is.null(base.model), var, paste0("Relative ", var, " vs ", base.model)),
                     ": Last ", use.n.years, " Years")) +
      ylab(ifelse(is.null(base.model), var, paste0("Relative ", var, " Difference"))) +
      xlab("Model") + 
      theme_bw()
  } else if (plot.style == "median_iqr") {
    # Compute summary statistics with 1.5x IQR whiskers
    res_summary <- res %>%
      group_by(Model, Label) %>%
      summarise(
        q1 = quantile(!!sym(var), 0.25),
        med = median(!!sym(var)),
        q3 = quantile(!!sym(var), 0.75),
        iqr = q3 - q1,
        .groups = "drop"
      ) %>%
      mutate(
        x = as.numeric(factor(Model)),
        ymin = if (show.whisker) q1 - 1.5 * iqr else NA_real_,
        ymax = if (show.whisker) q3 + 1.5 * iqr else NA_real_
      )
    
    # Clip whiskers to the observed range
    res_limits <- res %>%
      group_by(Model, Label) %>%
      summarise(
        min_val = min(!!sym(var), na.rm = TRUE),
        max_val = max(!!sym(var), na.rm = TRUE),
        .groups = "drop"
      )
    
    res_summary <- left_join(res_summary, res_limits, by = c("Model", "Label")) %>%
      mutate(
        ymin = pmax(ymin, min_val),
        ymax = pmin(ymax, max_val)
      )
    
    p <- ggplot(res_summary, aes(x = x, color = Model)) +
      # Whiskers (1.5 x IQR, clipped to observed range)
      {if (show.whisker) geom_segment(aes(x = x, xend = x, y = ymin, yend = q1)) } +
      {if (show.whisker) geom_segment(aes(x = x, xend = x, y = q3, yend = ymax)) } +
      
      # IQR box (no fill)
      geom_rect(aes(xmin = x - 0.3, xmax = x + 0.3, ymin = q1, ymax = q3, col = Model),
                fill = NA, linewidth = 0.8) +
      
      # Median line
      geom_segment(aes(x = x - 0.3, xend = x + 0.3, y = med, yend = med, color = Model),
                   linewidth = 0.8) +
      
      scale_x_continuous(breaks = res_summary$x, labels = res_summary$Model) +
      facet_grid(Label ~ ., scales = "free") +
      scale_color_viridis_d(option = col.opt) +
      ggtitle(paste0(ifelse(is.null(base.model), var, paste0("Relative ", var, " vs ", base.model)),
                     ": Last ", use.n.years, " Years")) +
      ylab(ifelse(is.null(base.model), var, paste0("Relative ", var, " Difference"))) +
      xlab("Model") + 
      theme_bw()
  } else {
    stop("Unknown plot.style. Choose 'boxplot' or 'median_iqr'.")
  }
  
  plot_name <- paste0(var, ifelse(is.null(base.model), "", "_Relative"), "_last_", use.n.years, "_years.PNG")
  
  # Create the new subfolder if it doesn't exist
  new_sub_dir <- file.path(main.dir, sub.dir, "Performance_Boxplot")
  
  if (!file.exists(new_sub_dir)){
    dir.create(new_sub_dir)
  }
  
  # Save the figure inside the new subfolder
  ggsave(file.path(new_sub_dir, plot_name), plot = p, width = width, height = height, dpi = dpi)
  
  return(p)
}


plot_catch_performance2 <- function(mods, is.nsim, main.dir, sub.dir, var = "Catch",
                                    width = 10, height = 7, dpi = 300, col.opt = "D",
                                    method = NULL,
                                    outlier.opt = NA,
                                    plot.style = "median_iqr", 
                                    show.whisker = TRUE,
                                    use.n.years = NULL,
                                    start.years = NULL,
                                    new_model_names = NULL,
                                    base.model = NULL
                                    ) {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(rlang)
  
  if (is.null(use.n.years)) {
    cat("\nuse.n.years is not specified, so default (5 years) is used here!\n")
    use.n.years <- 5
  }
  
  if (is.null(start.years)) {
    cat("\nstart.years is not specified, so default (1st year in historical period) is used here!\n")
    start.years <- 1
  }
  
  res <- NULL
  
  if (!is.nsim) {
    Years <- mods[[1]]$om$years
    res <- lapply(seq_along(mods), function(i) {
      tmp <- data.frame(
        Catch = mods[[i]]$om$rep$pred_catch,
        Model = paste0("Model", i),
        Year = Years,
        Realization = 1
      )
      start_idx <- start.years
      end_idx <- min(start.years + use.n.years - 1, nrow(tmp))
      tmp[start_idx:end_idx, ]
    }) %>% bind_rows()
  } else {
    Years <- mods[[1]][[1]]$om$years
    res <- lapply(seq_along(mods), function(r) {
      lapply(seq_along(mods[[r]]), function(m) {
        tmp <- data.frame(
          Catch = mods[[r]][[m]]$om$rep$pred_catch,
          Model = paste0("Model", m),
          Year = Years,
          Realization = r
        )
        start_idx <- start.years
        end_idx <- min(start.years + use.n.years - 1, nrow(tmp))
        tmp[start_idx:end_idx, ]
      }) %>% bind_rows()
    }) %>% bind_rows()
  }
  
  # Add the last column: total Catch by Year and Realization
  res <- res %>%
    rowwise() %>%
    mutate(Catch_Global = sum(c_across(starts_with("Catch")), na.rm = TRUE)) %>%
    ungroup()
  
  # Rename models if requested
  if (!is.null(new_model_names)) {
    if (length(new_model_names) != length(unique(res$Model))) {
      stop("Length of new_model_names must match the number of models.")
    }
    res$Model <- factor(res$Model,
                        levels = paste0("Model", seq_along(new_model_names)),
                        labels = new_model_names)
    # if (!is.null(base.model)) base.model <- new_model_names[base.model]
    if (!is.null(base.model)) {
      if (!(base.model %in% new_model_names)) {
        warning("base.model does not match any of the new_model_names.")
      }
    }
  }
  
  # Pivot longer if multiple types
  res <- pivot_longer(res, cols = starts_with(var), names_to = "Label", values_to = var)
  
  # Calculate relative difference if base.model is specified
  if (!is.null(base.model)) {
    base_df <- res %>% filter(Model == base.model) %>%
      rename(base_val = !!sym(var)) %>%
      select(Realization, Year, Label, base_val)
    
    res <- left_join(res, base_df, by = c("Realization", "Year", "Label")) %>%
      mutate(!!var := (!!sym(var)) / base_val - 1)
  }
  
  # Apply mean or median summarization if method is specified
  if (!is.null(method)) {
    res <- res %>%
      group_by(Model, Realization, Label) %>%
      summarise(!!var := if (method == "mean") mean(!!sym(var)) else median(!!sym(var)),
                .groups = "drop")
  }
  
  # Plot
  if (plot.style == "boxplot") {
    p <- ggplot(res, aes(x = Model, y = Catch, color = Model)) +
      geom_boxplot(lwd = 0.8, outlier.shape = outlier.opt) +
      facet_grid(Label ~ ., scales = "free") +
      scale_color_viridis_d(option = col.opt) +
      ggtitle(paste0(ifelse(is.null(base.model), var, paste0("Relative ", var, " vs ", base.model)),
                     ": Years ", start.years, " to ", start.years + use.n.years - 1)) +
      ylab(ifelse(is.null(base.model), var, paste0("Relative ", var, " Difference"))) +
      xlab("Model") + 
      theme_bw()
  } else if (plot.style == "median_iqr") {
    # Compute summary statistics with 1.5x IQR whiskers
    res_summary <- res %>%
      group_by(Model, Label) %>%
      summarise(
        q1 = quantile(!!sym(var), 0.25),
        med = median(!!sym(var)),
        q3 = quantile(!!sym(var), 0.75),
        iqr = q3 - q1,
        .groups = "drop"
      ) %>%
      mutate(
        x = as.numeric(factor(Model)),
        ymin = if (show.whisker) q1 - 1.5 * iqr else NA_real_,
        ymax = if (show.whisker) q3 + 1.5 * iqr else NA_real_
      )
    
    # Clip whiskers to the observed range
    res_limits <- res %>%
      group_by(Model, Label) %>%
      summarise(
        min_val = min(!!sym(var), na.rm = TRUE),
        max_val = max(!!sym(var), na.rm = TRUE),
        .groups = "drop"
      )
    
    res_summary <- left_join(res_summary, res_limits, by = c("Model", "Label")) %>%
      mutate(
        ymin = pmax(ymin, min_val),
        ymax = pmin(ymax, max_val)
      )
    
    p <- ggplot(res_summary, aes(x = x, color = Model)) +
      # Whiskers (1.5 x IQR, clipped to observed range)
      {if (show.whisker) geom_segment(aes(x = x, xend = x, y = ymin, yend = q1)) } +
      {if (show.whisker) geom_segment(aes(x = x, xend = x, y = q3, yend = ymax)) } +
      
      # IQR box (no fill)
      geom_rect(aes(xmin = x - 0.3, xmax = x + 0.3, ymin = q1, ymax = q3, col = Model),
                fill = NA, linewidth = 0.8) +
      
      # Median line
      geom_segment(aes(x = x - 0.3, xend = x + 0.3, y = med, yend = med, color = Model),
                   linewidth = 0.8) +
      
      scale_x_continuous(breaks = res_summary$x, labels = res_summary$Model) +
      facet_grid(Label ~ ., scales = "free") +
      scale_color_viridis_d(option = col.opt) +
      ggtitle(paste0(ifelse(is.null(base.model), var, paste0("Relative ", var, " vs ", base.model)),
                     ": Years ", start.years, " to ", start.years + use.n.years - 1)) +
      ylab(ifelse(is.null(base.model), var, paste0("Relative ", var, " Difference"))) +
      xlab("Model") + 
      theme_bw()
  } else {
    stop("Unknown plot.style. Choose 'boxplot' or 'median_iqr'.")
  }
  
  # Save the plot
  plot_name <- paste0(var, ifelse(is.null(base.model), "", "_Relative"),"_first_", use.n.years, "_years.PNG")
  
  # Create the new subfolder if it doesn't exist
  new_sub_dir <- file.path(main.dir, sub.dir, "Performance_Boxplot")
  
  if (!file.exists(new_sub_dir)){
    dir.create(new_sub_dir)
  }
  
  # Save the figure inside the new subfolder
  ggsave(file.path(new_sub_dir, plot_name), plot = p, width = width, height = height, dpi = dpi)
  
  return(p)
}

plot_ssb_status <- function(mods, is.nsim, main.dir, sub.dir, var = "SSB_status",
                            width = 10, height = 7, dpi = 300, col.opt = "D",
                            method = NULL,
                            outlier.opt = NA,
                            plot.style = "median_iqr", 
                            show.whisker = TRUE,
                            use.n.years = NULL,
                            new_model_names = NULL,
                            base.model = NULL,
                            plot_prob = TRUE
                            ) {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(rlang)
  
  if (is.null(use.n.years)) {
    cat("\nuse.n.years is not specified, so default (5 years) is used here!\n")
    use.n.years <- 5
  }
  
  if (!is.nsim) {
    if (is.null(mods[[1]]$om$rep$log_SSB_FXSPR)) {
      message("Biological Reference Point has not been calculated internally!")
      return(invisible(NULL))
    }
  } else {
    if (is.null(mods[[1]][[1]]$om$rep$log_SSB_FXSPR)) {
      message("Biological Reference Point has not been calculated internally!")
      return(invisible(NULL))
    }
  }

  if (!is.nsim) {
    Years <- mods[[1]]$om$years
    title_main <- paste0("Probability SSB/SSB", mods[[1]]$om$input$data$percentSPR, "% < 0.5: Last ", use.n.years, " Years")
    res_list <- lapply(seq_along(mods), function(i) {
      tmp <- mods[[i]]$om$rep$SSB
      tmp <- cbind(tmp, rowSums(tmp))
      tmp <- tmp / exp(mods[[i]]$om$rep$log_SSB_FXSPR)
      tmp <- as.data.frame(tmp)
      name_tmp <- paste0("SSB/SSB", mods[[i]]$om$input$data$percentSPR, "%")
      names(tmp) <- paste0(name_tmp, ".s", 1:ncol(tmp))
      names(tmp)[ncol(tmp)] <- name_tmp
      tmp$Model <- paste0("Model", i)
      tmp$Year <- Years
      tmp$Realization <- 1
      tail(tmp, use.n.years)
    })
    res <- bind_rows(res_list)
    prob <- lapply(seq_along(mods), function(i) {
      x <- res_list[[i]]
      y <- x[, grepl("^SSB/SSB", names(x)), drop = FALSE]
      out <- as.data.frame(t(colMeans(y < 0.5, na.rm = TRUE)))
      out$Model <- paste0("Model", i)
      out
    }) %>% bind_rows()
  } else {
    Years <- mods[[1]][[1]]$om$years
    title_main <- paste0("Probability SSB/SSB", mods[[1]][[1]]$om$input$data$percentSPR, "% < 0.5: Last ", use.n.years, " Years")
    res_list <- lapply(seq_along(mods), function(r) {
      lapply(seq_along(mods[[r]]), function(m) {
        tmp <- mods[[r]][[m]]$om$rep$SSB
        tmp <- cbind(tmp, rowSums(tmp))
        tmp <- tmp / exp(mods[[r]][[m]]$om$rep$log_SSB_FXSPR)
        tmp <- as.data.frame(tmp)
        name_tmp <- paste0("SSB/SSB", mods[[r]][[m]]$om$input$data$percentSPR, "%")
        names(tmp) <- paste0(name_tmp, ".s", 1:ncol(tmp))
        names(tmp)[ncol(tmp)] <- name_tmp
        tmp$Model <- paste0("Model", m)
        tmp$Year <- Years
        tmp$Realization <- r
        tail(tmp, use.n.years)
      })
    })
    res <- bind_rows(res_list)
    prob <- lapply(seq_along(mods), function(r) {
      lapply(seq_along(mods[[r]]), function(m) {
        x <- res_list[[r]][[m]]
        y <- x[, grepl("^SSB/SSB", names(x)), drop = FALSE]
        out <- as.data.frame(t(colMeans(y < 0.5, na.rm = TRUE)))
        out$Model <- paste0("Model", m)
        out
      }) %>% bind_rows()
    }) %>% bind_rows()
  }
  
  # Rename models
  rename_models <- function(df) {
    if (!is.null(new_model_names)) {
      if (length(new_model_names) != length(unique(df$Model))) {
        stop("Length of new_model_names must match the number of models.")
      }
      df$Model <- factor(df$Model,
                         levels = paste0("Model", seq_along(new_model_names)),
                         labels = new_model_names)
      # if (!is.null(base.model)) base.model <<- new_model_names[base.model]
      if (!is.null(base.model)) {
        if (!(base.model %in% new_model_names)) {
          warning("base.model does not match any of the new_model_names.")
        }
      }
    }
    df
  }
  
  res <- rename_models(res)
  prob <- rename_models(prob)
  
  # Boxplot
  
  # Get the variable name pattern
  var_name <- unique(gsub("\\.s\\d+", "", grep("^SSB[./]SSB", names(res), value = TRUE)))
  
  # Pivot to 'value' column to avoid dynamic naming issues
  res_long <- pivot_longer(
    res,
    cols = matches("^SSB[./]SSB"),
    names_to = "Label",
    values_to = "value"
  )
  
  # Relative difference if base.model provided
  if (!is.null(base.model)) {
    base_df <- res_long %>%
      filter(Model == base.model) %>%
      rename(base_val = value) %>%
      select(Realization, Year, Label, base_val)
    
    res_long <- left_join(res_long, base_df, by = c("Realization", "Year", "Label")) %>%
      mutate(value = value / base_val - 1)
  }
  
  # Apply mean or median summarization if method is specified
  if (!is.null(method)) {
    var = "value"
    res_long <- res_long %>%
      group_by(Model, Realization, Label) %>%
      summarise(!!var := if (method == "mean") mean(!!sym(var)) else median(!!sym(var)),
                .groups = "drop")
  }
  
  # Plot
  if (plot.style == "boxplot") {
    p1 <- ggplot(res_long, aes(x = Model, y = value, color = Model)) +
      geom_boxplot(lwd = 0.8, outlier.shape = outlier.opt) +
      facet_grid(Label ~ ., scales = "free") +
      scale_color_viridis_d(option = col.opt) +
      ggtitle(ifelse(is.null(base.model),
                     paste0(var_name, ": Last ", use.n.years, " Years"),
                     paste0("Relative ", var_name, " vs ", base.model, ": Last ", use.n.years, " Years"))) +
      ylab(ifelse(is.null(base.model),
                  var_name,
                  paste0("Relative ", var_name, " Difference"))) +
      xlab("Model") + 
      theme_bw()
  } else if (plot.style == "median_iqr") {
    # Compute summary statistics with 1.5x IQR whiskers
    res_summary <- res_long %>%
      group_by(Model, Label) %>%
      summarise(
        q1 = quantile(value, 0.25, na.rm = TRUE),
        med = median(value, na.rm = TRUE),
        q3 = quantile(value, 0.75, na.rm = TRUE),
        iqr = q3 - q1,
        .groups = "drop"
      ) %>%
      mutate(
        x = as.numeric(factor(Model)),
        ymin = if (show.whisker) q1 - 1.5 * iqr else NA_real_,
        ymax = if (show.whisker) q3 + 1.5 * iqr else NA_real_
      )
    
    # Clip whiskers to observed range
    res_limits <- res_long %>%
      group_by(Model, Label) %>%
      summarise(
        min_val = min(value, na.rm = TRUE),
        max_val = max(value, na.rm = TRUE),
        .groups = "drop"
      )
    
    res_summary <- left_join(res_summary, res_limits, by = c("Model", "Label")) %>%
      mutate(
        ymin = pmax(ymin, min_val),
        ymax = pmin(ymax, max_val)
      )
    
    # Plot
    p1 <- ggplot(res_summary, aes(x = x, color = Model)) +
      {if (show.whisker) geom_segment(aes(x = x, xend = x, y = ymin, yend = q1))} +
      {if (show.whisker) geom_segment(aes(x = x, xend = x, y = q3, yend = ymax))} +
      # IQR box (no fill)
      geom_rect(aes(xmin = x - 0.3, xmax = x + 0.3, ymin = q1, ymax = q3, col = Model),
                fill = NA, linewidth = 0.8) +
      # Median line
      geom_segment(aes(x = x - 0.3, xend = x + 0.3, y = med, yend = med, color = Model),
                   linewidth = 0.8) +
      scale_x_continuous(breaks = res_summary$x, labels = res_summary$Model) +
      facet_grid(Label ~ ., scales = "free") +
      scale_color_viridis_d(option = col.opt) +
      ggtitle(ifelse(is.null(base.model),
                     paste0(var_name, ": Last ", use.n.years, " Years"),
                     paste0("Relative ", var_name, " vs ", base.model, ": Last ", use.n.years, " Years"))) +
      ylab(ifelse(is.null(base.model),
                  var_name,
                  paste0("Relative ", var_name, " Difference"))) +
      xlab("Model") + 
      theme_bw()
  } else {
    stop("Unknown plot.style. Choose 'boxplot' or 'median_iqr'.")
  }
  
  plot_name = paste0("SSB_status", ifelse(is.null(base.model), "", "_Relative"), "_last_", use.n.years, "_years.PNG")
  
  # Create the new subfolder if it doesn't exist
  new_sub_dir <- file.path(main.dir, sub.dir, "Status_Boxplot")
  
  if (!file.exists(new_sub_dir)){
    dir.create(new_sub_dir)
  }
  
  # Save plot
  ggsave(file.path(main.dir, sub.dir, "Status_Boxplot", plot_name), p1, width = width, height = height, dpi = dpi)
  
  # Probability plot (not affected by base.model)
  prob_long <- pivot_longer(prob, cols = matches("^SSB[./]SSB"), names_to = "Label", values_to = "Prob")
  
  if(plot_prob) {
    p2 <- ggplot(prob_long, aes(x = Model, y = Prob, color = Model)) +
      geom_boxplot(lwd = 0.8, outlier.shape = outlier.opt) +
      facet_grid(Label ~ ., scales = "free") +
      scale_color_viridis_d(option = col.opt) +
      ggtitle(title_main) +
      ylab("Probability") +
      theme_bw()
    
    plot_name = paste0("SSB_status_overfished_prob_last_", use.n.years, "_years.PNG")
    
    # Create the new subfolder if it doesn't exist
    new_sub_dir <- file.path(main.dir, sub.dir, "Status_Boxplot")
    
    if (!file.exists(new_sub_dir)){
      dir.create(new_sub_dir)
    }
    
    # Save plot
    ggsave(file.path(main.dir, sub.dir, "Status_Boxplot", plot_name), p2, width = width, height = height, dpi = dpi)
    
  } else {
    p2 <- NULL
  }
  
  return(list(
    boxplot = p1,
    probplot = p2
  ))
}

plot_ssb_status2 <- function(mods, is.nsim, main.dir, sub.dir, var = "SSB_status",
                             width = 10, height = 7, dpi = 300, col.opt = "D",
                             method = NULL,
                             outlier.opt = NA,
                             plot.style = "median_iqr", 
                             show.whisker = TRUE,
                             use.n.years = NULL,
                             start.years = NULL,
                             new_model_names = NULL,
                             base.model = NULL,
                             plot_prob = TRUE
                             ) {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(rlang)
  
  if (is.null(use.n.years)) {
    cat("\nuse.n.years is not specified, so default (5 years) is used here!\n")
    use.n.years <- 5
  }
  
  if (is.null(start.years)) {
    cat("\nstart.years is not specified, so default (1st year in historical period) is used here!\n")
    start.years <- 1
  }
  
  if (!is.nsim) {
    if (is.null(mods[[1]]$om$rep$log_SSB_FXSPR)) {
      message("Biological Reference Point has not been calculated internally!")
      return(invisible(NULL))
    }
  } else {
    if (is.null(mods[[1]][[1]]$om$rep$log_SSB_FXSPR)) {
      message("Biological Reference Point has not been calculated internally!")
      return(invisible(NULL))
    }
  }
  
  if (!is.nsim) {
    Years <- mods[[1]]$om$years
    title_main <- paste0("Probability SSB/SSB", mods[[1]]$om$input$data$percentSPR, "% < 0.5: Years ", start.years, " to ", start.years + use.n.years - 1) 
    res_list <- lapply(seq_along(mods), function(i) {
      tmp <- mods[[i]]$om$rep$SSB
      tmp <- cbind(tmp, rowSums(tmp))
      tmp <- tmp / exp(mods[[i]]$om$rep$log_SSB_FXSPR)
      tmp <- as.data.frame(tmp)
      name_tmp <- paste0("SSB/SSB", mods[[i]]$om$input$data$percentSPR, "%")
      names(tmp) <- paste0(name_tmp, ".s", 1:ncol(tmp))
      names(tmp)[ncol(tmp)] <- name_tmp
      tmp$Model <- paste0("Model", i)
      tmp$Year <- Years
      tmp$Realization <- 1
      tmp[start.years:min(start.years + use.n.years - 1, nrow(tmp)), ]
    })
    res <- bind_rows(res_list)
    prob <- lapply(seq_along(mods), function(i) {
      x <- res_list[[i]]
      y <- x[, grepl("^SSB/SSB", names(x)), drop = FALSE]
      out <- as.data.frame(t(colMeans(y < 0.5, na.rm = TRUE)))
      out$Model <- paste0("Model", i)
      out
    }) %>% bind_rows()
  } else {
    Years <- mods[[1]][[1]]$om$years
    title_main <- paste0("Probability SSB/SSB", mods[[1]][[1]]$om$input$data$percentSPR, "% < 0.5: Years ", start.years, " to ", start.years + use.n.years - 1)
    res_list <- lapply(seq_along(mods), function(r) {
      lapply(seq_along(mods[[r]]), function(m) {
        tmp <- mods[[r]][[m]]$om$rep$SSB
        tmp <- cbind(tmp, rowSums(tmp))
        tmp <- tmp / exp(mods[[r]][[m]]$om$rep$log_SSB_FXSPR)
        tmp <- as.data.frame(tmp)
        name_tmp <- paste0("SSB/SSB", mods[[r]][[m]]$om$input$data$percentSPR, "%")
        names(tmp) <- paste0(name_tmp, ".s", 1:ncol(tmp))
        names(tmp)[ncol(tmp)] <- name_tmp
        tmp$Model <- paste0("Model", m)
        tmp$Year <- Years
        tmp$Realization <- r
        tmp[start.years:min(start.years + use.n.years - 1, nrow(tmp)), ]
      })
    })
    res <- bind_rows(res_list)
    prob <- lapply(seq_along(mods), function(r) {
      lapply(seq_along(mods[[r]]), function(m) {
        x <- res_list[[r]][[m]]
        y <- x[, grepl("^SSB/SSB", names(x)), drop = FALSE]
        out <- as.data.frame(t(colMeans(y < 0.5, na.rm = TRUE)))
        out$Model <- paste0("Model", m)
        out
      }) %>% bind_rows()
    }) %>% bind_rows()
  }
  
  # --- Rename models if needed ---
  rename_models <- function(df) {
    if (!is.null(new_model_names)) {
      if (length(new_model_names) != length(unique(df$Model))) {
        stop("Length of new_model_names must match the number of models.")
      }
      df$Model <- factor(df$Model,
                         levels = paste0("Model", seq_along(new_model_names)),
                         labels = new_model_names)
      # if (!is.null(base.model)) base.model <<- new_model_names[base.model]
      if (!is.null(base.model)) {
        if (!(base.model %in% new_model_names)) {
          warning("base.model does not match any of the new_model_names.")
        }
      }
    }
    df
  }
  
  res <- rename_models(res)
  prob <- rename_models(prob)
  
  # Get the variable name pattern
  var_name <- unique(gsub("\\.s\\d+", "", grep("^SSB[./]SSB", names(res), value = TRUE)))
  
  # Pivot to 'value' column to avoid dynamic naming issues
  res_long <- pivot_longer(
    res,
    cols = matches("^SSB[./]SSB"),
    names_to = "Label",
    values_to = "value"
  )
  
  # Relative difference if base.model provided
  if (!is.null(base.model)) {
    base_df <- res_long %>%
      filter(Model == base.model) %>%
      rename(base_val = value) %>%
      select(Realization, Year, Label, base_val)
    
    res_long <- left_join(res_long, base_df, by = c("Realization", "Year", "Label")) %>%
      mutate(value = value / base_val - 1)
  }
  
  # Apply mean or median summarization if method is specified
  if (!is.null(method)) {
    var = "value"
    res_long <- res_long %>%
      group_by(Model, Realization, Label) %>%
      summarise(!!var := if (method == "mean") mean(value) else median(value),
                .groups = "drop")
  }
  
  # Plot
  if (plot.style == "boxplot") {
    p1 <- ggplot(res_long, aes(x = Model, y = value, color = Model)) +
      geom_boxplot(lwd = 0.8, outlier.shape = outlier.opt) +
      facet_grid(Label ~ ., scales = "free") +
      scale_color_viridis_d(option = col.opt) +
      ggtitle(paste0(ifelse(is.null(base.model), var_name, paste0("Relative ", var_name, " vs ", base.model)),
                     ": Years ", start.years, " to ", start.years + use.n.years - 1)) +
      ylab(ifelse(is.null(base.model),
                  var_name,
                  paste0("Relative ", var_name, " Difference"))) +
      xlab("Model") + 
      theme_bw()
  } else if (plot.style == "median_iqr") {
    # Compute summary statistics with 1.5x IQR whiskers
    res_summary <- res_long %>%
      group_by(Model, Label) %>%
      summarise(
        q1 = quantile(value, 0.25, na.rm = TRUE),
        med = median(value, na.rm = TRUE),
        q3 = quantile(value, 0.75, na.rm = TRUE),
        iqr = q3 - q1,
        .groups = "drop"
      ) %>%
      mutate(
        x = as.numeric(factor(Model)),
        ymin = if (show.whisker) q1 - 1.5 * iqr else NA_real_,
        ymax = if (show.whisker) q3 + 1.5 * iqr else NA_real_
      )
    
    # Clip whiskers to observed range
    res_limits <- res_long %>%
      group_by(Model, Label) %>%
      summarise(
        min_val = min(value, na.rm = TRUE),
        max_val = max(value, na.rm = TRUE),
        .groups = "drop"
      )
    
    res_summary <- left_join(res_summary, res_limits, by = c("Model", "Label")) %>%
      mutate(
        ymin = pmax(ymin, min_val),
        ymax = pmin(ymax, max_val)
      )
    
    # Plot
    p1 <- ggplot(res_summary, aes(x = x, color = Model)) +
      {if (show.whisker) geom_segment(aes(x = x, xend = x, y = ymin, yend = q1))} +
      {if (show.whisker) geom_segment(aes(x = x, xend = x, y = q3, yend = ymax))} +
      # IQR box (no fill)
      geom_rect(aes(xmin = x - 0.3, xmax = x + 0.3, ymin = q1, ymax = q3, col = Model),
                fill = NA, linewidth = 0.8) +
      # Median line
      geom_segment(aes(x = x - 0.3, xend = x + 0.3, y = med, yend = med, color = Model),
                   linewidth = 0.8) +
      scale_x_continuous(breaks = res_summary$x, labels = res_summary$Model) +
      facet_grid(Label ~ ., scales = "free") +
      scale_color_viridis_d(option = col.opt) +
      ggtitle(paste0(ifelse(is.null(base.model), var_name, paste0("Relative ", var_name, " vs ", base.model)),
                     ": Years ", start.years, " to ", start.years + use.n.years - 1)) +
      ylab(ifelse(is.null(base.model),
                  var_name,
                  paste0("Relative ", var_name, " Difference"))) +
      xlab("Model") + 
      theme_bw()
  } else {
    stop("Unknown plot.style. Choose 'boxplot' or 'median_iqr'.")
  }
  
  # Save plot
  plot_name = paste0("SSB_status", ifelse(is.null(base.model), "", "_Relative"),"_first_", use.n.years, "_years.PNG")
  
  # Create the new subfolder if it doesn't exist
  new_sub_dir <- file.path(main.dir, sub.dir, "Status_Boxplot")
  
  if (!file.exists(new_sub_dir)){
    dir.create(new_sub_dir)
  }
  
  # Save plot
  ggsave(file.path(main.dir, sub.dir, "Status_Boxplot", plot_name), p1, width = width, height = height, dpi = dpi)
  
  # Point plot of probability SSB/SSBSPR < 0.5 (unchanged)
  prob_long <- pivot_longer(prob, cols = matches("^SSB[./]SSB"), names_to = "Label", values_to = "Prob")
  
  if(plot_prob) {
    
    p2 <- ggplot(prob_long, aes(x = Model, y = Prob, color = Model)) +
      geom_boxplot(lwd = 0.8, outlier.shape = outlier.opt) +
      facet_grid(Label ~ ., scales = "free") +
      scale_color_viridis_d(option = col.opt) +
      ggtitle(title_main) +
      ylab("Probability") +
      theme_bw()
    
    # Save plot
    plot_name = paste0("SSB_status_overfished_prob_first_", use.n.years, "_years.PNG")
    
    # Create the new subfolder if it doesn't exist
    new_sub_dir <- file.path(main.dir, sub.dir, "Status_Boxplot")
    
    if (!file.exists(new_sub_dir)){
      dir.create(new_sub_dir)
    }
    
    # Save plot
    ggsave(file.path(main.dir, sub.dir, "Status_Boxplot", plot_name), p2, width = width, height = height, dpi = dpi)
    
  } else {
    p2 = NULL
  }
  
  return(list(
    boxplot = p1,
    probplot = p2
  ))
}

plot_fbar_status <- function(mods, is.nsim, main.dir, sub.dir, var = "Fbar_status",
                             width = 10, height = 7, dpi = 300, col.opt = "D",
                             method = NULL,
                             outlier.opt = NA,
                             plot.style = "median_iqr", 
                             show.whisker = TRUE,
                             f.ymin = NULL, 
                             f.ymax = NULL, 
                             use.n.years = NULL,
                             new_model_names = NULL,
                             base.model = NULL,
                             plot_prob = TRUE
                             ) {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  if (is.null(use.n.years)) {
    cat("\nuse.n.years is not specified, so default (5 years) is used here!\n")
    use.n.years <- 5
  }
  
  if (!is.nsim) {
    if (is.null(mods[[1]]$om$rep$log_SSB_FXSPR)) {
      message("Biological Reference Point has not been calculated internally!")
      return(invisible(NULL))
    }
  } else {
    if (is.null(mods[[1]][[1]]$om$rep$log_SSB_FXSPR)) {
      message("Biological Reference Point has not been calculated internally!")
      return(invisible(NULL))
    }
  }
  
  make_plot_data <- function(index_range, label_prefix) {
    if (!is.nsim) {
      Years = mods[[1]]$om$years
      res_list <- lapply(seq_along(mods), function(i) {
        tmp1 <- mods[[i]]$om$rep$Fbar[, index_range, drop = FALSE]
        tmp2 <- exp(mods[[i]]$om$rep$log_Fbar_XSPR[, index_range, drop = FALSE])
        tmp <- as.data.frame(tmp1 / tmp2)
        names(tmp) <- paste0(label_prefix, seq_along(index_range))
        tmp$Model <- paste0("Model", i)
        tmp$Year <- Years
        tmp$Realization <- 1
        tmp <- tail(tmp, use.n.years)
      })
      res <- bind_rows(res_list)
      prob <- lapply(res_list, function(x) {
        y <- x[, 1:length(index_range), drop = FALSE]
        data.frame(t(colMeans(y > 1, na.rm = TRUE)))
      }) %>% bind_rows() %>% mutate(Model = paste0("Model", seq_along(mods)))
      return(list(data = res, prob = prob))
    } else {
      Years = mods[[1]][[1]]$om$years
      res_list <- lapply(seq_along(mods), function(r) {
        lapply(seq_along(mods[[r]]), function(m) {
          tmp1 <- mods[[r]][[m]]$om$rep$Fbar[, index_range, drop = FALSE]
          tmp2 <- exp(mods[[r]][[m]]$om$rep$log_Fbar_XSPR[, index_range, drop = FALSE])
          tmp <- as.data.frame(tmp1 / tmp2)
          names(tmp) <- paste0(label_prefix, seq_along(index_range))
          tmp$Model <- paste0("Model", m)
          tmp$Year <- Years
          tmp$Realization <- r
          tmp <- tail(tmp, use.n.years)
        })
      })
      res <- bind_rows(res_list)
      prob <- lapply(seq_along(mods), function(r) {
        lapply(seq_along(mods[[r]]), function(m) {
          x <- res_list[[r]][[m]]
          y <- x[, 1:length(index_range), drop = FALSE]
          df <- data.frame(t(colMeans(y > 1, na.rm = TRUE)))
          df$Model <- paste0("Model", m)
          df
        }) %>% bind_rows()
      }) %>% bind_rows()
      return(list(data = res, prob = prob))
    }
  }
  
  if (!is.nsim) {
    title_main <- paste0("F/F", mods[[1]]$om$input$data$percentSPR, "%")
    n_fleets <- mods[[1]]$om$input$data$n_fleets[1]
    n_regions <- mods[[1]]$om$input$data$n_regions[1]
  } else {
    title_main <- paste0("F/F", mods[[1]][[1]]$om$input$data$percentSPR, "%")
    n_fleets <- mods[[1]][[1]]$om$input$data$n_fleets[1]
    n_regions <- mods[[1]][[1]]$om$input$data$n_regions[1]
  }
  
  res_fleet_all <- make_plot_data(1:n_fleets, "Fleet_")
  res_region_all <- make_plot_data((n_fleets + 1):(n_fleets + n_regions), "Region_")
  res_global_all <- make_plot_data(n_fleets + n_regions + 1, "Global")
  
  res_fleet <- res_fleet_all$data
  res_region <- res_region_all$data
  res_global <- res_global_all$data
  prob_fleet <- res_fleet_all$prob
  prob_region <- res_region_all$prob
  prob_global <- res_global_all$prob
  
  rename_models <- function(res_or_prob) {
    if (!is.null(new_model_names)) {
      if (length(new_model_names) != length(unique(res_or_prob$Model))) {
        stop("Length of new_model_names must match the number of models.")
      }
      res_or_prob$Model <- factor(res_or_prob$Model,
                                  levels = paste0("Model", seq_along(new_model_names)),
                                  labels = new_model_names)
    }
    return(res_or_prob)
  }
  
  res_fleet <- rename_models(res_fleet)
  res_region <- rename_models(res_region)
  res_global <- rename_models(res_global)
  prob_fleet <- rename_models(prob_fleet)
  prob_region <- rename_models(prob_region)
  prob_global <- rename_models(prob_global)
  
  if (!is.null(base.model) && !is.null(new_model_names)) {
    if (!(base.model %in% new_model_names)) {
      warning("base.model does not match any of the new_model_names.")
    }
  }
  
  plot_boxplot <- function(res, title, ylab_text, filename) {
    res_long <- pivot_longer(res, cols = starts_with(c("Fleet_", "Region_", "Global")),
                             names_to = "Label", values_to = "Fbar")
    if (!is.null(base.model)) {
      base_df <- res_long %>%
        filter(Model == base.model) %>%
        rename(base_val = Fbar) %>%
        select(Realization, Year, Label, base_val)
      res_long <- left_join(res_long, base_df, by = c("Realization", "Year", "Label")) %>%
        mutate(Fbar = Fbar / base_val - 1)
    }
    
    if (!is.null(base.model)) {
      if(!is.null(f.ymin)) y1 = f.ymin else y1 = -1
      if(!is.null(f.ymax)) y2 = f.ymax else y2 = 2
    } else {
      if(!is.null(f.ymin)) y1 = f.ymin else y1 = 0
      if(!is.null(f.ymax)) y2 = f.ymax else y2 = 2
    }
    
    # Apply mean or median summarization if method is specified
    if (!is.null(method)) {
      var = "Fbar"
      res_long <- res_long %>%
        group_by(Model, Realization, Label) %>%
        summarise(!!var := if (method == "mean") mean(!!sym(var)) else median(!!sym(var)),
                  .groups = "drop")
    }
    
    # Plot
    if (plot.style == "boxplot") {
      p1 <- ggplot(res_long, aes(x = Model, y = Fbar, color = Model)) +
        geom_boxplot(lwd = 0.8, outlier.shape = outlier.opt) +
        coord_cartesian(ylim = c(y1, y2)) + 
        facet_grid(Label ~ ., scales = "free") +
        scale_color_viridis_d(option = col.opt) +
        ggtitle(ifelse(is.null(base.model),
                       paste0(title, ": Last ", use.n.years, " Years"),
                       paste0("Relative ", title, " vs ", base.model, ": Last ", use.n.years, " Years"))) +
        ylab(ifelse(is.null(base.model), ylab_text, "Relative Difference")) +
        xlab("Model") + 
        theme_bw()
    } else if (plot.style == "median_iqr") {
      # Compute summary statistics with 1.5x IQR whiskers
      res_long$value <- res_long$Fbar 
      res_summary <- res_long %>%
        group_by(Model, Label) %>%
        summarise(
          q1 = quantile(value, 0.25, na.rm = TRUE),
          med = median(value, na.rm = TRUE),
          q3 = quantile(value, 0.75, na.rm = TRUE),
          iqr = q3 - q1,
          .groups = "drop"
        ) %>%
        mutate(
          x = as.numeric(factor(Model)),
          ymin = if (show.whisker) q1 - 1.5 * iqr else NA_real_,
          ymax = if (show.whisker) q3 + 1.5 * iqr else NA_real_
        )
      
      # Clip whiskers to observed range
      res_limits <- res_long %>%
        group_by(Model, Label) %>%
        summarise(
          min_val = min(value, na.rm = TRUE),
          max_val = max(value, na.rm = TRUE),
          .groups = "drop"
        )
      
      res_summary <- left_join(res_summary, res_limits, by = c("Model", "Label")) %>%
        mutate(
          ymin = pmax(ymin, min_val),
          ymax = pmin(ymax, max_val)
        )
      
      # Plot
      p1 <- ggplot(res_summary, aes(x = x, color = Model)) +
        {if (show.whisker) geom_segment(aes(x = x, xend = x, y = ymin, yend = q1))} +
        {if (show.whisker) geom_segment(aes(x = x, xend = x, y = q3, yend = ymax))} +
        # IQR box (no fill)
        geom_rect(aes(xmin = x - 0.3, xmax = x + 0.3, ymin = q1, ymax = q3, col = Model),
                  fill = NA, linewidth = 0.8) +
        # Median line
        geom_segment(aes(x = x - 0.3, xend = x + 0.3, y = med, yend = med, color = Model),
                     linewidth = 0.8) +
        scale_x_continuous(breaks = res_summary$x, labels = res_summary$Model) +
        facet_grid(Label ~ ., scales = "free") +
        scale_color_viridis_d(option = col.opt) +
        ggtitle(ifelse(is.null(base.model),
                       paste0(title, ": Last ", use.n.years, " Years"),
                       paste0("Relative ", title, " vs ", base.model, ": Last ", use.n.years, " Years"))) +
        ylab(ifelse(is.null(base.model), ylab_text, "Relative Difference")) +
        xlab("Model") + 
        theme_bw()
      } else {
        stop("Unknown plot.style. Choose 'boxplot' or 'median_iqr'.")
      }
    
    # Create the new subfolder if it doesn't exist
    new_sub_dir <- file.path(main.dir, sub.dir, "Status_Boxplot")
    
    if (!file.exists(new_sub_dir)){
      dir.create(new_sub_dir)
    }
    
    ggsave(file.path(main.dir, sub.dir, "Status_Boxplot", paste0(filename, ifelse(is.null(base.model), "", "_Relative"), ".PNG")),
           p1, width = width, height = height, dpi = dpi)
    
    return(p1)
  }
  
  plot_pointplot <- function(prob, title, ylab_text, filename) {
    prob_long <- pivot_longer(prob, cols = starts_with(c("Fleet_", "Region_", "Global")),
                              names_to = "Label", values_to = "Prob")
    p2 <- ggplot(prob_long, aes(x = Model, y = Prob, color = Model)) +
      geom_boxplot(lwd = 0.8, outlier.shape = outlier.opt) +
      facet_grid(Label ~ ., scales = "free") +
      scale_color_viridis_d(option = col.opt) +
      ggtitle(title) +
      ylab(ylab_text) +
      xlab("Model") + 
      theme_bw()
    
    # Create the new subfolder if it doesn't exist
    new_sub_dir <- file.path(main.dir, sub.dir, "Status_Boxplot")
    
    if (!file.exists(new_sub_dir)){
      dir.create(new_sub_dir)
    }
    
    ggsave(file.path(main.dir, sub.dir, "Status_Boxplot", paste0(filename, ".PNG")), p2, width = width, height = height, dpi = dpi)
    
    return(p2)
  }
  
  p_fleet_box <- plot_boxplot(res_fleet, paste(title_main, "by Fleet"), title_main, paste0(var, "_fleet_last_", use.n.years, "_years"))
  p_region_box <- plot_boxplot(res_region, paste(title_main, "by Region"), title_main, paste0(var, "_region_last_", use.n.years, "_years"))
  p_global_box <- plot_boxplot(res_global, paste(title_main, "Global"), title_main, paste0(var, "_global_last_", use.n.years, "_years"))
  
  if (plot_prob) {
    p_fleet_point <- plot_pointplot(prob_fleet, paste0("Probability (", title_main, " > 1) - Fleet: Last ", use.n.years, " Years"), "Probability", paste0(var, "_fleet_overfishing_prob_fleet_last_", use.n.years, "_years"))
    p_region_point <- plot_pointplot(prob_region, paste0("Probability (", title_main, " > 1) - Region: Last ", use.n.years, " Years"), "Probability", paste0(var, "_region_overfishing_prob_fleet_last_", use.n.years, "_years"))
    p_global_point <- plot_pointplot(prob_global, paste0("Probability (", title_main, " > 1) - Global: Last ", use.n.years, " Years"), "Probability", paste0(var, "_global_overfishing_prob_fleet_last_", use.n.years, "_years"))
  } else {
    p_fleet_point <- NULL
    p_region_point <- NULL
    p_global_point <- NULL
  }

  return(list(
    fleet_box = p_fleet_box,
    region_box = p_region_box,
    global_box = p_global_box,
    fleet_prob = p_fleet_point,
    region_prob = p_region_point,
    global_prob = p_global_point
  ))
}

plot_fbar_status2 <- function(mods, is.nsim, main.dir, sub.dir, var = "Fbar_status",
                              width = 10, height = 7, dpi = 300, col.opt = "D",
                              method = NULL,
                              outlier.opt = NA,
                              plot.style = "median_iqr", 
                              show.whisker = TRUE,
                              f.ymin = NULL, 
                              f.ymax = NULL, 
                              use.n.years = NULL,
                              start.years = NULL,
                              new_model_names = NULL,
                              base.model = NULL,
                              plot_prob = TRUE
                              ) {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  
  if (is.null(use.n.years)) {
    cat("\nuse.n.years is not specified, so default (5 years) is used here!\n")
    use.n.years <- 5
  }
  
  if (is.null(start.years)) {
    cat("\nstart.years is not specified, so default (1st year in historical period) is used here!\n")
    start.years <- 1
  }
  
  if (!is.nsim) {
    if (is.null(mods[[1]]$om$rep$log_SSB_FXSPR)) {
      message("Biological Reference Point has not been calculated internally!")
      return(invisible(NULL))
    }
  } else {
    if (is.null(mods[[1]][[1]]$om$rep$log_SSB_FXSPR)) {
      message("Biological Reference Point has not been calculated internally!")
      return(invisible(NULL))
    }
  }
  
  make_plot_data <- function(index_range, label_prefix) {
    if (!is.nsim) {
      Years = mods[[1]]$om$years
      res_list <- lapply(seq_along(mods), function(i) {
        tmp1 <- mods[[i]]$om$rep$Fbar[, index_range, drop = FALSE]
        tmp2 <- exp(mods[[i]]$om$rep$log_Fbar_XSPR[, index_range, drop = FALSE])
        tmp <- as.data.frame(tmp1 / tmp2)
        names(tmp) <- paste0(label_prefix, seq_along(index_range))
        tmp$Model <- paste0("Model", i)
        tmp$Year <- Years
        tmp$Realization <- 1
        start_idx <- start.years
        end_idx <- min(start.years + use.n.years - 1, nrow(tmp))
        tmp[start_idx:end_idx, ]
      })
      res <- bind_rows(res_list)
      prob <- lapply(res_list, function(x) {
        y <- x[, 1:length(index_range), drop = FALSE]
        data.frame(t(colMeans(y > 1, na.rm = TRUE)))
      }) %>% bind_rows() %>% mutate(Model = paste0("Model", seq_along(mods)))
      return(list(data = res, prob = prob))
    } else {
      Years = mods[[1]][[1]]$om$years
      res_list <- lapply(seq_along(mods), function(r) {
        lapply(seq_along(mods[[r]]), function(m) {
          tmp1 <- mods[[r]][[m]]$om$rep$Fbar[, index_range, drop = FALSE]
          tmp2 <- exp(mods[[r]][[m]]$om$rep$log_Fbar_XSPR[, index_range, drop = FALSE])
          tmp <- as.data.frame(tmp1 / tmp2)
          names(tmp) <- paste0(label_prefix, seq_along(index_range))
          tmp$Model <- paste0("Model", m)
          tmp$Year <- Years
          tmp$Realization <- r
          start_idx <- start.years
          end_idx <- min(start.years + use.n.years - 1, nrow(tmp))
          tmp[start_idx:end_idx, ]
        })
      })
      res <- bind_rows(res_list)
      prob <- lapply(seq_along(mods), function(r) {
        lapply(seq_along(mods[[r]]), function(m) {
          x <- res_list[[r]][[m]]
          y <- x[, 1:length(index_range), drop = FALSE]
          df <- data.frame(t(colMeans(y > 1, na.rm = TRUE)))
          df$Model <- paste0("Model", m)
          df
        }) %>% bind_rows()
      }) %>% bind_rows()
      return(list(data = res, prob = prob))
    }
  }
  
  if (!is.nsim) {
    n_fleets <- mods[[1]]$om$input$data$n_fleets[1]
    n_regions <- mods[[1]]$om$input$data$n_regions[1]
  } else {
    n_fleets <- mods[[1]][[1]]$om$input$data$n_fleets[1]
    n_regions <- mods[[1]][[1]]$om$input$data$n_regions[1]
  }
  
  res_fleet_all <- make_plot_data(1:n_fleets, "Fleet_")
  res_region_all <- make_plot_data((n_fleets + 1):(n_fleets + n_regions), "Region_")
  res_global_all <- make_plot_data(n_fleets + n_regions + 1, "Global")
  res_fleet <- res_fleet_all$data
  res_region <- res_region_all$data
  res_global <- res_global_all$data
  prob_fleet <- res_fleet_all$prob
  prob_region <- res_region_all$prob
  prob_global <- res_global_all$prob
  
  rename_models <- function(res_or_prob) {
    if (!is.null(new_model_names)) {
      if (length(new_model_names) != length(unique(res_or_prob$Model))) {
        stop("Length of new_model_names must match the number of models.")
      }
      res_or_prob$Model <- factor(res_or_prob$Model,
                                  levels = paste0("Model", seq_along(new_model_names)),
                                  labels = new_model_names)
    }
    return(res_or_prob)
  }
  
  res_fleet <- rename_models(res_fleet)
  res_region <- rename_models(res_region)
  res_global <- rename_models(res_global)
  prob_fleet <- rename_models(prob_fleet)
  prob_region <- rename_models(prob_region)
  prob_global <- rename_models(prob_global)
  
  if (!is.null(base.model) && !is.null(new_model_names)) {
    if (!(base.model %in% new_model_names)) {
      warning("base.model does not match any of the new_model_names.")
    }
  }
  
  plot_boxplot <- function(res, title, ylab_text, filename) {
    res_long <- pivot_longer(res, cols = starts_with(c("Fleet_", "Region_", "Global")),
                             names_to = "Label", values_to = "Fbar")
    if (!is.null(base.model)) {
      base_df <- res_long %>%
        filter(Model == base.model) %>%
        rename(base_val = Fbar) %>%
        select(Realization, Year, Label, base_val)
      res_long <- left_join(res_long, base_df, by = c("Realization", "Year", "Label")) %>%
        mutate(Fbar = Fbar / base_val - 1)
    }
    
    if (!is.null(base.model)) {
      if(!is.null(f.ymin)) y1 = f.ymin else y1 = -1
      if(!is.null(f.ymax)) y2 = f.ymax else y2 = 2
    } else {
      if(!is.null(f.ymin)) y1 = f.ymin else y1 = 0
      if(!is.null(f.ymax)) y2 = f.ymax else y2 = 2
    }
    
    # Apply mean or median summarization if method is specified
    if (!is.null(method)) {
      var = "Fbar"
      res_long <- res_long %>%
        group_by(Model, Realization, Label) %>%
        summarise(!!var := if (method == "mean") mean(!!sym(var)) else median(!!sym(var)),
                  .groups = "drop")
    }
    
    # Plot
    if (plot.style == "boxplot") {
      p1 <- ggplot(res_long, aes(x = Model, y = Fbar, color = Model)) +
        geom_boxplot(lwd = 0.8, outlier.shape = outlier.opt) +
        coord_cartesian(ylim = c(y1, y2)) + 
        facet_grid(Label ~ ., scales = "free") +
        scale_color_viridis_d(option = col.opt) +
        ggtitle(paste0(ifelse(is.null(base.model), title, paste0("Relative ", title, " vs ", base.model)),
                       ": Years ", start.years, " to ", start.years + use.n.years - 1)) +
        ylab(ifelse(is.null(base.model), ylab_text, "Relative Difference")) +
        theme_bw()
    } else if (plot.style == "median_iqr") {
      # Compute summary statistics with 1.5x IQR whiskers
      res_long$value <- res_long$Fbar 
      res_summary <- res_long %>%
        group_by(Model, Label) %>%
        summarise(
          q1 = quantile(value, 0.25, na.rm = TRUE),
          med = median(value, na.rm = TRUE),
          q3 = quantile(value, 0.75, na.rm = TRUE),
          iqr = q3 - q1,
          .groups = "drop"
        ) %>%
        mutate(
          x = as.numeric(factor(Model)),
          ymin = if (show.whisker) q1 - 1.5 * iqr else NA_real_,
          ymax = if (show.whisker) q3 + 1.5 * iqr else NA_real_
        )
      
      # Clip whiskers to observed range
      res_limits <- res_long %>%
        group_by(Model, Label) %>%
        summarise(
          min_val = min(value, na.rm = TRUE),
          max_val = max(value, na.rm = TRUE),
          .groups = "drop"
        )
      
      res_summary <- left_join(res_summary, res_limits, by = c("Model", "Label")) %>%
        mutate(
          ymin = pmax(ymin, min_val),
          ymax = pmin(ymax, max_val)
        )
      
      # Plot
      p1 <- ggplot(res_summary, aes(x = x, color = Model)) +
        {if (show.whisker) geom_segment(aes(x = x, xend = x, y = ymin, yend = q1))} +
        {if (show.whisker) geom_segment(aes(x = x, xend = x, y = q3, yend = ymax))} +
        # IQR box (no fill)
        geom_rect(aes(xmin = x - 0.3, xmax = x + 0.3, ymin = q1, ymax = q3, col = Model),
                  fill = NA, linewidth = 0.8) +
        # Median line
        geom_segment(aes(x = x - 0.3, xend = x + 0.3, y = med, yend = med, color = Model),
                     linewidth = 0.8) +
        scale_x_continuous(breaks = res_summary$x, labels = res_summary$Model) +
        facet_grid(Label ~ ., scales = "free") +
        scale_color_viridis_d(option = col.opt) +
        ggtitle(paste0(ifelse(is.null(base.model), title, paste0("Relative ", title, " vs ", base.model)),
                       ": Years ", start.years, " to ", start.years + use.n.years - 1)) +
        ylab(ifelse(is.null(base.model), ylab_text, "Relative Difference")) +
        xlab("Model") + 
        theme_bw()
    } else {
      stop("Unknown plot.style. Choose 'boxplot' or 'median_iqr'.")
    }
    
    # Create the new subfolder if it doesn't exist
    new_sub_dir <- file.path(main.dir, sub.dir, "Status_Boxplot")
    
    if (!file.exists(new_sub_dir)){
      dir.create(new_sub_dir)
    }
    
    ggsave(file.path(main.dir, sub.dir, "Status_Boxplot", paste0(filename, ifelse(is.null(base.model), "", "_Relative"), ".PNG")),
           p1, width = width, height = height, dpi = dpi)
    return(p1)
  }
  
  plot_pointplot <- function(prob, title, ylab_text, filename) {
    prob_long <- pivot_longer(prob, cols = starts_with(c("Fleet_", "Region_", "Global")),
                              names_to = "Label", values_to = "Prob")
    p2 <- ggplot(prob_long, aes(x = Model, y = Prob, color = Model)) +
      geom_boxplot(lwd = 0.8, outlier.shape = outlier.opt) +
      facet_grid(Label ~ ., scales = "free") +
      scale_color_viridis_d(option = col.opt) +
      ggtitle(title) +
      ylab(ylab_text) +
      xlab("Model") + 
      theme_bw()
    
    # Create the new subfolder if it doesn't exist
    new_sub_dir <- file.path(main.dir, sub.dir, "Status_Boxplot")
    
    if (!file.exists(new_sub_dir)){
      dir.create(new_sub_dir)
    }
    
    ggsave(file.path(main.dir, sub.dir, "Status_Boxplot", paste0(filename, ".PNG")), p2, width = width, height = height, dpi = dpi)
    return(p2)
  }
  
  if (!is.nsim) {
    title_main <- paste0("F/F", mods[[1]]$om$input$data$percentSPR, "%")
  } else {
    title_main <- paste0("F/F", mods[[1]][[1]]$om$input$data$percentSPR, "%")
  }
    
  p_fleet_box <- plot_boxplot(res_fleet, paste(title_main, "by Fleet"), title_main, paste0(var, "_fleet_first_", use.n.years, "_years"))
  p_region_box <- plot_boxplot(res_region, paste(title_main, "by Region"), title_main, paste0(var, "_region_first_", use.n.years, "_years"))
  p_global_box <- plot_boxplot(res_global, paste(title_main, "Global"), title_main, paste0(var, "_global_first_", use.n.years, "_years"))
  
  if (plot_prob) {
    p_fleet_point <- plot_pointplot(prob_fleet, paste0("Probability (", title_main, " > 1) - Fleet: Years ", start.years, " to ", start.years + use.n.years - 1), "Probability", paste0(var, "_fleet_overfishing_prob_fleet_first_", use.n.years, "_years"))
    p_region_point <- plot_pointplot(prob_region, paste0("Probability (", title_main, " > 1) - Region: Years ", start.years, " to ", start.years + use.n.years - 1), "Probability", paste0(var, "_region_overfishing_prob_fleet_first_", use.n.years, "_years"))
    p_global_point <- plot_pointplot(prob_global, paste0("Probability (", title_main, " > 1) - Global: Years ", start.years, " to ", start.years + use.n.years - 1), "Probability", paste0(var, "_global_overfishing_prob_fleet_first_", use.n.years, "_years"))
  } else {
    p_fleet_point <- NULL
    p_region_point <- NULL
    p_global_point <- NULL
  }

  return(list(
    fleet_box = p_fleet_box,
    region_box = p_region_box,
    global_box = p_global_box,
    fleet_prob = p_fleet_point,
    region_prob = p_region_point,
    global_prob = p_global_point
  ))
}

# plot_kobe_status <- function(mods, is.nsim, main.dir, sub.dir, 
#                              width = 10, height = 7, dpi = 300, col.opt = "D",
#                              new_model_names = NULL,
#                              use.n.years = NULL,
#                              show_density = FALSE) {
#   library(dplyr)
#   library(tidyr)
#   library(ggplot2)
#   library(viridis)
#   
#   if (is.null(use.n.years)) {
#     cat("\nuse.n.years is not specified, so default (terminal year) is used here!\n")
#     use.n.years <- 1
#   }
#   
#   # === Extract SSB/SSBxx% ===
#   if (!is.nsim) {
#     Years <- mods[[1]]$om$years
#     ssb_list <- lapply(seq_along(mods), function(i) {
#       tmp <- mods[[i]]$om$rep$SSB
#       tmp <- cbind(tmp, rowSums(tmp))
#       tmp <- tmp / exp(mods[[i]]$om$rep$log_SSB_FXSPR)
#       ssb <- tail(tmp[, ncol(tmp)], use.n.years)
#       data.frame(Model = paste0("Model", i),
#                  Realization = 1,
#                  Year = tail(Years, use.n.years),
#                  Overfished = ssb)
#     })
#     ssb_df <- bind_rows(ssb_list)
#   } else {
#     Years <- mods[[1]][[1]]$om$years
#     ssb_list <- lapply(seq_along(mods), function(r) {
#       lapply(seq_along(mods[[r]]), function(m) {
#         tmp <- mods[[r]][[m]]$om$rep$SSB
#         tmp <- cbind(tmp, rowSums(tmp))
#         tmp <- tmp / exp(mods[[r]][[m]]$om$rep$log_SSB_FXSPR)
#         ssb <- tail(tmp[, ncol(tmp)], use.n.years)
#         data.frame(Model = paste0("Model", m),
#                    Realization = r,
#                    Year = tail(Years, use.n.years),
#                    Overfished = ssb)
#       })
#     })
#     ssb_df <- bind_rows(ssb_list)
#   }
#   
#   # === Extract F/Fxx% ===
#   if (!is.nsim) {
#     fbar_list <- lapply(seq_along(mods), function(i) {
#       fbar <- mods[[i]]$om$rep$Fbar[, ncol(mods[[i]]$om$rep$Fbar)]
#       fbar_ref <- exp(mods[[i]]$om$rep$log_Fbar_XSPR[, ncol(mods[[i]]$om$rep$log_Fbar_XSPR)])
#       ff <- tail(fbar / fbar_ref, use.n.years)
#       data.frame(Model = paste0("Model", i),
#                  Realization = 1,
#                  Year = tail(Years, use.n.years),
#                  Overfishing = ff)
#     })
#     fbar_df <- bind_rows(fbar_list)
#   } else {
#     fbar_list <- lapply(seq_along(mods), function(r) {
#       lapply(seq_along(mods[[r]]), function(m) {
#         fbar <- mods[[r]][[m]]$om$rep$Fbar[, ncol(mods[[r]][[m]]$om$rep$Fbar)]
#         fbar_ref <- exp(mods[[r]][[m]]$om$rep$log_Fbar_XSPR[, ncol(mods[[r]][[m]]$om$rep$log_Fbar_XSPR)])
#         ff <- tail(fbar / fbar_ref, use.n.years)
#         data.frame(Model = paste0("Model", m),
#                    Realization = r,
#                    Year = tail(Years, use.n.years),
#                    Overfishing = ff)
#       })
#     })
#     fbar_df <- bind_rows(fbar_list)
#   }
#   
#   # === Merge ===
#   temp <- left_join(ssb_df, fbar_df, by = c("Model", "Realization", "Year"))
#   temp$Index <- paste0("Year ", temp$Year)
#   
#   # === Rename Models if needed ===
#   if (!is.null(new_model_names)) {
#     if (length(new_model_names) != length(unique(temp$Model))) {
#       stop("Length of new_model_names must match number of models.")
#     }
#     temp$Model <- factor(temp$Model,
#                          levels = paste0("Model", seq_along(new_model_names)),
#                          labels = new_model_names)
#   }
#   
#   percentSPR <- if (!is.nsim) mods[[1]]$om$input$data$percentSPR else mods[[1]][[1]]$om$input$data$percentSPR
#   
#   # === Plot
#   p <- ggplot(temp, aes(x = Overfished, y = Overfishing)) +
#     facet_wrap(~ Model) +
#     annotate('rect', xmin = 0.5, xmax = Inf, ymin = -Inf, ymax = 1, alpha = 0.2, fill = "yellow")
#   
#   # === Add smoothed 2D density layer
#   if (show_density) {
#     p <- p + geom_density_2d_filled(aes(fill = after_stat(level)), alpha = 0.4, contour_var = "density", bins = 100) +
#       scale_fill_viridis_d(option = col.opt) + guides(fill = "none") 
#   }
#   
#   # === Add points and reference lines
#   if (show_density) {
#     p <- p +
#       geom_point(aes(color = Model), size = 1.5, alpha = 0.3) +
#       scale_color_viridis_d(option = col.opt) +
#       geom_vline(xintercept = 0.5, linetype = "dashed", color = "red", size = 1) +
#       geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 1) +
#       xlab(bquote(paste("SSB/", SSB[.(percentSPR)*"%"]))) +
#       ylab(bquote(paste("F/", F[.(percentSPR)*"%"]))) +
#       ggtitle(paste0("Stock Status in the Last ", use.n.years, " Year(s)")) +
#       theme_bw() +
#       theme(axis.text = element_text(size = 12),
#             axis.title = element_text(size = 20),
#             plot.title = element_text(size = 12),
#             strip.text = element_text(size = 12),
#             legend.text = element_text(size = 12),
#             legend.title = element_text(size = 12),
#             aspect.ratio = 1,
#             panel.grid.major = element_blank(),
#             panel.grid.minor = element_blank(),
#             panel.background = element_blank(),
#             axis.line = element_line(colour = "black")) +
#       coord_cartesian(xlim = c(0, quantile(temp$Overfished, 0.95, na.rm = TRUE)),
#                       ylim = c(0, quantile(temp$Overfishing, 0.95, na.rm = TRUE)))
#   } else {
#     p <- p +
#       geom_point(aes(color = Model), size = 1.5, alpha = 0.8) +
#       scale_color_viridis_d(option = col.opt) +
#       geom_vline(xintercept = 0.5, linetype = "dashed", color = "red", size = 1) +
#       geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 1) +
#       xlab(bquote(paste("SSB/", SSB[.(percentSPR)*"%"]))) +
#       ylab(bquote(paste("F/", F[.(percentSPR)*"%"]))) +
#       ggtitle(paste0("Stock Status in the Last ", use.n.years, " Year(s)")) +
#       theme_bw() +
#       theme(axis.text = element_text(size = 12),
#             axis.title = element_text(size = 20),
#             plot.title = element_text(size = 12),
#             strip.text = element_text(size = 12),
#             legend.text = element_text(size = 12),
#             legend.title = element_text(size = 12),
#             aspect.ratio = 1,
#             panel.grid.major = element_blank(),
#             panel.grid.minor = element_blank(),
#             panel.background = element_blank(),
#             axis.line = element_line(colour = "black")) +
#       coord_cartesian(xlim = c(0, quantile(temp$Overfished, 0.95, na.rm = TRUE)),
#                       ylim = c(0, quantile(temp$Overfishing, 0.95, na.rm = TRUE)))
#   }
#   
#   # Save plot
#   plot_name = paste0("Kobe_Plot_KDE_", use.n.years, "_Year.png")
#   
#   # Create the new subfolder if it doesn't exist
#   new_sub_dir <- file.path(main.dir, sub.dir, "KOBE_Plot")
#   
#   if (!file.exists(new_sub_dir)){
#     dir.create(new_sub_dir)
#   }
#   
#   # Save the figure inside the new subfolder
#   ggsave(file.path(new_sub_dir, plot_name), plot = p, width = width, height = height, dpi = dpi)
# 
#   return(p)
# }

# plot_kobe_status <- function(mods, is.nsim, main.dir, sub.dir, 
#                              width = 10, height = 7, dpi = 300, col.opt = "D",
#                              new_model_names = NULL,
#                              use.n.years = NULL,
#                              show_density = FALSE) {
#   library(dplyr)
#   library(tidyr)
#   library(ggplot2)
#   library(viridis)
#   
#   if (is.null(use.n.years)) {
#     cat("\nuse.n.years is not specified, so default (terminal year) is used here!\n")
#     use.n.years <- 1
#   }
#   
#   if (!is.nsim) {
#     if (is.null(mods[[1]]$om$rep$log_SSB_FXSPR)) {
#       message("Biological Reference Point has not been calculated internally!")
#       return(invisible(NULL))
#     }
#   } else {
#     if (is.null(mods[[1]][[1]]$om$rep$log_SSB_FXSPR)) {
#       message("Biological Reference Point has not been calculated internally!")
#       return(invisible(NULL))
#     }
#   }
#   
#   # === Extract SSB/SSBxx% ===
#   if (!is.nsim) {
#     Years <- mods[[1]]$om$years
#     ssb_list <- lapply(seq_along(mods), function(i) {
#       tmp <- mods[[i]]$om$rep$SSB
#       tmp <- cbind(tmp, rowSums(tmp))
#       tmp <- tmp / exp(mods[[i]]$om$rep$log_SSB_FXSPR)
#       ssb <- tail(tmp, use.n.years)
#       data.frame(Model = paste0("Model", i),
#                  Realization = 1,
#                  Year = tail(Years, use.n.years),
#                  Overfished = ssb)
#     })
#     ssb_df <- bind_rows(ssb_list)
#   } else {
#     Years <- mods[[1]][[1]]$om$years
#     ssb_list <- lapply(seq_along(mods), function(r) {
#       lapply(seq_along(mods[[r]]), function(m) {
#         tmp <- mods[[r]][[m]]$om$rep$SSB
#         tmp <- cbind(tmp, rowSums(tmp))
#         tmp <- tmp / exp(mods[[r]][[m]]$om$rep$log_SSB_FXSPR)
#         ssb <- tail(tmp, use.n.years)
#         data.frame(Model = paste0("Model", m),
#                    Realization = r,
#                    Year = tail(Years, use.n.years),
#                    Overfished = ssb)
#       })
#     })
#     ssb_df <- bind_rows(ssb_list)
#   }
#   
#   ncol = ncol(ssb_df) - 3
#   ssb <- list()
#   for (i in 1:ncol) {
#     ssb[[i]] <- ssb_df[,c(1:3,3+i)]
#     names(ssb[[i]])[4] <- "Overfished"
#   }
#   
#   # === Extract F/Fxx% ===
#   if (!is.nsim) {
#     fbar_list <- lapply(seq_along(mods), function(i) {
#       n_fleets <- mods[[i]]$om$input$data$n_fleets
#       fbar <- mods[[i]]$om$rep$Fbar[, -c(1:n_fleets)]
#       fbar_ref <- exp(mods[[i]]$om$rep$log_Fbar_XSPR[, -c(1:n_fleets)])
#       ff <- tail(fbar / fbar_ref, use.n.years)
#       data.frame(Model = paste0("Model", i),
#                  Realization = 1,
#                  Year = tail(Years, use.n.years),
#                  Overfishing = ff)
#     })
#     fbar_df <- bind_rows(fbar_list)
#   } else {
#     fbar_list <- lapply(seq_along(mods), function(r) {
#       lapply(seq_along(mods[[r]]), function(m) {
#         n_fleets <- mods[[r]][[m]]$om$input$data$n_fleets
#         fbar <- mods[[r]][[m]]$om$rep$Fbar[, -c(1:n_fleets)]
#         fbar_ref <- exp(mods[[r]][[m]]$om$rep$log_Fbar_XSPR[, -c(1:n_fleets)])
#         ff <- tail(fbar / fbar_ref, use.n.years)
#         data.frame(Model = paste0("Model", m),
#                    Realization = r,
#                    Year = tail(Years, use.n.years),
#                    Overfishing = ff)
#       })
#     })
#     fbar_df <- bind_rows(fbar_list)
#   }
#   
#   ncol = ncol(fbar_df) - 3
#   fbar <- list()
#   for (i in 1:ncol) {
#     fbar[[i]] <- fbar_df[,c(1:3,3+i)]
#     names(fbar[[i]])[4] <- "Overfishing"
#   }
#   
#   p <- list()
#   # === Merge ===
#   for (i in 1:ncol) {
#     temp <- left_join(ssb[[i]], fbar[[i]], by = c("Model", "Realization", "Year"))
#     temp$Index <- paste0("Year ", temp$Year)
#     
#     # === Rename Models if needed ===
#     if (!is.null(new_model_names)) {
#       if (length(new_model_names) != length(unique(temp$Model))) {
#         stop("Length of new_model_names must match number of models.")
#       }
#       temp$Model <- factor(temp$Model,
#                            levels = paste0("Model", seq_along(new_model_names)),
#                            labels = new_model_names)
#     }
#     
#     percentSPR <- if (!is.nsim) mods[[1]]$om$input$data$percentSPR else mods[[1]][[1]]$om$input$data$percentSPR
#     
#     # Check if any rows contain non-finite values in the specified columns
#     if (any(!is.finite(temp$Overfished) | !is.finite(temp$Overfishing))) {
#       # If non-finite values exist, print this message
#       message("Non-finite values detected. Removing them before plotting.")
#     }
#     
#     # This filtering step runs regardless of whether the message was printed
#     temp <- temp %>% 
#       filter(is.finite(Overfished) & is.finite(Overfishing))
#     
#     # === Plot
#     p[[i]] <- ggplot(temp, aes(x = Overfished, y = Overfishing)) +
#       facet_wrap(~ Model) +
#       annotate('rect', xmin = 0.5, xmax = Inf, ymin = -Inf, ymax = 1, alpha = 0.2, fill = "yellow")
#     
#     # === Add smoothed 2D density layer
#     if (show_density) {
#       p[[i]] <- p[[i]] + geom_density_2d_filled(aes(fill = after_stat(level)), alpha = 0.4, contour_var = "density", bins = 100) +
#         scale_fill_viridis_d(option = col.opt) + guides(fill = "none") 
#     }
#     
#     # === Add points and reference lines
#     if (show_density) {
#       
#       plot_title <- if (i < ncol) {
#         paste0("Stock Status in the Last ", use.n.years, " Year(s)\nRegion ", i)
#       } else {
#         paste0("Stock Status in the Last ", use.n.years, " Year(s)\nGlobal")
#       }
#       
#       p[[i]] <- p[[i]] +
#         geom_point(aes(color = Model), size = 1.5, alpha = 0.3) +
#         scale_color_viridis_d(option = col.opt) +
#         geom_vline(xintercept = 0.5, linetype = "dashed", color = "red", size = 1) +
#         geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 1) +
#         xlab(bquote(paste("SSB/", SSB[.(percentSPR)*"%"]))) +
#         ylab(bquote(paste("F/", F[.(percentSPR)*"%"]))) +
#         ggtitle(plot_title) +
#         theme_bw() +
#         theme(axis.text = element_text(size = 12),
#               axis.title = element_text(size = 20),
#               plot.title = element_text(size = 12),
#               strip.text = element_text(size = 12),
#               legend.text = element_text(size = 12),
#               legend.title = element_text(size = 12),
#               aspect.ratio = 1,
#               panel.grid.major = element_blank(),
#               panel.grid.minor = element_blank(),
#               panel.background = element_blank(),
#               axis.line = element_line(colour = "black")) +
#         coord_cartesian(xlim = c(0, quantile(temp$Overfished, 0.95, na.rm = TRUE)),
#                         ylim = c(0, quantile(temp$Overfishing, 0.95, na.rm = TRUE)))
#     } else {
#       
#       plot_title <- if (i < ncol) {
#         paste0("Stock Status in the Last ", use.n.years, " Year(s)\nRegion ", i)
#       } else {
#         paste0("Stock Status in the Last ", use.n.years, " Year(s)\nGlobal")
#       }
#       
#       p[[i]] <- p[[i]] +
#         geom_point(aes(color = Model), size = 1.5, alpha = 0.8) +
#         scale_color_viridis_d(option = col.opt) +
#         geom_vline(xintercept = 0.5, linetype = "dashed", color = "red", size = 1) +
#         geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 1) +
#         xlab(bquote(paste("SSB/", SSB[.(percentSPR)*"%"]))) +
#         ylab(bquote(paste("F/", F[.(percentSPR)*"%"]))) +
#         ggtitle(plot_title) +
#         theme_bw() +
#         theme(axis.text = element_text(size = 12),
#               axis.title = element_text(size = 20),
#               plot.title = element_text(size = 12),
#               strip.text = element_text(size = 12),
#               legend.text = element_text(size = 12),
#               legend.title = element_text(size = 12),
#               aspect.ratio = 1,
#               panel.grid.major = element_blank(),
#               panel.grid.minor = element_blank(),
#               panel.background = element_blank(),
#               axis.line = element_line(colour = "black")) +
#         coord_cartesian(xlim = c(0, quantile(temp$Overfished, 0.95, na.rm = TRUE)),
#                         ylim = c(0, quantile(temp$Overfishing, 0.95, na.rm = TRUE)))
#     }
#     
#     # Save plot
#     plot_name = paste0("Kobe_Plot_KDE_", use.n.years, "_Year_", i, ".png")
#     
#     # Create the new subfolder if it doesn't exist
#     new_sub_dir <- file.path(main.dir, sub.dir, "KOBE_Plot")
#     
#     if (!file.exists(new_sub_dir)){
#       dir.create(new_sub_dir)
#     }
#     
#     # Save the figure inside the new subfolder
#     ggsave(file.path(new_sub_dir, plot_name), plot = p[[i]], width = width, height = height, dpi = dpi)
#     
#   }
#   
#   return(p)
# }

plot_kobe_status <- function(mods, is.nsim, main.dir, sub.dir, 
                             width = 10, height = 7, dpi = 300, col.opt = "D",
                             new_model_names = NULL,
                             use.n.years = NULL,
                             show_density = FALSE,
                             nbins = 20) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(viridis)
  
  if (is.null(use.n.years)) {
    cat("\nuse.n.years is not specified, so default (terminal year) is used here!\n")
    use.n.years <- 1
  }
  
  if (!is.nsim) {
    if (is.null(mods[[1]]$om$rep$log_SSB_FXSPR)) {
      message("Biological Reference Point has not been calculated internally!")
      return(invisible(NULL))
    }
  } else {
    if (is.null(mods[[1]][[1]]$om$rep$log_SSB_FXSPR)) {
      message("Biological Reference Point has not been calculated internally!")
      return(invisible(NULL))
    }
  }
  
  # === Extract SSB/SSBxx% ===
  if (!is.nsim) {
    Years <- mods[[1]]$om$years
    ssb_list <- lapply(seq_along(mods), function(i) {
      tmp <- mods[[i]]$om$rep$SSB
      tmp <- cbind(tmp, rowSums(tmp))
      tmp <- tmp / exp(mods[[i]]$om$rep$log_SSB_FXSPR)
      ssb <- tail(tmp, use.n.years)
      data.frame(Model = paste0("Model", i),
                 Realization = 1,
                 Year = tail(Years, use.n.years),
                 Overfished = ssb)
    })
    ssb_df <- bind_rows(ssb_list)
  } else {
    Years <- mods[[1]][[1]]$om$years
    ssb_list <- lapply(seq_along(mods), function(r) {
      lapply(seq_along(mods[[r]]), function(m) {
        tmp <- mods[[r]][[m]]$om$rep$SSB
        tmp <- cbind(tmp, rowSums(tmp))
        tmp <- tmp / exp(mods[[r]][[m]]$om$rep$log_SSB_FXSPR)
        ssb <- tail(tmp, use.n.years)
        data.frame(Model = paste0("Model", m),
                   Realization = r,
                   Year = tail(Years, use.n.years),
                   Overfished = ssb)
      })
    })
    ssb_df <- bind_rows(ssb_list)
  }
  
  ncol = ncol(ssb_df) - 3
  ssb <- list()
  for (i in 1:ncol) {
    ssb[[i]] <- ssb_df[,c(1:3,3+i)]
    names(ssb[[i]])[4] <- "Overfished"
  }
  
  # === Extract F/Fxx% ===
  if (!is.nsim) {
    fbar_list <- lapply(seq_along(mods), function(i) {
      n_fleets <- mods[[i]]$om$input$data$n_fleets
      fbar <- mods[[i]]$om$rep$Fbar[, -c(1:n_fleets)]
      fbar_ref <- exp(mods[[i]]$om$rep$log_Fbar_XSPR[, -c(1:n_fleets)])
      ff <- tail(fbar / fbar_ref, use.n.years)
      data.frame(Model = paste0("Model", i),
                 Realization = 1,
                 Year = tail(Years, use.n.years),
                 Overfishing = ff)
    })
    fbar_df <- bind_rows(fbar_list)
  } else {
    fbar_list <- lapply(seq_along(mods), function(r) {
      lapply(seq_along(mods[[r]]), function(m) {
        n_fleets <- mods[[r]][[m]]$om$input$data$n_fleets
        fbar <- mods[[r]][[m]]$om$rep$Fbar[, -c(1:n_fleets)]
        fbar_ref <- exp(mods[[r]][[m]]$om$rep$log_Fbar_XSPR[, -c(1:n_fleets)])
        ff <- tail(fbar / fbar_ref, use.n.years)
        data.frame(Model = paste0("Model", m),
                   Realization = r,
                   Year = tail(Years, use.n.years),
                   Overfishing = ff)
      })
    })
    fbar_df <- bind_rows(fbar_list)
  }
  
  ncol = ncol(fbar_df) - 3
  fbar <- list()
  for (i in 1:ncol) {
    fbar[[i]] <- fbar_df[,c(1:3,3+i)]
    names(fbar[[i]])[4] <- "Overfishing"
  }
  
  p <- list()
  # === Merge ===
  for (i in 1:ncol) {
    temp <- left_join(ssb[[i]], fbar[[i]], by = c("Model", "Realization", "Year"))
    temp$Index <- paste0("Year ", temp$Year)
    
    # === Rename Models if needed ===
    if (!is.null(new_model_names)) {
      if (length(new_model_names) != length(unique(temp$Model))) {
        stop("Length of new_model_names must match number of models.")
      }
      temp$Model <- factor(temp$Model,
                           levels = paste0("Model", seq_along(new_model_names)),
                           labels = new_model_names)
    }
    
    percentSPR <- if (!is.nsim) mods[[1]]$om$input$data$percentSPR else mods[[1]][[1]]$om$input$data$percentSPR
    
    # Check if any rows contain non-finite values in the specified columns
    if (any(!is.finite(temp$Overfished) | !is.finite(temp$Overfishing))) {
      # If non-finite values exist, print this message
      message("Non-finite values detected. Removing them before plotting.")
    }
    
    # This filtering step runs regardless of whether the message was printed
    temp <- temp %>% 
      filter(is.finite(Overfished) & is.finite(Overfishing))
    
    # === Plot
    p[[i]] <- ggplot(temp, aes(x = Overfished, y = Overfishing)) +
      facet_wrap(~ Model) 
    
    # === Add points and reference lines
    if (show_density) {
      
      plot_title <- if (i < ncol) {
        paste0("Stock Status in the Last ", use.n.years, " Year(s)\nRegion ", i)
      } else {
        paste0("Stock Status in the Last ", use.n.years, " Year(s)\nGlobal")
      }
      
      p[[i]] <- p[[i]] +
        geom_point(size = 1, alpha = 0.1) +
        scale_color_viridis_d(option = col.opt) +
        geom_vline(xintercept = 0.5, linetype = "dashed", color = "red", size = 1) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 1) +
        xlab(bquote(paste("SSB/", SSB[.(percentSPR)*"%"]))) +
        ylab(bquote(paste("F/", F[.(percentSPR)*"%"]))) +
        ggtitle(plot_title) +
        theme_bw() +
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 20),
              plot.title = element_text(size = 12),
              strip.text = element_text(size = 12),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 12),
              aspect.ratio = 1,
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black")) +
        coord_cartesian(xlim = c(0, quantile(temp$Overfished, 0.95, na.rm = TRUE)),
                        ylim = c(0, quantile(temp$Overfishing, 0.95, na.rm = TRUE)))
      
      # === Add smoothed 2D density layer
      p[[i]] <- p[[i]] + geom_density_2d_filled(aes(fill = after_stat(level)), alpha = 0.9, contour_var = "density", bins = nbins) +
        scale_fill_viridis_d(option = col.opt) + guides(fill = "none") 
      
    } else {
      
      plot_title <- if (i < ncol) {
        paste0("Stock Status in the Last ", use.n.years, " Year(s)\nRegion ", i)
      } else {
        paste0("Stock Status in the Last ", use.n.years, " Year(s)\nGlobal")
      }
      
      p[[i]] <- p[[i]] +
        geom_point(aes(color = Model), size = 1.5, alpha = 0.9) +
        scale_color_viridis_d(option = col.opt) +
        annotate('rect', xmin = 0.5, xmax = Inf, ymin = -Inf, ymax = 1, alpha = 0.2, fill = "yellow") + 
        geom_vline(xintercept = 0.5, linetype = "dashed", color = "red", size = 1) +
        geom_hline(yintercept = 1, linetype = "dashed", color = "red", size = 1) +
        xlab(bquote(paste("SSB/", SSB[.(percentSPR)*"%"]))) +
        ylab(bquote(paste("F/", F[.(percentSPR)*"%"]))) +
        ggtitle(plot_title) +
        theme_bw() +
        theme(axis.text = element_text(size = 12),
              axis.title = element_text(size = 20),
              plot.title = element_text(size = 12),
              strip.text = element_text(size = 12),
              legend.text = element_text(size = 12),
              legend.title = element_text(size = 12),
              aspect.ratio = 1,
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.background = element_blank(),
              axis.line = element_line(colour = "black")) +
        coord_cartesian(xlim = c(0, quantile(temp$Overfished, 0.95, na.rm = TRUE)),
                        ylim = c(0, quantile(temp$Overfishing, 0.95, na.rm = TRUE)))
      
    }
    
    # Save plot
    plot_name = paste0("Kobe_Plot_KDE_", use.n.years, "_Year_", i, ".png")
    
    # Create the new subfolder if it doesn't exist
    new_sub_dir <- file.path(main.dir, sub.dir, "KOBE_Plot")
    
    if (!file.exists(new_sub_dir)){
      dir.create(new_sub_dir)
    }
    
    if (show_density) plot_name = paste0("Kobe_Plot_KDE_", use.n.years, "_Year_", i, "_density.png")
    
    # Save the figure inside the new subfolder
    ggsave(file.path(new_sub_dir, plot_name), plot = p[[i]], width = width, height = height, dpi = dpi)
    
  }
  
  return(p)
}

# plot_model_performance_radar <- function(mods, is.nsim, main.dir, sub.dir, 
#                                          width = 10, height = 10, dpi = 300, col.opt = "D",
#                                          method = NULL,
#                                          use.n.years.first = 5,
#                                          use.n.years.last = 5,
#                                          start.years = 1, 
#                                          new_model_names = NULL) {
#   
#   library(dplyr)
#   library(tidyr)
#   library(fmsb)
#   library(viridis)
#   
#   if (is.nsim) {
#     
#     n_models <- length(mods[[1]])
#     n_reps <- length(mods)
#     results <- list()
#     
#     for (r in seq_len(n_reps)) {
#       tmp <- data.frame(Model = paste0("Model", seq_len(n_models)))
#       
#       for (m in seq_len(n_models)) {
#         rep <- mods[[r]][[m]]$om$rep
#         
#         # Catch and SSB as rowSums
#         catch_ts <- rowSums(rep$pred_catch)
#         ssb_ts <- rowSums(rep$SSB)
#         
#         # Global Fbar
#         n_fleets <- mods[[r]][[m]]$om$input$data$n_fleets[1]
#         n_regions <- mods[[r]][[m]]$om$input$data$n_regions[1]
#         fbar_ts <- rep$Fbar[, n_fleets + n_regions + 1]
#         
#         # First and last period means
#         if(is.null(method)) method = "median"
#         if (method == "median") {
#           tmp$Catch_first[m] <- median(catch_ts[start.years:(start.years + use.n.years.first - 1)])
#           tmp$SSB_first[m] <- median(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
#           tmp$Fbar_first[m] <- median(fbar_ts[start.years:(start.years + use.n.years.first - 1)])
#           tmp$Catch_last[m] <- median(tail(catch_ts, use.n.years.last))
#           tmp$SSB_last[m] <- median(tail(ssb_ts, use.n.years.last))
#           tmp$Fbar_last[m] <- median(tail(fbar_ts, use.n.years.last))
#         } else if (method == "mean") {
#           tmp$Catch_first[m] <- mean(catch_ts[start.years:(start.years + use.n.years.first - 1)])
#           tmp$SSB_first[m] <- mean(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
#           tmp$Fbar_first[m] <- mean(fbar_ts[start.years:(start.years + use.n.years.first - 1)])
#           tmp$Catch_last[m] <- mean(tail(catch_ts, use.n.years.last))
#           tmp$SSB_last[m] <- mean(tail(ssb_ts, use.n.years.last))
#           tmp$Fbar_last[m] <- mean(tail(fbar_ts, use.n.years.last))
#         } 
#       }
#       
#       # Normalize: higher is better except for Fbar
#       for (v in c("Catch_first", "SSB_first", "Catch_last", "SSB_last")) {
#         range_val <- max(tmp[[v]]) - min(tmp[[v]])
#         if (range_val == 0) {
#           tmp[[v]] <- 100
#         } else {
#           tmp[[v]] <- 100 * (tmp[[v]] - min(tmp[[v]])) / range_val
#         }
#       }
#       for (v in c("Fbar_first", "Fbar_last")) {
#         range_val <- max(tmp[[v]]) - min(tmp[[v]])
#         if (range_val == 0) {
#           tmp[[v]] <- 100
#         } else {
#           norm_f <- (tmp[[v]] - min(tmp[[v]])) / range_val
#           tmp[[v]] <- 100 * (1 - norm_f)
#         }
#       }
#       
#       results[[r]] <- tmp
#     }
#     
#     # Combine and compute median across realizations
#     combined <- bind_rows(results, .id = "Realization")
#     scores_median <- combined %>%
#       group_by(Model) %>%
#       summarise(across(-Realization, median), .groups = "drop")
#     
#     # Optional renaming
#     if (!is.null(new_model_names)) {
#       if (length(new_model_names) != length(unique(scores_median$Model))) {
#         stop("Length of new_model_names must match the number of models.")
#       }
#       scores_median$Model <- factor(scores_median$Model,
#                                     levels = paste0("Model", seq_along(new_model_names)),
#                                     labels = new_model_names)
#     }
#     
#     # Radar chart prep
#     plot_df <- as.data.frame(scores_median)         # Ensure it's a base data.frame
#     rownames(plot_df) <- plot_df$Model              # Set row names safely
#     plot_df$Model <- NULL                           # Drop Model column
#     plot_df <- as.data.frame(t(plot_df))            # Transpose for radar chart
#     
#     # Add min/max rows for radar scale
#     plot_df <- rbind(rep(100, ncol(plot_df)), rep(0, ncol(plot_df)), plot_df)
#     
#   } else {
#     
#     n_models <- length(mods)
#     results <- list()
#     
#     tmp <- data.frame(Model = paste0("Model", seq_len(n_models)))
#     
#     for (m in seq_len(n_models)) {
#       
#       rep <- mods[[m]]$om$rep
#       
#       # Catch and SSB as rowSums
#       catch_ts <- rowSums(rep$pred_catch)
#       ssb_ts <- rowSums(rep$SSB)
#       
#       # Global Fbar
#       n_fleets <- mods[[m]]$om$input$data$n_fleets[1]
#       n_regions <- mods[[m]]$om$input$data$n_regions[1]
#       fbar_ts <- rep$Fbar[, n_fleets + n_regions + 1]
#       
#       if(is.null(method)) method = "median"
#       # First and last period means
#       if (method == "median") {
#         tmp$Catch_first[m] <- median(catch_ts[start.years:(start.years + use.n.years.first - 1)])
#         tmp$SSB_first[m] <- median(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
#         tmp$Fbar_first[m] <- median(fbar_ts[start.years:(start.years + use.n.years.first - 1)])
#         tmp$Catch_last[m] <- median(tail(catch_ts, use.n.years.last))
#         tmp$SSB_last[m] <- median(tail(ssb_ts, use.n.years.last))
#         tmp$Fbar_last[m] <- median(tail(fbar_ts, use.n.years.last))
#       } else if (method == "mean") {
#         tmp$Catch_first[m] <- mean(catch_ts[start.years:(start.years + use.n.years.first - 1)])
#         tmp$SSB_first[m] <- mean(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
#         tmp$Fbar_first[m] <- mean(fbar_ts[start.years:(start.years + use.n.years.first - 1)])
#         tmp$Catch_last[m] <- mean(tail(catch_ts, use.n.years.last))
#         tmp$SSB_last[m] <- mean(tail(ssb_ts, use.n.years.last))
#         tmp$Fbar_last[m] <- mean(tail(fbar_ts, use.n.years.last))
#       }
#     }
#     
#     # Normalize: higher is better except for Fbar
#     for (v in c("Catch_first", "SSB_first", "Catch_last", "SSB_last")) {
#       range_val <- max(tmp[[v]]) - min(tmp[[v]])
#       if (range_val == 0) {
#         tmp[[v]] <- 100  # All values are identical; assign full score
#       } else {
#         tmp[[v]] <- 100 * (tmp[[v]] - min(tmp[[v]])) / range_val
#       }
#     }
#     
#     for (v in c("Fbar_first", "Fbar_last")) {
#       range_val <- max(tmp[[v]]) - min(tmp[[v]])
#       if (range_val == 0) {
#         tmp[[v]] <- 100  # All values identical; assign full score
#       } else {
#         norm_f <- (tmp[[v]] - min(tmp[[v]])) / range_val
#         tmp[[v]] <- 100 * (1 - norm_f)
#       }
#     }
#     
#     results[[1]] <- tmp
#     
#     # Combine and compute median across realizations
#     combined <- bind_rows(results, .id = "Realization")
#     scores_median <- combined %>%
#       group_by(Model) %>%
#       summarise(across(-Realization, median), .groups = "drop")
#     
#     # Optional renaming
#     if (!is.null(new_model_names)) {
#       if (length(new_model_names) != length(unique(scores_median$Model))) {
#         stop("Length of new_model_names must match the number of models.")
#       }
#       scores_median$Model <- factor(scores_median$Model,
#                                     levels = paste0("Model", seq_along(new_model_names)),
#                                     labels = new_model_names)
#     }
#     
#     # Radar chart prep
#     plot_df <- as.data.frame(scores_median)         # Ensure it's a base data.frame
#     rownames(plot_df) <- plot_df$Model              # Set row names safely
#     plot_df$Model <- NULL                           # Drop Model column
#     plot_df <- as.data.frame(t(plot_df))            # Transpose for radar chart
#     
#     # Add min/max rows for radar scale
#     plot_df <- rbind(rep(100, ncol(plot_df)), rep(0, ncol(plot_df)), plot_df)
#     
#   }
#   
#   if (ncol(plot_df) < 3) {
#     message("Radar chart needs at least 3 models. Showing barplot instead.")
#     return(invisible(NULL))
#   }
#   
#   # Plot without any further changes
#   colors <- viridisLite::viridis(n = nrow(plot_df) - 2, option = col.opt) # -2 for min/max rows
#   
#   # Create the new subfolder if it doesn't exist
#   new_sub_dir <- file.path(main.dir, sub.dir, "Radar_Holistic_Plot")
#   
#   if (!file.exists(new_sub_dir)){
#     dir.create(new_sub_dir)
#   }
#   
#   # Correct way to assign expressions to rownames
#   my_legend_labels <- c(expression(Catch[ST]),
#                         expression(SSB[ST]),
#                         expression(F[ST]),
#                         expression(Catch[LT]),
#                         expression(SSB[LT]),
#                         expression(F[LT]))
#   
#   # Save to PNG
#   output_file <- file.path(file.path(main.dir, sub.dir, "Radar_Holistic_Plot", "model_performance_radar.png"))
#   png(filename = output_file, width = width, height = height, units = "in", res = dpi)
#   radarchart(plot_df,
#              axistype = 4,
#              pcol = colors,
#              plwd = 3,
#              # pty = 1:length(colors),
#              plty = 1:length(colors),
#              cglcol = "grey80",
#              cglty = 1,
#              axislabcol = "grey30",
#              vlcex = 0.9)
#   legend(x = 1.0, y = 0.3, legend = my_legend_labels, col = colors,
#          lty = 1:length(colors), 
#          # pch = 1:length(colors),
#          lwd = 3, cex = 0.9, y.intersp = 1.5)
#   dev.off()
#   
#   # Also plot to screen (inline)
#   op <- par(mar = c(1, 1, 3, 1))  # NEW: minimize margins
#   radarchart(plot_df,
#              axistype = 4,
#              pcol = colors,
#              plwd = 3,
#              # pty = 1:length(colors),
#              plty = 1:length(colors),
#              cglcol = "grey80",
#              cglty = 1,
#              axislabcol = "grey30",
#              vlcex = 0.9)
#   legend(x = 1.0, y = 0.3, legend = my_legend_labels, col = colors,
#          lty = 1:length(colors), 
#          # pch = 1:length(colors),
#          lwd = 3, cex = 0.9, y.intersp = 1.5)
#   on.exit(par(op))  # restore after plotting
# }

# plot_model_performance_radar2 <- function(mods, is.nsim, main.dir, sub.dir, 
#                                           width = 10, height = 10, dpi = 300, col.opt = "D",
#                                           method = NULL,
#                                           use.n.years.first = 5,
#                                           use.n.years.last = 5,
#                                           start.years = 1, 
#                                           new_model_names = NULL) {
#   
#   library(dplyr)
#   library(tidyr)
#   library(fmsb)
#   library(viridis)
#   library(viridisLite)
#   
#   # Helper: calculate Average Annual Catch Variation (AACV)
#   calculate_aacv <- function(values) {
#     if (!is.numeric(values)) {
#       stop("Input must be a numeric vector.")
#     }
#     diffs <- abs(diff(values))
#     aacv <- sum(diffs) / sum(values[-length(values)])
#     return(aacv)
#   }
#   
#   if (is.nsim) {
#     
#     n_models <- length(mods[[1]])
#     n_reps <- length(mods)
#     
#     results <- list()
#     
#     for (r in seq_len(n_reps)) {
#       tmp <- data.frame(Model = paste0("Model", seq_len(n_models)))
#       
#       for (m in seq_len(n_models)) {
#         rep <- mods[[r]][[m]]$om$rep
#         
#         # Catch and SSB as rowSums
#         catch_ts <- rowSums(rep$pred_catch)
#         ssb_ts <- rowSums(rep$SSB)
#         
#         # Global Fbar
#         n_fleets <- mods[[r]][[m]]$om$input$data$n_fleets[1]
#         n_regions <- mods[[r]][[m]]$om$input$data$n_regions[1]
#         
#         # Catch variation
#         catch_aacv <- calculate_aacv(rowSums(rep$pred_catch))
#         
#         # SSB variation
#         ssb_aacv <- calculate_aacv(rowSums(rep$SSB))
#         
#         # Fbar variation (global is last col)
#         fbar_aacv <- calculate_aacv(rep$Fbar[, ncol(rep$Fbar)])
#         
#         # First and last period means
#         if(is.null(method)) method = "median"
#         if (method == "median") {
#           tmp$Catch_first[m] <- median(catch_ts[start.years:(start.years + use.n.years.first - 1)])
#           tmp$SSB_first[m] <- median(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
#           tmp$Catch_last[m] <- median(tail(catch_ts, use.n.years.last))
#           tmp$SSB_last[m] <- median(tail(ssb_ts, use.n.years.last))
#           tmp$catch_aacv[m] <- catch_aacv
#           tmp$ssb_aacv[m] <- ssb_aacv
#           tmp$fbar_aacv[m] <- fbar_aacv
#         } else if (method == "mean") {
#           tmp$Catch_first[m] <- mean(catch_ts[start.years:(start.years + use.n.years.first - 1)])
#           tmp$SSB_first[m] <- mean(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
#           tmp$Catch_last[m] <- mean(tail(catch_ts, use.n.years.last))
#           tmp$SSB_last[m] <- mean(tail(ssb_ts, use.n.years.last))
#           tmp$catch_aacv[m] <- catch_aacv
#           tmp$ssb_aacv[m] <- ssb_aacv
#           tmp$fbar_aacv[m] <- fbar_aacv
#         } 
#       }
#       
#       # Normalize: higher is better except for Fbar
#       for (v in c("Catch_first", "SSB_first", "Catch_last", "SSB_last")) {
#         range_val <- max(tmp[[v]]) - min(tmp[[v]])
#         if (range_val == 0) {
#           tmp[[v]] <- 100
#         } else {
#           tmp[[v]] <- 100 * (tmp[[v]] - min(tmp[[v]])) / range_val
#         }
#       }
#       for (v in c("catch_aacv", "ssb_aacv", "fbar_aacv")) {
#         range_val <- max(tmp[[v]]) - min(tmp[[v]])
#         if (range_val == 0) {
#           tmp[[v]] <- 100
#         } else {
#           norm_f <- (tmp[[v]] - min(tmp[[v]])) / range_val
#           tmp[[v]] <- 100 * (1 - norm_f)
#         }
#       }
#       
#       results[[r]] <- tmp
#     }
#     
#     # Combine and compute median across realizations
#     combined <- bind_rows(results, .id = "Realization")
#     scores_median <- combined %>%
#       group_by(Model) %>%
#       summarise(across(-Realization, median), .groups = "drop")
#     
#     # Optional renaming
#     if (!is.null(new_model_names)) {
#       if (length(new_model_names) != length(unique(scores_median$Model))) {
#         stop("Length of new_model_names must match the number of models.")
#       }
#       scores_median$Model <- factor(scores_median$Model,
#                                     levels = paste0("Model", seq_along(new_model_names)),
#                                     labels = new_model_names)
#     }
#     
#     # Radar chart prep
#     plot_df <- as.data.frame(scores_median)         # Ensure it's a base data.frame
#     rownames(plot_df) <- plot_df$Model              # Set row names safely
#     plot_df$Model <- NULL                           # Drop Model column
#     plot_df <- as.data.frame(t(plot_df))            # Transpose for radar chart
#     
#     # Add min/max rows for radar scale
#     plot_df <- rbind(rep(100, ncol(plot_df)), rep(0, ncol(plot_df)), plot_df)
#     
#   } else {
#     
#     n_models <- length(mods)
#     results <- list()
#     
#     tmp <- data.frame(Model = paste0("Model", seq_len(n_models)))
#     
#     for (m in seq_len(n_models)) {
#       
#       rep <- mods[[m]]$om$rep
#       
#       # Catch and SSB as rowSums
#       catch_ts <- rowSums(rep$pred_catch)
#       ssb_ts <- rowSums(rep$SSB)
#       
#       # Global Fbar
#       n_fleets <- mods[[m]]$om$input$data$n_fleets[1]
#       n_regions <- mods[[m]]$om$input$data$n_regions[1]
#       
#       # Catch variation
#       catch_aacv <- calculate_aacv(rowSums(rep$pred_catch))
#       
#       # SSB variation
#       ssb_aacv <- calculate_aacv(rowSums(rep$SSB))
#       
#       # Fbar variation (global is last col)
#       fbar_aacv <- calculate_aacv(rep$Fbar[, ncol(rep$Fbar)])
#       
#       if(is.null(method)) method = "median"
#       # First and last period means
#       if (method == "median") {
#         tmp$Catch_first[m] <- median(catch_ts[start.years:(start.years + use.n.years.first - 1)])
#         tmp$SSB_first[m] <- median(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
#         tmp$Catch_last[m] <- median(tail(catch_ts, use.n.years.last))
#         tmp$SSB_last[m] <- median(tail(ssb_ts, use.n.years.last))
#         tmp$catch_aacv[m] <- catch_aacv
#         tmp$ssb_aacv[m] <- ssb_aacv
#         tmp$fbar_aacv[m] <- fbar_aacv
#       } else if (method == "mean") {
#         tmp$Catch_first[m] <- mean(catch_ts[start.years:(start.years + use.n.years.first - 1)])
#         tmp$SSB_first[m] <- mean(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
#         tmp$Catch_last[m] <- mean(tail(catch_ts, use.n.years.last))
#         tmp$SSB_last[m] <- mean(tail(ssb_ts, use.n.years.last))
#         tmp$catch_aacv[m] <- catch_aacv
#         tmp$ssb_aacv[m] <- ssb_aacv
#         tmp$fbar_aacv[m] <- fbar_aacv
#       }
#     }
#     
#     # Normalize: higher is better except for Fbar
#     for (v in c("Catch_first", "SSB_first", "Catch_last", "SSB_last")) {
#       range_val <- max(tmp[[v]]) - min(tmp[[v]])
#       if (range_val == 0) {
#         tmp[[v]] <- 100  # All values are identical; assign full score
#       } else {
#         tmp[[v]] <- 100 * (tmp[[v]] - min(tmp[[v]])) / range_val
#       }
#     }
#     
#     for (v in c("catch_aacv", "ssb_aacv", "fbar_aacv")) {
#       range_val <- max(tmp[[v]]) - min(tmp[[v]])
#       if (range_val == 0) {
#         tmp[[v]] <- 100  # All values identical; assign full score
#       } else {
#         norm_f <- (tmp[[v]] - min(tmp[[v]])) / range_val
#         tmp[[v]] <- 100 * (1 - norm_f)
#       }
#     }
#     
#     results[[1]] <- tmp
#     
#     # Combine and compute median across realizations
#     combined <- bind_rows(results, .id = "Realization")
#     scores_median <- combined %>%
#       group_by(Model) %>%
#       summarise(across(-Realization, median), .groups = "drop")
#     
#     # Optional renaming
#     if (!is.null(new_model_names)) {
#       if (length(new_model_names) != length(unique(scores_median$Model))) {
#         stop("Length of new_model_names must match the number of models.")
#       }
#       scores_median$Model <- factor(scores_median$Model,
#                                     levels = paste0("Model", seq_along(new_model_names)),
#                                     labels = new_model_names)
#     }
#     
#     # Radar chart prep
#     plot_df <- as.data.frame(scores_median)         # Ensure it's a base data.frame
#     rownames(plot_df) <- plot_df$Model              # Set row names safely
#     plot_df$Model <- NULL                           # Drop Model column
#     plot_df <- as.data.frame(t(plot_df))            # Transpose for radar chart
#     
#     # Add min/max rows for radar scale
#     plot_df <- rbind(rep(100, ncol(plot_df)), rep(0, ncol(plot_df)), plot_df)
#     
#   }
#   
#   if (ncol(plot_df) < 3) {
#     message("Radar chart needs at least 3 models. Showing barplot instead.")
#     return(invisible(NULL))
#   }
#   
#   # Plot without any further changes
#   colors <- viridisLite::viridis(n = nrow(plot_df) - 2, option = col.opt) # -2 for min/max rows
#   
#   # Create the new subfolder if it doesn't exist
#   new_sub_dir <- file.path(main.dir, sub.dir, "Radar_Holistic_Plot")
#   
#   if (!file.exists(new_sub_dir)){
#     dir.create(new_sub_dir)
#   }
#   
#   rownames(plot_df)
#   # Correct way to assign expressions to rownames
#   my_legend_labels <- c(expression(Catch[ST]),
#                         expression(SSB[ST]),
#                         expression(Catch[LT]),
#                         expression(SSB[LT]),
#                         expression(Catch[AAV]),
#                         expression(SSB[AAV]),
#                         expression(F[AAV]))
#   
#   # Save to PNG
#   output_file <- file.path(file.path(main.dir, sub.dir, "Radar_Holistic_Plot", "model_performance_radar2.png"))
#   png(filename = output_file, width = width, height = height, units = "in", res = dpi)
#   radarchart(plot_df,
#              axistype = 4,
#              pcol = colors,
#              plwd = 3,
#              # pty = 1:length(colors),
#              plty = 1:length(colors),
#              cglcol = "grey80",
#              cglty = 1,
#              axislabcol = "grey30",
#              vlcex = 0.9)
#   legend(x = 1.0, y = 0.3, legend = my_legend_labels, col = colors,
#          lty = 1:length(colors), 
#          # pch = 1:length(colors),
#          lwd = 3, cex = 0.9, y.intersp = 1.5)
#   dev.off()
#   
#   # Also plot to screen (inline)
#   op <- par(mar = c(1, 1, 3, 1))  # NEW: minimize margins
#   radarchart(plot_df,
#              axistype = 4,
#              pcol = colors,
#              plwd = 3,
#              # pty = 1:length(colors),
#              plty = 1:length(colors),
#              cglcol = "grey80",
#              cglty = 1,
#              axislabcol = "grey30",
#              vlcex = 0.9)
#   legend(x = 1.0, y = 0.3, legend = my_legend_labels, col = colors,
#          lty = 1:length(colors), 
#          # pch = 1:length(colors),
#          lwd = 3, cex = 0.9, y.intersp = 1.5)
#   on.exit(par(op))  # restore after plotting
# }

# plot_model_performance_radar <- function(mods, is.nsim, main.dir, sub.dir, 
#                                          width = 10, height = 10, dpi = 300, col.opt = "D",
#                                          method = NULL,
#                                          use.n.years.first = 5,
#                                          use.n.years.last = 5,
#                                          start.years = 1, 
#                                          new_model_names = NULL) {
#   
#   library(dplyr)
#   library(tidyr)
#   library(fmsb)
#   library(viridis)
#   
#   if (is.nsim) {
#     
#     n_models <- length(mods[[1]])
#     n_reps <- length(mods)
#     results <- list()
#     
#     for (r in seq_len(n_reps)) {
#       tmp <- data.frame(Model = paste0("Model", seq_len(n_models)))
#       
#       for (m in seq_len(n_models)) {
#         rep <- mods[[r]][[m]]$om$rep
#         
#         # Catch and SSB as rowSums
#         catch_ts <- rowSums(rep$pred_catch)
#         ssb_ts <- rowSums(rep$SSB)
#         
#         # Global Fbar
#         n_fleets <- mods[[r]][[m]]$om$input$data$n_fleets[1]
#         n_regions <- mods[[r]][[m]]$om$input$data$n_regions[1]
#         fbar_ts <- rep$Fbar[, n_fleets + n_regions + 1]
#         
#         # First and last period means
#         if(is.null(method)) method = "median"
#         if (method == "median") {
#           tmp$Catch_first[m] <- median(catch_ts[start.years:(start.years + use.n.years.first - 1)])
#           tmp$SSB_first[m] <- median(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
#           tmp$Fbar_first[m] <- median(fbar_ts[start.years:(start.years + use.n.years.first - 1)])
#           tmp$Catch_last[m] <- median(tail(catch_ts, use.n.years.last))
#           tmp$SSB_last[m] <- median(tail(ssb_ts, use.n.years.last))
#           tmp$Fbar_last[m] <- median(tail(fbar_ts, use.n.years.last))
#         } else if (method == "mean") {
#           tmp$Catch_first[m] <- mean(catch_ts[start.years:(start.years + use.n.years.first - 1)])
#           tmp$SSB_first[m] <- mean(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
#           tmp$Fbar_first[m] <- mean(fbar_ts[start.years:(start.years + use.n.years.first - 1)])
#           tmp$Catch_last[m] <- mean(tail(catch_ts, use.n.years.last))
#           tmp$SSB_last[m] <- mean(tail(ssb_ts, use.n.years.last))
#           tmp$Fbar_last[m] <- mean(tail(fbar_ts, use.n.years.last))
#         } 
#       }
#       
#       # Normalize: higher is better except for Fbar
#       for (v in c("Catch_first", "SSB_first", "Catch_last", "SSB_last")) {
#         range_val <- max(tmp[[v]]) - min(tmp[[v]])
#         if (range_val == 0) {
#           tmp[[v]] <- 100
#         } else {
#           tmp[[v]] <- 100 * (tmp[[v]] - min(tmp[[v]])) / range_val
#         }
#       }
#       for (v in c("Fbar_first", "Fbar_last")) {
#         range_val <- max(tmp[[v]]) - min(tmp[[v]])
#         if (range_val == 0) {
#           tmp[[v]] <- 100
#         } else {
#           norm_f <- (tmp[[v]] - min(tmp[[v]])) / range_val
#           tmp[[v]] <- 100 * (1 - norm_f)
#         }
#       }
#       
#       results[[r]] <- tmp
#     }
#     
#     # Combine and compute median across realizations
#     combined <- bind_rows(results, .id = "Realization")
#     scores_median <- combined %>%
#       group_by(Model) %>%
#       summarise(across(-Realization, median), .groups = "drop")
#     
#     # Optional renaming
#     if (!is.null(new_model_names)) {
#       if (length(new_model_names) != length(unique(scores_median$Model))) {
#         stop("Length of new_model_names must match the number of models.")
#       }
#       scores_median$Model <- factor(scores_median$Model,
#                                     levels = paste0("Model", seq_along(new_model_names)),
#                                     labels = new_model_names)
#     }
#     
#     # Radar chart prep
#     plot_df <- as.data.frame(scores_median)         # Ensure it's a base data.frame
#     rownames(plot_df) <- plot_df$Model              # Set row names safely
#     plot_df$Model <- NULL                           # Drop Model column
#     plot_df <- as.data.frame(t(plot_df))            # Transpose for radar chart
#     
#     # Add min/max rows for radar scale
#     plot_df <- rbind(rep(100, ncol(plot_df)), rep(0, ncol(plot_df)), plot_df)
#     
#   } else {
#     
#     n_models <- length(mods)
#     results <- list()
#     
#     tmp <- data.frame(Model = paste0("Model", seq_len(n_models)))
#     
#     for (m in seq_len(n_models)) {
#       
#       rep <- mods[[m]]$om$rep
#       
#       # Catch and SSB as rowSums
#       catch_ts <- rowSums(rep$pred_catch)
#       ssb_ts <- rowSums(rep$SSB)
#       
#       # Global Fbar
#       n_fleets <- mods[[m]]$om$input$data$n_fleets[1]
#       n_regions <- mods[[m]]$om$input$data$n_regions[1]
#       fbar_ts <- rep$Fbar[, n_fleets + n_regions + 1]
#       
#       if(is.null(method)) method = "median"
#       # First and last period means
#       if (method == "median") {
#         tmp$Catch_first[m] <- median(catch_ts[start.years:(start.years + use.n.years.first - 1)])
#         tmp$SSB_first[m] <- median(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
#         tmp$Fbar_first[m] <- median(fbar_ts[start.years:(start.years + use.n.years.first - 1)])
#         tmp$Catch_last[m] <- median(tail(catch_ts, use.n.years.last))
#         tmp$SSB_last[m] <- median(tail(ssb_ts, use.n.years.last))
#         tmp$Fbar_last[m] <- median(tail(fbar_ts, use.n.years.last))
#       } else if (method == "mean") {
#         tmp$Catch_first[m] <- mean(catch_ts[start.years:(start.years + use.n.years.first - 1)])
#         tmp$SSB_first[m] <- mean(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
#         tmp$Fbar_first[m] <- mean(fbar_ts[start.years:(start.years + use.n.years.first - 1)])
#         tmp$Catch_last[m] <- mean(tail(catch_ts, use.n.years.last))
#         tmp$SSB_last[m] <- mean(tail(ssb_ts, use.n.years.last))
#         tmp$Fbar_last[m] <- mean(tail(fbar_ts, use.n.years.last))
#       }
#     }
#     
#     # Normalize: higher is better except for Fbar
#     for (v in c("Catch_first", "SSB_first", "Catch_last", "SSB_last")) {
#       range_val <- max(tmp[[v]]) - min(tmp[[v]])
#       if (range_val == 0) {
#         tmp[[v]] <- 100  # All values are identical; assign full score
#       } else {
#         tmp[[v]] <- 100 * (tmp[[v]] - min(tmp[[v]])) / range_val
#       }
#     }
#     
#     for (v in c("Fbar_first", "Fbar_last")) {
#       range_val <- max(tmp[[v]]) - min(tmp[[v]])
#       if (range_val == 0) {
#         tmp[[v]] <- 100  # All values identical; assign full score
#       } else {
#         norm_f <- (tmp[[v]] - min(tmp[[v]])) / range_val
#         tmp[[v]] <- 100 * (1 - norm_f)
#       }
#     }
#     
#     results[[1]] <- tmp
#     
#     # Combine and compute median across realizations
#     combined <- bind_rows(results, .id = "Realization")
#     scores_median <- combined %>%
#       group_by(Model) %>%
#       summarise(across(-Realization, median), .groups = "drop")
#     
#     # Optional renaming
#     if (!is.null(new_model_names)) {
#       if (length(new_model_names) != length(unique(scores_median$Model))) {
#         stop("Length of new_model_names must match the number of models.")
#       }
#       scores_median$Model <- factor(scores_median$Model,
#                                     levels = paste0("Model", seq_along(new_model_names)),
#                                     labels = new_model_names)
#     }
#     
#     # Radar chart prep
#     plot_df <- as.data.frame(scores_median)         # Ensure it's a base data.frame
#     rownames(plot_df) <- plot_df$Model              # Set row names safely
#     plot_df$Model <- NULL                           # Drop Model column
#     plot_df <- as.data.frame(t(plot_df))            # Transpose for radar chart
#     
#     # Add min/max rows for radar scale
#     plot_df <- rbind(rep(100, ncol(plot_df)), rep(0, ncol(plot_df)), plot_df)
#     
#   }
#   
#   if (ncol(plot_df) < 3) {
#     warning("Radar chart needs at least 3 models. Showing barplot instead.")
#     return(invisible(NULL))
#   }
#   
#   # Plot without any further changes
#   colors <- viridisLite::viridis(n = nrow(plot_df) - 2, option = col.opt) # -2 for min/max rows
#   
#   # Create the new subfolder if it doesn't exist
#   new_sub_dir <- file.path(main.dir, sub.dir, "Radar_Holistic_Plot")
#   
#   if (!file.exists(new_sub_dir)){
#     dir.create(new_sub_dir)
#   }
#   
#   # Save to PNG
#   output_file <- file.path(file.path(main.dir, sub.dir, "Radar_Holistic_Plot", "model_performance_radar.png"))
#   png(filename = output_file, width = width, height = height, units = "in", res = dpi)
#   radarchart(plot_df,
#              axistype = 4,
#              pcol = colors,
#              plwd = 2,
#              plty = 1,
#              cglcol = "grey80",
#              cglty = 1,
#              axislabcol = "grey30",
#              vlcex = 0.9,
#              title = paste("Model Performance"))
#   legend(x = 1.25, y = 0.3, legend = rownames(plot_df)[-c(1,2)], col = colors,
#          lty = 1, lwd = 2, cex = 0.8, bty = "n")
#   dev.off()
#   
#   # Also plot to screen (inline)
#   op <- par(mar = c(1, 1, 3, 1))  # NEW: minimize margins
#   radarchart(plot_df,
#              axistype = 4,
#              pcol = colors,
#              plwd = 2,
#              plty = 1,
#              cglcol = "grey80",
#              cglty = 1,
#              axislabcol = "grey30",
#              vlcex = 0.9,
#              title = paste("Model Performance"))
#   legend(x = 1.25, y = 0.3, legend = rownames(plot_df)[-c(1,2)], col = colors,
#          lty = 1, lwd = 2, cex = 0.8, bty = "n")
#   on.exit(par(op))  # restore after plotting
# }

plot_mean_rec_par <- function(mods, is.nsim, main.dir, sub.dir, 
                              width = 10, height = 7, dpi = 300, col.opt = "D",
                              outlier.opt = NA,
                              new_model_names = NULL) {
  library(dplyr)
  library(ggplot2)
  library(viridisLite)
  
  # Assume all models have the same OM (e.g., from the first realization)
  om_tmp <- if (is.nsim) mods[[1]][[1]] else mods[[1]]
  mean_rec_true <- exp(om_tmp$om$parList$mean_rec_pars[, 1])
  is.n.regions <- if(om_tmp$om$input$data$n_regions>1) TRUE else FALSE
  
  res <- NULL
  
  if (is.nsim) {
    for (i in seq_along(mods)) {
      for (j in seq_along(mods[[1]])) {
        tmp <- mods[[i]][[j]]
        k <- length(tmp$par.est)
        
        if (any(names(tmp$par.est[[k]]) == "mean_rec_pars")) {
          temp <- exp(tmp$par.est[[k]]$mean_rec_pars[, 1])
        } else {
          m <- length(tmp$par.est[[k]])
          temp <- NULL
          for (n in seq_len(m)) {
            temp1 <- exp(tmp$par.est[[k]][[n]]$mean_rec_pars[, 1])
            temp <- c(temp, temp1)
          }
        }
        
        res1 <- data.frame(
          Model = j,
          nsim = i,
          Value = temp,
          Var = if (length(temp) == 1) "Mean_Rec" else paste0("Mean_Rec_", seq_along(temp))
        )
        
        res <- rbind(res, res1)
      }
    }
  } else {
    for (j in seq_along(mods)) {
      tmp <- mods[[j]]
      k <- length(tmp$par.est)
      
      if (any(names(tmp$par.est[[k]]) == "mean_rec_pars")) {
        temp <- exp(tmp$par.est[[k]]$mean_rec_pars[, 1])
      } else {
        m <- length(tmp$par.est[[k]])
        temp <- NULL
        for (n in seq_len(m)) {
          temp1 <- exp(tmp$par.est[[k]][[n]]$mean_rec_pars[, 1])
          temp <- c(temp, temp1)
        }
      }
      
      res1 <- data.frame(
        Model = j,
        Value = temp,
        Var = if (length(temp) == 1) "Mean_Rec" else paste0("Mean_Rec_", seq_along(temp))
      )
      
      res <- rbind(res, res1)
    }
  }
  
  res$Model <- as.factor(paste0("Model",res$Model))
  
  if (!is.null(new_model_names)) {
    if (length(new_model_names) != length(unique(res$Model))) {
      stop("Length of new_model_names must match the number of models.")
    }
    res$Model <- factor(res$Model,
                        levels = paste0("Model", seq_along(new_model_names)),
                        labels = new_model_names)
  }
  
  res$True_Value <- NA
  
  library(stringr)
  
  # Match entries like "Mean_Rec_1", "Mean_Rec_2", etc.
  if(length(unique(res$Var)) > 1) {
    
    is_indexed <- grepl("^Mean_Rec_\\d+$", res$Var)
    
    # Extract index
    rec_idx <- as.numeric(str_extract(res$Var[is_indexed], "\\d+"))
    
    # Assign based on extracted index
    res$True_Value[is_indexed] <- mean_rec_true[rec_idx]
    
    # Handle non-indexed single value
    res$True_Value[res$Var == "Mean_Rec"] <- mean_rec_true[1]
    
    if(length(res$Var) > 1) res$True_Value[res$Var == "Mean_Rec"] = NA
    
  } else {
    
    if (is.n.regions) {
      
      is_indexed <- grepl("Mean_Rec", res$Var)
      
      # Extract index
      rec_idx <- as.numeric(str_extract(res$Var[is_indexed], "\\d+"))
      
      # Assign based on extracted index
      res$True_Value[is_indexed] <- sum(mean_rec_true)
      
    } else {
      
      is_indexed <- grepl("Mean_Rec", res$Var)
      
      # Extract index
      rec_idx <- as.numeric(str_extract(res$Var[is_indexed], "\\d+"))
      
      # Assign based on extracted index
      res$True_Value[is_indexed] <- mean_rec_true[rec_idx]
      
      # Handle non-indexed single value
      res$True_Value[res$Var == "Mean_Rec"] <- mean_rec_true[1]
    }

  }
  
  
  p1 <- ggplot(res, aes(x = Model, y = Value, col = Model)) +
    geom_boxplot(lwd = 0.8, outlier.shape = outlier.opt) +
    geom_hline(aes(yintercept = True_Value), col = "red", linetype = "dashed") +
    facet_grid(Var ~ ., scales = "free") +
    scale_color_viridis_d(option = col.opt) +
    ggtitle("Mean Recruitment from the Last EM") +
    ylab("") +
    theme_bw() +
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 15),
      plot.title = element_text(size = 12),
      strip.text = element_text(size = 10, color = "black"),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10)
    )
  
  # Create the new subfolder if it doesn't exist
  new_sub_dir <- file.path(main.dir, sub.dir, "Diagnostic_Results")
  
  if (!file.exists(new_sub_dir)){
    dir.create(new_sub_dir)
  }
  
  ggsave(file.path(main.dir, sub.dir, "Diagnostic_Results", "Mean_rec_par.png"), p1, width = 10, height = 10, dpi = dpi)
  
  return(p1)
}

plot_NAA_sigma_par <- function(mods, is.nsim, main.dir, sub.dir,
                               width = 10, height = 15, dpi = 300, col.opt = "D",
                               outlier.opt = NA,
                               new_model_names = NULL) {
  library(dplyr)
  library(ggplot2)
  library(viridisLite)
  
  res <- NULL
  
  # Assume all models have the same OM (e.g., from the first realization)
  om_tmp <- if (is.nsim) mods[[1]][[1]] else mods[[1]]
  rec_sig_true <- exp(om_tmp$om$parList$log_NAA_sigma[, 1, 1])
  naa_sig_true <- exp(om_tmp$om$parList$log_NAA_sigma[, 1, 2])
  
  if (is.nsim) {
    for (i in seq_along(mods)) {
      for (j in seq_along(mods[[1]])) {
        tmp <- mods[[i]][[j]]
        k <- length(tmp$par.est)
        
        if (any(names(tmp$par.est[[k]]) == "log_NAA_sigma")) {
          rec_sig <- exp(tmp$par.est[[k]]$log_NAA_sigma[, 1, 1])
          naa_sig <- exp(tmp$par.est[[k]]$log_NAA_sigma[, 1, 2])
          temp <- c(rec_sig, naa_sig)
        } else {
          m <- length(tmp$par.est[[k]])
          temp <- NULL
          for (n in seq_len(m)) {
            rec_sig <- exp(tmp$par.est[[k]][[n]]$log_NAA_sigma[, 1, 1])
            naa_sig <- exp(tmp$par.est[[k]][[n]]$log_NAA_sigma[, 1, 2])
            temp <- c(temp, rec_sig, naa_sig)
          }
        }
        
        res1 <- data.frame(
          Model = j,
          nsim = i,
          Value = temp,
          Var = if (length(temp) == 2) {
            c("Rec_sigma", "NAA_sigma")
          } else {
            c(paste0("Rec_sigma", seq_along(rec_sig)), paste0("NAA_sigma", seq_along(naa_sig)))
          }
        )
        res <- rbind(res, res1)
      }
    }
  } else {
    for (j in seq_along(mods)) {
      tmp <- mods[[j]]
      k <- length(tmp$par.est)
      
      if (any(names(tmp$par.est[[k]]) == "log_NAA_sigma")) {
        rec_sig <- exp(tmp$par.est[[k]]$log_NAA_sigma[, 1, 1])
        naa_sig <- exp(tmp$par.est[[k]]$log_NAA_sigma[, 1, 2])
        temp <- c(rec_sig, naa_sig)
      } else {
        m <- length(tmp$par.est[[k]])
        temp <- NULL
        for (n in seq_len(m)) {
          rec_sig <- exp(tmp$par.est[[k]][[n]]$log_NAA_sigma[, 1, 1])
          naa_sig <- exp(tmp$par.est[[k]][[n]]$log_NAA_sigma[, 1, 2])
          temp <- c(temp, rec_sig, naa_sig)
        }
      }
      
      res1 <- data.frame(
        Model = j,
        Value = temp,
        Var = if (length(temp) == 2) {
          c("Rec_sigma", "NAA_sigma")
        } else {
          c(paste0("Rec_sigma", seq_along(rec_sig)), paste0("NAA_sigma", seq_along(naa_sig)))
        }
      )
      
      res <- rbind(res, res1)
    }
  }
  
  res$Model <- as.factor(paste0("Model", res$Model))
  
  if (!is.null(new_model_names)) {
    if (length(new_model_names) != length(unique(res$Model))) {
      stop("Length of new_model_names must match the number of models.")
    }
    res$Model <- factor(res$Model,
                        levels = paste0("Model", seq_along(new_model_names)),
                        labels = new_model_names)
  }
  
  res$True_Value <- NA
  
  is.n.regions <- if(om_tmp$om$input$data$n_regions>1) TRUE else FALSE
  
  if (!is.n.regions) {
    if (length(unique(res$Var)) == 2) {
      # Simple case: one Rec_sigma and one NAA_sigma
      res$True_Value[res$Var == "Rec_sigma"] <- rec_sig_true
      res$True_Value[res$Var == "NAA_sigma"] <- naa_sig_true
    } else {
      # Multiple sigmas, match by index
      res$True_Value[grepl("Rec_sigma", res$Var)] <- rec_sig_true[as.numeric(gsub("Rec_sigma", "", res$Var[grepl("Rec_sigma", res$Var)]))]
      res$True_Value[grepl("NAA_sigma", res$Var)] <- naa_sig_true[as.numeric(gsub("NAA_sigma", "", res$Var[grepl("NAA_sigma", res$Var)]))]
    }
  }
  
  p2 <- ggplot(res, aes(x = Model, y = Value, col = Model)) +
    geom_boxplot(lwd = 0.8, outlier.shape = outlier.opt) +
    {if (!is.n.regions) geom_hline(aes(yintercept = True_Value), col = "red", linetype = "dashed")} +
    facet_grid(Var ~ ., scales = "free") +
    scale_color_viridis_d(option = col.opt) +
    ggtitle("Standard Deviation of NAA from the Last EM") +
    ylab("") +
    theme_bw() +
    theme(
      axis.text = element_text(size = 10),
      axis.title = element_text(size = 15),
      plot.title = element_text(size = 12),
      strip.text = element_text(size = 10, color = "black"),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 10)
    )
  
  # Create the new subfolder if it doesn't exist
  new_sub_dir <- file.path(main.dir, sub.dir, "Diagnostic_Results")
  
  if (!file.exists(new_sub_dir)){
    dir.create(new_sub_dir)
  }
  
  ggsave(file.path(main.dir, sub.dir, "Diagnostic_Results", "Variance_Para_NAA.png"), p2, width = width, height = height, dpi = dpi)
  
  return(p2)
}

plot_model_performance_bar <- function(mods, is.nsim,
                                       main.dir = ".",
                                       sub.dir = ".",
                                       new_model_names = NULL,
                                       width = 12, height = 8, dpi = 300,
                                       col.opt = "D",
                                       method = NULL,
                                       use.n.years.first = 5,
                                       use.n.years.last = 5,
                                       start.years = 1) {

  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(viridisLite)

  results_list <- list()

  if (is.nsim) {
    n_models <- length(mods[[1]])
    n_reps <- length(mods)

    for (r in seq_len(n_reps)) {
      tmp <- data.frame(Model = paste0("Model", seq_len(n_models)))

      for (m in seq_len(n_models)) {
        rep <- mods[[r]][[m]]$om$rep
        input <- mods[[r]][[m]]$om$input$data
        n_fleets <- input$n_fleets[1]
        n_regions <- input$n_regions[1]

        # Collect timeseries
        catch_ts <- rowSums(rep$pred_catch)
        ssb_ts <- rowSums(rep$SSB)
        fbar_ts <- rep$Fbar[, n_fleets + n_regions + 1]

        if(is.null(method)) method = "median"
        # First and last period means
        if (method == "median") {
          tmp$Catch_first[m] <- median(catch_ts[start.years:(start.years + use.n.years.first - 1)])
          tmp$SSB_first[m] <- median(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
          tmp$Fbar_first[m] <- median(fbar_ts[start.years:(start.years + use.n.years.first - 1)])
          tmp$Catch_last[m] <- median(tail(catch_ts, use.n.years.last))
          tmp$SSB_last[m] <- median(tail(ssb_ts, use.n.years.last))
          tmp$Fbar_last[m] <- median(tail(fbar_ts, use.n.years.last))
        } else if (method == "mean") {
          tmp$Catch_first[m] <- mean(catch_ts[start.years:(start.years + use.n.years.first - 1)])
          tmp$SSB_first[m] <- mean(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
          tmp$Fbar_first[m] <- mean(fbar_ts[start.years:(start.years + use.n.years.first - 1)])
          tmp$Catch_last[m] <- mean(tail(catch_ts, use.n.years.last))
          tmp$SSB_last[m] <- mean(tail(ssb_ts, use.n.years.last))
          tmp$Fbar_last[m] <- mean(tail(fbar_ts, use.n.years.last))
        }
      }

      # Normalize within realization (safe check)
      for (v in c("Catch_first", "SSB_first", "Catch_last", "SSB_last")) {
        range_val <- max(tmp[[v]]) - min(tmp[[v]])
        if (range_val == 0) {
          tmp[[v]] <- 100
        } else {
          tmp[[v]] <- 100 * (tmp[[v]] - min(tmp[[v]])) / range_val
        }
      }
      for (v in c("Fbar_first", "Fbar_last")) {
        range_val <- max(tmp[[v]]) - min(tmp[[v]])
        if (range_val == 0) {
          tmp[[v]] <- 100
        } else {
          norm_f <- (tmp[[v]] - min(tmp[[v]])) / range_val
          tmp[[v]] <- 100 * (1 - norm_f)
        }
      }

      tmp$Realization <- r
      results_list[[r]] <- tmp
    }

    combined <- bind_rows(results_list)

    # Summarize across realizations
    summary_data <- combined %>%
      pivot_longer(-c(Model, Realization), names_to = "Metric", values_to = "Score") %>%
      group_by(Model, Metric) %>%
      summarize(
        Median = median(Score, na.rm = TRUE),
        Q1 = quantile(Score, 0.25, na.rm = TRUE),
        Q3 = quantile(Score, 0.75, na.rm = TRUE),
        .groups = "drop"
      )

  } else {
    n_models <- length(mods)
    tmp <- data.frame(Model = paste0("Model", seq_len(n_models)))

    for (m in seq_len(n_models)) {
      rep <- mods[[m]]$om$rep
      input <- mods[[m]]$om$input$data
      n_fleets <- input$n_fleets[1]
      n_regions <- input$n_regions[1]

      # Collect timeseries
      catch_ts <- rowSums(rep$pred_catch)
      ssb_ts <- rowSums(rep$SSB)
      fbar_ts <- rep$Fbar[, n_fleets + n_regions + 1]

      if(is.null(method)) method = "median"
      # First and last period means
      if (method == "median") {
        tmp$Catch_first[m] <- median(catch_ts[start.years:(start.years + use.n.years.first - 1)])
        tmp$SSB_first[m] <- median(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
        tmp$Fbar_first[m] <- median(fbar_ts[start.years:(start.years + use.n.years.first - 1)])
        tmp$Catch_last[m] <- median(tail(catch_ts, use.n.years.last))
        tmp$SSB_last[m] <- median(tail(ssb_ts, use.n.years.last))
        tmp$Fbar_last[m] <- median(tail(fbar_ts, use.n.years.last))
      } else if (method == "mean") {
        tmp$Catch_first[m] <- mean(catch_ts[start.years:(start.years + use.n.years.first - 1)])
        tmp$SSB_first[m] <- mean(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
        tmp$Fbar_first[m] <- mean(fbar_ts[start.years:(start.years + use.n.years.first - 1)])
        tmp$Catch_last[m] <- mean(tail(catch_ts, use.n.years.last))
        tmp$SSB_last[m] <- mean(tail(ssb_ts, use.n.years.last))
        tmp$Fbar_last[m] <- mean(tail(fbar_ts, use.n.years.last))
      }
    }

    # Normalize across models (safe check)
    for (v in c("Catch_first", "SSB_first", "Catch_last", "SSB_last")) {
      range_val <- max(tmp[[v]]) - min(tmp[[v]])
      if (range_val == 0) {
        tmp[[v]] <- 100
      } else {
        tmp[[v]] <- 100 * (tmp[[v]] - min(tmp[[v]])) / range_val
      }
    }
    for (v in c("Fbar_first", "Fbar_last")) {
      range_val <- max(tmp[[v]]) - min(tmp[[v]])
      if (range_val == 0) {
        tmp[[v]] <- 100
      } else {
        norm_f <- (tmp[[v]] - min(tmp[[v]])) / range_val
        tmp[[v]] <- 100 * (1 - norm_f)
      }
    }

    tmp$Realization <- 1
    combined <- tmp

    summary_data <- tmp %>%
      pivot_longer(-c(Model, Realization), names_to = "Metric", values_to = "Median") %>%
      mutate(Q1 = Median, Q3 = Median)
  }

  # Rename if needed
  if (!is.null(new_model_names)) {
    if (length(new_model_names) != length(unique(summary_data$Model))) {
      stop("Length of new_model_names must match number of models.")
    }
    summary_data$Model <- factor(summary_data$Model,
                                 levels = paste0("Model", seq_along(new_model_names)),
                                 labels = new_model_names)
  }

  # Colors
  n_metrics <- length(unique(summary_data$Metric))
  my_colors <- viridisLite::viridis(n = n_metrics, option = col.opt)

  my_metric_labels <- c(
    expression(Catch[ST]), # Catch will be bold
    expression(SSB[ST]),   # SSB will be bold
    expression(F[ST]),     # F will be bold
    expression(Catch[LT]), # Catch will be bold
    expression(SSB[LT]),   # SSB will be bold
    expression(F[LT])      # F will be bold
  )

  metric_levels <- c(
    "Catch_first",
    "SSB_first",
    "Fbar_first", # Note: Fbar_first maps to F[ST]
    "Catch_last",
    "SSB_last",
    "Fbar_last"   # Note: Fbar_last maps to F[LT]
  )

  summary_data$Metric <- factor(summary_data$Metric,
                                levels = metric_levels,
                                labels = my_metric_labels)

  # Plot
  plot <- ggplot(summary_data, aes(x = Median, y = Model, fill = Metric, color = Metric)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.7) +
    geom_errorbar(aes(xmin = Q1, xmax = Q3), col = "black",
                  position = position_dodge(width = 0.9), width = 0.3, alpha = 0.3) +
    scale_fill_manual(values = my_colors,
                      labels = my_metric_labels) + # Add labels argument for fill legend
    scale_color_manual(values = my_colors,
                       labels = my_metric_labels) + # Add labels argument for color legend
    theme_bw() +
    labs(title = paste("Holistic Model Performance"),
         x = "Score (0-100)",
         y = "Estimation Model",
         fill = "Metric",
         color = "Metric") +
    theme(
      axis.text = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      plot.title = element_text(size = 18, face = "bold")
    )

  # Print and save
  print(plot)

  # Create the new subfolder if it doesn't exist
  new_sub_dir <- file.path(main.dir, sub.dir, "Holistic_Bar_Plot")

  if (!file.exists(new_sub_dir)){
    dir.create(new_sub_dir)
  }

  ggsave(filename = file.path(new_sub_dir, "Overall_Performance.png"),
         plot = plot, width = width, height = height, dpi = dpi)
}

plot_model_performance_bar2 <- function(mods, is.nsim,
                                        main.dir = ".",
                                        sub.dir = ".",
                                        new_model_names = NULL,
                                        width = 12, height = 8, dpi = 300,
                                        col.opt = "D",
                                        method = NULL,
                                        use.n.years.first = 5,
                                        use.n.years.last = 5,
                                        start.years = 1) {

  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(viridisLite)

  results_list <- list()

  if (is.nsim) {
    n_models <- length(mods[[1]])
    n_reps <- length(mods)

    for (r in seq_len(n_reps)) {
      tmp <- data.frame(Model = paste0("Model", seq_len(n_models)))

      for (m in seq_len(n_models)) {
        rep <- mods[[r]][[m]]$om$rep
        input <- mods[[r]][[m]]$om$input$data
        n_fleets <- input$n_fleets[1]
        n_regions <- input$n_regions[1]

        # Collect timeseries
        catch_ts <- rowSums(rep$pred_catch)
        ssb_ts <- rowSums(rep$SSB)
        fbar_ts <- rep$Fbar[, n_fleets + n_regions + 1]

        if(is.null(method)) method = "median"
        # First and last period means
        if (method == "median") {
          tmp$Catch_first[m] <- median(catch_ts[start.years:(start.years + use.n.years.first - 1)])
          tmp$SSB_first[m] <- median(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
          tmp$Fbar_first[m] <- median(fbar_ts[start.years:(start.years + use.n.years.first - 1)])
          tmp$Catch_last[m] <- median(tail(catch_ts, use.n.years.last))
          tmp$SSB_last[m] <- median(tail(ssb_ts, use.n.years.last))
          tmp$Fbar_last[m] <- median(tail(fbar_ts, use.n.years.last))
        } else if (method == "mean") {
          tmp$Catch_first[m] <- mean(catch_ts[start.years:(start.years + use.n.years.first - 1)])
          tmp$SSB_first[m] <- mean(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
          tmp$Fbar_first[m] <- mean(fbar_ts[start.years:(start.years + use.n.years.first - 1)])
          tmp$Catch_last[m] <- mean(tail(catch_ts, use.n.years.last))
          tmp$SSB_last[m] <- mean(tail(ssb_ts, use.n.years.last))
          tmp$Fbar_last[m] <- mean(tail(fbar_ts, use.n.years.last))
        }
      }

      # Normalize within realization (safe check)
      for (v in c("Catch_first", "SSB_first", "Catch_last", "SSB_last")) {
        range_val <- max(tmp[[v]]) - min(tmp[[v]])
        if (range_val == 0) {
          tmp[[v]] <- 100
        } else {
          tmp[[v]] <- 100 * (tmp[[v]] - min(tmp[[v]])) / range_val
        }
      }
      for (v in c("Fbar_first", "Fbar_last")) {
        range_val <- max(tmp[[v]]) - min(tmp[[v]])
        if (range_val == 0) {
          tmp[[v]] <- 100
        } else {
          norm_f <- (tmp[[v]] - min(tmp[[v]])) / range_val
          tmp[[v]] <- 100 * (1 - norm_f)
        }
      }

      tmp$Realization <- r
      results_list[[r]] <- tmp
    }

    combined <- bind_rows(results_list)

    # Summarize across realizations
    summary_data <- combined %>%
      pivot_longer(-c(Model, Realization), names_to = "Metric", values_to = "Score") %>%
      group_by(Model, Metric) %>%
      summarize(
        Median = median(Score, na.rm = TRUE),
        Q1 = quantile(Score, 0.25, na.rm = TRUE),
        Q3 = quantile(Score, 0.75, na.rm = TRUE),
        .groups = "drop"
      )

  } else {
    n_models <- length(mods)
    tmp <- data.frame(Model = paste0("Model", seq_len(n_models)))

    for (m in seq_len(n_models)) {
      rep <- mods[[m]]$om$rep
      input <- mods[[m]]$om$input$data
      n_fleets <- input$n_fleets[1]
      n_regions <- input$n_regions[1]

      # Collect timeseries
      catch_ts <- rowSums(rep$pred_catch)
      ssb_ts <- rowSums(rep$SSB)
      fbar_ts <- rep$Fbar[, n_fleets + n_regions + 1]

      if(is.null(method)) method = "median"
      # First and last period means
      if (method == "median") {
        tmp$Catch_first[m] <- median(catch_ts[start.years:(start.years + use.n.years.first - 1)])
        tmp$SSB_first[m] <- median(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
        tmp$Fbar_first[m] <- median(fbar_ts[start.years:(start.years + use.n.years.first - 1)])
        tmp$Catch_last[m] <- median(tail(catch_ts, use.n.years.last))
        tmp$SSB_last[m] <- median(tail(ssb_ts, use.n.years.last))
        tmp$Fbar_last[m] <- median(tail(fbar_ts, use.n.years.last))
      } else if (method == "mean") {
        tmp$Catch_first[m] <- mean(catch_ts[start.years:(start.years + use.n.years.first - 1)])
        tmp$SSB_first[m] <- mean(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
        tmp$Fbar_first[m] <- mean(fbar_ts[start.years:(start.years + use.n.years.first - 1)])
        tmp$Catch_last[m] <- mean(tail(catch_ts, use.n.years.last))
        tmp$SSB_last[m] <- mean(tail(ssb_ts, use.n.years.last))
        tmp$Fbar_last[m] <- mean(tail(fbar_ts, use.n.years.last))
      }
    }

    # Normalize across models (safe check)
    for (v in c("Catch_first", "SSB_first", "Catch_last", "SSB_last")) {
      range_val <- max(tmp[[v]]) - min(tmp[[v]])
      if (range_val == 0) {
        tmp[[v]] <- 100
      } else {
        tmp[[v]] <- 100 * (tmp[[v]] - min(tmp[[v]])) / range_val
      }
    }
    for (v in c("Fbar_first", "Fbar_last")) {
      range_val <- max(tmp[[v]]) - min(tmp[[v]])
      if (range_val == 0) {
        tmp[[v]] <- 100
      } else {
        norm_f <- (tmp[[v]] - min(tmp[[v]])) / range_val
        tmp[[v]] <- 100 * (1 - norm_f)
      }
    }

    tmp$Realization <- 1
    combined <- tmp

    summary_data <- tmp %>%
      pivot_longer(-c(Model, Realization), names_to = "Metric", values_to = "Median") %>%
      mutate(Q1 = Median, Q3 = Median)
  }

  # Rename if needed
  if (!is.null(new_model_names)) {
    if (length(new_model_names) != length(unique(summary_data$Model))) {
      stop("Length of new_model_names must match number of models.")
    }
    summary_data$Model <- factor(summary_data$Model,
                                 levels = paste0("Model", seq_along(new_model_names)),
                                 labels = new_model_names)
  }

  # Colors
  n_metrics <- length(unique(summary_data$Model))
  my_colors <- viridisLite::viridis(n = n_metrics, option = col.opt)

  my_metric_labels <- c(
    expression(bold(Catch[ST])), # Catch will be bold
    expression(bold(SSB[ST])),   # SSB will be bold
    expression(bold(F[ST])),     # F will be bold
    expression(bold(Catch[LT])), # Catch will be bold
    expression(bold(SSB[LT])),   # SSB will be bold
    expression(bold(F[LT]))      # F will be bold
  )

  metric_levels <- c(
    "Catch_first",
    "SSB_first",
    "Fbar_first", # Note: Fbar_first maps to F[ST]
    "Catch_last",
    "SSB_last",
    "Fbar_last"   # Note: Fbar_last maps to F[LT]
  )

  summary_data$Metric <- factor(summary_data$Metric,
                                levels = metric_levels,
                                labels = my_metric_labels)
  # Plot
  plot <- ggplot(summary_data, aes(x = Median, y = Metric, fill = Model, color = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.7) +
    geom_errorbar(aes(xmin = Q1, xmax = Q3), col = "black",
                  position = position_dodge(width = 0.9), width = 0.3, alpha = 0.3) +
    scale_fill_manual(values = my_colors) +
    scale_color_manual(values = my_colors) +
    theme_bw() +
    labs(title = paste("Holistic Model Performance"),
         x = "Score (0-100)",
         y = "Metric",
         fill = "Model",
         color = "Model") +
    scale_y_discrete(labels = my_metric_labels,
                     limits = rev(levels(summary_data$Metric))) +
    theme(
      axis.text = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      plot.title = element_text(size = 18, face = "bold")
    )

  # Print and save
  print(plot)

  # Create the new subfolder if it doesn't exist
  new_sub_dir <- file.path(main.dir, sub.dir, "Holistic_Bar_Plot")

  if (!file.exists(new_sub_dir)){
    dir.create(new_sub_dir)
  }

  ggsave(filename = file.path(new_sub_dir, "Overall_Performance2.png"),
         plot = plot, width = width, height = height, dpi = dpi)
}

plot_model_performance_triangle <- function(mods, is.nsim,
                                            main.dir, sub.dir,
                                            width = 8, height = 7, dpi = 300,
                                            col.opt = "D",
                                            method = NULL,
                                            new_model_names = NULL,
                                            use.n.years.first = 5,
                                            use.n.years.last = 5,
                                            start.years = 1) {
  library(dplyr)
  library(ggtern)
  library(viridisLite)
  
  # Safer warning
  if (missing(start.years)) {
    warning("start.years not provided; defaulting to 1 (first historical year).")
  }
  
  # Function to extract mean Catch, SSB, and Fbar
  process_scores <- function(rep, n_fleets, n_regions, use.n.years, start.years = NULL, type = c("short", "long")) {
    catch_ts <- rowSums(rep$pred_catch)
    ssb_ts <- rowSums(rep$SSB)
    fbar_ts <- rep$Fbar[, ncol(rep$Fbar)] # SAFER: use last column (global Fbar)
    
    if (type == "short") {
      idx <- start.years:(start.years + use.n.years - 1)
    } else {
      idx <- (length(catch_ts) - use.n.years + 1):length(catch_ts)
    }
    
    if(is.null(method)) method = "median"
    if(method == "median") {
      return(c(median(catch_ts[idx]), median(ssb_ts[idx]), median(fbar_ts[idx])))
    }
    if(method == "mean") {
      return(c(mean(catch_ts[idx]), mean(ssb_ts[idx]), mean(fbar_ts[idx])))
    }
  }
  
  results_short <- list()
  results_long <- list()
  
  if (is.nsim) {
    n_models <- length(mods[[1]])
    n_reps <- length(mods)
    
    for (r in seq_len(n_reps)) {
      df_short <- df_long <- data.frame(Model = paste0("Model", seq_len(n_models)))
      for (m in seq_len(n_models)) {
        rep <- mods[[r]][[m]]$om$rep
        input <- mods[[r]][[m]]$om$input$data
        n_fleets <- input$n_fleets[1]
        n_regions <- input$n_regions[1]
        
        short_vals <- process_scores(rep, n_fleets, n_regions, use.n.years.first, start.years, "short")
        long_vals <- process_scores(rep, n_fleets, n_regions, use.n.years.last, NULL, "long")
        
        df_short[m, c("Catch", "SSB", "Fbar")] <- short_vals
        df_long[m, c("Catch", "SSB", "Fbar")] <- long_vals
      }
      
      # Normalize within realization (Catch & SSB up, Fbar inverted)
      for (v in c("Catch", "SSB")) {
        range_v <- max(df_short[[v]]) - min(df_short[[v]])
        df_short[[v]] <- if (range_v == 0) 1 else (df_short[[v]] - min(df_short[[v]])) / range_v
        
        range_v_long <- max(df_long[[v]]) - min(df_long[[v]])
        df_long[[v]] <- if (range_v_long == 0) 1 else (df_long[[v]] - min(df_long[[v]])) / range_v_long
      }
      for (v in c("Fbar")) {
        range_v <- max(df_short[[v]]) - min(df_short[[v]])
        norm_f <- if (range_v == 0) 0 else (df_short[[v]] - min(df_short[[v]])) / range_v
        df_short[[v]] <- 1 - norm_f # inverted
        
        range_v_long <- max(df_long[[v]]) - min(df_long[[v]])
        norm_f_long <- if (range_v_long == 0) 0 else (df_long[[v]] - min(df_long[[v]])) / range_v_long
        df_long[[v]] <- 1 - norm_f_long
      }
      
      # Now we already have values 0–1: just renormalize to sum = 1
      df_short <- df_short %>%
        mutate(total = Catch + SSB + Fbar,
               Catch = ifelse(total == 0, 0, Catch / total),
               SSB = ifelse(total == 0, 0, SSB / total),
               Fbar = ifelse(total == 0, 0, Fbar / total),
               Realization = r)
      
      df_long <- df_long %>%
        mutate(total = Catch + SSB + Fbar,
               Catch = ifelse(total == 0, 0, Catch / total),
               SSB = ifelse(total == 0, 0, SSB / total),
               Fbar = ifelse(total == 0, 0, Fbar / total),
               Realization = r)
      
      results_short[[r]] <- df_short
      results_long[[r]] <- df_long
    }
  } else {
    n_models <- length(mods)
    df_short <- df_long <- data.frame(Model = paste0("Model", seq_len(n_models)))
    
    for (m in seq_len(n_models)) {
      rep <- mods[[m]]$om$rep
      input <- mods[[m]]$om$input$data
      n_fleets <- input$n_fleets[1]
      n_regions <- input$n_regions[1]
      
      short_vals <- process_scores(rep, n_fleets, n_regions, use.n.years.first, start.years, "short")
      long_vals <- process_scores(rep, n_fleets, n_regions, use.n.years.last, NULL, "long")
      
      df_short[m, c("Catch", "SSB", "Fbar")] <- short_vals
      df_long[m, c("Catch", "SSB", "Fbar")] <- long_vals
    }
    
    for (v in c("Catch", "SSB")) {
      range_v <- max(df_short[[v]]) - min(df_short[[v]])
      df_short[[v]] <- if (range_v == 0) 1 else (df_short[[v]] - min(df_short[[v]])) / range_v
      
      range_v_long <- max(df_long[[v]]) - min(df_long[[v]])
      df_long[[v]] <- if (range_v_long == 0) 1 else (df_long[[v]] - min(df_long[[v]])) / range_v_long
    }
    for (v in c("Fbar")) {
      range_v <- max(df_short[[v]]) - min(df_short[[v]])
      norm_f <- if (range_v == 0) 0 else (df_short[[v]] - min(df_short[[v]])) / range_v
      df_short[[v]] <- 1 - norm_f
      
      range_v_long <- max(df_long[[v]]) - min(df_long[[v]])
      norm_f_long <- if (range_v_long == 0) 0 else (df_long[[v]] - min(df_long[[v]])) / range_v_long
      df_long[[v]] <- 1 - norm_f_long
    }
    
    df_short <- df_short %>%
      mutate(total = Catch + SSB + Fbar,
             Catch = ifelse(total == 0, 0, Catch / total),
             SSB = ifelse(total == 0, 0, SSB / total),
             Fbar = ifelse(total == 0, 0, Fbar / total),
             Realization = 1)
    
    df_long <- df_long %>%
      mutate(total = Catch + SSB + Fbar,
             Catch = ifelse(total == 0, 0, Catch / total),
             SSB = ifelse(total == 0, 0, SSB / total),
             Fbar = ifelse(total == 0, 0, Fbar / total),
             Realization = 1)
    
    results_short[[1]] <- df_short
    results_long[[1]] <- df_long
  }
  
  df_short_all <- bind_rows(results_short)
  df_long_all <- bind_rows(results_long)
  
  df_short_all <- bind_rows(results_short) %>%
    group_by(Model) %>%
    summarise(
      Catch = median(Catch),
      SSB = median(SSB),
      Fbar = median(Fbar),
      .groups = "drop"
    )
  
  df_long_all <- bind_rows(results_long) %>%
    group_by(Model) %>%
    summarise(
      Catch = median(Catch),
      SSB = median(SSB),
      Fbar = median(Fbar),
      .groups = "drop"
    )
  
  if (!is.null(new_model_names)) {
    if (length(new_model_names) != length(unique(df_short_all$Model))) {
      stop("Length of new_model_names must match number of models.")
    }
    df_short_all$Model <- factor(df_short_all$Model,
                                 levels = paste0("Model", seq_along(new_model_names)),
                                 labels = new_model_names)
    df_long_all$Model <- factor(df_long_all$Model,
                                levels = paste0("Model", seq_along(new_model_names)),
                                labels = new_model_names)
  }
  
  # === Plot function ===
  plot_and_save <- function(df, title, file) {
    colors <- viridisLite::viridis(n = length(unique(df$Model)), option = col.opt)
    p <- ggtern(df, aes(x = Catch, y = SSB, z = Fbar, color = Model)) +
      geom_point(alpha = 0.8, size = 8) +
      scale_color_manual(values = colors) +
      labs(title = title, T = "Fbar", L = "SSB", R = "Catch") +
      theme_rgbw() +
      theme(plot.title = element_text(hjust = 0.5))
    
    # Create the new subfolder if it doesn't exist
    new_sub_dir <- file.path(main.dir, sub.dir, "Holistic_Triangle_Plot")
    
    if (!file.exists(new_sub_dir)){
      dir.create(new_sub_dir)
    }
    
    # Save the figure inside the new subfolder
    ggsave(filename = file.path(new_sub_dir, file), plot = p, width = width, height = height, dpi = dpi)
    
    print(p)
  }
  
  plot_and_save2 <- function(df, title, file) {
    colors <- viridisLite::viridis(n = length(unique(df$Model)), option = col.opt)
    p <- ggtern(df, aes(x = Catch, y = SSB, z = Fbar, color = Model)) +
      geom_point(size = 1) +
      scale_color_manual(values = colors) +
      labs(title = title, T = "Fbar", L = "SSB", R = "Catch") +
      theme_rgbw() +
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_confidence_tern(breaks = 0.95)
    
    # Create the new subfolder if it doesn't exist
    new_sub_dir <- file.path(main.dir, sub.dir, "Holistic_Triangle_Plot")
    
    if (!file.exists(new_sub_dir)){
      dir.create(new_sub_dir)
    } 
    
    # Save the figure inside the new subfolder
    ggsave(filename = file.path(new_sub_dir, file), plot = p, width = width, height = height, dpi = dpi)
    
    print(p)
  }
  
  plot_and_save(df_short_all, paste0("Short-term Performance (Normalized): Years ", start.years, " to ", start.years + use.n.years.first - 1), "Performance_Short(Norm).png")
  plot_and_save(df_long_all, paste0("Long-term Performance (Normalized): Last ", use.n.years.last, " Years"), "Performance_Long(Norm).png")
  
  results_short2 <- bind_rows(results_short)
  
  results_long2 <- bind_rows(results_long)
  
  if (!is.null(new_model_names)) {
    if (length(new_model_names) != length(unique(results_short2$Model))) {
      stop("Length of new_model_names must match number of models.")
    }
    results_short2$Model <- factor(results_short2$Model,
                                   levels = paste0("Model", seq_along(new_model_names)),
                                   labels = new_model_names)
    results_long2$Model <- factor(results_long2$Model,
                                  levels = paste0("Model", seq_along(new_model_names)),
                                  labels = new_model_names)
  }
  
  plot_and_save2(results_short2, paste0("Short-term Performance: Years ", start.years, " to ", start.years + use.n.years.first - 1), "Performance_Short(Raw).png")
  plot_and_save2(results_long2, paste0("Long-term Performance: Last ", use.n.years.last, " Years"), "Performance_Long(Raw).png")
}

calculate_aacv <- function(catch_values) {
  # Ensure catch_values is a numeric vector
  if (!is.numeric(catch_values)) {
    stop("Input catch_values must be a numeric vector.")
  }
  
  if (length(catch_values) < 2) {
    warning("Need at least 2 years of catch to calculate AACV.")
    return(NA)
  }
  
  catch_diff <- abs(diff(catch_values))
  denom <- sum(catch_values[-length(catch_values)])
  
  if (denom == 0) {
    warning("Total catch is zero; returning NA for AACV.")
    return(NA)
  }
  
  aacv <- sum(catch_diff) / denom
  return(aacv)
}

plot_catch_variation <- function(mods, is.nsim, main.dir, sub.dir, var = "Catch",
                                 width = 10, height = 7, dpi = 300, col.opt = "D",
                                 outlier.opt = NA,
                                 new_model_names = NULL,
                                 base.model = NULL) {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(rlang)
  
  # Helper to calculate AACV (Average Annual Catch Variation)
  calculate_aacv <- function(catch_values) {
    if (!is.numeric(catch_values)) {
      stop("Input catch_values must be a numeric vector.")
    }
    catch_diff <- abs(diff(catch_values))
    aacv <- sum(catch_diff) / sum(catch_values[-length(catch_values)])
    return(aacv)
  }
  
  res <- NULL
  
  if (!is.nsim) {
    n_fleets <- ncol(mods[[1]]$om$rep$pred_catch)
    Years <- mods[[1]]$om$years
    
    res <- lapply(seq_along(mods), function(i) {
      catch_mat <- mods[[i]]$om$rep$pred_catch # n_years × n_fleets
      aacv_list <- lapply(seq_len(n_fleets), function(f) {
        calculate_aacv(catch_mat[, f])
      })
      # Global catch AACV
      global_aacv <- calculate_aacv(rowSums(catch_mat))
      
      tmp <- data.frame(
        Local = t(unlist(aacv_list)),
        Global = global_aacv,
        Model = paste0("Model", i),
        Realization = 1
      )
      colnames(tmp)[1:n_fleets] <- paste0(var, "_Fleet", seq_len(n_fleets))
      colnames(tmp)[n_fleets + 1] <- paste0(var, "_Global")
      tmp
    }) %>% bind_rows()
    
  } else {
    n_fleets <- ncol(mods[[1]][[1]]$om$rep$pred_catch)
    Years <- mods[[1]][[1]]$om$years
    
    res <- lapply(seq_along(mods), function(r) {
      lapply(seq_along(mods[[r]]), function(m) {
        catch_mat <- mods[[r]][[m]]$om$rep$pred_catch
        aacv_list <- lapply(seq_len(n_fleets), function(f) {
          calculate_aacv(catch_mat[, f])
        })
        global_aacv <- calculate_aacv(rowSums(catch_mat))
        
        tmp <- data.frame(
          Local = t(unlist(aacv_list)),
          Global = global_aacv,
          Model = paste0("Model", m),
          Realization = r
        )
        colnames(tmp)[1:n_fleets] <- paste0(var, "_Fleet", seq_len(n_fleets))
        colnames(tmp)[n_fleets + 1] <- paste0(var, "_Global")
        tmp
      }) %>% bind_rows()
    }) %>% bind_rows()
  }
  
  # Rename models if specified
  if (!is.null(new_model_names)) {
    if (length(new_model_names) != length(unique(res$Model))) {
      stop("Length of new_model_names must match the number of models.")
    }
    res$Model <- factor(res$Model,
                        levels = paste0("Model", seq_along(new_model_names)),
                        labels = new_model_names)
    # if (!is.null(base.model)) base.model <- new_model_names[base.model]
    if (!is.null(base.model)) {
      if (!(base.model %in% new_model_names)) {
        warning("base.model does not match any of the new_model_names.")
      }
    }
  }
  
  # Pivot longer
  res <- pivot_longer(res, cols = starts_with(var), names_to = "Label", values_to = "AACV")
  
  # Relative difference if base.model specified
  if (!is.null(base.model)) {
    base_df <- res %>% filter(Model == base.model) %>%
      rename(base_val = AACV) %>%
      select(Realization, Label, base_val)
    
    res <- left_join(res, base_df, by = c("Realization", "Label")) %>%
      mutate(AACV = AACV / base_val - 1)
  }
  
  # Plot
  p1 <- ggplot(res, aes(x = Model, y = AACV, color = Model)) +
    geom_boxplot(lwd = 0.8, outlier.shape = outlier.opt) +
    facet_grid(Label ~ ., scales = "free") +
    scale_color_viridis_d(option = col.opt) +
    ggtitle(paste0("Average Annual Catch Variation",
                   if (!is.null(base.model)) paste0(" (Relative to ", base.model, ")"))) +
    ylab(ifelse(is.null(base.model), "AACV", "Relative AACV Difference")) +
    theme_bw() 
  
  # Save plot
  plot_name <- paste0(var, "_variation", ifelse(is.null(base.model), "", "_Relative"), ".png")
  
  # Create the new subfolder if it doesn't exist
  new_sub_dir <- file.path(main.dir, sub.dir, "Annual_Variation_Boxplot")
  
  if (!file.exists(new_sub_dir)){
    dir.create(new_sub_dir)
  }
  
  ggsave(file.path(main.dir, sub.dir, "Annual_Variation_Boxplot", plot_name), p1, width = width, height = height, dpi = dpi)
  
  return(p1)
}

plot_relative_trajectories <- function(mods, is.nsim,
                                       main.dir, sub.dir,
                                       base.model = NULL,
                                       new_model_names = NULL,
                                       width = 10, height = 7, dpi = 300,
                                       col.opt = "D",
                                       user.Q1 = NULL,
                                       user.Q3 = NULL) {

  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(rlang)

  # ---- Helper: extract relative trajectories ----
  extract_relative_trajectories <- function(mods, base_em_idx) {

    data <- NULL
    Years <- if (!is.nsim) mods[[1]]$om$years else mods[[1]][[1]]$om$years

    if (!is.nsim) {
      n_stocks <- ncol(mods[[1]]$om$rep$SSB)
      n_fleets <- ncol(mods[[1]]$om$rep$pred_catch)

      for (m in seq_along(mods)) {
        ssb_mat <- mods[[m]]$om$rep$SSB
        ssb_base <- mods[[base_em_idx]]$om$rep$SSB
        catch_mat <- mods[[m]]$om$rep$pred_catch
        catch_base <- mods[[base_em_idx]]$om$rep$pred_catch

        # SSB per stock
        for (s in seq_len(n_stocks)) {
          tmp <- ssb_mat[, s] / ssb_base[, s] - 1
          res <- data.frame(Realization = 1, EM = m, Years = Years,
                            Index = paste0("SSB_", s), Value = tmp)
          data <- bind_rows(data, res)
        }

        # SSB global
        tmp <- rowSums(ssb_mat) / rowSums(ssb_base) - 1
        res <- data.frame(Realization = 1, EM = m, Years = Years,
                          Index = "SSB_Global", Value = tmp)
        data <- bind_rows(data, res)

        # Catch per fleet
        for (f in seq_len(n_fleets)) {
          tmp <- catch_mat[, f] / catch_base[, f] - 1
          res <- data.frame(Realization = 1, EM = m, Years = Years,
                            Index = paste0("Catch_", f), Value = tmp)
          data <- bind_rows(data, res)
        }

        # Catch global
        tmp <- rowSums(catch_mat) / rowSums(catch_base) - 1
        res <- data.frame(Realization = 1, EM = m, Years = Years,
                          Index = "Catch_Global", Value = tmp)
        data <- bind_rows(data, res)
      }

    } else {
      n_stocks <- ncol(mods[[1]][[1]]$om$rep$SSB)
      n_fleets <- ncol(mods[[1]][[1]]$om$rep$pred_catch)

      for (r in seq_along(mods)) {
        for (m in seq_along(mods[[r]])) {
          ssb_mat <- mods[[r]][[m]]$om$rep$SSB
          ssb_base <- mods[[r]][[base_em_idx]]$om$rep$SSB
          catch_mat <- mods[[r]][[m]]$om$rep$pred_catch
          catch_base <- mods[[r]][[base_em_idx]]$om$rep$pred_catch

          # SSB per stock
          for (s in seq_len(n_stocks)) {
            tmp <- ssb_mat[, s] / ssb_base[, s] - 1
            res <- data.frame(Realization = r, EM = m, Years = Years,
                              Index = paste0("SSB_", s), Value = tmp)
            data <- bind_rows(data, res)
          }

          # SSB global
          tmp <- rowSums(ssb_mat) / rowSums(ssb_base) - 1
          res <- data.frame(Realization = r, EM = m, Years = Years,
                            Index = "SSB_Global", Value = tmp)
          data <- bind_rows(data, res)

          # Catch per fleet
          for (f in seq_len(n_fleets)) {
            tmp <- catch_mat[, f] / catch_base[, f] - 1
            res <- data.frame(Realization = r, EM = m, Years = Years,
                              Index = paste0("Catch_", f), Value = tmp)
            data <- bind_rows(data, res)
          }

          # Catch global
          tmp <- rowSums(catch_mat) / rowSums(catch_base) - 1
          res <- data.frame(Realization = r, EM = m, Years = Years,
                            Index = "Catch_Global", Value = tmp)
          data <- bind_rows(data, res)
        }
      }
    }

    return(data)
  }

  # --- Determine base_em_idx ---
  n_models <- if (!is.nsim) length(mods) else length(mods[[1]])

  if (is.null(new_model_names)) {
    new_model_names <- paste0("Model", seq_len(n_models))
  }

  if (is.character(base.model)) {
    base_em_idx <- which(new_model_names == base.model)
    if (length(base_em_idx) == 0) stop("base.model not found in new_model_names.")
  } else {
    stop("base.model must be specified as a model name (e.g., 'Model1').")
  }
  
  # --- Extract data ---
  data <- extract_relative_trajectories(mods, base_em_idx = base_em_idx)

  # --- Rename EMs ---
  data$EM <- factor(data$EM,
                    levels = seq_along(new_model_names),
                    labels = new_model_names)

  # --- Summarize ---
  if(is.null(user.Q1)) user.Q1 = 0.4
  if(is.null(user.Q3)) user.Q3 = 0.6
  sum_data <- data %>%
    group_by(EM, Index, Years) %>%
    summarize(
      Median = median(Value, na.rm = TRUE),
      Q1 = quantile(Value, user.Q1, na.rm = TRUE),
      Q3 = quantile(Value, user.Q3, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Create the new subfolder if it doesn't exist
  new_sub_dir <- file.path(main.dir, sub.dir, "Relative_Trajectory")
  
  if (!file.exists(new_sub_dir)){
    dir.create(new_sub_dir)
  }
  
  plot = list()

  t = 0
  # --- Plot loop ---
  for (name in unique(sum_data$Index)) {

    subset_data <- filter(sum_data, Index == name)

    t = t + 1

    plot[[t]] <- ggplot(subset_data, aes(x = Years, y = Median, color = EM)) +
      geom_line(size = 1) +
      geom_ribbon(aes(ymin = Q1, ymax = Q3, fill = EM), alpha = 0.1, color = NA) +
      scale_color_viridis_d(option = col.opt) +
      scale_fill_viridis_d(option = col.opt) +
      labs(x = "Years", y = "Relative Difference", color = "EM", fill = "EM") +
      ggtitle(name) +
      theme_bw() +
      geom_hline(yintercept = 0, col = "red", linetype = "dashed")

    plot_name <- paste0(name,"_trajectories.png")
    
    ggsave(file.path(new_sub_dir, plot_name), plot[[t]], width = width, height = height, dpi = dpi)

    # return(list(plot[[t]]))
    # print(plot[[t]])
  }
  print(plot[[t]])
}


plot_relative_trajectories1 <- function(mods, is.nsim,
                                        main.dir, sub.dir,
                                        base.model = "Model1",
                                        new_model_names = NULL,
                                        width = 10, height = 7, dpi = 300,
                                        col.opt = "D",
                                        user.Q1 = NULL,
                                        user.Q3 = NULL) {
  
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(rlang)
  
  # ---- Helper: extract relative trajectories ----
  extract_relative_trajectories <- function(mods, base_em_idx) {
    
    data <- NULL
    Years <- if (!is.nsim) mods[[1]]$om$years else mods[[1]][[1]]$om$years
    
    if (!is.nsim) {
      n_stocks <- ncol(mods[[1]]$om$rep$SSB)
      
      for (m in seq_along(mods)) {
        ssb_mat <- mods[[m]]$om$rep$SSB
        ssb_base <- mods[[base_em_idx]]$om$rep$SSB
        
        # SSB per stock
        for (s in seq_len(n_stocks)) {
          tmp <- ssb_mat[, s] / ssb_base[, s] - 1
          res <- data.frame(Realization = 1, EM = m, Years = Years,
                            Index = paste0("SSB_", s), Value = tmp)
          data <- bind_rows(data, res)
        }
        
        # SSB global
        tmp <- rowSums(ssb_mat) / rowSums(ssb_base) - 1
        res <- data.frame(Realization = 1, EM = m, Years = Years,
                          Index = "SSB_Global", Value = tmp)
        data <- bind_rows(data, res)
        
      }
      
    } else {
      n_stocks <- ncol(mods[[1]][[1]]$om$rep$SSB)
      n_fleets <- ncol(mods[[1]][[1]]$om$rep$pred_catch)
      
      for (r in seq_along(mods)) {
        for (m in seq_along(mods[[r]])) {
          ssb_mat <- mods[[r]][[m]]$om$rep$SSB
          ssb_base <- mods[[r]][[base_em_idx]]$om$rep$SSB
          
          # SSB per stock
          for (s in seq_len(n_stocks)) {
            tmp <- ssb_mat[, s] / ssb_base[, s] - 1
            res <- data.frame(Realization = r, EM = m, Years = Years,
                              Index = paste0("SSB_", s), Value = tmp)
            data <- bind_rows(data, res)
          }
          
          # SSB global
          tmp <- rowSums(ssb_mat) / rowSums(ssb_base) - 1
          res <- data.frame(Realization = r, EM = m, Years = Years,
                            Index = "SSB_Global", Value = tmp)
          data <- bind_rows(data, res)
          
        }
      }
    }
    
    return(data)
  }
  
  # --- Determine base_em_idx ---
  n_models <- if (!is.nsim) length(mods) else length(mods[[1]])
  
  if (is.null(new_model_names)) {
    new_model_names <- paste0("Model", seq_len(n_models))
  }
  
  if (is.character(base.model)) {
    base_em_idx <- which(new_model_names == base.model)
    if (length(base_em_idx) == 0) stop("base.model not found in new_model_names.")
  } else {
    stop("base.model must be specified as a model name (e.g., 'Model1').")
  }
  
  # --- Extract data ---
  data <- extract_relative_trajectories(mods, base_em_idx = base_em_idx)
  
  # --- Rename EMs ---
  data$EM <- factor(data$EM,
                    levels = seq_along(new_model_names),
                    labels = new_model_names)
  
  # --- Summarize ---
  if(is.null(user.Q1)) user.Q1 = 0.4
  if(is.null(user.Q3)) user.Q3 = 0.6
  sum_data <- data %>%
    group_by(EM, Index, Years) %>%
    summarize(
      Median = median(Value, na.rm = TRUE),
      Q1 = quantile(Value, user.Q1, na.rm = TRUE),
      Q3 = quantile(Value, user.Q3, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Create the new subfolder if it doesn't exist
  new_sub_dir <- file.path(main.dir, sub.dir, "Relative_Trajectory")
  
  if (!file.exists(new_sub_dir)){
    dir.create(new_sub_dir)
  }
  
  # --- Plot loop ---
  for (name in unique(sum_data$Index)) {
    
    subset_data <- filter(sum_data, Index == name)
    
    plot1 <- ggplot(subset_data, aes(x = Years, y = Median, color = EM)) +
      geom_line(size = 1) +
      geom_ribbon(aes(ymin = Q1, ymax = Q3, fill = EM), alpha = 0.1, color = NA) +
      scale_color_viridis_d(option = col.opt) +
      scale_fill_viridis_d(option = col.opt) +
      labs(x = "Years", y = "Relative Difference", color = "EM", fill = "EM") +
      ggtitle(name) +
      theme_bw() +
      geom_hline(yintercept = 0, col = "red", linetype = "dashed")
    
    plot_name <- paste0(name,"_trajectories.png")
    ggsave(file.path(new_sub_dir, plot_name), plot1, width = width, height = height, dpi = dpi)
  }
  return(plot1)
}

plot_relative_trajectories2 <- function(mods, is.nsim,
                                        main.dir, sub.dir,
                                        base.model = "Model1",
                                        new_model_names = NULL,
                                        width = 10, height = 7, dpi = 300,
                                        col.opt = "D",
                                        user.Q1 = NULL,
                                        user.Q3 = NULL) {
  
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(rlang)
  
  # ---- Helper: extract relative trajectories ----
  extract_relative_trajectories <- function(mods, base_em_idx) {
    
    data <- NULL
    Years <- if (!is.nsim) mods[[1]]$om$years else mods[[1]][[1]]$om$years
    
    if (!is.nsim) {
      n_stocks <- ncol(mods[[1]]$om$rep$SSB)
      n_fleets <- ncol(mods[[1]]$om$rep$pred_catch)
      
      for (m in seq_along(mods)) {
        catch_mat <- mods[[m]]$om$rep$pred_catch
        catch_base <- mods[[base_em_idx]]$om$rep$pred_catch
        
        # Catch per fleet
        for (f in seq_len(n_fleets)) {
          tmp <- catch_mat[, f] / catch_base[, f] - 1
          res <- data.frame(Realization = 1, EM = m, Years = Years,
                            Index = paste0("Catch_", f), Value = tmp)
          data <- bind_rows(data, res)
        }
        
        # Catch global
        tmp <- rowSums(catch_mat) / rowSums(catch_base) - 1
        res <- data.frame(Realization = 1, EM = m, Years = Years,
                          Index = "Catch_Global", Value = tmp)
        data <- bind_rows(data, res)
      }
      
    } else {
      n_stocks <- ncol(mods[[1]][[1]]$om$rep$SSB)
      n_fleets <- ncol(mods[[1]][[1]]$om$rep$pred_catch)
      
      for (r in seq_along(mods)) {
        for (m in seq_along(mods[[r]])) {
          catch_mat <- mods[[r]][[m]]$om$rep$pred_catch
          catch_base <- mods[[r]][[base_em_idx]]$om$rep$pred_catch
          
          # Catch per fleet
          for (f in seq_len(n_fleets)) {
            tmp <- catch_mat[, f] / catch_base[, f] - 1
            res <- data.frame(Realization = r, EM = m, Years = Years,
                              Index = paste0("Catch_", f), Value = tmp)
            data <- bind_rows(data, res)
          }
          
          # Catch global
          tmp <- rowSums(catch_mat) / rowSums(catch_base) - 1
          res <- data.frame(Realization = r, EM = m, Years = Years,
                            Index = "Catch_Global", Value = tmp)
          data <- bind_rows(data, res)
        }
      }
    }
    
    return(data)
  }
  
  # --- Determine base_em_idx ---
  n_models <- if (!is.nsim) length(mods) else length(mods[[1]])
  
  if (is.null(new_model_names)) {
    new_model_names <- paste0("Model", seq_len(n_models))
  }
  
  if (is.character(base.model)) {
    base_em_idx <- which(new_model_names == base.model)
    if (length(base_em_idx) == 0) stop("base.model not found in new_model_names.")
  } else {
    stop("base.model must be specified as a model name (e.g., 'Model1').")
  }
  
  # --- Extract data ---
  data <- extract_relative_trajectories(mods, base_em_idx = base_em_idx)
  
  # --- Rename EMs ---
  data$EM <- factor(data$EM,
                    levels = seq_along(new_model_names),
                    labels = new_model_names)
  
  # --- Summarize ---
  if(is.null(user.Q1)) user.Q1 = 0.4
  if(is.null(user.Q3)) user.Q3 = 0.6
  sum_data <- data %>%
    group_by(EM, Index, Years) %>%
    summarize(
      Median = median(Value, na.rm = TRUE),
      Q1 = quantile(Value, user.Q1, na.rm = TRUE),
      Q3 = quantile(Value, user.Q3, na.rm = TRUE),
      .groups = "drop"
    )
  
  # Create the new subfolder if it doesn't exist
  new_sub_dir <- file.path(main.dir, sub.dir, "Relative_Trajectory")
  
  if (!file.exists(new_sub_dir)){
    dir.create(new_sub_dir)
  }
  
  # --- Plot loop ---
  for (name in unique(sum_data$Index)) {
    
    subset_data <- filter(sum_data, Index == name)
    
    plot2 <- ggplot(subset_data, aes(x = Years, y = Median, color = EM)) +
      geom_line(size = 1) +
      geom_ribbon(aes(ymin = Q1, ymax = Q3, fill = EM), alpha = 0.1, color = NA) +
      scale_color_viridis_d(option = col.opt) +
      scale_fill_viridis_d(option = col.opt) +
      labs(x = "Years", y = "Relative Difference", color = "EM", fill = "EM") +
      ggtitle(name) +
      theme_bw() +
      geom_hline(yintercept = 0, col = "red", linetype = "dashed")
    
    plot_name <- paste0(name,"_trajectories.png")
    ggsave(file.path(new_sub_dir, plot_name), plot2, width = width, height = height, dpi = dpi)
    
  }
  return(plot2)
}

plot_ssb_variation <- function(mods, is.nsim, main.dir, sub.dir, var = "SSB",
                               width = 10, height = 7, dpi = 300, col.opt = "D",
                               outlier.opt = NA,
                               new_model_names = NULL,
                               base.model = NULL) {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(rlang)
  
  # Helper to calculate AACV (Average Annual Catch Variation)
  calculate_aacv <- function(catch_values) {
    if (!is.numeric(catch_values)) {
      stop("Input catch_values must be a numeric vector.")
    }
    catch_diff <- abs(diff(catch_values))
    aacv <- sum(catch_diff) / sum(catch_values[-length(catch_values)])
    return(aacv)
  }
  
  res <- NULL
  
  if (!is.nsim) {
    n_regions <- ncol(mods[[1]]$om$rep$SSB)
    Years <- mods[[1]]$om$years
    
    res <- lapply(seq_along(mods), function(i) {
      catch_mat <- mods[[i]]$om$rep$SSB # n_years × n_regions
      aacv_list <- lapply(seq_len(n_regions), function(f) {
        calculate_aacv(catch_mat[, f])
      })
      # Global catch AACV
      global_aacv <- calculate_aacv(rowSums(catch_mat))
      
      tmp <- data.frame(
        Local = t(unlist(aacv_list)),
        Global = global_aacv,
        Model = paste0("Model", i),
        Realization = 1
      )
      colnames(tmp)[1:n_regions] <- paste0(var, "_Region", seq_len(n_regions))
      colnames(tmp)[n_regions + 1] <- paste0(var, "_Global")
      tmp
    }) %>% bind_rows()
    
  } else {
    n_regions <- ncol(mods[[1]][[1]]$om$rep$SSB)
    Years <- mods[[1]][[1]]$om$years
    
    res <- lapply(seq_along(mods), function(r) {
      lapply(seq_along(mods[[r]]), function(m) {
        catch_mat <- mods[[r]][[m]]$om$rep$SSB
        aacv_list <- lapply(seq_len(n_regions), function(f) {
          calculate_aacv(catch_mat[, f])
        })
        global_aacv <- calculate_aacv(rowSums(catch_mat))
        
        tmp <- data.frame(
          Local = t(unlist(aacv_list)),
          Global = global_aacv,
          Model = paste0("Model", m),
          Realization = r
        )
        colnames(tmp)[1:n_regions] <- paste0(var, "_Region", seq_len(n_regions))
        colnames(tmp)[n_regions + 1] <- paste0(var, "_Global")
        tmp
      }) %>% bind_rows()
    }) %>% bind_rows()
  }
  
  # Rename models if specified
  if (!is.null(new_model_names)) {
    if (length(new_model_names) != length(unique(res$Model))) {
      stop("Length of new_model_names must match the number of models.")
    }
    res$Model <- factor(res$Model,
                        levels = paste0("Model", seq_along(new_model_names)),
                        labels = new_model_names)
    # if (!is.null(base.model)) base.model <- new_model_names[base.model]
    if (!is.null(base.model)) {
      if (!(base.model %in% new_model_names)) {
        warning("base.model does not match any of the new_model_names.")
      }
    }
  }
  
  # Pivot longer
  res <- pivot_longer(res, cols = starts_with(var), names_to = "Label", values_to = "AACV")
  
  # Ensure Global is always the last row in facet
  labels_order <- res %>%
    distinct(Label) %>%
    arrange(if_else(grepl("Global", Label), 2, 1), Label) %>%
    pull(Label)
  
  res$Label <- factor(res$Label, levels = labels_order)
  
  # Relative difference if base.model specified
  if (!is.null(base.model)) {
    base_df <- res %>% filter(Model == base.model) %>%
      rename(base_val = AACV) %>%
      select(Realization, Label, base_val)
    
    res <- left_join(res, base_df, by = c("Realization", "Label")) %>%
      mutate(AACV = AACV / base_val - 1)
  }
  
  # Plot
  p1 <- ggplot(res, aes(x = Model, y = AACV, color = Model)) +
    geom_boxplot(lwd = 0.8, outlier.shape = outlier.opt) +
    facet_grid(Label ~ ., scales = "free") +
    scale_color_viridis_d(option = col.opt) +
    ggtitle(paste0("Average Annual SSB Variation",
                   if (!is.null(base.model)) paste0(" (Relative to ", base.model, ")"))) +
    ylab(ifelse(is.null(base.model), "AASV", "Relative AASV Difference")) +
    theme_bw() 
  
  # Save plot
  plot_name <- paste0(var, "_variation", ifelse(is.null(base.model), "", "_Relative"), ".png")
  
  # Create the new subfolder if it doesn't exist
  new_sub_dir <- file.path(main.dir, sub.dir, "Annual_Variation_Boxplot")
  
  if (!file.exists(new_sub_dir)){
    dir.create(new_sub_dir)
  }
  
  ggsave(file.path(main.dir, sub.dir, "Annual_Variation_Boxplot", plot_name), p1, width = width, height = height, dpi = dpi)
  
  return(p1)
}


plot_fbar_variation <- function(mods, is.nsim, main.dir, sub.dir, var = "Fbar",
                                width = 10, height = 7, dpi = 300, col.opt = "D",
                                outlier.opt = NA,
                                new_model_names = NULL,
                                base.model = NULL) {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(rlang)
  
  calculate_aacv <- function(values) {
    if (!is.numeric(values)) stop("Input must be numeric.")
    diff_abs <- abs(diff(values))
    aacv <- sum(diff_abs) / sum(values[-length(values)])
    return(aacv)
  }
  
  res <- NULL
  
  if (!is.nsim) {
    n_fleets <- mods[[1]]$om$input$data$n_fleets[1]
    n_regions <- mods[[1]]$om$input$data$n_regions[1]
    Years <- mods[[1]]$om$years
    
    res <- lapply(seq_along(mods), function(i) {
      fbar_mat <- mods[[i]]$om$rep$Fbar # Fbar matrix
      
      # Fleets
      aacv_fleet <- lapply(seq_len(n_fleets), function(f) {
        calculate_aacv(fbar_mat[, f])
      })
      
      # Regions
      aacv_region <- lapply(seq_len(n_regions), function(r) {
        calculate_aacv(fbar_mat[, n_fleets + r])
      })
      
      # Global
      aacv_global <- calculate_aacv(fbar_mat[, n_fleets + n_regions + 1])
      
      tmp <- data.frame(
        t(unlist(aacv_fleet)),
        t(unlist(aacv_region)),
        Global = aacv_global,
        Model = paste0("Model", i),
        Realization = 1
      )
      colnames(tmp)[1:n_fleets] <- paste0(var, "_Fleet", seq_len(n_fleets))
      colnames(tmp)[(n_fleets + 1):(n_fleets + n_regions)] <- paste0(var, "_Region", seq_len(n_regions))
      colnames(tmp)[n_fleets + n_regions + 1] <- paste0(var, "_Global")
      tmp
    }) %>% bind_rows()
    
  } else {
    n_fleets <- mods[[1]][[1]]$om$input$data$n_fleets[1]
    n_regions <- mods[[1]][[1]]$om$input$data$n_regions[1]
    Years <- mods[[1]][[1]]$om$years
    
    res <- lapply(seq_along(mods), function(r) {
      lapply(seq_along(mods[[r]]), function(m) {
        fbar_mat <- mods[[r]][[m]]$om$rep$Fbar
        
        # Fleets
        aacv_fleet <- lapply(seq_len(n_fleets), function(f) {
          calculate_aacv(fbar_mat[, f])
        })
        
        # Regions
        aacv_region <- lapply(seq_len(n_regions), function(rr) {
          calculate_aacv(fbar_mat[, n_fleets + rr])
        })
        
        # Global
        aacv_global <- calculate_aacv(fbar_mat[, n_fleets + n_regions + 1])
        
        tmp <- data.frame(
          t(unlist(aacv_fleet)),
          t(unlist(aacv_region)),
          Global = aacv_global,
          Model = paste0("Model", m),
          Realization = r
        )
        colnames(tmp)[1:n_fleets] <- paste0(var, "_Fleet", seq_len(n_fleets))
        colnames(tmp)[(n_fleets + 1):(n_fleets + n_regions)] <- paste0(var, "_Region", seq_len(n_regions))
        colnames(tmp)[n_fleets + n_regions + 1] <- paste0(var, "_Global")
        tmp
      }) %>% bind_rows()
    }) %>% bind_rows()
  }
  
  # Rename models if needed
  if (!is.null(new_model_names)) {
    if (length(new_model_names) != length(unique(res$Model))) {
      stop("Length of new_model_names must match the number of models.")
    }
    res$Model <- factor(res$Model,
                        levels = paste0("Model", seq_along(new_model_names)),
                        labels = new_model_names)
    if (!is.null(base.model)) {
      if (!(base.model %in% new_model_names)) {
        warning("base.model does not match any of the new_model_names.")
      }
    }
  }
  
  res <- pivot_longer(res, cols = starts_with(var), names_to = "Label", values_to = "AACV")
  
  # Ensure Global is always the last row in facet
  labels_order <- res %>%
    distinct(Label) %>%
    arrange(if_else(grepl("Global", Label), 2, 1), Label) %>%
    pull(Label)
  
  res$Label <- factor(res$Label, levels = labels_order)
  
  
  if (!is.null(base.model)) {
    base_df <- res %>% filter(Model == base.model) %>%
      rename(base_val = AACV) %>%
      select(Realization, Label, base_val)
    
    res <- left_join(res, base_df, by = c("Realization", "Label")) %>%
      mutate(AACV = AACV / base_val - 1)
  }
  
  p1 <- ggplot(res, aes(x = Model, y = AACV, color = Model)) +
    geom_boxplot(lwd = 0.8, outlier.shape = outlier.opt) +
    facet_grid(Label ~ ., scales = "free") +
    scale_color_viridis_d(option = col.opt) +
    ggtitle(paste0("Average Annual Fbar Variation",
                   if (!is.null(base.model)) paste0(" (Relative to ", base.model, ")"))) +
    ylab(ifelse(is.null(base.model), "AAFV", "Relative AAFV Difference")) +
    theme_bw()
  
  # Save
  out_dir <- file.path(main.dir, sub.dir)
  if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)
  
  # Create the new subfolder if it doesn't exist
  new_sub_dir <- file.path(main.dir, sub.dir, "Annual_Variation_Boxplot")
  
  if (!file.exists(new_sub_dir)){
    dir.create(new_sub_dir)
  }
  
  plot_name <- paste0(var, "_variation", ifelse(is.null(base.model), "", "_Relative"), ".png")
  ggsave(file.path(out_dir, "Annual_Variation_Boxplot", plot_name), p1, width = width, height = height, dpi = dpi)
  
  return(p1)
}


plot_AAV_performance <- function(mods, is.nsim,
                                 main.dir, sub.dir,
                                 width = 8, height = 7, dpi = 300,
                                 new_model_names = NULL,
                                 col.opt = "D") {
  library(dplyr)
  library(ggtern)
  library(viridisLite)
  
  # Helper: calculate Average Annual Catch Variation (AACV)
  calculate_aacv <- function(values) {
    if (!is.numeric(values)) {
      stop("Input must be a numeric vector.")
    }
    diffs <- abs(diff(values))
    aacv <- sum(diffs) / sum(values[-length(values)])
    return(aacv)
  }
  
  # --- Loop through each mod and collect variation measures ---
  res_list <- list()
  
  if (!is.nsim) {
    n_fleets <- ncol(mods[[1]]$om$rep$pred_catch)
    n_regions <- ncol(mods[[1]]$om$rep$SSB)
    
    for (i in seq_along(mods)) {
      om_rep <- mods[[i]]$om$rep
      
      # Catch variation
      catch_aacv <- calculate_aacv(rowSums(om_rep$pred_catch))
      
      # SSB variation
      ssb_aacv <- calculate_aacv(rowSums(om_rep$SSB))
      
      # Fbar variation (global is last col)
      fbar_aacv <- calculate_aacv(om_rep$Fbar[, ncol(om_rep$Fbar)])
      
      res_list[[i]] <- data.frame(
        Model = paste0("Model", i),
        Catch = catch_aacv,
        SSB = ssb_aacv,
        Fbar = fbar_aacv
      )
    }
    
  } else {
    n_fleets <- ncol(mods[[1]][[1]]$om$rep$pred_catch)
    n_regions <- ncol(mods[[1]][[1]]$om$rep$SSB)
    
    for (r in seq_along(mods)) {
      for (m in seq_along(mods[[r]])) {
        om_rep <- mods[[r]][[m]]$om$rep
        
        # Catch variation
        catch_aacv <- calculate_aacv(rowSums(om_rep$pred_catch))
        
        # SSB variation
        ssb_aacv <- calculate_aacv(rowSums(om_rep$SSB))
        
        # Fbar variation (global is last col)
        fbar_aacv <- calculate_aacv(om_rep$Fbar[, ncol(om_rep$Fbar)])
        
        res_list[[length(res_list) + 1]] <- data.frame(
          Model = paste0("Model", m),
          Realization = r,
          Catch = catch_aacv,
          SSB = ssb_aacv,
          Fbar = fbar_aacv
        )
      }
    }
  }
  
  res_all <- bind_rows(res_list)
  
  # --- Take median across realizations if needed ---
  df_all <- res_all %>%
    group_by(Model) %>%
    summarise(
      Catch = median(Catch),
      SSB = median(SSB),
      Fbar = median(Fbar),
      .groups = "drop"
    )
  
  # Rename models if user provided new names
  if (!is.null(new_model_names)) {
    if (length(new_model_names) != length(unique(df_all$Model))) {
      stop("Length of new_model_names must match the number of models.")
    }
    df_all$Model <- factor(df_all$Model,
                           levels = paste0("Model", seq_along(new_model_names)),
                           labels = new_model_names)
  }
  
  # --- Normalize & invert scores (higher variation is BAD) ---
  for (v in c("Catch", "SSB", "Fbar")) {
    range_v <- max(df_all[[v]]) - min(df_all[[v]])
    if (range_v == 0) {
      df_all[[paste0(v, "_score")]] <- 100 # if no variation, give full score
    } else {
      norm_v <- (df_all[[v]] - min(df_all[[v]])) / range_v
      df_all[[paste0(v, "_score")]] <- 100 * (1 - norm_v)
    }
  }
  
  # --- Plot triangle using inverted scores ---
  colors <- viridisLite::viridis(n = length(unique(df_all$Model)), option = col.opt)
  
  p <- ggtern(df_all, aes(x = Catch_score, y = SSB_score, z = Fbar_score, color = Model)) +
    geom_point(alpha = 0.8, size = 8) +
    scale_color_manual(values = colors) +
    labs(title = "Average Annual Variation in Catch, SSB, and F\nHigher Score = Lower Variation",
         T = "Catch",
         L = "SSB",
         R = "Fbar") +
    theme_rgbw() +
    theme(plot.title = element_text(hjust = 0.5))
  
  # Create the new subfolder if it doesn't exist
  new_sub_dir <- file.path(main.dir, sub.dir, "Holistic_Annual_Variation")
  
  if (!file.exists(new_sub_dir)){
    dir.create(new_sub_dir)
  }
  
  # Save & print
  ggsave(filename = file.path(main.dir, sub.dir, "Holistic_Annual_Variation", "Overall_AAV_Median.png"),
         plot = p, width = width, height = height, dpi = dpi)
  
  print(p)
  
  # --- Normalize & invert scores (higher variation is BAD) ---
  for (v in c("Catch", "SSB", "Fbar")) {
    range_v <- max(res_all[[v]]) - min(res_all[[v]])
    if (range_v == 0) {
      res_all[[paste0(v, "_score")]] <- 100 # if no variation, give full score
    } else {
      norm_v <- (res_all[[v]] - min(res_all[[v]])) / range_v
      res_all[[paste0(v, "_score")]] <- 100 * (1 - norm_v)
    }
  }
  
  # Rename models if user provided new names
  if (!is.null(new_model_names)) {
    if (length(new_model_names) != length(unique(res_all$Model))) {
      stop("Length of new_model_names must match the number of models.")
    }
    res_all$Model <- factor(res_all$Model,
                           levels = paste0("Model", seq_along(new_model_names)),
                           labels = new_model_names)
  }
  
  p2 <- ggtern(res_all, aes(x = Catch_score, y = SSB_score, z = Fbar_score, color = Model)) +
    geom_point(size = 1) +
    scale_color_manual(values = colors) +
    labs(title = "Average Annual Variation in Catch, SSB, and F\nHigher Score = Lower Variation",
         T = "Catch",
         L = "SSB",
         R = "Fbar") +
    theme_rgbw() +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_confidence_tern(breaks = 0.95)
  
  print(p2)
  
  # Save & print
  ggsave(filename = file.path(main.dir, sub.dir, "Holistic_Annual_Variation", "Overall_AAV_Distribution.png"),
         plot = p2, width = width, height = height, dpi = dpi)
  
  # Raw values
  res_long <- res_all %>%
    pivot_longer(cols = c(Catch, SSB, Fbar),
                 names_to = "Variable", values_to = "Score") %>%
    mutate(Variable = gsub("_score", "", Variable))  # clean name
  
  # if (!is.null(new_model_names)) {
  #   if (length(new_model_names) != length(unique(res_long$Model))) {
  #     stop("Length of new_model_names must match the number of models.")
  #   }
  #   res_long$Model <- factor(res_long$Model,
  #                            levels = paste0("Model", seq_along(new_model_names)),
  #                            labels = new_model_names)
  # }
  
  p3 <- ggplot(res_long, aes(x = Score, y = Model, fill = Model)) +
    ggridges::geom_density_ridges(alpha = 0.8, scale = 1.0, 
                        quantile_lines = TRUE, 
                        quantiles = 0.5,
                        linetype = 2) +
    theme_bw() +
    scale_fill_viridis_d(option = col.opt) +
    labs(
      title = "Average Annual Variation (Raw) Distributions (Lower = Better)",
      x = "Value",
      y = "Model"
    ) +
    facet_wrap(~ Variable)
  
  print(p3)
  
  # Save & print
  ggsave(filename = file.path(main.dir, sub.dir, "Holistic_Annual_Variation", "Overall_AAV_Density(Raw).png"),
         plot = p3, width = width, height = height, dpi = dpi)
  
  # Normalized
  res_norm <- res_all
  
  # Normalize and invert across all realizations
  for (v in c("Catch", "SSB", "Fbar")) {
    range_v <- max(res_norm[[v]]) - min(res_norm[[v]])
    if (range_v == 0) {
      res_norm[[paste0(v, "_score")]] <- 100
    } else {
      norm_v <- (res_norm[[v]] - min(res_norm[[v]])) / range_v
      res_norm[[paste0(v, "_score")]] <- 100 * (1 - norm_v)
    }
  }
  
  library(tidyr)
  library(ggridges)
  library(ggplot2)
  
  # Reshape normalized data
  res_long <- res_norm %>%
    pivot_longer(cols = c(Catch_score, SSB_score, Fbar_score),
                 names_to = "Variable", values_to = "Score") %>%
    mutate(Variable = gsub("_score", "", Variable))  # clean name
  
  # Optional relabeling
  # if (!is.null(new_model_names)) {
  #   if (length(new_model_names) != length(unique(res_long$Model))) {
  #     stop("Length of new_model_names must match the number of models.")
  #   }
  #   res_long$Model <- factor(res_long$Model,
  #                            levels = paste0("Model", seq_along(new_model_names)),
  #                            labels = new_model_names)
  # }
  
  # Plot ridge
  # ggplot(res_long, aes(x = Score, y = Variable, fill = Variable)) +
  #   geom_density_ridges(alpha = 0.8, scale = 1.0) +
  #   theme_bw() +
  #   scale_fill_viridis_d(option = "plasma") +
  #   labs(
  #     title = "Performance Score Distributions (Higher = Better)",
  #     x = "Score (Normalized, Inverted)",
  #     y = "Metric"
  #   ) +
  #   facet_wrap(~ Model)
  
  p4 <- ggplot(res_long, aes(x = Score, y = Model, fill = Model)) +
    ggridges::geom_density_ridges(alpha = 0.8, scale = 1.0, 
                        quantile_lines = TRUE, 
                        quantiles = 0.5,
                        linetype = 2) +
    # geom_vline(data = medians, aes(xintercept = Median, color = Variable), linetype = "dashed",
    #            linewidth = 0.8, show.legend = FALSE) + 
    theme_bw() +
    scale_fill_viridis_d(option = col.opt) +
    labs(
      title = "Average Annual Variation (Normalized) Distributions (Higher = Better)",
      x = "Value",
      y = "Model"
    ) +
    facet_wrap(~ Variable)
  
  print(p4)
  # Save & print
  ggsave(filename = file.path(main.dir, sub.dir, "Holistic_Annual_Variation", "Overall_AAV_Density(Norm).png"),
         plot = p4, width = width, height = height, dpi = dpi)
}


# plot_status_triangle <- function(mods, is.nsim,
#                                  main.dir, sub.dir,
#                                  width = 8, height = 7, dpi = 300,
#                                  col.opt = "D",
#                                  method = NULL,
#                                  new_model_names = NULL,
#                                  use.n.years.first = 5,
#                                  use.n.years.last = 5,
#                                  start.years = 1) {
#   library(dplyr)
#   library(ggtern)
#   library(viridisLite)
#   
#   # Safer warning
#   if (missing(start.years)) {
#     warning("start.years not provided; defaulting to 1 (first historical year).")
#   }
#   
#   # Function to extract mean Catch, SSB, and Fbar
#   process_scores <- function(rep, n_fleets, n_regions, use.n.years, start.years = NULL, type = c("short", "long")) {
#     catch_ts <- rowSums(rep$pred_catch)
#     catch_ts_brp <- exp(rep$log_Y_FXSPR[, ncol(rep$Fbar)])
#     catch_ts <- catch_ts/catch_ts_brp
#     
#     ssb_ts <- rowSums(rep$SSB)
#     ssb_ts_brp <- exp(rep$log_SSB_FXSPR[,n_regions+1])
#     ssb_ts <- ssb_ts/ssb_ts_brp
#     
#     fbar_ts <- rep$Fbar[, ncol(rep$Fbar)] # SAFER: use last column (global Fbar)
#     fbar_ts_brp <- exp(rep$log_Fbar_XSPR[, ncol(rep$Fbar)])
#     fbar_ts <- fbar_ts/fbar_ts_brp
#     
#     if (type == "short") {
#       idx <- start.years:(start.years + use.n.years - 1)
#     } else {
#       idx <- (length(catch_ts) - use.n.years + 1):length(catch_ts)
#     }
#     
#     if(is.null(method)) method = "median"
#     if(method == "median") {
#       return(c(median(catch_ts[idx]), median(ssb_ts[idx]), median(fbar_ts[idx])))
#     }
#     if(method == "mean") {
#       return(c(mean(catch_ts[idx]), mean(ssb_ts[idx]), mean(fbar_ts[idx])))
#     }
#   }
#   
#   results_short <- list()
#   results_long <- list()
#   
#   if (is.nsim) {
#     n_models <- length(mods[[1]])
#     n_reps <- length(mods)
#     
#     for (r in seq_len(n_reps)) {
#       df_short <- df_long <- data.frame(Model = paste0("Model", seq_len(n_models)))
#       for (m in seq_len(n_models)) {
#         rep <- mods[[r]][[m]]$om$rep
#         input <- mods[[r]][[m]]$om$input$data
#         n_fleets <- input$n_fleets[1]
#         n_regions <- input$n_regions[1]
#         
#         short_vals <- process_scores(rep, n_fleets, n_regions, use.n.years.first, start.years, "short")
#         long_vals <- process_scores(rep, n_fleets, n_regions, use.n.years.last, NULL, "long")
#         
#         df_short[m, c("Catch", "SSB", "Fbar")] <- short_vals
#         df_long[m, c("Catch", "SSB", "Fbar")] <- long_vals
#       }
#       
#       # Normalize within realization (Catch & SSB up, Fbar inverted)
#       for (v in c("Catch", "SSB")) {
#         range_v <- max(df_short[[v]]) - min(df_short[[v]])
#         df_short[[v]] <- if (range_v == 0) 1 else (df_short[[v]] - min(df_short[[v]])) / range_v
#         
#         range_v_long <- max(df_long[[v]]) - min(df_long[[v]])
#         df_long[[v]] <- if (range_v_long == 0) 1 else (df_long[[v]] - min(df_long[[v]])) / range_v_long
#       }
#       for (v in c("Fbar")) {
#         range_v <- max(df_short[[v]]) - min(df_short[[v]])
#         norm_f <- if (range_v == 0) 0 else (df_short[[v]] - min(df_short[[v]])) / range_v
#         df_short[[v]] <- 1 - norm_f # inverted
#         
#         range_v_long <- max(df_long[[v]]) - min(df_long[[v]])
#         norm_f_long <- if (range_v_long == 0) 0 else (df_long[[v]] - min(df_long[[v]])) / range_v_long
#         df_long[[v]] <- 1 - norm_f_long
#       }
#       
#       # Now we already have values 0–1: just renormalize to sum = 1
#       df_short <- df_short %>%
#         mutate(total = Catch + SSB + Fbar,
#                Catch = ifelse(total == 0, 0, Catch / total),
#                SSB = ifelse(total == 0, 0, SSB / total),
#                Fbar = ifelse(total == 0, 0, Fbar / total),
#                Realization = r)
#       
#       df_long <- df_long %>%
#         mutate(total = Catch + SSB + Fbar,
#                Catch = ifelse(total == 0, 0, Catch / total),
#                SSB = ifelse(total == 0, 0, SSB / total),
#                Fbar = ifelse(total == 0, 0, Fbar / total),
#                Realization = r)
#       
#       results_short[[r]] <- df_short
#       results_long[[r]] <- df_long
#     }
#   } else {
#     n_models <- length(mods)
#     df_short <- df_long <- data.frame(Model = paste0("Model", seq_len(n_models)))
#     
#     for (m in seq_len(n_models)) {
#       rep <- mods[[m]]$om$rep
#       input <- mods[[m]]$om$input$data
#       n_fleets <- input$n_fleets[1]
#       n_regions <- input$n_regions[1]
#       
#       short_vals <- process_scores(rep, n_fleets, n_regions, use.n.years.first, start.years, "short")
#       long_vals <- process_scores(rep, n_fleets, n_regions, use.n.years.last, NULL, "long")
#       
#       df_short[m, c("Catch", "SSB", "Fbar")] <- short_vals
#       df_long[m, c("Catch", "SSB", "Fbar")] <- long_vals
#     }
#     
#     for (v in c("Catch", "SSB")) {
#       range_v <- max(df_short[[v]]) - min(df_short[[v]])
#       df_short[[v]] <- if (range_v == 0) 1 else (df_short[[v]] - min(df_short[[v]])) / range_v
#       
#       range_v_long <- max(df_long[[v]]) - min(df_long[[v]])
#       df_long[[v]] <- if (range_v_long == 0) 1 else (df_long[[v]] - min(df_long[[v]])) / range_v_long
#     }
#     for (v in c("Fbar")) {
#       range_v <- max(df_short[[v]]) - min(df_short[[v]])
#       norm_f <- if (range_v == 0) 0 else (df_short[[v]] - min(df_short[[v]])) / range_v
#       df_short[[v]] <- 1 - norm_f
#       
#       range_v_long <- max(df_long[[v]]) - min(df_long[[v]])
#       norm_f_long <- if (range_v_long == 0) 0 else (df_long[[v]] - min(df_long[[v]])) / range_v_long
#       df_long[[v]] <- 1 - norm_f_long
#     }
#     
#     df_short <- df_short %>%
#       mutate(total = Catch + SSB + Fbar,
#              Catch = ifelse(total == 0, 0, Catch / total),
#              SSB = ifelse(total == 0, 0, SSB / total),
#              Fbar = ifelse(total == 0, 0, Fbar / total),
#              Realization = 1)
#     
#     df_long <- df_long %>%
#       mutate(total = Catch + SSB + Fbar,
#              Catch = ifelse(total == 0, 0, Catch / total),
#              SSB = ifelse(total == 0, 0, SSB / total),
#              Fbar = ifelse(total == 0, 0, Fbar / total),
#              Realization = 1)
#     
#     results_short[[1]] <- df_short
#     results_long[[1]] <- df_long
#   }
#   
#   df_short_all <- bind_rows(results_short)
#   df_long_all <- bind_rows(results_long)
#   
#   df_short_all <- bind_rows(results_short) %>%
#     group_by(Model) %>%
#     summarise(
#       Catch = median(Catch),
#       SSB = median(SSB),
#       Fbar = median(Fbar),
#       .groups = "drop"
#     )
#   
#   df_long_all <- bind_rows(results_long) %>%
#     group_by(Model) %>%
#     summarise(
#       Catch = median(Catch),
#       SSB = median(SSB),
#       Fbar = median(Fbar),
#       .groups = "drop"
#     )
#   
#   if (!is.null(new_model_names)) {
#     if (length(new_model_names) != length(unique(df_short_all$Model))) {
#       stop("Length of new_model_names must match number of models.")
#     }
#     df_short_all$Model <- factor(df_short_all$Model,
#                                  levels = paste0("Model", seq_along(new_model_names)),
#                                  labels = new_model_names)
#     df_long_all$Model <- factor(df_long_all$Model,
#                                 levels = paste0("Model", seq_along(new_model_names)),
#                                 labels = new_model_names)
#   }
#   
#   # === Plot function ===
#   plot_and_save <- function(df, title, file) {
#     colors <- viridisLite::viridis(n = length(unique(df$Model)), option = col.opt)
#     p <- ggtern(df, aes(x = Catch, y = SSB, z = Fbar, color = Model)) +
#       geom_point(alpha = 0.8, size = 8) +
#       scale_color_manual(values = colors) +
#       labs(title = title, T = "Fbar\nstatus", L = "SSB\nstatus", R = "Catch\nstatus") +
#       theme_rgbw() +
#       theme(plot.title = element_text(hjust = 0.5))
#     
#     ggsave(filename = file.path(main.dir, sub.dir, file), plot = p,
#            width = width, height = height, dpi = dpi)
#     
#     print(p)
#   }
#   
#   plot_and_save2 <- function(df, title, file) {
#     colors <- viridisLite::viridis(n = length(unique(df$Model)), option = col.opt)
#     p <- ggtern(df, aes(x = Catch, y = SSB, z = Fbar, color = Model)) +
#       geom_point(size = 1) +
#       scale_color_manual(values = colors) +
#       labs(title = title, T = "Fbar\nstatus", L = "SSB\nstatus", R = "Catch\nstatus") +
#       theme_rgbw() +
#       theme(plot.title = element_text(hjust = 0.5)) +
#       geom_confidence_tern(breaks = 0.95)
#     
#     ggsave(filename = file.path(main.dir, sub.dir, file), plot = p,
#            width = width, height = height, dpi = dpi)
#     
#     print(p)
#   }
#   
#   plot_and_save(df_short_all, paste0("Short-term Relative Status (Normalized): Years ", start.years, " to ", start.years + use.n.years.first - 1), "Relative_status_triangle_short.png")
#   plot_and_save(df_long_all, paste0("Long-term Relative Status (Normalized): Last ", use.n.years.last, " Years"), "Relative_status_triangle_long.png")
#   
#   results_short2 <- bind_rows(results_short)
#   results_long2 <- bind_rows(results_long)
#   
#   if (!is.null(new_model_names)) {
#     if (length(new_model_names) != length(unique(results_short2$Model))) {
#       stop("Length of new_model_names must match number of models.")
#     }
#     results_short2$Model <- factor(results_short2$Model,
#                                  levels = paste0("Model", seq_along(new_model_names)),
#                                  labels = new_model_names)
#     results_long2$Model <- factor(results_long2$Model,
#                                 levels = paste0("Model", seq_along(new_model_names)),
#                                 labels = new_model_names)
#   }
#   
#   plot_and_save2(results_short2, paste0("Short-term Relative Status: Years ", start.years, " to ", start.years + use.n.years.first - 1), "Relative_status_triangle_short_raw.png")
#   plot_and_save2(results_long2, paste0("Long-term Relative Status: Last ", use.n.years.last, " Years"), "Relative_status_triangle_long_raw.png")
# }


plot_status_triangle <- function(mods, is.nsim,
                                 main.dir, sub.dir,
                                 width = 8, height = 7, dpi = 300,
                                 col.opt = "D",
                                 method = NULL,
                                 new_model_names = NULL,
                                 use.n.years.first = 5,
                                 use.n.years.last = 5,
                                 start.years = 1,
                                 density_ridges = FALSE,
                                 show.whisker = FALSE) {
  library(dplyr)
  library(ggtern)
  library(viridisLite)
  
  # Safer warning
  if (missing(start.years)) {
    warning("start.years not provided; defaulting to 1 (first historical year).")
  }
  
  if (!is.nsim) {
    if (is.null(mods[[1]]$om$rep$log_SSB_FXSPR)) {
      message("Biological Reference Point has not been calculated internally!")
      return(invisible(NULL))
    }
  } else {
    if (is.null(mods[[1]][[1]]$om$rep$log_SSB_FXSPR)) {
      message("Biological Reference Point has not been calculated internally!")
      return(invisible(NULL))
    }
  }
  
  process_scores <- function(rep, n_fleets, n_regions, use.n.years, start.years = NULL, type = c("short", "long")) {
    # --- Calculate time series ratios ---
    catch_ts <- rowSums(rep$pred_catch) / exp(rep$log_Y_FXSPR[, ncol(rep$Fbar)])
    ssb_ts <- rowSums(rep$SSB) / exp(rep$log_SSB_FXSPR[, n_regions + 1])
    fbar_ts <- rep$Fbar[, ncol(rep$Fbar)] / exp(rep$log_Fbar_XSPR[, ncol(rep$Fbar)])
    
    # --- Define time window index ---
    if (type == "short") {
      idx <- start.years:(start.years + use.n.years - 1)
    } else {
      idx <- (length(catch_ts) - use.n.years + 1):length(catch_ts)
    }
    
    # --- Summarize time series, ignoring NA values ---
    if (is.null(method)) method <- "median"
    
    if (method == "median") {
      return(c(
        median(catch_ts[idx], na.rm = TRUE),
        median(ssb_ts[idx], na.rm = TRUE),
        median(fbar_ts[idx], na.rm = TRUE)
      ))
    }
    
    if (method == "mean") {
      return(c(
        mean(catch_ts[idx], na.rm = TRUE),
        mean(ssb_ts[idx], na.rm = TRUE),
        mean(fbar_ts[idx], na.rm = TRUE)
      ))
    }
  }
  
  res_short_raw <- list()
  res_long_raw <- list()
  res_short <- list()
  res_long <- list()
  results_short <- list()
  results_long <- list()
  
  if (is.nsim) {
    n_models <- length(mods[[1]])
    n_reps <- length(mods)
    
    for (r in seq_len(n_reps)) {
      df_short <- df_long <- data.frame(Model = paste0("Model", seq_len(n_models)))
      for (m in seq_len(n_models)) {
        rep <- mods[[r]][[m]]$om$rep
        input <- mods[[r]][[m]]$om$input$data
        n_fleets <- input$n_fleets[1]
        n_regions <- input$n_regions[1]
        
        short_vals <- process_scores(rep, n_fleets, n_regions, use.n.years.first, start.years, "short")
        long_vals <- process_scores(rep, n_fleets, n_regions, use.n.years.last, NULL, "long")
        
        df_short[m, c("Catch", "SSB", "Fbar")] <- short_vals
        df_long[m, c("Catch", "SSB", "Fbar")] <- long_vals
      }
      
      res_short_raw[[r]] <- df_short
      res_long_raw[[r]] <- df_long
      
      # Normalize within realization (Catch & SSB up, Fbar inverted)
      for (v in c("Catch", "SSB")) {
        range_v <- max(df_short[[v]]) - min(df_short[[v]])
        df_short[[v]] <- if (range_v == 0) 1 else (df_short[[v]] - min(df_short[[v]])) / range_v
        
        range_v_long <- max(df_long[[v]]) - min(df_long[[v]])
        df_long[[v]] <- if (range_v_long == 0) 1 else (df_long[[v]] - min(df_long[[v]])) / range_v_long
      }
      for (v in c("Fbar")) {
        range_v <- max(df_short[[v]]) - min(df_short[[v]])
        norm_f <- if (range_v == 0) 0 else (df_short[[v]] - min(df_short[[v]])) / range_v
        df_short[[v]] <- 1 - norm_f # inverted
        
        range_v_long <- max(df_long[[v]]) - min(df_long[[v]])
        norm_f_long <- if (range_v_long == 0) 0 else (df_long[[v]] - min(df_long[[v]])) / range_v_long
        df_long[[v]] <- 1 - norm_f_long
      }
      
      res_short[[r]] <- df_short
      res_long[[r]] <- df_long
      
      # Now we already have values 0–1: just renormalize to sum = 1
      df_short <- df_short %>%
        mutate(total = Catch + SSB + Fbar,
               Catch = ifelse(total == 0, 0, Catch / total),
               SSB = ifelse(total == 0, 0, SSB / total),
               Fbar = ifelse(total == 0, 0, Fbar / total),
               Realization = r)
      
      df_long <- df_long %>%
        mutate(total = Catch + SSB + Fbar,
               Catch = ifelse(total == 0, 0, Catch / total),
               SSB = ifelse(total == 0, 0, SSB / total),
               Fbar = ifelse(total == 0, 0, Fbar / total),
               Realization = r)
      
      results_short[[r]] <- df_short
      results_long[[r]] <- df_long
      
    }
  } else {
    n_models <- length(mods)
    df_short <- df_long <- data.frame(Model = paste0("Model", seq_len(n_models)))
    
    for (m in seq_len(n_models)) {
      rep <- mods[[m]]$om$rep
      input <- mods[[m]]$om$input$data
      n_fleets <- input$n_fleets[1]
      n_regions <- input$n_regions[1]
      
      short_vals <- process_scores(rep, n_fleets, n_regions, use.n.years.first, start.years, "short")
      long_vals <- process_scores(rep, n_fleets, n_regions, use.n.years.last, NULL, "long")
      
      df_short[m, c("Catch", "SSB", "Fbar")] <- short_vals
      df_long[m, c("Catch", "SSB", "Fbar")] <- long_vals
    }
    
    res_short_raw[[1]] <- df_short
    res_long_raw[[1]] <- df_long
    
    for (v in c("Catch", "SSB")) {
      range_v <- max(df_short[[v]]) - min(df_short[[v]])
      df_short[[v]] <- if (range_v == 0) 1 else (df_short[[v]] - min(df_short[[v]])) / range_v
      
      range_v_long <- max(df_long[[v]]) - min(df_long[[v]])
      df_long[[v]] <- if (range_v_long == 0) 1 else (df_long[[v]] - min(df_long[[v]])) / range_v_long
    }
    for (v in c("Fbar")) {
      range_v <- max(df_short[[v]]) - min(df_short[[v]])
      norm_f <- if (range_v == 0) 0 else (df_short[[v]] - min(df_short[[v]])) / range_v
      df_short[[v]] <- 1 - norm_f
      
      range_v_long <- max(df_long[[v]]) - min(df_long[[v]])
      norm_f_long <- if (range_v_long == 0) 0 else (df_long[[v]] - min(df_long[[v]])) / range_v_long
      df_long[[v]] <- 1 - norm_f_long
    }
    
    res_short[[1]] <- df_short
    res_long[[1]] <- df_long
    
    df_short <- df_short %>%
      mutate(total = Catch + SSB + Fbar,
             Catch = ifelse(total == 0, 0, Catch / total),
             SSB = ifelse(total == 0, 0, SSB / total),
             Fbar = ifelse(total == 0, 0, Fbar / total),
             Realization = 1)
    
    df_long <- df_long %>%
      mutate(total = Catch + SSB + Fbar,
             Catch = ifelse(total == 0, 0, Catch / total),
             SSB = ifelse(total == 0, 0, SSB / total),
             Fbar = ifelse(total == 0, 0, Fbar / total),
             Realization = 1)
    
    results_short[[1]] <- df_short
    results_long[[1]] <- df_long
  }
  
  df_short_all <- bind_rows(results_short)
  df_long_all <- bind_rows(results_long)
  
  df_short_all <- bind_rows(results_short) %>%
    group_by(Model) %>%
    summarise(
      Catch = median(Catch),
      SSB = median(SSB),
      Fbar = median(Fbar),
      .groups = "drop"
    )
  
  df_long_all <- bind_rows(results_long) %>%
    group_by(Model) %>%
    summarise(
      Catch = median(Catch),
      SSB = median(SSB),
      Fbar = median(Fbar),
      .groups = "drop"
    )
  
  if (!is.null(new_model_names)) {
    if (length(new_model_names) != length(unique(df_short_all$Model))) {
      stop("Length of new_model_names must match number of models.")
    }
    df_short_all$Model <- factor(df_short_all$Model,
                                 levels = paste0("Model", seq_along(new_model_names)),
                                 labels = new_model_names)
    df_long_all$Model <- factor(df_long_all$Model,
                                levels = paste0("Model", seq_along(new_model_names)),
                                labels = new_model_names)
  }
  
  # === Plot function ===
  plot_and_save <- function(df, title, file) {
    colors <- viridisLite::viridis(n = length(unique(df$Model)), option = col.opt)
    p <- ggtern(df, aes(x = Catch, y = SSB, z = Fbar, color = Model)) +
      geom_point(alpha = 0.8, size = 8) +
      scale_color_manual(values = colors) +
      labs(title = title, T = "Fbar\nstatus", L = "SSB\nstatus", R = "Catch\nstatus") +
      theme_rgbw() +
      theme(plot.title = element_text(hjust = 0.5))
    
    # Create the new subfolder if it doesn't exist
    new_sub_dir <- file.path(main.dir, sub.dir, "Status_Triangle_Plot")
    
    if (!file.exists(new_sub_dir)){
      dir.create(new_sub_dir)
    }
    
    ggsave(filename = file.path(main.dir, sub.dir, "Status_Triangle_Plot", file), plot = p,
           width = width, height = height, dpi = dpi)
    
    print(p)
  }
  
  plot_and_save2 <- function(df, title, file) {
    colors <- viridisLite::viridis(n = length(unique(df$Model)), option = col.opt)
    p <- ggtern(df, aes(x = Catch, y = SSB, z = Fbar, color = Model)) +
      geom_point(size = 1) +
      scale_color_manual(values = colors) +
      labs(title = title, T = "Fbar\nstatus", L = "SSB\nstatus", R = "Catch\nstatus") +
      theme_rgbw() +
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_confidence_tern(breaks = 0.95)
    
    # Create the new subfolder if it doesn't exist
    new_sub_dir <- file.path(main.dir, sub.dir, "Status_Triangle_Plot")
    
    if (!file.exists(new_sub_dir)){
      dir.create(new_sub_dir)
    }
    
    ggsave(filename = file.path(main.dir, sub.dir, "Status_Triangle_Plot", file), plot = p,
           width = width, height = height, dpi = dpi)
    
    print(p)
  }
  
  plot_and_save(df_short_all, paste0("Short-term Relative Status (Normalized): Years ", start.years, " to ", start.years + use.n.years.first - 1), "Status_Short.png")
  plot_and_save(df_long_all, paste0("Long-term Relative Status (Normalized): Last ", use.n.years.last, " Years"), "Status_Long.png")
  
  results_short2 <- bind_rows(results_short)
  results_long2 <- bind_rows(results_long)
  
  if (!is.null(new_model_names)) {
    if (length(new_model_names) != length(unique(results_short2$Model))) {
      stop("Length of new_model_names must match number of models.")
    }
    results_short2$Model <- factor(results_short2$Model,
                                   levels = paste0("Model", seq_along(new_model_names)),
                                   labels = new_model_names)
    results_long2$Model <- factor(results_long2$Model,
                                  levels = paste0("Model", seq_along(new_model_names)),
                                  labels = new_model_names)
  }
  
  plot_and_save2(results_short2, paste0("Short-term Relative Status: Years ", start.years, " to ", start.years + use.n.years.first - 1), "Status_Short(Raw).png")
  plot_and_save2(results_long2, paste0("Long-term Relative Status: Last ", use.n.years.last, " Years"), "Status_Long(Raw).png")
  
  res_short_all_raw <- bind_rows(res_short_raw)
  res_long_all_raw <- bind_rows(res_long_raw)
  res_short_all <- bind_rows(res_short)
  res_long_all <- bind_rows(res_long)
  
  if (!is.null(new_model_names)) {
    if (length(new_model_names) != length(unique(res_short_all$Model))) {
      stop("Length of new_model_names must match number of models.")
    }
    res_short_all$Model <- factor(res_short_all$Model,
                                  levels = paste0("Model", seq_along(new_model_names)),
                                  labels = new_model_names)
    res_long_all$Model <- factor(res_long_all$Model,
                                 levels = paste0("Model", seq_along(new_model_names)),
                                 labels = new_model_names)
    res_short_all_raw$Model <- factor(res_short_all_raw$Model,
                                      levels = paste0("Model", seq_along(new_model_names)),
                                      labels = new_model_names)
    res_long_all_raw$Model <- factor(res_long_all_raw$Model,
                                     levels = paste0("Model", seq_along(new_model_names)),
                                     labels = new_model_names)
  }
  
  # Raw values
  res_short_all2 <- res_short_all %>%
    pivot_longer(cols = c(Catch, SSB, Fbar),
                 names_to = "Variable", values_to = "Score") %>%
    mutate(Variable = gsub("_score", "", Variable))  # clean name
  
  if (density_ridges) {
    
    # Create the new subfolder if it doesn't exist
    new_sub_dir <- file.path(main.dir, sub.dir, "Status_Density_Plot")
    
    if (!file.exists(new_sub_dir)){
      dir.create(new_sub_dir)
    }
    
    p3 <- ggplot(res_short_all2, aes(x = Score, y = Model, fill = Model)) +
      ggridges::geom_density_ridges(alpha = 0.8, scale = 1.0, 
                          quantile_lines = TRUE, 
                          quantiles = 0.5,
                          linetype = 2) +
      theme_bw() +
      scale_fill_viridis_d(option = col.opt) +
      labs(
        title = paste0("Short-term Relative Status (Normalized): Years ", start.years, " to ", start.years + use.n.years.first - 1), 
        x = "Value",
        y = "Model"
      ) +
      facet_wrap(~ Variable)
    
  } else {
    
    # Create the new subfolder if it doesn't exist
    new_sub_dir <- file.path(main.dir, sub.dir, "Status_Bar_Plot")
    
    if (!file.exists(new_sub_dir)){
      dir.create(new_sub_dir)
    }
    
    p3 <- ggplot(res_short_all2, aes(x = Score, y = Model, fill = Model)) +
      geom_boxplot(width = 0.6, col = "red") +
      stat_summary(fun = median, geom = "point", shape = 21, fill = "white", size = 2) +
      theme_bw() +
      scale_fill_viridis_d(option = col.opt) +
      labs(
        title = paste0("Short-term Performance (Normalized): Years ", start.years, " to ", start.years + use.n.years.first - 1), 
        x = "Value",
        y = "Model"
      ) +
      facet_wrap(~ Variable)
  }

  
  print(p3)
  # Save & print
  ggsave(filename = file.path(new_sub_dir, "Status_Short2(Norm).png"),
         plot = p3, width = width, height = height, dpi = dpi)
  
  # Normalized
  library(tidyr)
  library(ggridges)
  library(ggplot2)
  
  # Raw values
  res_long_all2 <- res_long_all %>%
    pivot_longer(cols = c(Catch, SSB, Fbar),
                 names_to = "Variable", values_to = "Score") %>%
    mutate(Variable = gsub("_score", "", Variable))  # clean name
  
  if (density_ridges) {
    
    # Create the new subfolder if it doesn't exist
    new_sub_dir <- file.path(main.dir, sub.dir, "Status_Density_Plot")
    
    if (!file.exists(new_sub_dir)){
      dir.create(new_sub_dir)
    }
    
    p4 <- ggplot(res_long_all2, aes(x = Score, y = Model, fill = Model)) +
      ggridges::geom_density_ridges(alpha = 0.8, scale = 1.0, 
                          quantile_lines = TRUE, 
                          quantiles = 0.5,
                          linetype = 2) +
      theme_bw() +
      scale_fill_viridis_d(option = col.opt) +
      labs(
        title = paste0("Long-term Relative Status (Normalized): Last ", use.n.years.last, " Years"), 
        x = "Value",
        y = "Model"
      ) +
      facet_wrap(~ Variable)
    
  } else {
    
    # Create the new subfolder if it doesn't exist
    new_sub_dir <- file.path(main.dir, sub.dir, "Status_Bar_Plot")
    
    if (!file.exists(new_sub_dir)){
      dir.create(new_sub_dir)
    }
    
    p4 <- ggplot(res_long_all2, aes(x = Score, y = Model, fill = Model)) +
      geom_boxplot(width = 0.6, col = "red") +
      stat_summary(fun = median, geom = "point", shape = 21, fill = "white", size = 2) +
      theme_bw() +
      scale_fill_viridis_d(option = col.opt) +
      labs(
        title = paste0("Long-term Performance (Normalized): Last ", use.n.years.last, " Years"), 
        x = "Value",
        y = "Model"
      ) +
      facet_wrap(~ Variable)
  }
  
  print(p4)
  
  # Save & print
  ggsave(filename = file.path(new_sub_dir, "Status_Long2(Norm).png"),
         plot = p4, width = width, height = height, dpi = dpi)
  
  # Raw values
  res_short_all_raw2 <- res_short_all_raw %>%
    pivot_longer(cols = c(Catch, SSB, Fbar),
                 names_to = "Variable", values_to = "Score") %>%
    mutate(Variable = gsub("_score", "", Variable))  # clean name
  
  # Raw values
  res_long_all_raw2 <- res_long_all_raw %>%
    pivot_longer(cols = c(Catch, SSB, Fbar),
                 names_to = "Variable", values_to = "Score") %>%
    mutate(Variable = gsub("_score", "", Variable))  # clean name
  
  if (density_ridges) {
    
    # Create the new subfolder if it doesn't exist
    new_sub_dir <- file.path(main.dir, sub.dir, "Status_Density_Plot")
    
    if (!file.exists(new_sub_dir)){
      dir.create(new_sub_dir)
    }
    
    p5 <- ggplot(res_short_all_raw2, aes(x = Score, y = Model, fill = Model)) +
      ggridges::geom_density_ridges(alpha = 0.8, scale = 1.0, 
                          quantile_lines = TRUE, 
                          quantiles = 0.5,
                          linetype = 2) +
      theme_bw() +
      scale_fill_viridis_d(option = col.opt) +
      labs(
        title = paste0("Short-term Relative Status (Raw): Years ", start.years, " to ", start.years + use.n.years.first - 1), 
        x = "Value",
        y = "Model"
      ) +
      facet_wrap(~ Variable)
    
  } else {
    
    # Create the new subfolder if it doesn't exist
    new_sub_dir <- file.path(main.dir, sub.dir, "Status_Bar_Plot")
    
    if (!file.exists(new_sub_dir)){
      dir.create(new_sub_dir)
    }
    
    # Step 1: Summarize data
    res_summary <- res_short_all_raw2 %>%
      group_by(Model, Variable) %>%
      summarise(
        q1 = quantile(Score, 0.25, na.rm = TRUE),
        med = median(Score, na.rm = TRUE),
        q3 = quantile(Score, 0.75, na.rm = TRUE),
        iqr = q3 - q1,
        .groups = "drop"
      ) %>%
      mutate(
        y = as.numeric(factor(Model)),
        ymin = if (show.whisker) q1 - 1.5 * iqr else NA_real_,
        ymax = if (show.whisker) q3 + 1.5 * iqr else NA_real_
      )
    
    p5 <- ggplot(res_summary, aes(y = y, x = med, fill = Model)) +
      # Box body from Q1 to Q3
      geom_rect(aes(xmin = q1, xmax = q3, ymin = y - 0.3, ymax = y + 0.3, color = Model),
                fill = NA, linewidth = 0.8) +
      # Median line
      geom_segment(aes(x = med, xend = med, y = y - 0.3, yend = y + 0.3, color = Model),
                   linewidth = 0.8) +
      # Optional whiskers
      {if (show.whisker) geom_segment(aes(x = ymin, xend = q1, y = y, yend = y), linetype = "dashed")} +
      {if (show.whisker) geom_segment(aes(x = q3, xend = ymax, y = y, yend = y), linetype = "dashed")} +
      scale_y_continuous(breaks = res_summary$y, labels = res_summary$Model) +
      scale_color_viridis_d(option = col.opt) +
      theme_bw() +
      labs(
        title = paste0("Short-term Performance (Raw): Years ", start.years, " to ", start.years + use.n.years.first - 1), 
        x = "Value",
        y = "Model"
      ) +
      facet_wrap(~ Variable, scales = "free_x")
  }
  
  print(p5)
  # Save & print
  ggsave(filename = file.path(new_sub_dir, "Status_Short2(Raw).png"),
         plot = p5, width = width, height = height, dpi = dpi)
  
  if (density_ridges) {
    
    # Create the new subfolder if it doesn't exist
    new_sub_dir <- file.path(main.dir, sub.dir, "Status_Density_Plot")
    
    if (!file.exists(new_sub_dir)){
      dir.create(new_sub_dir)
    }
    
    p6 <- ggplot(res_long_all_raw2, aes(x = Score, y = Model, fill = Model)) +
      ggridges::geom_density_ridges(alpha = 0.8, scale = 1.0, 
                          quantile_lines = TRUE, 
                          quantiles = 0.5,
                          linetype = 2) +
      theme_bw() +
      scale_fill_viridis_d(option = col.opt) +
      labs(
        title = paste0("Long-term Relative Status (Raw): Last ", use.n.years.last, " Years"), 
        x = "Value",
        y = "Model"
      ) +
      facet_wrap(~ Variable)
    
  } else {
    
    # Create the new subfolder if it doesn't exist
    new_sub_dir <- file.path(main.dir, sub.dir, "Status_Bar_Plot")
    
    if (!file.exists(new_sub_dir)){
      dir.create(new_sub_dir)
    }
    
    # Step 1: Summarize data
    res_summary <- res_long_all_raw2 %>%
      group_by(Model, Variable) %>%
      summarise(
        q1 = quantile(Score, 0.25, na.rm = TRUE),
        med = median(Score, na.rm = TRUE),
        q3 = quantile(Score, 0.75, na.rm = TRUE),
        iqr = q3 - q1,
        .groups = "drop"
      ) %>%
      mutate(
        y = as.numeric(factor(Model)),
        ymin = if (show.whisker) q1 - 1.5 * iqr else NA_real_,
        ymax = if (show.whisker) q3 + 1.5 * iqr else NA_real_
      )
    
    # Step 2: Plot
    p6 <- ggplot(res_summary, aes(y = y, x = med, fill = Model)) +
      # Box body from Q1 to Q3
      geom_rect(aes(xmin = q1, xmax = q3, ymin = y - 0.3, ymax = y + 0.3, color = Model),
                fill = NA, linewidth = 0.8) +
      # Median line
      geom_segment(aes(x = med, xend = med, y = y - 0.3, yend = y + 0.3, color = Model),
                   linewidth = 0.8) +
      # Optional whiskers
      {if (show.whisker) geom_segment(aes(x = ymin, xend = q1, y = y, yend = y), linetype = "dashed")} +
      {if (show.whisker) geom_segment(aes(x = q3, xend = ymax, y = y, yend = y), linetype = "dashed")} +
      scale_y_continuous(breaks = res_summary$y, labels = res_summary$Model) +
      scale_color_viridis_d(option = col.opt) +
      theme_bw() +
      labs(
        title = paste0("Long-term Performance (Raw): Last ", use.n.years.last, " Years"),
        x = "Score",
        y = "Model"
      ) +
      facet_wrap(~ Variable, scales = "free_x")
  }
  
  print(p6)
  
  # Save & print
  ggsave(filename = file.path(new_sub_dir, "Status_Long2(Raw).png"),
         plot = p6, width = width, height = height, dpi = dpi)
  
}

# plot_model_performance_triangle2 <- function(mods, is.nsim,
#                                              main.dir, sub.dir,
#                                              width = 8, height = 7, dpi = 300,
#                                              col.opt = "D",
#                                              method = NULL,
#                                              new_model_names = NULL,
#                                              use.n.years.first = 5,
#                                              use.n.years.last = 5,
#                                              start.years = 1) {
#   library(dplyr)
#   library(ggtern)
#   library(viridisLite)
#   library(tidyr)
#   library(ggplot2)
#   library(rlang)
#   
#   # Helper to calculate AACV (Average Annual Catch Variation)
#   calculate_aacv <- function(catch_values) {
#     if (!is.numeric(catch_values)) {
#       stop("Input catch_values must be a numeric vector.")
#     }
#     catch_diff <- abs(diff(catch_values))
#     aacv <- sum(catch_diff) / sum(catch_values[-length(catch_values)])
#     return(aacv)
#   }
#   
#   # Safer warning
#   if (missing(start.years)) {
#     warning("start.years not provided; defaulting to 1 (first historical year).")
#   }
#   
#   # Function to extract mean Catch, SSB, and Fbar
#   process_scores <- function(rep, n_fleets, n_regions, use.n.years, start.years = NULL, type = c("short", "long")) {
#     catch_ts <- rowSums(rep$pred_catch)
#     ssb_ts <- rowSums(rep$SSB)
#     fbar_ts <- rep$Fbar[, ncol(rep$Fbar)] # SAFER: use last column (global Fbar)
#     
#     if (type == "short") {
#       idx <- start.years:(start.years + use.n.years - 1)
#     } else {
#       idx <- (length(catch_ts) - use.n.years + 1):length(catch_ts)
#     }
#     
#     if(is.null(method)) method = "median"
#     if(method == "median") {
#       return(c(calculate_aacv(catch_ts[idx]), median(ssb_ts[idx]), median(catch_ts[idx])))
#     }
#     if(method == "mean") {
#       return(c(calculate_aacv(catch_ts[idx]), mean(ssb_ts[idx]), mean(catch_ts[idx])))
#     }
#   }
#   
#   results_short <- list()
#   results_long <- list()
#   
#   if (is.nsim) {
#     n_models <- length(mods[[1]])
#     n_reps <- length(mods)
#     
#     for (r in seq_len(n_reps)) {
#       df_short <- df_long <- data.frame(Model = paste0("Model", seq_len(n_models)))
#       for (m in seq_len(n_models)) {
#         rep <- mods[[r]][[m]]$om$rep
#         input <- mods[[r]][[m]]$om$input$data
#         n_fleets <- input$n_fleets[1]
#         n_regions <- input$n_regions[1]
#         
#         short_vals <- process_scores(rep, n_fleets, n_regions, use.n.years.first, start.years, "short")
#         long_vals <- process_scores(rep, n_fleets, n_regions, use.n.years.last, NULL, "long")
#         
#         df_short[m, c("Catch", "SSB", "Fbar")] <- short_vals
#         df_long[m, c("Catch", "SSB", "Fbar")] <- long_vals
#       }
#       
#       # Normalize within realization (Catch & SSB up, Fbar inverted)
#       for (v in c("Catch", "SSB")) {
#         range_v <- max(df_short[[v]]) - min(df_short[[v]])
#         df_short[[v]] <- if (range_v == 0) 1 else (df_short[[v]] - min(df_short[[v]])) / range_v
#         
#         range_v_long <- max(df_long[[v]]) - min(df_long[[v]])
#         df_long[[v]] <- if (range_v_long == 0) 1 else (df_long[[v]] - min(df_long[[v]])) / range_v_long
#       }
#       for (v in c("Fbar")) {
#         range_v <- max(df_short[[v]]) - min(df_short[[v]])
#         norm_f <- if (range_v == 0) 0 else (df_short[[v]] - min(df_short[[v]])) / range_v
#         df_short[[v]] <- 1 - norm_f # inverted
#         
#         range_v_long <- max(df_long[[v]]) - min(df_long[[v]])
#         norm_f_long <- if (range_v_long == 0) 0 else (df_long[[v]] - min(df_long[[v]])) / range_v_long
#         df_long[[v]] <- 1 - norm_f_long
#       }
#       
#       # Now we already have values 0–1: just renormalize to sum = 1
#       df_short <- df_short %>%
#         mutate(total = Catch + SSB + Fbar,
#                Catch = ifelse(total == 0, 0, Catch / total),
#                SSB = ifelse(total == 0, 0, SSB / total),
#                Fbar = ifelse(total == 0, 0, Fbar / total),
#                Realization = r)
#       
#       df_long <- df_long %>%
#         mutate(total = Catch + SSB + Fbar,
#                Catch = ifelse(total == 0, 0, Catch / total),
#                SSB = ifelse(total == 0, 0, SSB / total),
#                Fbar = ifelse(total == 0, 0, Fbar / total),
#                Realization = r)
#       
#       results_short[[r]] <- df_short
#       results_long[[r]] <- df_long
#     }
#   } else {
#     n_models <- length(mods)
#     df_short <- df_long <- data.frame(Model = paste0("Model", seq_len(n_models)))
#     
#     for (m in seq_len(n_models)) {
#       rep <- mods[[m]]$om$rep
#       input <- mods[[m]]$om$input$data
#       n_fleets <- input$n_fleets[1]
#       n_regions <- input$n_regions[1]
#       
#       short_vals <- process_scores(rep, n_fleets, n_regions, use.n.years.first, start.years, "short")
#       long_vals <- process_scores(rep, n_fleets, n_regions, use.n.years.last, NULL, "long")
#       
#       df_short[m, c("Catch", "SSB", "Fbar")] <- short_vals
#       df_long[m, c("Catch", "SSB", "Fbar")] <- long_vals
#     }
#     
#     for (v in c("Catch", "SSB")) {
#       range_v <- max(df_short[[v]]) - min(df_short[[v]])
#       df_short[[v]] <- if (range_v == 0) 1 else (df_short[[v]] - min(df_short[[v]])) / range_v
#       
#       range_v_long <- max(df_long[[v]]) - min(df_long[[v]])
#       df_long[[v]] <- if (range_v_long == 0) 1 else (df_long[[v]] - min(df_long[[v]])) / range_v_long
#     }
#     for (v in c("Fbar")) {
#       range_v <- max(df_short[[v]]) - min(df_short[[v]])
#       norm_f <- if (range_v == 0) 0 else (df_short[[v]] - min(df_short[[v]])) / range_v
#       df_short[[v]] <- 1 - norm_f
#       
#       range_v_long <- max(df_long[[v]]) - min(df_long[[v]])
#       norm_f_long <- if (range_v_long == 0) 0 else (df_long[[v]] - min(df_long[[v]])) / range_v_long
#       df_long[[v]] <- 1 - norm_f_long
#     }
#     
#     df_short <- df_short %>%
#       mutate(total = Catch + SSB + Fbar,
#              Catch = ifelse(total == 0, 0, Catch / total),
#              SSB = ifelse(total == 0, 0, SSB / total),
#              Fbar = ifelse(total == 0, 0, Fbar / total),
#              Realization = 1)
#     
#     df_long <- df_long %>%
#       mutate(total = Catch + SSB + Fbar,
#              Catch = ifelse(total == 0, 0, Catch / total),
#              SSB = ifelse(total == 0, 0, SSB / total),
#              Fbar = ifelse(total == 0, 0, Fbar / total),
#              Realization = 1)
#     
#     results_short[[1]] <- df_short
#     results_long[[1]] <- df_long
#   }
#   
#   df_short_all <- bind_rows(results_short)
#   df_long_all <- bind_rows(results_long)
#   
#   df_short_all <- bind_rows(results_short) %>%
#     group_by(Model) %>%
#     summarise(
#       Catch = median(Catch),
#       SSB = median(SSB),
#       Fbar = median(Fbar),
#       .groups = "drop"
#     )
#   
#   df_long_all <- bind_rows(results_long) %>%
#     group_by(Model) %>%
#     summarise(
#       Catch = median(Catch),
#       SSB = median(SSB),
#       Fbar = median(Fbar),
#       .groups = "drop"
#     )
#   
#   if (!is.null(new_model_names)) {
#     if (length(new_model_names) != length(unique(df_short_all$Model))) {
#       stop("Length of new_model_names must match number of models.")
#     }
#     df_short_all$Model <- factor(df_short_all$Model,
#                                  levels = paste0("Model", seq_along(new_model_names)),
#                                  labels = new_model_names)
#     df_long_all$Model <- factor(df_long_all$Model,
#                                 levels = paste0("Model", seq_along(new_model_names)),
#                                 labels = new_model_names)
#   }
#   
#   # === Plot function ===
#   plot_and_save <- function(df, title, file) {
#     colors <- viridisLite::viridis(n = length(unique(df$Model)), option = col.opt)
#     p <- ggtern(df, aes(x = Catch, y = SSB, z = Fbar, color = Model)) +
#       geom_point(alpha = 0.8, size = 8) +
#       scale_color_manual(values = colors) +
#       labs(title = title, T = "AACV", L = "SSB", R = "Catch") +
#       theme_rgbw() +
#       theme(plot.title = element_text(hjust = 0.5))
#     
#     ggsave(filename = file.path(main.dir, sub.dir, file), plot = p,
#            width = width, height = height, dpi = dpi)
#     
#     print(p)
#   }
#   
#   plot_and_save2 <- function(df, title, file) {
#     colors <- viridisLite::viridis(n = length(unique(df$Model)), option = col.opt)
#     p <- ggtern(df, aes(x = Catch, y = SSB, z = Fbar, color = Model)) +
#       geom_point(size = 1) +
#       scale_color_manual(values = colors) +
#       labs(title = title, T = "AACV", L = "SSB", R = "Catch") +
#       theme_rgbw() +
#       theme(plot.title = element_text(hjust = 0.5)) +
#       geom_confidence_tern(breaks = 0.95)
#     
#     ggsave(filename = file.path(main.dir, sub.dir, file), plot = p,
#            width = width, height = height, dpi = dpi)
#     
#     print(p)
#   }
#   
#   plot_and_save(df_short_all, paste0("Short-term Performance (Normalized): Years ", start.years, " to ", start.years + use.n.years.first - 1), "model_performance_triangle_short2.png")
#   plot_and_save(df_long_all, paste0("Long-term Performance (Normalized): Last ", use.n.years.last, " Years"), "model_performance_triangle_long2.png")
#   
#   results_short2 <- bind_rows(results_short)
#   results_long2 <- bind_rows(results_long)
#   
#   if (!is.null(new_model_names)) {
#     if (length(new_model_names) != length(unique(results_short2$Model))) {
#       stop("Length of new_model_names must match number of models.")
#     }
#     results_short2$Model <- factor(results_short2$Model,
#                                    levels = paste0("Model", seq_along(new_model_names)),
#                                    labels = new_model_names)
#     results_long2$Model <- factor(results_long2$Model,
#                                   levels = paste0("Model", seq_along(new_model_names)),
#                                   labels = new_model_names)
#   }
#   
#   plot_and_save2(results_short2, paste0("Short-term Performance: Years ", start.years, " to ", start.years + use.n.years.first - 1), "model_performance_triangle_short_raw2.png")
#   plot_and_save2(results_long2, paste0("Long-term Performance: Last ", use.n.years.last, " Years"), "model_performance_triangle_long_raw2.png")
# }

plot_model_performance_triangle2 <- function(mods, is.nsim,
                                             main.dir, sub.dir,
                                             width = 8, height = 7, dpi = 300,
                                             col.opt = "D",
                                             method = NULL,
                                             new_model_names = NULL,
                                             use.n.years.first = 5,
                                             use.n.years.last = 5,
                                             start.years = 1,
                                             density_ridges = FALSE,
                                             show.whisker = FALSE) {
  library(dplyr)
  library(ggtern)
  library(viridisLite)
  library(tidyr)
  library(ggplot2)
  library(rlang)
  
  # Helper to calculate AACV (Average Annual Catch Variation)
  calculate_aacv <- function(catch_values) {
    if (!is.numeric(catch_values)) {
      stop("Input catch_values must be a numeric vector.")
    }
    catch_diff <- abs(diff(catch_values))
    aacv <- sum(catch_diff) / sum(catch_values[-length(catch_values)])
    return(aacv)
  }
  
  # Safer warning
  if (missing(start.years)) {
    warning("start.years not provided; defaulting to 1 (first historical year).")
  }
  
  # Function to extract mean Catch, SSB, and Fbar
  process_scores <- function(rep, n_fleets, n_regions, use.n.years, start.years = NULL, type = c("short", "long")) {
    catch_ts <- rowSums(rep$pred_catch)
    ssb_ts <- rowSums(rep$SSB)
    
    if (type == "short") {
      idx <- start.years:(start.years + use.n.years - 1)
    } else {
      idx <- (length(catch_ts) - use.n.years + 1):length(catch_ts)
    }
    
    if(is.null(method)) method = "median"
    if(method == "median") {
      return(c(calculate_aacv(catch_ts[idx]), median(ssb_ts[idx]), median(catch_ts[idx])))
    }
    if(method == "mean") {
      return(c(calculate_aacv(catch_ts[idx]), mean(ssb_ts[idx]), mean(catch_ts[idx])))
    }
  }
  
  results_short <- list()
  results_long <- list()
  res_short_raw <- list()
  res_long_raw <- list()
  res_short <- list()
  res_long <- list()
  results_short <- list()
  results_long <- list()
  
  if (is.nsim) {
    n_models <- length(mods[[1]])
    n_reps <- length(mods)
    
    for (r in seq_len(n_reps)) {
      df_short <- df_long <- data.frame(Model = paste0("Model", seq_len(n_models)))
      for (m in seq_len(n_models)) {
        rep <- mods[[r]][[m]]$om$rep
        input <- mods[[r]][[m]]$om$input$data
        n_fleets <- input$n_fleets[1]
        n_regions <- input$n_regions[1]
        
        short_vals <- process_scores(rep, n_fleets, n_regions, use.n.years.first, start.years, "short")
        long_vals <- process_scores(rep, n_fleets, n_regions, use.n.years.last, NULL, "long")
        
        df_short[m, c("AACV", "SSB", "Catch")] <- short_vals
        df_long[m, c("AACV", "SSB", "Catch")] <- long_vals
      }
      
      res_short_raw[[r]] <- df_short
      res_long_raw[[r]] <- df_long
      
      # Normalize
      for (v in c("SSB", "Catch")) {
        range_v <- max(df_short[[v]]) - min(df_short[[v]])
        df_short[[v]] <- if (range_v == 0) 1 else (df_short[[v]] - min(df_short[[v]])) / range_v
        
        range_v_long <- max(df_long[[v]]) - min(df_long[[v]])
        df_long[[v]] <- if (range_v_long == 0) 1 else (df_long[[v]] - min(df_long[[v]])) / range_v_long
      }
      for (v in c("AACV")) {
        range_v <- max(df_short[[v]]) - min(df_short[[v]])
        norm_f <- if (range_v == 0) 0 else (df_short[[v]] - min(df_short[[v]])) / range_v
        df_short[[v]] <- 1 - norm_f # inverted
        
        range_v_long <- max(df_long[[v]]) - min(df_long[[v]])
        norm_f_long <- if (range_v_long == 0) 0 else (df_long[[v]] - min(df_long[[v]])) / range_v_long
        df_long[[v]] <- 1 - norm_f_long
      }
      
      res_short[[r]] <- df_short
      res_long[[r]] <- df_long
      
      # Now we already have values 0–1: just renormalize to sum = 1
      df_short <- df_short %>%
        mutate(total = AACV + SSB + Catch,
               AACV = ifelse(total == 0, 0, AACV / total),
               SSB = ifelse(total == 0, 0, SSB / total),
               Catch = ifelse(total == 0, 0, Catch / total),
               Realization = r)
      
      df_long <- df_long %>%
        mutate(total = AACV + SSB + Catch,
               AACV = ifelse(total == 0, 0, AACV / total),
               SSB = ifelse(total == 0, 0, SSB / total),
               Catch = ifelse(total == 0, 0, Catch / total),
               Realization = r)
      
      results_short[[r]] <- df_short
      results_long[[r]] <- df_long
    }
  } else {
    n_models <- length(mods)
    df_short <- df_long <- data.frame(Model = paste0("Model", seq_len(n_models)))
    
    for (m in seq_len(n_models)) {
      rep <- mods[[m]]$om$rep
      input <- mods[[m]]$om$input$data
      n_fleets <- input$n_fleets[1]
      n_regions <- input$n_regions[1]
      
      short_vals <- process_scores(rep, n_fleets, n_regions, use.n.years.first, start.years, "short")
      long_vals <- process_scores(rep, n_fleets, n_regions, use.n.years.last, NULL, "long")
      
      df_short[m, c("AACV", "SSB", "Catch")] <- short_vals
      df_long[m, c("AACV", "SSB", "Catch")] <- long_vals
    }
    
    res_short_raw[[1]] <- df_short
    res_long_raw[[1]] <- df_long
    
    for (v in c("SSB", "Catch")) {
      range_v <- max(df_short[[v]]) - min(df_short[[v]])
      df_short[[v]] <- if (range_v == 0) 1 else (df_short[[v]] - min(df_short[[v]])) / range_v
      
      range_v_long <- max(df_long[[v]]) - min(df_long[[v]])
      df_long[[v]] <- if (range_v_long == 0) 1 else (df_long[[v]] - min(df_long[[v]])) / range_v_long
    }
    for (v in c("AACV")) {
      range_v <- max(df_short[[v]]) - min(df_short[[v]])
      norm_f <- if (range_v == 0) 0 else (df_short[[v]] - min(df_short[[v]])) / range_v
      df_short[[v]] <- 1 - norm_f
      
      range_v_long <- max(df_long[[v]]) - min(df_long[[v]])
      norm_f_long <- if (range_v_long == 0) 0 else (df_long[[v]] - min(df_long[[v]])) / range_v_long
      df_long[[v]] <- 1 - norm_f_long
    }
    
    res_short[[1]] <- df_short
    res_long[[1]] <- df_long
    
    df_short <- df_short %>%
      mutate(total = AACV + SSB + Catch,
             AACV = ifelse(total == 0, 0, AACV / total),
             SSB = ifelse(total == 0, 0, SSB / total),
             Catch = ifelse(total == 0, 0, Catch / total),
             Realization = 1)
    
    df_long <- df_long %>%
      mutate(total = AACV + SSB + Catch,
             AACV = ifelse(total == 0, 0, AACV / total),
             SSB = ifelse(total == 0, 0, SSB / total),
             Catch = ifelse(total == 0, 0, Catch / total),
             Realization = 1)
    
    results_short[[1]] <- df_short
    results_long[[1]] <- df_long
  }
  
  df_short_all <- bind_rows(results_short)
  df_long_all <- bind_rows(results_long)
  
  df_short_all <- bind_rows(results_short) %>%
    group_by(Model) %>%
    summarise(
      AACV = median(AACV),
      SSB = median(SSB),
      Catch = median(Catch),
      .groups = "drop"
    )
  
  df_long_all <- bind_rows(results_long) %>%
    group_by(Model) %>%
    summarise(
      AACV = median(AACV),
      SSB = median(SSB),
      Catch = median(Catch),
      .groups = "drop"
    )
  
  if (!is.null(new_model_names)) {
    if (length(new_model_names) != length(unique(df_short_all$Model))) {
      stop("Length of new_model_names must match number of models.")
    }
    df_short_all$Model <- factor(df_short_all$Model,
                                 levels = paste0("Model", seq_along(new_model_names)),
                                 labels = new_model_names)
    df_long_all$Model <- factor(df_long_all$Model,
                                levels = paste0("Model", seq_along(new_model_names)),
                                labels = new_model_names)
  }
  
  # === Plot function ===
  plot_and_save <- function(df, title, file) {
    colors <- viridisLite::viridis(n = length(unique(df$Model)), option = col.opt)
    p <- ggtern(df, aes(x = AACV, y = SSB, z = Catch, color = Model)) +
      geom_point(alpha = 0.8, size = 8) +
      scale_color_manual(values = colors) +
      labs(title = title, T = "AACV", L = "SSB", R = "Catch") +
      theme_rgbw() +
      theme(plot.title = element_text(hjust = 0.5))
    
    # Create the new subfolder if it doesn't exist
    new_sub_dir <- file.path(main.dir, sub.dir, "Holistic_Triangle_Plot")
    
    if (!file.exists(new_sub_dir)){
      dir.create(new_sub_dir)
    }
    
    # Save the figure inside the new subfolder
    ggsave(filename = file.path(new_sub_dir, file), plot = p, width = width, height = height, dpi = dpi)
    
    print(p)
  }
  
  plot_and_save2 <- function(df, title, file) {
    colors <- viridisLite::viridis(n = length(unique(df$Model)), option = col.opt)
    p <- ggtern(df, aes(x = AACV, y = SSB, z = Catch, color = Model)) +
      geom_point(size = 1) +
      scale_color_manual(values = colors) +
      labs(title = title, T = "AACV", L = "SSB", R = "Catch") +
      theme_rgbw() +
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_confidence_tern(breaks = 0.95)
    
    # Create the new subfolder if it doesn't exist
    new_sub_dir <- file.path(main.dir, sub.dir, "Holistic_Triangle_Plot")
    
    if (!file.exists(new_sub_dir)){
      dir.create(new_sub_dir)
    }
    
    # Save the figure inside the new subfolder
    ggsave(filename = file.path(new_sub_dir, file), plot = p, width = width, height = height, dpi = dpi)
    
    print(p)
  }
  
  plot_and_save(df_short_all, paste0("Short-term Performance (Normalized): Years ", start.years, " to ", start.years + use.n.years.first - 1), "Performance_Short2(Norm).png")
  plot_and_save(df_long_all, paste0("Long-term Performance (Normalized): Last ", use.n.years.last, " Years"), "Performance_Long2(Norm).png")
  
  results_short2 <- bind_rows(results_short)
  results_long2 <- bind_rows(results_long)
  
  if (!is.null(new_model_names)) {
    if (length(new_model_names) != length(unique(results_short2$Model))) {
      stop("Length of new_model_names must match number of models.")
    }
    results_short2$Model <- factor(results_short2$Model,
                                   levels = paste0("Model", seq_along(new_model_names)),
                                   labels = new_model_names)
    results_long2$Model <- factor(results_long2$Model,
                                  levels = paste0("Model", seq_along(new_model_names)),
                                  labels = new_model_names)
  }
  
  plot_and_save2(results_short2, paste0("Short-term Performance: Years ", start.years, " to ", start.years + use.n.years.first - 1), "Performance_Short2(Raw).png")
  plot_and_save2(results_long2, paste0("Long-term Performance: Last ", use.n.years.last, " Years"), "Performance_Long2(Raw).png")
  
  res_short_all_raw <- bind_rows(res_short_raw)
  res_long_all_raw <- bind_rows(res_long_raw)
  res_short_all <- bind_rows(res_short)
  res_long_all <- bind_rows(res_long)
  
  if (!is.null(new_model_names)) {
    if (length(new_model_names) != length(unique(res_short_all$Model))) {
      stop("Length of new_model_names must match number of models.")
    }
    res_short_all$Model <- factor(res_short_all$Model,
                                  levels = paste0("Model", seq_along(new_model_names)),
                                  labels = new_model_names)
    res_long_all$Model <- factor(res_long_all$Model,
                                 levels = paste0("Model", seq_along(new_model_names)),
                                 labels = new_model_names)
    res_short_all_raw$Model <- factor(res_short_all_raw$Model,
                                      levels = paste0("Model", seq_along(new_model_names)),
                                      labels = new_model_names)
    res_long_all_raw$Model <- factor(res_long_all_raw$Model,
                                     levels = paste0("Model", seq_along(new_model_names)),
                                     labels = new_model_names)
  }
  
  # Raw values
  res_short_all2 <- res_short_all %>%
    pivot_longer(cols = c(AACV, SSB, Catch),
                 names_to = "Variable", values_to = "Score") %>%
    mutate(Variable = gsub("_score", "", Variable))  # clean name
  
  if (density_ridges) {
    # Create the new subfolder if it doesn't exist
    new_sub_dir <- file.path(main.dir, sub.dir, "Holistic_Density_Plot")
    
    if (!file.exists(new_sub_dir)){
      dir.create(new_sub_dir)
    }
    
    p3 <- ggplot(res_short_all2, aes(x = Score, y = Model, fill = Model)) +
      ggridges::geom_density_ridges(alpha = 0.8, scale = 1.0, 
                          quantile_lines = TRUE, 
                          quantiles = 0.5,
                          linetype = 2) +
      theme_bw() +
      scale_fill_viridis_d(option = col.opt) +
      labs(
        title = paste0("Short-term Performance (Normalized): Years ", start.years, " to ", start.years + use.n.years.first - 1), 
        x = "Value",
        y = "Model"
      ) +
      facet_wrap(~ Variable)
    
  } else {
    
    # Create the new subfolder if it doesn't exist
    new_sub_dir <- file.path(main.dir, sub.dir, "Holistic_Bar_Plot")
    
    if (!file.exists(new_sub_dir)){
      dir.create(new_sub_dir)
    }
    
    p3 <- ggplot(res_short_all2, aes(x = Score, y = Model, fill = Model)) +
      geom_boxplot(width = 0.6, col = "red") +
      stat_summary(fun = median, geom = "point", shape = 21, fill = "white", size = 2) +
      theme_bw() +
      scale_fill_viridis_d(option = col.opt) +
      labs(
        title = paste0("Short-term Performance (Normalized): Years ", start.years, " to ", start.years + use.n.years.first - 1), 
        x = "Value",
        y = "Model"
      ) +
      facet_wrap(~ Variable)
  }
  
  print(p3)
  # Save & print
  ggsave(filename = file.path(new_sub_dir, "Performance_Short2(Norm).png"), plot = p3, width = width, height = height, dpi = dpi)
  
  # Normalized
  library(tidyr)
  library(ggridges)
  library(ggplot2)
  
  # Raw values
  res_long_all2 <- res_long_all %>%
    pivot_longer(cols = c(AACV, SSB, Catch),
                 names_to = "Variable", values_to = "Score") %>%
    mutate(Variable = gsub("_score", "", Variable))  # clean name
  
  if (density_ridges) {
    
    # Create the new subfolder if it doesn't exist
    new_sub_dir <- file.path(main.dir, sub.dir, "Holistic_Density_Plot")
    
    if (!file.exists(new_sub_dir)){
      dir.create(new_sub_dir)
    }
    
    p4 <- ggplot(res_long_all2, aes(x = Score, y = Model, fill = Model)) +
      ggridges::geom_density_ridges(alpha = 0.8, scale = 1.0, 
                          quantile_lines = TRUE, 
                          quantiles = 0.5,
                          linetype = 2) +
      theme_bw() +
      scale_fill_viridis_d(option = col.opt) +
      labs(
        title = paste0("Long-term Performance (Normalized): Last ", use.n.years.last, " Years"), 
        x = "Value",
        y = "Model"
      ) +
      facet_wrap(~ Variable)
    
  } else {
    
    # Create the new subfolder if it doesn't exist
    new_sub_dir <- file.path(main.dir, sub.dir, "Holistic_Bar_Plot")
    
    if (!file.exists(new_sub_dir)){
      dir.create(new_sub_dir)
    }
    
    p4 <- ggplot(res_long_all2, aes(x = Score, y = Model, fill = Model)) +
      geom_boxplot(width = 0.6, col = "red") +
      stat_summary(fun = median, geom = "point", shape = 21, fill = "white", size = 2) +
      theme_bw() +
      scale_fill_viridis_d(option = col.opt) +
      labs(
        title = paste0("Long-term Performance (Normalized): Last ", use.n.years.last, " Years"), 
        x = "Value",
        y = "Model"
      ) +
      facet_wrap(~ Variable)
  }
  
  print(p4)
  
  # Save & print
  ggsave(filename = file.path(new_sub_dir, "Performance_Long2(Norm).png"), plot = p4, width = width, height = height, dpi = dpi)
  
  # Raw values
  res_short_all_raw2 <- res_short_all_raw %>%
    pivot_longer(cols = c(AACV, SSB, Catch),
                 names_to = "Variable", values_to = "Score") %>%
    mutate(Variable = gsub("_score", "", Variable))  # clean name
  
  # Raw values
  res_long_all_raw2 <- res_long_all_raw %>%
    pivot_longer(cols = c(AACV, SSB, Catch),
                 names_to = "Variable", values_to = "Score") %>%
    mutate(Variable = gsub("_score", "", Variable))  # clean name
  
  if (density_ridges) {
    
    # Create the new subfolder if it doesn't exist
    new_sub_dir <- file.path(main.dir, sub.dir, "Holistic_Density_Plot")
    
    if (!file.exists(new_sub_dir)){
      dir.create(new_sub_dir)
    }
    
    p5 <- ggplot(res_short_all_raw2, aes(x = Score, y = Model, fill = Model)) +
      ggridges::geom_density_ridges(alpha = 0.8, scale = 1.0, 
                          quantile_lines = TRUE, 
                          quantiles = 0.5,
                          linetype = 2
                          ) +
      theme_bw() +
      scale_fill_viridis_d(option = col.opt) +
      labs(
        title = paste0("Short-term Performance (Raw): Years ", start.years, " to ", start.years + use.n.years.first - 1), 
        x = "Value",
        y = "Model"
      ) +
      facet_wrap(~ Variable)
    
  } else {
    
    # Create the new subfolder if it doesn't exist
    new_sub_dir <- file.path(main.dir, sub.dir, "Holistic_Bar_Plot")
    
    if (!file.exists(new_sub_dir)){
      dir.create(new_sub_dir)
    }
    
    # Step 1: Summarize data
    res_summary <- res_short_all_raw2 %>%
      group_by(Model, Variable) %>%
      summarise(
        q1 = quantile(Score, 0.25, na.rm = TRUE),
        med = median(Score, na.rm = TRUE),
        q3 = quantile(Score, 0.75, na.rm = TRUE),
        iqr = q3 - q1,
        .groups = "drop"
      ) %>%
      mutate(
        y = as.numeric(factor(Model)),
        ymin = if (show.whisker) q1 - 1.5 * iqr else NA_real_,
        ymax = if (show.whisker) q3 + 1.5 * iqr else NA_real_
      )
    
    p5 <- ggplot(res_summary, aes(y = y, x = med, fill = Model)) +
      # Box body from Q1 to Q3
      geom_rect(aes(xmin = q1, xmax = q3, ymin = y - 0.3, ymax = y + 0.3, color = Model),
                fill = NA, linewidth = 0.8) +
      # Median line
      geom_segment(aes(x = med, xend = med, y = y - 0.3, yend = y + 0.3, color = Model),
                   linewidth = 0.8) +
      # Optional whiskers
      {if (show.whisker) geom_segment(aes(x = ymin, xend = q1, y = y, yend = y), linetype = "dashed")} +
      {if (show.whisker) geom_segment(aes(x = q3, xend = ymax, y = y, yend = y), linetype = "dashed")} +
      scale_y_continuous(breaks = res_summary$y, labels = res_summary$Model) +
      scale_color_viridis_d(option = col.opt) +
      theme_bw() +
      labs(
        title = paste0("Short-term Performance (Raw): Years ", start.years, " to ", start.years + use.n.years.first - 1), 
        x = "Value",
        y = "Model"
      ) +
      facet_wrap(~ Variable, scales = "free_x")
  }
  
  print(p5)
  # Save & print
  ggsave(filename = file.path(new_sub_dir, "Performance_Short2(Raw).png"), plot = p5, width = width, height = height, dpi = dpi)
  
  
  if (density_ridges) {
    
    # Create the new subfolder if it doesn't exist
    new_sub_dir <- file.path(main.dir, sub.dir, "Holistic_Density_Plot")
    
    if (!file.exists(new_sub_dir)){
      dir.create(new_sub_dir)
    }
    
    p6 <- ggplot(res_long_all_raw2, aes(x = Score, y = Model, fill = Model)) +
      ggridges::geom_density_ridges(alpha = 0.8, scale = 1.0, 
                          quantile_lines = TRUE, 
                          quantiles = 0.5,
                          linetype = 2) +
      theme_bw() +
      scale_fill_viridis_d(option = col.opt) +
      labs(
        title = paste0("Long-term Performance (Raw): Last ", use.n.years.last, " Years"), 
        x = "Value",
        y = "Model"
      ) +
      facet_wrap(~ Variable)
    
  } else {
    
    # Create the new subfolder if it doesn't exist
    new_sub_dir <- file.path(main.dir, sub.dir, "Holistic_Bar_Plot")
    
    if (!file.exists(new_sub_dir)){
      dir.create(new_sub_dir)
    }
    
    # Step 1: Summarize data
    res_summary <- res_long_all_raw2 %>%
      group_by(Model, Variable) %>%
      summarise(
        q1 = quantile(Score, 0.25, na.rm = TRUE),
        med = median(Score, na.rm = TRUE),
        q3 = quantile(Score, 0.75, na.rm = TRUE),
        iqr = q3 - q1,
        .groups = "drop"
      ) %>%
      mutate(
        y = as.numeric(factor(Model)),
        ymin = if (show.whisker) q1 - 1.5 * iqr else NA_real_,
        ymax = if (show.whisker) q3 + 1.5 * iqr else NA_real_
      )
    
    # Step 2: Plot
    p6 <- ggplot(res_summary, aes(y = y, x = med, fill = Model)) +
      # Box body from Q1 to Q3
      geom_rect(aes(xmin = q1, xmax = q3, ymin = y - 0.3, ymax = y + 0.3, color = Model),
                fill = NA, linewidth = 0.8) +
      # Median line
      geom_segment(aes(x = med, xend = med, y = y - 0.3, yend = y + 0.3, color = Model),
                   linewidth = 0.8) +
      # Optional whiskers
      {if (show.whisker) geom_segment(aes(x = ymin, xend = q1, y = y, yend = y), linetype = "dashed")} +
      {if (show.whisker) geom_segment(aes(x = q3, xend = ymax, y = y, yend = y), linetype = "dashed")} +
      scale_y_continuous(breaks = res_summary$y, labels = res_summary$Model) +
      scale_color_viridis_d(option = col.opt) +
      theme_bw() +
      labs(
        title = paste0("Long-term Performance (Raw): Last ", use.n.years.last, " Years"),
        x = "Score",
        y = "Model"
      ) +
      facet_wrap(~ Variable, scales = "free_x")
    
  }
  
  print(p6)
  
  # Save & print
  ggsave(filename = file.path(new_sub_dir, "Performance_Long2(Raw).png"), plot = p6, width = width, height = height, dpi = dpi)
  
}


plot_total_catch_performance <- function(mods, is.nsim, main.dir, sub.dir, var = "Catch",
                                         width = 10, height = 7, dpi = 300, col.opt = "D",
                                         method = NULL,
                                         outlier.opt = NA,
                                         plot.style = "median_iqr", 
                                         show.whisker = TRUE,
                                         use.n.years = NULL,
                                         new_model_names = NULL,
                                         base.model = NULL) {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(rlang)
  
  if (is.null(use.n.years)) {
    cat("\nuse.n.years is not specified, so default (5 years) is used here!\n")
    use.n.years <- 5
  }
  
  res <- NULL
  
  # Prepare data
  if (!is.nsim) {
    Years <- mods[[1]]$om$years
    res <- lapply(seq_along(mods), function(i) {
      data.frame(
        Catch = mods[[i]]$om$rep$pred_catch,
        Model = paste0("Model", i),
        Year = Years,
        Realization = 1
      ) %>% tail(use.n.years)
    }) %>% bind_rows()
  } else {
    Years <- mods[[1]][[1]]$om$years
    res <- lapply(seq_along(mods), function(r) {
      lapply(seq_along(mods[[r]]), function(m) {
        data.frame(
          Catch = mods[[r]][[m]]$om$rep$pred_catch,
          Model = paste0("Model", m),
          Year = Years,
          Realization = r
        ) %>% tail(use.n.years)
      }) %>% bind_rows()
    }) %>% bind_rows()
  }
  
  # Add the last column: total Catch by Year and Realization
  res <- res %>%
    rowwise() %>%
    mutate(Catch_Global = sum(c_across(starts_with("Catch")), na.rm = TRUE)) %>%
    ungroup()
  
  # Rename models if specified
  if (!is.null(new_model_names)) {
    if (length(new_model_names) != length(unique(res$Model))) {
      stop("Length of new_model_names must match the number of models.")
    }
    res$Model <- factor(res$Model,
                        levels = paste0("Model", seq_along(new_model_names)),
                        labels = new_model_names)
    # if (!is.null(base.model)) base.model <- new_model_names[base.model]
    if (!is.null(base.model)) {
      if (!(base.model %in% new_model_names)) {
        warning("base.model does not match any of the new_model_names.")
      }
    }
  }
  
  # Pivot longer if needed (for multiple Catch columns)
  res <- pivot_longer(res, cols = starts_with(var), names_to = "Label", values_to = var)
  
  # Then sum over years (use.n.years) by Model, Realization, and Label
  res <- res %>%
    filter(Year %in% tail(Years, use.n.years)) %>%
    group_by(Model, Realization, Label) %>%
    summarise(!!var := sum(!!sym(var), na.rm = TRUE), .groups = "drop")
  
  # Relative difference from base.model
  if (!is.null(base.model)) {
    base_df <- res %>%
      filter(Model == base.model) %>%
      rename(base_val = !!sym(var)) %>%
      select(Realization, Label, base_val)
    
    res <- left_join(res, base_df, by = c("Realization", "Label")) %>%
      mutate(!!var := (!!sym(var)) / base_val - 1)
  }
  
  # Apply mean or median summarization if method is specified
  if (!is.null(method)) {
    res <- res %>%
      group_by(Model, Realization, Label) %>%
      summarise(!!var := if (method == "mean") mean(!!sym(var)) else median(!!sym(var)),
                .groups = "drop")
  }
  
  # Plot
  if (plot.style == "boxplot") {
    p <- ggplot(res, aes(x = Model, y = Catch, color = Model)) +
      geom_boxplot(lwd = 0.8, outlier.shape = outlier.opt) +
      facet_grid(Label ~ ., scales = "free") +
      scale_color_viridis_d(option = col.opt) +
      ggtitle(paste0(ifelse(is.null(base.model), paste0("Total ", var), paste0("Relative Total ", var, " vs ", base.model)),
                     ": Last ", use.n.years, " Years")) +
      ylab(ifelse(is.null(base.model), paste0("Total ", var), paste0("Relative Total ", var, " Difference"))) +
      xlab("Model") + 
      theme_bw()
  } else if (plot.style == "median_iqr") {
    # Compute summary statistics with 1.5x IQR whiskers
    res_summary <- res %>%
      group_by(Model, Label) %>%
      summarise(
        q1 = quantile(!!sym(var), 0.25),
        med = median(!!sym(var)),
        q3 = quantile(!!sym(var), 0.75),
        iqr = q3 - q1,
        .groups = "drop"
      ) %>%
      mutate(
        x = as.numeric(factor(Model)),
        ymin = if (show.whisker) q1 - 1.5 * iqr else NA_real_,
        ymax = if (show.whisker) q3 + 1.5 * iqr else NA_real_
      )
    
    # Clip whiskers to the observed range
    res_limits <- res %>%
      group_by(Model, Label) %>%
      summarise(
        min_val = min(!!sym(var), na.rm = TRUE),
        max_val = max(!!sym(var), na.rm = TRUE),
        .groups = "drop"
      )
    
    res_summary <- left_join(res_summary, res_limits, by = c("Model", "Label")) %>%
      mutate(
        ymin = pmax(ymin, min_val),
        ymax = pmin(ymax, max_val)
      )
    
    p <- ggplot(res_summary, aes(x = x, color = Model)) +
      # Whiskers (1.5 x IQR, clipped to observed range)
      {if (show.whisker) geom_segment(aes(x = x, xend = x, y = ymin, yend = q1)) } +
      {if (show.whisker) geom_segment(aes(x = x, xend = x, y = q3, yend = ymax)) } +
      
      # IQR box (no fill)
      geom_rect(aes(xmin = x - 0.3, xmax = x + 0.3, ymin = q1, ymax = q3, col = Model),
                fill = NA, linewidth = 0.8) +
      
      # Median line
      geom_segment(aes(x = x - 0.3, xend = x + 0.3, y = med, yend = med, color = Model),
                   linewidth = 0.8) +
      
      scale_x_continuous(breaks = res_summary$x, labels = res_summary$Model) +
      facet_grid(Label ~ ., scales = "free") +
      scale_color_viridis_d(option = col.opt) +
      ggtitle(paste0(ifelse(is.null(base.model), paste0("Total ", var), paste0("Relative Total ", var, " vs ", base.model)),
                     ": Last ", use.n.years, " Years")) +
      ylab(ifelse(is.null(base.model), paste0("Total ", var), paste0("Relative Total ", var, " Difference"))) +
      xlab("Model") + 
      theme_bw()
  } else {
    stop("Unknown plot.style. Choose 'boxplot' or 'median_iqr'.")
  }
  
  plot_name <- paste0("Total_",var, ifelse(is.null(base.model), "", "_Relative"), "_last_", use.n.years, "_years.PNG")
  
  # Create the new subfolder if it doesn't exist
  new_sub_dir <- file.path(main.dir, sub.dir, "Performance_Boxplot")
  
  if (!file.exists(new_sub_dir)){
    dir.create(new_sub_dir)
  }
  
  # Save the figure inside the new subfolder
  ggsave(file.path(new_sub_dir, plot_name), plot = p, width = width, height = height, dpi = dpi)
  
  return(p)
}


plot_total_catch_performance2 <- function(mods, is.nsim, main.dir, sub.dir, var = "Catch",
                                          width = 10, height = 7, dpi = 300, col.opt = "D",
                                          method = NULL,
                                          outlier.opt = NA,
                                          plot.style = "median_iqr", 
                                          show.whisker = TRUE,
                                          use.n.years = NULL,
                                          start.years = NULL,
                                          new_model_names = NULL,
                                          base.model = NULL) {
  
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(rlang)
  
  if (is.null(use.n.years)) {
    cat("\nuse.n.years is not specified, so default (5 years) is used here!\n")
    use.n.years <- 5
  }
  
  if (is.null(start.years)) {
    cat("\nstart.years is not specified, so default (1st year in historical period) is used here!\n")
    start.years <- 1
  }
  
  res <- NULL
  
  if (!is.nsim) {
    Years <- mods[[1]]$om$years
    res <- lapply(seq_along(mods), function(i) {
      tmp <- data.frame(
        Catch = mods[[i]]$om$rep$pred_catch,
        Model = paste0("Model", i),
        Year = Years,
        Realization = 1
      )
      start_idx <- start.years
      end_idx <- min(start.years + use.n.years - 1, nrow(tmp))
      tmp[start_idx:end_idx, ]
    }) %>% bind_rows()
  } else {
    Years <- mods[[1]][[1]]$om$years
    res <- lapply(seq_along(mods), function(r) {
      lapply(seq_along(mods[[r]]), function(m) {
        tmp <- data.frame(
          Catch = mods[[r]][[m]]$om$rep$pred_catch,
          Model = paste0("Model", m),
          Year = Years,
          Realization = r
        )
        start_idx <- start.years
        end_idx <- min(start.years + use.n.years - 1, nrow(tmp))
        tmp[start_idx:end_idx, ]
      }) %>% bind_rows()
    }) %>% bind_rows()
  }
  
  # Add the last column: total Catch by Year and Realization
  res <- res %>%
    rowwise() %>%
    mutate(Catch_Global = sum(c_across(starts_with("Catch")), na.rm = TRUE)) %>%
    ungroup()
  
  # Rename models if requested
  if (!is.null(new_model_names)) {
    if (length(new_model_names) != length(unique(res$Model))) {
      stop("Length of new_model_names must match the number of models.")
    }
    res$Model <- factor(res$Model,
                        levels = paste0("Model", seq_along(new_model_names)),
                        labels = new_model_names)
    # if (!is.null(base.model)) base.model <- new_model_names[base.model]
    if (!is.null(base.model)) {
      if (!(base.model %in% new_model_names)) {
        warning("base.model does not match any of the new_model_names.")
      }
    }
  }
  
  # Pivot longer if multiple types
  res <- pivot_longer(res, cols = starts_with(var), names_to = "Label", values_to = var)
  
  # Then sum over years (use.n.years) by Model, Realization, and Label
  res <- res %>%
    group_by(Model, Realization, Label) %>%
    summarise(!!var := sum(!!sym(var), na.rm = TRUE), .groups = "drop")
  
  # Relative difference from base.model
  if (!is.null(base.model)) {
    base_df <- res %>%
      filter(Model == base.model) %>%
      rename(base_val = !!sym(var)) %>%
      select(Realization, Label, base_val)
    
    res <- left_join(res, base_df, by = c("Realization", "Label")) %>%
      mutate(!!var := (!!sym(var)) / base_val - 1)
  }
  
  # Apply mean or median summarization if method is specified
  if (!is.null(method)) {
    res <- res %>%
      group_by(Model, Realization, Label) %>%
      summarise(!!var := if (method == "mean") mean(!!sym(var)) else median(!!sym(var)),
                .groups = "drop")
  }
  
  # Plot
  if (plot.style == "boxplot") {
    p <- ggplot(res, aes(x = Model, y = Catch, color = Model)) +
      geom_boxplot(lwd = 0.8, outlier.shape = outlier.opt) +
      facet_grid(Label ~ ., scales = "free") +
      scale_color_viridis_d(option = col.opt) +
      ggtitle(paste0(ifelse(is.null(base.model), paste0("Total ", var), paste0("Relative Total ", var, " vs ", base.model)),
                     ": Years ", start.years, " to ", start.years + use.n.years - 1)) +
      ylab(ifelse(is.null(base.model), paste0("Total ", var), paste0("Relative Total ", var, " Difference"))) +
      xlab("Model") + 
      theme_bw()
  } else if (plot.style == "median_iqr") {
    # Compute summary statistics with 1.5x IQR whiskers
    res_summary <- res %>%
      group_by(Model, Label) %>%
      summarise(
        q1 = quantile(!!sym(var), 0.25),
        med = median(!!sym(var)),
        q3 = quantile(!!sym(var), 0.75),
        iqr = q3 - q1,
        .groups = "drop"
      ) %>%
      mutate(
        x = as.numeric(factor(Model)),
        ymin = if (show.whisker) q1 - 1.5 * iqr else NA_real_,
        ymax = if (show.whisker) q3 + 1.5 * iqr else NA_real_
      )
    
    # Clip whiskers to the observed range
    res_limits <- res %>%
      group_by(Model, Label) %>%
      summarise(
        min_val = min(!!sym(var), na.rm = TRUE),
        max_val = max(!!sym(var), na.rm = TRUE),
        .groups = "drop"
      )
    
    res_summary <- left_join(res_summary, res_limits, by = c("Model", "Label")) %>%
      mutate(
        ymin = pmax(ymin, min_val),
        ymax = pmin(ymax, max_val)
      )
    
    p <- ggplot(res_summary, aes(x = x, color = Model)) +
      # Whiskers (1.5 x IQR, clipped to observed range)
      {if (show.whisker) geom_segment(aes(x = x, xend = x, y = ymin, yend = q1)) } +
      {if (show.whisker) geom_segment(aes(x = x, xend = x, y = q3, yend = ymax)) } +
      
      # IQR box (no fill)
      geom_rect(aes(xmin = x - 0.3, xmax = x + 0.3, ymin = q1, ymax = q3, col = Model),
                fill = NA, linewidth = 0.8) +
      
      # Median line
      geom_segment(aes(x = x - 0.3, xend = x + 0.3, y = med, yend = med, color = Model),
                   linewidth = 0.8) +
      
      scale_x_continuous(breaks = res_summary$x, labels = res_summary$Model) +
      facet_grid(Label ~ ., scales = "free") +
      scale_color_viridis_d(option = col.opt) +
      ggtitle(paste0(ifelse(is.null(base.model), paste0("Total ", var), paste0("Relative Total_", var, " vs ", base.model)),
                     ": Years ", start.years, " to ", start.years + use.n.years - 1)) +
      ylab(ifelse(is.null(base.model), paste0("Total ", var), paste0("Relative Total ", var, " Difference"))) +
      xlab("Model") + 
      theme_bw()
  } else {
    stop("Unknown plot.style. Choose 'boxplot' or 'median_iqr'.")
  }
  
  # Save the plot
  plot_name <- paste0("Total_", var, ifelse(is.null(base.model), "", "_Relative"),"_first_", use.n.years, "_years.PNG")
  
  # Create the new subfolder if it doesn't exist
  new_sub_dir <- file.path(main.dir, sub.dir, "Performance_Boxplot")
  
  if (!file.exists(new_sub_dir)){
    dir.create(new_sub_dir)
  }
  
  # Save the figure inside the new subfolder
  ggsave(file.path(new_sub_dir, plot_name), plot = p, width = width, height = height, dpi = dpi)
  
  return(p)
}

# plot_model_performance_radar3 <- function(mods, is.nsim, main.dir, sub.dir, 
#                                           width = 10, height = 10, dpi = 300, col.opt = "D",
#                                           method = NULL,
#                                           use.n.years.first = 5,
#                                           use.n.years.last = 5,
#                                           start.years = 1, 
#                                           new_model_names = NULL) {
#   
#   library(dplyr)
#   library(tidyr)
#   library(fmsb)
#   library(viridis)
#   library(viridisLite)
#   
#   # Helper: calculate Average Annual Catch Variation (AACV)
#   calculate_aacv <- function(values) {
#     if (!is.numeric(values)) {
#       stop("Input must be a numeric vector.")
#     }
#     diffs <- abs(diff(values))
#     aacv <- sum(diffs) / sum(values[-length(values)])
#     return(aacv)
#   }
#   
#   if (is.nsim) {
#     
#     n_models <- length(mods[[1]])
#     n_reps <- length(mods)
#     
#     results <- list()
#     
#     for (r in seq_len(n_reps)) {
#       tmp <- data.frame(Model = paste0("Model", seq_len(n_models)))
#       
#       for (m in seq_len(n_models)) {
#         rep <- mods[[r]][[m]]$om$rep
#         
#         # Catch and SSB as rowSums
#         catch_ts <- rowSums(rep$pred_catch)
#         ssb_ts <- rowSums(rep$SSB)
#         
#         # Global Fbar
#         n_fleets <- mods[[r]][[m]]$om$input$data$n_fleets[1]
#         n_regions <- mods[[r]][[m]]$om$input$data$n_regions[1]
#         
#         # Catch variation
#         catch_aacv <- calculate_aacv(rowSums(rep$pred_catch))
#         
#         # SSB variation
#         ssb_aacv <- calculate_aacv(rowSums(rep$SSB))
#         
#         # Fbar variation (global is last col)
#         fbar_aacv <- calculate_aacv(rep$Fbar[, ncol(rep$Fbar)])
#         
#         # First and last period means
#         if(is.null(method)) method = "median"
#         if (method == "median") {
#           tmp$Catch_first[m] <- median(catch_ts[start.years:(start.years + use.n.years.first - 1)])
#           tmp$SSB_first[m] <- median(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
#           tmp$Catch_last[m] <- median(tail(catch_ts, use.n.years.last))
#           tmp$SSB_last[m] <- median(tail(ssb_ts, use.n.years.last))
#           tmp$catch_aacv[m] <- catch_aacv
#           tmp$ssb_aacv[m] <- ssb_aacv
#           tmp$fbar_aacv[m] <- fbar_aacv
#         } else if (method == "mean") {
#           tmp$Catch_first[m] <- mean(catch_ts[start.years:(start.years + use.n.years.first - 1)])
#           tmp$SSB_first[m] <- mean(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
#           tmp$Catch_last[m] <- mean(tail(catch_ts, use.n.years.last))
#           tmp$SSB_last[m] <- mean(tail(ssb_ts, use.n.years.last))
#           tmp$catch_aacv[m] <- catch_aacv
#           tmp$ssb_aacv[m] <- ssb_aacv
#           tmp$fbar_aacv[m] <- fbar_aacv
#         }
#         
#         if (is.null(mods[[1]][[1]]$om$rep$log_SSB_FXSPR)) {
#           message("Biological Reference Point has not been calculated internally!")
#           return(invisible(NULL))
#         }
#         
#         tmp1 <- rep$Fbar[, ncol(rep$Fbar)]
#         tmp2 <- exp(rep$log_Fbar_XSPR[, ncol(rep$Fbar)])
#         temp <- tmp1/tmp2
#         temp1 <- temp[start.years:(start.years + use.n.years.first - 1)]
#         temp2 <- tail(temp, use.n.years.last)
#         tmp$prob_first[m] <- mean(temp1 > 1, na.rm = TRUE)
#         tmp$prob_last[m] <- mean(temp2 > 1, na.rm = TRUE)
#       }
#       
#       # Normalize: higher is better except for Fbar
#       for (v in c("Catch_first", "SSB_first", "Catch_last", "SSB_last")) {
#         range_val <- max(tmp[[v]]) - min(tmp[[v]])
#         if (range_val == 0) {
#           tmp[[v]] <- 100
#         } else {
#           tmp[[v]] <- 100 * (tmp[[v]] - min(tmp[[v]])) / range_val
#         }
#       }
#       for (v in c("catch_aacv", "ssb_aacv", "fbar_aacv","prob_first","prob_last")) {
#         range_val <- max(tmp[[v]]) - min(tmp[[v]])
#         if (range_val == 0) {
#           tmp[[v]] <- 100
#         } else {
#           norm_f <- (tmp[[v]] - min(tmp[[v]])) / range_val
#           tmp[[v]] <- 100 * (1 - norm_f)
#         }
#       }
#       
#       results[[r]] <- tmp
#     }
#     
#     # Combine and compute median across realizations
#     combined <- bind_rows(results, .id = "Realization")
#     scores_median <- combined %>%
#       group_by(Model) %>%
#       summarise(across(-Realization, median), .groups = "drop")
#     
#     # Optional renaming
#     if (!is.null(new_model_names)) {
#       if (length(new_model_names) != length(unique(scores_median$Model))) {
#         stop("Length of new_model_names must match the number of models.")
#       }
#       scores_median$Model <- factor(scores_median$Model,
#                                     levels = paste0("Model", seq_along(new_model_names)),
#                                     labels = new_model_names)
#     }
#     
#     # Radar chart prep
#     plot_df <- as.data.frame(scores_median)         # Ensure it's a base data.frame
#     rownames(plot_df) <- plot_df$Model              # Set row names safely
#     plot_df$Model <- NULL                           # Drop Model column
#     plot_df <- as.data.frame(t(plot_df))            # Transpose for radar chart
#     
#     # Add min/max rows for radar scale
#     plot_df <- rbind(rep(100, ncol(plot_df)), rep(0, ncol(plot_df)), plot_df)
#     
#   } else {
#     
#     n_models <- length(mods)
#     results <- list()
#     
#     tmp <- data.frame(Model = paste0("Model", seq_len(n_models)))
#     
#     for (m in seq_len(n_models)) {
#       
#       rep <- mods[[m]]$om$rep
#       
#       # Catch and SSB as rowSums
#       catch_ts <- rowSums(rep$pred_catch)
#       ssb_ts <- rowSums(rep$SSB)
#       
#       # Global Fbar
#       n_fleets <- mods[[m]]$om$input$data$n_fleets[1]
#       n_regions <- mods[[m]]$om$input$data$n_regions[1]
#       
#       # Catch variation
#       catch_aacv <- calculate_aacv(rowSums(rep$pred_catch))
#       
#       # SSB variation
#       ssb_aacv <- calculate_aacv(rowSums(rep$SSB))
#       
#       # Fbar variation (global is last col)
#       fbar_aacv <- calculate_aacv(rep$Fbar[, ncol(rep$Fbar)])
#       
#       if(is.null(method)) method = "median"
#       # First and last period means
#       if (method == "median") {
#         tmp$Catch_first[m] <- median(catch_ts[start.years:(start.years + use.n.years.first - 1)])
#         tmp$SSB_first[m] <- median(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
#         tmp$Catch_last[m] <- median(tail(catch_ts, use.n.years.last))
#         tmp$SSB_last[m] <- median(tail(ssb_ts, use.n.years.last))
#         tmp$catch_aacv[m] <- catch_aacv
#         tmp$ssb_aacv[m] <- ssb_aacv
#         tmp$fbar_aacv[m] <- fbar_aacv
#         
#       } else if (method == "mean") {
#         tmp$Catch_first[m] <- mean(catch_ts[start.years:(start.years + use.n.years.first - 1)])
#         tmp$SSB_first[m] <- mean(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
#         tmp$Catch_last[m] <- mean(tail(catch_ts, use.n.years.last))
#         tmp$SSB_last[m] <- mean(tail(ssb_ts, use.n.years.last))
#         tmp$catch_aacv[m] <- catch_aacv
#         tmp$ssb_aacv[m] <- ssb_aacv
#         tmp$fbar_aacv[m] <- fbar_aacv
#       }
#     }
#     
#     if (is.null(rep$log_SSB_FXSPR)) {
#       message("Biological Reference Point has not been calculated internally!")
#       return(invisible(NULL))
#     }
#     
#     tmp1 <- rep$Fbar[, ncol(rep$Fbar)]
#     tmp2 <- exp(rep$log_Fbar_XSPR[, ncol(rep$Fbar)])
#     temp <- tmp1/tmp2
#     temp1 <- temp[start.years:(start.years + use.n.years.first - 1)]
#     temp2 <- tail(temp, use.n.years.last)
#     tmp$prob_first[m] <- mean(temp1 > 1, na.rm = TRUE)
#     tmp$prob_last[m] <- mean(temp2 > 1, na.rm = TRUE)
#     
#     # Normalize: higher is better except for Fbar
#     for (v in c("Catch_first", "SSB_first", "Catch_last", "SSB_last")) {
#       range_val <- max(tmp[[v]]) - min(tmp[[v]])
#       if (range_val == 0) {
#         tmp[[v]] <- 100  # All values are identical; assign full score
#       } else {
#         tmp[[v]] <- 100 * (tmp[[v]] - min(tmp[[v]])) / range_val
#       }
#     }
#     
#     for (v in c("catch_aacv", "ssb_aacv", "fbar_aacv","prob_first","prob_last")) {
#       range_val <- max(tmp[[v]]) - min(tmp[[v]])
#       if (range_val == 0) {
#         tmp[[v]] <- 100  # All values identical; assign full score
#       } else {
#         norm_f <- (tmp[[v]] - min(tmp[[v]])) / range_val
#         tmp[[v]] <- 100 * (1 - norm_f)
#       }
#     }
#     
#     results[[1]] <- tmp
#     
#     # Combine and compute median across realizations
#     combined <- bind_rows(results, .id = "Realization")
#     scores_median <- combined %>%
#       group_by(Model) %>%
#       summarise(across(-Realization, median), .groups = "drop")
#     
#     # Optional renaming
#     if (!is.null(new_model_names)) {
#       if (length(new_model_names) != length(unique(scores_median$Model))) {
#         stop("Length of new_model_names must match the number of models.")
#       }
#       scores_median$Model <- factor(scores_median$Model,
#                                     levels = paste0("Model", seq_along(new_model_names)),
#                                     labels = new_model_names)
#     }
#     
#     # Radar chart prep
#     plot_df <- as.data.frame(scores_median)         # Ensure it's a base data.frame
#     rownames(plot_df) <- plot_df$Model              # Set row names safely
#     plot_df$Model <- NULL                           # Drop Model column
#     plot_df <- as.data.frame(t(plot_df))            # Transpose for radar chart
#     
#     # Add min/max rows for radar scale
#     plot_df <- rbind(rep(100, ncol(plot_df)), rep(0, ncol(plot_df)), plot_df)
#     
#   }
#   
#   if (ncol(plot_df) < 3) {
#     message("Radar chart needs at least 3 models. Showing barplot instead.")
#     return(invisible(NULL))
#   }
#   
#   # Plot without any further changes
#   colors <- viridisLite::viridis(n = nrow(plot_df) - 2, option = col.opt) # -2 for min/max rows
#   
#   # Create the new subfolder if it doesn't exist
#   new_sub_dir <- file.path(main.dir, sub.dir, "Radar_Holistic_Plot")
#   
#   if (!file.exists(new_sub_dir)){
#     dir.create(new_sub_dir)
#   }
#   
#   rownames(plot_df)
#   # Correct way to assign expressions to rownames
#   my_legend_labels <- c(expression(Catch[ST]),
#                         expression(SSB[ST]),
#                         expression(Catch[LT]),
#                         expression(SSB[LT]),
#                         expression(Catch[AAV]),
#                         expression(SSB[AAV]),
#                         expression(F[AAV]),
#                         expression(P(F[ST] > F[MSY])),
#                         expression(P(F[LT] > F[MSY])))
#   
#   # Save to PNG
#   output_file <- file.path(file.path(main.dir, sub.dir, "Radar_Holistic_Plot", "model_performance_radar3.png"))
#   png(filename = output_file, width = width, height = height, units = "in", res = dpi)
#   radarchart(plot_df,
#              axistype = 4,
#              pcol = colors,
#              plwd = 3,
#              # pty = 1:length(colors),
#              plty = 1:length(colors),
#              cglcol = "grey80",
#              cglty = 1,
#              axislabcol = "grey30",
#              vlcex = 0.9)
#   legend(x = 0.9, y = 0.3, legend = my_legend_labels, col = colors,
#          lty = 1:length(colors), 
#          # pch = 1:length(colors),
#          lwd = 3, cex = 0.9, y.intersp = 1.5)
#   dev.off()
#   
#   # Also plot to screen (inline)
#   op <- par(mar = c(1, 1, 3, 1))  # NEW: minimize margins
#   radarchart(plot_df,
#              axistype = 4,
#              pcol = colors,
#              plwd = 3,
#              # pty = 1:length(colors),
#              plty = 1:length(colors),
#              cglcol = "grey80",
#              cglty = 1,
#              axislabcol = "grey30",
#              vlcex = 0.9)
#   legend(x = 0.9, y = 0.3, legend = my_legend_labels, col = colors,
#          lty = 1:length(colors), 
#          # pch = 1:length(colors),
#          lwd = 3, cex = 0.9, y.intersp = 1.5)
#   on.exit(par(op))  # restore after plotting
# }


plot_model_performance_bar3 <- function(mods, is.nsim,
                                        main.dir = ".",
                                        sub.dir = ".",
                                        new_model_names = NULL,
                                        width = 12, height = 8, dpi = 300,
                                        col.opt = "D",
                                        method = NULL,
                                        use.n.years.first = 5,
                                        use.n.years.last = 5,
                                        start.years = 1) {

  library(dplyr)
  library(tidyr)
  library(fmsb)
  library(viridis)
  library(viridisLite)

  # Helper: calculate Average Annual Catch Variation (AACV)
  calculate_aacv <- function(values) {
    if (!is.numeric(values)) {
      stop("Input must be a numeric vector.")
    }
    diffs <- abs(diff(values))
    aacv <- sum(diffs) / sum(values[-length(values)])
    return(aacv)
  }

  if (is.nsim) {

    n_models <- length(mods[[1]])
    n_reps <- length(mods)

    results <- list()

    for (r in seq_len(n_reps)) {
      tmp <- data.frame(Model = paste0("Model", seq_len(n_models)))

      for (m in seq_len(n_models)) {
        rep <- mods[[r]][[m]]$om$rep

        # Catch and SSB as rowSums
        catch_ts <- rowSums(rep$pred_catch)
        ssb_ts <- rowSums(rep$SSB)

        # Global Fbar
        n_fleets <- mods[[r]][[m]]$om$input$data$n_fleets[1]
        n_regions <- mods[[r]][[m]]$om$input$data$n_regions[1]

        # Catch variation
        catch_aacv <- calculate_aacv(rowSums(rep$pred_catch))

        # SSB variation
        ssb_aacv <- calculate_aacv(rowSums(rep$SSB))

        # Fbar variation (global is last col)
        fbar_aacv <- calculate_aacv(rep$Fbar[, ncol(rep$Fbar)])

        # First and last period means
        if(is.null(method)) method = "median"
        if (method == "median") {
          tmp$Catch_first[m] <- median(catch_ts[start.years:(start.years + use.n.years.first - 1)])
          tmp$SSB_first[m] <- median(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
          tmp$Catch_last[m] <- median(tail(catch_ts, use.n.years.last))
          tmp$SSB_last[m] <- median(tail(ssb_ts, use.n.years.last))
          tmp$catch_aacv[m] <- catch_aacv
          tmp$ssb_aacv[m] <- ssb_aacv
          tmp$fbar_aacv[m] <- fbar_aacv
        } else if (method == "mean") {
          tmp$Catch_first[m] <- mean(catch_ts[start.years:(start.years + use.n.years.first - 1)])
          tmp$SSB_first[m] <- mean(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
          tmp$Catch_last[m] <- mean(tail(catch_ts, use.n.years.last))
          tmp$SSB_last[m] <- mean(tail(ssb_ts, use.n.years.last))
          tmp$catch_aacv[m] <- catch_aacv
          tmp$ssb_aacv[m] <- ssb_aacv
          tmp$fbar_aacv[m] <- fbar_aacv
        }

        if (is.null(mods[[1]][[1]]$om$rep$log_SSB_FXSPR)) {
          message("Biological Reference Point has not been calculated internally!")
          return(invisible(NULL))
        }

        tmp1 <- rep$Fbar[, ncol(rep$Fbar)]
        tmp2 <- exp(rep$log_Fbar_XSPR[, ncol(rep$Fbar)])
        temp <- tmp1/tmp2
        temp1 <- temp[start.years:(start.years + use.n.years.first - 1)]
        temp2 <- tail(temp, use.n.years.last)
        tmp$prob_first[m] <- mean(temp1 > 1, na.rm = TRUE)
        tmp$prob_last[m] <- mean(temp2 > 1, na.rm = TRUE)
      }

      # Normalize: higher is better except for Fbar
      for (v in c("Catch_first", "SSB_first", "Catch_last", "SSB_last")) {
        range_val <- max(tmp[[v]]) - min(tmp[[v]])
        if (range_val == 0) {
          tmp[[v]] <- 100
        } else {
          tmp[[v]] <- 100 * (tmp[[v]] - min(tmp[[v]])) / range_val
        }
      }
      for (v in c("catch_aacv", "ssb_aacv", "fbar_aacv","prob_first","prob_last")) {
        range_val <- max(tmp[[v]]) - min(tmp[[v]])
        if (range_val == 0) {
          tmp[[v]] <- 100
        } else {
          norm_f <- (tmp[[v]] - min(tmp[[v]])) / range_val
          tmp[[v]] <- 100 * (1 - norm_f)
        }
      }

      tmp$Realization <- r
      results[[r]] <- tmp
    }

    # Combine and compute median across realizations
    combined <- bind_rows(results, .id = "Realization")

    # Summarize across realizations
    scores_median <- combined %>%
      pivot_longer(-c(Model, Realization), names_to = "Metric", values_to = "Score") %>%
      group_by(Model, Metric) %>%
      summarize(
        Median = median(Score, na.rm = TRUE),
        Q1 = quantile(Score, 0.25, na.rm = TRUE),
        Q3 = quantile(Score, 0.75, na.rm = TRUE),
        .groups = "drop"
      )

    # Optional renaming
    if (!is.null(new_model_names)) {
      if (length(new_model_names) != length(unique(scores_median$Model))) {
        stop("Length of new_model_names must match the number of models.")
      }
      scores_median$Model <- factor(scores_median$Model,
                                    levels = paste0("Model", seq_along(new_model_names)),
                                    labels = new_model_names)
    }

  } else {

    n_models <- length(mods)
    results <- list()

    tmp <- data.frame(Model = paste0("Model", seq_len(n_models)))

    for (m in seq_len(n_models)) {

      rep <- mods[[m]]$om$rep

      # Catch and SSB as rowSums
      catch_ts <- rowSums(rep$pred_catch)
      ssb_ts <- rowSums(rep$SSB)

      # Global Fbar
      n_fleets <- mods[[m]]$om$input$data$n_fleets[1]
      n_regions <- mods[[m]]$om$input$data$n_regions[1]

      # Catch variation
      catch_aacv <- calculate_aacv(rowSums(rep$pred_catch))

      # SSB variation
      ssb_aacv <- calculate_aacv(rowSums(rep$SSB))

      # Fbar variation (global is last col)
      fbar_aacv <- calculate_aacv(rep$Fbar[, ncol(rep$Fbar)])

      if(is.null(method)) method = "median"
      # First and last period means
      if (method == "median") {
        tmp$Catch_first[m] <- median(catch_ts[start.years:(start.years + use.n.years.first - 1)])
        tmp$SSB_first[m] <- median(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
        tmp$Catch_last[m] <- median(tail(catch_ts, use.n.years.last))
        tmp$SSB_last[m] <- median(tail(ssb_ts, use.n.years.last))
        tmp$catch_aacv[m] <- catch_aacv
        tmp$ssb_aacv[m] <- ssb_aacv
        tmp$fbar_aacv[m] <- fbar_aacv

      } else if (method == "mean") {
        tmp$Catch_first[m] <- mean(catch_ts[start.years:(start.years + use.n.years.first - 1)])
        tmp$SSB_first[m] <- mean(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
        tmp$Catch_last[m] <- mean(tail(catch_ts, use.n.years.last))
        tmp$SSB_last[m] <- mean(tail(ssb_ts, use.n.years.last))
        tmp$catch_aacv[m] <- catch_aacv
        tmp$ssb_aacv[m] <- ssb_aacv
        tmp$fbar_aacv[m] <- fbar_aacv
      }
    }

    if (is.null(rep$log_SSB_FXSPR)) {
      message("Biological Reference Point has not been calculated internally!")
      return(invisible(NULL))
    }

    tmp1 <- rep$Fbar[, ncol(rep$Fbar)]
    tmp2 <- exp(rep$log_Fbar_XSPR[, ncol(rep$Fbar)])
    temp <- tmp1/tmp2
    temp1 <- temp[start.years:(start.years + use.n.years.first - 1)]
    temp2 <- tail(temp, use.n.years.last)
    tmp$prob_first[m] <- mean(temp1 > 1, na.rm = TRUE)
    tmp$prob_last[m] <- mean(temp2 > 1, na.rm = TRUE)

    # Normalize: higher is better except for Fbar
    for (v in c("Catch_first", "SSB_first", "Catch_last", "SSB_last")) {
      range_val <- max(tmp[[v]]) - min(tmp[[v]])
      if (range_val == 0) {
        tmp[[v]] <- 100  # All values are identical; assign full score
      } else {
        tmp[[v]] <- 100 * (tmp[[v]] - min(tmp[[v]])) / range_val
      }
    }

    for (v in c("catch_aacv", "ssb_aacv", "fbar_aacv","prob_first","prob_last")) {
      range_val <- max(tmp[[v]]) - min(tmp[[v]])
      if (range_val == 0) {
        tmp[[v]] <- 100  # All values identical; assign full score
      } else {
        norm_f <- (tmp[[v]] - min(tmp[[v]])) / range_val
        tmp[[v]] <- 100 * (1 - norm_f)
      }
    }

    tmp$Realization <- 1
    combined <- tmp

    scores_median <- tmp %>%
      pivot_longer(-c(Model, Realization), names_to = "Metric", values_to = "Median") %>%
      mutate(Q1 = Median, Q3 = Median)

    # Optional renaming
    if (!is.null(new_model_names)) {
      if (length(new_model_names) != length(unique(scores_median$Model))) {
        stop("Length of new_model_names must match the number of models.")
      }
      scores_median$Model <- factor(scores_median$Model,
                                    levels = paste0("Model", seq_along(new_model_names)),
                                    labels = new_model_names)
    }
  }

  # Colors
  n_metrics <- length(unique(scores_median$Metric))
  # Colors
  # Plot without any further changes
  colors <- viridisLite::viridis(n = n_metrics, option = col.opt) # -2 for min/max rows

  # Correct way to assign expressions to rownames
  my_legend_labels <- c(expression(Catch[ST]),
                        expression(SSB[ST]),
                        expression(Catch[LT]),
                        expression(SSB[LT]),
                        expression(Catch[AAV]),
                        expression(SSB[AAV]),
                        expression(F[AAV]),
                        expression(P(F[ST] > F[MSY])),
                        expression(P(F[LT] > F[MSY])))

  metric_levels <- c(
    "Catch_first",
    "SSB_first",
    "Catch_last",
    "SSB_last",
    "catch_aacv",
    "ssb_aacv",
    "fbar_aacv",
    "prob_first",
    "prob_last"
  )

  scores_median$Metric <- factor(scores_median$Metric,
                                 levels = metric_levels,
                                 labels = my_legend_labels)

  # Plot
  plot <- ggplot(scores_median, aes(x = Median, y = Model, fill = Metric, color = Metric)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.7) +
    geom_errorbar(aes(xmin = Q1, xmax = Q3), col = "black",
                  position = position_dodge(width = 0.9), width = 0.3, alpha = 0.3) +
    scale_fill_manual(values = colors,
                      labels = my_legend_labels) + # Add labels argument for fill legend
    scale_color_manual(values = colors,
                       labels = my_legend_labels) + # Add labels argument for color legend
    theme_bw() +
    labs(title = paste("Holistic Model Performance"),
         x = "Score (0-100)",
         y = "Estimation Model",
         fill = "Metric",
         color = "Metric") +
    theme(
      axis.text = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      plot.title = element_text(size = 18, face = "bold")
    )

  # Print and save
  print(plot)

  # Create the new subfolder if it doesn't exist
  new_sub_dir <- file.path(main.dir, sub.dir, "Holistic_Bar_Plot")

  if (!file.exists(new_sub_dir)){
    dir.create(new_sub_dir)
  }

  ggsave(filename = file.path(new_sub_dir, "Overall_Performance3.png"),
         plot = plot, width = width, height = height, dpi = dpi)
}


plot_model_performance_bar4 <- function(mods, is.nsim,
                                        main.dir = ".",
                                        sub.dir = ".",
                                        new_model_names = NULL,
                                        width = 12, height = 8, dpi = 300,
                                        col.opt = "D",
                                        method = NULL,
                                        use.n.years.first = 5,
                                        use.n.years.last = 5,
                                        start.years = 1) {
  
  library(dplyr)
  library(tidyr)
  library(fmsb)
  library(viridis)
  library(viridisLite)
  
  # Helper: calculate Average Annual Catch Variation (AACV)
  calculate_aacv <- function(values) {
    if (!is.numeric(values)) {
      stop("Input must be a numeric vector.")
    }
    diffs <- abs(diff(values))
    aacv <- sum(diffs) / sum(values[-length(values)])
    return(aacv)
  }
  
  if (is.nsim) {
    
    n_models <- length(mods[[1]])
    n_reps <- length(mods)
    
    results <- list()
    
    for (r in seq_len(n_reps)) {
      tmp <- data.frame(Model = paste0("Model", seq_len(n_models)))
      
      for (m in seq_len(n_models)) {
        rep <- mods[[r]][[m]]$om$rep
        
        # Catch and SSB as rowSums
        catch_ts <- rowSums(rep$pred_catch)
        ssb_ts <- rowSums(rep$SSB)
        
        # Global Fbar
        n_fleets <- mods[[r]][[m]]$om$input$data$n_fleets[1]
        n_regions <- mods[[r]][[m]]$om$input$data$n_regions[1]
        
        # Catch variation
        catch_aacv <- calculate_aacv(rowSums(rep$pred_catch))
        
        # SSB variation
        ssb_aacv <- calculate_aacv(rowSums(rep$SSB))
        
        # Fbar variation (global is last col)
        fbar_aacv <- calculate_aacv(rep$Fbar[, ncol(rep$Fbar)])
        
        # First and last period means
        if(is.null(method)) method = "median"
        if (method == "median") {
          tmp$Catch_first[m] <- median(catch_ts[start.years:(start.years + use.n.years.first - 1)])
          tmp$SSB_first[m] <- median(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
          tmp$Catch_last[m] <- median(tail(catch_ts, use.n.years.last))
          tmp$SSB_last[m] <- median(tail(ssb_ts, use.n.years.last))
          tmp$catch_aacv[m] <- catch_aacv
          tmp$ssb_aacv[m] <- ssb_aacv
          tmp$fbar_aacv[m] <- fbar_aacv
        } else if (method == "mean") {
          tmp$Catch_first[m] <- mean(catch_ts[start.years:(start.years + use.n.years.first - 1)])
          tmp$SSB_first[m] <- mean(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
          tmp$Catch_last[m] <- mean(tail(catch_ts, use.n.years.last))
          tmp$SSB_last[m] <- mean(tail(ssb_ts, use.n.years.last))
          tmp$catch_aacv[m] <- catch_aacv
          tmp$ssb_aacv[m] <- ssb_aacv
          tmp$fbar_aacv[m] <- fbar_aacv
        }
        
        if (is.null(mods[[1]][[1]]$om$rep$log_SSB_FXSPR)) {
          message("Biological Reference Point has not been calculated internally!")
          return(invisible(NULL))
        }
        
        tmp1 <- rep$Fbar[, ncol(rep$Fbar)]
        tmp2 <- exp(rep$log_Fbar_XSPR[, ncol(rep$Fbar)])
        temp <- tmp1/tmp2
        temp1 <- temp[start.years:(start.years + use.n.years.first - 1)]
        temp2 <- tail(temp, use.n.years.last)
        tmp$prob_first[m] <- mean(temp1 > 1, na.rm = TRUE)
        tmp$prob_last[m] <- mean(temp2 > 1, na.rm = TRUE)
      }
      
      # Normalize: higher is better except for Fbar
      for (v in c("Catch_first", "SSB_first", "Catch_last", "SSB_last")) {
        range_val <- max(tmp[[v]]) - min(tmp[[v]])
        if (range_val == 0) {
          tmp[[v]] <- 100
        } else {
          tmp[[v]] <- 100 * (tmp[[v]] - min(tmp[[v]])) / range_val
        }
      }
      for (v in c("catch_aacv", "ssb_aacv", "fbar_aacv","prob_first","prob_last")) {
        range_val <- max(tmp[[v]]) - min(tmp[[v]])
        if (range_val == 0) {
          tmp[[v]] <- 100
        } else {
          norm_f <- (tmp[[v]] - min(tmp[[v]])) / range_val
          tmp[[v]] <- 100 * (1 - norm_f)
        }
      }
      
      results[[r]] <- tmp
    }
    
    # Combine and compute median across realizations
    combined <- bind_rows(results, .id = "Realization")
    
    # Summarize across realizations
    scores_median <- combined %>%
      pivot_longer(-c(Model, Realization), names_to = "Metric", values_to = "Score") %>%
      group_by(Model, Metric) %>%
      summarize(
        Median = median(Score, na.rm = TRUE),
        Q1 = quantile(Score, 0.25, na.rm = TRUE),
        Q3 = quantile(Score, 0.75, na.rm = TRUE),
        .groups = "drop"
      )
    
    # Optional renaming
    if (!is.null(new_model_names)) {
      if (length(new_model_names) != length(unique(scores_median$Model))) {
        stop("Length of new_model_names must match the number of models.")
      }
      scores_median$Model <- factor(scores_median$Model,
                                    levels = paste0("Model", seq_along(new_model_names)),
                                    labels = new_model_names)
    }
    
  } 
  else {
    
    n_models <- length(mods)
    results <- list()
    
    tmp <- data.frame(Model = paste0("Model", seq_len(n_models)))
    
    for (m in seq_len(n_models)) {
      
      rep <- mods[[m]]$om$rep
      
      # Catch and SSB as rowSums
      catch_ts <- rowSums(rep$pred_catch)
      ssb_ts <- rowSums(rep$SSB)
      
      # Global Fbar
      n_fleets <- mods[[m]]$om$input$data$n_fleets[1]
      n_regions <- mods[[m]]$om$input$data$n_regions[1]
      
      # Catch variation
      catch_aacv <- calculate_aacv(rowSums(rep$pred_catch))
      
      # SSB variation
      ssb_aacv <- calculate_aacv(rowSums(rep$SSB))
      
      # Fbar variation (global is last col)
      fbar_aacv <- calculate_aacv(rep$Fbar[, ncol(rep$Fbar)])
      
      if(is.null(method)) method = "median"
      # First and last period means
      if (method == "median") {
        tmp$Catch_first[m] <- median(catch_ts[start.years:(start.years + use.n.years.first - 1)])
        tmp$SSB_first[m] <- median(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
        tmp$Catch_last[m] <- median(tail(catch_ts, use.n.years.last))
        tmp$SSB_last[m] <- median(tail(ssb_ts, use.n.years.last))
        tmp$catch_aacv[m] <- catch_aacv
        tmp$ssb_aacv[m] <- ssb_aacv
        tmp$fbar_aacv[m] <- fbar_aacv
        
      } else if (method == "mean") {
        tmp$Catch_first[m] <- mean(catch_ts[start.years:(start.years + use.n.years.first - 1)])
        tmp$SSB_first[m] <- mean(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
        tmp$Catch_last[m] <- mean(tail(catch_ts, use.n.years.last))
        tmp$SSB_last[m] <- mean(tail(ssb_ts, use.n.years.last))
        tmp$catch_aacv[m] <- catch_aacv
        tmp$ssb_aacv[m] <- ssb_aacv
        tmp$fbar_aacv[m] <- fbar_aacv
      }
    }
    
    if (is.null(rep$log_SSB_FXSPR)) {
      message("Biological Reference Point has not been calculated internally!")
      return(invisible(NULL))
    }
    
    tmp1 <- rep$Fbar[, ncol(rep$Fbar)]
    tmp2 <- exp(rep$log_Fbar_XSPR[, ncol(rep$Fbar)])
    temp <- tmp1/tmp2
    temp1 <- temp[start.years:(start.years + use.n.years.first - 1)]
    temp2 <- tail(temp, use.n.years.last)
    tmp$prob_first[m] <- mean(temp1 > 1, na.rm = TRUE)
    tmp$prob_last[m] <- mean(temp2 > 1, na.rm = TRUE)
    
    # Normalize: higher is better except for Fbar
    for (v in c("Catch_first", "SSB_first", "Catch_last", "SSB_last")) {
      range_val <- max(tmp[[v]]) - min(tmp[[v]])
      if (range_val == 0) {
        tmp[[v]] <- 100  # All values are identical; assign full score
      } else {
        tmp[[v]] <- 100 * (tmp[[v]] - min(tmp[[v]])) / range_val
      }
    }
    
    for (v in c("catch_aacv", "ssb_aacv", "fbar_aacv","prob_first","prob_last")) {
      range_val <- max(tmp[[v]]) - min(tmp[[v]])
      if (range_val == 0) {
        tmp[[v]] <- 100  # All values identical; assign full score
      } else {
        norm_f <- (tmp[[v]] - min(tmp[[v]])) / range_val
        tmp[[v]] <- 100 * (1 - norm_f)
      }
    }
    
    tmp$Realization <- 1
    combined <- tmp
    
    scores_median <- tmp %>%
      pivot_longer(-c(Model, Realization), names_to = "Metric", values_to = "Median") %>%
      mutate(Q1 = Median, Q3 = Median)
    
    # Optional renaming
    if (!is.null(new_model_names)) {
      if (length(new_model_names) != length(unique(scores_median$Model))) {
        stop("Length of new_model_names must match the number of models.")
      }
      scores_median$Model <- factor(scores_median$Model,
                                    levels = paste0("Model", seq_along(new_model_names)),
                                    labels = new_model_names)
    }
    
  }
  
  # Colors
  n_metrics <- length(unique(scores_median$Model))
  
  # Colors
  # Plot without any further changes
  colors <- viridisLite::viridis(n = n_metrics, option = col.opt) 

  # Correct way to assign expressions to rownames
  my_legend_labels <- c(expression(bold(Catch[ST])),
                        expression(bold(SSB[ST])),
                        expression(bold(Catch[LT])),
                        expression(bold(SSB[LT])),
                        expression(bold(Catch[AAV])),
                        expression(bold(SSB[AAV])),
                        expression(bold(F[AAV])),
                        expression(bold(P(F[ST] > F[MSY]))),
                        expression(bold(P(F[LT] > F[MSY]))))
  
  metric_levels <- c(
    "Catch_first",
    "SSB_first",
    "Catch_last",
    "SSB_last",
    "catch_aacv",
    "ssb_aacv",
    "fbar_aacv",
    "prob_first",
    "prob_last"
  )
  
  scores_median$Metric <- factor(scores_median$Metric,
                                 levels = metric_levels,
                                 labels = my_legend_labels)
  
  # Plot
  plot <- ggplot(scores_median, aes(x = Median, y = Metric, fill = Model, color = Model)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.7) +
    geom_errorbar(aes(xmin = Q1, xmax = Q3), col = "black",
                  position = position_dodge(width = 0.9), width = 0.3, alpha = 0.3) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    theme_bw() +
    labs(title = paste("Holistic Model Performance"),
         x = "Score (0-100)",
         y = "Metric",
         fill = "Model",
         color = "Model") +
    scale_y_discrete(labels = my_legend_labels,
                     limits = rev(levels(scores_median$Metric))) +
    theme(
      axis.text = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 16, face = "bold"),
      legend.text = element_text(size = 14),
      legend.title = element_text(size = 16),
      plot.title = element_text(size = 18, face = "bold")
    )
  
  # Print and save
  print(plot)
  
  # Create the new subfolder if it doesn't exist
  new_sub_dir <- file.path(main.dir, sub.dir, "Holistic_Bar_Plot")
  
  if (!file.exists(new_sub_dir)){
    dir.create(new_sub_dir)
  }
  
  ggsave(filename = file.path(new_sub_dir, "Overall_Performance4.png"),
         plot = plot, width = width, height = height, dpi = dpi)
}


plot_model_performance_radar <- function(mods, is.nsim, main.dir, sub.dir, 
                                         width = 10, height = 10, dpi = 300, col.opt = "D",
                                         method = NULL,
                                         use.n.years.first = 5,
                                         use.n.years.last = 5,
                                         start.years = 1, 
                                         new_model_names = NULL) {
  
  library(dplyr)
  library(tidyr)
  library(fmsb)
  library(viridis)
  
  if (is.nsim) {
    
    n_models <- length(mods[[1]])
    n_reps <- length(mods)
    results <- list()
    
    for (r in seq_len(n_reps)) {
      tmp <- data.frame(Model = paste0("Model", seq_len(n_models)))
      
      for (m in seq_len(n_models)) {
        rep <- mods[[r]][[m]]$om$rep
        
        # Catch and SSB as rowSums
        catch_ts <- rowSums(rep$pred_catch)
        ssb_ts <- rowSums(rep$SSB)
        
        # Global Fbar
        n_fleets <- mods[[r]][[m]]$om$input$data$n_fleets[1]
        n_regions <- mods[[r]][[m]]$om$input$data$n_regions[1]
        fbar_ts <- rep$Fbar[, n_fleets + n_regions + 1]
        
        # First and last period means
        if(is.null(method)) method = "median"
        if (method == "median") {
          tmp$Catch_first[m] <- median(catch_ts[start.years:(start.years + use.n.years.first - 1)])
          tmp$SSB_first[m] <- median(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
          tmp$Fbar_first[m] <- median(fbar_ts[start.years:(start.years + use.n.years.first - 1)])
          tmp$Catch_last[m] <- median(tail(catch_ts, use.n.years.last))
          tmp$SSB_last[m] <- median(tail(ssb_ts, use.n.years.last))
          tmp$Fbar_last[m] <- median(tail(fbar_ts, use.n.years.last))
        } else if (method == "mean") {
          tmp$Catch_first[m] <- mean(catch_ts[start.years:(start.years + use.n.years.first - 1)])
          tmp$SSB_first[m] <- mean(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
          tmp$Fbar_first[m] <- mean(fbar_ts[start.years:(start.years + use.n.years.first - 1)])
          tmp$Catch_last[m] <- mean(tail(catch_ts, use.n.years.last))
          tmp$SSB_last[m] <- mean(tail(ssb_ts, use.n.years.last))
          tmp$Fbar_last[m] <- mean(tail(fbar_ts, use.n.years.last))
        } 
      }
      
      # Normalize: higher is better except for Fbar
      for (v in c("Catch_first", "SSB_first", "Catch_last", "SSB_last")) {
        range_val <- max(tmp[[v]]) - min(tmp[[v]])
        if (range_val == 0) {
          tmp[[v]] <- 100
        } else {
          tmp[[v]] <- 100 * (tmp[[v]] - min(tmp[[v]])) / range_val
        }
      }
      for (v in c("Fbar_first", "Fbar_last")) {
        range_val <- max(tmp[[v]]) - min(tmp[[v]])
        if (range_val == 0) {
          tmp[[v]] <- 100
        } else {
          norm_f <- (tmp[[v]] - min(tmp[[v]])) / range_val
          tmp[[v]] <- 100 * (1 - norm_f)
        }
      }
      
      results[[r]] <- tmp
    }
    
    # Combine and compute median across realizations
    combined <- bind_rows(results, .id = "Realization")
    scores_median <- combined %>%
      group_by(Model) %>%
      summarise(across(-Realization, median), .groups = "drop")
    
    # Optional renaming
    if (!is.null(new_model_names)) {
      if (length(new_model_names) != length(unique(scores_median$Model))) {
        stop("Length of new_model_names must match the number of models.")
      }
      scores_median$Model <- factor(scores_median$Model,
                                    levels = paste0("Model", seq_along(new_model_names)),
                                    labels = new_model_names)
    }
    
    # Radar chart prep
    plot_df <- as.data.frame(scores_median)         # Ensure it's a base data.frame
    rownames(plot_df) <- plot_df$Model              # Set row names safely
    plot_df$Model <- NULL                           # Drop Model column
    plot_df <- as.data.frame(t(plot_df))            # Transpose for radar chart
    
    # Add min/max rows for radar scale
    plot_df <- rbind(rep(100, ncol(plot_df)), rep(0, ncol(plot_df)), plot_df)
    
  } else {
    
    n_models <- length(mods)
    results <- list()
    
    tmp <- data.frame(Model = paste0("Model", seq_len(n_models)))
    
    for (m in seq_len(n_models)) {
      
      rep <- mods[[m]]$om$rep
      
      # Catch and SSB as rowSums
      catch_ts <- rowSums(rep$pred_catch)
      ssb_ts <- rowSums(rep$SSB)
      
      # Global Fbar
      n_fleets <- mods[[m]]$om$input$data$n_fleets[1]
      n_regions <- mods[[m]]$om$input$data$n_regions[1]
      fbar_ts <- rep$Fbar[, n_fleets + n_regions + 1]
      
      if(is.null(method)) method = "median"
      # First and last period means
      if (method == "median") {
        tmp$Catch_first[m] <- median(catch_ts[start.years:(start.years + use.n.years.first - 1)])
        tmp$SSB_first[m] <- median(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
        tmp$Fbar_first[m] <- median(fbar_ts[start.years:(start.years + use.n.years.first - 1)])
        tmp$Catch_last[m] <- median(tail(catch_ts, use.n.years.last))
        tmp$SSB_last[m] <- median(tail(ssb_ts, use.n.years.last))
        tmp$Fbar_last[m] <- median(tail(fbar_ts, use.n.years.last))
      } else if (method == "mean") {
        tmp$Catch_first[m] <- mean(catch_ts[start.years:(start.years + use.n.years.first - 1)])
        tmp$SSB_first[m] <- mean(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
        tmp$Fbar_first[m] <- mean(fbar_ts[start.years:(start.years + use.n.years.first - 1)])
        tmp$Catch_last[m] <- mean(tail(catch_ts, use.n.years.last))
        tmp$SSB_last[m] <- mean(tail(ssb_ts, use.n.years.last))
        tmp$Fbar_last[m] <- mean(tail(fbar_ts, use.n.years.last))
      }
    }
    
    # Normalize: higher is better except for Fbar
    for (v in c("Catch_first", "SSB_first", "Catch_last", "SSB_last")) {
      range_val <- max(tmp[[v]]) - min(tmp[[v]])
      if (range_val == 0) {
        tmp[[v]] <- 100  # All values are identical; assign full score
      } else {
        tmp[[v]] <- 100 * (tmp[[v]] - min(tmp[[v]])) / range_val
      }
    }
    
    for (v in c("Fbar_first", "Fbar_last")) {
      range_val <- max(tmp[[v]]) - min(tmp[[v]])
      if (range_val == 0) {
        tmp[[v]] <- 100  # All values identical; assign full score
      } else {
        norm_f <- (tmp[[v]] - min(tmp[[v]])) / range_val
        tmp[[v]] <- 100 * (1 - norm_f)
      }
    }
    
    results[[1]] <- tmp
    
    # Combine and compute median across realizations
    combined <- bind_rows(results, .id = "Realization")
    scores_median <- combined %>%
      group_by(Model) %>%
      summarise(across(-Realization, median), .groups = "drop")
    
    # Optional renaming
    if (!is.null(new_model_names)) {
      if (length(new_model_names) != length(unique(scores_median$Model))) {
        stop("Length of new_model_names must match the number of models.")
      }
      scores_median$Model <- factor(scores_median$Model,
                                    levels = paste0("Model", seq_along(new_model_names)),
                                    labels = new_model_names)
    }
    
    # Radar chart prep
    plot_df <- as.data.frame(scores_median)         # Ensure it's a base data.frame
    rownames(plot_df) <- plot_df$Model              # Set row names safely
    plot_df$Model <- NULL                           # Drop Model column
    plot_df <- as.data.frame(t(plot_df))            # Transpose for radar chart
    
    # Add min/max rows for radar scale
    plot_df <- rbind(rep(100, ncol(plot_df)), rep(0, ncol(plot_df)), plot_df)
    
  }
  
  if (ncol(plot_df) < 3) {
    message("Radar chart needs at least 3 models. Showing barplot instead.")
    return(invisible(NULL))
  }
  
  # Plot without any further changes
  colors <- viridisLite::viridis(n = nrow(plot_df) - 2, option = col.opt) # -2 for min/max rows
  
  # Create the new subfolder if it doesn't exist
  new_sub_dir <- file.path(main.dir, sub.dir, "Radar_Holistic_Plot")
  
  if (!file.exists(new_sub_dir)){
    dir.create(new_sub_dir)
  }
  
  # Correct way to assign expressions to rownames
  my_legend_labels <- c(expression(Catch[ST]),
                        expression(SSB[ST]),
                        expression(F[ST]),
                        expression(Catch[LT]),
                        expression(SSB[LT]),
                        expression(F[LT]))
  
  # --- Start of Modified Plotting Block ---
  
  # Save to PNG with a two-column layout
  output_file <- file.path(main.dir, sub.dir, "Radar_Holistic_Plot", "model_performance_radar.png")
  png(filename = output_file, width = width, height = height, units = "in", res = dpi)
  
  # 1. SET UP LAYOUT: Create a 2-column layout (left is 3x wider than right)
  layout(matrix(c(1, 2), nrow = 1), widths = c(3, 1))
  
  # 2. PLOT CHART: Plot radarchart in the left column (area 1)
  par(mar = c(1, 1, 1, 1)) # Set margins for the chart
  radarchart(plot_df,
             axistype = 4,
             pcol = colors,
             plwd = 3,
             plty = 1:length(colors),
             cglcol = "grey80",
             cglty = 1,
             axislabcol = "grey30",
             vlcex = 1.2)
  
  # 3. PLOT LEGEND: Plot legend in the right column (area 2)
  par(mar = c(1, 1, 1, 1))  # Set margins for the legend area
  plot.new()               # Create an empty plot area
  legend("center",         # Place legend in the center of its area
         legend = my_legend_labels,
         col = colors,
         lty = 1:length(colors),
         lwd = 3,
         cex = 0.9,
         y.intersp = 1.5)
  
  dev.off()
  
  # --- Also plot to screen (inline) using the same layout ---
  op <- par(no.readonly = TRUE) # Save original graphics settings
  on.exit(par(op))             # Restore settings when function exits
  
  layout(matrix(c(1, 2), nrow = 1), widths = c(3, 1))
  
  par(mar = c(1, 1, 2, 1))
  radarchart(plot_df,
             axistype = 4, pcol = colors, plwd = 3, plty = 1:length(colors),
             cglcol = "grey80", cglty = 1, axislabcol = "grey30", vlcex = 1.2)
  
  par(mar = c(1, 1, 2, 1))
  plot.new()
  legend("center",         # Place legend in the center of its area
         legend = my_legend_labels,
         col = colors,
         lty = 1:length(colors),
         lwd = 3,
         cex = 0.9,
         y.intersp = 2, 
         bty = "n")
  
  # --- End of Modified Plotting Block ---
}


plot_model_performance_radar2 <- function(mods, is.nsim, main.dir, sub.dir, 
                                          width = 10, height = 10, dpi = 300, col.opt = "D",
                                          method = NULL,
                                          use.n.years.first = 5,
                                          use.n.years.last = 5,
                                          start.years = 1, 
                                          new_model_names = NULL) {
  
  library(dplyr)
  library(tidyr)
  library(fmsb)
  library(viridis)
  library(viridisLite)
  
  # Helper: calculate Average Annual Catch Variation (AACV)
  calculate_aacv <- function(values) {
    if (!is.numeric(values)) {
      stop("Input must be a numeric vector.")
    }
    diffs <- abs(diff(values))
    aacv <- sum(diffs) / sum(values[-length(values)])
    return(aacv)
  }
  
  if (is.nsim) {
    
    n_models <- length(mods[[1]])
    n_reps <- length(mods)
    
    results <- list()
    
    for (r in seq_len(n_reps)) {
      tmp <- data.frame(Model = paste0("Model", seq_len(n_models)))
      
      for (m in seq_len(n_models)) {
        rep <- mods[[r]][[m]]$om$rep
        
        # Catch and SSB as rowSums
        catch_ts <- rowSums(rep$pred_catch)
        ssb_ts <- rowSums(rep$SSB)
        
        # Global Fbar
        n_fleets <- mods[[r]][[m]]$om$input$data$n_fleets[1]
        n_regions <- mods[[r]][[m]]$om$input$data$n_regions[1]
        
        # Catch variation
        catch_aacv <- calculate_aacv(rowSums(rep$pred_catch))
        
        # SSB variation
        ssb_aacv <- calculate_aacv(rowSums(rep$SSB))
        
        # Fbar variation (global is last col)
        fbar_aacv <- calculate_aacv(rep$Fbar[, ncol(rep$Fbar)])
        
        # First and last period means
        if(is.null(method)) method = "median"
        if (method == "median") {
          tmp$Catch_first[m] <- median(catch_ts[start.years:(start.years + use.n.years.first - 1)])
          tmp$SSB_first[m] <- median(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
          tmp$Catch_last[m] <- median(tail(catch_ts, use.n.years.last))
          tmp$SSB_last[m] <- median(tail(ssb_ts, use.n.years.last))
          tmp$catch_aacv[m] <- catch_aacv
          tmp$ssb_aacv[m] <- ssb_aacv
          tmp$fbar_aacv[m] <- fbar_aacv
        } else if (method == "mean") {
          tmp$Catch_first[m] <- mean(catch_ts[start.years:(start.years + use.n.years.first - 1)])
          tmp$SSB_first[m] <- mean(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
          tmp$Catch_last[m] <- mean(tail(catch_ts, use.n.years.last))
          tmp$SSB_last[m] <- mean(tail(ssb_ts, use.n.years.last))
          tmp$catch_aacv[m] <- catch_aacv
          tmp$ssb_aacv[m] <- ssb_aacv
          tmp$fbar_aacv[m] <- fbar_aacv
        } 
      }
      
      # Normalize: higher is better except for Fbar
      for (v in c("Catch_first", "SSB_first", "Catch_last", "SSB_last")) {
        range_val <- max(tmp[[v]]) - min(tmp[[v]])
        if (range_val == 0) {
          tmp[[v]] <- 100
        } else {
          tmp[[v]] <- 100 * (tmp[[v]] - min(tmp[[v]])) / range_val
        }
      }
      for (v in c("catch_aacv", "ssb_aacv", "fbar_aacv")) {
        range_val <- max(tmp[[v]]) - min(tmp[[v]])
        if (range_val == 0) {
          tmp[[v]] <- 100
        } else {
          norm_f <- (tmp[[v]] - min(tmp[[v]])) / range_val
          tmp[[v]] <- 100 * (1 - norm_f)
        }
      }
      
      results[[r]] <- tmp
    }
    
    # Combine and compute median across realizations
    combined <- bind_rows(results, .id = "Realization")
    scores_median <- combined %>%
      group_by(Model) %>%
      summarise(across(-Realization, median), .groups = "drop")
    
    # Optional renaming
    if (!is.null(new_model_names)) {
      if (length(new_model_names) != length(unique(scores_median$Model))) {
        stop("Length of new_model_names must match the number of models.")
      }
      scores_median$Model <- factor(scores_median$Model,
                                    levels = paste0("Model", seq_along(new_model_names)),
                                    labels = new_model_names)
    }
    
    # Radar chart prep
    plot_df <- as.data.frame(scores_median)         # Ensure it's a base data.frame
    rownames(plot_df) <- plot_df$Model              # Set row names safely
    plot_df$Model <- NULL                           # Drop Model column
    plot_df <- as.data.frame(t(plot_df))            # Transpose for radar chart
    
    # Add min/max rows for radar scale
    plot_df <- rbind(rep(100, ncol(plot_df)), rep(0, ncol(plot_df)), plot_df)
    
  } else {
    
    n_models <- length(mods)
    results <- list()
    
    tmp <- data.frame(Model = paste0("Model", seq_len(n_models)))
    
    for (m in seq_len(n_models)) {
      
      rep <- mods[[m]]$om$rep
      
      # Catch and SSB as rowSums
      catch_ts <- rowSums(rep$pred_catch)
      ssb_ts <- rowSums(rep$SSB)
      
      # Global Fbar
      n_fleets <- mods[[m]]$om$input$data$n_fleets[1]
      n_regions <- mods[[m]]$om$input$data$n_regions[1]
      
      # Catch variation
      catch_aacv <- calculate_aacv(rowSums(rep$pred_catch))
      
      # SSB variation
      ssb_aacv <- calculate_aacv(rowSums(rep$SSB))
      
      # Fbar variation (global is last col)
      fbar_aacv <- calculate_aacv(rep$Fbar[, ncol(rep$Fbar)])
      
      if(is.null(method)) method = "median"
      # First and last period means
      if (method == "median") {
        tmp$Catch_first[m] <- median(catch_ts[start.years:(start.years + use.n.years.first - 1)])
        tmp$SSB_first[m] <- median(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
        tmp$Catch_last[m] <- median(tail(catch_ts, use.n.years.last))
        tmp$SSB_last[m] <- median(tail(ssb_ts, use.n.years.last))
        tmp$catch_aacv[m] <- catch_aacv
        tmp$ssb_aacv[m] <- ssb_aacv
        tmp$fbar_aacv[m] <- fbar_aacv
      } else if (method == "mean") {
        tmp$Catch_first[m] <- mean(catch_ts[start.years:(start.years + use.n.years.first - 1)])
        tmp$SSB_first[m] <- mean(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
        tmp$Catch_last[m] <- mean(tail(catch_ts, use.n.years.last))
        tmp$SSB_last[m] <- mean(tail(ssb_ts, use.n.years.last))
        tmp$catch_aacv[m] <- catch_aacv
        tmp$ssb_aacv[m] <- ssb_aacv
        tmp$fbar_aacv[m] <- fbar_aacv
      }
    }
    
    # Normalize: higher is better except for Fbar
    for (v in c("Catch_first", "SSB_first", "Catch_last", "SSB_last")) {
      range_val <- max(tmp[[v]]) - min(tmp[[v]])
      if (range_val == 0) {
        tmp[[v]] <- 100  # All values are identical; assign full score
      } else {
        tmp[[v]] <- 100 * (tmp[[v]] - min(tmp[[v]])) / range_val
      }
    }
    
    for (v in c("catch_aacv", "ssb_aacv", "fbar_aacv")) {
      range_val <- max(tmp[[v]]) - min(tmp[[v]])
      if (range_val == 0) {
        tmp[[v]] <- 100  # All values identical; assign full score
      } else {
        norm_f <- (tmp[[v]] - min(tmp[[v]])) / range_val
        tmp[[v]] <- 100 * (1 - norm_f)
      }
    }
    
    results[[1]] <- tmp
    
    # Combine and compute median across realizations
    combined <- bind_rows(results, .id = "Realization")
    scores_median <- combined %>%
      group_by(Model) %>%
      summarise(across(-Realization, median), .groups = "drop")
    
    # Optional renaming
    if (!is.null(new_model_names)) {
      if (length(new_model_names) != length(unique(scores_median$Model))) {
        stop("Length of new_model_names must match the number of models.")
      }
      scores_median$Model <- factor(scores_median$Model,
                                    levels = paste0("Model", seq_along(new_model_names)),
                                    labels = new_model_names)
    }
    
    # Radar chart prep
    plot_df <- as.data.frame(scores_median)         # Ensure it's a base data.frame
    rownames(plot_df) <- plot_df$Model              # Set row names safely
    plot_df$Model <- NULL                           # Drop Model column
    plot_df <- as.data.frame(t(plot_df))            # Transpose for radar chart
    
    # Add min/max rows for radar scale
    plot_df <- rbind(rep(100, ncol(plot_df)), rep(0, ncol(plot_df)), plot_df)
    
  }
  
  if (ncol(plot_df) < 3) {
    message("Radar chart needs at least 3 models. Showing barplot instead.")
    return(invisible(NULL))
  }
  
  # Plot without any further changes
  colors <- viridisLite::viridis(n = nrow(plot_df) - 2, option = col.opt) # -2 for min/max rows
  
  # Create the new subfolder if it doesn't exist
  new_sub_dir <- file.path(main.dir, sub.dir, "Radar_Holistic_Plot")
  
  if (!file.exists(new_sub_dir)){
    dir.create(new_sub_dir)
  }
  
  rownames(plot_df)
  # Correct way to assign expressions to rownames
  my_legend_labels <- c(expression(Catch[ST]),
                        expression(SSB[ST]),
                        expression(Catch[LT]),
                        expression(SSB[LT]),
                        expression(Catch[AAV]),
                        expression(SSB[AAV]),
                        expression(F[AAV]))
  
  # Save to PNG
  output_file <- file.path(file.path(main.dir, sub.dir, "Radar_Holistic_Plot", "model_performance_radar2.png"))
  png(filename = output_file, width = width, height = height, units = "in", res = dpi)
  
  # 1. SET UP LAYOUT: Create a 2-column layout (left is 3x wider than right)
  layout(matrix(c(1, 2), nrow = 1), widths = c(3, 1))
  
  # 2. PLOT CHART: Plot radarchart in the left column (area 1)
  par(mar = c(1, 1, 1, 1)) # Set margins for the chart
  radarchart(plot_df,
             axistype = 4,
             pcol = colors,
             plwd = 3,
             plty = 1:length(colors),
             cglcol = "grey80",
             cglty = 1,
             axislabcol = "grey30",
             vlcex = 1.2)
  
  # 3. PLOT LEGEND: Plot legend in the right column (area 2)
  par(mar = c(1, 1, 1, 1))  # Set margins for the legend area
  plot.new()               # Create an empty plot area
  legend("center",         # Place legend in the center of its area
         legend = my_legend_labels,
         col = colors,
         lty = 1:length(colors),
         lwd = 3,
         cex = 0.9,
         y.intersp = 1.5)
  
  dev.off()
  
  # --- Also plot to screen (inline) using the same layout ---
  op <- par(no.readonly = TRUE) # Save original graphics settings
  on.exit(par(op))             # Restore settings when function exits
  
  layout(matrix(c(1, 2), nrow = 1), widths = c(3, 1))
  
  par(mar = c(1, 1, 2, 1))
  radarchart(plot_df,
             axistype = 4, pcol = colors, plwd = 3, plty = 1:length(colors),
             cglcol = "grey80", cglty = 1, axislabcol = "grey30", vlcex = 1.2)
  
  par(mar = c(1, 1, 2, 1))
  plot.new()
  legend("center",         # Place legend in the center of its area
         legend = my_legend_labels,
         col = colors,
         lty = 1:length(colors),
         lwd = 3,
         cex = 0.9,
         y.intersp = 2, 
         bty = "n")
  
  # --- End of Modified Plotting Block ---
}

plot_model_performance_radar3 <- function(mods, is.nsim, main.dir, sub.dir, 
                                          width = 10, height = 10, dpi = 300, col.opt = "D",
                                          method = NULL,
                                          use.n.years.first = 5,
                                          use.n.years.last = 5,
                                          start.years = 1, 
                                          new_model_names = NULL) {
  
  library(dplyr)
  library(tidyr)
  library(fmsb)
  library(viridis)
  library(viridisLite)
  
  # Helper: calculate Average Annual Catch Variation (AACV)
  calculate_aacv <- function(values) {
    if (!is.numeric(values)) {
      stop("Input must be a numeric vector.")
    }
    diffs <- abs(diff(values))
    aacv <- sum(diffs) / sum(values[-length(values)])
    return(aacv)
  }
  
  if (is.nsim) {
    
    n_models <- length(mods[[1]])
    n_reps <- length(mods)
    
    results <- list()
    
    for (r in seq_len(n_reps)) {
      tmp <- data.frame(Model = paste0("Model", seq_len(n_models)))
      
      for (m in seq_len(n_models)) {
        rep <- mods[[r]][[m]]$om$rep
        
        # Catch and SSB as rowSums
        catch_ts <- rowSums(rep$pred_catch)
        ssb_ts <- rowSums(rep$SSB)
        
        # Global Fbar
        n_fleets <- mods[[r]][[m]]$om$input$data$n_fleets[1]
        n_regions <- mods[[r]][[m]]$om$input$data$n_regions[1]
        
        # Catch variation
        catch_aacv <- calculate_aacv(rowSums(rep$pred_catch))
        
        # SSB variation
        ssb_aacv <- calculate_aacv(rowSums(rep$SSB))
        
        # Fbar variation (global is last col)
        fbar_aacv <- calculate_aacv(rep$Fbar[, ncol(rep$Fbar)])
        
        # First and last period means
        if(is.null(method)) method = "median"
        if (method == "median") {
          tmp$Catch_first[m] <- median(catch_ts[start.years:(start.years + use.n.years.first - 1)])
          tmp$SSB_first[m] <- median(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
          tmp$Catch_last[m] <- median(tail(catch_ts, use.n.years.last))
          tmp$SSB_last[m] <- median(tail(ssb_ts, use.n.years.last))
          tmp$catch_aacv[m] <- catch_aacv
          tmp$ssb_aacv[m] <- ssb_aacv
          tmp$fbar_aacv[m] <- fbar_aacv
        } else if (method == "mean") {
          tmp$Catch_first[m] <- mean(catch_ts[start.years:(start.years + use.n.years.first - 1)])
          tmp$SSB_first[m] <- mean(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
          tmp$Catch_last[m] <- mean(tail(catch_ts, use.n.years.last))
          tmp$SSB_last[m] <- mean(tail(ssb_ts, use.n.years.last))
          tmp$catch_aacv[m] <- catch_aacv
          tmp$ssb_aacv[m] <- ssb_aacv
          tmp$fbar_aacv[m] <- fbar_aacv
        }
        
        if (is.null(mods[[1]][[1]]$om$rep$log_SSB_FXSPR)) {
          message("Biological Reference Point has not been calculated internally!")
          return(invisible(NULL))
        }
        
        tmp1 <- rep$Fbar[, ncol(rep$Fbar)]
        tmp2 <- exp(rep$log_Fbar_XSPR[, ncol(rep$Fbar)])
        temp <- tmp1/tmp2
        temp1 <- temp[start.years:(start.years + use.n.years.first - 1)]
        temp2 <- tail(temp, use.n.years.last)
        tmp$prob_first[m] <- mean(temp1 > 1, na.rm = TRUE)
        tmp$prob_last[m] <- mean(temp2 > 1, na.rm = TRUE)
      }
      
      # Normalize: higher is better except for Fbar
      for (v in c("Catch_first", "SSB_first", "Catch_last", "SSB_last")) {
        range_val <- max(tmp[[v]]) - min(tmp[[v]])
        if (range_val == 0) {
          tmp[[v]] <- 100
        } else {
          tmp[[v]] <- 100 * (tmp[[v]] - min(tmp[[v]])) / range_val
        }
      }
      for (v in c("catch_aacv", "ssb_aacv", "fbar_aacv","prob_first","prob_last")) {
        range_val <- max(tmp[[v]]) - min(tmp[[v]])
        if (range_val == 0) {
          tmp[[v]] <- 100
        } else {
          norm_f <- (tmp[[v]] - min(tmp[[v]])) / range_val
          tmp[[v]] <- 100 * (1 - norm_f)
        }
      }
      
      results[[r]] <- tmp
    }
    
    # Combine and compute median across realizations
    combined <- bind_rows(results, .id = "Realization")
    scores_median <- combined %>%
      group_by(Model) %>%
      summarise(across(-Realization, median), .groups = "drop")
    
    # Optional renaming
    if (!is.null(new_model_names)) {
      if (length(new_model_names) != length(unique(scores_median$Model))) {
        stop("Length of new_model_names must match the number of models.")
      }
      scores_median$Model <- factor(scores_median$Model,
                                    levels = paste0("Model", seq_along(new_model_names)),
                                    labels = new_model_names)
    }
    
    # Radar chart prep
    plot_df <- as.data.frame(scores_median)         # Ensure it's a base data.frame
    rownames(plot_df) <- plot_df$Model              # Set row names safely
    plot_df$Model <- NULL                           # Drop Model column
    plot_df <- as.data.frame(t(plot_df))            # Transpose for radar chart
    
    # Add min/max rows for radar scale
    plot_df <- rbind(rep(100, ncol(plot_df)), rep(0, ncol(plot_df)), plot_df)
    
  } else {
    
    n_models <- length(mods)
    results <- list()
    
    tmp <- data.frame(Model = paste0("Model", seq_len(n_models)))
    
    for (m in seq_len(n_models)) {
      
      rep <- mods[[m]]$om$rep
      
      # Catch and SSB as rowSums
      catch_ts <- rowSums(rep$pred_catch)
      ssb_ts <- rowSums(rep$SSB)
      
      # Global Fbar
      n_fleets <- mods[[m]]$om$input$data$n_fleets[1]
      n_regions <- mods[[m]]$om$input$data$n_regions[1]
      
      # Catch variation
      catch_aacv <- calculate_aacv(rowSums(rep$pred_catch))
      
      # SSB variation
      ssb_aacv <- calculate_aacv(rowSums(rep$SSB))
      
      # Fbar variation (global is last col)
      fbar_aacv <- calculate_aacv(rep$Fbar[, ncol(rep$Fbar)])
      
      if(is.null(method)) method = "median"
      # First and last period means
      if (method == "median") {
        tmp$Catch_first[m] <- median(catch_ts[start.years:(start.years + use.n.years.first - 1)])
        tmp$SSB_first[m] <- median(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
        tmp$Catch_last[m] <- median(tail(catch_ts, use.n.years.last))
        tmp$SSB_last[m] <- median(tail(ssb_ts, use.n.years.last))
        tmp$catch_aacv[m] <- catch_aacv
        tmp$ssb_aacv[m] <- ssb_aacv
        tmp$fbar_aacv[m] <- fbar_aacv
        
      } else if (method == "mean") {
        tmp$Catch_first[m] <- mean(catch_ts[start.years:(start.years + use.n.years.first - 1)])
        tmp$SSB_first[m] <- mean(ssb_ts[start.years:(start.years + use.n.years.first - 1)])
        tmp$Catch_last[m] <- mean(tail(catch_ts, use.n.years.last))
        tmp$SSB_last[m] <- mean(tail(ssb_ts, use.n.years.last))
        tmp$catch_aacv[m] <- catch_aacv
        tmp$ssb_aacv[m] <- ssb_aacv
        tmp$fbar_aacv[m] <- fbar_aacv
      }
    }
    
    if (is.null(rep$log_SSB_FXSPR)) {
      message("Biological Reference Point has not been calculated internally!")
      return(invisible(NULL))
    }
    
    tmp1 <- rep$Fbar[, ncol(rep$Fbar)]
    tmp2 <- exp(rep$log_Fbar_XSPR[, ncol(rep$Fbar)])
    temp <- tmp1/tmp2
    temp1 <- temp[start.years:(start.years + use.n.years.first - 1)]
    temp2 <- tail(temp, use.n.years.last)
    tmp$prob_first[m] <- mean(temp1 > 1, na.rm = TRUE)
    tmp$prob_last[m] <- mean(temp2 > 1, na.rm = TRUE)
    
    # Normalize: higher is better except for Fbar
    for (v in c("Catch_first", "SSB_first", "Catch_last", "SSB_last")) {
      range_val <- max(tmp[[v]]) - min(tmp[[v]])
      if (range_val == 0) {
        tmp[[v]] <- 100  # All values are identical; assign full score
      } else {
        tmp[[v]] <- 100 * (tmp[[v]] - min(tmp[[v]])) / range_val
      }
    }
    
    for (v in c("catch_aacv", "ssb_aacv", "fbar_aacv","prob_first","prob_last")) {
      range_val <- max(tmp[[v]]) - min(tmp[[v]])
      if (range_val == 0) {
        tmp[[v]] <- 100  # All values identical; assign full score
      } else {
        norm_f <- (tmp[[v]] - min(tmp[[v]])) / range_val
        tmp[[v]] <- 100 * (1 - norm_f)
      }
    }
    
    results[[1]] <- tmp
    
    # Combine and compute median across realizations
    combined <- bind_rows(results, .id = "Realization")
    scores_median <- combined %>%
      group_by(Model) %>%
      summarise(across(-Realization, median), .groups = "drop")
    
    # Optional renaming
    if (!is.null(new_model_names)) {
      if (length(new_model_names) != length(unique(scores_median$Model))) {
        stop("Length of new_model_names must match the number of models.")
      }
      scores_median$Model <- factor(scores_median$Model,
                                    levels = paste0("Model", seq_along(new_model_names)),
                                    labels = new_model_names)
    }
    
    # Radar chart prep
    plot_df <- as.data.frame(scores_median)         # Ensure it's a base data.frame
    rownames(plot_df) <- plot_df$Model              # Set row names safely
    plot_df$Model <- NULL                           # Drop Model column
    plot_df <- as.data.frame(t(plot_df))            # Transpose for radar chart
    
    # Add min/max rows for radar scale
    plot_df <- rbind(rep(100, ncol(plot_df)), rep(0, ncol(plot_df)), plot_df)
    
  }
  
  if (ncol(plot_df) < 3) {
    message("Radar chart needs at least 3 models. Showing barplot instead.")
    return(invisible(NULL))
  }
  
  # Plot without any further changes
  colors <- viridisLite::viridis(n = nrow(plot_df) - 2, option = col.opt) # -2 for min/max rows
  
  # Create the new subfolder if it doesn't exist
  new_sub_dir <- file.path(main.dir, sub.dir, "Radar_Holistic_Plot")
  
  if (!file.exists(new_sub_dir)){
    dir.create(new_sub_dir)
  }
  
  rownames(plot_df)
  # Correct way to assign expressions to rownames
  my_legend_labels <- c(expression(Catch[ST]),
                        expression(SSB[ST]),
                        expression(Catch[LT]),
                        expression(SSB[LT]),
                        expression(Catch[AAV]),
                        expression(SSB[AAV]),
                        expression(F[AAV]),
                        expression(P(F[ST] > F[MSY])),
                        expression(P(F[LT] > F[MSY])))
  
  # Save to PNG
  output_file <- file.path(file.path(main.dir, sub.dir, "Radar_Holistic_Plot", "model_performance_radar3.png"))
  png(filename = output_file, width = width, height = height, units = "in", res = dpi)
  
  # 1. SET UP LAYOUT: Create a 2-column layout (left is 3x wider than right)
  layout(matrix(c(1, 2), nrow = 1), widths = c(3, 1))
  
  # 2. PLOT CHART: Plot radarchart in the left column (area 1)
  par(mar = c(1, 1, 1, 1)) # Set margins for the chart
  radarchart(plot_df,
             axistype = 4,
             pcol = colors,
             plwd = 3,
             plty = 1:length(colors),
             cglcol = "grey80",
             cglty = 1,
             axislabcol = "grey30",
             vlcex = 1.2)
  
  # 3. PLOT LEGEND: Plot legend in the right column (area 2)
  par(mar = c(1, 1, 1, 1))  # Set margins for the legend area
  plot.new()               # Create an empty plot area
  legend("center",         # Place legend in the center of its area
         legend = my_legend_labels,
         col = colors,
         lty = 1:length(colors),
         lwd = 3,
         cex = 0.9,
         y.intersp = 1.5)
  
  dev.off()
  
  # --- Also plot to screen (inline) using the same layout ---
  op <- par(no.readonly = TRUE) # Save original graphics settings
  on.exit(par(op))             # Restore settings when function exits
  
  layout(matrix(c(1, 2), nrow = 1), widths = c(3, 1))
  
  par(mar = c(1, 1, 2, 1))
  radarchart(plot_df,
             axistype = 4, pcol = colors, plwd = 3, plty = 1:length(colors),
             cglcol = "grey80", cglty = 1, axislabcol = "grey30", vlcex = 1.2)
  
  par(mar = c(1, 1, 2, 1))
  plot.new()
  legend("center",         # Place legend in the center of its area
         legend = my_legend_labels,
         col = colors,
         lty = 1:length(colors),
         lwd = 3,
         cex = 0.9,
         y.intersp = 2, 
         bty = "n")
  
  # --- End of Modified Plotting Block ---
}
