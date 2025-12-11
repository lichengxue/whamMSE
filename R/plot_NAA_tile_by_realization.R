#' Tile plots of NAA deviations by realization (OM or EM)
#'
#' Create heatmap (tile) plots of numbers-at-age deviations (`NAA_devs`) across
#' years and ages, for each realization in an MSE run. The function can plot
#' either the operating model (OM) deviations or the estimation model (EM)
#' deviations, optionally looping over all assessments and models.
#'
#' For `type = "OM"`, the function:
#' \itemize{
#'   \item assumes \code{mods} is a list of realizations, each containing at least
#'         one model with an \code{$om} component,
#'   \item extracts \code{NAA_devs} from the OM report object
#'         (\code{mod$rep$NAA_devs}) for each stock and region, and
#'   \item plots deviations as tiles over \code{Year} (x-axis) and \code{Age} (y-axis),
#'         faceted by \code{Stock} and \code{Region} with one row per realization.
#' }
#'
#' For `type = "EM"`, the function:
#' \itemize{
#'   \item assumes \code{mods} is a list of realizations, each being a list of models
#'         with elements \code{$em_list} (assessment fits) and \code{$em_input},
#'   \item loops over assessments within each realization and across all models,
#'   \item extracts \code{NAA_devs} for each stock and region from the EM report object,
#'         and
#'   \item plots deviations as tiles over \code{Year} and \code{Age}, faceted by
#'         a combined label of model, stock, and region.
#' }
#'
#' All plots are written to a single multi-page PDF in
#' \code{file.path(main.dir, sub.dir, pdf_name)}.
#'
#' @param mods A list of MSE results. For OM plots, this must be a list of
#'   realizations, each containing at least one model with an \code{$om}
#'   component (with \code{$env$data}, \code{$rep}, \code{$years}, etc.).
#'   For EM plots, each realization must be a list of models, each having
#'   \code{$em_list} (a list of assessment fits with \code{NAA_devs}) and
#'   \code{$em_input} (input lists with \code{years_full}).
#' @param type Character string indicating which component to plot:
#'   \code{"OM"} for operating model \code{NAA_devs}, or \code{"EM"} for
#'   estimation model \code{NAA_devs}.
#' @param em_index Optional integer index of the estimation model to focus on.
#'   Currently not used; reserved for potential future extensions.
#' @param main.dir Character path to the main output directory.
#' @param sub.dir Character folder name inside \code{main.dir} where the PDF
#'   will be written. The directory is created if it does not exist.
#' @param pdf_name Name of the output PDF file (default:
#'   \code{"NAA_tile_by_realization.pdf"}).
#' @param width,height Numeric width and height (in inches) for the PDF device.
#'   If both are \code{NULL}, defaults depend on \code{type}:
#'   \itemize{
#'     \item \code{type == "OM"}: \code{width = 15}, \code{height = 5}
#'     \item \code{type == "EM"}: \code{width = 10}, \code{height = 10}
#'   }
#' @param fontfam Optional font family name passed to the PDF device. Currently
#'   not explicitly used in the plotting theme but can be set via global theme
#'   modifications if desired.
#' @param realizations_to_plot Optional integer vector of realization indices
#'   to include. If \code{NULL}, all realizations (\code{seq_len(length(mods))})
#'   are plotted.
#' @param new_model_names Optional character vector of new model names for EM
#'   plots. If provided, must have length equal to the number of models in each
#'   realization, and will be used to relabel the \code{Model} factor.
#'
#' @details
#' The function checks for a "multi-simulation" structure in \code{mods} by
#' testing if \code{mods[[1]]} is a list and if \code{mods[[1]][[1]]$om} is
#' non-\code{NULL}. If this is not the case, the function stops with an error,
#' as it expects \code{mods} to be a list of realizations.
#'
#' For OM plots (`type = "OM"`), the code uses:
#' \itemize{
#'   \item \code{mod$env$data} to obtain dimensions (\code{n_ages}, \code{n_stocks},
#'         \code{n_regions}, \code{n_years_proj}),
#'   \item \code{mod$rep$NAA_devs[s, rn, , ]} as the age-by-year deviations for
#'         stock \code{s} and region \code{rn}, and
#'   \item \code{mod$years_full} when projection years exist, otherwise
#'         \code{mod$years}.
#' }
#'
#' For EM plots (`type = "EM"`), the function loops over assessments within
#' each model and realization:
#' \itemize{
#'   \item \code{em_list[[x]]$NAA_devs[s, rn, , ]} is used as the age-by-year
#'         deviations for assessment \code{x},
#'   \item \code{em_input_list[[x]]$years_full} provides the full year vector,
#'   \item a combined facet label \code{"ModelX_stock_s_region_r"} is built for
#'         each model/stock/region combination.
#' }
#'
#' The color scale for deviations uses \code{viridis::scale_fill_viridis()}.
#'
#' @return
#' Invisibly returns \code{NULL}. The primary side effect is the creation of a
#' multi-page PDF of NAA deviation tiles in \code{file.path(main.dir, sub.dir, pdf_name)}.
#'
#' @seealso
#' Other diagnostic plots in the package (e.g. terminal SSB/F bias plots) for
#' assessing estimation performance across realizations and assessments.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Suppose `mods` is a list of realizations produced by whamMSE:
#' # Plot OM NAA devs for all realizations:
#' plot_NAA_tile_by_realization(
#'   mods      = mods,
#'   type      = "OM",
#'   main.dir  = "Results",
#'   sub.dir   = "Diagnostics"
#' )
#'
#' # Plot EM NAA devs for realizations 1 and 2 only, relabeling models:
#' plot_NAA_tile_by_realization(
#'   mods                 = mods,
#'   type                 = "EM",
#'   main.dir             = "Results",
#'   sub.dir              = "Diagnostics",
#'   realizations_to_plot = c(1, 2),
#'   new_model_names      = c("Equal", "FleetCatch", "IndexEqual")
#' )
#' }
plot_NAA_tile_by_realization <- function(mods, type = "OM", em_index = NULL,
                                         main.dir, sub.dir,
                                         pdf_name = "NAA_tile_by_realization.pdf",
                                         width = NULL, height = NULL,
                                         fontfam = "", realizations_to_plot = NULL,
                                         new_model_names = NULL) {
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(grDevices)
  
  is.nsim <- is.list(mods[[1]]) && !is.null(mods[[1]][[1]]$om)
  
  if (!is.nsim) stop("mods must be a list of realizations, each with multiple models.")
  
  n_realizations <- length(mods)
  n_models <- length(mods[[1]])
  
  if (is.null(realizations_to_plot)) {
    realizations_to_plot <- seq_len(n_realizations)
  }
  
  dir.create(file.path(main.dir, sub.dir), recursive = TRUE, showWarnings = FALSE)
  
  if (is.null(width) && is.null(height)) {
    if (type == "OM") {
      width  <- 15
      height <- 5
    } else {
      width  <- 10
      height <- 10
    }
  }
  pdf(file = file.path(main.dir, sub.dir, pdf_name), width = width, height = height, family = fontfam)
  
  for (r in realizations_to_plot) {
    df.plot.all <- data.frame()
    
    if (type == "OM") {
      mod  <- mods[[r]][[1]]$om
      dat  <- mod$env$data
      rep  <- mod$rep
      years <- mod$years
      n_ages <- dat$n_ages
      ages.lab <- if (!is.null(mod$ages.lab)) mod$ages.lab else 1:n_ages
      years_full <- if (dat$n_years_proj > 0) mod$years_full else years
      
      for (s in 1:dat$n_stocks) for (rn in 1:dat$n_regions) {
        df.tmp <- as.data.frame(rep$NAA_devs[s, rn, , ])
        colnames(df.tmp) <- paste0("Age_", 1:n_ages)
        df.tmp <- cbind.data.frame(
          Realization = paste0("R", r),
          Stock  = mod$input$stock_names[s],
          Region = mod$input$region_names[rn],
          Year   = years_full,
          df.tmp
        )
        df.tmp <- df.tmp %>%
          pivot_longer(starts_with("Age"),
                       names_to = "Age",
                       names_prefix = "Age_",
                       values_to = "Deviation") %>%
          mutate(Age = factor(as.integer(Age), labels = ages.lab))
        df.plot.all <- rbind(df.plot.all, df.tmp)
      }
      
      p <- ggplot(df.plot.all, aes(x = Year, y = Age, fill = Deviation)) +
        geom_tile() +
        scale_x_continuous(expand = c(0, 0)) +
        scale_y_discrete(expand = c(0, 0)) +
        facet_grid(Realization ~ Stock + Region, labeller = label_both) +
        theme_bw() +
        ggtitle(paste0("OM NAA_devs for realization ", r)) +
        theme(plot.title = element_text(hjust = 0.5)) +
        viridis::scale_fill_viridis()
      
      print(p)
    }
    
    if (type == "EM") {
      cat("Processing EM plots...\n")
      
      em_list_first_model <- mods[[r]][[1]]$em_list
      n_ems <- length(em_list_first_model)
      
      for (x in seq_len(n_ems)) {
        cat("  Realization", r, "- Assessment", x, "\n")
        df.plot.all <- data.frame()
        
        for (m in 1:n_models) {
          em_list       <- mods[[r]][[m]]$em_list
          em_input_list <- mods[[r]][[m]]$em_input
          
          rep        <- em_list[[x]]
          years_full <- em_input_list[[x]]$years_full
          n_ages     <- dim(rep$NAA)[4]
          ages.lab   <- 1:n_ages
          
          for (s in 1:dim(rep$NAA)[1]) for (rn in 1:dim(rep$NAA)[2]) {
            df.tmp <- as.data.frame(rep$NAA_devs[s, rn, , ])
            colnames(df.tmp) <- paste0("Age_", 1:n_ages)
            df.tmp <- cbind.data.frame(
              Realization = paste0("R", r),
              Model       = paste0("Model", m),
              Assessment  = paste0("Assessment", x),
              Stock       = paste0("stock_", s),
              Region      = paste0("region_", rn),
              Year        = years_full,
              df.tmp
            )
            df.tmp <- df.tmp %>%
              pivot_longer(starts_with("Age"),
                           names_to = "Age",
                           names_prefix = "Age_",
                           values_to = "Deviation") %>%
              mutate(Age = factor(as.integer(Age), labels = ages.lab)) %>%
              mutate(FacetLabel = paste0("Model", m, "_", Stock, "_", Region))
            df.plot.all <- rbind(df.plot.all, df.tmp)
          }
        }
        
        if (!is.null(new_model_names)) {
          df.plot.all$Model <- factor(df.plot.all$Model,
                                      levels = paste0("Model", seq_len(n_models)),
                                      labels = new_model_names)
        }
        
        p <- ggplot(df.plot.all, aes(x = Year, y = Age, fill = Deviation)) +
          geom_tile() +
          scale_x_continuous(expand = c(0, 0)) +
          scale_y_discrete(expand = c(0, 0)) +
          facet_wrap(~ FacetLabel, ncol = 2) +
          theme_bw() +
          ggtitle(paste0("Realization R", r, ", Assessment ", x, "\nModels ",
                         paste(unique(df.plot.all$Model), collapse = ", "))) +
          theme(plot.title = element_text(hjust = 0.5)) +
          viridis::scale_fill_viridis()
        
        print(p)
      }
    }
  }
  
  dev.off()
}
