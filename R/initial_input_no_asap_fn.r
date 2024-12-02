#' Prepare wham input without ASAP file
#'
#' Function to prepare wham input without ASAP file
#' 
#' @param input Initial wham input
#' @param basic_info A list of basic information 
#' 
#' @return wham input file
#' 
#' @export
initial_input_no_asap_fn <- function(input, basic_info) {
  
  input$years <- 1975:2014
  data <- input$data
  
  if (!is.null(basic_info$years)) {
    if (!is.integer(basic_info$years)) stop("basic_info$years has been specified, but it is not an integer vector")
    input$years <- basic_info$years
  }
  data$n_years_model <- length(input$years)
  
  data$n_stocks <- 1
  if (!is.null(basic_info$n_stocks)) {
    if (!is.integer(basic_info$n_stocks)) stop("basic_info$n_stocks has been specified, but it is not an integer")
    data$n_stocks <- basic_info$n_stocks
  }
  
  data$n_regions <- 1
  if (!is.null(basic_info$n_regions)) {
    if (!is.integer(basic_info$n_regions)) stop("basic_info$n_regions has been specified, but it is not an integer.")
    data$n_regions <- basic_info$n_regions
  }
  
  data$n_seasons <- 1
  if (!is.null(basic_info$n_seasons)) {
    if (!is.integer(basic_info$n_seasons)) stop("basic_info$n_seasons has been specified, but it is not an integer.")
    data$n_seasons <- basic_info$n_seasons
  }
  
  data$fracyr_seasons <- 1
  if (!is.null(basic_info$fracyr_seasons)) {
    if (!is.numeric(basic_info$fracyr_seasons)) stop("basic_info$fracyr_seasons has been specified, but it is not numeric.")
    if (length(basic_info$fracyr_seasons) != data$n_seasons) stop("length of basic_info$fracyr_seasons not equal to n_seasons.")
    data$fracyr_seasons <- basic_info$fracyr_seasons
    if (length(basic_info$fracyr_seasons) != 1) {
      input$log$misc <- c(input$log$misc, "basic_info$fracyr_seasons implies more than one season. Ensure that
          input$data$fracyr_SSB, input$data$spawn_seasons, input$data$fracyr_indices, input$data$index_seasons, are specified appropriately.\n")
    }
  }
  
  # Handling movement-related parameters
  if (!is.null(basic_info$onto_move)) {
    data$onto_move <- if (basic_info$onto_move == 1) {
      matrix(1, data$n_stocks, data$n_regions)
    } else {
      matrix(0, data$n_stocks, data$n_regions)
    }
    
    # Handling onto_move_pars
    data$onto_move_pars <- if (basic_info$onto_move == 1) {
      if (is.null(basic_info$onto_move_pars)) {
        # Default: set to ones if not provided
        array(1, dim = c(data$n_stocks, data$n_regions, 4))
      } else if (is.vector(basic_info$onto_move_pars)) {
        # If it's a vector, expand it to an array with the correct dimensions
        if (length(basic_info$onto_move_pars) == 4) {
          array(rep(basic_info$onto_move_pars, each = data$n_stocks * data$n_regions), dim = c(data$n_stocks, data$n_regions, 4))
        } else {
          stop("onto_move_pars must be a vector of length 4 or an array with dimensions (n_stocks, n_regions, 4).")
        }
      } else if (is.array(basic_info$onto_move_pars)) {
        # Check if the provided array matches the expected dimensions
        if (!all(dim(basic_info$onto_move_pars) == c(data$n_stocks, data$n_regions, 4))) {
          stop("onto_move_pars array must have dimensions (n_stocks, n_regions, 4).")
        }
        basic_info$onto_move_pars
      } else {
        stop("onto_move_pars must be either a vector of length 4 or an array with dimensions (n_stocks, n_regions, 4).")
      }
    } else {
      NULL  # No movement parameters needed if onto_move != 1
    }
  }
  
  # Handling meta-population indicator
  data$move_dyn <- if (!is.null(basic_info$move_dyn)) {
    basic_info$move_dyn
  } else {
    0 # Default to no meta-population if not specified
  }
  
  # Handling age_mu_devs based on onto_move
  data$age_mu_devs <- if (!is.null(basic_info$onto_move) && basic_info$onto_move == 5) {
    if (!is.null(basic_info$age_mu_devs)) {
      basic_info$age_mu_devs
    } else {
      stop("When onto_move is set to 5, age_mu_devs must be defined by the user.")
    }
  } else {
    array(0, dim = c(data$n_stocks, data$n_regions, data$n_ages))
  }
  
  data$spawn_seasons <- rep(1, data$n_stocks)
  if (!is.null(basic_info$spawn_seasons)) {
    if (length(basic_info$spawn_seasons) != data$n_stocks) stop("length of basic_info$spawn_seasons not equal to n_stocks")
    if (any(!(basic_info$spawn_seasons %in% 1:data$n_seasons))) stop("one or more specified spawning seasons are not in 1:n_seasons")
    data$spawn_seasons <- basic_info$spawn_seasons
  }
  
  data$spawn_regions <- rep(1, data$n_stocks)
  if (!is.null(basic_info$spawn_regions)) {
    if (length(basic_info$spawn_regions) != data$n_stocks) stop("length of basic_info$spawn_regions not equal to n_stocks")
    if (any(!(basic_info$spawn_regions %in% 1:data$n_regions))) stop("one or more specified spawning regions are not in 1:n_regions")
    data$spawn_regions <- basic_info$spawn_regions
  }
  
  # Handling recruitment and movement trends
  data$apply_re_trend <- if (!is.null(basic_info$apply_re_trend)) {
    basic_info$apply_re_trend
  } else {
    0 # Default to no recruitment trend if not specified
  }
  
  data$trend_re_rate <- if (!is.null(basic_info$trend_re_rate)) {
    basic_info$trend_re_rate
  } else {
    0 # Default to no recruitment trend rate if not specified
  }
  
  data$apply_mu_trend <- if (!is.null(basic_info$apply_mu_trend)) {
    basic_info$apply_mu_trend
  } else {
    0 # Default to no movement trend if not specified
  }
  
  data$trend_mu_rate <- if (!is.null(basic_info$trend_mu_rate)) {
    basic_info$trend_mu_rate
  } else {
    0 # Default to no movement trend rate if not specified
  }
  
  # Additional logic for other parameters as needed...
  input$data <- data
  return(input)
}
