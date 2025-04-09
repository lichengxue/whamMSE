#' Generate Basic Information for Separate Independent Estimation Models
#'
#' This helper function generates basic biological and fishery information
#' for simulation-estimation and management strategy evaluation (MSE).
#' It builds on `generate_basic_info` by allowing parameter values
#' to be inherited from an existing estimation model (EM) setup.
#'
#' @param em_info List. Contains parameter inputs from a previously defined estimation model.
#' @param em_years Numeric vector. Years to be used for the estimation model.
#' @param n_stocks Integer. (Optional) Number of stocks; defaults to `em_info$par_inputs$n_stocks`.
#' @param n_regions Integer. (Optional) Number of regions; defaults to `em_info$par_inputs$n_regions`.
#' @param n_fleets Integer. (Optional) Number of fleets; defaults to `em_info$par_inputs$n_fleets`.
#' @param n_indices Integer. (Optional) Number of indices; defaults to `em_info$par_inputs$n_indices`.
#' @param n_seasons Integer. (Optional) Number of seasons; defaults to `em_info$par_inputs$n_seasons`.
#' @param filter_indices Integer (0/1) vector (optional). User-specified which indices are excluded from the assessment model 
#' 
#' @return A list containing the generated biological and fishery information,
#'         formatted for input into WHAM or other stock assessment models.
#'
#' @export
#'
#' @seealso `generate_basic_info`, `prepare_wham_input`
#'
generate_basic_info_em_SEP <- function(em_info, 
                                   em_years, 
                                   n_stocks = NULL, 
                                   n_regions = NULL, 
                                   n_fleets = NULL, 
                                   n_indices = NULL, 
                                   filter_indices = NULL) {
  
  if (is.null(filter_indices) || all(filter_indices != 0)) {
    # Default to values from em_info if not provided
    if (is.null(n_stocks)) n_stocks <- em_info$par_inputs$n_stocks
    if (is.null(n_regions)) n_regions <- em_info$par_inputs$n_regions
    if (is.null(n_fleets)) n_fleets <- em_info$par_inputs$n_fleets
    if (is.null(n_indices)) n_indices <- em_info$par_inputs$n_indices
    
    # Generate basic info using provided or default values from em_info
    basic_info <- generate_basic_info(
      n_stocks = n_stocks,
      n_regions = n_regions,
      n_indices = n_indices,
      n_fleets = n_fleets,
      n_seasons = em_info$par_inputs$n_seasons,
      base.years = em_years,
      life_history = em_info$par_inputs$life_history,
      n_ages = em_info$par_inputs$n_ages,
      
      # Using F_info from em_info to construct F_info list
      F_info = list(
        F.year1 = em_info$par_inputs$F.year1,
        Fhist = em_info$par_inputs$Fhist,
        Fmax = em_info$par_inputs$Fmax,
        Fmin = em_info$par_inputs$Fmin,
        change_time = em_info$par_inputs$change_time,
        user_F = NULL
      ),
      
      # Construct catch_info list
      catch_info = list(
        n_fleets = em_info$par_inputs$n_fleets,
        agg_catch = em_info$par_inputs$agg_catch,
        catch_paa = em_info$par_inputs$catch_paa,
        catch_cv = em_info$par_inputs$catch_cv,
        catch_Neff = em_info$par_inputs$catch_Neff,
        use_agg_catch = em_info$par_inputs$use_agg_catch, 
        use_catch_paa = em_info$par_inputs$use_catch_paa
      ),
      
      # Construct index_info list
      index_info = list(
        n_indices <- em_info$par_inputs$n_indices,
        agg_indices = em_info$par_inputs$agg_indices, 
        index_paa = em_info$par_inputs$index_paa,
        index_cv = em_info$par_inputs$index_cv,
        index_Neff = em_info$par_inputs$index_Neff,
        fracyr_indices = em_info$par_inputs$fracyr_indices,
        q = em_info$par_inputs$q,
        use_indices = em_info$par_inputs$use_indices, 
        use_index_paa = em_info$par_inputs$use_index_paa, 
        units_indices = em_info$par_inputs$units_indices, 
        units_index_paa = em_info$par_inputs$units_index_paa
      ),
      
      fracyr_spawn = em_info$par_inputs$fracyr_spawn,
      fracyr_seasons = em_info$par_inputs$fracyr_seasons, # Pass user-defined season fractions if available
      fleet_regions = em_info$par_inputs$fleet_regions, # Fleet region allocation
      index_regions = em_info$par_inputs$index_regions, # Index region allocation
      user_waa = em_info$par_inputs$user_waa, # Use user-defined weight-at-age if provided
      user_maturity = em_info$par_inputs$user_maturity, # Use user-defined maturity-at-age if provided
      bias.correct.process = em_info$par_inputs$bias_correct_process,
      bias.correct.observation = em_info$par_inputs$bias_correct_observation,
      bias.correct.BRPs = em_info$par_inputs$bias_correct_BRPs,
      mig_type = em_info$par_inputs$mig_type,
      XSPR_R_opt = em_info$par_inputs$XSPR_R_opt,
      
      # Add meta-population and movement parameters if specified
      move_dyn = em_info$par_inputs$move_dyn,
      onto_move = em_info$par_inputs$onto_move,
      onto_move_pars = em_info$par_inputs$onto_move_pars,
      
      # Add trend options if specified
      apply_re_trend = em_info$par_inputs$apply_re_trend,
      trend_re_rate = em_info$par_inputs$trend_re_rate,
      apply_mu_trend = em_info$par_inputs$apply_mu_trend,
      trend_mu_rate = em_info$par_inputs$trend_mu_rate,
      age_mu_devs = em_info$par_inputs$age_mu_devs
    ) 
    
    basic_info$fleets_to_remove = em_info$par_inputs$fleets_to_remove
    basic_info$indices_to_remove = em_info$par_inputs$indices_to_remove 
    
  } else {
    
    if (any(filter_indices == 0)) {
      
      # Default to values from em_info if not provided
      if (is.null(n_stocks)) n_stocks <- em_info$par_inputs$n_stocks
      if (is.null(n_regions)) n_regions <- em_info$par_inputs$n_regions
      if (is.null(n_fleets)) n_fleets <- em_info$par_inputs$n_fleets
      
      n_indices_om = em_info$par_inputs$n_indices
      
      # Generate basic info using provided or default values from em_info
      basic_info <- generate_basic_info(
        n_stocks = n_stocks,
        n_regions = n_regions,
        n_indices = em_info$par_inputs$n_indices,
        n_fleets = n_fleets,
        n_seasons = em_info$par_inputs$n_seasons,
        base.years = em_years,
        life_history = em_info$par_inputs$life_history,
        n_ages = em_info$par_inputs$n_ages,
        
        # Using F_info from em_info to construct F_info list
        F_info = list(
          F.year1 = em_info$par_inputs$F.year1,
          Fhist = em_info$par_inputs$Fhist,
          Fmax = em_info$par_inputs$Fmax,
          Fmin = em_info$par_inputs$Fmin,
          change_time = em_info$par_inputs$change_time,
          user_F = NULL
        ),
        
        catch_info = list(
          n_fleets = em_info$par_inputs$n_fleets,
          agg_catch = em_info$par_inputs$agg_catch,
          catch_paa = em_info$par_inputs$catch_paa,
          catch_cv = em_info$par_inputs$catch_cv,
          catch_Neff = em_info$par_inputs$catch_Neff,
          use_agg_catch = em_info$par_inputs$use_agg_catch, 
          use_catch_paa = em_info$par_inputs$use_catch_paa
        ),
        
        # Construct index_info list
        index_info = list(
          n_indices <- em_info$par_inputs$n_indices,
          agg_indices = em_info$par_inputs$agg_indices, 
          index_paa = em_info$par_inputs$index_paa,
          index_cv = em_info$par_inputs$index_cv,
          index_Neff = em_info$par_inputs$index_Neff,
          fracyr_indices = em_info$par_inputs$fracyr_indices,
          q = em_info$par_inputs$q,
          use_indices = em_info$par_inputs$use_indices, 
          use_index_paa = em_info$par_inputs$use_index_paa, 
          units_indices = em_info$par_inputs$units_indices, 
          units_index_paa = em_info$par_inputs$units_index_paa
        ),
        
        fracyr_spawn = em_info$par_inputs$fracyr_spawn,
        fracyr_seasons = em_info$par_inputs$fracyr_seasons, # Pass user-defined season fractions if available
        fleet_regions = em_info$par_inputs$fleet_regions, # Fleet region allocation
        index_regions = em_info$par_inputs$index_regions, # Index region allocation
        user_waa = em_info$par_inputs$user_waa, # Use user-defined weight-at-age if provided
        user_maturity = em_info$par_inputs$user_maturity, # Use user-defined maturity-at-age if provided
        bias.correct.process = em_info$par_inputs$bias_correct_process,
        bias.correct.observation = em_info$par_inputs$bias_correct_observation,
        bias.correct.BRPs = em_info$par_inputs$bias_correct_BRPs,
        mig_type = em_info$par_inputs$mig_type,
        XSPR_R_opt = em_info$par_inputs$XSPR_R_opt,
        
        # Add meta-population and movement parameters if specified
        move_dyn = em_info$par_inputs$move_dyn,
        onto_move = em_info$par_inputs$onto_move,
        onto_move_pars = em_info$par_inputs$onto_move_pars,
        
        # Add trend options if specified
        apply_re_trend = em_info$par_inputs$apply_re_trend,
        trend_re_rate = em_info$par_inputs$trend_re_rate,
        apply_mu_trend = em_info$par_inputs$apply_mu_trend,
        trend_mu_rate = em_info$par_inputs$trend_mu_rate,
        age_mu_devs = em_info$par_inputs$age_mu_devs
      )
      
      # Double check this
      basic_info$fleets_to_remove = em_info$par_inputs$fleets_to_remove
      basic_info$indices_to_remove = em_info$par_inputs$indices_to_remove 
      
    } else {
      stop("filter_indices is not specified correctly!")
    }
  }
  
  return(basic_info)
}
