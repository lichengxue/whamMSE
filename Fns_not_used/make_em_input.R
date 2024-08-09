#' Generate input data for the assessment model
#' 
#' A function to generate the input for the estimation model for management strategy evaluation. 
#' 
#' @param om Operating model 
#' @param info A list of information used to generate the operating model
#' @param M_em Natural mortality random effects
#' @param sel_em Selectivity random effects   
#' @param NAA_re_em Numbers-at-age random effects 
#' @param move_em Movement random effects
#' @param em.opt Movement random effects
#'   \itemize{
#'     \item \code{$separate.em} TRUE = No Global SPR, FALSE = Global SPR
#'     \item \code{$separate.em.type} only if separate.em = TRUE \cr
#'     {=1} panmictic (spatially-aggregated) \cr
#'     {=2} fleets-as-areas \cr
#'     {=3} n single assessment models (n = n_regions) \cr
#'     \item \code{$do.move} T/F movement is included (use if separate.em = FALSE)
#'     \item \code{$est.move} T/F movement rate is estimated (use if separate.em = FALSE)
#'     }
#' @param em_years Years used in the assessment model
#' @param year.use Number of years used in the assessment model
#' @param age_comp_em Likelihood distribution of age composition data in the assessment model
#'   \itemize{
#'     \item \code{"multinomial"} Default
#'     \item \code{"dir-mult"}
#'     \item \code{"dirichlet-miss0"}
#'     \item \code{"dirichlet-pool0"}
#'     \item \code{"dir-mult"}
#'     \item \code{"logistic-normal-miss0"}    
#'     \item \code{"logistic-normal-ar1-miss0"}
#'     \item \code{"logistic-normal-pool0"}
#'     \item \code{"logistic-normal-01-infl"}
#'     \item \code{"logistic-normal-pool0"}
#'     \item \code{"logistic-normal-01-infl-2par"}
#'     \item \code{"mvtweedie"}
#'     \item \code{"dir-mult-linear"}
#'     }
#'     
#' @return a wham input 
#' 
#' @export
#' 
#' @seealso \code{\link{loop_through_fn}}
#' @examples
#' \dontrun{
#' data <- generate_basic_info()
#' input <- prepare_wham_input(basic_info = data)
#' mod <- fit_wham(input, do.fit = FALSE)
#' data <- make_em_input(mod, em.opt = list(separate.em = TRUE, separate.em.type = 1, do.move = FALSE, est.move = FALSE), year.use = 10, em_years = 1973:2022)
#' }

make_em_input <- function(om, 
                          info,
                          M_em, 
                          sel_em, 
                          NAA_re_em, 
                          move_em,
                          em.opt,
                          em_years,
                          year.use,
                          age_comp_em = "multinomial") {
  
  if (is.null(em.opt)) stop("em.opt must be specified!")
  if (em.opt$separate.em) em.opt$do.move <- FALSE 
  if (!em.opt$separate.em & !em.opt$do.move) move.type <- 3 # no movement
  if (!em.opt$separate.em & em.opt$do.move & all(move_em$stock_move)) move.type <- 2 # bidirectional
  if (!em.opt$separate.em & em.opt$do.move & !all(move_em$stock_move)) move.type <- 1 # unidirectional
  
  data <- om$input$data
  
  if (!is.null(year.use)) {
    if (year.use > length(em_years)) {
      cat("year.use must be <= em_years!\nyear.use is set to the length of em_years\n")
      year.use <- length(em_years)
    }
    ind_em <- (length(em_years) - year.use + 1):length(em_years)
    em_years <- tail(em_years, year.use)
  } else {
    year.use <- length(em_years)
    ind_em <- (length(em_years) - year.use + 1):length(em_years)
  }
  
  if (em.opt$separate.em) {
    
    if (em.opt$separate.em.type == 1) { 
      
      info <- generate_basic_info(n_stocks = 1, 
                                  n_regions = 1, 
                                  n_indices = 1, 
                                  n_fleets = 1, 
                                  base.years = em_years,
                                  life_history = info$par_inputs$life_history, 
                                  n_ages= info$par_inputs$n_ages, 
                                  Fbar_ages = info$par_inputs$Fbar_ages, 
                                  recruit_model = info$par_inputs$recruit_model, 
                                  F_info = list(F.year1 = info$par_inputs$F.year1, Fhist = info$par_inputs$Fhist, Fmax = info$par_inputs$Fmax, Fmin = info$par_inputs$Fmin, change_time = info$par_inputs$change_time),
                                  catch_info = list(catch_cv = info$par_inputs$catch_cv,catch_Neff = info$par_inputs$catch_Neff), 
                                  index_info = list(index_cv = info$par_inputs$index_cv,index_Neff = info$par_inputs$index_Neff, fracyr_indices = info$par_inputs$fracyr_indices, q = info$par_inputs$q), 
                                  fracyr_spawn = info$par_inputs$fracyr_spawn, 
                                  bias.correct.process = info$par_inputs$bias.correct.process, 
                                  bias.correct.observation = info$par_inputs$bias.correct.observation, 
                                  bias.correct.BRPs= info$par_inputs$bias.correct.BRPs, 
                                  mig_type = info$par_inputs$mig_type)
      
      basic_info <- info$basic_info
      
      # Fill in the data from the operating model simulation
      info$catch_info$agg_catch <- data$agg_catch[ind_em, , drop = FALSE]
      info$catch_info$agg_catch <- matrix(rowSums(info$catch_info$agg_catch), ncol = 1)
      info$index_info$agg_indices <- data$agg_indices[ind_em, , drop = FALSE]
      info$index_info$agg_indices <- matrix(rowSums(info$index_info$agg_indices), ncol = 1)
      
      info$catch_info$catch_paa <- data$catch_paa[, ind_em, , drop = FALSE]
      catch <- data$agg_catch[ind_em, , drop = FALSE]
      ratio <- data$catch_paa[, ind_em, ]
      result <- 0
      for (i in 1:dim(ratio)[1]) {
        tmp <- ratio[i, , ] * catch[, i]
        result <- result + tmp
      }
      result <- t(apply(result, 1, function(row) row / sum(row)))
      info$catch_info$catch_paa <- array(result, dim = c(1, nrow(result), ncol(result)))
      
      info$index_info$index_paa <- data$index_paa[, ind_em, , drop = FALSE]
      catch <- data$agg_indices[ind_em, , drop = FALSE]
      ratio <- data$index_paa[, ind_em, ]
      result <- 0
      for (i in 1:dim(ratio)[1]) {
        tmp <- ratio[i, , ] * catch[, i]
        result <- result + tmp
      }
      result <- t(apply(result, 1, function(row) row / sum(row)))
      info$index_info$index_paa <- array(result, dim = c(1, nrow(result), ncol(result)))
      
      em_input <- prepare_wham_input(basic_info = basic_info, 
                                     selectivity = sel_em, 
                                     M = M_em, 
                                     NAA_re = NAA_re_em, 
                                     move = NULL,
                                     age_comp = age_comp_em,
                                     catch_info = info$catch_info, 
                                     index_info = info$index_info,
                                     F = info$F)
    }
    
    if (em.opt$separate.em.type == 2) { 
      
      n_fleets <- data$n_fleets
      n_indices <- data$n_indices
      
      info <- generate_basic_info(n_stocks = 1, 
                                  n_regions = 1, 
                                  n_indices = n_indices, 
                                  n_fleets = n_fleets, 
                                  base.years = em_years,
                                  life_history = info$par_inputs$life_history, 
                                  n_ages= info$par_inputs$n_ages, 
                                  Fbar_ages = info$par_inputs$Fbar_ages, 
                                  recruit_model = info$par_inputs$recruit_model, 
                                  F_info = list(F.year1 = info$par_inputs$F.year1, Fhist = info$par_inputs$Fhist, Fmax = info$par_inputs$Fmax, Fmin = info$par_inputs$Fmin, change_time = info$par_inputs$change_time),
                                  catch_info = list(catch_cv = info$par_inputs$catch_cv,catch_Neff = info$par_inputs$catch_Neff), 
                                  index_info = list(index_cv = info$par_inputs$index_cv,index_Neff = info$par_inputs$index_Neff, fracyr_indices = info$par_inputs$fracyr_indices, q = info$par_inputs$q), 
                                  fracyr_spawn = info$par_inputs$fracyr_spawn, 
                                  bias.correct.process = info$par_inputs$bias.correct.process, 
                                  bias.correct.observation = info$par_inputs$bias.correct.observation, 
                                  bias.correct.BRPs= info$par_inputs$bias.correct.BRPs, 
                                  mig_type = info$par_inputs$mig_type)
      
      basic_info <- info$basic_info
      
      # Fill in the data from the operating model simulation
      info$catch_info$agg_catch <- data$agg_catch[ind_em, , drop = FALSE]
      info$index_info$agg_indices <- data$agg_indices[ind_em, , drop = FALSE]
      info$catch_info$catch_paa <- data$catch_paa[, ind_em, , drop = FALSE]
      info$index_info$index_paa <- data$index_paa[, ind_em, , drop = FALSE]
      
      em_input <- prepare_wham_input(basic_info = basic_info, 
                                     selectivity = sel_em, 
                                     M = M_em, 
                                     NAA_re = NAA_re_em, 
                                     move = NULL,
                                     age_comp = age_comp_em,
                                     catch_info = info$catch_info, 
                                     index_info = info$index_info,
                                     F = info$F)
    }
    
    if (em.opt$separate.em.type == 3) {
      
      n_stocks <- data$n_stocks
      em_input <- list()
      
      for (s in 1:n_stocks) {
        
        info <- generate_basic_info(n_stocks = 1, 
                                    n_regions = 1, 
                                    n_indices = 1, 
                                    n_fleets = 1, 
                                    base.years = em_years,
                                    life_history = info$par_inputs$life_history, 
                                    n_ages= info$par_inputs$n_ages, 
                                    Fbar_ages = info$par_inputs$Fbar_ages, 
                                    recruit_model = info$par_inputs$recruit_model, 
                                    F_info = list(F.year1 = info$par_inputs$F.year1, Fhist = info$par_inputs$Fhist, Fmax = info$par_inputs$Fmax, Fmin = info$par_inputs$Fmin, change_time = info$par_inputs$change_time),
                                    catch_info = list(catch_cv = info$par_inputs$catch_cv,catch_Neff = info$par_inputs$catch_Neff), 
                                    index_info = list(index_cv = info$par_inputs$index_cv,index_Neff = info$par_inputs$index_Neff, fracyr_indices = info$par_inputs$fracyr_indices, q = info$par_inputs$q), 
                                    fracyr_spawn = info$par_inputs$fracyr_spawn, 
                                    bias.correct.process = info$par_inputs$bias.correct.process, 
                                    bias.correct.observation = info$par_inputs$bias.correct.observation, 
                                    bias.correct.BRPs= info$par_inputs$bias.correct.BRPs, 
                                    mig_type = info$par_inputs$mig_type)
        
        basic_info <- info$basic_info
        
        # Fill in the data from the operating model simulation
        info$catch_info$agg_catch <- data$agg_catch[ind_em, s, drop = FALSE]
        info$index_info$agg_indices <- data$agg_indices[ind_em, s, drop = FALSE]
        info$catch_info$catch_paa <- data$catch_paa[s, ind_em, , drop = FALSE]
        info$index_info$index_paa <- data$index_paa[s, ind_em, , drop = FALSE]
        
        em_input[[s]] <- prepare_wham_input(basic_info = basic_info, 
                                            selectivity = sel_em, 
                                            M = M_em, 
                                            NAA_re = NAA_re_em, 
                                            move = NULL,
                                            age_comp = age_comp_em,
                                            catch_info = info$catch_info, 
                                            index_info = info$index_info,
                                            F = info$F)
      }
    }
    
  } else {
    
    info <- generate_basic_info(n_stocks = info$par_inputs$n_stocks, 
                                n_regions = info$par_inputs$n_regions, 
                                n_indices = info$par_inputs$n_indices, 
                                n_fleets = info$par_inputs$n_fleets, 
                                base.years = em_years,
                                life_history = info$par_inputs$life_history, 
                                n_ages= info$par_inputs$n_ages, 
                                Fbar_ages = info$par_inputs$Fbar_ages, 
                                recruit_model = info$par_inputs$recruit_model, 
                                F_info = list(F.year1 = info$par_inputs$F.year1, Fhist = info$par_inputs$Fhist, Fmax = info$par_inputs$Fmax, Fmin = info$par_inputs$Fmin, change_time = info$par_inputs$change_time),
                                catch_info = list(catch_cv = info$par_inputs$catch_cv,catch_Neff = info$par_inputs$catch_Neff), 
                                index_info = list(index_cv = info$par_inputs$index_cv,index_Neff = info$par_inputs$index_Neff, fracyr_indices = info$par_inputs$fracyr_indices, q = info$par_inputs$q), 
                                fracyr_spawn = info$par_inputs$fracyr_spawn, 
                                bias.correct.process = info$par_inputs$bias.correct.process, 
                                bias.correct.observation = info$par_inputs$bias.correct.observation, 
                                bias.correct.BRPs= info$par_inputs$bias.correct.BRPs, 
                                mig_type = info$par_inputs$mig_type)
    
    basic_info <- info$basic_info
    
    basic_info <- generate_NAA_where(basic_info = basic_info, move.type = move.type)
    
    # Fill in the data from the operating model simulation
    info$catch_info$agg_catch <- data$agg_catch[ind_em, , drop = FALSE]
    info$index_info$agg_indices <- data$agg_indices[ind_em, , drop = FALSE]
    info$catch_info$catch_paa <- data$catch_paa[, ind_em, , drop = FALSE]
    info$index_info$index_paa <- data$index_paa[, ind_em, , drop = FALSE]
    
    if (em.opt$do.move) {
      
      NAA_re_em$NAA_where = basic_info$NAA_where
      # basic_info$NAA_where = NULL
      
      em_input <- prepare_wham_input(basic_info = basic_info,
                                     selectivity = sel_em,
                                     M = M_em,
                                     NAA_re = NAA_re_em,
                                     move = move_em,
                                     age_comp = age_comp_em,
                                     catch_info = info$catch_info,
                                     index_info = info$index_info,
                                     F = info$F)
      
      # if (sum(em_input$data$NAA_where) == basic_info$n_stocks*basic_info$n_ages) stop("NAA_where is not specified correctly!")
      
      # if (!em.opt$est.move) em_input <- fix_move(em_input)
      
    } else {
      
      basic_info$NAA_where = NULL
      
      NAA_re_em$NAA_where = NULL
      
      em_input <- prepare_wham_input(basic_info = basic_info, 
                                     selectivity = sel_em, 
                                     M = M_em, 
                                     NAA_re = NAA_re_em, 
                                     move = NULL, 
                                     age_comp = age_comp_em,
                                     catch_info = info$catch_info, 
                                     index_info = info$index_info,
                                     F = info$F)
      
      # if (sum(em_input$data$NAA_where) == basic_info$n_stocks*basic_info$n_ages)  stop("NAA_where is not specified correctly!")
      
    }
  }
  
  return(em_input)
}
