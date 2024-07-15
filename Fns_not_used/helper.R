#' Specify catch selectivity blocks and aggregate and age composition observations for catch
#'
#' @param input list containing data, parameters, map, and random elements (output from \code{\link{wham::prepare_wham_input}})
#' @param catch_info (optional) list specifying various aspects about catch by fleet (see details)
#' 
#' \code{catch_info} specifies observations, and various configuration options for fleet-specific catch observations and will overwrite attributes specified in the ASAP data file.
#' If \code{NULL}, all settings from the ASAP data file or basic_info are used.
#' \code{catch_info} is a list with any of the following entries:
#'   \describe{
#'     \item{$n_fleets}{number of fleets}
#'     \item{$fleet_regions}{vector (n_fleets) of regions where each fleet operates.}
#'     \item{$fleet_seasons}{matrix (n_fleets x n_seasons) of 0/1 values flagging which seasons each fleet operates.}
#'     \item{$agg_catch}{matrix (n_years_model x n_fleets) of annual aggregate catches by fleet.}
#'     \item{$agg_catch_cv}{matrix (n_years_model x n_fleets) of CVs for annual aggregate catches by fleet.}
#'     \item{$catch_paa}{array (n_fleets x n_years_model x n_ages) of annual catch proportions at age by fleet.}
#'     \item{$use_catch_paa}{matrix (n_years_model x n_fleets) of 0/1 values flagging whether to use proportions at age observations.}
#'     \item{$catch_Neff}{matrix (n_years_model x n_fleets) of effective sample sizes for proportions at age observations.}
#'     \item{$selblock_pointer_fleets}{matrix (n_years_model x n_fleets) of itegers indicated selblocks to use.}
#'   }
#'
#' @return a named list with same elements as the input provided with catch observations and fleet options modified.
#'
#' @seealso \code{\link{prepare_wham_input}} 
#'
#' @examples
#' \dontrun{
#' wham.dir <- find.package("wham")
#' path_to_examples <- system.file("extdata", package="wham")
#' asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))
#' input <- prepare_wham_input(asap3)
#` newcatch <- matrix(500, input$data$n_years_model, input$data$n_fleets)
#' input <- set_catch(input, catch_info = list(agg_catch = newcatch)) #constant catch of 500 mt
#' }
#'
#' @export
set_catch = function(input, catch_info= NULL) {
  data = input$data
  asap3 = input$asap3
  if(is.null(asap3)){
    data$n_fleets = 1
  } else {
    input$fleet_names <- NULL
    for(i in 1:length(asap3)) asap3[[i]]$use_catch_acomp <- rep(1,asap3[[i]]$n_fleets) #default is to use age comp for catch
    n_fleets_per_region = sapply(asap3, function(x) x$n_fleets)
    data$n_fleets = sum(n_fleets_per_region)
    for(i in 1:length(asap3)) if(!is.null(asap3[[i]]$fleet.names)) input$fleet_names <- c(input$fleet_names, asap3[[i]]$fleet.names)
  }
  if(!is.null(catch_info$n_fleets)) data$n_fleets = catch_info$n_fleets 
  if(is.null(input$fleet_names)) input$fleet_names <- paste0("fleet_", 1:data$n_fleets)
  
  data$fleet_regions = rep(1, data$n_fleets)
  data$fleet_seasons = matrix(1, data$n_fleets, data$n_seasons)
  
  data$agg_catch = data$agg_catch_sigma = data$catch_Neff = matrix(NA, data$n_years_model, data$n_fleets)
  data$catch_paa = array(NA, dim = c(data$n_fleets, data$n_years_model, data$n_ages))
  data$use_agg_catch = matrix(1, data$n_years_model, data$n_fleets)
  data$use_catch_paa = matrix(1, data$n_years_model, data$n_fleets)
  data$selblock_pointer_fleets = matrix(0, data$n_years_model, data$n_fleets)
  
  if(!is.null(asap3)) {
    #data$n_fleets = asap3$n_fleets
    k <- 1
    for(i in 1:length(asap3)) {
      asap3[[i]]$use_catch_acomp <- rep(1,asap3[[i]]$n_fleets) #default is to use age comp for catch
      for(j in 1:asap3[[i]]$n_fleets) {
        data$fleet_regions[k] = i #each asap file is a separate region
        data$agg_catch[,k] = asap3[[i]]$CAA_mats[[j]][,data$n_ages + 1]
        data$agg_catch_sigma[,k] = asap3[[i]]$catch_cv[,j]
        temp = asap3[[i]]$CAA_mats[[j]][,1:data$n_ages]
        temp[which(is.na(temp))] = 0
        temp[which(temp<0)] = 0
        data$catch_paa[k,,] = temp/apply(temp,1,sum)
        for(y in 1:data$n_years_model) if(asap3[[i]]$CAA_mats[[j]][y,data$n_ages+1] < 1e-15) data$use_agg_catch[y,k] = 0
        if(asap3[[i]]$use_catch_acomp[j] != 1){
          data$use_catch_paa[,k] = 0
        } else { # use catch paa in at least some years - not necessarily all, have to go through year by year
          for(y in 1:data$n_years_model){
            if(is.na(sum(data$catch_paa[k,y,] > 1e-15))){ # handle negative or NA paa
              data$use_catch_paa[y,k] = 0
            } else {
              if(asap3[[i]]$catch_Neff[y,j] < 1e-15 | sum(data$catch_paa[k,y,] > 1e-15)<2) data$use_catch_paa[y,k] = 0
            }
          } 
        }
        data$catch_Neff[,k] = asap3[[i]]$catch_Neff[,j]
        temp = asap3[[i]]$sel_block_assign[[j]]
        temp = match(temp,sort(unique(temp))) #index unique values
        data$selblock_pointer_fleets[,k] = max(data$selblock_pointer_fleets) + temp #max grows each time
        k <- k + 1
      }
    }
    data$agg_catch_sigma[which(data$agg_catch_sigma < 1e-15)] = 100
    data$agg_catch_sigma = sqrt(log(data$agg_catch_sigma^2 + 1))
  }
  else
  {
    #data$n_fleets = 1
    data$agg_catch[] = 1000  	
    data$catch_paa[] = 1/data$n_ages
    data$agg_catch_sigma[] = sqrt(log(0.1^2 + 1))
    data$catch_Neff[] = 200	  
    for(i in 1:data$n_fleets) for(y in 1:data$n_years_model){ 
      if(data$catch_Neff[y,i] < 1e-15 | sum(data$catch_paa[i,y,] > 1e-15)<2 | any(is.na(data$catch_paa[i,y,]))) data$use_catch_paa[y,i] = 0
    }
    data$selblock_pointer_fleets[] = rep(1:data$n_fleets, each = data$n_years_model)
    input$fleet_names <- paste0("Fleet ", 1:data$n_fleets)
  }
  
  if(!is.null(catch_info$fleet_regions)) data$fleet_regions[] = catch_info$fleet_regions
  if(!is.null(catch_info$agg_catch)) data$agg_catch[] = catch_info$agg_catch
  if(!is.null(catch_info$catch_paa)) data$catch_paa[] = catch_info$catch_paa
  if(!is.null(catch_info$catch_cv)) data$agg_catch_sigma[] = sqrt((log(catch_info$catch_cv^2 + 1)))
  if(!is.null(catch_info$catch_Neff)) data$catch_Neff[] = catch_info$catch_Neff
  if(!is.null(catch_info$use_catch_paa)) data$use_catch_paa[] = catch_info$use_catch_paa
  if(!is.null(catch_info$selblock_pointer_fleets)) data$selblock_pointer_fleets[] = catch_info$selblock_pointer_fleets
  
  data$catch_paa[is.na(data$catch_paa)] = 0
  
  input$par$log_catch_sig_scale = rep(0, data$n_fleets)
  input$map$log_catch_sig_scale = factor(rep(NA, data$n_fleets))
  
  input$data = data
  input$asap3 <- asap3
  input$options$catch <- catch_info
  if(!is.null(input$par$logit_selpars)) input <- set_selectivity(input, input$options$selectivity)
  if(!is.null(input$data$obsvec)) input <- set_osa_obs(input)
  
  return(input)
}

#' Specify configuration for environmental covariates, effects on the population, and parameter values
#'
#' @param input list containing data, parameters, map, and random elements (output from \code{\link{wham::prepare_wham_input}})
#' @param ecov (optional) named list of environmental covariate data and parameters (see details)
#' 
#' \code{ecov} specifies any environmental covariate data and model as well as the effect on the population. Environmental covariate data need not span
#' the same years as the fisheries data. It can be \code{NULL} if no environmental data are to be fit.
#' Otherwise, it must be a named list with the following components:
#'   \describe{
#'     \item{$label}{Name(s) of the environmental covariate(s). Used in printing.}
#'     \item{$mean}{Mean observations (matrix). number of years x number of covariates. Missing values = NA.}
#'     \item{$logsigma}{Configure observation standard errors. Options:
#'       \describe{
#'         \item{Matrix of \eqn{log} standard errors with same dimensions as \code{$mean}}{Specified values for each time step }
#'         \item{log standard errors for each covariate, numeric vector or matrix w/ dim 1 x n.ecov}{Specified value the same for all time steps}
#'         \item{estimation option (for all covariates). character string:}{
#'           \code{"est_1"}: Estimated, one value shared among time steps.
#'           \code{"est_re"}: Estimated value for each time step as random effects with two parameters (mean, var)}
#'         \item{list of two elements.}{
#'           First is the matrix of log standard errors or the vector of single values for each covariate as above. 
#'           Second is a character vector of estimation options (\code{NA}, \code{"est_1"},\code{"est_re"}) for each covariate. 
#'           For covariates with non-NA values, values in the first element are ignored.}
#'       }
#'     }
#'     \item{$year}{Years corresponding to observations (vector of same length as \code{$mean} and \code{$logsigma})}
#'     \item{$use_obs}{T/F (or 0/1) vector/matrix of the same dimension as \code{$mean} and \code{$logsigma}.
#'     Use the observation? Can be used to ignore subsets of the ecov without changing data files.}
#'     \item{$process_model}{Process model for the ecov time-series. \code{"rw"} = random walk, \code{"ar1"} = 1st order autoregressive, 
#'        \code{NA} = do not fit}
#'     \item{$process_mean_vals}{vector of (initial) mean values for the ecov time-series.}
#'     \item{$process_sig_vals}{vector of (initial) standard deviation values for the ecov time-series.}
#'     \item{$process_cor_vals}{vector of (initial) correlation values for the ecov time-series.}
#'     \item{$recruitment_how}{character matrix (n_Ecov x n_stocks) indicating how each ecov affects recruitment for each stock. 
#'        Options are based on (see \href{https://www.sciencedirect.com/science/article/pii/S1385110197000221}{Iles & Beverton (1998)}) 
#'        combined with the order of orthogonal polynomial of the covariate and has the form "type-lag-order". "type" can be:
#'         \describe{
#'           \item{= "none"}{no effect.}
#'           \item{= "controlling"}{pre-recruit density-independent mortality.}
#'           \item{= "limiting"}{ maximum recruitment, e.g. ecov determines amount of suitable habitat)}
#'           \item{= "lethal"}{threshold, i.e. R --> 0 at some ecov value.}
#'           \item{= "masking"}{metabolic/growth, decreases dR/dS}
#'           \item{= "directive"}{e.g. behavioral}
#'          }
#'         for type other than "none", "lag" can be:
#'         \describe{
#'           \item{= "lag-n"}{lag = n which can be 0,1,2,.... lag-1 implies the covariate in year y affects recruitment in year y+1.}
#'         }
#'         for "type" being other than "none", "order" can be:
#'         \describe{
#'           \item{= "linear"}{the covariate effect is linear on the transformed recruitment parameter (e.g., log).}
#'           \item{= "poly-n"}{orthogonal polynomial where n = 1 (same as "linear"),2,...}
#'         }
#'         so "limiting-lag-1-poly-2" would model the covariate affecting recruitment the next year (lag = 1) as a second order orthogonal 
#'          polynomial (\eqn{b_0 + b_1*ecov + b_2*ecov^2 + ...}) limiting effect.
#'      }
#'     \item{$M_how}{character array (n_Ecov x n_stocks x n_ages x n_regions) indicating how each ecov affects M by age,stock,region and 
#'        has the form "lag-order". "lag" can be:
#'         \describe{
#'           \item{= "none"}{no effect.}
#'           \item{= "lag-n"}{lag = n which can be 0,1,2,.... lag-1 implies the covariate in year y affects M in year y+1.}
#'         }
#'         for "lag" being other than "none", "order" can be:
#'         \describe{
#'           \item{= "linear"}{the covariate effect is linear on the transformed M parameter (e.g., log).}
#'           \item{= "poly-n"}{orthogonal polynomial where n = 1 (same as "linear"),2,...}
#'          }
#'       }
#'     \item{$M_effect_map}{integer array (n_stocks x n_ages x n_regions x n_Ecov) indicating which estimated effects are common by age,stock,region.
#'       If not specified there the same effect is estimated for all M where $M_how is other than "none" for each covariate.}
#'     \item{$q_how}{character matrix (n_Ecov x n_indices) indicating whether each ecov affects catchability for each index. and has 
#'      the form "lag-order". "lag" can be:
#'         \describe{
#'           \item{= "none"}{no effect.}
#'           \item{= "lag-n"}{lag = n which can be 0,1,2,.... lag-1 implies the covariate in year y affects catchability in year y+1.}
#'         }
#'         for "lag" being other than "none", "order" can be:
#'         \describe{
#'           \item{= "linear"}{the covariate effect is linear on the transformed catchability parameter (e.g., log).}
#'           \item{= "poly-n"}{orthogonal polynomial where n = 1 (same as "linear"),2,...}
#'          }
#'       }
#'     \item{$move_how}{character array (n_Ecov x n_stocks x n_ages x n_seasons x n_regions x n_regions - 1) indicating whether each ecov 
#'        affects movement from one region to the others by stock,age,season. and has the form "lag-order". "lag" can be:
#'         \describe{
#'           \item{= "none"}{no effect.}
#'           \item{= "lag-n"}{lag = n which can be 0,1,2,.... lag-1 implies the covariate in year y affects a movement parameter in year y+1.}
#'         }
#'         for "lag" being other than "none", "order" can be:
#'         \describe{
#'           \item{= "linear"}{the covariate effect is linear on the transformed movement parameter (e.g., log).}
#'           \item{= "poly-n"}{orthogonal polynomial where n = 1 (same as "linear"),2,...}
#'          }
#'        }
#'     \item{$move_effect_map}{integer array (n_stocks x n_ages x n_seasons x n_regions x n_regions-1 x n_Ecov) indicating which estimated 
#'        effects are common by age,stock,region, season etc. If not specified the same effect is estimated for all movement parameters
#'        where $move_how is other than "none" for each covariate.}
#'     \item{$beta_R_vals}{n_stocks x n_ecov x max(n_poly_R) array of initial values for effects on recruitment.}
#'     \item{$beta_M_vals}{n_stocks x n_ages x n_regions x n_ecov x max(n_poly_M) array of initial values for effects on natural mortality.}
#'     \item{$beta_q_vals}{n_indices x n_ecov x max(n_poly_q) array of initial values for effects on catchability.}
#'     \item{$beta_mu_vals}{n_stocks x n_ages x n_seasons x n_regions x n_regions - 1 x n_ecov x max(n_poly_move) array of initial values for effects on movement parameters.}
#'   }
#'
#' @return a named list with same elements as the input provided with environmental covariate observations, effects, and model options modified.
#'
#' @seealso \code{\link{prepare_wham_input}} 
#'
#' @examples
#' \dontrun{
#' wham.dir <- find.package("wham")
#' path_to_examples <- system.file("extdata", package="wham")
#' asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))
#' env.dat <- read.csv(file.path(path_to_examples,"GSI.csv"), header=T)
#' input <- prepare_wham_input(asap3, NAA_re = list(sigma = "rec"))
#' ecov <- list(
#'  label = "GSI",
#'  mean = as.matrix(env.dat$GSI),
#'  logsigma = 'est_1', # estimate obs sigma, 1 value shared across years
#'  year = env.dat$year,
#'  use_obs = matrix(1, ncol=1, nrow=dim(env.dat)[1]), # use all obs (=1)
#'  process_model = 'ar1') # "rw" or "ar1"
#' input <- set_ecov(input, ecov = ecov) #GSI in the model without any effects
#' }
#'
#' @export
set_ecov = function(input, ecov) {
  data = input$data
  par = input$par
  map = input$map
  input$log$ecov <- list()
  #clear any map definitions that may exist. necessary because some configurations may not define map elements.
  ecov_pars = c("Ecov_re", "Ecov_beta_R", "Ecov_beta_M", "Ecov_beta_mu", "Ecov_beta_q", 
                "Ecov_process_pars", "Ecov_obs_log_sigma", "Ecov_obs_logsigma_re", "Ecov_obs_sigma_par")
  map <- map[(!names(map) %in% ecov_pars)]
  
  #define new dimensions for all effects a given ecov can have on assessment model
  #currently the order is recruitment, M, index 1, ..., n_indices
  #n_effects = 2 + data$n_indices
  #index_effects = 2+1:data$n_indices
  
  # --------------------------------------------------------------------------------
  # Environmental covariate data
  #set up for ecov == NULL
  data$Ecov_obs <- matrix(1, nrow=1, ncol=1)
  data$Ecov_obs_sigma_opt <- 1
  data$n_Ecov <- 1
  data$Ecov_model <- rep(0, data$n_Ecov)
  input$years_Ecov <- input$years[1]
  data$n_years_Ecov <- 1
  data$years_use_Ecov <- 0
  data$Ecov_use_obs <- matrix(0, nrow=1, ncol=1)
  
  data$Ecov_how_R <- matrix(0, data$n_Ecov, data$n_stocks)
  data$Ecov_how_q <- matrix(0, data$n_Ecov, data$n_indices)
  data$Ecov_how_M <- array(0, dim = c(data$n_Ecov, data$n_stocks, data$n_ages,data$n_regions))
  data$Ecov_how_mu <- array(0, dim = c(data$n_Ecov, data$n_stocks, data$n_ages, data$n_seasons, data$n_regions, data$n_regions-1))
  data$ind_Ecov_out_start_R <- data$ind_Ecov_out_end_R <- matrix(0, data$n_Ecov, data$n_stocks)
  data$ind_Ecov_out_start_q <- data$ind_Ecov_out_end_q <- matrix(0, data$n_Ecov, data$n_indices)
  data$ind_Ecov_out_start_M <- data$ind_Ecov_out_end_M <- array(0, dim = c(data$n_Ecov, data$n_stocks, data$n_ages, data$n_regions))
  data$ind_Ecov_out_start_mu <- data$ind_Ecov_out_end_mu <- array(0, dim = c(data$n_Ecov, data$n_stocks, data$n_ages, data$n_seasons, 
                                                                             data$n_regions, data$n_regions-1))
  
  input$Ecov_names <- "none"
  data$Ecov_use_re <- rep(0, data$n_Ecov)
  
  data$n_poly_Ecov_R <- matrix(1,data$n_Ecov, data$n_stocks)
  data$n_poly_Ecov_M <- array(1,dim = c(data$n_Ecov, data$n_stocks, data$n_ages, data$n_regions))
  data$n_poly_Ecov_mu <- array(1,dim = c(data$n_Ecov, data$n_stocks, data$n_ages, data$n_seasons, data$n_regions, data$n_regions-1))
  #print(data$n_Ecov)
  #print(dim(data$n_poly_Ecov_mu))
  #print(max(data$n_poly_Ecov_mu))
  #print(data$n_poly_Ecov_mu)
  data$n_poly_Ecov_q <- matrix(1,data$n_Ecov, data$n_indices)
  
  par$Ecov_obs_logsigma <- matrix(-2.3, nrow=1, ncol=1)
  par$Ecov_obs_logsigma_re <- matrix(-2.3, nrow=1, ncol=1)
  par$Ecov_obs_sigma_par <- matrix(0, nrow=1, ncol=1)
  par$Ecov_process_pars <- matrix(0, 3, data$n_Ecov)  
  par$Ecov_re <- matrix(0, data$n_years_Ecov, data$n_Ecov)
  par$Ecov_beta_R <- array(0, dim = c(data$n_stocks, data$n_Ecov, max(data$n_poly_Ecov_R)))
  par$Ecov_beta_q <- array(0, dim = c(data$n_indices, data$n_Ecov, max(data$n_poly_Ecov_q)))
  par$Ecov_beta_M <- array(0, dim = c(data$n_stocks, data$n_ages, data$n_regions, data$n_Ecov, max(data$n_poly_Ecov_M)))
  par$Ecov_beta_mu <- array(0, dim = c(data$n_stocks, data$n_ages, data$n_seasons, data$n_regions, data$n_regions-1, data$n_Ecov, max(data$n_poly_Ecov_mu,0)))
  map$Ecov_obs_logsigma <- rep(NA, length(par$Ecov_obs_logsigma))
  map$Ecov_obs_logsigma_re <- rep(NA, length(par$Ecov_obs_logsigma_re))
  map$Ecov_obs_sigma_par <- rep(NA, length(par$Ecov_obs_sigma_par))
  map$Ecov_process_pars <- rep(NA, length(par$Ecov_process_pars))
  map$Ecov_re <- rep(NA, length(par$Ecov_re))
  map$Ecov_beta_R <- rep(NA, length(par$Ecov_beta_R))
  map$Ecov_beta_q <- rep(NA, length(par$Ecov_beta_q))
  map$Ecov_beta_M <- rep(NA, length(par$Ecov_beta_M))
  map$Ecov_beta_mu <- rep(NA, length(par$Ecov_beta_mu))
  
  if(!is.null(ecov)){
    
    #FIRST: configure observation models, likelihoods, and process models, likelihoods for ecovs
    if(class(ecov$mean)[1] == "matrix") {data$Ecov_obs <- ecov$mean} else{
      input$log$ecov <- c(input$log$ecov, "NOTE: ecov$mean is not a matrix. Coercing to a matrix...")
      data$Ecov_obs <- as.matrix(ecov$mean)
    }
    if(class(ecov$use_obs)[1] == "matrix"){
      data$Ecov_use_obs <- ecov$use_obs
    } else{
      input$log$ecov <- c(input$log$ecov, "NOTE: ecov$use_obs is not a matrix with same dimensions as ecov$mean. Coercing to a matrix...")
      data$Ecov_use_obs <- as.matrix(as.integer(ecov$use_obs))
    }
    if(!identical(dim(data$Ecov_use_obs), dim(data$Ecov_obs))) stop("Dimensions of ecov$use_obs != dimensions of ecov$mean")
    
    # Handle Ecov sigma options
    data$n_Ecov <- dim(data$Ecov_obs)[2] # num Ecovs
    n_Ecov_obs <- dim(data$Ecov_obs)[1] # num Ecov obs
    data$Ecov_obs_sigma_opt = rep(1, data$n_Ecov) # Ecov sigma given, initialized at given values, not estimated by default
    
    if(length(ecov$year) != n_Ecov_obs) stop("ecov$year is not the same length as # rows in ecov$mean")
    #data$Ecov_year <- as.numeric(ecov$year)
    input$years_Ecov <- as.numeric(ecov$year)
    end_model <- tail(input$years,1)
    end_Ecov <- tail(ecov$year,1)
    
    #define labels and report info about effects of each
    if(length(ecov$label) == data$n_Ecov){
      input$Ecov_names <- ecov$label
    } else {
      input$log$ecov <- c(input$log$ecov, "NOTE: Number of Ecov labels not equal to number of Ecovs
              Setting Ecov labels = 'Ecov 1', 'Ecov 2', ...")
      input$Ecov_names = paste0("Ecov ",1:data$n_Ecov)
    }    
    
    
    #default NAs for parameter matrix of observation standard deviations
    par$Ecov_obs_logsigma = matrix(NA, n_Ecov_obs, data$n_Ecov)
    map$Ecov_obs_logsigma <- matrix(NA, nrow=n_Ecov_obs, ncol=data$n_Ecov) # turn off estimation
    
    logsig_more = list()
    if(is.list(ecov$logsigma)){
      n = length(ecov$logsigma)
      if(n>1) logsig_more = ecov$logsigma[[2]] #further elements than the second are ignored.
      ecov$logsigma = ecov$logsigma[[1]]
    }
    
    if(class(ecov$logsigma)[1] == "matrix"){
      par$Ecov_obs_logsigma <- ecov$logsigma
      if(!identical(dim(par$Ecov_obs_logsigma), dim(data$Ecov_obs))) stop("Dimensions of ecov$mean != dimensions of ecov$logsigma")
    }
    if(class(ecov$logsigma)[1] == 'numeric'){
      #data$Ecov_obs_sigma_opt[] = 1 #defined above
      input$log$ecov <- c(input$log$ecov, "ecov$logsigma is numeric. Coercing to a matrix... \n")
      if(length(ecov$logsigma) == data$n_Ecov) par$Ecov_obs_logsigma <- matrix(rep(ecov$logsigma, each=n_Ecov_obs), ncol=data$n_Ecov)
      if(length(ecov$logsigma) == n_Ecov_obs & data$n_Ecov == 1) par$Ecov_obs_logsigma <- matrix(ecov$logsigma, ncol=1)
      if(length(ecov$logsigma) != data$n_Ecov & length(ecov$logsigma) != n_Ecov_obs) stop("ecov$logsigma is numeric but length is not equal to # of ecovs or ecov observations")
    }
    
    #set up and check length of ecov$process_model
    if(length(ecov$process_model) == 1) ecov$process_model = rep(ecov$process_model, data$n_Ecov) #use the single value for all Ecovs
    if(length(ecov$process_model) != data$n_Ecov) stop("length of ecov$process_model must be either 1 or the number of Ecovs")
    
    #now over write ecov$logsigma with logsig_more if available because the fixed obs var matrices have been defined
    if(length(logsig_more)) ecov$logsigma = logsig_more
    
    
    if(class(ecov$logsigma)[1] == 'character'){
      #check that estimation options are right
      if(!all(ecov$logsigma %in% c("est_1", "est_re", NA))){
        stop("ecov$logsigma or ecov$logsigma[[2]] is character and must be NA (do not estimate), 'est_1' (single variance parameter), or 'est_re' (iid re annual variance parameters)")
      }
      if(length(ecov$logsigma) == 1) ecov$logsigma = rep(ecov$logsigma, data$n_Ecov) #use the single value for all Ecovs
      #check length of estimation options
      if(length(ecov$logsigma) != data$n_Ecov) stop("length of ecov$logsigma when character must be either 1 or the number of Ecovs")
      
      
      for(i in 1:data$n_Ecov) {
        if(!is.na(ecov$logsigma[i])) if(ecov$logsigma[i] == 'est_1'){ # estimate 1 Ecov obs sigma for each Ecov
          data$Ecov_obs_sigma_opt[i] = 2
          par$Ecov_obs_logsigma[,i] <- -1.3 
          map$Ecov_obs_logsigma[,i] <- i 
        }
        if(!is.na(ecov$logsigma[i])) if(ecov$logsigma[i] == 'est_re'){
          data$Ecov_obs_sigma_opt[i] = 4
          map$Ecov_obs_logsigma[,i] <- NA # turn off estimation of fixed effects
          par$Ecov_obs_sigma_par[,i] <- c(-1.3, -2.3) # random effect pars
          map$Ecov_obs_sigma_par[,i] <- max(0, map$Ecov_obs_sigma_par, na.rm =T) + 1:2 
        }
      }
      
      #default values for: 
      #map of observation error variance parameters
      map$Ecov_obs_logsigma_re <- matrix(NA, nrow=n_Ecov_obs, ncol=data$n_Ecov) # turn off estimation
      #initial values of random effects
      par$Ecov_obs_logsigma_re = matrix(0, n_Ecov_obs, data$n_Ecov)
      
      for(i in 1:data$n_Ecov) if(!is.na(ecov$logsigma[i])) if(ecov$logsigma[i] == 'est_re') {
        map$Ecov_obs_logsigma_re[,i] = max(0, map$Ecov_obs_logsigma_re, na.rm=T) + 1:data$n_years_Ecov
        par$Ecov_obs_logsigma_re[,i] <- par$Ecov_obs_logsigma[,i] # random effect initialize at values in matrix provided
      }
      
    }
    
    
    if(!all(ecov$process_model %in% c(NA,"rw", "ar1"))){
      stop("ecov$process_model must be 'rw' (random walk), 'ar1', or NA (do not fit)")
    }
    data$Ecov_model <- sapply(ecov$process_model, match, c("rw", "ar1"))
    data$Ecov_model[is.na(data$Ecov_model)] <- 0 #don't fit for NA
    
    
    if(is.null(ecov$recruitment_how)) ecov$recruitment_how <- matrix("none", data$n_Ecov, data$n_stocks)
    if(is.null(ecov$q_how)) ecov$q_how <- matrix("none", data$n_Ecov, data$n_indices)
    if(is.null(ecov$M_how)) ecov$M_how <- array("none", dim = c(data$n_Ecov, data$n_stocks, data$n_ages, data$n_regions))
    if(is.null(ecov$M_effect_map)) ecov$M_effect_map <- array(NA, dim = c(data$n_stocks, data$n_ages, data$n_regions, data$n_Ecov))
    if(is.null(ecov$move_how)) ecov$move_how <- array("none", dim = c(data$n_Ecov, data$n_stocks, data$n_ages, data$n_seasons, data$n_regions, data$n_regions-1))
    if(is.null(ecov$move_effect_map)) ecov$move_effect_map <- array(NA, dim = c(data$n_stocks, data$n_ages, data$n_seasons, data$n_regions, data$n_regions-1, data$n_Ecov))
    
    for(i in 1:data$n_Ecov) {
      ecov_used <- any(c(ecov$recruitment_how[i,],ecov$q_how[i,],ecov$M_how[i,,,],ecov$move_how[i,,,,,]) != "none")
      if(data$Ecov_model[i] == 0 & ecov_used) {
        stop(paste0("ecov$process_model ", i, " is turned off (NA) but an effect is specified in ecov$recruitment_how, ecov$q_how, ecov$M_how, and/or ecov$move_how.
       Either 1) choose an ecov process model ('rw' or 'ar1'),
              2) turn off ecov (set the ecov$type_how = 'none' and ecov$process_model = NA),
           or 3) fit ecov but with no effect on population (set the ecov$type_how = 'none' and ecov$process_model[i] = 'rw' or 'ar1')."))
      }
    }
    
    #make R
    ecov$lag_R <- matrix(1, data$n_Ecov, data$n_stocks)
    data$n_poly_Ecov_R <- matrix(1, data$n_Ecov, data$n_stocks)
    data$Ecov_how_R <- matrix(0, data$n_Ecov, data$n_stocks)
    if(!is.null(ecov$recruitment_how)) {
      for(s in 1:data$n_stocks) for(i in 1:data$n_Ecov) {
        tmp <- strsplit(ecov$recruitment_how[i,s], split = "-")[[1]]
        if(!(length(tmp) %in% c(1,4:5))) stop(paste0("Form of ecov$recruitment_how[",i,",", s, "] for ecov, ", i, " and stock , ", s, 
                                                     " should be 'none', 'type-lag-l-poly-p' or 'type-lag-l-linear'"))
        if(tmp[1] != "none") ecov$lag_R[i,s] <- as.integer(tmp[3])
        if(length(tmp) == 5) data$n_poly_Ecov_R[i,s] <- as.integer(tmp[5])
        if(tmp[1] != "none" & data$recruit_model[s] == 1) stop(paste0("Random walk recruitment for stock ", s, " cannot have an ecov effect on recruitment.
            Either choose a different recruit_model (2, 3, or 4), or remove the Ecov effect."))
        if(data$recruit_model[s] == 2 & !(tmp[1] %in% c("none", "controlling"))) stop(paste0("Random recruitment about mean for stock ", s, 
                                                                                             "only allows effect type to be 'none' or 'controlling'."))
        if(data$recruit_model[s] == 4 & tmp[1] == "limiting") stop(paste0("'Limiting' ecov effect on Ricker recruitment for stock ", s, " not implemented.
          Either set ecov$how_R = 0 (no effect), 1 (controlling), or 4 (masking)...
          Or set recruit_model = 3 (Bev-Holt)."))
        if(tmp[1] != "none" & data$NAA_re_model[s] == 0) stop(paste0("Cannot estimate ecov effect on recruitment for stock ", s, " when
          recruitment in each year is estimated freely as fixed effect parameters.
          Either remove ecov-recruit effect or estimate recruitment
          (or all numbers-at-age) as random effects."))
        if(!(tmp[1] %in% c("none", "controlling", "limiting", "lethal", "masking", "directive"))) stop(paste0("The first component of each
          element of ecov$recruitment_how must be one of the following: 'none', 'controlling', 'limiting', 'lethal', 'masking', 'directive'"))
        if(tmp[1] == "controlling") data$Ecov_how_R[i,s] <- 1
        if(tmp[1] == "limiting") data$Ecov_how_R[i,s] <- 2
        if(tmp[1] == "lethal") data$Ecov_how_R[i,s] <- 3
        if(tmp[1] == "masking") data$Ecov_how_R[i,s] <- 4
        if(tmp[1] == "directive") data$Ecov_how_R[i,s] <- 5
        #it is already set to 1
        #if(length(tmp)== 4 & tmp[4] == "linear") ecov_poly_R[i,s] <- 1 
      }
    } else {
      input$log$ecov <- c(input$log$ecov, "no ecov$recruitment_how is provided so no effects on recruitment will be estimated. \n")
    }
    
    par$Ecov_beta_R <- array(0, dim = c(data$n_stocks, data$n_Ecov, max(data$n_poly_Ecov_R)))
    map$Ecov_beta_R <- array(NA, dim = dim(par$Ecov_beta_R))
    k <- 0
    for(s in 1:data$n_stocks) for(i in 1:data$n_Ecov) if(data$Ecov_how_R[i,s]>0){
      n <- data$n_poly_Ecov_R[i,s]
      map$Ecov_beta_R[s,i,1:n] <- k + 1:n
      k <- k + n
    }
    
    #make M
    ecov$lag_M <- array(0, dim = c(data$n_Ecov, data$n_stocks, data$n_ages, data$n_regions))
    data$n_poly_Ecov_M <- array(1, dim = c(data$n_Ecov,data$n_stocks, data$n_ages, data$n_regions))
    data$Ecov_how_M <- array(0, dim = c(data$n_Ecov,data$n_stocks, data$n_ages, data$n_regions))
    if(!is.null(ecov$M_how)) {
      for(s in 1:data$n_stocks) for(r in 1:data$n_regions) for(a in 1:data$n_ages) for(i in 1:data$n_Ecov) {
        tmp <- strsplit(ecov$M_how[i,s,a,r], split = "-")[[1]]
        if(!(length(tmp) %in% c(1,3:4))) stop(paste0("Form of ecov$M_how[",i,",", s, ",", a, ",", r, "] for ecov, ", i, ", stock , ", s, 
                                                     ", age , ", a, ", region , ", r," should be 'none', 'lag-l-poly-p' or 'lag-l-linear'"))
        if(tmp[1] != "none"){
          data$Ecov_how_M[i,s,a,r] <- 1
          ecov$lag_M[i,s,a,r] <- as.integer(tmp[2])
          if(length(tmp) == 4) {
            #print(tmp)
            data$n_poly_Ecov_M[i,s,a,r] <- as.integer(tmp[4])
          }
        }
        #it is already set to 1
        #if(length(tmp)== 3 & tmp[4] == "linear") ecov_poly_R[i,s] <- 1 
      }
    } else {
      input$log$ecov <- c(input$log$ecov, "no ecov$M_how is provided so no effects on natural mortality will be estimated. \n")
    }
    #print(c(data$n_stocks, data$n_ages, data$n_regions, data$n_Ecov, max(data$n_poly_Ecov_M)))
    # stop()
    par$Ecov_beta_M <- array(0, dim = c(data$n_stocks, data$n_ages, data$n_regions, data$n_Ecov, max(data$n_poly_Ecov_M)))
    map$Ecov_beta_M <- array(NA, dim = dim(par$Ecov_beta_M))
    
    if(all(is.na(ecov$M_effect_map))){ #wasn't provided or some non-NA values were provided
      k <- 0 #same beta for all effects of this covariate
      for(i in 1:data$n_Ecov) {
        for(s in 1:data$n_stocks) for(a in 1:data$n_ages) for(r in 1:data$n_regions) {
          if(data$Ecov_how_M[i,s,a,r] != 0) {
            map$Ecov_beta_M[s,a,r,i,] <- k + 1:data$n_poly_Ecov_M[i,s,a,r]
            #ecov$M_effect_map[s,a,r,i] <- i #same beta for all M for each Ecov
          }
        }
        k <- max(map$Ecov_beta_M, 0, na.rm = TRUE)
      }
    } else { #user provided 
      ecov$M_effect_map[] <- as.integer(factor(ecov$M_effect_map)) #reset to 1:n_unique_effects_M
      n_eff_tot <- max(ecov$M_effect_map, na.rm = TRUE)
      for(i in 1:data$n_Ecov) {
        for(s in 1:data$n_stocks) for(a in 1:data$n_ages) for(r in 1:data$n_regions) {
          if(data$Ecov_how_M[i,s,a,r] != 0) {
            n_eff <- data$n_poly_Ecov_M[i,s,a,r]
            k <- ecov$M_effect_map[s,a,r,i]
            if(n_eff>1) {
              k <- c(k, neff_tot+1:(neff-1))
              n_eff_tot <- max(k)
            }
            map$Ecov_beta_M[s,a,r,i,1:neff] <- k
          }
        }
      }
    }
    # unique_effects_M <- unique(ecov$M_effect_map[which(!is.na(ecov$M_effect_map))]) #will be at least 1 value
    # n_unique_effects_M <- length(unique_effects_M)
    # for(s in 1:data$n_stocks) for(a in 1:data$n_ages) for(r in 1:data$n_regions) for(i in 1:data$n_Ecov) if(data$Ecov_how_M[i,s,a,r]>0){
    #   n <- data$n_poly_Ecov_M[i,s,a,r]
    #   if(n>1) { #need unique parameters for higher order polynomials
    #     for(eff in 2:n) map$Ecov_beta_M[s,a,r,i,eff] <- map$Ecov_beta_M[s,a,r,i,1] + n_unique_effects_M * (eff-1)
    #   }
    # }
    
    #make q
    ecov$lag_q <- matrix(0, data$n_Ecov, data$n_indices)
    data$n_poly_Ecov_q <- matrix(1, data$n_Ecov, data$n_indices)
    data$Ecov_how_q <- matrix(0, data$n_Ecov, data$n_indices)
    if(!is.null(ecov$q_how)) {
      for(ind in 1:data$n_indices) for(i in 1:data$n_Ecov) {
        tmp <- strsplit(ecov$q_how[i,ind], split = "-")[[1]]
        if(!(length(tmp) %in% c(1,3:4))) stop(paste0("Form of ecov$q_how[",i,",", ind, "] for ecov, ", i, ", index , ", ind,
                                                     " should be 'none', 'lag-l-poly-p' or 'lag-l-linear'"))
        if(tmp[1] != "none"){
          data$Ecov_how_q[i,ind] <- 1
          ecov$lag_q[i,ind] <- as.integer(tmp[2])
          if(length(tmp) == 4) data$n_poly_Ecov_q[i,ind] <- as.integer(tmp[4])
        }
      }
    } else {
      input$log$ecov <- c(input$log$ecov, "no ecov$q_how is provided so no effects on catchability will be estimated. \n")
    }
    
    par$Ecov_beta_q <- array(0, dim = c(data$n_indices, data$n_Ecov, max(data$n_poly_Ecov_q)))
    map$Ecov_beta_q <- array(NA, dim = dim(par$Ecov_beta_q))
    k <- 1
    for(ind in 1:data$n_indices) for(i in 1:data$n_Ecov) if(data$Ecov_how_q[i,ind]>0){
      n <- data$n_poly_Ecov_q[i,ind]
      map$Ecov_beta_q[ind,i,1:n] <- k + 1:n
      k <- k + n
    }
    
    #make mu
    ecov$lag_mu <- array(0, dim = c(data$n_Ecov, data$n_stocks, data$n_ages, data$n_seasons, data$n_regions, data$n_regions-1))
    
    data$n_poly_Ecov_mu <- array(1, dim = c(data$n_Ecov,data$n_stocks, data$n_ages, data$n_seasons, data$n_regions, data$n_regions-1))
    data$Ecov_how_mu <- array(0, dim = c(data$n_Ecov,data$n_stocks, data$n_ages, data$n_seasons, data$n_regions, data$n_regions-1))
    if(data$n_regions > 1 & !is.null(ecov$move_how)) {
      for(s in 1:data$n_stocks) for(a in 1:data$n_ages) for(t in 1:data$n_seasons) for(r in 1:data$n_regions) for(rr in 1:(data$n_regions-1)) for(i in 1:data$n_Ecov) {
        tmp <- strsplit(ecov$move_how[i,s,a,t,r,rr], split = "-")[[1]]
        if(!(length(tmp) %in% c(1,3:4))) stop(paste0("Form of ecov$move_how[" ,i, ",", s, ",", a, ",", t, ",", r, ",", rr, "] for ecov, ", i, 
                                                     ", stock , ", s, ", age , ", a, ", season , ", t, ", from region , ", r, " to region , ", rr, " should be 'none', 'lag-l-poly-p' or 'lag-l-linear'"))
        if(tmp[1] != "none"){
          data$Ecov_how_mu[i,s,a,t,r,rr] <- 1
          ecov$lag_mu[i,s,a,t,r,rr] <- as.integer(tmp[2])
          if(length(tmp) == 4) data$n_poly_Ecov_mu[i,s,a,t,r,rr] <- as.integer(tmp[4])
        }
        #it is already set to 1
        #if(length(tmp)== 3 & tmp[4] == "linear") ecov_poly_R[i,s] <- 1 
      }
    } else {
      if(data$n_regions> 1) input$log$ecov <- c(input$log$ecov, "There is more than 1 region, but no ecov$move_how is provided so no effects on movement parameters will be estimated.\n")
    }
    par$Ecov_beta_mu <- array(0, dim = c(data$n_stocks, data$n_ages, data$n_seasons, data$n_regions, data$n_regions-1, data$n_Ecov, 
                                         max(data$n_poly_Ecov_mu,0)))
    map$Ecov_beta_mu <- array(NA, dim = dim(par$Ecov_beta_mu))
    
    if(data$n_regions>1) {
      if(all(is.na(ecov$move_effect_map))){ #wasn't provided or some non-NA values were provided
        k <- 0 #same beta for all effects of this covariate
        for(i in 1:data$n_Ecov) {
          for(s in 1:data$n_stocks) for(a in 1:data$n_ages) for(t in 1:data$n_seasons) for(r in 1:data$n_regions) for(rr in 1:(data$n_regions-1)){
            if(data$Ecov_how_mu[i,s,a,t,r,rr] != 0) {
              #ecov$move_effect_map[s,a,t,r,rr,i] <- NA #no effect so map it so.
              map$Ecov_beta_mu[s,a,t,r,rr,i,] <- k + 1:data$n_poly_Ecov_mu[i,s,a,t,r,rr]
            }
          }
          k <- max(map$Ecov_beta_mu, 0, na.rm = TRUE)
        }
      } else { #user provided 
        ecov$move_effect_map[] <- as.integer(factor(ecov$move_effect_map))
        n_eff_tot <- max(ecov$move_effect_map, na.rm = TRUE)
        for(i in 1:data$n_Ecov) {
          for(s in 1:data$n_stocks) for(a in 1:data$n_ages) for(t in 1:data$n_seasons) for(r in 1:data$n_regions) for(rr in 1:(data$n_regions-1)){
            if(data$Ecov_how_mu[i,s,a,t,r,rr] != 0) {
              n_eff <- data$n_poly_Ecov_mu[i,s,a,t,r,rr]
              k <- ecov$move_effect_map[s,a,t,r,rr,i]
              if(n_eff>1) {
                k <- c(k, neff_tot+1:(neff-1))
                n_eff_tot <- max(k)
              }
              map$Ecov_beta_mu[s,a,t,r,rr,i,] <- k
            }
          }
        }
      }
    }
    
    #   map$Ecov_beta_mu[,,,,,,1] <- ecov$move_effect_map# k + 1:n
    #   unique_effects_mu <- unique(ecov$move_effect_map[which(!is.na(ecov$move_effect_map))]) #will be at least 1 value
    #   n_unique_effects_mu <- length(unique_effects_mu)
    #   for(s in 1:data$n_stocks) for(a in 1:data$n_ages) for(t in 1:data$n_seasons) for(r in 1:data$n_regions) for(rr in 1:(data$n_regions-1)) for(i in 1:data$n_Ecov) if(data$Ecov_how_M[i,s,a,r]>0){
    #     n <- data$n_poly_Ecov_mu[i,s,a,t,r,rr]
    #     if(n>1) { #need unique parameters for higher order polynomials
    #       for(eff in 2:n) map$Ecov_beta_mu[s,a,t,r,rr,i,eff] <- map$Ecov_beta_mu[s,a,t,r,rr,i,1] + n_unique_effects_mu * (eff-1)
    #     }
    #   }
    # }
    
    # # check that Ecov year vector doesn't have missing gaps
    # pad Ecov if it starts after model year1 - max(lag)
    max.lag = max(c(ecov$lag_R,ecov$lag_M,ecov$lag_mu,ecov$lag_q))
    #if(is.null(ecov$lag)) stop("ecov$lag needs to be provided for each ecov")
    #if(!is.list(ecov$lag)) ecov$lag = lapply(ecov$lag, function(x) rep(x,n_effects))
    # print(max.lag)
    if(input$years_Ecov[1] > input$years[1] - max.lag){
      input$log$ecov <- c(input$log$ecov, "one or more ecov does not start by model year 1 - max(lag). Padding ecov... \n")
      data$Ecov_obs <- rbind(matrix(0, nrow = input$years_Ecov[1]-(input$years[1]-max.lag), ncol = data$n_Ecov), data$Ecov_obs)
      par$Ecov_obs_logsigma <- rbind(matrix(par$Ecov_obs_logsigma[1,], nrow = input$years_Ecov[1]-(input$years[1]-max.lag), ncol = data$n_Ecov, byrow=T), par$Ecov_obs_logsigma)
      map$Ecov_obs_logsigma <- rbind(matrix(NA, nrow = input$years_Ecov[1]-(input$years[1]-max.lag), ncol = data$n_Ecov), map$Ecov_obs_logsigma)
      data$Ecov_use_obs <- rbind(matrix(0, nrow = input$years_Ecov[1]-(input$years[1]-max.lag), ncol = data$n_Ecov), data$Ecov_use_obs)
      input$years_Ecov <- c(seq(input$years[1] - max.lag, input$years_Ecov[1]-1), input$years_Ecov)
    }
    
    # pad Ecov if it ends before last model year
    if(end_Ecov < end_model){
      input$log$ecov <- c(input$log$ecov, "Ecov last year is before model last year. Padding Ecov... \n")
      data$Ecov_obs <- rbind(data$Ecov_obs, matrix(0, nrow = end_model-end_Ecov, ncol = data$n_Ecov))
      par$Ecov_obs_logsigma <- rbind(par$Ecov_obs_logsigma, matrix(par$Ecov_obs_logsigma[NROW(par$Ecov_obs_logsigma),], nrow = end_model-end_Ecov, ncol = data$n_Ecov, byrow=T))
      map$Ecov_obs_logsigma <- rbind(map$Ecov_obs_logsigma, matrix(NA, nrow = end_model-end_Ecov, ncol = data$n_Ecov))
      data$Ecov_use_obs <- rbind(data$Ecov_use_obs, matrix(0, nrow = end_model-end_Ecov, ncol = data$n_Ecov))
      input$years_Ecov <- c(input$years_Ecov, seq(end_Ecov+1, end_model))
      end_Ecov <- end_model
    }
    data$n_years_Ecov <- dim(data$Ecov_obs)[1] # num years Ecov to model (padded)
    data$years_use_Ecov <- 1:data$n_years_Ecov - 1
    
    data$Ecov_use_re <- rep(1, data$n_Ecov)
    for(i in 1:data$n_Ecov){
      if(all(data$Ecov_use_obs[,i]==0)) {
        data$Ecov_use_re[i] <- data$Ecov_model[i] <- 0
        input$log$ecov <- c(input$log$ecov, paste0("No observations for Ecov", i, " so no latent process will be estimated.\n"))
      }
    }
    
    
    #set up Ecov_re with padded dimensions
    #par$Ecov_re = matrix(rnorm(data$n_years_Ecov*data$n_Ecov), data$n_years_Ecov, data$n_Ecov)
    par$Ecov_re = matrix(0, data$n_years_Ecov, data$n_Ecov)
    map$Ecov_re <- matrix(1:length(par$Ecov_re), data$n_years_Ecov, data$n_Ecov, byrow=FALSE)
    for(i in 1:data$n_Ecov){
      #tmp.pars[,i] <- if(data$Ecov_model[i]==0) rep(NA,3) else tmp.pars[,i]
      map$Ecov_re[,i] <- if(data$Ecov_model[i]==0) rep(NA,data$n_years_Ecov) else map$Ecov_re[,i]
      if(data$Ecov_model[i]==1) map$Ecov_re[1,i] <- NA # if Ecov is a rw, first year of Ecov_re is not used bc Ecov_x[1] uses Ecov1 (fixed effect)
    }
    ind.notNA <- which(!is.na(map$Ecov_re))
    map$Ecov_re[ind.notNA] <- 1:length(ind.notNA)
    
    
    # get index of Ecov_x to use for Ecov_out (Ecovs can have diff lag)
    # print(end_model)
    # print(ecov$lag_R)
    #stop()
    data$ind_Ecov_out_start_R <- which(input$years_Ecov == input$years[1]) - ecov$lag_R - 1
    data$ind_Ecov_out_end_R <- which(input$years_Ecov==end_model)- ecov$lag_R - 1 # -1 is for cpp indexing
    data$ind_Ecov_out_start_M <- which(input$years_Ecov == input$years[1]) - ecov$lag_M - 1
    data$ind_Ecov_out_end_M <- which(input$years_Ecov==end_model)- ecov$lag_M - 1 # -1 is for cpp indexing
    data$ind_Ecov_out_start_q <- which(input$years_Ecov == input$years[1]) - ecov$lag_q - 1
    data$ind_Ecov_out_end_q <- which(input$years_Ecov==end_model)- ecov$lag_q - 1 # -1 is for cpp indexing
    data$ind_Ecov_out_start_mu <- which(input$years_Ecov == input$years[1]) - ecov$lag_mu - 1
    data$ind_Ecov_out_end_mu <- which(input$years_Ecov==end_model)- ecov$lag_mu - 1 # -1 is for cpp indexing
    
    input$log$ecov <- c(input$log$ecov, paste0("Please check that the environmental covariates have been loaded and interpreted correctly.

      Model years: ", input$years[1], " to ", end_model,"
      Ecov years: ", input$years_Ecov[1], " to ", end_Ecov,"

    "))
    
    #if(!identical(length(ecov$lag), length(ecov$label), data$n_Ecov)) stop("Length of ecov$lag and ecov$label not equal to # Ecov")
    for(i in 1:data$n_Ecov){
      years <- input$years_Ecov[as.logical(data$Ecov_use_obs[,i])]
      lastyr <- tail(years,1)
      
      # recruitment
      for(s in 1:data$n_stocks) if(data$Ecov_how_R[i,s] > 0){ 
        input$log$ecov <- c(input$log$ecov, paste0("Ecov ",i,": ",ecov$label[i],"
          ",c('*NO*','Controlling','Limiting','Lethal','Masking','Directive')[data$Ecov_how_R[i,s]+1]," (",
                                                   ifelse(data$n_poly_Ecov_R[i,s] == 1, "linear", paste0("polynomial order = ", data$n_poly_Ecov_R[i,s])),
                                                   ") effect on: recruitment for stock ", s, "

          Model years:
        "))
        
        input$log$ecov <- c(input$log$ecov, paste0(years, collapse=', '))
        
        input$log$ecov <- c(input$log$ecov, paste0("Lag: ",ecov$lag_R[i,s],"
        Ex: ",ecov$label[i]," in ",years[1]," affects recruitment in ",years[1+ecov$lag_R[i,s]],"
            ",ecov$label[i]," in ",lastyr," affects recruitment in ",lastyr+ecov$lag_R[i,s],"

        \n"))
      }
      
      #M
      for(s in 1:data$n_stocks) for(a in 1:data$n_ages) for(r in 1:data$n_regions) if(data$Ecov_how_M[i,s,a,r] == 1){
        input$log$ecov <- c(input$log$ecov, paste0("Ecov ",i,": ",ecov$label[i]," effect (", 
                                                   ifelse(data$n_poly_Ecov_M[i,s,a,r] == 1, "linear", paste0("polynomial order = ", data$n_poly_Ecov_M[i,s,a,r])), 
                                                   ") on: M for stock ", s, " at age ", a, " in region ", r,
                                                   "

          Model years:
        "))
        input$log$ecov <- c(input$log$ecov, paste0(years, collapse= ", "))
        
        input$log$ecov <- c(input$log$ecov, paste0("Lag: ",ecov$lag_M[i,s,a,r],"
        Ex: ",ecov$label[i]," in ",years[1]," affects M in ",years[1+ecov$lag_M[i,s,a,r]],"
            ",ecov$label[i]," in ",lastyr," affects M in ",lastyr+ecov$lag_M[i,s,a,r],"

        \n"))
      }
      
      # q
      for(j in 1:data$n_indices) if(data$Ecov_how_q[i,j] == 1){
        input$log$ecov <- c(input$log$ecov, paste0("Ecov ",i,": ",ecov$label[i]," effect (", 
                                                   ifelse(data$n_poly_Ecov_q[i,j] == 1, "linear", paste0("polynomial order = ", data$n_poly_Ecov_q[i,j])), 
                                                   ") on: q for index ", j, 
                                                   " 

          Model years:
        "))
        input$log$ecov <- c(input$log$ecov, paste0(years, collapse = ", "))
        
        input$log$ecov <- c(input$log$ecov, paste0("Lag: ",ecov$lag_q[i,j],"
        Ex: ",ecov$label[i]," in ",years[1]," affects index ", j, " in ",years[1+ecov$lag_q[i,j]],"
            ",ecov$label[i]," in ",lastyr," affects index ", j, " in ",lastyr+ecov$lag_q[i,j],"

        \n"))
      }
      
      #movement
      if(data$n_regions>1) for(s in 1:data$n_stocks) for(a in 1:data$n_ages) for(t in 1:data$n_seasons) if(data$n_regions>1) {
        for(r in 1:data$n_regions) for(rr in 1:(data$n_regions-1)) if(data$Ecov_how_mu[i,s,a,t,r,rr] == 1){
          input$log$ecov <- c(input$log$ecov, paste0("Ecov ",i,": ",ecov$label[i]," effect (", 
                                                     ifelse(data$n_poly_Ecov_mu[i,s,a,t,r,rr] == 1, "linear", paste0("polynomial order = ", data$n_poly_Ecov_mu[i,s,a,t,r,rr])), 
                                                     ") on: movement for stock ", s, " at age ", a, "in season ", t, "from region ", r, "to ", rr, "of the other regions
            
            Model years:
          "))
          input$log$ecov <- c(input$log$ecov, paste0(years, collapse=", "))
          
          input$log$ecov <- c(input$log$ecov, paste0("Lag: ",ecov$lag_mu[i,s,a,t,r,rr],"
          Ex: ",ecov$label[i]," in ",years[1]," affects movement in ",years[1+ecov$lag_mu[i,s,a,t,r,rr]],"
              ",ecov$label[i]," in ",lastyr," affects movement in ",lastyr+lag_mu[i,s,a,t,r,rr],"

          \n"))
        }
      }
    }
    #input$Ecov_names <- list(input$Ecov_names)
    
    # Ecov process pars
    par$Ecov_process_pars = matrix(0, 3, data$n_Ecov) # nrows = RW: 2 par (Ecov1, log_sig), AR1: 3 par (mu, log_sig, phi); ncol = N_ecov
    #this row is for the mean not the sd of the process
    par$Ecov_process_pars[1,] = -1.3 # start sig_ecov at 0.27
    #changing the initial value for sig_ecov to the right place actually causes tests to not pass!
    #par$Ecov_process_pars[2,] = -1.3 # start sig_ecov at 0.27
    if(!is.null(ecov$process_mean_vals)){
      for(i in 1:data$n_Ecov) if(data$Ecov_model[i]==2) par$Ecov_process_pars[1,i] = ecov$process_mean_vals[i]
    }
    if(!is.null(ecov$process_sig_vals)){
      for(i in 1:data$n_Ecov) par$Ecov_process_pars[2,i] = log(ecov$process_sig_vals[i])
    }
    if(!is.null(ecov$process_cor_vals)){
      inv_trans_rho <- function(rho) log(rho+1) - log(1-rho) 
      for(i in 1:data$n_Ecov)  if(data$Ecov_model[i]==2) par$Ecov_process_pars[3,i] = inv_trans_rho(ecov$process_cor_vals[i])
    }
    
    # turn off Ecov pars if no Ecov (re, process)
    # turn off 3rd Ecov par if it's a RW
    map$Ecov_process_pars <- par$Ecov_process_pars
    for(i in 1:data$n_Ecov) {
      map$Ecov_process_pars[,i] <- if(data$Ecov_model[i]==0) rep(NA,3) else map$Ecov_process_pars[,i]
      map$Ecov_process_pars[3,i] <- ifelse(data$Ecov_model[i]==1, NA, 0)
    }
    ind.notNA <- which(!is.na(map$Ecov_process_pars))
    map$Ecov_process_pars[ind.notNA] <- 1:length(ind.notNA)
    
    #fill in initial values for any covariate effects
    for(i in c("beta_q", "beta_R", "beta_M", "beta_mu")){
      inp <- paste0(i, "_vals")
      if(!is.null(ecov[[inp]])){
        if(!length(dim(ecov[[inp]])) == length(dim(par[[paste0("Ecov_",i)]])) | !all(dim(ecov[[inp]]) == dim(par[[paste0("Ecov_",i)]]))){
          stop(paste0("ecov$",i, "_vals is not an array of the correct dimensions."))
        }
        par[[paste0("Ecov_",i)]] <- ecov[[inp]]
      }
    }
    
  } # end load Ecov
  
  
  map$Ecov_obs_logsigma <- factor(map$Ecov_obs_logsigma)
  map$Ecov_obs_sigma_par <- factor(map$Ecov_obs_sigma_par)
  map$Ecov_obs_logsigma_re <- factor(map$Ecov_obs_logsigma_re)
  map$Ecov_process_pars = factor(map$Ecov_process_pars)
  map$Ecov_re <- factor(map$Ecov_re)
  map$Ecov_beta_R <- factor(map$Ecov_beta_R)
  map$Ecov_beta_M <- factor(map$Ecov_beta_M)
  map$Ecov_beta_q <- factor(map$Ecov_beta_q)
  map$Ecov_beta_mu <- factor(map$Ecov_beta_mu)
  
  input$data = data
  input$par = par
  input$map = map
  if(length(input$log$ecov)) input$log$ecov <- c("Ecov: \n", input$log$ecov)
  # add vector of all observations for one step ahead residuals ==========================
  if(!is.null(input$data$obsvec)) input <- set_osa_obs(input)
  #print("osa_obs")
  
  # projection data will always be modified by 'prepare_projection'
  input = set_proj(input, proj.opts = NULL) #proj options are used later after model fit, right?
  #print("proj")
  
  #set any parameters as random effects
  input$random = NULL
  input = set_random(input)
  #print("random")
  input$options$ecov <- ecov
  
  return(input)
}

#' Specify configuration for fully-selected fishing mortality
#'
#' @param input list containing data, parameters, map, and random elements (output from \code{\link{wham::prepare_wham_input}})
#' @param F_opts (optional) named list of initial values for annual fully-selected fishing mortality and configuration method for estimation.
#' 
#' \code{F_opts} specifies a few as well as the effect on the population. Environmental covariate data need not span
#' the same years as the fisheries data. It can be \code{NULL} if no environmental data are to be fit.
#' Otherwise, it must be a named list with the following components:
#'   \describe{
#'     \item{$F}{matrix (n_years x n_fleets) of (initial) values for fully-selected fishing morality.}
#'     \item{$F_config}{integer 1: (default) configure F parameters (on log scale) as an F in the initial year and then deviations from one year to the next,
#'        or 2: configure F parameters as (log) annual values.}
#'  }
#'
#' @export
set_F = function(input, F_opts = NULL)
{
  asap3 = input$asap3
  input$data$F_config <- 1 #log_F1, F_devs
  input$par$F_pars <- matrix(0,input$data$n_years_model, input$data$n_fleets)
  #F_devs = matrix(0, input$data$n_years_model-1, input$data$n_fleets)
  if(!is.null(asap3)) {
    k = 1
    for(i in 1:length(asap3))for(j in 1:length(asap3[[i]]$F1_ini)){
      input$par$F_pars[1,k] = log(asap3[[i]]$F1_ini[j]) # use F1_ini values from asap3 file  
      k = k + 1
    }
  } else {
    input$par$F_pars[1,] = log(0.2) # old
  }
  if(!is.null(F_opts$F_config)) input$data$F_config <- F_opts$F_config
  if(!is.null(F_opts$F)) {
    if(input$data$F_config == 2) {
      input$par$F_pars[] = log(F_opts$F[])
    } else { # F_config = 1
      input$par$F_pars[1,] <- log(F_opts$F[1,])
      for(f in 1:input$data$n_fleets) input$par$F_pars[-1,f] <- diff(log(F_opts$F[,f]))
    }
  }
  input$options$F <- F_opts
  return(input)
}

#' Make one or more selectivity blocks with age-specific parameters
#'
#' @param input list containing data and parameters (output from \code{\link{prepare_wham_input}})
#' @param selblocks numeric, number of age-specific selectivity blocks
#'
#' @return a modified list of data and parameters
#'
#' @examples
#' \dontrun{
#' asap3 = read_asap3_dat("ASAP_SNEMAYT.dat")
#' input = prepare_wham_input(asap3)
#' input = set_age_sel0(input, selblocks = 1:3)
#' mod = fit_wham(input)
#' }
#'
#' @export
set_age_sel0 <- function(input, selblocks){
  temp = input
  temp$map$logit_selpars = as.integer(as.character(temp$map$logit_selpars))
  temp$map$logit_selpars = matrix(temp$map$logit_selpars, temp$data$n_selblocks, temp$data$n_ages + 6)
  temp$map$logit_selpars[selblocks,] = NA
  n.par  = sum(!is.na(temp$map$logit_selpars))
  if(n.par) temp$map$logit_selpars[which(!is.na(temp$map$logit_selpars))] = 1:sum(!is.na(temp$map$logit_selpars))
  temp$data$selblock_models[selblocks] = 1
  for(x in selblocks){
    ind = list(fleets = which(apply(temp$data$selblock_pointer_fleets == x,2,sum) > 0))
    ind$indices = which(apply(temp$data$selblock_pointer_indices == x,2,sum) > 0)
    paa = matrix(nrow = 0, ncol = temp$data$n_ages)
    if(length(ind$fleets)) for(f in ind$fleets)
    {
      y = temp$data$catch_paa[f,which(temp$data$selblock_pointer_fleets[,f] == x & temp$data$use_catch_paa[,f] == 1),]
      paa = rbind(paa,y)
    }
    if(length(ind$indices)) for(i in ind$indices)
    {
      y = temp$data$index_paa[i,which(temp$data$selblock_pointer_indices[,i] == x & temp$data$use_index_paa[,i] == 1),]
      paa = rbind(paa,y)
    }
    y = apply(paa,2,sum)
    temp$par$logit_selpars[x,temp$data$n_ages + 1:6] = Inf
    temp$par$logit_selpars[x,which(y < 1e-5)] = -Inf
    temp$par$logit_selpars[x,which(y >= 1e-5)] = 0
    if(sum(y >= 1e-5)) temp$map$logit_selpars[x,which(y >= 1e-5)] = n.par + 1:sum(y >= 1e-5)
    n.par = n.par + sum(y >= 1e-5)
  }
  temp$map$logit_selpars = factor(temp$map$logit_selpars)
  return(temp)
}

#'  Specify the age composition models for fleet(s) and indices.
#'
#' @param input list containing data, parameters, map, and random elements (output from \code{\link{wham::prepare_wham_input}})
#' @param age_comp specifies the age composition models for fleet(s) and indices. If \code{NULL}, the multinomial is used because this was the only option in ASAP. 
#'
#' The age composition models available are:
#'   \describe{
#'     \item{\code{"multinomial"}}{Multinomial. This is the default because it was the only option in ASAP. 0 parameters.}
#'     \item{\code{"dir-mult"}}{Saturating Dirichlet-multinomial, parameterized such that effective-sample-size is a nonlinear and saturating function with respect to input-sample-size. 1 parameter. Effective sample size is estimated by the model (\href{https://www.ccamlr.org/es/system/files/science_journal_papers/07candy.pdf}{Candy 2008})}
#'     \item{\code{"dirichlet-pool0"}}{Dirichlet, pooling zero observations with adjacent age classes. 1. parameter. See \href{https://www.sciencedirect.com/science/article/abs/pii/S0165783613003093}{Francis 2014} and \href{https://cdnsciencepub.com/doi/abs/10.1139/cjfas-2015-0532}{Albertsen et al. 2016}}
#'     \item{\code{"dirichlet-miss0"}}{}{Dirichlet, treating zero observations as missing. 1 parameter.}
#'     \item{\code{"logistic-normal-miss0"}}{Logistic normal, treating zero observations as missing. 1 parameter.}
#'     \item{\code{"logistic-normal-ar1-miss0"}}{Logistic normal, treating zero observations as missing. 1 parameter.}
#'     \item{\code{"logistic-normal-pool0"}}{Logistic normal, pooling zero observations with adjacent age classes. 1 parameter. See \href{https://doi.org/10.1093/icesjms/fsl024}{Schnute and Haigh (2007)} and \href{https://doi.org/10.1016/j.fishres.2013.12.015}{Francis (2014)}}.
#'     \item{\code{"logistic-normal-01-infl"}}{Zero-or-one inflated logistic normal. Inspired by zero-one inflated beta in \href{https://www.sciencedirect.com/science/article/abs/pii/S0167947311003628}{Ospina and Ferrari (2012)}. 3 parameters. . No OSA residuals.}
#'     \item{\code{"logistic-normal-01-infl-2par"}}{Zero-one inflated logistic normal where p0 is a function of binomial sample size. 2 parameters. No OSA residuals.}
#'     \item{\code{"mvtweedie"}}{Multivariate-tweedie, where the product of composition proportions and input sample sizes follows a distribution with mean equal to the product of predicted proportions and input sample size, and other parameters define the ratio of effective to input sample size (with is bounded 0 to Inf) and the probability of zeros. 2 parameters. No OSA residuals.}
#'     \item{\code{"dir-mult-linear"}}{Linear Dirichlet-multinomial, parameterized such that effective-sample-size is a linear function with respect to input-sample-size, estimating 1 parameter, \eqn{log(\theta)}, where the ratio of effective and input sample size is approximately \eqn{\theta / (1+\theta)}, i.e., the logistic transformation of the estimated parameter \eqn{log(\theta)}.  (\href{https://doi.org/10.1016/j.fishres.2016.06.005}{Thorson et al. 2017}) }
#'   }
#' The two Dirichlet-multinomial options will only differ when input-sample-size differs among years.  In these cases, the linear-Dirichlet multinomial is designed to decrease the effective sample size in each year by approximately the same proportion, while the saturating-Dirichlet multinomial will decrease the years with highest input-sample-size much more than those with lower input-sample-size.
#' One-step-ahead residuals will be calculated for all but options 8-10 when \code{do.osa=TRUE} (Nielsen et al. in prep.). An age composition model needs
#' to be specified for each fleet and index. If you would like all fleets and indices to use the same age composition likelihood, you 
#' can simply specify one of the strings above, i.e. \code{age_comp = "logistic-normal-miss0"}. If you do not want the same
#' age composition model to be used for all fleets and indices, you must specify a named list with the following entries:
#'   \describe{
#'     \item{$fleets}{A vector of the above strings with length = the number of fleets.}
#'     \item{$indices}{A vector of the above strings with length = the number of indices.}
#'   }
#'
#' @return a named list with same elements as the input provided with age composition likelihood options modified.
#'
#' @seealso \code{\link{prepare_wham_input}} 
#'
#' @examples
#' \dontrun{
#' wham.dir <- find.package("wham")
#' path_to_examples <- system.file("extdata", package="wham")
#' asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))
#' input <- prepare_wham_input(asap3)
#' input <- set_age_comp(input, age_comp = "logistic-normal-miss0") #no longer multinomial
#' }
#'
#' @export
#'
set_age_comp = function(input, age_comp)
{
  data = input$data
  par = input$par
  map = input$map
  all_models <- c( "multinomial",
                   "dir-mult",
                   "dirichlet-miss0",
                   "dirichlet-pool0",
                   "logistic-normal-miss0",
                   "logistic-normal-ar1-miss0",
                   "logistic-normal-pool0",
                   "logistic-normal-01-infl",
                   "logistic-normal-01-infl-2par",
                   "mvtweedie",
                   "dir-mult-linear" )
  input$log$age_comp <- list()
  n_pars <- c(0,1,1,1,1,2,1,3,2,2,1)
  if(is.null(age_comp)){
    data$age_comp_model_fleets = rep(1, data$n_fleets) # multinomial by default
    data$age_comp_model_indices = rep(1, data$n_indices) # multinomial by default
  } else {
    if(is.character(age_comp) & length(age_comp)==1){ # all use the same model
      name_change = age_comp == "dirichlet"
      if(any(name_change)){
        input$log$age_comp <- c(input$log$age_comp, "'dirichlet' is no longer an option and the old option is equivalent to 'dirichlet-pool0' so using that.\n")
        age_comp[which(name_change)] = "dirichlet-pool0"
      }
      themod <- match(age_comp, all_models)
      if(is.na(themod)) stop("age_comp option not recognized. See ?prepare_wham_input.")
      data$age_comp_model_fleets = rep(themod, data$n_fleets)
      data$age_comp_model_indices = rep(themod, data$n_indices)
    } else {
      if(all(names(age_comp) %in% c("fleets","indices"))){
        name_change = list(fleets = age_comp$fleets == "dirichlet", indices = age_comp$indices == "dirichlet")
        if(any(unlist(name_change))){
          input$log$age_comp <- c(input$log$age_comp, "'dirichlet' is no longer an option and the old option is equivalent to 'dirichlet-pool0' so using that.\n")
          age_comp$fleets[which(name_change$fleets)] = "dirichlet-pool0"
          age_comp$indices[which(name_change$indices)] = "dirichlet-pool0"
        }
        themods <- match(age_comp$fleets, all_models)
        if(any(is.na(themods))) stop("age_comp$fleets option not recognized. See ?prepare_wham_input for available options.")
        if(length(themods) != data$n_fleets) stop("age_comp$fleets must have length = the number of fleets")
        data$age_comp_model_fleets = themods
        
        themods <- match(age_comp$indices, all_models)
        if(any(is.na(themods))) stop("age_comp$indices option not recognized. See ?prepare_wham_input for available options.")
        if(length(themods) != data$n_indices) stop("age_comp$indices must have length = the number of indices")
        data$age_comp_model_indices = themods
      } else {
        stop("age_comp must either be a character or a named ('fleets','indices') list. See ?prepare_wham_input.")
      }
    }
  }
  
  # age comp pars
  par$catch_paa_pars = matrix(0,data$n_fleets, 3) 
  par$index_paa_pars = matrix(0,data$n_indices, 3) 
  neff <- data$catch_Neff
  neff[neff <= 0] <- 1
  catch_neff <- apply(neff,2,mean, na.rm=TRUE)
  ind = which(data$age_comp_model_fleets %in% 5:7)
  par$catch_paa_pars[ind,1] <- 0.5*log(catch_neff[ind])
  neff <- data$index_Neff
  neff[neff <= 0] <- 1
  index_neff <- apply(neff,2,mean, na.rm=TRUE)#[which(apply(data$use_index_paa,2,sum)>0)]
  ind = which(data$age_comp_model_indices %in% 5:7)
  par$index_paa_pars[ind,1] <- 0.5*log(index_neff[ind])
  
  map$index_paa_pars = matrix(NA,data$n_indices, 3)
  for(i in 1:data$n_indices) if(sum(data$use_index_paa[,i])){
    if(data$age_comp_model_indices[i] %in% c(2:5,7,11))  map$index_paa_pars[i,1] = 1
    if(data$age_comp_model_indices[i] %in% c(6,9,10))  map$index_paa_pars[i,1:2] = 1
    if(data$age_comp_model_indices[i] %in% 8)  map$index_paa_pars[i,1:3] = 1
  }
  map$catch_paa_pars = matrix(NA,data$n_fleets, 3)
  for(i in 1:data$n_fleets) if(sum(data$use_catch_paa[,i])){
    if(data$age_comp_model_fleets[i] %in% c(2:5,7,11))  map$catch_paa_pars[i,1] = 1
    if(data$age_comp_model_fleets[i] %in% c(6,9,10))  map$catch_paa_pars[i,1:2] = 1
    if(data$age_comp_model_fleets[i] %in% 8)  map$catch_paa_pars[i,1:3] = 1
  }
  nest = sum(map$index_paa_pars,na.rm=TRUE)
  if(nest) map$index_paa_pars[which(!is.na(map$index_paa_pars))] = 1:nest
  map$index_paa_pars = factor(map$index_paa_pars)
  nest = sum(map$catch_paa_pars,na.rm=TRUE)
  if(nest) map$catch_paa_pars[which(!is.na(map$catch_paa_pars))] = 1:nest
  map$catch_paa_pars = factor(map$catch_paa_pars)
  
  input$data = data
  input$par = par
  input$map = map
  if(length(input$log$age_comp))	input$log$age_comp <- c("Age composition: \n", input$log$age_comp)
  input$options$age_comp <- age_comp
  
  if(!is.null(input$data$obsvec)) input <- set_osa_obs(input)
  
  return(input)
}

#' Specify index selectivity blocks and aggregate and age composition observations for indices
#'
#' @param input list containing data, parameters, map, and random elements (output from \code{\link{wham::prepare_wham_input}})
#' @param index_info (optional) list specifying various aspects about catch by indices (see details)
#' 
#' \code{index_info} specifies observations, and various configuration options for index-specific catch observations and will overwrite attributes specified in the ASAP data file.
#' If \code{NULL}, all settings from the ASAP data file or basic_info are used.
#' \code{index_info} is a list with any of the following entries:
#'   \describe{
#'     \item{$n_indices}{number of indices}
#'     \item{$index_regions}{vector (n_indices) of regions where each fleet operates.}
#'     \item{$index_seasons}{vector (n_indices) of 0/1 values flagging which seasons each index occurs.}
#'     \item{$agg_indices}{matrix (n_years_model x n_indices) of annual aggregate index catches.}
#'     \item{$agg_index_cv}{matrix (n_years_model x n_indices) of CVs for annual aggregate index catches.}
#'     \item{$fracyy_indices}{matrix (n_years_model x n_indices) of fractions of year at which index occurs within the season (difference between time of survey and time at start of season).}
#'     \item{$use_indices}{matrix (n_years_model x n_indices) of 0/1 values flagging whether to use aggregate observations.}
#'     \item{$units_indices}{matrix (n_years_model x n_indices) of 1/2 values flagging whether aggregate observations are biomass (1) or numbers (2).}
#'     \item{$index_paa}{array (n_indices x n_years_model x n_ages) of annual catch proportions at age by index.}
#'     \item{$use_index_paa}{matrix (n_years_model x n_indices) of 0/1 values flagging whether to use proportions at age observations.}
#'     \item{$units_index_paa}{matrix (n_years_model x n_indices) of 1/2 values flagging whether composition observations are biomass (1) or numbers (2).}
#'     \item{$index_Neff}{matrix (n_years_model x n_indices) of effective sample sizes for proportions at age observations.}
#'     \item{$selblock_pointer_indices}{matrix (n_years_model x n_indices) of itegers indicated selblocks to use.}
#'   }
#'
#' @export
set_indices = function(input, index_info=NULL) {
  data = input$data
  asap3 = input$asap3
  if(is.null(asap3)) {
    data$n_indices = 1
  } else {
    input$index_names <- NULL
    for(i in 1:length(asap3)) {
      which_indices <- which(asap3[[i]]$use_index ==1)
      asap3[[i]]$n_indices = length(which_indices)
      asap3[[i]]$survey_index_units <- asap3[[i]]$index_units[which_indices]
      asap3[[i]]$survey_acomp_units <- asap3[[i]]$index_acomp_units[which_indices]
      asap3[[i]]$survey_WAA_pointers <- asap3[[i]]$index_WAA_pointers[which_indices]
      asap3[[i]]$survey_month <- asap3[[i]]$index_month[which_indices]
      #asap3[[i]]$survey_month <- matrix(asap3[[i]]$index_month[which_indices], asap3[[i]]$n_years, asap3[[i]]$n_indices, byrow = TRUE)
      asap3[[i]]$use_survey_acomp <- asap3[[i]]$use_index_acomp[which_indices]
      asap3[[i]]$index_WAA_pointers = asap3[[i]]$index_WAA_pointers[which_indices]
      asap3[[i]]$IAA_mats <- asap3[[i]]$IAA_mats[which_indices]
      asap3[[i]]$use_survey <- asap3[[i]]$use_index[which_indices]
      if(!is.null(asap3[[i]]$index.names)) input$index_names <- c(input$index_names, asap3[[i]]$index.names[which_indices])
    }
    n_indices_per_region = sapply(asap3, function(x) x$n_indices)
    data$n_indices = sum(n_indices_per_region)
  } 
  if(!is.null(index_info$n_indices)) data$n_indices = index_info$n_indices
  if(is.null(input$index_names)) input$index_names <- paste0("index_", 1:data$n_indices)
  
  data$index_regions = rep(1, data$n_indices)
  data$index_seasons = rep(1, data$n_indices)
  
  data$agg_indices = data$agg_index_sigma = data$index_Neff = matrix(NA, data$n_years_model, data$n_indices)
  data$index_paa = array(NA, dim = c(data$n_indices, data$n_years_model, data$n_ages))
  data$use_indices = matrix(1, data$n_years_model, data$n_indices)
  data$use_index_paa = matrix(1, data$n_years_model, data$n_indices)
  data$selblock_pointer_indices = matrix(0, data$n_years_model, data$n_indices)
  data$units_indices = rep(2,data$n_indices)
  data$units_index_paa = rep(2,data$n_indices)
  data$fracyr_indices = matrix(data$fracyr_seasons[1]*0.5, data$n_years_model, data$n_indices)
  
  if(!is.null(asap3)) {
    k <- 1
    for(i in 1:length(asap3)) {
      for(j in 1:asap3[[i]]$n_indices) {
        data$index_regions[k] = i #each asap file is a separate region
        data$units_indices[k] <- asap3[[i]]$survey_index_units[j]
        tmp = (asap3[[i]]$survey_month[j]-1)/12 #make sure that this is right
        int_starts <- cumsum(c(0,data$fracyr_seasons))
        ind = max(which(int_starts <= tmp))
        data$index_seasons[k] = ind
        data$fracyr_indices[,k] = tmp - int_starts[ind]
        
        data$agg_indices[,k] = asap3[[i]]$IAA_mats[[j]][,2]
        for(y in 1:data$n_years_model) if(asap3[[i]]$IAA_mats[[j]][y,2] < 1e-15) data$use_indices[y,k] = 0
        data$agg_index_sigma[,k] = asap3[[i]]$IAA_mats[[j]][,3]
        
        temp = asap3[[i]]$IAA_mats[[j]][,3 + 1:data$n_ages]
        temp[which(is.na(temp))] = 0
        temp[which(temp<0)] = 0
        data$index_paa[k,,] = temp/apply(temp,1,sum) #all 0s will make NaN
        if(asap3[[i]]$use_survey_acomp[j] != 1){
          data$use_index_paa[,k] = 0
        } else {
          for(y in 1:data$n_years_model) {
            flag <- asap3[[i]]$IAA_mats[[j]][y,4 + data$n_ages] < 1e-15 | sum(data$index_paa[k,y,] > 1e-15) < 2 | any(is.na(data$index_paa[k,y,]))
            if(flag) data$use_index_paa[y,k] = 0
          }
        }
        data$units_index_paa[k] <- asap3[[i]]$survey_acomp_units[j]
        data$index_Neff[,k] = asap3[[i]]$IAA_mats[[j]][,4 + data$n_ages]
        data$selblock_pointer_indices[,k] = max(data$selblock_pointer_fleets) + k #set_catch already called
        k <- k + 1
      }
    }
  }
  else {
    data$agg_indices[] = 10
    data$agg_index_sigma[] = 0.3
    data$index_paa[] = 1/data$n_ages
    data$index_Neff[] = 100
    data$selblock_pointer_indices[] = rep(1:data$n_indices, each = data$n_years_model) + max(data$selblock_pointer_fleets)
    input$index_names <- paste0("Index ", 1:data$n_indices)
  }
  
  if(!is.null(index_info$use_indices)) data$use_indices[] = index_info$use_indices
  if(!is.null(index_info$use_index_paa)) data$use_index_paa[] = index_info$use_index_paa
  if(!is.null(index_info$units_indices)) data$units_indices[] = index_info$units_indices
  if(!is.null(index_info$fracyr_indices)) data$fracyr_indices[] = index_info$fracyr_indices
  if(!is.null(index_info$agg_indices)) data$agg_indices[] = index_info$agg_indices
  if(!is.null(index_info$index_cv)) data$agg_index_sigma[] = index_info$index_cv
  if(!is.null(index_info$index_paa)) data$index_paa[] = index_info$index_paa
  if(!is.null(index_info$units_index_paa)) data$units_index_paa[] = index_info$units_index_paa
  if(!is.null(index_info$index_Neff)) data$index_Neff[] = index_info$index_Neff
  if(!is.null(index_info$selblock_pointer_indices)) data$selblock_pointer_indices[] = index_info$selblock_pointer_indices
  if(!is.null(index_info$index_seasons)) data$index_seasons[] = index_info$index_seasons
  if(!is.null(index_info$index_regions)) data$index_regions[] = index_info$index_regions
  
  ################################################################################
  # for plotting, in years where index is not used set sigma = avg of used years
  tmp <- data$agg_index_sigma
  tmp[data$use_indices == 0] = NA
  mean_agg_ind_sigma <- apply(tmp, 2, mean, na.rm=T)
  for(i in 1:data$n_indices) data$agg_index_sigma[data$use_indices[,i] == 0,i] = mean_agg_ind_sigma[i]
  ################################################################################
  
  data$index_paa[is.na(data$index_paa)] = 0
  data$agg_index_sigma[which(data$agg_index_sigma < 1e-15)] = 100  
  data$agg_index_sigma = sqrt(log(data$agg_index_sigma^2 + 1))
  
  input$par$log_index_sig_scale = rep(0, data$n_indices)
  input$map$log_index_sig_scale = factor(rep(NA, data$n_indices))
  input$asap3 <- asap3
  
  input$data = data
  input$options$index <- index_info
  if(!is.null(input$par$logit_selpars)) input <- set_selectivity(input, input$options$selectivity)
  if(!is.null(input$data$obsvec)) input <- set_osa_obs(input)
  return(input)
}

#' Specify model and parameter configuration for "extra" mortality not directly attributed to natural mortality
#'
#' @param input list containing data, parameters, map, and random elements (output from \code{\link{wham::prepare_wham_input}})
#' @param L (optional) list specifying "extra" mortality options: model, random effects, initial values, and parameters to fix (see details)
#' 
#' \code{L} specifies estimation options for "extra" mortality.
#' If \code{NULL}, This mortality source is not used. 
#' \code{L} is a list with the following entries:
#'   \describe{
#'     \item{$model}{length = n_regions. "extra" mortality model options are:
#'                    \describe{
#'                      \item{"none"}{(default) no extra mortality for this region.}
#'                      \item{"constant"}{estimate a single mean mortality for the region shared across all ages}
#'                      \item{"iid_re"}{estimate independent random effects over years, for the region}
#'                      \item{"ar1_re"}{estimate random effect correlated over years, for the region}
#'                    }
#'                  }
#'     \item{$initial_means}{Initial/mean L-at-region}
#'     \item{$sigma_vals}{Initial standard deviation by region value to use for the L random effects. Values are not used if \code{L$model} = "none".}
#'     \item{$cor_vals}{Initial correlation values to use for the L random effects. If unspecified all initial values are 0}
#'   }
#'
#' @export
set_L = function(input, L)
{
  data = input$data
  par = input$par
  map = input$map
  #asap3 = input$asap3
  
  #clear any map definitions that may exist. necessary because some configurations may not define map elements.
  map <- map[(!names(map) %in% c("L_repars", "L_re"))]
  
  #L$model length is n_regions
  data$L_model = rep(0, data$n_regions)
  par$L_re = matrix(0, data$n_years_model, data$n_regions)
  par$L_repars = matrix(0,data$n_regions, 3) #mean, sig, rho
  map$L_re = matrix(NA, data$n_years_model, data$n_regions)
  map$L_repars = matrix(NA,data$n_regions, 3) #mean, sig, rho
  
  # natural mortality options, default = use values from ASAP file, no estimation
  if(!is.null(L)){
    if(!is.null(L$model)){ # L model options
      L_mods = c("none","constant","iid_re","ar1_re")
      if(!(L$model %in% L_mods)) stop(paste0("L$model must be one of these: ", paste0(L_mods, collapse=",")))
      data$L_model[] = match(L$model, L_mods) - 1
    }
  }
  inv_trans_rho <- function(rho, s = 1) (log(rho+1) - log(1-rho))/s
  k = 1
  for(r in 1:data$n_regions) {
    if(data$L_model[r] >0) {
      map$L_repars[r,1] <- k
      k <- k + 1
      if(!is.null(L$initial_means)){
        par$L_repars[r,1] = log(L$initial_means[r])
      }
    }
    if(data$L_model[r] > 1){
      map$L_repars[r,2] <- k
      k <- k + 1
      map$L_re[,r] <- 1
      if(!is.null(L$sigma_vals)){
        par$L_repars[r,2] = log(L$sigma_vals[r])
      }
    }
    if(data$L_model[r] > 2){
      map$L_repars[r,3] <- k
      k <- k + 1
      if(!is.null(L$cor_vals)){
        par$L_repars[r,3] = inv_trans_rho(L$cor_vals[r])
      }
    }
  }
  map$L_re[which(map$L_re==1)] <- 1:sum(map$L_re==1, na.rm = TRUE)
  map$L_re = factor(map$L_re)
  map$L_repars = factor(map$L_repars)
  
  input$data = data
  input$par = par
  input$map = map
  
  #may need to update these 
  # projection data will always be modified by 'prepare_projection'
  input = set_proj(input, proj.opts = NULL) #proj options are used later after model fit, right?
  
  #set any parameters as random effects
  input$random = NULL
  input = set_random(input)
  input$options$L <- L
  return(input)
  
}

#' Specify model and parameter configuration for natural mortality
#'
#' @param input list containing data, parameters, map, and random elements (output from \code{\link{wham::prepare_wham_input}})
#' @param M (optional) list specifying natural mortality options: model, random effects, initial values, and parameters to fix (see details)
#' 
#' \code{M} specifies estimation options for natural mortality and can overwrite M-at-age values specified in the ASAP data file.
#' If \code{NULL}, the M-at-age matrix from the ASAP data file is used (M fixed, not estimated). To estimate M-at-age
#' shared/mirrored among some but not all ages, modify \code{M$means_map} (see vignette for more details). \code{M} is a list 
#' with the following entries:
#'   \describe{
#'     \item{$mean_model}{Character describing the type of model for M stock and regional models for natural mortality. Options are:
#'       \describe{
#'         \item{"fixed-M"}{Use initial values from ASAP3 dat files or \code{$initial_means} for (mean) M as fixed values. If no ASAP3 files
#'           and \code{$initial_means} is not provided, default is M = 0.2 for all stocks, regions and ages}
#'         \item{"estimate-M"}{estimate one or more (mean) M parameters. Default is to estimate a single M shared across all stocks and ages, but
#'           use \code{$means_map} to fix or estimate parameters for specific stocks, regions, ages.}
#'         \item{"weight-at-age"}{specifies M as a function of weight-at-age, \eqn{M_y,a = exp(b0 + b1*log(W_y,a))}, as in
#'           \href{https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1095-8649.1996.tb00060.x}{Lorenzen (1996)} and
#'           \href{https://www.nrcresearchpress.com/doi/10.1139/cjfas-2017-0035}{Miller & Hyun (2018)}.
#'           Default is to estimate a single model shared across all stocks and regions, but
#'           use \code{$means_map[s,r,1]} to fix or estimate the intercept for specific stocks, regions. See also \code{$logb_prior}
#'           and \code{$initial_b} configuring the slope on log scale.}
#'       }
#'     }
#'     \item{$initial_means}{array (n_stocks x n_regions x n_ages) of initial/mean M by stock, region and age. If \code{NULL}, initial 
#'       mean M-at-age values for a given stock and region are taken from the first row of the MAA matrix in the ASAP data file. If no
#'        ASAP data file, M = 0.2 is the default. If \code{$mean_model} is "weight-at-age" only 
#'       elements for the first age (\code{$initial_means[,,1]}) are used (for the intercept of log(M)).}
#'     \item{$means_map}{array (n_stocks x n_regions x n_ages) of NA or integers ( 0 <= max <= n_stocks * n_regions * n_ages) indicating 
#'       which ages to estimate (mean) M and whether to set any ages to be identical. E.g. in a model with 2 stock, 2 regions 
#'       and 6 ages with constant M estimated for each stock across regions and ages  \code{$M_ages_map[1,,] = 1} 
#'       and \code{$M_ages_map[2,,] = 2}. \code{$M_ages_map[1,1,] = c(NA,1,1,2,2,3)} will fix M for age 1 at the initial value, 
#'       and estimates for ages 2 and 3 are identical as are those for ages 4 and 5 and different from age 6+ for stock 1 and 
#'       region 1. If \code{NULL}, specifies all ages fixed at \code{M$initial_means}.  If \code{$mean_model} is "weight-at-age"
#'       these are used for all stocks and regions and only the elements for the first age (\code{$M_ages_map[,,1]}) 
#'       are used (for the intercept of log(M)).}
#'     \item{$b_model}{"constant","stock","region", "stock-region" defining whether parameter is constant, stock-specific, region-specific, 
#'       stock- and region-specific. Only used if \code{M$mean_model} = "weight-at-age".}
#'     \item{$b_prior}{T/F, should a N(mu, 0.08) prior (where mu = log(0.305) by default) be used on log_b? Based on Fig. 1 and Table 1 
#'       (marine fish) in  \href{https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1095-8649.1996.tb00060.x}{Lorenzen (1996)}. (Only used if 
#'       \code{$mean_model} is "weight-at-age").}
#'     \item{$intial_b}{if any elements of \code{$mean_model} is "weight-at-age", initial value for mean b for weight-at-age
#'       model.}
#'     \item{$re_model}{Character matrix (n_stocks x n_regions) of options for time- and age-varying (random effects) on M by stock and region.
#'       Possible values are:
#'       \describe{
#'         \item{"none"}{(default) No random effects by age or year.}
#'         \item{"iid_a"}{uncorrelated M by age, constant in time.}
#'         \item{"iid_y"}{uncorrelated M by year, constant all ages.}
#'         \item{"ar1_a"}{M correlated by age (AR1), constant in time.}
#'         \item{"ar1_y"}{M correlated by year (AR1), constant all ages.}
#'         \item{"iid_ay"}{M uncorrelated by year and age (2D).}
#'         \item{"ar1_ay"}{M correlated by year and age (2D AR1), as in \href{https://www.nrcresearchpress.com/doi/10.1139/cjfas-2015-0047}{Cadigan
#'           (2016)}.}
#'       }
#'     }
#'     \item{$re_map}{array (n_stocks x n_regions x n_ages) of integers ( 0 <= max <= n_ages) indicating which ages, for a given 
#'       stock and region, have random effects and whether to set RE for any ages to be identical. E.g. in a model with 2 stock, 
#'       2 regions and 6 ages, \code{$re_map[2,1,] = c(0,1,1,2,2,3)} will not estimate RE for age 1, and those for ages 
#'       2 and 3 are identical as are those for ages 4 and 5 and different from age 6+ for stock 2 and region 1. If \code{NULL}, 
#'       and \code{$re_model} specifies M random effects at age, at least two ages must be 
#'       specified for correlation among ages to be estimated.}
#'     \item{$sigma_vals}{n_stocks x n_regions matrix Initial standard deviation value to use for the M random effects. Values are not used 
#'       if \code{M$re_model} = "none". Otherwise, a single value. If unspecified all values are 0.1.}
#'     \item{$cor_vals}{n_stocks x n_regions x 2 array of initial correlation values to use for the M deviations. If unspecified all initial 
#'       values are 0. When \code{M$re_model} = 
#'       \describe{
#'         \item{"iid_a", "iid_y", "iid_ay" or "none"}{values are not used.}
#'         \item{"ar1_a" }{first value cor_vals[s,r,1] is used.}
#'         \item{"ar1_y" }{second value cor_vals[s,r,2] is used.}
#'         \item{"ar1_ay"}{First is for "age", second is for "year".}
#'       }
#'     }
#'     \item{$sigma_map}{n_stocks x n_region matrix of NA or integers indicating which random effects sd is estimated and whether to 
#'       set any to be identical. If not supplied a single sd will be estimated for any stock and region where $re_model is other than "none".}
#'     \item{$cor_map}{n_stocks x n_region matrix x 2 array of NA or integers indicating which random effects correlation parameters are estimated
#'       and whether to set any to be identical. If not supplied a single value for age and/or year will be estimated for any stock and region where 
#'       $re_model is other than "none", "iid_a", "iid_y".}
#'   }
#'
#' @return a named list with same elements as the input provided with natural mortality options modified.
#'
#' @seealso \code{\link{prepare_wham_input}} 
#'
#' @examples
#' \dontrun{
#' wham.dir <- find.package("wham")
#' path_to_examples <- system.file("extdata", package="wham")
#' asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))
#' input <- prepare_wham_input(asap3)
#' M = list(mean_model = "estimate-M")
#' input <- set_q(input, M = M) #estimate a constant M parameters
#' }
#'
#' @export
set_M = function(input, M)
{
  data = input$data
  par = input$par
  map = input$map
  asap3 = input$asap3
  input$log$M <- list()
  
  # elements of M: model, initial_means, means_map, logb_prior, intial_b, re_model, re_map, sigma_vals, cor_vals
  # data elements are: n_M_re, M_re_index, M_model, M_re_model, use_b_prior, log_b_model
  # data$waa_pointer_M is supplied in waa_info argument to set_WAA.
  # par elements are: Mpars, M_re, M_repars, log_b
  re_mods <- c("none", "iid_a", "iid_y", "ar1_a", "ar1_y", "iid_ay", "ar1_ay")
  #clear any map definitions that may exist. necessary because some configurations may not define map elements.
  map <- map[(!names(map) %in% c("log_b", "M_repars", "Mpars","M_re"))]
  
  #M_model length is n_regions; 1: Estimate-M, 2: f(WAA), 3: f(WAA) by stock
  data$M_model = 1 #Fixed-M
  data$log_b_model = 1 #constant if estimated
  data$use_b_prior = 0
  data$n_M_re <- matrix(1, data$n_stocks, data$n_regions)
  data$M_re_index <- array(1, dim = c(data$n_stocks, data$n_regions, data$n_ages))
  #for(s in 1:data$n_stocks) for(r in 1:data$n_regions) data$M_re_index[s,r,] <- 1:data$n_ages
  #M_re_model length = n_regions; 1 = none, 2 = IID, 3 = ar1_a, 4 = ar1_y, 5 = 2dar1
  data$M_re_model = matrix(1,data$n_stocks,data$n_regions) # default = no RE / 'none'
  
  par$Mpars <- array(NA, dim = c(data$n_stocks, data$n_regions, data$n_ages))
  map$Mpars <- par$Mpars
  par$Mpars[] <- log(0.2)
  
  par$M_re <- array(NA, dim = c(data$n_stocks, data$n_regions, data$n_years_model, data$n_ages))
  map$M_re <- par$M_re
  par$M_re[] <- 0# if estimating mean M for any ages, initialize yearly deviations at 0
  
  par$M_repars <- array(NA, dim = c(data$n_stocks, data$n_regions, 3))
  map$M_repars <- par$M_repars
  par$M_repars[] <- 0
  par$M_repars[,,1] <- log(0.1)
  
  if(is.null(M$mean_model)) M$mean_model <- "fixed-M"
  if(is.null(M$re_model)) M$re_model <- matrix("none", data$n_stocks, data$n_regions)
  if(is.null(M$b_prior)) M$logb_prior = FALSE
  
  if(is.null(M$intial_b)) M$intial_b = 0.305
  if(is.null(M$initial_mean)) {
    if(!is.null(asap3)){
      M$initial_means <- array(NA, dim = c(data$n_stocks, data$n_regions, data$n_ages))
      for(i in 1:length(asap3)) for(r in 1:data$n_regions) {
        M$initial_means[i,r,] <- asap3[[i]]$M[1,]
        # print(M$initial_means[i,r,])
        # print(log(asap3[[i]]$M))
        for(y in 1:data$n_years_model) {
          par$M_re[i,r,y,] <- log(asap3[[i]]$M[y,]) - log(M$initial_means[i,r,])
        }
      }
    }
    else M$intitial_means <- array(0.2, dim = c(data$n_stocks, data$n_regions, data$n_ages))
  }
  
  par$log_b <- matrix(log(M$intial_b), data$n_stocks, data$n_regions)
  map$log_b <- matrix(NA, data$n_stocks, data$n_regions)
  
  if(length(M$mean_model) != 1 & !is.character(M$mean_model)) stop("M$mean_model must be 'fixed-M', 'estimate-M', or 'weight-at-age'")
  M_mods = c("fixed-M", "estimate-M", "weight-at-age")
  if(!(M$mean_model %in% M_mods)) stop(paste0("M$mean_model must be one of these: ", paste0(M_mods, collapse=",")))
  data$M_model[] <- ifelse(M$mean_model %in% M_mods[1:2], 1, 2) #only needs to know f(age) or f(waa)
  
  if(!is.null(M$initial_means)){
    if(!is.array(M$initial_means)) stop("M$initial_means must now be an array with dimensions = c(n_stocks,n_regions,n_ages)") 
    dimsM = dim(M$initial_means)
    if(length(dimsM) != 3) stop("dimensions of M$initial_means must be c(n_stocks,n_regions,n_ages)")
    if(!all(dimsM == dim(par$Mpars))) stop("dimensions of M$initial_means must be c(n_stocks,n_regions,n_ages)")  
    par$Mpars[] <- log(M$initial_means)
  }
  
  if(!is.null(M$re_model)){
    if(!is.matrix(M$re_model)) stop("M$re_model must be a character n_stocks x n_regions matrix.")
    dimsM = dim(M$re_model)
    if(length(dimsM) != 2) stop("dimensions of M$re_model must be n_stocks x n_regions.")
    if(!all(dimsM == dim(par$Mpars)[1:2])) stop("dimensions of M$re_model must be n_stocks x n_regions.")  
    if(any(!(M$re_model %in% re_mods))) {
      stop(paste0("Each M$re_model must be one of the following: ", paste(re_mods, collapse = ",")))
    }
  }
  
  if(!is.null(M$means_map)){
    if(!is.array(M$means_map)) stop("M$means_map must be an array with dimensions = c(n_stocks,n_regions,n_ages)")
    dimsM = dim(M$means_map)
    if(length(dimsM) != 3) stop("dimensions of M$means_map must be c(n_stocks,n_regions,n_ages)")
    if(!all(dimsM == dim(par$Mpars))) stop("dimensions of M$means_map must be c(n_stocks,n_regions,n_ages)")
  }
  
  if(!is.null(M$re_map)){
    if(!is.array(M$re_map)) stop("M$re_map must be an array with dimensions = c(n_stocks,n_regions,n_ages)")
    dimsM = dim(M$re_map)
    if(length(dimsM) != 3) stop("dimensions of M$re_map must be c(n_stocks,n_regions,n_ages)")
    if(!all(dimsM == dim(par$Mpars))) stop("dimensions of M$re_map must be c(n_stocks,n_regions,n_ages)")
    if(any(!(M$re_map %in% 0:data$n_ages))) stop("Entries in M$re_map must be between 0 and n_ages.")
    k <- 0
    if(M$mean_model != "weight-at-age" ) for(s in 1:data$n_stocks) for(r in 1:data$n_regions) {
      if(M$re_model[s,r] %in% c("ar1_a","iid_a", "iid_ay","ar1_ay")){
        if(all(M$re_map[s,r,] == 0)) stop("random effects by age have been specified, but M$re_map are all 0.")
        temp <- M$re_map[s,r,]
        temp[which(temp==0)] <- NA
        data$n_M_re[s,r] <- length(unique(temp[!is.na(temp)]))
        data$M_re_index[s,r,] <- unclass(factor(M$re_map[s,r,]))
        ind <- unclass(factor(temp))
        if(M$re_model[s,r] %in% c("ar1_a","iid_a")) {
          for(y in 1:dim(map$M_re)[3]) map$M_re[s,r,y,] <- ind + k 
          k <- k + data$n_M_re[s,r]
        }
        if(M$re_model[s,r] %in% c("ar1_ay","iid_ay")) {
          for(y in 1:dim(map$M_re)[3]) {
            map$M_re[s,r,y,] <- ind + k 
            k <- k + data$n_M_re[s,r]
          }
        }
      } #for M$mean_model == "weight-at-age" (M = f(WAA)), These should already be set up right.
      if(M$re_model[s,r] %in% c("ar1_y","iid_y")){
        data$n_M_re[s,r] <- 1
        data$M_re_index[s,r,] <- 1
        for(a in 1:dim(map$M_re)[4]) map$M_re[s,r,,a] <- 1:dim(map$M_re)[3] + k 
        k <- k + dim(map$M_re)[3]
      } #for M$mean_model == "weight-at-age" (M = f(WAA)), These should already be set up right.
      k <- k + 1
    }
  } else{
    k <- 0
    if(M$mean_model != "weight-at-age") for(s in 1:data$n_stocks) for(r in 1:data$n_regions) {
      if(M$re_model[s,r] %in% c("ar1_a","iid_a")){
        data$n_M_re[s,r] <- data$n_ages
        data$M_re_index[s,r,] <- 1:data$n_ages
        for(y in 1:dim(map$M_re)[3]) map$M_re[s,r,y,] <- 1:data$n_ages + k
        k <- k +  data$n_ages
      } #for M$mean_model == "weight-at-age" (M = f(WAA)), These should already be set up right.
      if(M$re_model[s,r] %in% c("iid_ay","ar1_ay")){
        data$n_M_re[s,r] <- data$n_ages
        data$M_re_index[s,r,] <- 1:data$n_ages
        map$M_re[s,r,,] <- 1:(dim(map$M_re)[3]*data$n_ages) + k
        k <- k + dim(map$M_re)[3]*data$n_ages
      } #for M$mean_model == "weight-at-age" (M = f(WAA)), These should already be set up right.
      if(M$re_model[s,r] %in% c("ar1_y","iid_y")){
        data$n_M_re[s,r] <- 1
        data$M_re_index[s,r,] <- 1
        for(a in 1:dim(map$M_re)[4]) map$M_re[s,r,,a] <- 1:dim(map$M_re)[3] + k
        k <- k + dim(map$M_re)[3]
      } #for M$mean_model == "weight-at-age" (M = f(WAA)), These should already be set up right.
    }
    
  }
  
  if(is.null(M$means_map)){ #constant M or weight at age being used.
    if(M$mean_model != "fixed-M") { #Estimating constant M or WAA intercept
      map$Mpars[] = 1
      #only use first value of M from asap files if estimating M for constant or f(WAA)
      for(s in 1:data$n_stocks) for(r in 1:data$n_regions) {
        if(length(unique(par$Mpars[s,r,])) > 1) {
          input$log$M <- c(input$log$M, 
                           "M is estimated and no M$means_map is specified so only 1 mean M parameter), but for some stock/region, MAA has > 1 unique value.
          Initializing M at M$initial_means[1,1,1]. To avoid this warning without changing ASAP file, specify M$initial_means appropriately.\n")
        }
      }
      par$Mpars[] = par$Mpars[1,1,1]
    }
  }
  
  if((M$mean_model == "weight-at-age") & any(M$re_model == "ar1_a")) stop("Cannot estimate M as a function of weight at age and \n 
    random effects just on age. If you want random effects M parameters for on weight-at-age, set M$re_model = 'iid_y' or 'ar1_y'.")    
  if(any(M$re_model %in% c("ar1_a")) & !is.null(M$means_map)) { #this flag may need to be stock and region specific
    for(i in 1:data$n_stocks) for(j in 1:data$n_regions){
      if(length(unique(M$means_map[i,j,])) != 1) stop("If M$re_model is 'ar1_a', all M$means_map must be the same across ages \n
        (NA = fixed mean or an integer = estimated mean).")
    }
  }
  
  temp <- c("none", "ar1_a","ar1_y", "ar1_a","ar1_y", "ar1_ay", "ar1_ay")
  temp <- temp[match(M$re_model, re_mods)]
  data$M_re_model[] <- match(temp, c("none","ar1_a","ar1_y","ar1_ay"))
  
  if(M$mean_model == "weight-at-age"){
    if(is.null(M$b_model)){
      M$b_model <- "constant"
      input$log$M <- c(input$log$M, "M$b_model was not specified, so M as a function of weight at age will be used for all stocks and regions.\n")
    }
    if(length(M$b_model) != 1 | (!M$b_model %in% c("constant", "stock", "regions", "stock_region"))){
      stop("M$b_model must be a single value: 'constant', 'stock', 'region', or 'stock_region'.")
    }
    map$log_b[] <- 1 #log_b_model = 1: constant across stock, region
    map$Mpars[] <- 1 #the same "intercept": constant across stock, region
    if(M$b_model == "constant") par$Mpars[] <- par$Mpars[1,1,1]
    if(M$b_model == "stock") {
      data$log_b_model <- 2
      for(s in 1:data$n_stocks) {
        map$log_b[s,] <- s
        map$Mpars[s,,] <- s
        par$Mpars[s,,] <- par$Mpars[s,1,1]
      }
    }
    if(M$b_model == "region") {
      data$log_b_model <- 3
      for(r in 1:data$n_regions) {
        map$log_b[,r] <- r
        map$Mpars[,r,] <- r
        par$Mpars[,r,] <- par$Mpars[1,r,1]
      }
    }
    if(M$b_model == "stock_region") {
      data$log_b_model <- 4
      for(s in 1:data$n_stocks) for(r in 1:data$n_regions){
        map$log_b[s,r] <- r + (s-1) * data$n_regions
        map$Mpars[s,r,] <- r + (s-1) * data$n_regions
        par$Mpars[s,r,] <- par$Mpars[s,r,1]
      }
    }
    if(is.null(M$b_prior)){
      M$b_prior = FALSE
      input$log$M <- c(input$log$M, "M$b_prior was not specified, so prior for the b parameter of M_a = aW_a^b will not be used.\n")
    }
    if(length(M$b_prior) != 1 | !is.logical(M$b_prior)) stop("M$b_prior must be single value: TRUE or FALSE")
    if(M$b_prior) data$use_b_prior = 1
  }
  
  if(!is.null(M$sigma_vals)){
    if(!is.matrix(M$sigma_vals)) stop("M$sigma_vals must be a n_stocks x n_regions matrix.")
    dimsM = dim(M$sigma_vals)
    if(length(dimsM) != 2) stop("dimensions of M$sigma_vals must be n_stocks x n_regions.")
    if(!all(dimsM == dim(par$Mpars)[1:2])) stop("dimensions of M$sigma_vals must be n_stocks x n_regions.")
    if(any(M$sigma_vals<0)) stop("M$sigma_vals must be > 0.")
    for(s in 1:data$n_stocks)for(r in 1:data$n_regions) if(M$re_model[s,r] != "none") par$M_repars[s,r,1] <- log(M$sigma_vals[s,r])
  }
  if(!is.null(M$sigma_map)){
    if(!is.matrix(M$sigma_map)) stop("M$sigma_map must be a n_stocks x n_regions matrix.")
    dimsM = dim(M$sigma_map)
    if(length(dimsM) != 2) stop("dimensions of M$sigma_map must be n_stocks x n_regions.")
    if(!all(dimsM == dim(par$Mpars)[1:2])) stop("dimensions of M$sigma_map must be n_stocks x n_regions.")
  }
  
  #inv_trans_rho <- function(rho) 0.5 * (log(rho+1) - log(1-rho)) # 0.5 because needed transformation on cpp side is unusual.
  inv_trans_rho <- function(rho) (log(rho+1) - log(1-rho)) # 0.5 because needed transformation on cpp side is unusual.
  if(!is.null(M$cor_vals)){
    if(!is.array(M$cor_vals)) stop("M$cor_vals must be an array with dimensions = c(n_stocks,n_regions,2)")
    dimsM = dim(M$cor_vals)
    if(length(dimsM) != 3) stop("dimensions of M$cor_vals must be c(n_stocks,n_regions,2)")
    if(!all(dimsM == c(data$n_stocks,data$n_regions,2))) stop("dimensions of M$cor_vals must be c(n_stocks,n_regions,2)")
    if(any(abs(M$cor_vals) > 1)) stop("M$cor_vals must be > -1 and < 1.")
    for(s in 1:data$n_stocks)for(r in 1:data$n_regions) {
      if(M$re_model[s,r] %in% c("ar1_ay","ar1_y")) par$M_repars[s,r,3] <- inv_trans_rho(M$cor_vals[s,r,2])
      if(M$re_model[s,r] %in% c("ar1_ay","ar1_a")) par$M_repars[s,r,2] <- inv_trans_rho(M$cor_vals[s,r,1])
    }
  }
  if(!is.null(M$cor_map)){
    if(!is.array(M$cor_map)) stop("M$cor_map must be an array with dimensions = c(n_stocks,n_regions,2)")
    dimsM = dim(M$cor_map)
    if(length(dimsM) != 3) stop("dimensions of M$cor_map must be c(n_stocks,n_regions,2)")
    if(!all(dimsM == c(data$n_stocks,data$n_regions,2))) stop("dimensions of M$cor_map must be c(n_stocks,n_regions,2)")
  }
  
  # map mean M 
  if(is.null(M$means_map)){
    #estimate a single M for all ages or a single intercept for M = f(WAA) 
    if(M$mean_model != "fixed-M") map$Mpars[] = 1
  } else{ #use M$means_map
    if(M$mean_model == "weight-at-age"){ #only use the first age for the intercept
      map$Mpars[,,1] <- unclass(factor(M$means_map[,,1]))
    } else{
      #print(M$means_map)
      map$Mpars[] <- unclass(factor(M$means_map)) 
    }
  }
  
  #############
  #map M_repars
  if(!is.null(M$sigma_map)){
    M$sigma_map[which(M$re_model=="none")] <- NA
    M$sigma_map[] <- unclass(factor(M$sigma_map))
    for(s in 1:data$n_stocks)for(r in 1:data$n_regions) if(M$re_model[s,r] != "none") map$M_repars[s,r,,1] <- M$sigma_map[s,r]
  } else{ #all the same value
    for(s in 1:data$n_stocks)for(r in 1:data$n_regions) if(M$re_model[s,r] != "none") map$M_repars[s,r,1] <- 1
  }
  k <- 0
  #print(map$M_repars[1,1,])
  if(any(!is.na(map$M_repars))) k <- max(map$M_repars, na.rm= T)
  #print(k)
  
  if(!is.null(M$cor_map)){
    for(s in 1:data$n_stocks)for(r in 1:data$n_regions) {
      if(M$re_model[s,r] %in% c("none","iid_y","ar1_y", "iid_ay")) M$cor_map[s,r,1] <- NA
      if(M$re_model[s,r] %in% c("none","iid_a","ar1_a", "iid_ay")) M$cor_map[s,r,2] <- NA
    }
    M$cor_map[] <- k + unclass(factor(M$cor_map))
  } else { #same values for age and same values for year
    M$cor_map <- array(NA, dim = c(data$n_stocks, data$n_regions, 2))
    if(any(M$re_model %in% c("ar1_a", "ar1_ay"))) {
      M$cor_map[,,1] <- k + 1
      k <- k + 1
    }
    if(any(M$re_model %in% c("ar1_y", "ar1_ay"))) {
      M$cor_map[,,2] <- k + 1
      k <- k + 1
    }
  }
  #print(M$cor_map[1,1,])
  #print(map$M_repars[1,1,])
  for(s in 1:data$n_stocks)for(r in 1:data$n_regions) map$M_repars[s,r,2:3] <- M$cor_map[s,r,]
  #print(map$M_repars[1,1,])
  #############
  
  map$M_repars = factor(map$M_repars)
  map$M_re <- factor(map$M_re)
  map$Mpars = factor(map$Mpars)
  map$log_b = factor(map$log_b)
  
  input$data = data
  input$par = par
  input$map = map
  if(length(input$log$M)) input$log$M <- c("Natural Mortality: \n", input$log$M)
  #may need to update these 
  # projection data will always be modified by 'prepare_projection'
  input = set_proj(input, proj.opts = NULL) #proj options are used later after model fit, right?
  
  #set any parameters as random effects
  input$random = NULL
  input = set_random(input)
  input$options$M <- M
  return(input)
  
}

#' Specify model and parameter configuration for movement when input$data$n_regions > 1
#'
#' @param input list containing data, parameters, map, and random elements (output from \code{\link{wham::prepare_wham_input}})
#' @param move (optional) list specifying movement options: model, random effects, initial values, and parameters to fix (see details)
#' 
#' \code{move} specifies estimation options for movement.
#' If \code{NULL}, no movement will occur. If there are multiple regions, each stock will be modeled separately in different regions without movement. 
#' \code{move} is a list with the following entries:
#'   \describe{
#'     \item{$stock_move}{length = n_stocks, T/F whether each stock can move. If not provided then movement will be defined below for all stocks.}
#'     \item{$separable}{length = n_stocks, T/F whether movement should be modeled separably from mortality or both occuring simultaneously.}
#'     \item{$mean_model}{matrix (n_regions x (n_regions-1)): model options for fixed effects (mean and possibly variance) for each movement parameter are:
#'                    \describe{
#'                      \item{"none"}{(default) no movement between regions.}
#'                      \item{"constant"}{estimate a single movement rate to each region shared across all stocks, seasons, ages, years}
#'                      \item{"season"}{estimate movement rates  to each region for each season shared across all stocks, ages, years}
#'                      \item{"stock_constant"}{estimate a movement rate for each stock to each region shared across all seasons, ages, years}
#'                      \item{"stock_season"}{estimate a movement rate for each stock each season to each region shared across all ages, years}
#'                    }
#'     }
#'     \item{$age_re}{matrix (n_regions x (n_regions-1)): options for age random effects (for each mean parameter defined in \code{move$mean_model}):
#'                    \describe{
#'                      \item{"none"}{(default) no movement rate random effects by age.}
#'                      \item{"iid"}{independent movement rate random effects by age.}
#'                      \item{"ar1"}{allow first order autoregressive correlation of movement rate random effects by age.}
#'                    }
#'     }
#'     \item{$year_re}{matrix (n_regions x (n_regions-1)): options for yearly random effects (for each mean parameter defined in \code{move$mean_model}):
#'                    \describe{
#'                      \item{"none"}{(default) no movement rate random effects by year.}
#'                      \item{"iid"}{independent movement rate random effects by year.}
#'                      \item{"ar1"}{allow first order autoregressive correlation of movement rate random effects by year.}
#'                    }
#'     }
#'     \item{$prior_sigma}{array (n_stocks x n_seasons x n_regions x n_regions - 1) of sd parameters for normal priors on mean movement parameters on transformed scale (-Inf,Inf)}
#'     \item{$use_prior}{array (n_stocks x n_seasons x n_regions x n_regions - 1) 0/1 indicator whether to include prior for mean movement parameters in joint log-likelihood.}
#'     \item{$can_move}{array (n_stocks x n_seasons x n_regions x n_regions) 0/1 indicator whether movement can occur from one region to another.}
#'     \item{$must_move}{array (n_stocks x n_seasons x n_regions) 0/1 indicator whether movement from region must occur.}
#'     \item{$mean_vals}{array (n_stocks x n_seasons x n_regions x n_regions-1) of initial movement rate parameters *from* each region. Usage depends on \code{move$mean_model}.}
#'     \item{$sigma_vals}{array (n_stocks x n_seasons x n_regions x n_regions -1) of initial standard deviations to use for random effects. Usage depends on \code{move$age_re} and \code{move$year_re}.}
#'     \item{$cor_vals}{array (n_stocks x n_seasons x n_regions x n_regions - 1x 2) of initial correlation values to use for random effects. Usage depends on \code{move$age_re} and \code{move$year_re}.
#'        cor_vals[,,,,1] is for correlation with age, and cor_vals[,,,,2] is for correlation with year.}
#'   }
#'
#' @export
set_move = function(input, move)
{
  data = input$data
  par = input$par
  map = input$map
  input$log$move <- list()
  
  #clear any map definitions that may exist. necessary because some configurations may not define map elements.
  map <- map[(!names(map) %in% c("mu_repars", "mu_re"))]
  
  data$can_move = array(0, dim = c(data$n_stocks, data$n_seasons, data$n_regions, data$n_regions))
  if(!is.null(move$can_move)) {
    data$can_move[] = move$can_move
    if(sum(data$can_move) == 0) input$log$move <- c(input$log$move, "All of move$can_move = 0, so model assumes no movement for any stocks.\n")
  }
  
  data$use_mu_prior = array(0, dim = c(data$n_stocks, data$n_seasons, data$n_regions, data$n_regions-1))
  data$mig_type = rep(0,data$n_stocks)
  data$mu_model <- matrix(1, data$n_regions, data$n_regions-1)
  data$trans_mu_prior_sigma = array(0.1, dim = c(data$n_stocks, data$n_seasons, data$n_regions, data$n_regions-1))
  data$must_move = array(0, dim = c(data$n_stocks, data$n_seasons, data$n_regions))
  par$mu_prior_re = array(0, dim = c(data$n_stocks, data$n_seasons, data$n_regions, data$n_regions-1))
  par$trans_mu = array(0, dim = c(data$n_stocks, data$n_seasons, data$n_regions, data$n_regions-1))
  par$mu_re = array(0, dim = c(data$n_stocks, data$n_ages, data$n_seasons, data$n_years_model, data$n_regions, data$n_regions-1))
  par$mu_repars = array(0, dim = c(data$n_stocks, data$n_seasons, data$n_regions, data$n_regions-1, 3))
  map$mu_prior_re = array(NA, dim = dim(par$mu_prior_re))
  map$trans_mu = array(NA, dim = dim(par$trans_mu))
  map$mu_re = array(NA, dim = dim(par$mu_re))
  map$mu_repars = array(NA, dim = dim(par$mu_repars))
  
  if(!is.null(move$must_move)) data$must_move[] = move$must_move
  if(data$n_stocks>1) for(s in 1:data$n_stocks) if(sum(data$can_move[s,,,]) == 0) input$log$move <- c(input$log$move, paste0("Model assumes no movement for stock ", s, ".\n"))
  
  if(data$n_regions ==1 | sum(data$can_move)==0){
    map$trans_mu = factor(map$trans_mu)
    map$mu_repars = factor(map$mu_repars)
    map$mu_prior_re = factor(map$mu_prior_re)
    map$mu_re = factor(map$mu_re)
    input$data = data
    input$par = par
    input$map = map
    return(input)
  }
  
  if(!is.null(move$use_prior)) data$use_mu_prior[] = move$use_prior
  mean_mods <- c("none","constant", "season")
  mean_mods <- c(mean_mods, paste0("stock_", mean_mods[2:3]))
  if(is.null(move$mean_model)){
    if(is.null(move$can_move) | data$n_regions == 1) {
      move$mean_model <- matrix("none", data$n_regions, data$n_regions-1)
    } else{
      if(sum(move$can_move)>0) move$mean_model <- matrix("constant", data$n_regions, data$n_regions-1)
      else {
        if(data$n_regions==1) move$mean_model <- matrix("none", data$n_regions, data$n_regions-1)
        else move$mean_model <- matrix("constant", data$n_regions, data$n_regions-1)
      }
    }
    input$log$move <- c(input$log$move, paste0("\n move$mean_model was not specified and set to ", move$mean_model[1], " for all movement parameters based on data$n_regions and move$can_move if provided. \n"))
  }
  if(!is.matrix(move$mean_model)) stop(paste0("move$mean_model must be a n_regions x n_region-1 matrix filled with one of ", paste0(mean_mods, collapse = ", "), " when provided."))
  if(any(!move$mean_model %in% mean_mods)) stop(paste0("elements of move$mean_model must each be one of ", paste0(mean_mods, collapse = ", "), " when provided."))
  if(all(move$mean_model == "none")) {
    if(data$n_regions>1) input$log$move <- c(input$log$move, "all move$mean_model = 'none', so model assumes no movement for any stocks.\n")
    data$can_move[] <- 0
    map$trans_mu = factor(map$trans_mu)
    map$mu_repars = factor(map$mu_repars)
    map$mu_prior_re = factor(map$mu_prior_re)
    map$mu_re = factor(map$mu_re)
    input$data = data
    input$par = par
    input$map = map
    return(input)
  } else {
    if(sum(data$can_move)==0){ 
      input$log$move <- c(input$log$move, "move$model is not 'none', but all data$can_move = 0, so model assumes no movement for any stock.\n")
      map$trans_mu = factor(map$trans_mu)
      map$mu_repars = factor(map$mu_repars)
      map$mu_prior_re = factor(map$mu_prior_re)
      map$mu_re = factor(map$mu_re)
      input$data = data
      input$par = par
      input$map = map
      return(input)
    }
  }
  
  #if we've gotten this far then some parameters will be estimated
  #define data$mu_model and map for mu_re and mu_repars
  data$mu_model[] = 1 #this value needs to be between 1-8 after below
  data$mu_model[which(move$mean_model == "constant")] = 1 #already defined
  data$mu_model[which(move$mean_model == "stock_constant")] = 5 #already defined
  data$mu_model[which(move$mean_model == "season")] = 9 #already defined
  data$mu_model[which(move$mean_model == "stock_season")] = 13 #already defined
  re_mods <- c("none","iid", "ar1")
  if(is.null(move$age_re)) move$age_re <- matrix("none", data$n_regions, data$n_regions -1)
  if(!is.matrix(move$age_re)) stop(paste0("move$age_re must be a n_regions x n_region-1 matrix filled with one of ", paste0(re_mods, collapse = ", "), " when provided."))
  if(!all(dim(move$age_re) == data$n_regions - 0:1)) stop(paste0("move$age_re must be a n_regions x n_region-1 matrix filled with one of ", paste0(re_mods, collapse = ", "), " when provided."))
  if(is.null(move$year_re)) move$year_re = matrix("none", data$n_regions, data$n_regions -1)
  if(!is.matrix(move$year_re)) stop(paste0("move$year_re must be a n_regions x n_region-1 matrix filled with one of ", paste0(re_mods, collapse = ", "), " when provided."))
  if(!all(dim(move$year_re) == data$n_regions - 0:1)) stop(paste0("move$year_re must be a n_regions x n_region-1 matrix filled with one of ", paste0(re_mods, collapse = ", "), " when provided."))
  ind <- which(move$age_re != "none" & move$year_re == "none")
  data$mu_model[ind] <- data$mu_model[ind] + 1
  ind <- which(move$age_re == "none" & move$year_re != "none")
  data$mu_model[ind] <- data$mu_model[ind] + 2
  ind <- which(move$age_re != "none" & move$year_re != "none")
  data$mu_model[ind] <- data$mu_model[ind] + 3
  
  if(is.null(move$can_move)) for(r in 1:data$n_regions) {
    for(rr in 1:data$n_regions)
      rrless <- rr
    if(rr>r) rrless <- rr-1
    if(move$mean_model[r,rr] != "none") if(r!=rr){
      data$can_move[,,r,rr] <- 1
    }
  }
  
  use_mu_prior <- data$use_mu_prior
  map$trans_mu = array(NA, dim = dim(par$trans_mu)) 
  if(any(data$use_mu_prior>0)) {
    map$mu_prior_re <- array(NA, dim = dim(par$mu_prior_re))
    use_mu_prior[] = 0 #need to set use_mu_prior based on mean_model
  }
  
  #define map for mu_re and mu_repars
  i <- re_i <- 1
  for(r in 1:data$n_regions) {
    k <- 1
    for(rr in 1:data$n_regions) if(rr!= r) {
      if(data$mu_model[r,k] %in% 2:4) {
        if(sum(data$can_move[,,r,rr])>0) {
          map$mu_repars[,,r,k,1] <- i
          i <- i + 1
          if(data$mu_model[r,k] %in% c(2,4) & move$age_re[r,k] == "ar1") {
            map$mu_repars[,,r,k,2] <- i
            i <- i + 1
          }
          if(data$mu_model[r,k] %in% c(3,4) & move$year_re[r,k] == "ar1") {
            map$mu_repars[,,r,k,3] <- i
            i <- i + 1
          }
          for(s in 1:data$n_stocks) for(t in 1:data$n_seasons) { #mu_re are constant across stock, season
            if(data$mu_model[r,k] == 2) for(y in 1:data$n_years_model) map$mu_re[s,,t,y,r,k] <- (re_i-1)*data$n_ages + 1:data$n_ages
            if(data$mu_model[r,k] == 3) for(a in 1:data$n_ages) map$mu_re[s,a,t,,r,k] <- (re_i-1)*data$n_years_model + 1:data$n_years_model
            if(data$mu_model[r,k] == 4) for(a in 1:data$n_ages) map$mu_re[s,,t,,r,k] <- (re_i-1)*data$n_years_model*data$n_ages + 1:data$n_years_model*data$n_ages
          }
          re_i <- re_i + 1
        }
      }
      if(data$mu_model[r,k] %in% 6:8) {
        for(s in 1:data$n_stocks) {
          if(sum(data$can_move[s,,r,rr])>0) {
            map$mu_repars[s,,r,k,1] <- i
            i <- i + 1
            if(data$mu_model[r,k] %in% c(6,8) & move$age_re[r,k] == "ar1") {
              map$mu_repars[s,,r,k,2] <- i
              i <- i + 1
            }
            if(data$mu_model[r,k] %in% c(7,8) & move$year_re[r,k] == "ar1") {
              map$mu_repars[s,,r,k,3] <- i
              i <- i + 1
            }
            for(t in 1:data$n_seasons) { #mu_re are constant across season
              if(data$mu_model[r,k] == 6) for(y in 1:data$n_years_model) map$mu_re[s,,t,y,r,k] <- (re_i-1)*data$n_ages + 1:data$n_ages
              if(data$mu_model[r,k] == 7) for(a in 1:data$n_ages) map$mu_re[s,a,t,,r,k] <- (re_i-1)*data$n_years_model + 1:data$n_years_model
              if(data$mu_model[r,k] == 8) for(a in 1:data$n_ages) map$mu_re[s,,t,,r,k] <- (re_i-1)*data$n_years_model*data$n_ages + 1:data$n_years_model*data$n_ages
            }
            re_i <- re_i + 1
          }
        }
      }
      if(data$mu_model[r,k] %in% 10:12){ #season
        for(t in 1:data$n_seasons) {
          if(sum(data$can_move[,t,r,rr])>0) {
            map$mu_repars[,t,r,k,1] <- i
            i <- i + 1
            if(data$mu_model[r,k] %in% c(10,12) & move$age_re[r,k] == "ar1") {
              map$mu_repars[,t,r,k,2] <- i
              i <- i + 1
            }
            if(data$mu_model[r,k] %in% c(10,12) & move$year_re[r,k] == "ar1") {
              map$mu_repars[,t,r,k,3] <- i
              i <- i + 1
            }
            for(s in 1:data$n_stocks) { #mu_re are constant across stock
              if(data$mu_model[r,k] == 10) for(y in 1:data$n_years_model) map$mu_re[s,,t,y,r,k] <- (re_i-1)*data$n_ages + 1:data$n_ages
              if(data$mu_model[r,k] == 11) for(a in 1:data$n_ages) map$mu_re[s,a,t,,r,k] <- (re_i-1)*data$n_years_model + 1:data$n_years_model
              if(data$mu_model[r,k] == 12) for(a in 1:data$n_ages) map$mu_re[s,,t,,r,k] <- (re_i-1)*data$n_years_model*data$n_ages + 1:data$n_years_model*data$n_ages
            }
            re_i <- re_i + 1
          }
        }
      }
      if(data$mu_model[r,k] %in% 14:16){ #stock,season
        for(s in 1:data$n_stocks) for(t in 1:data$n_seasons) {
          if(data$can_move[s,t,r,rr]) {
            map$mu_repars[s,t,r,k,1] <- i
            i <- i + 1
            if(data$mu_model[r,k] %in% c(14,16) & move$age_re[r,k] == "ar1") {
              map$mu_repars[s,t,r,k,2] <- i
              i <- i + 1
            }
            if(data$mu_model[r,k] %in% c(15,16) & move$year_re[r,k] == "ar1") {
              map$mu_repars[s,t,r,k,3] <- i
              i <- i + 1
            } #mu_re are different for each stock, season
            if(data$mu_model[r,k] == 14) for(y in 1:data$n_years_model) map$mu_re[s,,t,y,r,k] <- (re_i-1)*data$n_ages + 1:data$n_ages
            if(data$mu_model[r,k] == 15) for(a in 1:data$n_ages) map$mu_re[s,a,t,,r,k] <- (re_i-1)*data$n_years_model + 1:data$n_years_model
            if(data$mu_model[r,k] == 16) for(a in 1:data$n_ages) map$mu_re[s,,t,,r,k] <- (re_i-1)*data$n_years_model*data$n_ages + 1:data$n_years_model*data$n_ages
            re_i <- re_i + 1
          }
        }
      }
      k <- k + 1
    }
  }
  # if(data$mu_model %in% 6:8){ #stock
  #   i <- re_i <- 1
  #   for(s in 1:data$n_stocks) for(r in 1:data$n_regions) {
  #     k <- 1
  #     for(rr in 1:data$n_regions) if(rr!= r) {
  #       if(sum(data$can_move[s,,r,rr])>0) {
  #         map$mu_repars[s,,r,k,1] <- i
  #         i <- i + 1
  #         if(data$mu_model %in% c(6,8)) {
  #           map$mu_repars[s,,r,k,2] <- i
  #           i <- i + 1
  #         }
  #         if(data$mu_model %in% c(7,8)) {
  #           map$mu_repars[s,,r,k,3] <- i
  #           i <- i + 1
  #         }
  #         for(t in 1:data$n_seasons) {
  #           if(data$mu_model == 6) for(y in 1:data$n_years_model) map$mu_re[s,,t,y,r,k] <- (re_i-1)*data$n_ages + 1:data$n_ages
  #           if(data$mu_model == 7) for(a in 1:data$n_ages) map$mu_re[s,a,t,,r,k] <- (re_i-1)*data$n_years_model + 1:data$n_years_model
  #           if(data$mu_model == 8) for(a in 1:data$n_ages) map$mu_re[s,,t,,r,k] <- (re_i-1)*data$n_years_model*data$n_ages + 1:data$n_years_model*data$n_ages
  #         }
  #         re_i <- re_i + 1
  #       }
  #       k <- k + 1
  #     }
  #   }
  # }
  # if(data$mu_model %in% 10:12){ #season
  #   i <- re_i <- 1
  #   for(t in 1:data$n_seasons) for(r in 1:data$n_regions) {
  #     k <- 1
  #     for(rr in 1:data$n_regions) if(rr!= r) {
  #       if(sum(data$can_move[,t,r,rr])>0) {
  #         map$mu_repars[,t,r,k,1] <- i
  #         i <- i + 1
  #         if(data$mu_model %in% c(10,12)) {
  #           map$mu_repars[,t,r,k,2] <- i
  #           i <- i + 1
  #         }
  #         if(data$mu_model %in% c(10,12)) {
  #           map$mu_repars[,t,r,k,3] <- i
  #           i <- i + 1
  #         }
  #         for(s in 1:data$n_stocks) {
  #           if(data$mu_model == 10) for(y in 1:data$n_years_model) map$mu_re[s,,t,y,r,k] <- (re_i-1)*data$n_ages + 1:data$n_ages
  #           if(data$mu_model == 11) for(a in 1:data$n_ages) map$mu_re[s,a,t,,r,k] <- (re_i-1)*data$n_years_model + 1:data$n_years_model
  #           if(data$mu_model == 12) for(a in 1:data$n_ages) map$mu_re[s,,t,,r,k] <- (re_i-1)*data$n_years_model*data$n_ages + 1:data$n_years_model*data$n_ages
  #         }
  #         re_i <- re_i + 1
  #       }
  #       k <- k + 1
  #     }
  #   }
  # }
  # if(data$mu_model %in% 14:16){ #stock,season
  #   i <- re_i <- 1
  #   for(s in 1:data$n_stocks) for(t in 1:data$n_seasons) for(r in 1:data$n_regions) {
  #     k <- 1
  #     for(rr in 1:data$n_regions) if(rr!= r) {
  #       if(data$can_move[s,t,r,rr]) {
  #         map$mu_repars[s,t,r,k,1] <- i
  #         i <- i + 1
  #         if(data$mu_model %in% c(14,16)) {
  #           map$mu_repars[s,t,r,k,2] <- i
  #           i <- i + 1
  #         }
  #         if(data$mu_model %in% c(15,16)) {
  #           map$mu_repars[s,t,r,k,3] <- i
  #           i <- i + 1
  #         }
  #         if(data$mu_model == 14) for(y in 1:data$n_years_model) map$mu_re[s,,t,y,r,k] <- (re_i-1)*data$n_ages + 1:data$n_ages
  #         if(data$mu_model == 15) for(a in 1:data$n_ages) map$mu_re[s,a,t,,r,k] <- (re_i-1)*data$n_years_model + 1:data$n_years_model
  #         if(data$mu_model == 16) for(a in 1:data$n_ages) map$mu_re[s,,t,,r,k] <- (re_i-1)*data$n_years_model*data$n_ages + 1:data$n_years_model*data$n_ages
  #         re_i <- re_i + 1
  #       }
  #       k <- k + 1
  #     }
  #   }
  # }
  
  #define map and "use" elements for mean mu (trans_mu) parameters and prior-based RE (mu_prior_re) parameters 
  i <- 1
  for(r in 1:data$n_regions) {
    k <- 1
    for(rr in 1:data$n_regions) if(rr!= r) {
      if(data$mu_model[r,k] %in% 1:4) { #constant
        if(sum(data$can_move[,,r,rr])>0) {
          if(sum(data$use_mu_prior[,,r,k])>0) { #mean mu is a random effect with defined prior
            use_mu_prior[1,1,r,k] <- 1 #only evaluate the likelihood once
            map$mu_prior_re[,,r,k] <- i #all values the same
          } else {
            map$trans_mu[,,r,k] <- i #all values the same
          }
          i <- i + 1
        }
      }
      if(data$mu_model[r,k] %in% 5:8) { #stock
        for(s in 1:data$n_stocks) {
          if(sum(data$can_move[s,,r,rr])>0) {
            #mean mu is a random effect with defined prior
            if(sum(data$use_mu_prior[s,,r,k])>0) {
              use_mu_prior[s,1,r,k] <- 1 #only evaluate the likelihood once
              map$mu_prior_re[s,,r,k] <- i #all values the same
            } else {
              map$trans_mu[s,,r,k] <- i #all values the same
            }
            i <- i + 1
          }
        }
      }
      if(data$mu_model[r,k] %in% 9:12) { #season
        for(t in 1:data$n_seasons) {
          if(sum(data$can_move[,t,r,rr])>0) {
            if(sum(data$use_mu_prior[,t,r,k])>0) { #mean mu is a random effect with defined prior
              use_mu_prior[1,t,r,k] <- 1 #only evaluate the likelihood once
              map$mu_prior_re[,t,r,k] <- i #all values the same
            } else {
              map$trans_mu[,t,r,k] <- i #all values the same
            }
            i <- i + 1
          }
        }
      }
      if(data$mu_model[r,k] %in% 13:16) { #stock_season
        for(s in 1:data$n_stocks) for(t in 1:data$n_seasons) {
          if(sum(data$can_move[s,t,r,rr])>0) {
            if(sum(data$use_mu_prior[s,t,r,k])>0) { #mean mu is a random effect with defined prior
              use_mu_prior[s,t,r,k] <- 1 #only evaluate the likelihood once
              map$mu_prior_re[s,t,r,k] <- i #all values the same
            } else {
              map$trans_mu[s,t,r,k] <- i #all values the same
            }
            i <- i + 1
          }
        }
      }
      k <- k + 1
    }
  }
  # if(data$mu_model %in% 5:8) { #stock
  #   i <- 1
  #   for(s in 1:data$n_stocks) for(r in 1:data$n_regions) {
  #     k <- 1
  #     for(rr in 1:data$n_regions) if(rr!= r) {
  #       if(sum(data$can_move[s,,r,rr])>0) {
  #         #mean mu is a random effect with defined prior
  #         if(sum(data$use_mu_prior[s,,r,k])>0) {
  #           use_mu_prior[s,1,r,k] <- 1 #only evaluate the likelihood once
  #           map$mu_prior_re[s,,r,k] <- i #all values the same
  #         } else {
  #           map$trans_mu[s,,r,k] <- i #all values the same
  #         }
  #         i <- i + 1
  #       }
  #       k <- k + 1
  #     }
  #   }
  # }
  # if(data$mu_model %in% 9:12) { #season
  #   i <- 1
  #   for(t in 1:data$n_seasons) for(r in 1:data$n_regions) {
  #     k <- 1
  #     for(rr in 1:data$n_regions) if(rr!= r) {
  #       if(sum(data$can_move[,t,r,rr])>0) {
  #         if(sum(data$use_mu_prior[,t,r,k])>0) { #mean mu is a random effect with defined prior
  #           use_mu_prior[1,t,r,k] <- 1 #only evaluate the likelihood once
  #           map$mu_prior_re[,t,r,k] <- i #all values the same
  #         } else {
  #           map$trans_mu[,t,r,k] <- i #all values the same
  #         }
  #         i <- i + 1
  #       }
  #       k <- k + 1
  #     }
  #   }
  # }
  # if(data$mu_model %in% 13:16) { #stock_season
  #   i <- 1
  #   for(s in 1:data$n_stocks) for(t in 1:data$n_seasons) for(r in 1:data$n_regions) {
  #     k <- 1
  #     for(rr in 1:data$n_regions) if(rr!= r) {
  #       if(sum(data$can_move[s,t,r,rr])>0) {
  #         if(sum(data$use_mu_prior[s,t,r,k])>0) { #mean mu is a random effect with defined prior
  #           use_mu_prior[s,t,r,k] <- 1 #only evaluate the likelihood once
  #           map$mu_prior_re[s,t,r,k] <- i #all values the same
  #         } else {
  #           map$trans_mu[s,t,r,k] <- i #all values the same
  #         }
  #         i <- i + 1
  #       }
  #       k <- k + 1
  #     }
  #   }
  # }
  data$use_mu_prior <- use_mu_prior 
  
  if(!is.null(move$separable)) {
    if(any(!move$separable)) {
      input$log$move <- c(input$log$move, "NOTE: movement and mortality must be assumed to occur seequentially for now.")
      move$separable[] <- TRUE
    }
    data$mig_type[] = as.integer(!move$separable)
  }
  if(data$n_regions>1) {
    if(length(unique(move$separable))==1) {
      if(move$separable[1] == 0) input$log$move <- c(input$log$move, "movement and mortality will be assumed to occur simultaneously within each season.\n")
      if(move$separable[1] == 1) input$log$move <- c(input$log$move, "movement and mortality will be assumed to occur sequentially within each season.\n")
    }
    else{
      for(s in 1:data$n_stocks){
        if(data$mig_type[s] == 1) input$log$move <- c(input$log$move, paste0("for stock ", s, ", movement and mortality will be assumed to occur simultaneously within each season.\n"))
        if(data$mig_type[s] == 0) input$log$move <- c(input$log$move, paste0("for stock ", s, ", movement and mortality will be assumed to occur sequentially within each season.\n"))
      }
    }
  }
  if(!is.null(move$prior_sigma)) data$trans_mu_prior_sigma[] = move$prior_sigma
  
  inv_trans_rho <- function(rho, s = 1) (log(rho+1) - log(1-rho))/s
  k <- 1
  if(!is.null(move$mean_vals)){
    for(s in 1:data$n_stocks) {
      if(data$mig_type[s]) { #simultaneous
        for(t in 1:data$n_seasons) for(r in 1:data$n_regions) {
          #ind <- which(data$can_move[s,t,r,-r]==1) #not necessary because trans_mu only used if can_move == 1
          par$trans_mu[s,,t,r,]<- log(move$mean_vals[s,t,r,])
        }
      } else { #separable
        for(t in 1:data$n_seasons) for(r in 1:data$n_regions) {
          #print(data$can_move[s,t,r,-r])
          #ind <- which(data$can_move[s,t,r,-r]==1) #not necessary because trans_mu only used if can_move == 1
          p <- move$mean_vals[s,t,r,] #n_r - 1 long!
          par$trans_mu[s,t,r,] <- log(p) - log(1-sum(p)) #additive
        }
      }
    }
  }
  if(!is.null(move$sigma_vals)){
    par$mu_repars[,,,1] = log(move$sigma_vals)
  }
  if(!is.null(move$cor_vals)){
    par$mu_repars[,,,,2] = inv_trans_rho(move$cor_vals[,,,,1])
    par$mu_repars[,,,,3] = inv_trans_rho(move$cor_vals[,,,,2])
  }
  map$trans_mu = factor(map$trans_mu)
  map$mu_repars = factor(map$mu_repars)
  map$mu_prior_re = factor(map$mu_prior_re)
  map$mu_re = factor(map$mu_re)
  
  input$data = data
  input$par = par
  input$map = map
  if(length(input$log$move)) input$log$move <- c("Movement: \n", input$log$move)
  #may need to update these 
  # projection data will always be modified by 'prepare_projection'
  input = set_proj(input, proj.opts = NULL) #proj options are used later after model fit, right?
  
  #set any parameters as random effects
  input$random = NULL
  input = set_random(input)
  input$options$move <- move
  return(input)
  
}

#' Specify model and parameter configuration for numbers at age
#'
#' @param input list containing data, parameters, map, and random elements (output from \code{\link{prepare_wham_input}})
#' @param NAA_re (optional) list specifying options for numbers-at-age random effects, initial parameter values, and recruitment model (see details)
#' 
#' If \code{NAA_re = NULL}, a traditional statistical catch-at-age model is fit (NAA = pred_NAA for all ages, deterministic). Otherwise,
#' \code{NAA_re} specifies numbers-at-age configuration. It is a list with the following possible entries:
#'   \describe{
#'     \item{$sigma}{Which ages allow deviations from the predicted NAA given NAA from previous time step? Must be a single character string described below or a vector
#'                    of length n_stocks. If length = 1, assumptions will be applied to all stocks. Options are specified with the strings:
#'                    \describe{
#'                      \item{"rec"}{Random effects on recruitment (deviations), all other ages deterministic}
#'                      \item{"rec+1"}{"Full state space" model with 2 estimated \code{sigma_a}, one for recruitment and one shared among other ages}
#'                    }
#'                  }
#'     \item{$sigma_map}{You can specify a more complex parameter configuration by entering an integer array (nstocks x n_regions x n_ages), where each entry points to the
#'                   NAA_sigma to use for that stock and age. E.g., for 2 stocks, 2 regions, and 6 ages, array(rep(c(1,2,2,3,3,3), each = 4),c(2,2,6)) will estimate 3 
#'                   sigmas, with recruitment (age-1) deviations having their own \code{sigma_R}, ages 2-3 sharing \code{sigma_2}, and ages 4-6 sharing \code{sigma_3}. 
#'                   All parameters being the same for both stocks and across regions. The user must be sure that a compatible \code{NAA_re$sigma} configuration is defined.
#'                   Values are not used if recruit_model = 1 and \code{NAA_re$sigma} is not specified.
#'                  }
#'     \item{$sigma_vals}{Initial standard deviation values to use for the NAA deviations. Values are not used if recruit_model = 1 and \code{NAA_re$sigma} is
#'                  not specified. Otherwise when \code{NAA_re$sigma} =
#'                  \describe{
#'                    \item{"rec"}{must be a list (length = n_stocks) of single values .}
#'                    \item{"rec+1"}{Either 1) a list (length = n_stocks) of 2 values must be specified. First is for the first age class (recruits), second is for all other ages,
#'                      or 2) an array (nstocks x nregions x nages)}{only the value for age 1 in the spawning region is used.}
#'                  }
#'                  If \code{NAA_re$sigma_map} is defined, the user must ensure that the configuration is compatible with \code{NAA_re$sigma_vals}
#'                }
#'     \item{$cor}{Correlation structure for the NAA deviations. Must be a single character string described below or a vector
#'                    of length n_stocks. If length = 1, assumptions will be applied to all stocks. Options are:
#'                  \describe{
#'                    \item{"iid"}{NAA deviations vary by year and age, but uncorrelated.}
#'                    \item{"ar1_a"}{NAA deviations correlated by age (AR1).}
#'                    \item{"ar1_y"}{NAA deviations correlated by year (AR1).}
#'                    \item{"2dar1"}{NAA deviations correlated by year and age (2D AR1).}
#'                  }
#'                }
#'     \item{$cor_vals}{Initial correlation values to use for the NAA deviations. If unspecified all initial values are 0. When \code{NAA_re$cor} = 
#'                  \describe{
#'                    \item{"iid"}{values are not used.}
#'                    \item{"ar1_a" or "ar1_y"}{must be a list (length = n_stocks) each with a single value.}
#'                    \item{"2dar1"}{must be a list (length = n_stocks) each with 2 values. First is for "age", second is for "year".}
#'                  }
#'                }
#'     \item{$decouple_recruitment}{T/F determining whether correlation structure of recruitment is independent of RE deviations for older ages 
#'        (default = FALSE). Only applicable for \code{NAA_re$sigma = "rec+1"} and correlation across ages is specified. If TRUE and \code{NAA_re$cor = "ar1_a"}, only deviations for ages>1 
#'        have the correlation structure. If TRUE and NAA_re$cor is "ar1_y" or "2dar1" separate year correlation parameters are estimated for recruitment and older
#'        ages.
#'     }
#'     \item{$N1_model}{Integer vector (n_stocks) determining which way to model the initial numbers at age:
#'       \describe{
#'          \item{"age-specific-fe"}{(default) age- and region-specific fixed effects parameters}
#'          \item{"equilibrium"}{2 fixed effects parameters: an initial recruitment and an instantaneous fishing mortality rate to generate an equilibruim abundance at age.}
#'          \item{"iid-re"}{(default) age- and region-specific iid random effects parameters. 2 parameters: mean and sd for log NAA}
#'          \item{"ar1-re"}{(default) age- and region-specific random effects parameters. 3 parameters: mean and sd, and cor for log NAA}
#'       }
#'     }
#'     \item{$N1_pars}{An (n_stocks x n_regions x n_ages) array. If N1_model = 0, then this should be filled with the initial values to use for abundance at age by stock and region in the first year. 
#'        If N1_model = 1 (equilibrium assumption), only the first two values in the ages dimension are used: the (s,r,1) value is recruitment for stock (and region) and (s,r,2) is the fully-selected 
#'        equilibrium fishing mortality rate generating the rest of the numbers at age in the first year.
#'     }
#'     \item{$recruit_model}{Integer vector (n_stocks) determining how to model recruitment. Overrides \code{recruit_model} argument to \code{prepare_wham_input}. Must make sure \code{NAA_re$sigma}, \code{NAA_re$cor}
#'        and \code{ecov} are properly specified.
#'       \describe{
#'           \item{1}{SCAA, estimating annual recruitments as fixed effects or a random walk if NAA_re$sigma specified}
#'           \item{2}{estimating a mean recruitment with annual recruitments as random effects}
#'           \item{3}{Beverton-Holt stock-recruitment with annual recruitments as random effects}
#'           \item{4}{Ricker stock-recruitment with annual recruitments as random effects}
#'       }
#'     }
#'     \item{$recruit_pars}{list (length = n_stocks) of vectors of initial parameters for recruitment model. If $recruit_model is 3 or 4, parameters are "alpha" and "beta".
#'     }
#'   }
#'
#' @return a named list with same elements as the input provided with abundance modeling options modified.
#'
#' @seealso \code{\link{prepare_wham_input}} 
#'
#' @examples
#' \dontrun{
#' wham.dir <- find.package("wham")
#' path_to_examples <- system.file("extdata", package="wham")
#' asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))
#' input <- prepare_wham_input(asap3)
#' NAA = list(sigma = "rec")
#' input <- set_q(input, NAA_re = NAA) #estimate recruitment as random effects
#' }
#'
#' @export
set_NAA = function(input, NAA_re=NULL)
{
  
  data = input$data
  par = input$par
  map = input$map
  asap3 = input$asap3
  inv_trans_rho <- function(rho, s = 1) (log(rho+1) - log(1-rho))/s 
  input$log$NAA <- list()
  #clear any map definitions that may exist. necessary because some configurations may not define map elements.
  map <- map[(!names(map) %in% c("mean_rec_pars", "log_N1", "log_NAA_sigma", "trans_NAA_rho","logR_proj", "log_NAA"))]
  
  #if(is.null(input$asap3)) asap3 = NULL
  #else asap3 = input$asap3
  
  
  #set up initial NAA
  #0: just age-specific numbers at age
  data$N1_model = rep(0, data$n_stocks)
  
  data$decouple_recruitment <- 0 #until all examples, tests, vignettes are changed
  # data$decouple_recruitment <- 1 #decouple is default now!
  if(!is.null(NAA_re$decouple_recruitment)) data$decouple_recruitment <- as.integer(NAA_re$decouple_recruitment)
  
  par$log_N1 = array(0,dim = c(data$n_stocks,data$n_regions,data$n_ages))
  map$log_N1 = array(NA,dim = c(data$n_stocks,data$n_regions,data$n_ages))
  par$N1_repars = array(0,dim = c(data$n_stocks,data$n_regions,3))
  map$N1_repars = array(NA,dim = c(data$n_stocks,data$n_regions,3))
  init_NAA = log(exp(10)*exp(-(0:(data$n_ages-1))*0.2))
  if(!is.null(NAA_re$N1_model)) {
    options = c("age-specific-fe", "equilibrium","iid-re", "ar1-re")
    k = 1
    for(s in 1:data$n_stocks) {
      if(!(NAA_re$N1_model[s] %in% options)) stop("NAA_re$N1_model must all be 'age-specific-fe', 'equilibrium', 'iid-re' or 'ar1-re'.")
      if(NAA_re$N1_model[s] == options[1]) data$N1_model[s] = 0
      if(NAA_re$N1_model[s] == options[2]) data$N1_model[s] = 1
      if(NAA_re$N1_model[s] %in% options[3:4]) {
        data$N1_model[s] = 2
        map$N1_repars[s,,1] = k
        map$N1_repars[s,,2] = k + 1
        k = k + 2
        if(NAA_re$N1_model[s] == options[4]){ #ar1 with age
          map$N1_repars[s,,2] = k
          k = k + 1
        }
      }
    }
  }
  map$N1_repars = factor(map$N1_repars)
  k = 1
  if(any(data$N1_model == 2) & any(data$N1_model != 2)) stop("If any initial numbers at age are treated as RE, then all must.")
  for(s in 1:data$n_stocks) {
    if(data$N1_model[s] == 0){
      for(r in 1:data$n_regions) for(a in 1:data$n_ages) {
        if(data$NAA_where[s,r,a] == 1) {
          if(!is.null(asap3)) {
            par$log_N1[s,r,a] = log(asap3[[s]]$N1_ini[a]) # use N1_ini values from asap3 file
          } else{
            par$log_N1[s,r,a] = init_NAA[a]
          }
          map$log_N1[s,r,a] = k
          k <- k + 1
        }
      }
    }
    if(data$N1_model[s] == 1) { #equilibrium assumption, 2 pars per stock
      par$log_N1[s,data$spawn_regions[s],1:2] = c(10,log(0.1)) # allowed in wham.cpp but no option to set here (must be hard-coded after calling prepare_wham_input)
      map$log_N1[s,data$spawn_regions[s],1:2] = k + 0:1
      k = k + 2
    }
    if(data$N1_model[s] == 2) { #RE
      for(r in 1:data$n_regions) for(a in 1:data$n_ages) {
        if(data$NAA_where[s,r,a] ==1) {
          par$log_N1[s,r,a] = init_NAA[a]
          map$log_N1[s,r,a] = k
          k <- k + 1
        }
      }
    }
  }
  if(!is.null(NAA_re[["N1_pars"]])){
    par$log_N1[] = log(NAA_re$N1_pars)
  }
  map$log_N1 = factor(map$log_N1)
  
  # NAA_re options for beyond year 1
  # default = SCAA for each stock
  data$NAA_re_model = rep(0, data$n_stocks)
  par$log_NAA_sigma = array(0, c(data$n_stocks, data$n_regions, data$n_ages))
  map$log_NAA_sigma = array(NA, c(data$n_stocks, data$n_regions, data$n_ages))
  par$trans_NAA_rho = array(0,c(data$n_stocks, data$n_regions, 3))
  map$trans_NAA_rho = array(NA,c(data$n_stocks, data$n_regions, 3))
  par$log_NAA = array(10,dim = c(data$n_stocks, data$n_regions, data$n_years_model-1, data$n_ages))
  map$log_NAA = array(NA,dim = c(data$n_stocks, data$n_regions, data$n_years_model-1, data$n_ages))
  for(s in 1:data$n_stocks){
    map$log_NAA[s,data$spawn_regions[s],,1] <- 1 #change to unique values later
  }
  
  #if(is.null(NAA_re$sigma)){ 
  #  data$n_NAA_sigma <- 0
  #  data$NAA_sigma_pointers <- rep(1,data$n_ages)
  #}
  if(!is.null(NAA_re$sigma)){
    k <- 1
    if(!length(NAA_re$sigma) %in% c(1,data$n_stocks)) stop("NAA_re$sigma length must be 1 or equal to the number of stocks.")
    if(length(NAA_re$sigma) == 1) {
      input$log$NAA <- c(input$log$NAA, paste0("\n Same NAA_re$sigma being used for all stocks (", NAA_re$sigma[[1]][1], ").\n"))
      #NAA_re$sigma = rep(list(NAA_re$sigma), data$n_stocks)
      map$log_NAA_sigma[,,1] <- k
      data$NAA_re_model[] <- 1
      #for(s in 1:data$n_stocks) for(r in 1:data$n_regions) if(data$NAA_where[s,r,1]==1) map$log_NAA[s,r,,1] <- 1 #change to unique values later
      if(NAA_re$sigma[[1]][1] == "rec+1"){ # default state-space model with two NAA_sigma (one for recruits, one for ages > 1)
        map$log_NAA_sigma[,,-1] <- k+1
        data$NAA_re_model[] <- 2
        #also RE for ages>1
        for(s in 1:data$n_stocks) for(r in 1:data$n_regions) for(a in 2:data$n_ages) if(data$NAA_where[s,r,a]==1) map$log_NAA[s,r,,a] <- 1 #change to unique values later
      }
    } else {
      for(s in 1:data$n_stocks) {
        if(NAA_re$sigma[[s]][1] == "rec"){
          data$NAA_re_model[s] <- 1
          map$log_NAA_sigma[s,,1] <- k
          #below is already done above for SCAA
          #map$log_NAA[s,data$spawn_regions[s],,1] = 1 #change to unique values later
          k <- k + 1
          #data$n_NAA_sigma <- 1
          #data$NAA_sigma_pointers <- rep(1,data$n_ages)
        } else {
          if(NAA_re$sigma[[s]][1] == "rec+1"){ # default state-space model with two NAA_sigma (one for recruits, one for ages > 1)
            data$NAA_re_model[s] <- 2
            for(r in 1:data$n_regions) map$log_NAA_sigma[s,r,] <- c(k, rep(k+1, data$n_ages-1))
            for(r in 1:data$n_regions) for(a in 1:data$n_ages) if(data$NAA_where[s,r,a]==1) map$log_NAA[s,r,,a] <- 1 #change to unique values later
            k <- k + 2
            #data$n_NAA_sigma <- 2
            #data$NAA_sigma_pointers <- c(1,rep(2,data$n_ages-1))
          } else {
            if(length(NAA_re$sigma[[s]]) != data$n_ages) stop("each element of NAA_re$sigma must either be 'rec' (random effects on recruitment only), 
    'rec+1' (random effects on all NAA with ages > 1 sharing sigma_a,
    or a vector with length == n.ages specifying which sigma_a to use for each age.")
            #if(length(NAA_re$sigma[[s]]) == data$n_ages){
            if(any(is.na(unique(NAA_re$sigma[[s]])))){
              #use SCAA because of na values
              #data$n_NAA_sigma <- 0
              #data$NAA_sigma_pointers <- rep(1,data$n_ages)            
            } else {
              tmp <- unique(NAA_re_sigma[[s]])
              ind <- 1:length(tmp)
              ind <- ind[match(NAA_re_sigma[[s]],tmp)] - 1
              for(r in 1:data$n_regions) map$log_NAA_sigma[s,r,] <- k + ind
              for(r in 1:data$n_regions) for(a in 1:data$n_ages) if(data$NAA_where[s,r,a]==1) map$log_NAA[s,r,,a] <- 1 #change to unique values later
              k <- max(k + ind, na.rm=T)
              #data$n_NAA_sigma <- max(unique(NAA_re$sigma), na.rm=T)
              #data$NAA_sigma_pointers <- NAA_re$sigma
            }
            #}
          }
        }
      }
    }
    if(!is.null(NAA_re$sigma_vals)) {
      if(!is.array(NAA_re$sigma_vals)) {
        if(length(NAA_re$sigma_vals) != data$n_stocks) stop("NAA_re$sigma_vals must be a list with length = number of stocks")
        for(s in 1:data$n_stocks) {
          if(!(length(NAA_re$sigma_vals[[s]]) %in% c(1,data$n_ages))) stop(paste0("length of NAA_re$sigma_vals[[s]] must be 1 or ", data$n_ages, "."))
          par$log_NAA_sigma[s,,] <- log(NAA_re$sigma_vals[[s]])
        }
      } else {
        if(any(dim(NAA_re$sigma_vals) != c(data$n_stocks, data$n_regions, data$n_ages))) stop("dimensions of NAA_re$sigma_vals array are not c(nstocks,nregions,nages)")
        par$log_NAA_sigma[] <- log(NAA_re$sigma_vals)
      }
    }
    
    if(is.null(NAA_re$cor)) NAA_re$cor <- "iid"
    if(!(length(NAA_re$cor) %in% c(1,data$n_stocks))) stop("NAA_r$cor must have length 1 or n_stocks")
    k <- 0
    constant <- length(NAA_re$cor)==1
    for(s in 1:data$n_stocks) {
      ind <- ifelse(constant,1,s)
      if(!constant) k <- max(c(0, map$trans_NAA_rho), na.rm = TRUE)
      if(!is.null(NAA_re$cor[ind])){
        if(!NAA_re$cor[ind] %in% c("iid","ar1_a","ar1_y","2dar1")) stop("NAA_re$cor must be one of 'iid','ar1_a','ar1_y','2dar1'")
        if(NAA_re$cor[ind] %in% c("ar1_a")) map$trans_NAA_rho[s,,1] <- k+1
        if(NAA_re$cor[ind] %in% c("ar1_y")) {
          map$trans_NAA_rho[s,,2] <- k+1
          if(data$decouple_rec) map$trans_NAA_rho[s,,3] <- k+2
        }
        if(NAA_re$cor[ind] == "2dar1") for(r in 1:data$n_regions) {
          map$trans_NAA_rho[s,r,1:2] <- k + 1:2
          if(data$decouple_rec) map$trans_NAA_rho[s,3] <- k + 3
        } else {
          # NAA_re$cor[s] <- 'iid'
        }
      }
    }
    # k <- 1
    # if(length(NAA_re$cor) == 1) {
    #   input$log$NAA <- c(input$log$NAA, paste0("\n Same NAA_re$cor being used for all stocks (",NAA_re$cor[[1]][1],").\n"))
    #   #NAA_re$cor = rep(list(NAA_re$cor), data$n_stocks)
    #   if(!NAA_re$cor[[1]][1] %in% c("iid","ar1_a","ar1_y","2dar1")) stop("NAA_re$cor must be one of 'iid','ar1_a','ar1_y','2dar1'")
    #   if(NAA_re$cor[[1]][1] %in% c("ar1_a","2dar1")) map$trans_NAA_rho[,,1] <- k
    #   if(NAA_re$cor[[1]][1] == "ar1_y") map$trans_NAA_rho[,,2] <- k
    #   if(NAA_re$cor[[1]][1] == "2dar1") map$trans_NAA_rho[,,2] <- k + 1
    #   #if(NAA_re$cor[[1]][1] != "iid") k <- max(map$trans_NAA_rho, na.rm = TRUE)
    # } else {
    #   for(s in 1:data$n_stocks) {
    #     if(!is.null(NAA_re$cor[[s]])){
    #       if(!NAA_re$cor[[s]] %in% c("iid","ar1_a","ar1_y","2dar1")) stop("NAA_re$cor[[s]] must be one of 'iid','ar1_a','ar1_y','2dar1'")
    #       if(NAA_re$cor[[s]] == "ar1_a") map$trans_NAA_rho[s,,1] <- k
    #       if(NAA_re$cor[[s]] == "ar1_y") map$trans_NAA_rho[s,,2] <- k
    #       if(NAA_re$cor[[s]] == "2dar1") for(r in 1:data$n_regions) map$trans_NAA_rho[s,r,1:2] <- k + 0:1
    #       if(any(!is.na(map$trans_NAA_rho))) k <- max(map$trans_NAA_rho, na.rm = TRUE) + 1
    #       #if(NAA_re$cor[[s]] != "iid") k <- max(map$trans_NAA_rho, na.rm = TRUE) + 1
    #     } else {
    #       NAA_re$cor[[s]] <- 'iid'
    #     }
    #   }
    # }
    if(!is.null(NAA_re$cor_vals)) {
      if(length(NAA_re$cor_vals) != data$n_stocks) stop("NAA_re$cor_vals must be a list with length = number of stocks")
      for(s in 1:data$n_stocks) {
        if(!is.null(NAA_re$cor_vals[[s]])) {
          #FIXME: add in checks for decoupled recruitment, which will be default
          if(length(NAA_re$cor_vals[[s]]) > 3) stop(paste0("length of NAA_re$cor_vals[[s]] is not consistent with other elements of NAA_re$cor."))
          if(length(NAA_re$cor_vals[[s]]) == 2)  for(r in 1:data$n_regions) par$trans_NAA_rho[s,r,1:2] <- inv_trans_rho(NAA_re$cor_vals[[s]])
          if(length(NAA_re$cor_vals[[s]]) == 1) {
            if(NAA_re$cor[[s]] == "ar1_a") {
              par$trans_NAA_rho[s,,1] <- inv_trans_rho(NAA_re$cor_vals[[s]])
            }
            if(NAA_re$cor[[s]] == "ar1_y") {
              par$trans_NAA_rho[s,,2] <- inv_trans_rho(NAA_re$cor_vals[[s]])
            }
          }
        }
      }
    }
    if(!is.null(NAA_re$sigma_map)) {
      if(!is.array(NAA_re$sigma_map)) stop("NAA_re$sigma_map must be an array with dimensions = nstocks x nregions x nages")
      if(any(dim(NAA_re$sigma_map) != c(data$n_stocks, data$n_regions, data$n_ages))) stop("dimensions of NAA_re$sigma_vals array are not c(nstocks,nregions,nages)")
      map$log_NAA_sigma[] <- NAA_re$sigma_map
    }
  }
  #map$trans_NAA_rho[which(!is.na(map$trans_NAA_rho))] <- 1:sum(!is.na(map$trans_NAA_rho))
  map$trans_NAA_rho <- factor(map$trans_NAA_rho)
  map$log_NAA[which(!is.na(map$log_NAA))] <- 1:sum(!is.na(map$log_NAA))
  map$log_NAA <- factor(map$log_NAA)
  map$log_NAA_sigma <- factor(map$log_NAA_sigma)
  
  if(any(data$recruit_model > 2 & data$NAA_re_model == 0)) input$log$NAA <- c(input$log$NAA, "NOTE: SCAA model specified, yearly recruitment deviations estimated as fixed effects. Stock-recruit function also specified. WHAM will fit the SCAA model but without estimating a stock-recruit function.
    This message will not appear if you set recruit_model = 2 (random about mean).")
  
  #set up recruitment
  if(!is.null(NAA_re$recruit_model)) {
    data$recruit_model[] = NAA_re$recruit_model #overrides recruit_model argument to wham::prepare_wham_input
    for(s in 1:data$n_stocks) if(data$recruit_model[s] > 1 & data$NAA_re_model[s] == 0) {
      stop("NAA_re$recruit_model[s] > 1 has been specified, but NAA_re$sigma[[s]] must either be 'rec' (random effects on recruitment only), 
      'rec+1' (random effects on all NAA with ages > 1 sharing sigma_a,
      or a vector with length == n.ages specifying which sigma_a to use for each age.")
    }
  }
  par$mean_rec_pars = matrix(0, data$n_stocks, 2)
  map$mean_rec_pars = matrix(NA, data$n_stocks, 2)
  
  for(s in 1:data$n_stocks) {
    if(!is.null(NAA_re$recruit_pars[[s]])){
      if(data$recruit_model[s] == 2) par$mean_rec_pars[s,1] = log(NAA_re$recruit_pars[[s]][1])
      if(data$recruit_model[s] %in% 3:4) par$mean_rec_pars[s,1:2] = log(NAA_re$recruit_pars[[s]][1:2])
    } else{
      if(data$recruit_model[s]==2) {
        if(!is.null(asap3)) par$mean_rec_pars[s,1] = log(asap3[[s]]$N1_ini[1]) # initialize R0 at initial age-1
        else par$mean_rec_pars[s,1] = 10
      }  
      if(data$recruit_model[s]==4) par$mean_rec_pars[s,2] = -10
    }
    if(data$NAA_re_model[s] > 0){
      if(data$recruit_model[s]==2) map$mean_rec_pars[s,1] = s
      if(data$recruit_model[s]>2) map$mean_rec_pars[s,] = 1:2 + (s-1) * 2
    }
  }
  map$mean_rec_pars = factor(map$mean_rec_pars)
  
  input$data = data
  input$par = par
  input$map = map
  if(length(input$log$NAA))	input$log$NAA <- c("NAA: \n", input$log$NAA)
  #may need to update these 
  # projection data will always be modified by 'prepare_projection'
  #input = wham:::set_proj(input, proj.opts = NULL) #proj options are used later after model fit, right?
  
  #set any parameters as random effects
  #print(sort(names(input$data)))
  input$random = NULL
  input = set_random(input)
  input$options$NAA_re <- NAA_re
  return(input)
}


#' Set up observation vector that is used by the model for likelihood calculations and one-step-ahead residuals.
#'
#' @param input list containing data, parameters, map, and random elements (output from \code{\link{wham::prepare_wham_input}})
#'
#' @return the same input list as provided, but with $obs and $obsvec configured.
#' This is run after any changes have been made to the data
#' @export
set_osa_obs = function(input)
{
  data = input$data
  par = input$par
  map = input$map
  input$log$osa_obs <- list()
  obs.colnames <- c("year","fleet","age","type","val")
  obs <- data.frame(matrix(ncol = length(obs.colnames), nrow = 0))
  colnames(obs) <- obs.colnames
  
  # 5 components: Ecov, fleet catch (log), index catch (log), paa catch, paa index
  # do all Ecov first
  # 1. Ecov
  if(!all(data$Ecov_use_obs==0)){
    x <- as.data.frame(data$Ecov_obs)
    x[data$Ecov_use_obs==0] <- NA # only include index data to fit in obsvec
    colnames(x) <- paste0("Ecov_", 1:data$n_Ecov)
    x$year <- seq(from=input$years_Ecov[1]-input$years[1]+1, length.out=data$n_years_Ecov) # don't assume Ecov and model years are the same
    tmp <- tidyr::pivot_longer(x, cols = -year, values_to = 'val', names_to="fleet")
    tmp <- tmp[complete.cases(tmp),]
    tmp$age <- NA
    tmp$type <- "Ecov"
    obs <- rbind(obs, tmp[, obs.colnames])
  }
  
  
  #x <- data.frame(NA,nrow = data$n_years_model, ncol = data$n_indices)
  #x$year <- 1:data$n_years_indices # code assumes you have index and catch in all years - this will not work if we extend catch to 1930s
  
  #find out if there are any age-specific selectivity parameters set to 0 and not estimated
  #if the selblock is used, then we will reduce the observation vector to omit these ages. 
  #an OSA residual for ages with selectivity assumed to be zero will be useless and causes estimation problems with conditional construction of 
  #age comp likelihoods.
  ind <- which(data$selblock_models == 1) #age-specific
  ages_omit <- lapply(1:data$n_selblocks, function(x) integer(0))
  if(length(ind)){ #some
    if(is.null(map$logit_selpars)) {
      input$log$osa_obs <- c(input$log$osa_obs, "NOTE: input$map$logit_selapars has not been specified. This will only work for simulation.")
    } else {
      lsp_map <- matrix(as.integer(map$logit_selpars),nrow = data$n_selblocks, ncol = data$n_ages+6)
      for(i in ind){
        age_pars <- par$logit_selpars[i,1:data$n_ages]  
        age_map <- lsp_map[i,1:data$n_ages]
        ages_omit[[i]] <- which(age_pars == -Inf & is.na(age_map)) #if length then observations for these ages will be omitted.
        if(length(ages_omit[[i]])){
          if(i %in% data$selblock_pointer_indices) for(j in 1:data$n_indices){ #some indices affected
            bad_year = which(data$selblock_pointer_indices[,j] == i & data$use_index_paa[,j] == 1)
            if(length(bad_year)){
              input$log$osa_obs <- c(input$log$osa_obs, paste0(
                "Selectivity block ", i, " has parameters for ages ", paste0(input$ages.lab[ages_omit[[i]]], collapse = ", "), 
                " fixed at 0 and is used by index ", j, "\n",
                " in years ", paste0(input$years[bad_year], collapse = ", "), ".\n", 
                "Observations at these ages will be omitted and remaining observations will be rescaled.\n\n"))
            }
          }
          if(i %in% data$selblock_pointer_fleets) for(j in 1:data$n_fleets){ #some fleets affected
            bad_year = which(data$selblock_pointer_fleets[,j] == i & data$use_catch_paa[,j] == 1)
            if(length(bad_year)){
              input$log$osa_obs <- c(input$log$osa_obs, paste0(
                "Selectivity block ", i, " has parameters for ages ", paste0(input$ages.lab[ages_omit[[i]]], collapse = ", "), 
                " fixed at 0 and is used by fleet ", j, "\n",
                " in years ", paste0(input$years[bad_year], collapse = ", "), ".\n", 
                "Observations at these ages will be omitted and remaining observations will be rescaled.\n\n"))
            }
          }
        }
      }
    }
  }
  
  #step through time adding data
  for(y in 1:data$n_years_model){
    
    # 2. log index catch
    x <- as.data.frame(data$agg_indices)
    x[data$use_indices == 0] <- NA # only include index data to fit in obsvec
    colnames(x) <- paste0("index_", 1:data$n_indices)
    for(i in 1:data$n_indices) {
      if(!is.na(x[y,i])) {
        tmp = data.frame(year = y, fleet = colnames(x)[i], age = NA, type = 'logindex', val = log(x[y,i]))
        obs <- rbind(obs, tmp[, obs.colnames])
      }
    }
    
    # 3. log fleet catch
    x <- as.data.frame(data$agg_catch)
    x[data$use_agg_catch==0] <- NA # can't fit to fleets/years with 0 catch
    colnames(x) <- paste0("fleet_", 1:data$n_fleets)
    for(i in 1:data$n_fleets) {
      if(!is.na(x[y,i])) {
        tmp = data.frame(year = y, fleet = colnames(x)[i], age = NA, type = 'logcatch', val = log(x[y,i]))
        obs <- rbind(obs, tmp[, obs.colnames])
      }
    }
    
    # 4. paa index
    for(i in 1:data$n_indices){
      #obs_levels <- c(obs_levels, paste0("fleet_",i, "_paa"))
      x <- data$index_paa[i,,]
      x[which(data$use_index_paa[,i]==0),] <- NA # only include catch data to fit in obsvec
      indices = paste0("index_", 1:data$n_indices)
      if(data$use_index_paa[y,i] == 1)
      {
        obs_y = x[y,]
        tmp <- ages_omit[[data$selblock_pointer_indices[y,i]]]
        res = transform_paa_obs(obs_y, data$age_comp_model_indices[i], ages_omit = tmp)
        obs_y = res[[1]]
        ind = res[[2]] #now the ages to use is specified for all likelihods by transform_paa_obs
        #multinom, D-M, mvtweedie
        if(data$age_comp_model_indices[i] %in% c(1:2,10,11)) obs_y = obs_y * data$index_Neff[y,i]
        
        #if(data$age_comp_model_indices[i] %in% 3:7) {
        #  ind = res[[2]]
        #} else ind = 1:data$n_ages
        if(length(ind)) {
          tmp = data.frame(year = y, fleet = indices[i], age = (1:data$n_ages)[ind], type = 'indexpaa', val = obs_y[ind])
          #tmp = data.frame(year = y, fleet = indices[i], age = (1:data$n_ages), type = 'indexpaa', val = obs_y)
          obs <- rbind(obs, tmp[, obs.colnames])
        } else {
          data$use_index_paa[y,i] <- 0 #set to not use because there are not enough positive values 
          input$log$osa_obs <- c(input$log$osa_obs, paste0("Setting data$use_index_paa to 0 for index ", i, " and year ", y, "because not enough positive values. \n"))
        }
      }
    }
    
    
    # 5. paa catch
    for(i in 1:data$n_fleets){
      #obs_levels <- c(obs_levels, paste0("fleet_",i, "_paa"))
      x <- data$catch_paa[i,,]
      x[which(data$use_catch_paa[,i]==0),] <- NA # only include catch data to fit in obsvec
      #x = as.data.frame(x)
      fleets = paste0("fleet_", 1:data$n_fleets)
      if(data$use_catch_paa[y,i] == 1) {
        obs_y = x[y,]
        tmp <- ages_omit[[data$selblock_pointer_fleets[y,i]]]
        #multinom, D-M, mvtweedie
        res = transform_paa_obs(obs_y, data$age_comp_model_fleets[i], ages_omit = tmp)
        obs_y = res[[1]]
        ind = res[[2]] #now the ages to use is specified for all likelihods by transform_paa_obs
        if(data$age_comp_model_fleets[i] %in% c(1:2,10,11)) obs_y = obs_y * data$catch_Neff[y,i]
        #if(data$age_comp_model_fleets[i] %in% 3:7) {
        #  ind = res[[2]]
        #} else ind = 1:data$n_ages
        
        if(length(ind)) {
          tmp = data.frame(year = y, fleet = fleets[i], age = (1:data$n_ages)[ind], type = 'catchpaa', val = obs_y[ind])
          obs <- rbind(obs, tmp[, obs.colnames])
        } else {
          data$use_catch_paa[y,i] <- 0 #set to not use because there are not enough positive values
          input$log$osa_obs <- c(input$log$osa_obs, paste0("Setting data$use_catch_paa to 0 for fleet ", i, " and year ", y, "because not enough positive values. \n"))
        }
      }
    }
    
  }
  
  # calculate obsvec indices in keep arrays
  obs$ind <- 1:dim(obs)[1]
  
  
  data$keep_E <- matrix(NA, nrow=data$n_years_Ecov, ncol=data$n_Ecov)
  for(i in 1:data$n_Ecov){
    ind <- which(data$Ecov_use_obs[,i]==1) 
    data$keep_E[ind,i] <- subset(obs, type=='Ecov' & fleet == paste0("Ecov_",i))$ind
  }
  # subtract 1 bc TMB indexes from 0
  data$keep_E <- data$keep_E - 1
  
  data$keep_C <- matrix(NA, nrow=data$n_years_model, ncol=data$n_fleets)
  for(y in 1:data$n_years_model) for(i in 1:data$n_fleets){
    if(data$use_agg_catch[y,i]==1){ 
      data$keep_C[y,i] <- subset(obs, year == y & type=='logcatch' & fleet == paste0("fleet_",i))$ind
    }
  }
  # subtract 1 bc TMB indexes from 0
  data$keep_C <- data$keep_C - 1
  
  data$keep_I <- matrix(NA, nrow=data$n_years_model, ncol=data$n_indices)
  for(y in 1:data$n_years_model) for(i in 1:data$n_indices){
    if(data$use_indices[y,i]==1){ 
      data$keep_I[y,i] <- subset(obs, year == y & type=='logindex' & fleet == paste0("index_",i))$ind
    }
  }
  # subtract 1 bc TMB indexes from 0
  data$keep_I <- data$keep_I - 1
  
  data$condition_no_osa = NULL #to condition on age comps for likelihoods we don't have osa residuals set up for.
  data$subset_discrete_osa = NULL #age comp obs for likelihoods we need to specify as discrete obs.
  
  data$keep_Cpaa <- array(NA, dim=c(data$n_fleets, data$n_years_model, 2))
  for(i in 1:data$n_fleets) {
    for(y in 1:data$n_years_model) if(data$use_catch_paa[y,i]==1){
      tmp = subset(obs, year == y & type=='catchpaa' & fleet==paste0("fleet_",i))
      if(length(tmp$ind)) #should always be TRUE because use_paa changed above
      {
        data$keep_Cpaa[i,y,1:2] <- c(tmp$ind[1], length(tmp$ind))
        #if(data$age_comp_model_fleets[i] %in% 1:2) data$subset_discrete_osa = c(data$subset_discrete_osa, tmp$ind)
        #subset for oneStepPredict can't include these
        if(data$age_comp_model_fleets[i] %in% 8:10) data$condition_no_osa = c(data$condition_no_osa, tmp$ind)
      }
    }
  }
  # subtract 1 bc TMB indexes from 0
  data$keep_Cpaa[,,1] <- data$keep_Cpaa[,,1] - 1
  
  data$keep_Ipaa <- array(NA, dim=c(data$n_indices, data$n_years_model, 2))
  for(i in 1:data$n_indices) {
    for(y in 1:data$n_years_model) if(data$use_index_paa[y,i]==1){
      tmp = subset(obs, year == y & type=='indexpaa' & fleet==paste0("index_",i))
      if(length(tmp$ind)) #should always be TRUE because use_paa changed above
      {
        data$keep_Ipaa[i,y,1:2] <- c(tmp$ind[1], length(tmp$ind))
        #if(data$age_comp_model_indices[i] %in% 1:2) data$subset_discrete_osa = c(data$subset_discrete_osa, tmp$ind)
        #subset for oneStepPredict can't include these
        if(data$age_comp_model_indices[i] %in% 8:10) data$condition_no_osa = c(data$condition_no_osa, tmp$ind)
      }
    }
  }
  # subtract 1 bc TMB indexes from 0
  data$keep_Ipaa[,,1] <- data$keep_Ipaa[,,1] - 1
  
  obs$cohort = as.numeric(obs$year) - as.numeric(obs$age)   #make cohort column. could be useful for analyzing age comp OSA residuals
  data$obs <- obs
  data$obsvec <- obs$val
  data$agesvec <- obs$age #potentially needed for AR1 sigma correlation of logistic-normal paa obs. 
  data$do_osa = 0 #this will be changed when TMB::oneStepPredict is called by fit_wham
  attr(data, "check.passed") <- NULL #if data have been generated from obj$simulate(complete=T), this can be problematic
  
  #data$do_post_samp = rep(0,5) #this will be changed in fit_wham when a sample of posterior process residuals are to be calculated
  if(length(input$log$osa_obs)) input$log$osa_obs <- c("OSA obs: \n", input$log$osa_obs)
  input$data = data
  return(input)
}

transform_paa_obs = function(x, model, zero.criteria = 1e-15, do_mult = FALSE, ages_omit = integer(0)){
  #transforms paa obs for obsvec and appropriate for OSA residuals once model is fit.
  #remove zeros for dirichlet and logistic-normal
  #remove ages/classes where predicted probability = 0 for all likelihoods (i.e., surveys that only observe a subset of ages).
  #transform logistic-normal obs to MVN obs (analogous to log-catch and log-indices)
  all_models <- c("multinomial","dir-mult","dirichlet-miss0","dirichlet-pool0",
                  "logistic-normal-miss0", "logistic-normal-ar1-miss0", "logistic-normal-pool0",
                  "logistic-normal-01-infl","logistic-normal-01-infl-2par", "mvtweedie", "dir-mult-linear")
  # if model %in% 1:2 do nothing for multinomial and D-m
  is_pred_pos = !(1:length(x) %in% ages_omit)
  x[which(!is_pred_pos)] = 0 #if obs in omitted ages are zero this does nothing
  x =  x/sum(x) #rescale to a reduced set of ages, if obs in omitted ages are zero this does nothing
  if(model>2 & model<10){ #not multinom, D-M, mvtweedie or DM-linear
    is_pos = x> zero.criteria # looking for zeros here will also omit the ages with predicted probability = 0
    pos_ind = which(is_pos)
    npos = sum(is_pos)
    zero_ind = which(!is_pos)
    if(npos>1){ #need at least 2 categories
      if(length(zero_ind)){
        x_pos = x[pos_ind] # 0s missing or pooled
        # if(model %in% c(3,5:6)) x_pos = x[pos_ind] # 0s missing
        # if(model %in% c(4,7)) { # 0s pooled
        #   index = 0
        #   x_pos = rep(NA,npos)
        #   for(a in 1:length(x))
        #   {
        #     x_pos[index] = x_pos[index] + x[a]; #just adding zeros. Not necessary? Only necessary for expected proportions
        #     if((x[a] > zero.criteria) & index < npos) index = index+ 1
        #   }
        # }
        x[zero_ind] <- NA
      } else x_pos = x
      
      #transform to multivariate normal obs for logistic-normal.
      #Necessary because TMB::oneStepPredict needs observation transformations on R side.
      #Conditional distribution format for logistic-normal would be needed to avoid this. 
      if(model %in% 5:7){
        y = log(x_pos[-length(x_pos)])
        if(do_mult){ #multiplicative
          for(i in 1:length(y)) {
            y[i] = y[i] - log(1-sum(x_pos[1:i]))
          }
        } else { #additive
          y = y - log(x_pos[length(x_pos)])
        }
        y = c(y,NA) #NA for last category which = 1 - sum of the other categories
        x_pos = y
      }
      x[pos_ind] = x_pos
    } else { #only 1 positive category
      x[] = NA
    }
  } else { #multinom, D-M, mvtweedie
    x[which(!is_pred_pos)] = NA
    pos_ind = which(is_pred_pos)
  }
  return(list(x, pos_ind))
}


set_proj = function(input, proj.opts = NULL)
{
  data = input$data
  if(is.null(proj.opts))
  {
    input$par$logR_proj <- matrix(0, data$n_stocks, 1) # will be set by prepare_projection if SCAA
    input$map$logR_proj <- factor(rep(NA, data$n_stocks))
    
    data$do_proj <- 0
    data$n_years_proj <- 0
    #data$n_years_proj_Ecov <- 0
    data$avg_years_Ecov <- data$n_years_model - (5:1) # c++ indices start at 0
    data$proj_F_opt <- 1
    data$proj_Fcatch <- 0
    data$proj_M_opt <- 1
    data$proj_R_opt <- 1
    data$proj_L_opt <- 1
    data$proj_mu_opt <- 1
    data$mature_proj <- array(0,dim =c(1,1,1))
    data$waa_proj <- array(0,dim =c(1,1,1))
    if(!is.null(data$n_Ecov)){ #sometimes set_proj is called within other set_X functions before set_ecov has been called.
      data$proj_Ecov_opt <- rep(1, data$n_Ecov)
      data$Ecov_use_proj <- matrix(0, 1, data$n_Ecov)
    }
    data$logR_mean <- 0 # only used for SCAA projections
    data$logR_sd <- 0 # only used for SCAA projections
    data$FXSPR_init = rep(0.5, data$n_years_model + data$n_years_proj)
    data$FMSY_init = rep(0.5, data$n_years_model + data$n_years_proj)
    data$F_proj_init = 0.1
    data$percentFMSY = 100
    data$percentFXSPR = 100
  }
  else {
    
    #do nothing?
    #prepare_projection requires a model object returned by fit_wham
    
  }
  input$data = data
  input$options$proj <- proj.opts
  return(input)
}

#' Specify model and parameter configuration for catchability
#'
#' @param input list containing data, parameters, map, and random elements (output from \code{\link{prepare_wham_input}})
#' @param catchability (optional) list specifying options for numbers-at-age random effects, initial parameter values, and recruitment model (see details)
#' 
#' \code{catchability} specifies options for catchability. If \code{NULL} and \code{asap3} is not NULL, a single catchability parameter for each index is used with initial values specified in ASAP file. If both are NULL, initial catchabilities for all indices = 0.3.
#' Otherwise, it is a list with the following optional entries:
#'   \describe{
#'     \item{$re}{Time-varying (random effects) for each index. Vector with length = number of indices.
#'                  Each entry of \code{catchability$re} must be one of the following options:
#'                  \describe{
#'                    \item{"none"}{(default) are constant}
#'                    \item{"iid"}{vary by year and age/par, but uncorrelated}
#'                    \item{"ar1"}{correlated by year (AR1)}
#'                  }
#'                 }
#'     \item{$initial_q}{Initial catchabilities for each index. vector length = number of indices. Will override values provided in \code{basic_info$q}.
#'        If \code{basic_info$q} and \code{asap3} are not provided, default q values are 0.3.}
#'     \item{$q_lower}{Lower bound for catchabilities for each index. vector length = number of indices. For indices with NULL components default lower values are 0.}
#'     \item{$q_upper}{Upper bound for catchabilities for each index. vector length = number of indices. For indices with NULL components default lower values are 1000.}
#'     \item{$prior_sd}{vector of NA and standard deviations to use for gaussian prior on logit transform of catchability parameter. Length = number of indices.
#'       Indices with non-NA values will have mean logit q as a random effect with mean determined by logit transform of \code{catchability$initial_q} and
#'       sigma as standard error.}
#'     \item{$sigma_val}{Vector of initial standard deviation values to use for annual random effects for each index. Values are not used if \code{q$re} = "none". Otherwise, a single value for all indices.}
#'     \item{$sigma_map}{Specify which sigma parameters to fix for the random effect deviations. Must be a vector with length = number of indices. 
#'                Use \code{NA} to fix a parameter and integers to estimate. Use the same integer for multiple indices to share the same sigma parameter.
#'                Not used if \code{re = 'none'} for all indices.}
#'     \item{$cor_vals}{Vector of initial correlation values to use for annual random effects for each index. If unspecified all initial values are 0. Only used if \code{q$re} = "ar1"}
#'     \item{$cor_map}{Specify which ar1 correlation parameters to fix for the random effect deviations. Must be a vector with length = number of indices. 
#'                If \code{re = 'ar1'}, each element (index) must be a single value. Use \code{NA} to fix a parameter and integers to estimate. 
#'                Use the same integer for multiple indices to share the same correlation parameter. Not used if \code{re = 'none'} or \code{re = 'iid'} for all indices.}
#'   }
#'
#' @return a named list with same elements as the input provided with catchability options modified.
#'
#' @seealso \code{\link{prepare_wham_input}} 
#'
#' @examples
#' \dontrun{
#' wham.dir <- find.package("wham")
#' path_to_examples <- system.file("extdata", package="wham")
#' asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))
#' input <- prepare_wham_input(asap3)
#' catchability = list(re = c("iid", "none"))
#' input <- set_q(input, catchability = catchability) #independent time-varying random effects on q for first survey.
#' }
#'
#' @export
set_q = function(input, catchability = NULL){
  
  q_opts = catchability
  data = input$data
  par = input$par
  map = input$map
  
  #clear any map definitions that may exist. necessary because some configurations may not define map elements.
  map <- map[(!names(map) %in% c("q_prior_re", "q_re", "q_repars"))]
  input$log$q <- list()
  asap3 = input$asap3
  #if(is.null(input$asap3)) asap3 = NULL
  #else asap3 = input$asap3
  
  data$q_lower <- rep(0,data$n_indices)
  data$q_upper <- rep(1000,data$n_indices)
  if(!is.null(q_opts$q_lower)) {
    if(length(q_opts$q_lower) != data$n_indices) stop("the length of catchability$q_lower is not equal to the number of indices")
    data$q_lower = q_opts$q_lower
  }
  if(!is.null(q_opts$q_upper)) {
    if(length(q_opts$q_upper) != data$n_indices) stop("the length of catchability$q_upper is not equal to the number of indices")
    data$q_upper = q_opts$q_upper
  }
  par$logit_q = wham:::gen.logit(0.3, data$q_lower, data$q_upper)
  if(!is.null(asap3)) {
    k = 0
    for(i in 1:length(asap3)){
      ind = which(asap3[[i]]$use_index ==1)
      k_ind = k + 1:length(ind)
      par$logit_q[k_ind] = wham:::gen.logit(asap3[[i]]$q_ini[ind], data$q_lower[k_ind], data$q_upper[k_ind]) # use q_ini values from asap3 file
      k = max(k_ind)
    }
  }
  if(!is.null(q_opts$initial_q)) {
    if(length(q_opts$initial_q) != data$n_indices) stop("the length of catchability$initial_q is not equal to the number of indices")
    par$logit_q = gen.logit(q_opts$initial_q, data$q_lower, data$q_upper)
  }
  
  data$use_q_re = rep(0, data$n_indices)
  par$q_re = matrix(0, data$n_years_model, data$n_indices)
  map$q_re = matrix(NA, data$n_years_model, data$n_indices)
  par$q_repars = matrix(0, data$n_indices, 2)
  map$q_repars = matrix(NA, data$n_indices, 2)
  
  if(!is.null(q_opts$re)){
    ind = which(q_opts$re %in% c("iid", "ar1"))
    map$q_re[,ind] = 1:(length(ind) * (data$n_years_model)) #turn on appropriate columns of q_re
    if(!is.null(q_opts$sigma_val)){
      if(length(q_opts$sigma_val) != data$n_indices) stop("the length of catchability$sigma_val provided is not equal to the number of indices")
      if(q_opts$sigma_val[ind]<0) stop("catchability$sigma_val must be greater than 0")
      par$q_repars[ind,1] = log(q_opts$sigma_val[ind])
    }
    iids = q_opts$re == "iid"
    ar1s = q_opts$re == "ar1"
    if(!is.null(q_opts$cor_val)){
      if(length(q_opts$cor_val) != data$n_indices) stop("the length of catchability$cor_val provided is not equal to the number of indices")
      if(abs(q_opts$cor_val[ind])>1) stop("it must be that -1 < catchability$cor_val < 1 ")
      if(any(iids & abs(q_opts$cor_val)>1e-10)) input$log$q <- c(input$log$q, "certain indices have re='iid' and cor_val not = 0. Those will values will be ignored. \n")
      par$q_repars[ind,2] = gen.logit(q_opts$cor_val[ind], -1, 1, 1) 
      par$q_repars[which(iids),2] = 0 #iids must be set to 0.
      
    }
    if(!is.null(q_opts$sigma_map)){
      if(length(q_opts$sigma_map) != data$n_indices) stop("catchability$sigma_map must be a vector of length = number of indices")
      map$q_repars[which(!is.na(q_opts$sigma_map)),1] = 1:sum(!is.na(q_opts$sigma_map))
    } else map$q_repars[ind,1] = 1:length(ind)
    
    if(!is.null(q_opts$cor_map)){
      if(length(q_opts$cor_map) != data$n_indices) stop("catchability$cor_map must be a list of length = number of indices")
      if(any(iids & !is.na(q_opts$cor_map))) stop("certain indices have re='iid' and cor_map is not NA.")
      map$q_repars[which(!is.na(q_opts$cor_map)),2] = sum(!is.na(map$q_repars[,1])) + 1:sum(!is.na(q_opts$cor_map))
      #"iid" not necessary because already set to NA
    } else {
      if(sum(ar1s)) map$q_repars[ind,2] = sum(!is.na(map$q_repars[,1])) + 1:sum(ar1s)	  	
    }
    data$use_q_re[ind] = 1
  }
  map$q_repars = factor(map$q_repars)
  map$q_re = factor(map$q_re)
  
  data$use_q_prior = rep(0, data$n_indices)
  par$q_prior_re = rep(0, data$n_indices)
  map$q_prior_re = rep(NA, data$n_indices)
  data$logit_q_prior_sigma = rep(1, data$n_indices)
  map$logit_q = 1:data$n_indices
  if(!is.null(q_opts$prior)){
    ind = which(!is.na(q_opts$prior_sd))
    data$use_q_prior[ind] = 1
    data$logit_q_prior_sigma[ind] = q_opts$prior_sd[ind]
    par$q_prior_re[ind] = par$logit_q[ind]
    map$q_prior_re[ind] = 1:length(ind)
    map$logit_q[ind] = NA #turn off logit q (mean of the prior) when random effects are used
  }
  map$q_prior_re = factor(map$q_prior_re)
  map$logit_q = factor(map$logit_q)
  
  input$data = data
  input$par = par
  input$map = map
  if(length(input$log$q))	input$log$q <- c("Catchability: \n", input$log$q)
  
  #set any parameters as random effects
  input$random = NULL
  input = set_random(input)
  input$options$q <- q_opts
  return(input)
}

set_random <- function(input){
  random = NULL
  #if(input$data$Ecov_obs_sigma_opt == 4) random = "Ecov_obs_logsigma"
  if(any(input$data$Ecov_obs_sigma_opt==4)) random = c(random, "Ecov_obs_logsigma_re")
  if(any(input$data$selblock_models_re > 1)) random = c(random, "selpars_re")
  if(any(input$data$M_re_model > 1)) random = c(random, "M_re")
  if(any(input$data$Ecov_model > 0)) random = c(random, "Ecov_re")
  if(any(input$data$NAA_re_model > 0)) random = c(random, "log_NAA")
  if(any(input$data$N1_model == 2)) random = c(random, "log_N1")
  if(any(input$data$use_mu_prior > 0)) random = c(random, "mu_prior_re")
  #print("here")
  #print(input$data$use_mu_prior)
  #print(input$data$mu_model)
  if(any(input$data$mu_model %in% c(2:4,6:8,10:12,14:16))) random = c(random, "mu_re")
  if(sum(input$data$use_q_prior)) random = c(random, "q_prior_re")
  if(sum(input$data$use_q_re)) random = c(random, "q_re")
  input$random = random	
  return(input)
}

#' Specify model and parameter configuration for selectivity
#'
#' @param input list containing data, parameters, map, and random elements (output from \code{\link{prepare_wham_input}})
#' @param selectivity (optional) list specifying options for selectivity blocks, models, initial parameter values, parameter fixing/mapping, and random effects (see details)
#'
#' \code{set_selectivity} specifies options for selectivity and allows you to overwrite existing options 
#' in the \code{input} list or as specified in the ASAP data file. If \code{selectivity = NULL}, selectivity options from 
#' \code{input} are used. 
#'
#' \code{\link{prepare_wham_input}(..., selectivity=selectivity)} calls \code{set_selectivity(..., selectivity=selectivity)}.
#' If you already have created \code{input} with \code{prepare_wham_input}, you can also use \code{set_selectivity(input, selectivity=selectivity)}
#' to modify the selectivity specification.
#'
#' \code{selectivity} is a list with the following entries:
#'   \describe{
#'     \item{$model}{Selectivity model for each block. Vector with length = number of selectivity blocks. Each entry must be one of: "age-specific", "logistic", "double-logistic", or "decreasing-logistic".}
#'     \item{$re}{Time-varying (random effects) for each block. Vector with length = number of selectivity blocks.
#'                  If \code{NULL}, selectivity parameters in all blocks are constant over time and uncorrelated.
#'                  Each entry of \code{selectivity$re} must be one of the following options, where the selectivity parameters are:
#'                  \describe{
#'                    \item{"none"}{(default) are constant and uncorrelated}
#'                    \item{"iid"}{vary by year and age/par, but uncorrelated}
#'                    \item{"ar1"}{correlated by age/par (AR1), but not year}
#'                    \item{"ar1_y"}{correlated by year (AR1), but not age/par}
#'                    \item{"2dar1"}{correlated by year and age/par (2D AR1)}
#'                  }
#'                 }
#'     \item{$initial_pars}{Initial parameter values for each block. List of length = number of selectivity blocks. Each entry must be
#'                a vector of length # parameters in the block, i.e. \code{c(2,0.2)} for logistic (a50 and 1/slope) or \code{c(0.5,0.5,0.5,1,1,0.5)} for
#'                age-specific parameters when there are 6 ages. Default is to set at middle of parameter range. This is 0.5 for age-specific and n.ages/2 
#'                or logistic, double-logistic, and decreasing-logistic.}
#'     \item{$fix_pars}{Alternative to \code{$map_pars} for specifying which selectivity parameters (only fixed effects) to fix at initial values. 
#'                List of length = number of selectivity blocks. E.g. model with 3 age-specific blocks and 6 ages, 
#'                \code{list(4:5, 4, 2:4))} will fix ages 4 and 5 in block 1, age 4 in block 2, and ages 2, 3, and 4 in block 3.
#'                Use NULL to not fix any parameters for a block, e.g. list(NULL, 4, 2) does not fix any pars in block 1.}
#'     \item{$map_pars}{Alternative to \code{$fix_pars} for specifying how to fix selectivity parameters (only fixed effects), corresponds 
#'                to \code{map$logit_selpars}. List of length = number of selectivity blocks, where each item is a 
#'                vector of length = number of selectivity parameters (age-specific: n.ages, logistic: 2, 
#'                double-logistic: 4). Use \code{NA} to fix a parameter and integers to estimate. Use the same integer
#'                for multiple ages or fleets/indices to estimate a shared parameter. E.g. for a model with 3 age-specific 
#'                blocks (1 fleet, 2 indices) and 6 ages, \code{$map_pars = list(c(1,2,3,NA,NA,4), c(5,6,7,NA,8,8), c(1,2,3,NA,NA,4))}
#'                will estimate ages 1-3 and 6 in block 1 (fleet), ages 1-3 and 4-5 (shared) in block 2 (index 1), and then set the
#'                index 2 (block 3) selectivity equal to the fleet.}
#'     \item{$initial_sigma}{Initial standard deviation values to use for the random effect deviations. Must be a vector with length = number of blocks. 
#'                Use natural (not log) scale, must be positive. \code{par$sel_repars[,1]} will be estimated on log-scale. Not used if \code{re = 'none'} for all blocks.}
#'     \item{$map_sigma}{Specify which SD parameters to fix for the random effect deviations. Must be a vector with length = number of blocks. 
#'                Use \code{NA} to fix a parameter and integers to estimate. Use the same integer for multiple blocks to estimate a shared SD parameter.
#'                Not used if \code{re = 'none'} for all blocks.}
#'     \item{$initial_cor}{Initial correlation values to use for the random effect deviations. Must be a list with length = number of blocks. If \code{re = 'ar1'} or \code{ar1_y'}
#'                for a block, list element must be a single value. If \code{re = '2dar1'} for a block, list element must be a vector of length 2 (first is for "age", second is for "year") for each block. 
#'                Use natural scale, must be between -1 and 1. \code{par$sel_repars[,2:3]} will be estimated on a transform scale,  (2 / (1 + exp(-2x))) - 1. 
#'                Not used if \code{re = 'none'} or \code{re = 'iid'} for all blocks.}
#'     \item{$map_cor}{Specify which correlation parameters to fix for the random effect deviations. Must be a list with length = number of blocks. 
#'                If \code{re = 'ar1'} or \code{ar1_y'}, each list element (block) must be a single value. If \code{re = '2dar1'}, 
#'                must be a vector of length 2 (first is for "age", second is for "year"). Use \code{NA} to fix a parameter and integers to estimate. 
#'                Use the same integer for multiple blocks to estimate a shared correlation parameter. Not used if \code{re = 'none'} or \code{re = 'iid'} for all blocks.}
#'     \item{$n_selblocks}{How many selectivity blocks. Optional. If unspecified and no asap3 object, then this is set to the number 
#'                of fleets + indices. If specified, ensure other components of \code{selectivity} are consistent.}
#'   }
#'
#' @return a named list with same elements as the input provided with selectivity options modified.
#'
#' @seealso \code{\link{prepare_wham_input}} 
#'
#' @examples
#' \dontrun{
#' wham.dir <- find.package("wham")
#' path_to_examples <- system.file("extdata", package="wham")
#' asap3 <- read_asap3_dat(file.path(path_to_examples,"ex1_SNEMAYT.dat"))
#' input <- prepare_wham_input(asap3, NAA_re = list(sigma = "rec"))
#' sel <- list(model=rep("logistic",input$data$n_selblocks),
#'    initial_pars=rep(list(c(3,3)),input$data$n_selblocks),
#'    fix_pars=rep(list(NULL),input$data$n_selblocks),
#' input <- set_selectivity(input, selectivity = sel) #logistic selectivity for all selectivity blocks
#' }
#'
#' @export
set_selectivity = function(input, selectivity)
{
  data = input$data
  par = input$par
  map = input$map
  
  par_index = list(
    1:data$n_ages,
    data$n_ages + 1:2,
    data$n_ages + 3:6,
    data$n_ages + 1:2
  )
  asap3 = input$asap3
  input$log$selectivity <- list()
  if(is.null(asap3)) {
    if(is.null(selectivity$n_selblocks)) {
      #data$n_selblocks <- data$n_fleets+ data$n_indices
      #selblock pointers are defined upstream in set_indices and set_catch
      data$n_selblocks = max(input$data$selblock_pointer_fleets,input$data$selblock_pointer_indices)
      input$log$selectivity <- c(input$log$selectivity, paste0("number of selblocks, ", data$n_selblocks, 
                                                               ", is being determined by max(input$data$selblock_pointer_fleets,input$data$selblock_pointer_indices).\n"))
      #data$n_selblocks = data$n_fleets + data$n_indices  #1 for fleet, 1 for index
    } else data$n_selblocks = selectivity$n_selblocks
  } else {
    data$n_selblocks = 0
    data$selblock_models = integer(0)
    for(i in 1:length(asap3)){ #have to do all fleets first, then indices
      data$n_selblocks <- data$n_selblocks + asap3[[i]]$n_fleet_sel_blocks
      data$selblock_models <- c(data$selblock_models, asap3[[i]]$sel_block_option)
    }
    for(i in 1:length(asap3)){
      which_indices <- which(asap3[[i]]$use_index ==1) # length = data$n_indices
      asap3[[i]]$index_sel_option <- asap3[[i]]$index_sel_option[which_indices]
      asap3[[i]]$index_sel_ini = asap3[[i]]$index_sel_ini[which_indices]
      data$n_selblocks <- data$n_selblocks + length(which_indices)
      data$selblock_models <- c(data$selblock_models, asap3[[i]]$index_sel_option)
    }
    orig_sel_models <- data$selblock_models
  }
  no_asap = is.null(asap3)
  selopts <- c("age-specific","logistic","double-logistic","decreasing-logistic")
  
  if(is.null(selectivity$model)) {
    if(no_asap) data$selblock_models <- rep(2, data$n_selblocks)
  } 
  if(!is.null(selectivity$model)){
    if(length(selectivity$model) != data$n_selblocks) stop("Length of selectivity$model must equal number of selectivity blocks (e.g., asap3$n_fleet_sel_blocks + asap3$n_indices)")
    if(!all(selectivity$model %in% selopts)) stop("Each model entry must be one of the following: 'age-specific','logistic','double-logistic','decreasing-logistic'")
    data$selblock_models <- match(selectivity$model, selopts)
  }
  
  if(is.null(selectivity$re)) data$selblock_models_re <- rep(1, data$n_selblocks) # default: no RE on selectivity parameters
  if(!is.null(selectivity$re)){
    if(length(selectivity$re) != data$n_selblocks) stop("Length of selectivity$re must equal number of selectivity blocks (asap3$n_fleet_sel_blocks + asap3$n_indices)")
    if(!all(selectivity$re %in% c("none","iid","ar1","ar1_y","2dar1"))) stop("Each selectivity$re entry must be one of the following: 'none','iid','ar1','ar1_y','2dar1'")
    data$selblock_models_re <- match(selectivity$re, c("none","iid","ar1","ar1_y","2dar1"))
  }
  
  selblock_pointers <- cbind(data$selblock_pointer_fleets, data$selblock_pointer_indices)
  data$selblock_years <- matrix(0, nrow=data$n_years_model, ncol=data$n_selblocks)
  for(b in 1:data$n_selblocks) data$selblock_years[,b] <- apply(selblock_pointers, 1, function(x) b %in% x)
  data$n_years_selblocks <- apply(data$selblock_years, 2, sum)
  
  data$n_selpars <- c(data$n_ages,2,4,2)[data$selblock_models] # num selpars per block
  # Prep selectivity initial values  
  selpars_ini = matrix(NA, data$n_selblocks, data$n_ages + 6)
  # Prep selectivity map
  phase_selpars = matrix(-1, data$n_selblocks, data$n_ages + 6)
  for(b in 1:data$n_selblocks){
    phase_selpars[b,par_index[[data$selblock_models[b]]]] = 1
  }
  
  # initial values
  if(is.null(selectivity$initial_pars)) {
    if(!no_asap) {
      j = 1
      for(k in 1:length(asap3)) for(i in 1:asap3[[k]]$n_fleet_sel_blocks) {
        selpars_ini[j,] = asap3[[k]]$sel_ini[[i]][,1]
        j = j + 1
      }
      for(k in 1:length(asap3)) for(i in 1:length(asap3[[k]]$index_sel_ini)){
        selpars_ini[j,] = asap3[[k]]$index_sel_ini[[i]][,1]
        j = j + 1
      }
    }
    default_selpars <- list()
    dpars = c(0.5,data$n_ages/2)
    orig_selpars <- list()
    
    for(b in 1:data$n_selblocks){
      if(data$selblock_models[b] == 1) 
      {
        default_selpars[[b]] <- rep(0.5, data$n_ages) # default to middle of par range
      }
      if(data$selblock_models[b] %in% c(2,4)) {
        default_selpars[[b]] <- rep(data$n_ages/2, 2) # default to middle of par range
      }
      if(data$selblock_models[b] == 3) {
        default_selpars[[b]] <- rep(data$n_ages/2, 4) # default to middle of par range
      }
      if(!no_asap){
        orig_selpars[[b]] <- selpars_ini[b,par_index[[data$selblock_models[b]]]]
      }
    }
    if(no_asap) for(b in 1:data$n_selblocks){
      selpars_ini[b,] <- c(rep(0.5,data$n_ages), rep(data$n_ages/2, 6))#default_selpars[[b]] # default to middle of par range
    }
    if(!no_asap) {
      #defined above
      #orig_sel_models <- c(asap3[[i]]$sel_block_option, asap3[[i]]$index_sel_option)
      sel_mod_diff_warn <- NULL
      #for(b in 1:data$n_selblocks){
      for(b in 1:length(orig_sel_models)){
        if(data$selblock_models[b] != orig_sel_models[b]){
          selpars_ini[b,par_index[[data$selblock_models[b]]]] <- default_selpars[[b]] # default to middle of par range
          sel_mod_diff_warn <- paste(sel_mod_diff_warn, paste0("For Selectivity Block ",b,":"), 
                                     #paste0("  ASAP .dat file: ",selopts[orig_sel_models[b]]),
                                     paste0("  selectivity$model: ",selopts[data$selblock_model[b]]),
                                     paste0("  Changing values from ASAP3 .dat file ",
                                            #paste0(orig_selpars[[b]], collapse = ", "),
                                            " to inits in middle of par range ",
                                            paste0(default_selpars[[b]], collapse = ", ")), sep="\n")
        }
      }
      if(!is.null(sel_mod_diff_warn)){
        sel_mod_diff_warn <- paste("Selectivity models differ from ASAP .dat file but initial","parameter values not specified. Please check initial values","and specify with selectivity$initial_pars if desired.",sel_mod_diff_warn,sep="\n")
        input$log$selectivity <- c(input$log$selectivity, sel_mod_diff_warn, sep="\n")
      }
    }
  } else {
    if(!is.list(selectivity$initial_pars)) stop("selectivity$initial_pars must be a list")
    if(length(selectivity$initial_pars) != data$n_selblocks) stop("Length of selectivity$initial_pars must equal number of selectivity blocks (asap3$n_fleet_sel_blocks + asap3$n_indices)")
    for(b in 1:data$n_selblocks){
      if(length(selectivity$initial_pars[[b]]) != data$n_selpars[b]) stop(paste0("Length of vector ",b," in the selectivity$initial_pars list is not equal to the number of selectivity parameters for block ",b,": ",data$n_selpars[b]))
      selpars_ini[b,par_index[[data$selblock_models[b]]]] = selectivity$initial_pars[[b]]
    }
  }
  
  # which selpars to fix, either $fix_pars or $map_pars
  if(is.null(selectivity$fix_pars) & is.null(selectivity$map_pars)){
    if(!no_asap){
      j = 1
      for(k in 1:length(asap3)) for(i in 1:asap3[[k]]$n_fleet_sel_blocks) {
        phase_selpars[j,par_index[[asap3[[k]]$sel_block_option[i]]]] = asap3[[k]]$sel_ini[[i]][par_index[[asap3[[k]]$sel_block_option[i]]],2]
        j <- j + 1
      }
      for(k in 1:length(asap3)) for(i in 1:length(asap3[[k]]$index_sel_ini)) {
        phase_selpars[j,par_index[[asap3[[k]]$index_sel_option[i]]]] = asap3[[k]]$index_sel_ini[[i]][par_index[[asap3[[k]]$index_sel_option[i]]],2]
        j <- j + 1
      }
    } 
  } else {
    if(!is.null(selectivity$fix_pars) & !is.null(selectivity$map_pars)) stop("Cannot specify $fix_pars and $map_pars (both set which pars to estimate). Choose one or the other.")
    # fix_pars
    if(!is.null(selectivity$fix_pars)){
      if(!is.list(selectivity$fix_pars)) stop("selectivity$fix_pars must be a list")
      if(length(selectivity$fix_pars) != data$n_selblocks) stop("Length of selectivity$fix_pars must equal number of selectivity blocks (asap3$n_fleet_sel_blocks + asap3$n_indices).
        Use 'NULL' to not fix any parameters for a block, e.g. list(NULL,4,2) does not fix any pars in block 1")
      for(b in 1:data$n_selblocks){
        if(data$selblock_models[b] == 1) phase_selpars[b,selectivity$fix_pars[[b]]] = -1
        if(data$selblock_models[b] %in% c(2,4)) phase_selpars[b,data$n_ages+selectivity$fix_pars[[b]]] = -1
        if(data$selblock_models[b] == 3) phase_selpars[b,data$n_ages+2+selectivity$fix_pars[[b]]] = -1
      }
    }
    # map_pars
    if(!is.null(selectivity$map_pars)){
      if(!is.list(selectivity$map_pars)) stop("selectivity$map_pars must be a list")
      if(length(selectivity$map_pars) != data$n_selblocks) stop("Length of selectivity$map_pars must equal number of selectivity blocks (asap3$n_fleet_sel_blocks + asap3$n_indices)")
      for(b in 1:data$n_selblocks){
        naind <- which(is.na(selectivity$map_pars[[b]]))
        if(data$selblock_models[b] == 1) phase_selpars[b, naind] = -1
        if(data$selblock_models[b] %in% c(2,4)) phase_selpars[b, data$n_ages + naind] = -1
        if(data$selblock_models[b] == 3) phase_selpars[b,data$n_ages + 2 + naind] = -1
      }
    }
  }
  # For age-specific selectivity blocks, check for ages with ~zero catch and fix these at 0
  age_specific <- which(data$selblock_models==1)
  for(b in age_specific){
    if(all(phase_selpars[b,] < 0)){ # if no selpars estimated, keep fixed at specified initial values
      phase_selpars[b,] = -1
    } else {
      ind = list(fleets = which(apply(data$selblock_pointer_fleets == b,2,sum) > 0))
      ind$indices = which(apply(data$selblock_pointer_indices == b,2,sum) > 0)
      paa = matrix(nrow = 0, ncol = data$n_ages)
      if(length(ind$fleets)) for(f in ind$fleets) {
        y = data$catch_paa[f,which(data$selblock_pointer_fleets[,f] == b & data$use_catch_paa[,f] == 1),]
        paa = rbind(paa,y)
      }
      if(length(ind$indices)) for(i in ind$indices) {
        y = data$index_paa[i,which(data$selblock_pointer_indices[,i] == b & data$use_index_paa[,i] == 1),]
        paa = rbind(paa,y)
      }
      y = apply(paa,2,sum)
      ind <- which(y < 1e-5 & phase_selpars[b,1:data$n_ages] > 0) #if phase is set to -1 and 0s in all years for an age and selpars_ini >0, then keep selpars as is.
      selpars_ini[b, ind] = 0
      phase_selpars[b, ind] = -1
    }
  }
  temp <- matrix(NA, data$n_selblocks, data$n_ages + 6)
  # if(!is.null(selectivity$fix_pars)){ # use fix_pars
  temp[which(phase_selpars > 0)] = 1:sum(phase_selpars>0)
  # }
  if(!is.null(selectivity$map_pars)){ # use map_pars directly
    for(b in 1:data$n_selblocks) temp[b, par_index[[data$selblock_models[b]]]] = selectivity$map_pars[[b]]
  }
  for(b in 1:data$n_selblocks){ 
    if(data$selblock_models_re[b] == 3){
      bl <- temp[b, par_index[[data$selblock_models[b]]]]
      # print(b)
      # print(data$selblock_models_re[b])
      # print(data$selblock_models[b])
      # print(bl)
      if(sum(!is.na(bl)) < 3) {
        #data$selblock_models_re[b] <- 1 #no RE for this block
        stop(paste0("'ar1' (AR1(age)) selectivity random effects specified for block ",b,", but number of free mean parameters is <= 2, 
        which will not be identifiable. Use age-specific selectivity with more free mean selectivity parameters (par$logit_selpars) or use 
        a different random effects specification."))
      } else {
        # if re='ar1' (by age) with age-specific selectivity, warning that only estimate one mean shared across ages
        # but don't overwrite fixed pars (likely will be fixing one age at 1)
        input$log$selectivity <- c(input$log$selectivity, paste0("\nNOTE: 'ar1' (AR1(age)) with age-specific selectivity for block ",b,
                                                                 ". Only estimating one mean parameter to be shared across ages that are not fixed (there are at least 3 ages that are free here).\n"))
        # but don't overwrite fixed pars (likely will be fixing one age at 1)
        bl[!is.na(bl)] = min(bl, na.rm=TRUE)
        temp[b, par_index[[data$selblock_models[b]]]] = bl
      }
      # warning message if no mean sel pars (logit_selpars) are fixed
      # allow so user can fit model without fixing and then fix the age with highest sel at 1
      if(all(!is.na(bl))  & data$selblock_models[b] == 1) input$log$selectivity <- c(input$log$selectivity, paste0("\nNOTE: 'ar1' (AR1(age)) with age-specific selectivity for block ",b,
                                                                                                                   ", but no age fixed at 1. Advised to fit the current model and then fix the age with highest selectivity at 1. 
        Can use selectivity$fix_pars.\n"))
      # could add warning message if most ages are fixed, leaving less than xx ages to estimate with the AR1 re
    }
  }
  #data$selpars_est <- phase_selpars
  #data$selpars_est[data$selpars_est == -1] = 0
  data$selpars_est <- matrix(0, data$n_selblocks, data$n_ages + 6)
  data$selpars_est[which(!is.na(temp))] = 1
  # print(data$selpars_est)
  data$n_selpars_est <- apply(data$selpars_est > 0, 1, sum)
  map$logit_selpars = factor(temp)
  
  # initial values on logit scale, par$logit_selpars
  selpars_lo = selpars_hi = matrix(0, data$n_selblocks, data$n_ages + 6)
  selpars_hi[,1:data$n_ages] = 1
  selpars_hi[,data$n_ages + 1:6] = data$n_ages
  data$selpars_lower = selpars_lo #only need these for estimated parameters
  data$selpars_upper = selpars_hi
  
  selpars_ini[which(selpars_ini > selpars_hi)] <- selpars_hi[which(selpars_ini > selpars_hi)]
  selpars_ini[which(selpars_ini < selpars_lo)] <- selpars_lo[which(selpars_ini < selpars_lo)]
  par$logit_selpars = log(selpars_ini-selpars_lo) - log(selpars_hi - selpars_ini)
  par$logit_selpars[!is.na(map$logit_selpars) & is.infinite(par$logit_selpars) & par$logit_selpars<0] = -10
  par$logit_selpars[!is.na(map$logit_selpars) & is.infinite(par$logit_selpars) & par$logit_selpars>0] = 10
  
  # random effects, selpars_re
  # number of estimated selpars per block * number of years per block (only if that block has re)
  par$selpars_re <- array(0, c(data$n_selblocks, data$n_years_model, data$n_ages))
  map$selpars_re <- array(NA, c(data$n_selblocks, data$n_years_model, data$n_ages))
  if(any(data$selblock_models_re > 1)){
    #par$selpars_re <- rep(0, sum((data$selblock_models_re > 1)*data$n_selpars_est*data$n_years_selblocks))
    #tmp_vec <- c()
    ct <- 0
    for(b in 1:data$n_selblocks){
      if(data$selblock_models_re[b] > 1){
        #tmp <- matrix(0, nrow=data$n_years_selblocks[b], ncol=data$n_selpars_est[b])
        if(data$selblock_models_re[b] %in% c(2,5)){ # 2d ar1
          map$selpars_re[b,which(data$selblock_years[,b]==1),1:data$n_selpars_est[b]] <- ct + 1:(data$n_years_selblocks[b]*data$n_selpars_est[b])
          #tmp[] = 1:(dim(tmp)[1]*dim(tmp)[2]) + ct # all y,a estimated
        }
        if(data$selblock_models_re[b] == 3){ # ar1_a (devs by age, constant by year)
          #for(i in 1:dim(tmp)[2]) {
          for(i in 1:data$n_selpars_est[b]) {
            map$selpars_re[b,1:data$n_years_selblocks[b],i] <- ct + i
            #tmp[,i] = (i + ct)
          }
        }
        if(data$selblock_models_re[b] == 4){ # ar1_y (devs by year, constant by age)
          #for(i in 1:dim(tmp)[1]) {
          for(i in 1:data$n_years_selblocks[b]) {
            map$selpars_re[b,i,1:data$n_selpars_est[b]] <- ct + i
            #tmp[i,] = (i + ct)
          }
        }
        #if(length(tmp)){
        if(any(!is.na(map$selpars_re))){
          #ct = max(tmp)
          ct <- max(map$selpars_re, na.rm =T)
          #tmp_vec = c(tmp_vec, as.vector(tmp))
        } else stop(paste0("set_selectivity thinks you want to use random effects for selblock ", b, 
                           ", but either the selblock is not used, or there are no selectivity parameters being estimated."))
      }
    }
    #map$selpars_re <- factor(tmp_vec)
  } else {
    #par$selpars_re <- matrix(0)
    #map$selpars_re <- factor(NA)
  }
  map$selpars_re <- factor(map$selpars_re)
  
  # initial and map for parameters controlling selectivity RE
  trans <- function(x) return((2/(1 + exp(-x))) - 1) # transform for correlation par
  # default initial values: sigma = 0.1, rho = 0
  par$sel_repars <- matrix(0, nrow=data$n_selblocks, ncol=3)
  par$sel_repars[,1] <- log(0.1)
  for(b in 1:data$n_selblocks){
    if(data$selblock_models_re[b] == 3) par$sel_repars[b,3] <- 0 # if ar1 over ages only, fix rho_y = 0
    if(data$selblock_models_re[b] == 4) par$sel_repars[b,2] <- 0 # if ar1 over years only, fix rho = 0
    # check if only 1 estimated sel par (e.g. because all but 1 age is fixed), can't estimate rho
    if(data$n_selpars_est[b] < 2) par$sel_repars[b,2] <- 0
  }
  if(!is.null(selectivity$initial_sigma)){
    if(any(selectivity$initial_sigma < 0)) stop('Variance controlling selectivity random effects must be positive.') 
    par$sel_repars[,1] = log(selectivity$initial_sigma) # log scale
  }
  if(!is.null(selectivity$initial_cor)){
    if(any(sapply(selectivity$initial_cor, function(x) any(x < -1)))) stop('Correlation parameters controlling selectivity random effects must be between -1 and 1.')
    if(any(sapply(selectivity$initial_cor, function(x) any(x > 1)))) stop('Correlation parameters controlling selectivity random effects must be between -1 and 1.')
    for(b in 1:data$n_selblocks){
      if(data$selblock_models_re[b] == 3) par$sel_repars[b,2] <- trans(selectivity$initial_cor[[b]]) # if ar1 over ages, use specified initial
      if(data$selblock_models_re[b] == 4) par$sel_repars[b,3] <- trans(selectivity$initial_cor[[b]]) # if ar1 over years, use specified initial
      if(data$selblock_models_re[b] == 5) par$sel_repars[b,2:3] <- trans(selectivity$initial_cor[[b]]) # if 2dar1 over years, use both
    }
  }  
  
  # map selectivity RE
  tmp.sel.repars <- par$sel_repars
  for(b in 1:data$n_selblocks){
    if(data$selblock_models_re[b] == 1) tmp.sel.repars[b,] <- rep(NA,3) # no RE pars to estimate
    if(data$selblock_models_re[b] == 2) tmp.sel.repars[b,2:3] <- rep(NA,2) # estimate sigma
    if(data$selblock_models_re[b] == 3) tmp.sel.repars[b,3] <- NA # estimate sigma, rho
    if(data$selblock_models_re[b] == 4) tmp.sel.repars[b,2] <- NA # estimate sigma, rho_y
    if(data$n_selpars_est[b] < 2) tmp.sel.repars[b,2] <- NA # can't estimate rho if only 1 selpar estimated
  }
  if(!is.null(selectivity$map_sigma)){
    if(length(selectivity$map_sigma) != data$n_selblocks) stop("selectivity$map_sigma must be a vector of length = number of selectivity blocks")
    tmp.sel.repars[is.na(selectivity$map_sigma),1] = NA
  }
  if(!is.null(selectivity$map_cor)){
    if(length(selectivity$map_cor) != data$n_selblocks) stop("selectivity$map_cor must be a list of length = number of selectivity blocks")
    for(b in 1:data$n_selblocks){
      if(data$selblock_models_re[b] == 3) tmp.sel.repars[is.na(selectivity$map_cor[[b]]),2] = NA
      if(data$selblock_models_re[b] == 4) tmp.sel.repars[is.na(selectivity$map_cor[[b]]),3] = NA
      if(data$selblock_models_re[b] == 5) tmp.sel.repars[is.na(selectivity$map_cor[[b]]),2:3] = NA
    }
  }
  ind.notNA <- which(!is.na(tmp.sel.repars))
  tmp.sel.repars[ind.notNA] <- 1:length(ind.notNA)
  if(!is.null(selectivity$map_sigma)){
    tmp.sel.repars[,1] = selectivity$map_sigma
  }
  if(!is.null(selectivity$map_cor)){
    for(b in 1:data$n_selblocks){
      if(data$selblock_models_re[b] == 3) tmp.sel.repars[b,2] = selectivity$map_cor[[b]]
      if(data$selblock_models_re[b] == 4) tmp.sel.repars[b,3] = selectivity$map_cor[[b]]
      if(data$selblock_models_re[b] == 5) tmp.sel.repars[b,2:3] = selectivity$map_cor[[b]]
    }
    st <- ifelse(!all(is.na(tmp.sel.repars[,1])), max(tmp.sel.repars[,1], na.rm=T), 0)
    tmp.sel.repars[,2:3] = tmp.sel.repars[,2:3] + st
  }
  map$sel_repars = factor(tmp.sel.repars)
  
  input$data = data
  input$par = par
  input$map = map
  if(length(input$log$selectivity))	input$log$selectivity <- c("Selectivity: \n", input$log$selectivity)
  
  #set any parameters as random effects
  input$random = NULL
  input = set_random(input)
  input$options$selectivity <- selectivity
  return(input)
  
}

set_WAA <- function(input, waa_info = NULL) {
  data <- input$data
  asap3 <- input$asap3
  # Weight-at-age
  data$waa <- array(NA,dim = c(data$n_fleets + data$n_regions + data$n_indices + data$n_stocks, data$n_years_model, data$n_ages))
  if(!is.null(asap3)) {
    # waa_temp = list()
    data$waa_pointer_fleets <- integer()
    data$waa_pointer_ssb <- integer()
    data$waa_pointer_indices <- integer()
    waa_pointer_totcatch <- integer()
    #fill with fleet catch waa
    i <- 1
    for(k in 1:length(asap3)) {
      x <- asap3[[k]]
      for(f in 1:x$n_fleets){
        data$waa[i,,] <- x$WAA_mats[[x$WAA_pointers[2*f-1]]]
        data$waa_pointer_fleets[i] <- i
        i <- i + 1
      }
    }
    #fill with total catch waa (region)
    for(k in 1:length(asap3)) {
      x <- asap3[[k]]
      data$waa[data$n_fleets + k,,] <- x$WAA_mats[[x$WAA_pointers[2*x$n_fleets+1]]]
      waa_pointer_totcatch[k] <- data$n_fleets + k
    }
    #fill with index waa
    i <- 1
    for(k in 1:length(asap3)) {
      x <- asap3[[k]]
      for(f in 1:x$n_indices){
        data$waa[data$n_fleets + data$n_regions + i,,] <- x$WAA_mats[[x$index_WAA_pointers[f]]]
        data$waa_pointer_indices[i] <- data$n_fleets + data$n_regions + i
        i <- i + 1
      }
    }
    #fill with ssb waa (stocks)
    for(k in 1:length(asap3)) {
      x <- asap3[[k]]
      data$waa[data$n_fleets + data$n_regions + data$n_indices + k,,] <- x$WAA_mats[[x$WAA_pointers[2*x$n_fleets+3]]]
      data$waa_pointer_ssb[k] <- data$n_fleets + data$n_regions + data$n_indices + k
    }
    # n_waa_total = sapply(asap3, function(x){
    # 	i <- c(seq(1,(x$n_fleets+1)*2-1,2),(x$n_fleets+1)*2 + 1:2)
    #  	WAA_pointers <- x$WAA_pointers[i] #wham has no discard data, so remove those WAA matrices
    # 	WAA_pointers <- c(WAA_pointers,x$index_WAA_pointers) #need any for indices too
    # 	length(unique(WAA_pointers))
    # })
    #print(sum(n_waa_total))
    
    # data$waa_pointer_fleets <- integer()
    # data$waa_pointer_ssb <- integer()
    # data$waa_pointer_indices <- integer()
    # data$waa_pointer_totcatch <- integer()
    # data$waa <- array(NA,dim = c(sum(n_waa_total), data$n_years_model, data$n_ages))
    # j = 1
    # n_waa <- 0
    #   for(a in 1:length(asap3)) {
    # 	#fleet 1 catch, fleet 2 catch, ..., fleet n catch, totcatch, ssb
    # 	x <- asap3[[a]]
    # 	i <- c(seq(1,(x$n_fleets+1)*2-1,2),(x$n_fleets+1)*2 + 2)
    #  	WAA_pointers <- c(x$WAA_pointers[i],x$index_WAA_pointers) #wham has no discard data, so remove those WAA matrices
    # 	for(k in unique(WAA_pointers)){ #just retain the needed WAA matrices
    # 		data$waa[j,,] <- x$WAA_mats[[k]]  #note order is changing
    # 		j <- j + 1
    # 	}
    # 	#print(WAA_pointers)
    # 	new_pointer <- n_waa + (1:length(unique(WAA_pointers)))
    # 	#print(new_pointer)
    # 	new_pointer <- new_pointer[match(WAA_pointers,unique(WAA_pointers))]
    # 	#print(new_pointer)
    # 	data$waa_pointer_fleets <- c(data$waa_pointer_fleets,new_pointer[1:x$n_fleets])
    # 	#print(data$waa_pointer_fleets)
    # 	data$waa_pointer_totcatch = c(data$waa_pointer_totcatch,new_pointer[x$n_fleets+1])
    # 	#print(data$waa_pointer_totcatch)
    # 	data$waa_pointer_ssb = c(data$waa_pointer_ssb,new_pointer[x$n_fleets+2])
    # 	data$waa_pointer_indices = c(data$waa_pointer_indices,new_pointer[x$n_fleets+2 + 1:x$n_indices])
    # 	#print(data$waa_pointer_indices)
    # 	n_waa <- n_waa + n_waa_total[a]
    # }
  } else {
    L = 100*(1-exp(-0.3*(1:data$n_ages - 0)))
    W = rep(exp(-11)*L^3, each = data$n_years_model)
    data$waa = array(W, dim = c(1, data$n_years_model, data$n_ages))
    data$waa_pointer_fleets = rep(1,data$n_fleets)
    data$waa_pointer_ssb = 1
    data$waa_pointer_indices = rep(1,data$n_indices)
  }
  input$log$waa <- list()
  if(!is.null(waa_info$waa)){
    data$waa = waa_info$waa
    dim_waa = dim(data$waa)
    if(length(dim_waa) != 3) stop("waa_info$waa must be a 3d array. second index is number of years, third is number of ages.")
    if(is.null(waa_info$waa_pointer_fleets)){
      input$log$waa <- c(input$log$waa, "waa_info$waa is provided without waa_info$waa_pointer_fleets so the first waa matrix will be used for all fleets. \n")
      data$waa_pointer_fleets = rep(1,data$n_fleets)
    }
    if(is.null(waa_info$waa_pointer_indices)){
      input$log$waa <- c(input$log$waa, "waa_info$waa is provided without waa_info$waa_pointer_indices so the first waa matrix will be used. \n")
      data$waa_pointer_indices = rep(1,data$n_indices)
    }
    if(is.null(waa_info$waa_pointer_ssb)){
      input$log$waa <- c(input$log$waa, "waa_info$waa is provided without waa_info$waa_pointer_ssb so the first waa matrix will be used. \n")
      data$waa_pointer_ssb = rep(1,data$n_stocks)
    }
    if(is.null(waa_info$waa_pointer_M)){
      input$log$waa <- c(input$log$waa, "waa_info$waa is provided without waa_info$waa_pointer_M so the first waa matrix will be used. \n")
      data$waa_pointer_M = rep(1,data$n_stocks)
    }
  }
  
  if(!is.null(waa_info$waa_pointer_fleets)){
    if(any(!(waa_info$waa_pointer_fleets %in% 1:dim(data$waa)[1]))){
      stop("some waa_info$waa_pointer_fleets are outside the number of waa matrices.\n")
    }
    if(length(waa_info$waa_pointer_fleets) != data$n_fleets){
      stop("length of waa_info$waa_pointer_fleets is not equal to the number of fleets.\n")
    }
    data$waa_pointer_fleets = waa_info$waa_pointer_fleets
  }
  if(!is.null(waa_info$waa_pointer_indices)){
    if(any(!(waa_info$waa_pointer_indices %in% 1:dim(data$waa)[1]))){
      stop("some waa_info$waa_pointer_indices are outside the number of waa matrices.\n")
    }
    if(length(waa_info$waa_pointer_indices) != data$n_indices){
      stop("length of waa_info$waa_pointer_indices is not equal to the number of indices.\n")
    }
    data$waa_pointer_indices = waa_info$waa_pointer_indices
  }
  if(!is.null(waa_info$waa_pointer_ssb)){
    if(any(!(waa_info$waa_pointer_ssb %in% 1:dim(data$waa)[1]))){
      stop("some waa_info$waa_pointer_ssb are outside the number of waa matrices.\n")
    }
    if(length(waa_info$waa_pointer_ssb) != data$n_stocks){
      stop("length of waa_info$waa_pointer_ssb is not equal to the number of stocks.\n")
    }
    data$waa_pointer_ssb = waa_info$waa_pointer_ssb
  }
  if(!is.null(waa_info$waa_pointer_totcatch)){
    if(any(!(waa_info$waa_pointer_totcatch %in% 1:dim(data$waa)[1]))){
      stop("some waa_info$waa_pointer_totcatch are outside the number of waa matrices.\n")
    }
    if(length(waa_info$waa_pointer_totcatch) != data$n_regions){
      stop("length of waa_info$waa_pointer_totcatch is not equal to the number of regions.\n")
    }
    data$waa_pointer_totcatch = waa_info$waa_pointer_totcatch
  }
  data$waa_pointer_M <- data$waa_pointer_ssb
  if(!is.null(waa_info$waa_pointer_M)){
    if(any(!(waa_info$waa_pointer_M %in% 1:dim(data$waa)[1]))){
      stop("some waa_info$waa_pointer_M are outside the number of waa matrices.\n")
    }
    if(length(waa_info$waa_pointer_M) != data$n_stocks){
      stop("length of waa_info$waa_pointer_totcatch is not equal to the number of regions.\n")
    }
    data$waa_pointer_M[] = waa_info$waa_pointer_M
  }
  if(length(input$log$waa))	input$log$waa <- c("WAA: \n", input$log$waa)
  
  input$data = data
  input$options$waa <- waa_info
  return(input)
}
sim_fn <- function(om, self.fit = FALSE){
  input <- om$input
  input$data = om$simulate(complete=TRUE)
  if(self.fit) {
    fit <- fit_wham(input, do.osa = F, do.retro = F, do.brp = TRUE, MakeADFun.silent = FALSE)
    return(fit)
  } else return(input) 
}

Generate_Maturity <- function(life_history = NULL, na) {
  if (is.null(life_history)){
    warning("Life history is not specified and default is used!")
    maturity <- t(matrix(1/(1 + exp(-1*(1:na - na/2))), na))
  } else if (life_history == "short"){
    m50 = 1.75; mslope = 1; 
    maturity <- t(matrix(1/(1 + exp(-(1:na-m50)/mslope)), na))
  } else if (life_history == "medium"){
    m50 = 3.5; mslope = 1; 
    maturity <- t(matrix(1/(1 + exp(-(1:na-m50)/mslope)), na))
  } else if (life_history == "long"){
    m50 = 7; mslope = 1; 
    maturity <- t(matrix(1/(1 + exp(-(1:na-m50)/mslope)), na))
  } 
  return(maturity)
}

Generate_Len <- function(Linf,k,n_ages) {
  Len <- Linf*(1-exp(-k*1:n_ages))
  return(Len)
}

Generate_WAA <- function(life_history = NULL, na) {
  if (is.null(life_history)){
    warning("Life history is not specified and default is used!")
    Len <- 100*(1-exp(-0.3*(1:na - 0)))
    W <- 3e-6*Len^3
  } else if (life_history == "short"){
    k = 0.27; Linf = 90
    Len <- Generate_Len(Linf,k,na)
    LWexp = 3; LWscaler = 3e-6
    W <- LWscaler*Len^LWexp
  } else if (life_history == "medium"){
    k = 0.13; Linf = 90
    Len <- Generate_Len(Linf,k,na)
    LWexp = 3; LWscaler = 3e-6
    W <- LWscaler*Len^LWexp
  } else if (life_history == "long"){
    k = 0.07; Linf = 90
    Len <- Generate_Len(Linf,k,na)
    LWexp = 3; LWscaler = 3e-6
    W <- LWscaler*Len^LWexp
  } 
  return(W)
}
