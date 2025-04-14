#' Update the operating model and generate data
#' 
#' Function to use the operating model to generate data 
#' 
#' @param om Operating model 
#' @param seed Seed used to generate data 
#'     
#' @return an operating model with simulated data
#'   
#' @export
#' 
#' @examples
#' \dontrun{
#' basic_info <- generate_basic_info()
#' input <- prepare_wham_input(basic_info = basic_info)
#' mod <- fit_wham(input, do.fit = FALSE)
#' data <- generate_data(mod, seed = 123)
#' }
#' 
generate_data = function(om, seed = 123) {
  
  obs_names = c("agg_indices","agg_catch","catch_paa","index_paa", "Ecov_obs", "obsvec")
  
  set.seed(seed)
  om_sim = om$simulate(complete=TRUE) #resimulate the population and observations
  om$input$data[obs_names] = om_sim[obs_names] #update any simulated data
  om$input$par[om$input$random] = om_sim[om$input$random]
  om <- fit_wham(om$input, do.fit = FALSE, do.retro = FALSE, do.osa = FALSE, MakeADFun.silent = TRUE)
  
  return(om)
}
