#' Update Estimation Model Inputs with Aggregated Data
#'
#' Updates the estimation model (`em_info`) by incorporating aggregated data
#' from `aggregate_catch_info` and `aggregate_index_info`. This function adjusts
#' `par_inputs` and weight-at-age (`user_waa`) to match the new fleet, region,
#' index, and stock structure while handling fleet/index exclusions (`fleet_pointer == 0`).
#'
#' @param om List. An operating model containing observed data (`agg_catch`, `agg_indices`).
#' @param em_info List. The estimation model information containing:
#'   \itemize{
#'     \item `par_inputs` - Model parameters including fleet and index configurations.
#'     \item `user_waa` - Weight-at-age information for fleets, indices, and stocks.
#'   }
#' @param aggregate_catch_info List. Aggregated fleet catch data:
#'   \itemize{
#'     \item `fleet_pointer` - Vector mapping original fleets to aggregated fleets (`0` = exclude).
#'     \item `agg_catch` - Aggregated catch values per fleet.
#'     \item `catch_cv`, `catch_Neff`, `use_agg_catch`, `use_catch_paa`
#'   }
#' @param aggregate_index_info List. Aggregated survey index data:
#'   \itemize{
#'     \item `index_pointer` - Vector mapping original indices to aggregated indices (`0` = exclude).
#'     \item `agg_indices` - Aggregated index values.
#'     \item `index_cv`, `index_Neff`, `fracyr_indices`, `q`, `use_indices`, `use_index_paa`
#'   }
#' @param ind_em Vector. Indices specifying the years for the estimation model.
#'
#' @return List. Updated `em_info` with:
#'   \itemize{
#'     \item `par_inputs` - Adjusted estimation model parameters reflecting the new structure.
#'     \item `user_waa` - Updated weight-at-age (`n_fleets + n_regions + n_indices + n_stocks, n_ages`).
#'     \item `basic_info$waa` - Normalized weight-at-age data ensuring proper scaling.
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' updated_em_info <- update_em_with_basic_info(
#'   om = some_operating_model,
#'   em_info = some_em_info,
#'   aggregate_catch_info = list(fleet_pointer = c(1, 1, 2, 0), agg_catch = matrix(runif(40), 10, 4)),
#'   aggregate_index_info = list(index_pointer = c(1, 2, 0), agg_indices = matrix(runif(30), 10, 3)),
#'   ind_em = 1:10
#' )
#' }
#' 
update_em_with_basic_info <- function(om, em_info, aggregate_catch_info, aggregate_index_info, ind_em) {
  
  data <- om$input$data
  
  # ---- Update `par_inputs` from aggregate data ---- #
  em_info$par_inputs$catch_cv <- ifelse(is.null(aggregate_catch_info$catch_cv), em_info$par_inputs$catch_cv[1], aggregate_catch_info$catch_cv)
  em_info$par_inputs$catch_Neff <- ifelse(is.null(aggregate_catch_info$catch_Neff), em_info$par_inputs$catch_Neff[1], aggregate_catch_info$catch_Neff)
  em_info$par_inputs$use_agg_catch <- ifelse(is.null(aggregate_catch_info$use_agg_catch), em_info$par_inputs$use_agg_catch[1], aggregate_catch_info$use_agg_catch)
  em_info$par_inputs$use_catch_paa <- ifelse(is.null(aggregate_catch_info$use_catch_paa), em_info$par_inputs$use_catch_paa[1], aggregate_catch_info$use_catch_paa)
  
  em_info$par_inputs$index_cv <- ifelse(is.null(aggregate_index_info$index_cv), em_info$par_inputs$index_cv[1], aggregate_index_info$index_cv)
  em_info$par_inputs$index_Neff <- ifelse(is.null(aggregate_index_info$index_Neff), em_info$par_inputs$index_Neff[1], aggregate_index_info$index_Neff)
  em_info$par_inputs$fracyr_indices <- ifelse(is.null(aggregate_index_info$fracyr_indices), em_info$par_inputs$fracyr_indices[1], aggregate_index_info$fracyr_indices)
  em_info$par_inputs$q <- ifelse(is.null(aggregate_index_info$q), em_info$par_inputs$q[1], aggregate_index_info$q)
  em_info$par_inputs$use_indices <- ifelse(is.null(aggregate_index_info$use_indices), em_info$par_inputs$use_indices[1], aggregate_index_info$use_indices)
  em_info$par_inputs$use_index_paa <- ifelse(is.null(aggregate_index_info$use_index_paa), em_info$par_inputs$use_index_paa[1], aggregate_index_info$use_index_paa)
  em_info$par_inputs$units_indices <- ifelse(is.null(aggregate_index_info$units_indices), em_info$par_inputs$units_indices[1], aggregate_index_info$units_indices)
  em_info$par_inputs$units_index_paa <- ifelse(is.null(aggregate_index_info$units_index_paa), em_info$par_inputs$units_index_paa[1], aggregate_index_info$units_index_paa)
  
  # ---- Determine fleet, index, region, and stock numbers ---- #
  n_fleets <- ifelse(is.null(aggregate_catch_info$n_fleets), 1, aggregate_catch_info$n_fleets)
  n_indices <- ifelse(is.null(aggregate_index_info$n_indices), 1, aggregate_index_info$n_indices)
  em_info$par_inputs$n_fleets = em_info$basic_info$n_fleets = n_fleets
  em_info$par_inputs$n_indices = em_info$basic_info$n_indices = n_indices
  n_ages <- em_info$par_inputs$n_ages
  # Place holder for later change
  n_regions <- em_info$par_inputs$n_regions
  n_stocks <- em_info$par_inputs$n_stocks
  n_regions <- 1 # Now first assume 1 for panmictic assessment only
  n_stocks <- 1 # Now first assume 1 for panmictic assessment only
  
  # ---- Exclude fleets and indices marked with `0` in pointers ---- #
  if (is.null(aggregate_catch_info$fleet_pointer)) warnings("aggregate_catch_info$fleet_pointer is not specified!")
  if (is.null(aggregate_catch_info$fleet_pointer)) aggregate_catch_info$fleet_pointer = rep(1,aggregate_catch_info)
  if (length(aggregate_catch_info$fleet_pointer) != om$input$data$n_fleets) stop("Length of aggregate_catch_info$fleet_pointer is wrong!")
  fleet_pointer <- aggregate_catch_info$fleet_pointer
  valid_fleets <- unique(fleet_pointer[fleet_pointer > 0])
  
  if (is.null(aggregate_index_info$index_pointer)) stop("aggregate_index_info$index_pointer is not specified!")
  if (length(aggregate_index_info$index_pointer) != om$input$data$n_indices) stop("Length of aggregate_index_info$index_pointer is wrong!")
  index_pointer <- aggregate_index_info$index_pointer
  valid_indices <- unique(index_pointer[index_pointer > 0])
  
  # ---- Compute WAA Weights ---- #
  fleet_weights = list()
  for (f in unique(valid_fleets)) {
    fleet_weights[[f]] <- if (aggregate_catch_info$use_catch_weighted_waa) {
      agg_catch_filtered <- data$agg_catch[ind_em,which(fleet_pointer == f)]
      prop.table(ifelse(is.matrix(agg_catch_filtered),colMeans(agg_catch_filtered), mean(agg_catch_filtered)))
    } else {
      rep(1 / length(valid_fleets), length(valid_fleets))  # Equal weights
    }
  }
  
  index_weights = list()
  for (f in unique(valid_indices)) {
    index_weights[[f]] <- if (aggregate_index_info$use_catch_weighted_waa) {
      agg_indices_filtered <- data$agg_indices[ind_em,which(index_pointer == f)]
      prop.table(ifelse(is.matrix(agg_catch_filtered),colMeans(agg_catch_filtered), mean(agg_catch_filtered)))
    } else {
      rep(1 / length(valid_indices), length(valid_indices))  # Equal weights
    }
  }
  
  # ---- Aggregate WAA ---- #
  aggregated_fleet_waa = list()
  for (f in unique(valid_fleets)) {
    if(is.null(em_info$par_inputs$user_waa)) {
      aggregated_fleet_waa[[f]] <- em_info$basic_info$waa[which(fleet_pointer == f),1,]
      if (is.matrix(aggregated_fleet_waa[[f]])) aggregated_fleet_waa[[f]] <- colSums(em_info$basic_info$waa[which(fleet_pointer == f),1, ]*fleet_weights[[f]])
    } else {
      aggregated_fleet_waa[[f]] <- em_info$par_inputs$user_waa[which(fleet_pointer == f), ]
      if (is.matrix(aggregated_fleet_waa[[f]])) aggregated_fleet_waa[[f]] <- colSums(em_info$par_inputs$user_waa[which(fleet_pointer == f), ]*fleet_weights[[f]])
    }
  }
  
  aggregated_region_waa = list()
  for (f in 1) { # Assuming n_regions = 1
    ind = length(fleet_pointer) + 1:em_info$basic_info$n_regions
    if(is.null(em_info$par_inputs$user_waa)) {
      aggregated_region_waa[[f]] <- em_info$basic_info$waa[ind,1, ]
      if (is.matrix(aggregated_region_waa[[f]])) aggregated_region_waa[[f]] <- colMeans(as.matrix(aggregated_region_waa[[f]]))
    } else {
      aggregated_region_waa[[f]] <- em_info$par_inputs$user_waa[, ]
      if (is.matrix(aggregated_region_waa[[f]])) aggregated_region_waa[[f]] <- colMeans(as.matrix(aggregated_region_waa[[f]]))
    }
  }
  
  aggregated_index_waa = list()
  for (f in unique(valid_indices)) {
    ind = length(fleet_pointer) + 1:em_info$basic_info$n_regions + which(index_pointer == f)
    if(is.null(em_info$par_inputs$user_waa)) {
      aggregated_index_waa[[f]] <- em_info$basic_info$waa[ind,1, ]
      if (is.matrix(aggregated_index_waa[[f]])) aggregated_index_waa[[f]] <- colMeans(as.matrix(aggregated_index_waa[[f]]))
    } else {
      aggregated_index_waa[[f]] <- em_info$par_inputs$user_waa[ind, ]
      if (is.matrix(aggregated_index_waa[[f]])) aggregated_index_waa[[f]] <- colSums(em_info$par_inputs$user_waa[which(index_pointer == f), ]*index_weights[[f]])
    }
  }
  
  aggregated_stock_waa = list()
  for (f in 1) {
    ind = length(fleet_pointer) + em_info$basic_info$n_regions + length(index_pointer) + 1:em_info$basic_info$n_stocks
    if(is.null(em_info$par_inputs$user_waa)) {
      aggregated_stock_waa[[f]] <- em_info$basic_info$waa[ind,1, ]
      if (is.matrix(aggregated_stock_waa[[f]])) aggregated_stock_waa[[f]] <- colMeans(as.matrix(aggregated_stock_waa[[f]]))
    } else {
      aggregated_stock_waa[[f]] <- em_info$par_inputs$user_waa[ind, ]
      if (is.matrix(aggregated_stock_waa[[f]])) aggregated_stock_waa[[f]] <- colMeans(as.matrix(aggregated_stock_waa[[f]]))
    }
  }
  
  # Store the aggregated WAA back into `em_info`
  em_info$par_inputs$user_waa <- do.call(rbind,c(aggregated_fleet_waa,aggregated_region_waa,aggregated_index_waa,aggregated_stock_waa))
  em_info$basic_info$waa = em_info$par_inputs$user_waa
  
  em_info$par_inputs$fleet_regions = rep(1, n_fleets)
  em_info$par_inputs$index_regions = rep(1, n_indices)
  
  # ---- Return Updated Info ---- #
  return(em_info)
}
