#' Update fishing mortality in the operating model given the catch
#' 
#' Function to update F in the operating model given the catch advice 
#' 
#' @param om Operating model 
#' @param year Actual year(s) in the projection
#' @param catch Catch advice in the projection years 
#'     
#' @return an operating model with updated F time series
#'   
#' @export
#'
#' @seealso \code{\link{get_F_from_catch}}
#' 
update_om_F2 = function(om, year, catch){
  rep = om$rep #generate the reported values given the parameters
  #year = om$years[year]
  year_ind = which(om$years == year) #index corresponding to year
  Fsolve = get_F_from_catch(om, year = year_ind, catch) #find the F for the catch advice
  n_regions = length(om$input$region_names)
  sel = rep$FAA[,year_ind,]
  #sel_all = sel/max(sel)
  sel_all = t(apply(sel,1,function(row) row/max(row)))
  FAA_all = Fsolve * rbind(sel_all)
  F_region = apply(FAA_all, 1, max) #full F for each region
  om$input$par$F_pars[year_ind,] = log(F_region)
  om <- fit_wham(om$input, do.fit = FALSE, MakeADFun.silent = TRUE)
  return(om)
}

get_F_from_catch2 <- function(om, year, catch, Finit = 0.1, maxF = 10){ # here year has to be year starting from 1
  
  get_catch = function(log_F, naa, sel, waa, Maa){
    Faa = exp(log_F) * sel_tot
    Zaa = Maa + Faa
    Catch = 0
    for(a  in 1:length(naa)) Catch = Catch + waa[a] * naa[a] * Faa[a] *(1 - exp(-Zaa[a]))/Zaa[a];
    return(Catch)
  }
  
  rep = om$report()
  n_regions = length(om$input$region_names)
  Fsolve = NULL
  for (r in 1:n_regions) {
    naa = colSums(rep$NAA[,r,year,]) # n_stocks x n_regions x n_years x n_ages
    Maa = colMeans(rep$MAA[,r,year,]) # assume MAA is the same across stocks
    sel_tot <- rep$FAA[r,year,]/max(rep$FAA[r,year,])
    waa = om$input$data$waa[om$input$data$waa_pointer_totcatch, year,][r,] # assume MAA is the same across stocks
    
    obj = function(log_F) (catch[r] - get_catch(log_F, naa, sel_tot, waa, Maa))^2
    opt = try(nlminb(log(Finit), obj))
    if(!is.character(opt)) Fsolve.tmp = exp(opt$par)[1] else Fsolve.tmp = maxF
    if(Fsolve.tmp>10) Fsolve.tmp = maxF
    Fsolve[r] <- Fsolve.tmp
  }
  cat(paste0("\nFsolve: ", Fsolve,"\n"))
  return(Fsolve)
}
