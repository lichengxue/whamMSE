Get_F_catch <- function(om, year, catch, Finit = 0.1, maxF = 10){
  rep = om$rep
  naa = rep$NAA[1,1,year,]
  Maa = rep$MAA[1,1,year,]
  sel_tot = rbind(rep$FAA[,year,]/max(exp(rep$log_FAA_tot[year,]))) #matrix nfleets x n_ages
  waa = rbind(om$input$data$waa[om$input$data$waa_pointer_fleets, year,]) #matrix nfleets x n_ages
  nfleets <- length(om$input$data$waa_pointer_fleets)
  get_catch = function(log_F, naa, sel, waa, Maa){
    Faa = exp(log_F) * sel_tot
    Zaa = Maa + apply(Faa,2,sum)
    Catch = 0
    for(a  in 1:length(naa)) for(f in 1:nfleets) Catch = Catch + waa[f,a] * naa[a] * Faa[f,a] *(1 - exp(-Zaa[a]))/Zaa[a];
    return(Catch)
  }
  obj = function(log_F) (catch - get_catch(log_F, naa, sel_tot, waa, Maa))^2
  opt = try(nlminb(log(Finit), obj))
  if(!is.character(opt)) Fsolve = exp(opt$par)[1] else Fsolve = maxF
  if(Fsolve>10) Fsolve = maxF
  print(paste0("Fsolve: ", Fsolve))
  return(Fsolve)
}

Update_F <- function(om, year, catch){
  rep <- om$rep #generate the reported values given the parameters
  year_ind <- year #index corresponding to year
  Fsolve <- Get_F_catch(om, year, catch) #find the F for the catch advice
  
  #have to be careful if more than one fleet
  FAA <- rbind(rep$FAA[,year_ind,]) #n_fleets x n_ages
  age_ind <- om$env$data$which_F_age[year_ind] #which age is used by wham to define total "full F"
  old_max_F <- apply(FAA,2,sum)[age_ind] # n_ages
  # selAA[[om$env$data$selblock_pointer_fleets[year_ind]]][year_ind,] #this only works when there is a single fleet
  selAA <- FAA/old_max_F #sum(selAA[i,]) = 1
  new_FAA <- Fsolve * selAA #updated FAA
  F_fleet_y <- apply(new_FAA, 1, max) #full F for each fleet
  return(new_FAA)
}

check_om_feedback_F <- function(mod, Finit = 0.1, maxF = 10,
                                save_pdf = TRUE,
                                main.dir = getwd(),
                                sub.dir = "Report",
                                pdf_name = "F_summary_table.pdf") {
  library(gridExtra)
  library(grid)
  
  om <- mod$om
  catch.adv <- do.call(rbind, mod$catch_advice)
  n_feedback_years <- nrow(catch.adv)
  n_om_years <- nrow(om$rep$pred_catch)
  first_feedback_year <- n_om_years - n_feedback_years + 1
  
  results <- data.frame(
    Year = integer(),
    OM_Catch = numeric(),
    Catch_Advice = numeric(),
    Diff = numeric(),
    True_Fbar = numeric(),
    Est_Fbar = numeric()
  )
  
  for (i in 1:n_feedback_years) {
    om_year <- first_feedback_year + i - 1
    
    user.catch <- catch.adv[i, ]
    om.catch <- om$rep$pred_catch[om_year, ]
    
    F_est <- Get_F_catch(om, om_year, user.catch)
    FAA_est <- Update_F(om, om_year, user.catch)
    Fbar_est <- mean(FAA_est[, om$env$data$which_F_age[om_year]])
    Fbar_true <- om$rep$Fbar[om_year, 3]
    
    results <- rbind(results, data.frame(
      Year = om_year,
      OM_Catch = round(om.catch, 3),
      Catch_Advice = round(user.catch, 3),
      Diff = round(om.catch - user.catch, 4),
      True_Fbar = round(Fbar_true, 4),
      Est_Fbar = round(Fbar_est, 4)
    ))
  }
  
  cat("\n==== Summary Table ====\n")
  print(results, row.names = FALSE)
  
  if (save_pdf) {
    dir.create(file.path(main.dir, sub.dir), recursive = TRUE, showWarnings = FALSE)
    pdf_path <- file.path(main.dir, sub.dir, pdf_name)
    pdf(file = pdf_path, width = 8.5, height = 15)
    grid.newpage()
    grid.draw(tableGrob(results))
    dev.off()
    cat(paste0("Summary table saved to: ", pdf_path, "\n"))
  }
  
  return(results)
}



