devtools::install_github("timjmiller/wham", dependencies=FALSE, ref="devel")
# remove.packages("whamMSE")
main.dir = "C:/Users/chengxue.li/"
install.packages(file.path(main.dir,"whamMSE"), dependencies = TRUE, repos = NULL, type = "source")

library(whamMSE)
library(wham)

year_start  <- 1993  # starting year in the burn-in period
year_end    <- 2022  # end year in the burn-in period
MSE_years   <- 3     # number of years in the feedback loop

info <- generate_basic_info(n_stocks   = 2, 
                            n_regions  = 2, 
                            n_indices  = 2, 
                            n_fleets   = 2, 
                            n_seasons  = 4, 
                            base.years = year_start:year_end, 
                            n_feedback_years = MSE_years, 
                            life_history  = "medium", 
                            n_ages        = 12, 
                            Fbar_ages     = 12, 
                            recruit_model = 2, 
                            F_info     = list(F.year1 = 0.2, Fhist = "F-H-L", Fmax = 2, Fmin = 1, change_time = 0.5),
                            catch_info = list(catch_cv = 0.2, catch_Neff = 100), 
                            index_info = list(index_cv = 0.2, index_Neff = 100, fracyr_indices = 0.5, q = 0.2), 
                            fracyr_spawn = 0.5, 
                            bias.correct.process     = FALSE, 
                            bias.correct.observation = FALSE, 
                            bias.correct.BRPs        = FALSE, 
                            mig_type = 0) 

basic_info = info$basic_info

basic_info <- generate_NAA_where(basic_info = basic_info, move.type = 2)
# Note: default is move = 0.3 (constant) for stock1 and 0.1 (constant) for the other stocks
# move <- generate_move(basic_info = basic_info, move.type = 2, move.rate = 0.3, move.re = "iid_y", move.sigma = 0.2)
move <- generate_move(basic_info = basic_info, move.type = 2, move.rate = 0.3, move.re = "iid_y")

### 4. Configure selecitvity and natural mortality

n_stocks  <- as.integer(basic_info['n_stocks'])
n_regions <- as.integer(basic_info['n_regions'])
n_fleets  <- as.integer(basic_info['n_fleets'])
n_indices <- as.integer(basic_info['n_indices'])
n_ages    <- as.integer(basic_info['n_ages'])

# Selectivity Configuration
fleet_pars <- c(5,1)
index_pars <- c(2,1)
sel <- list(model=rep("logistic",n_fleets+n_indices),
            initial_pars=c(rep(list(fleet_pars),n_fleets),rep(list(index_pars),n_indices)))

# M Configuration
M <- list(model="constant",initial_means=array(0.2, dim = c(n_stocks,n_regions,n_ages)))

sigma        <- "rec+1"
re_cor       <- "iid"
ini.opt      <- "equilibrium" # option   <- c("age-specific-fe", "equilibrium")
Rec_sig      <- 0.2 # (sigma for recruitment)
NAA_sig      <- 0.2 # (sigma for NAA)

sigma_vals = array(NAA_sig, dim = c(n_stocks, n_regions, n_ages)) # n_stocks x n_regions x n_ages"
sigma_vals[,,1] = Rec_sig

# Set initial NAA for each stock
log_N1    <- rep(10,n_stocks)
log_N1[1] <- log(exp(10)*2) # Create difference between stocks
N1_pars   <- generate_ini_N1(log_N1,basic_info,ini.opt)

# Set mean recruitment para. for each stock
mean_rec_par <- list()
for (i in 1:n_stocks) mean_rec_par[[i]] <- exp(log_N1[i])

NAA_re <- list(N1_model=rep(ini.opt,n_stocks),
               sigma=rep(sigma,n_stocks),
               cor=rep(re_cor,n_stocks),
               recruit_model = 2,  # rec random around the mean
               recruit_pars = mean_rec_par, 
               sigma_vals = sigma_vals,  
               N1_pars = N1_pars,
               NAA_where = basic_info$NAA_where)

catch_info = info$catch_info
index_info = info$index_info
F_info = info$F

input <- prepare_wham_input(basic_info = basic_info, selectivity = sel, M = M, NAA_re = NAA_re, move = move,
                            catch_info = catch_info, index_info = index_info, F = F_info)

# input$random <- NULL #so inner optimization won't change simulated RE
random = input$random
input$random = NULL
om = fit_wham(input, do.fit = F, do.brps = T, MakeADFun.silent = T)
exp(om$rep$log_FXSPR_static)
exp(om$rep$log_FAA_XSPR_static)

### 8. Generate datasets
om_with_data <- update_om_fn(om, seed = 124, random = random)

assess.interval <- 3 # Assessment interval
base.years      <- year_start:year_end # Burn-in period
first.year      <- head(base.years,1)
terminal.year   <- tail(base.years,1)
assess.years    <- seq(terminal.year, tail(om$years,1)-assess.interval,by = assess.interval)

mods <- list()

n_stocks = n_regions = n_fleets = n_indices = 2

sel_em <- list(model=rep("logistic",n_fleets+n_indices),
               initial_pars=c(rep(list(fleet_pars),n_fleets),rep(list(index_pars),n_indices)))
NAA_re_em <- list(N1_model=rep("equilibrium",n_stocks),
                  sigma=rep("rec+1",n_stocks),
                  cor=rep("iid",n_stocks))
M_em <- list(model="constant",initial_means=array(0.2, dim = c(n_stocks,n_regions,n_ages)))

mods[[1]] = loop_through_fn(om = om_with_data, 
                            random = random,
                            M_om = M,
                            sel_om = sel, 
                            NAA_re_om = NAA_re, 
                            move_om = move,
                            M_em = M_em, 
                            sel_em = sel_em, 
                            NAA_re_em = NAA_re_em, 
                            move_em = NULL,
                            em.opt = list(separate.em = FALSE, separate.em.type = 3, do.move = FALSE, est.move = FALSE),
                            assess_years = assess.years, 
                            assess_interval = assess.interval, 
                            base_years = base.years,
                            year.use = 20, # number of years of data you want to use in the assessment model
                            seed = 124,
                            save.sdrep = FALSE)

NAA_re_em <- list(N1_model=rep("equilibrium",n_stocks),
                  sigma=rep("rec",n_stocks),
                  cor=rep("iid",n_stocks))
mods[[2]] = loop_through_fn(om = om_with_data, 
                            random = random,
                            M_om = M,
                            sel_om = sel, 
                            NAA_re_om = NAA_re, 
                            move_om = move,
                            M_em = M_em, 
                            sel_em = sel_em, 
                            NAA_re_em = NAA_re_em, 
                            move_em = NULL,
                            em.opt = list(separate.em = FALSE, separate.em.type = 3, do.move = FALSE, est.move = FALSE),
                            assess_years = assess.years, 
                            assess_interval = assess.interval, 
                            base_years = base.years,
                            year.use = 20, # number of years of data you want to use in the assessment model
                            seed = 124,
                            save.sdrep = FALSE)

mods[[2]]$om$rep$SSB - mods[[1]]$om$rep$SSB

mod_list = list()
mod_list[[1]] = mods
mod_list[[2]] = mods
plot_mse_output(mod_list)


stock_om_fb = om_with_data

res <- cbind(om_with_data$rep$SSB, mods[[1]]$om$rep$SSB,mods[[2]]$om$rep$SSB)
plot(om_with_data$years, om_with_data$rep$SSB[,2], type = "l", ylab = "SSB", xlab = "Year")
plo(mods[[1]]$om$rep$SSB[,2],col = "red",type = "l")
lines(mods[[2]]$om$rep$SSB[,2],col = "blue",type = "l")

sum(mods[[2]]$om$rep$NAA_index - mods[[1]]$om$rep$NAA_index)

plot_mse_output(mods, out.type = "pdf")

plot(res[,1], type = "l", ylab = "SSB", xlab = "Year")
lines(res[,3],col = "blue",type = "l")
lines(res[,5],col = "red",type = "l")

plot(res[,2], type = "l", ylab = "SSB", xlab = "Year")
lines(res[,4],col = "blue",type = "l")
lines(res[,6],col = "red",type = "l")

res[]