n_regions = n_stocks = 2
sigma_vals   <-  array(0.75, dim = c(n_stocks, n_regions, n_ages)) # NAA sigma 
sigma_vals[,,1] = 0.75 # What is the recruitment sigma also 0.75?

NAA_re_em <- list(N1_model=rep(ini.opt,n_stocks),
                  sigma=rep(sigma,n_stocks),
                  cor=rep(re_cor,n_stocks),
                  recruit_model = 2,
                  recruit_pars = list(16000,21000), # I assume rec in windfarm area is 3000 as an example
                  sigma_vals = sigma_vals) # NAA_where must be specified in basic_info!

NAA_where_em <- basic_info$NAA_where[-3,-3,]

move_em = list(stock_move = c(TRUE,FALSE), separable = TRUE) #north moves, south doesn't

move_em$must_move = array(0,dim = c(2,11,2))

move_em$must_move[1,5,2] <- 1 
move_em$can_move = array(0, dim = c(2,11,2,2))
move_em$can_move[1,1:4,2,1] = 1 #only north stock can move and in seasons prior to spawning and after spawning
move_em$can_move[1,7:11,1,2] = 1 #only north stock can move and in seasons prior to spawning and after spawning
move_em$can_move[1,5,2,] = 1

mean_vals <- array(0, dim = c(2,length(fracyr_seasons),2,1))
mean_vals[1,1:11,1,1] <- 0.02214863 #see here("2023.RT.Runs","transform_SS_move_rates.R") for how these numbers are derived.
mean_vals[1,1:11,2,1] <- 0.3130358
move_em$mean_vals = mean_vals #n_stocks x n_seasons x n_regions x n_regions - 1

move_em$mean_model = matrix("stock_constant", 2,1)

move_em$mean_vals <- array(0.3, dim = c(2,11,2,1)) # I'm not sure what movement rate we are going to assume here, we may create some movement scenarios here for example different magnitude of movement rate between north to south OR, between north and windfarm area, etc. 
move_em$use_prior <- array(0, dim = c(2,length(fracyr_seasons),2,1))
move_em$use_prior[1,1,1,1] <- 1
move_em$use_prior[1,1,2,1] <- 1
move_em$prior_sigma <- array(0, dim = c(2,length(fracyr_seasons),2,1))
move_em$prior_sigma[1,1,1,1] <- 0.2
move_em$prior_sigma[1,1,2,1] <- 0.2


em_info = NULL
random = NULL
M_em = NULL
sel_em = NULL
NAA_re_em = NULL
move_em = NULL
catchability_em = NULL
ecov_em = NULL
age_comp_em = "multinomial"
em.opt = list(separate.em = FALSE, separate.em.type = 1,
              do.move = FALSE, est.move = FALSE)
aggregate_catch_info = NULL
aggregate_index_info = NULL
aggregate_weights_info = NULL
reduce_region_info = NULL
filter_indices = NULL
update_catch_info = NULL
update_index_info = NULL
user_SPR_weights_info = NULL
assess_years = NULL
assess_interval = NULL
base_years = NULL
year.use = 20
add.years = FALSE
by_fleet = TRUE
FXSPR_init = NULL
hcr = list(hcr.type = 1, hcr.opts = NULL)
catch_alloc = list(weight_type = 1, method = "equal", user_weights = NULL, weight_years = 1)
implementation_error = NULL
do.retro = FALSE
do.osa = FALSE
do.brps = FALSE
seed = 123
save.sdrep = FALSE
save.last.em = FALSE

om = om_with_data
NAA_re_em = NAA_re
sel_em = sel
M_em = M
move_em = move
em_info = info
em.opt = list(separate.em = FALSE, # FALSE: Spatially disaggregate
              separate.em.type = 3, 
              do.move = TRUE, 
              est.move = TRUE) # If turn on, movement will be estimated using prior
assess_years = assess.years
assess_interval = assess.interval
base_years = base.years
year.use = 30
add.years = TRUE
hcr = list(hcr.type = 1, hcr.opts = NULL)
reduce_region_info = list(remove_regions = c(1,1,0), 
                          reassign = 1,
                          NAA_where_em = NAA_where_em,
                          sel_em = sel, # selectivity config should not change!
                          M_em = M_em,
                          NAA_re_em = NAA_re_em,
                          move_em = move_em
)
# for example we want to ignore the survey from the windfarm here.
save.last.em = T
seed = 123

y = 30

em_years = em.years

