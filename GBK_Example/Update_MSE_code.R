## ============================================================
## 0) Corrected Code
## ============================================================

# I found a bug in the UMassD-MSE branch so please reinstall
# I reorganized the code, and added explanation
# I made some changes in the code: 
#   1) change N1 model for the EM: NAA_re_em$N1_model = "equilibrium" # IMPORTANT! Turn to equilibrium
#   2) change initial numbers at age: om_input_P$par$log_N1[i, i, ] <- log(om_ecov$rep$N1)[i, i, ] # IMPORTANT!!!

devtools::install_github(
  "lichengxue/whamMSE",
  ref = "UMassD-MSE",           # <-- branch name
  dependencies = FALSE
  )
setwd("~/Desktop/UMassD-Cheng")                                 # <- ADJUST if needed

## Load packages (assumes installed)
library(wham)
library(whamMSE)
library(dplyr)
library(ggplot2)

## Input files (assumes they exist in the wd)
gb_dat   <- read_asap3_dat("GBK.DAT")                           # historical ASAP3 data (Georges Bank cod example)
om_ecov  <- readRDS("om_ecov.rds")                              # fitted WHAM OM with Ecov (your base OM)
env.dat_me <- read.csv("CI_indices.csv") %>% filter(Year > 1971)# environmental series with Year & bt_temp
input_hist <- prepare_wham_input(gb_dat)                        # WHAM input for historical-only data (for copying sigmas etc.)
Ecov_re <- readRDS("Ecov_re.RDS")                               # Ecov random effects estimated from a fitted OM (length 51 here)

## Sanity checks for columns we use
stopifnot(all(c("Year","bt_temp") %in% names(env.dat_me)))

## Base WHAM input (for dimensions if needed)
input <- prepare_wham_input(gb_dat)

## ============================================================
## 1) Global settings for OM/EM dimensions and MSE horizon
## ============================================================

n_stocks      <- 1L
n_regions     <- 1L
n_ages        <- 6L
n_indices     <- 3L
n_fleets      <- 1L
n_seasons     <- 1L

year_start    <- 1973L                                          # first assessment year in data
year_end      <- 2022L                                          # last assessment year in data (50 yrs total)
MSE_years     <- 6L                                             # projection/feedback years for MSE
n_proj_years  <- (year_end - year_start + 1L) + MSE_years       # 50 + 6 = 56, but Ecov uses lag year so 57 rows later

base.years    <- year_start:year_end                             # c(1973,...,2022)
assess.interval <- 3L                                            # assess every 3 yrs
assess.years  <- seq(tail(base.years,1), year_end + MSE_years - assess.interval, by = assess.interval)

## ============================================================
## 2) Build user-defined biology (maturity, WAA) over base+MSE
##    - Copy first 50 years from fitted OM, hold constant after
## ============================================================

## 2.1 Maturity-at-age array: [stock, year, age]
user_maturity <- array(NA_real_, dim = c(n_stocks, n_proj_years, n_ages))

## copy years 1:50 from the fitted OM (which has 50 base years)
user_maturity[, 1:50, ] <- om_ecov$input$data$mature[, 1:50, ]

## hold the last historical year (50) constant into projection years (51:56)
for (iy in 51:n_proj_years) {
  user_maturity[, iy, ] <- om_ecov$input$data$mature[, 50, , drop = FALSE]
}

## 2.2 Weight-at-age list structure expected by WHAM
user_waa <- list()
user_waa$waa <- array(NA_real_, dim = c(5, n_proj_years, n_ages))  # "5" is the WHAM WAA dimension slot (fleets/indices/SSB/M/etc.)

## copy historical 1:50
user_waa$waa[, 1:50, ] <- om_ecov$input$data$waa[, 1:50, ]

## extend to 51:56 by last historical year
for (iy in 51:n_proj_years) {
  user_waa$waa[, iy, ] <- om_ecov$input$data$waa[, 50, , drop = FALSE]
}

## copy WAA pointer metadata directly from base OM
user_waa$waa_pointer_fleets  <- om_ecov$input$data$waa_pointer_fleets
user_waa$waa_pointer_indices <- om_ecov$input$data$waa_pointer_indices
user_waa$waa_pointer_totcatch<- om_ecov$input$data$waa_pointer_ssb
user_waa$waa_pointer_ssb     <- om_ecov$input$data$waa_pointer_ssb
user_waa$waa_pointer_M       <- om_ecov$input$data$waa_pointer_M

## Spawning fraction by year from ASAP (WHAM wants this)
fracyr_spawn <- gb_dat[[1]]$dat$fracyr_spawn

## ============================================================
## 3) Observation model specs for catch and indices
## ============================================================

catch_info <- list(
  catch_cv       = 0.05,   # CV for aggregated catch
  catch_Neff     = 50,     # Neff for catch-at-age
  use_agg_catch  = 1,      # use aggregated catch
  use_catch_paa  = 1       # use catch-at-age
)

index_info <- list(
  index_cv        = rep(0.5, n_indices),      # CV for aggregated indices
  index_Neff      = rep(25, n_indices),       # Neff for index-at-age
  fracyr_indices  = c(0.25, 0.75, 0.25),      # timing of surveys
  q               = c(2.110e-4, 2.247e-4, 2.683e-4),
  use_indices     = rep(1, n_indices),        # use aggregated indices
  use_index_paa   = rep(1, n_indices),        # use index-at-age
  units_indices   = rep(2, n_indices),        # biomass/abundance units flags (as in your prior code)
  units_index_paa = rep(2, n_indices)
)

## ============================================================
## 4) Basic model info scaffold (whamMSE helper creates a big list)
##    Note: we’re calling a single helper to avoid hand-crafting dozens of slots
## ============================================================

info <- whamMSE::generate_basic_info(
  n_stocks = n_stocks, n_regions = n_regions, n_indices = n_indices,
  n_fleets = n_fleets, n_seasons = n_seasons,
  base.years = year_start:year_end,
  n_feedback_years = MSE_years,
  n_ages = n_ages,
  catch_info = catch_info, index_info = index_info,
  user_waa = user_waa, user_maturity = user_maturity,
  fracyr_spawn = fracyr_spawn
)

basic_info      <- info$basic_info
catch_info_use  <- info$catch_info
index_info_use  <- info$index_info

## ============================================================
## 5) Fishing mortality initializer: copy historical Fbar(ages 2–5) from base OM
## ============================================================

F_info <- info$F
F_info$F[1:50, ] <- om_ecov$rep$Fbar[1:50, 1, drop = FALSE]    # copy the historical Fbar trajectory
## leave projection years for MSE to overwrite

## ============================================================
## 6) Selectivity, Natural mortality, and NAA process
## ============================================================

## 6.1 Selectivity set: 1 fleet + 3 indices (as in your earlier config)
sel3 <- list(
  model = c("age-specific", "logistic", "logistic", "logistic"),
  re    = c("ar1_y", "none", "none", "none"),
  initial_pars = list(
    c(0.017, 0.249, 0.751, 1, 1, 1),   # fleet age-specific (fix older ages below)
    c(2.305, 0.327),                    # index 1 logistic
    c(1.611, 0.483),                    # index 2 logistic
    c(2.135, 0.217)                     # index 3 logistic
  ),
  fix_pars = list(c(4:6), NULL, NULL, NULL) # fix ages 4–6 at 1 for fleet selectivity
)

## 6.2 M (constant by age, from your provided vector)
M <- list(
  model = "constant",
  initial_means = array(c(0.57, 0.33, 0.26, 0.23, 0.22, 0.22), dim = c(n_stocks, n_regions, n_ages))
)

## 6.3 NAA process structure (rec+1 with age-specific sigma)
sigma_vals <- array(0.5, dim = c(n_stocks, n_regions, n_ages))
sigma_vals[,,1]    <- 0.655
sigma_vals[,,2:6]  <- 0.720

NAA_re <- list(
  N1_model     = rep("age-specific-fe", n_stocks),
  sigma        = rep("rec+1", n_stocks),
  cor          = rep("iid", n_stocks),
  recruit_model = 3,                                 # BH-like
  recruit_pars  = list(c(exp(8.082), 1.018e-4)),     # (logR0, a BH slope-ish param you used)
  sigma_vals    = sigma_vals
)

## ============================================================
## 7) Environmental covariate (Ecov) series construction
##    - years 1972–2028 (= 57 rows) with a 1-year lag
##    - we will build a PROJECTED trend scenario for 2023–2028
## ============================================================

## 7.1 Base Ecov list (template)
ecov_om <- list(
  label = "bt_temp",
  mean  = as.matrix(c(env.dat_me$bt_temp, rep(NA_real_, 7))),   # fill future later
  logsigma = "est_1",
  year  = 1972:2028,                                            # includes 1-lag year (1972) + 1973–2028
  use_obs = matrix(1, ncol = 1, nrow = 57),                     # observed everywhere (we may override in MSE)
  process_model = "ar1",
  recruitment_how = matrix("controlling-lag-1-linear")
)

## 7.2 PROJECTED trend for 2023–2028 by linear regression over historical (1972+)
lm_fit <- lm(bt_temp ~ Year, data = env.dat_me)
lm_sum <- summary(lm_fit)
proj_df <- data.frame(year_proj = 2022:2028)                    # we only need 2023–2028, but harmless to compute 2022 too
error_sd <- lm_sum[["sigma"]]                                   # residual SD
set.seed(42)
proj_df$SD <- rnorm(nrow(proj_df), 0, error_sd)
proj_df$Temp_proj <- lm_sum$coefficients[1,1] +
  lm_sum$coefficients[2,1]*proj_df$year_proj +
  proj_df$SD

## 7.3 Fill PROJECTED mean into 2023–2028 (rows 52:57 if year starts at 1972)
ecov_om_P <- ecov_om
ecov_om_P$mean[51:57, 1] <- proj_df$Temp_proj[proj_df$year_proj >= 2022]  # 51->2022, 52->2023, ..., 57->2028

## Optional CONSTANT future scenario (commented out)
# ecov_om_C <- ecov_om
# ecov_om_C$mean[51:57, 1] <- mean(env.dat_me$bt_temp, na.rm = TRUE)

## ============================================================
## 8) Build the WHAM OM input object using PROJECTED Ecov
##    - copy parameter guesses from fitted OM so it’s consistent
## ============================================================

om_input_P <- wham::prepare_wham_input(
  basic_info  = basic_info,
  selectivity = sel3,
  M           = M,
  NAA_re      = NAA_re,
  ecov        = ecov_om_P,
  catch_info  = catch_info_use,
  index_info  = index_info_use,
  F           = F_info,
  age_comp    = "logistic-normal-pool0"
)

## 8.1 Ensure WAA is applied via whamMSE helper (same as earlier but inline)
waa_info <- info$par_inputs$user_waa
om_input_P <- whamMSE::update_waa(om_input_P, waa_info = waa_info)

## 8.2 Copy parameter vectors from fitted OM for consistency
om_input_P$par$Ecov_process_pars <- om_ecov$parList$Ecov_process_pars
om_input_P$par$Ecov_beta_R       <- om_ecov$parList$Ecov_beta_R
om_input_P$par$catch_paa_pars    <- om_ecov$parList$catch_paa_pars
om_input_P$par$index_paa_pars    <- om_ecov$parList$index_paa_pars
om_input_P$par$sel_repars        <- om_ecov$parList$sel_repars

## 8.3 Initialize N1 for region i from base OM log_NAA (historical first-year)
for (i in 1:n_regions) {
  om_input_P$par$log_N1[i, i, ] <- log(om_ecov$rep$N1)[i, i, ] # IMPORTANT!!!
}
# om_ecov$rep$N1
## 8.4 Bring over historical index usage and sigma (rows 1:50) then extend flat
om_input_P$data$agg_index_sigma[1:50, ] <- input_hist$data$agg_index_sigma
om_input_P$data$use_indices[1:50, ]     <- input_hist$data$use_indices
om_input_P$data$use_index_paa[1:50, ]   <- input_hist$data$use_index_paa
for (iy in 51:n_proj_years) {
  om_input_P$data$agg_index_sigma[iy, ] <- input_hist$data$agg_index_sigma[50, ]
}

## 8.5 Compute remove-year matrices (zeros positions) for agg and paa (index)
##     We do it “inline” instead of a helper function.
loc_zeros_agg <- lapply(seq_len(ncol(input_hist$data$use_indices)), function(j)
  which(input_hist$data$use_indices[, j] == 0))
loc_zeros_paa <- lapply(seq_len(ncol(input_hist$data$use_index_paa)), function(j)
  which(input_hist$data$use_index_paa[, j] == 0))

max_agg <- max(sapply(loc_zeros_agg, length))
max_paa <- max(sapply(loc_zeros_paa, length))
remove_agg_years1 <- if (max_agg > 0) {
  out <- matrix(NA_integer_, nrow = max_agg, ncol = length(loc_zeros_agg))
  for (j in seq_along(loc_zeros_agg)) if (length(loc_zeros_agg[[j]])) out[1:length(loc_zeros_agg[[j]]), j] <- loc_zeros_agg[[j]]
  out
} else NULL
remove_paa_years1 <- if (max_paa > 0) {
  out <- matrix(NA_integer_, nrow = max_paa, ncol = length(loc_zeros_paa))
  for (j in seq_along(loc_zeros_paa)) if (length(loc_zeros_paa[[j]])) out[1:length(loc_zeros_paa[[j]]), j] <- loc_zeros_paa[[j]]
  out
} else NULL

## 8.6 Apply those removals to om_input_P
om_input_P <- whamMSE::update_input_index_info(
  om_input_P,
  agg_index_sigma = om_input_P$data$agg_index_sigma,
  index_Neff      = om_input_P$data$index_Neff,
  remove_agg = TRUE, remove_agg_pointer = 1:3, remove_agg_years = remove_agg_years1,
  remove_paa = TRUE, remove_paa_pointer = 1:3, remove_paa_years = remove_paa_years1
)

## ============================================================
## 9) Fix Ecov random effects path for simulation:
##    - set historical Ecov_re (1:51) from fitted OM
##    - set future Ecov_re (52:57) = observed minus long-term mean (for projected)
##    - IMPORTANT: remove Ecov_re from random effects (TMB won’t re-simulate)
## ============================================================

## 9.1 Long-term mean (your explicit value)
long_term_bt_mean <- 7.341484

## 9.2 Populate Ecov_re into parameter list
## (length alignment: 1:51 is 1972–2022 when lag=1; 52:57 is 2023–2028)
om_input_P$par$Ecov_re[1:51, 1] <- Ecov_re[1:51, 1]             # from fitted OM/your saved RDS

## For PROJECTED scenario: future “random effect” = obs - process mean
om_input_P$par$Ecov_re[52:57, 1] <- om_input_P$data$Ecov_obs[52:57, 1] - long_term_bt_mean

## Tell WHAM to not simulate Ecov_re (use what we set)
om_input_P$data$do_simulate_Ecov_re <- 0L

## Remove "Ecov_re" from random effects so TMB won’t re-estimate/simulate it
if ("Ecov_re" %in% om_input_P$random) {
  om_input_P$random <- om_input_P$random[om_input_P$random != "Ecov_re"]
}

## ============================================================
## 10) Build UNFITTED OM object (TMB compiled), but do not fit
## ============================================================

unfitted_om_P <- fit_wham(om_input_P, do.fit = FALSE, do.brps = FALSE, MakeADFun.silent = TRUE)

## Quick visual of Ecov obs vs latent x for this unfitted object
plot(unfitted_om_P$input$data$Ecov_obs, type = "l", col = "red",
     main = "PROJECTED scenario: Ecov obs vs latent",
     ylab = "Ecov", xlab = "row (1972..2028)")
lines(unfitted_om_P$rep$Ecov_x, type = "l", col = "blue")       # note: presence depends on model building; keep as per your script

## ============================================================
## 11) Simulate a concrete OM dataset from the unfitted OM
##     - We’ve removed Ecov_re from random list; keep others as-is
## ============================================================

random_P <- unfitted_om_P$input$random
random_P <- random_P[!random_P %in% "Ecov_re"]                   # be explicit

set.seed(123 + 1L)
om_with_data_P <- update_om_fn(unfitted_om_P, seed = 123 + 1L, random = random_P)

## Visual check again with simulated obs
plot(om_with_data_P$input$data$Ecov_obs, type = "l", col = "red",
     main = "Simulated OM: Ecov obs vs latent",
     ylab = "Ecov", xlab = "row (1972..2028)")
lines(om_with_data_P$rep$Ecov_x, type = "l", col = "blue")

## ============================================================
## 12) Build several Ecov “EM-side” projection scenarios (no function)
##     We create a small list with:
##       - AR1 (copy from OM input as-is)
##       - HistAvg (flat = mean of historical 1973–2022 for projection 2023–2028)
## ============================================================

## Years & sizes
years_all   <- 1972:2028
n_hist      <- length(base.years)          # 50
lag         <- 1L
n_proj      <- MSE_years                   # 6
n_total     <- n_hist + n_proj + lag       # 57

stopifnot(length(om_with_data_P$input$data$Ecov_obs) == n_total)

## Template for EM Ecov (copy from OM-simulated obs so EM "sees" this)
ecov_me <- list(
  label = "bt_temp",
  mean  = as.matrix(om_with_data_P$input$data$Ecov_obs),        # use observed series for EM template
  logsigma = "est_1",
  year  = years_all,
  use_obs = matrix(1, ncol = 1, nrow = n_total),
  process_model = "ar1",
  recruitment_how = matrix("controlling-lag-1-linear")
)

## 12.1 AR1 scenario for EM: identical to template (EM assumes AR1)
ecov_AR1 <- ecov_me

## 12.2 Historical Average scenario for EM: replace future mean with hist avg
hist_avg_val <- mean(ecov_me$mean[(1+lag):(n_hist+lag), 1], na.rm = TRUE)  # average over 1973–2022
ecov_HistAvg <- ecov_me
ecov_HistAvg$mean[(n_hist+lag+1):n_total, 1] <- hist_avg_val               # 2023–2028 rows

## Keep both in a named list
scenarios <- list(AR1 = ecov_AR1, HistAvg = ecov_HistAvg)

## Quick size checks
range(scenarios$AR1$year)   # 1972..2028
length(scenarios$AR1$year)  # 57

## ============================================================
## 13) Run two EM loops against the same OM:
##     - mod1: EM uses AR1 Ecov
##     - mod2: EM uses HistAvg Ecov in projection period (2023–2028)
## ============================================================

## Harvest control rule used in both runs
hcr <- list(hcr.type = 1, hcr.opts = list(use_FXSPR = TRUE, percentFXSPR = 75))

## For EM updates to index/catch info over full horizon,
## we reuse what we put into om_input_P (not input_P!)
update_index_info_list <- list(
  agg_index_sigma   = om_input_P$data$agg_index_sigma,
  index_Neff        = om_input_P$data$index_Neff,
  remove_agg        = TRUE,  remove_agg_pointer = 1:3, remove_agg_years = remove_agg_years1,
  remove_paa        = TRUE,  remove_paa_pointer = 1:3, remove_paa_years = remove_paa_years1
)

update_catch_info_list <- list(
  agg_catch_sigma = om_input_P$data$agg_catch_sigma,
  catch_Neff      = om_input_P$data$catch_Neff
)

## Common EM structural pieces (use same as OM to isolate Ecov differences)
M_em       <- M
sel_em     <- sel3
NAA_re_em  <- NAA_re
NAA_re_em$N1_model = "equilibrium" # IMPORTANT! Turn to equilibrium
age_comp_em <- "logistic-normal-pool0"

## ————— mod1: EM uses AR1 Ecov —————
set.seed(123 + 1L)
mod1 <- loop_through_fn(
  om = om_with_data_P,
  em_info = info,
  random  = random_P,
  M_em = M_em,
  sel_em = sel_em,
  NAA_re_em = NAA_re_em,
  ecov_em = scenarios$AR1,                       # <- AR1 (same structure as OM)
  age_comp_em = age_comp_em,
  em.opt = list(separate.em = FALSE, separate.em.type = 1, do.move = FALSE, est.move = FALSE),
  update_index_info = update_index_info_list,
  update_catch_info = update_catch_info_list,
  assess_years = assess.years,
  assess_interval = assess.interval,
  base_years = base.years,
  year.use = 50,
  add.years = TRUE,
  seed = 123 + 1L,
  hcr = hcr,
  save.sdrep = FALSE,
  save.last.em = TRUE,
  FXSPR_init = 0.5,
  do.brps = TRUE # IMPORTANT we calculate Reference point!
)
saveRDS(mod1, "mod1.RDS")

## ————— mod2: EM uses Historical Average in projection —————
## identify projection rows (2022..2028 indices in Ecov vector, we’ll pass period for EM option)
period_rows <- which(scenarios$AR1$year == 2022):which(scenarios$AR1$year == 2028)

set.seed(123 + 1L)
mod2 <- loop_through_fn(
  om = om_with_data_P,
  em_info = info,
  random  = random_P,
  M_em = M_em,
  sel_em = sel_em,
  NAA_re_em = NAA_re_em,
  ecov_em = scenarios$HistAvg,                   # <- EM will use HistAvg in proj
  ecov_em_opts = list(use_ecov_em = TRUE, period = period_rows, lag = 1),
  age_comp_em = age_comp_em,
  em.opt = list(separate.em = FALSE, separate.em.type = 1, do.move = FALSE, est.move = FALSE),
  update_index_info = update_index_info_list,
  update_catch_info = update_catch_info_list,
  assess_years = assess.years,
  assess_interval = assess.interval,
  base_years = base.years,
  year.use = 50,
  add.years = TRUE,
  seed = 123 + 1L,
  hcr = hcr,
  save.sdrep = FALSE,
  save.last.em = TRUE,
  FXSPR_init = 0.5,
  do.brps = TRUE # IMPORTANT we calculate Reference point!
)
saveRDS(mod2, "mod2.RDS")

## ============================================================
## 14) Small 2×2+legend Ecov diagnostic plots 
## ============================================================

layout(matrix(c(1,2,
                3,4,
                5,5), ncol = 2, byrow = TRUE), heights = c(1,1,0.3))
par(mar = c(4, 4, 2, 1))

## EM: observation used (HistAvg scenario)
plot(mod2$em_full[[1]]$input$data$Ecov_obs, type = "l", col = "orange",
     ylab = "Ecov", xlab = "Year", lwd = 2,
     main = "EM observation (input)")

## EM: true process trajectory (latent)
plot(mod2$em_full[[1]]$rep$Ecov_x, type = "l", col = "red", lwd = 2,
     ylab = "Ecov", xlab = "Year",
     main = "EM true process (estimated)")

## OM: observation (simulated with obs error)
plot(mod2$om$input$data$Ecov_obs, type = "o", col = "blue", lwd = 2, pch = 16,
     ylab = "Ecov", xlab = "Year",
     main = "OM observation (simulated)")

## OM: latent process (no obs error)
plot(mod2$om$rep$Ecov_x, type = "o", col = "green", lwd = 2, pch = 16,
     ylab = "Ecov", xlab = "Year",
     main = "OM true process")

## Legend
par(mar = c(0, 0, 0, 0))
plot.new()
legend("center",
       legend = c("EM observation (HistAvg for proj)",
                  "EM true process (given obs)",
                  "OM observation (AR1; with obs err)",
                  "OM true process (process, no obs err)"),
       col = c("orange", "red", "blue", "green"),
       lty = 1, lwd = 2, pch = c(NA, NA, 16, 16),
       horiz = FALSE, bty = "n", cex = 0.9)

## ---------- 14) Quick parameter checks & comparisons ----------
# Convenience back-transform for AR1 rho
back_trans_rho <- function(a) (exp(a) - 1) / (exp(a) + 1)

# Compare OM vs EM (mod1) parameter blocks
om <- mod1$om # OM that includes all "Truth"
em <- mod1$em_full[[1]] # Last EM in the feedback

# B-H SRR function has 2 pars.

# Rec par: mean log(B-H a) intercept
om$parList$mean_rec_pars[,1] # True value in the OM
em$parList$mean_rec_pars[,1] # Estimated value in the EM

# Rec par: B-H b
exp(om$parList$mean_rec_pars[,2]) # True value in the OM
exp(em$parList$mean_rec_pars[,2]) # Estimated value in the EM

# Catchability between true and estimate
om$parList$logit_q;        em$parList$logit_q

par(mar = c(2, 2, 2, 2), mfrow = c(1,1))
# N1 / NAA trajectories
plot(exp(om$parList$log_N1), main = "OM log_N1 (exp)") # OM True NAA in the first year
lines(em$rep$NAA[,,1,], col = "steelblue") # EM estimated NAA in the first year

# Fishing mortality pars
plot(exp(om$parList$F_pars), col = "red", main = "F params (exp)") # OM True F
lines(exp(em$parList$F_pars)) # EM est. F

# NAA sigmas
exp(om$parList$log_NAA_sigma) # OM true
exp(em$parList$log_NAA_sigma) # EM est

# Logistic-normal age-comp parameters
om$parList$catch_paa_pars # OM true
em$parList$catch_paa_pars # EM est

# Ecov (Temp) process pars: mean (state), logsigma, logit_rho
om$parList$Ecov_process_pars   # true OM
em$parList$Ecov_process_pars   # EM estimates

# True rho (AR1)
rho_om <- back_trans_rho(om$parList$Ecov_process_pars[3,])
rho_om

# Est rho (AR1)
rho_em <- back_trans_rho(em$parList$Ecov_process_pars[3,])
rho_em

# Observation error SD (Ecov)
exp(om$parList$Ecov_obs_logsigma)  # true
exp(em$parList$Ecov_obs_logsigma)  # EM estimate (often < true under ML)

# Recruitment Ecov beta
om$parList$Ecov_beta_R # OM true 
em$parList$Ecov_beta_R # EM est

# SSB (EM vs OM)
plot(em$rep$SSB, col = "red", main = "SSB: EM (red) vs OM (blue)")
lines(mod1$om$rep$SSB, col = "blue")

# Ecov states & RE
om$parList$Ecov_re
plot(em$parList$Ecov_re, type = "l", main = "EM Ecov RE")
plot(om$rep$Ecov_x, type = "o", main = "OM Ecov latent states")

# Numbers-at-age comparison for Age 2 & 3
plot(om$rep$NAA[,,,2], type = "l", col = "red", main = "NAA Age 2")
lines(em$rep$NAA[,,,2], col = "blue")
plot(om$rep$NAA[,,,3], type = "l", col = "red", main = "NAA Age 3")
lines(em$rep$NAA[,,,3], col = "blue")

# Catch 
plot(om$rep$pred_catch[1:53],type = "l") # Actual catch in the OM
lines(em$rep$pred_catch, col = "red") # EM est catch

# for the SAME OM, SAME SEED, IF you have 2 EMs (in this case, we have mod1, mod2), you can compare management performance between mod1 MP and mod2 MP, see below 

# If you want to compare Catch performance between 2 management strategies (only feedback years)
plot(mod1$om$rep$pred_catch[53:56], type = "l")
lines(mod2$om$rep$pred_catch[53:56], col = "red")

# If you want to compare SSB performance between 2 management strategies (only feedback years)
plot(mod1$om$rep$SSB[53:56,], type = "l")
lines(mod2$om$rep$SSB[53:56,], col = "red")

# If you want to compare F (at age max) performance between 2 management strategies (only feedback years)
plot(mod1$om$rep$Fbar[53:56,1], type = "l")
lines(mod2$om$rep$Fbar[53:56,1], col = "red")

# If you want to compare Overfished performance between 2 management strategies (only feedback years)
plot(mod1$om$rep$SSB/exp(mod1$om$rep$log_SSB_FXSPR_static[1]), type = "l") # Static BRP
plot(mod1$om$rep$SSB/exp(mod1$om$rep$log_SSB_FXSPR)[,1], type = "l") # Dynamic BRP

plot(mod2$om$rep$SSB/exp(mod2$om$rep$log_SSB_FXSPR_static[1]), type = "l") # Static BRP
plot(mod2$om$rep$SSB/exp(mod2$om$rep$log_SSB_FXSPR)[,1], type = "l") # Dynamic BRP

# If you want to compare Overfishing and Overfished performance between 2 management strategies (only feedback years)
plot(mod1$om$rep$Fbar[,1]/exp(mod1$om$rep$log_Fbar_XSPR)[,1], type = "l") # Static BRP
plot(mod1$om$rep$Fbar[,1]/exp(mod1$om$rep$log_Fbar_XSPR_static)[1], type = "l") # Dynamic BRP

plot(mod2$om$rep$Fbar[,1]/exp(mod2$om$rep$log_Fbar_XSPR)[,1], type = "l") # Static BRP
plot(mod2$om$rep$Fbar[,1]/exp(mod2$om$rep$log_Fbar_XSPR_static)[1], type = "l") # Dynamic BRP
