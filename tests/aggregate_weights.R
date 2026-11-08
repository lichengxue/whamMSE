library(whamMSE)

make_fixture <- function() {
  ny <- 4L
  na <- 3L
  catches <- rbind(c(1, 2, 3, 4), c(4, 3, 2, 1),
                   c(2, 1, 3, 4), c(5, 1, 1, 5))
  indices <- rbind(c(1, 3), c(4, 1), c(2, 2), c(3, 1))
  waa <- array(0, c(8, ny, na))
  for (block in 1:8) for (age in 1:na) {
    waa[block, , age] <- 10 * block + 100 * (age - 1)
  }
  maturity <- array(0, c(2, ny, na))
  maturity[1, , ] <- rep(c(0.2, 0.4, 1), each = ny)
  maturity[2, , ] <- rep(c(0.4, 0.8, 1), each = ny)
  par <- list(
    n_fleets = 4, n_indices = 2, n_regions = 2, n_stocks = 2,
    n_ages = na, catch_cv = rep(0.1, 4), catch_Neff = rep(50, 4),
    use_agg_catch = rep(1, 4), use_catch_paa = rep(0, 4),
    index_cv = rep(0.4, 2), index_Neff = rep(50, 2), q = rep(1, 2),
    use_indices = rep(1, 2), use_index_paa = rep(0, 2),
    units_indices = rep(2, 2), units_index_paa = rep(2, 2),
    fracyr_indices = rep(0.25, 2),
    user_waa = list(waa = waa, waa_pointer_fleets = 1:4,
                    waa_pointer_indices = 5:6, waa_pointer_ssb = 7:8,
                    waa_pointer_M = 7:8),
    user_maturity = maturity
  )
  list(
    om = list(input = list(data = list(
      n_fleets = 4, n_indices = 2, n_regions = 2, n_stocks = 2,
      agg_catch = catches, agg_indices = indices
    ))),
    em = list(par_inputs = par, basic_info = list())
  )
}

aggregate_fixture <- function(x, years = 1:4, pointer = 1, fleet = FALSE,
                              weighted = TRUE) {
  rule <- list(fleet = fleet, index = !fleet, pointer = pointer)
  whamMSE::make_aggregate_data(
    x$om, x$em, years,
    aggregate_catch_info = list(
      n_fleets = 2, fleet_pointer = c(1, 2, 1, 2),
      use_catch_weighted_waa = weighted
    ),
    aggregate_index_info = list(
      n_indices = 1, index_pointer = c(1, 1),
      use_catch_weighted_waa = weighted
    ),
    aggregate_weights_info = list(ssb_waa_weights = rule,
                                  maturity_weights = rule)
  )
}

equal_numeric <- function(x, y) {
  stopifnot(isTRUE(all.equal(unname(x), unname(y), tolerance = 1e-12)))
}

x <- make_fixture()
answer <- aggregate_fixture(x)
w <- answer$par_inputs$user_waa$waa
m <- answer$par_inputs$user_maturity

# Independent expected annual values for each fleet, survey, and population.
equal_numeric(w[1, , 1], c(25, 50/3, 22, 40/3))
equal_numeric(w[2, , 1], c(100/3, 25, 36, 110/3))
equal_numeric(w[3, , 1], c(57.5, 52, 55, 52.5))
equal_numeric(w[4, , 1], c(77.5, 72, 75, 72.5))
equal_numeric(m[1, , 2], c(0.7, 0.48, 0.6, 0.5))
equal_numeric(m[1, , 3], rep(1, 4))
stopifnot(all(m >= 0 & m <= 1))

# A complete vector of original members is equivalent to its aggregate group.
vector_answer <- aggregate_fixture(x, pointer = c(1, 2))
reverse_answer <- aggregate_fixture(x, pointer = c(2, 1))
stopifnot(identical(answer, vector_answer), identical(answer, reverse_answer))
# Scalar fleet groups and their complete original-member vectors both work.
commercial <- aggregate_fixture(x, fleet = TRUE)
commercial_members <- aggregate_fixture(x, fleet = TRUE, pointer = c(1, 3))
stopifnot(identical(commercial, commercial_members))
equal_numeric(commercial$par_inputs$user_waa$waa[4, , 1], c(77.5, 220/3, 76, 215/3))
equal_numeric(commercial$par_inputs$user_maturity[1, , 2],
              c(0.7, 8/15, 0.64, 7/15))

# Adding assessment years preserves earlier rows and uses new observations.
short <- aggregate_fixture(x, years = 1:2, pointer = c(1, 2))
equal_numeric(short$par_inputs$user_waa$waa, w[, 1:2, , drop = FALSE])
equal_numeric(short$par_inputs$user_maturity, m[, 1:2, , drop = FALSE])
shifted <- aggregate_fixture(x, years = 2:4)
equal_numeric(shifted$par_inputs$user_waa$waa, w[, 2:4, , drop = FALSE])
equal_numeric(shifted$par_inputs$user_maturity, m[, 2:4, , drop = FALSE])
changed <- x
changed$om$input$data$agg_indices[4, ] <- c(1, 3)
updated <- aggregate_fixture(changed, pointer = c(1, 2))
equal_numeric(updated$par_inputs$user_waa$waa[4, 4, 1], 77.5)
equal_numeric(updated$par_inputs$user_maturity[1, 4, 2], 0.7)
equal_numeric(updated$par_inputs$user_maturity[, 1:3, , drop = FALSE],
              m[, 1:3, , drop = FALSE])

# Ambiguous or invalid vectors are rejected instead of selecting a list entry.
for (pointer in list(c(1, 1), c(1, 3), c(1, NA), c(1, 1.5))) {
  result <- try(aggregate_fixture(x, pointer = pointer), silent = TRUE)
  stopifnot(inherits(result, "try-error"))
}
cross_group <- try(aggregate_fixture(x, fleet = TRUE, pointer = c(1, 2)),
                   silent = TRUE)
stopifnot(inherits(cross_group, "try-error"))

# Changes to recent catch observations also update catch-weighted inputs.
changed_catch <- x
changed_catch$om$input$data$agg_catch[4, c(1, 3)] <- c(1, 3)
updated_catch <- aggregate_fixture(changed_catch, fleet = TRUE, pointer = c(1, 3))
equal_numeric(updated_catch$par_inputs$user_waa$waa[1, 4, 1], 25)
equal_numeric(updated_catch$par_inputs$user_waa$waa[4, 4, 1], 77.5)
equal_numeric(updated_catch$par_inputs$user_maturity[1, 4, 2], 0.7)
equal_numeric(updated_catch$par_inputs$user_waa$waa[, 1:3, , drop = FALSE],
              commercial$par_inputs$user_waa$waa[, 1:3, , drop = FALSE])
