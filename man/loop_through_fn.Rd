% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loop_through_fn.R
\name{loop_through_fn}
\alias{loop_through_fn}
\title{Perform Management Strategy Evaluation (MSE)}
\usage{
loop_through_fn(
  om,
  em_info = NULL,
  random = NULL,
  M_em = NULL,
  sel_em = NULL,
  NAA_re_em = NULL,
  move_em = NULL,
  age_comp_em = NULL,
  em.opt = list(separate.em = TRUE, separate.em.type = 1, do.move = FALSE, est.move =
    FALSE),
  global_waa = TRUE,
  aggregate_catch_info = NULL,
  aggregate_index_info = NULL,
  filter_indices = NULL,
  assess_years = NULL,
  assess_interval = NULL,
  base_years = NULL,
  year.use = 20,
  add.years = FALSE,
  hcr = list(hcr.type = 1, hcr.opts = NULL),
  catch_alloc = list(weight_type = 1, method = "equal", user_weights = NULL, weight_years
    = 1),
  do.retro = FALSE,
  do.osa = FALSE,
  do.brps = FALSE,
  seed = 123,
  save.sdrep = FALSE,
  save.last.em = FALSE,
  by_fleet = TRUE,
  FXSPR_init = NULL
)
}
\arguments{
\item{om}{A fitted operating model containing pseudo data.}

\item{em_info}{A list containing information used to generate the estimation model (default = NULL).}

\item{random}{Character vector of processes treated as random effects in the operating model (default = "log_NAA").}

\item{M_em}{Configuration of natural mortality in the estimation model.}

\item{sel_em}{Configuration of selectivity in the estimation model.}

\item{NAA_re_em}{Configuration of numbers-at-age (NAA) in the estimation model.}

\item{move_em}{Configuration of movement in the estimation model.}

\item{age_comp_em}{Character. The likelihood distribution of age composition data in the estimation model. Options include:
\itemize{
  \item \code{"multinomial"} (default)
  \item \code{"dir-mult"}
  \item \code{"dirichlet-miss0"}
  \item \code{"dirichlet-pool0"}
  \item \code{"logistic-normal-miss0"}
  \item \code{"logistic-normal-ar1-miss0"}
  \item \code{"logistic-normal-pool0"}
  \item \code{"logistic-normal-01-infl"}
  \item \code{"logistic-normal-01-infl-2par"}
  \item \code{"mvtweedie"}
  \item \code{"dir-mult-linear"}
}}

\item{em.opt}{List of options for the assessment model, including:
\describe{
  \item{\code{separate.em}}{Logical. If TRUE, separate estimation models are used for different stocks (default = FALSE).}
  \item{\code{separate.em.type}}{Integer. Used only if \code{separate.em = TRUE}:}
    \itemize{
      \item 1 = Panmictic (spatially-aggregated model).
      \item 2 = Fleets-as-areas model.
      \item 3 = Separate assessment models for each region (n = number of regions).
    }
  \item{\code{do.move}}{Logical. If TRUE, movement is included (default = FALSE).}
  \item{\code{est.move}}{Logical. If TRUE, movement rates are estimated (default = FALSE).}
}}

\item{aggregate_catch_info}{List (optional). User-specified catch aggregation settings 
for panmictic models using aggregate catch.
\itemize{
  \item `$n_fleets` Integer. Number of fleets.
  \item `$catch_cv` Numeric vector (`n_fleets`). Catch CVs for each fleet.
  \item `$catch_Neff` Numeric vector (`n_fleets`). Effective sample sizes for catch.
  \item `$use_agg_catch` Integer vector (`n_fleets`). 0/1 flag for using aggregate catch.
  \item `$use_catch_paa` Integer vector (`n_fleets`). 0/1 flag for using age composition.
  \item `$fleet_pointer` Integer vector (`n_fleets`). Defines fleet grouping (0 = exclude).
  \item `$use_catch_weighted_waa` Logical. Whether to weight-at-age using catch.
}}

\item{aggregate_index_info}{List (optional). User-specified index aggregation settings for panmictic models using aggregate indices.
\itemize{
  \item `$n_indices` Integer. Number of indices.
  \item `$index_cv` Numeric vector (`n_indices`). CVs for survey indices.
  \item `$index_Neff` Numeric vector (`n_indices`). Effective sample sizes for indices.
  \item `$fracyr_indices` Numeric vector (`n_indices`). Fraction of the year for each survey.
  \item `$q` Numeric vector (`n_indices`). Survey catchability coefficients.
  \item `$use_indices` Integer vector (`n_indices`). 0/1 flag for using survey indices.
  \item `$use_index_paa` Integer vector (`n_indices`). 0/1 flag for using age composition.
  \item `$units_indices` Integer vector (`n_indices`). 1 = Biomass, 2 = Numbers.
  \item `$units_index_paa` Integer vector (`n_indices`). 1 = Biomass, 2 = Numbers.
  \item `$index_pointer` Integer vector (`n_indices`). Defines index grouping (0 = exclude).
  \item `$use_index_weighted_waa` Logical. Whether to weight-at-age using indices.
   }}

\item{filter_indices}{Integer (0/1) vector (optional). User-specified which indices are excluded from the assessment model. For example, c(1,0,1,1) indicates Index 1 (include) and 2 (exclude) in region 1, Index 3 (include) and 4 (include) in region 2}

\item{assess_years}{Vector of years when assessments are conducted.}

\item{assess_interval}{Integer. The interval between stock assessments in the MSE feedback loop.}

\item{base_years}{Vector of years used in the burn-in period.}

\item{year.use}{Integer. Number of years included in the estimation model (default = 20).}

\item{add.years}{Logical. Whether or not using entire time series of data in the assessment model (default = FALSE).}

\item{hcr}{List containing harvest control rule (HCR) settings:
\describe{
  \item{\code{hcr.type}}{Integer. The type of harvest control rule:}
    \itemize{
      \item \code{1} Annual projected catch based on 75\% of F40\% (default).
      \item \code{2} Constant catch based on 75\% of F40\%.
      \item \code{3} "Hockey stick" catch based on stock status.
    }
  \item{\code{hcr.opts}}{(Only for \code{hcr.type = 3}) A list containing additional HCR options:}
    \describe{
      \item{\code{max_percent}}{Maximum percentage of F\_XSPR used for catch projections (default = 75).}
      \item{\code{min_percent}}{Minimum percentage of F\_XSPR used for catch projections (default = 0.01).}
      \item{\code{BThresh_up}}{Upper bound of overfished biomass threshold (default = 0.5).}
      \item{\code{BThresh_low}}{Lower bound of overfished biomass threshold (default = 0.1).}
    }
}}

\item{catch_alloc}{List. Contains specifications for catch allocation:
\itemize{
  \item `$weight_type` - Integer (1–4) indicating the type of weighting used.
    \itemize{
      \item `1` - Equal weighting across fleets.
      \item `2` - Weighting based on past catch data.
      \item `3` - Weighting based on past survey index data.
      \item `4` - User-specified weights.
    }
  \item `$method` - String. Specifies the catch allocation method:
    \itemize{
      \item `"equal"` - Use equal weights, split total catch into fleet-specific catches (use when weight_type = 1).
      \item `"fleet_equal"` - Use equal weights, split total fleet catch into fleet-specific catches (use when weight_type = 1).
      \item `"fleet_region"` - Uses regional catch totals to compute weights (use when weight_type = 2).
      \item `"fleet_gear"` - Uses gear-type catch totals to compute weights (use when weight_type = 2).
      \item `"fleet_combined"` - First allocates to gear-specific total catch, then splits into regions (use when weight_type = 2).
      \item `"fleet_catch"` - Uses fleet-specific catch to compute weights (use when weight_type = 2).
      \item `"index_equal"` - Uses survey-based regional weighting, then assigns equally among fleets in the same region (use when weight_type = 3).
      \item `"index_gear"` - Uses survey-based regional weighting, then assigns based on gear-specific weights (use when weight_type = 3).
      \item `"index_multiple"` - Uses multiple survey-based regional weighting, then assigns equally among fleets in the same region (use when weight_type = 3).
      \item `"user_defined"` - Uses manually specified weights for each region or each fleet (use when weight_type = 4).
    }
  \item `$user_weights` - Numeric vector (`n_regions` or `n_fleets`). Optional. User-defined weights summing to 1 (use when weight_type = 4).
  \item `$weight_years` - Integer. Number of years to average for calculating historical catch weights.
  \item `$survey_pointer` - Integer/Vector. Specifies which survey index type to use for weighting (use when weight_type = 3).
}}

\item{do.retro}{Logical. If TRUE, performs retrospective analysis for the assessment model (default = TRUE).}

\item{do.osa}{Logical. If TRUE, calculates one-step-ahead (OSA) residuals for the assessment model (default = TRUE).}

\item{do.brps}{Logical. If TRUE, calculates reference points in the operating model (default = FALSE).}

\item{seed}{Integer. The random seed used for stochastic processes.}

\item{save.sdrep}{Logical. If TRUE, saves the results of every assessment model (memory-intensive).}

\item{save.last.em}{Logical. If TRUE, saves only the last estimation model (default = FALSE).}

\item{by_fleet}{Logical. If TRUE, estimates F separately for each fleet. If FALSE, estimates a single global F (default = TRUE).}

\item{FXSPR_init}{Numeric. change initial F for estimating reference point in the assessment model (default = NULL).}

\item{ind_em}{Vector. Indices specifying the years for which the estimation model should use data.}
}
\value{
A list containing:
  \describe{
    \item{\code{om}}{Updated operating model.}
    \item{\code{em_list}}{List of estimation model results across assessment years.}
    \item{\code{par.est}}{Estimated parameter values.}
    \item{\code{par.se}}{Standard errors of parameter estimates.}
    \item{\code{adrep.est}}{Estimated values from ADMB report.}
    \item{\code{adrep.se}}{Standard errors from ADMB report.}
    \item{\code{opt_list}}{Optimization results from WHAM.}
    \item{\code{converge_list}}{Convergence diagnostics.}
    \item{\code{catch_advice}}{Projected catch advice over assessment years.}
    \item{\code{em_full}}{List of full estimation model outputs.}
    \item{\code{runtime}}{Elapsed time for function execution.}
  }
}
\description{
A wrapper function to iterate through the feedback period in a management strategy evaluation (MSE), 
updating the operating model (OM), refitting the estimation model (EM), and generating catch advice.
}
\seealso{
\code{\link{make_em_input}}, \code{\link{update_om_fn}}, \code{\link{advice_fn}}
}
