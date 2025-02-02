% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loop_through_fn.R
\name{loop_through_fn}
\alias{loop_through_fn}
\title{Function to perform management strategy evaluation}
\usage{
loop_through_fn(
  om,
  em_info = NULL,
  random = "log_NAA",
  M_em,
  sel_em,
  NAA_re_em,
  move_em,
  age_comp_em = "multinomial",
  em.opt = NULL,
  aggregate_catch_info = NULL,
  aggregate_index_info = NULL,
  assess_years = NULL,
  assess_interval = NULL,
  base_years = NULL,
  year.use = 30,
  hcr.type = 1,
  hcr.opts = NULL,
  weight_type = 1,
  weight_years = 1,
  user_weights = NULL,
  do.retro = FALSE,
  do.osa = FALSE,
  seed = 123,
  save.sdrep = FALSE,
  save.last.em = FALSE
)
}
\arguments{
\item{om}{Operating model including pseudo data}

\item{em_info}{A list of information used to generate the estimation model (default = NULL)}

\item{random}{A vector of processes that are treated as random effects in the operating model (default = "log_NAA")}

\item{M_em}{Natural mortality configuration in the assessment model}

\item{sel_em}{Selectivity configuration in the assessment model}

\item{NAA_re_em}{Numbers-at-age (NAA) configuration in the assessment model}

\item{move_em}{Movement configuration in the assessment model}

\item{age_comp_em}{Likelihood distribution of age composition data in the assessment model
\itemize{
  \item \code{"multinomial"} Default
  \item \code{"dir-mult"}
  \item \code{"dirichlet-miss0"}
  \item \code{"dirichlet-pool0"}
  \item \code{"dir-mult"}
  \item \code{"logistic-normal-miss0"}
  \item \code{"logistic-normal-ar1-miss0"}
  \item \code{"logistic-normal-pool0"}
  \item \code{"logistic-normal-01-infl"}
  \item \code{"logistic-normal-01-infl-2par"}
  \item \code{"mvtweedie"}
  \item \code{"dir-mult-linear"}
}}

\item{em.opt}{List of options for the assessment model
\itemize{
  \item \code{$separate.em} TRUE = No Global SPR, FALSE = Global SPR
  \item \code{$separate.em.type} only if separate.em = TRUE \cr
  {=1} panmictic (spatially-aggregated) \cr
  {=2} fleets-as-areas \cr
  {=3} n single assessment models (n = n_regions) \cr
  \item \code{$do.move} T/F movement is included (use if separate.em = FALSE)
  \item \code{$est.move} T/F movement rate is estimated (use if separate.em = FALSE)
}}

\item{aggregate_catch_info}{(optional) User specified list of infomation for aggregate catch (default=NULL: using the infomation from the first fleet). Only when pamictic model with aggregate catch is used
\itemize{
  \item \code{$catch_cv} cv for the catch (can be a single value or a vector of length n_fleets)
  \item \code{$catch_Neff} effective sample size for the catch (can be a single value or a vector of length n_fleets)
  }}

\item{aggregate_index_info}{(optional) User specified list of infomation for aggregate index (default=NULL: using the infomation from the first index). Only when pamictic model with aggregate index is used
\itemize{
  \item \code{$index_cv} cv for the indices (can be a single value or a vector of length n_indices)
  \item \code{$index_Neff} effective sample size for the indices (can be a single value or a vector of length n_indices)
  \item \code{$fracyr_indices} fraction of the year when survey is conducted (can be a single value or a vector of length n_indices) 
  \item \code{$q} survey catchability (can be a single value or a vector of length n_indices)
  }}

\item{assess_years}{Year in which the assessment is conducted}

\item{assess_interval}{Assessment interval used in the MSE feedback loop}

\item{base_years}{Years used in the burn-in period}

\item{year.use}{Number of years used in the assessment model}

\item{hcr.type}{Type of harvest control rule
\itemize{
  \item \code{"1"} Annual projected catch based on 75% of F40% (default)
  \item \code{"2"} Constant catch based on 75% of F40%
  \item \code{"3"} "Hockey stick" catch based on stock status
}}

\item{hcr.opts}{only used if hcr.type = 3
\itemize{
  \item \code{"max_percent"} maximum percent of F_XSPR to use for calculating catch in projections (default = 75)
  \item \code{"min_percent"} minimum percent of F_XSPR to use for calculating catch in projections (default = 0.01)
  \item \code{"BThresh_up"} Upper bound of overfished level (default = 0.5)  
  \item \code{"BThresh_low"} Lower bound of overfished level (default = 0.1)
}}

\item{weight_type}{Type of weights to use for allocating catch
\itemize{
  \item \code{"1"} Uniform: Assign equal weights to all regions
  \item \code{"2"} User defined: User provides weights (provide weights as \code{user_weights})
  \item \code{"3"} Catch history: Use catch data averaged over specified number of years (provide years as \code{weight_years})
  \item \code{"4"} Survey average: Use survey data averaged over specified number of years (provide years as \code{weight_years})
  \item \code{"5"} Recruitment proportional: Use average recruitment estimates to allocate weights proportionally
}}

\item{weight_years}{Number of years to average survey or catch data for weights (only if \code{weight_type = 3 or 4})}

\item{user_weights}{User-defined weights for allocating catch (only if \code{weight_type = 2})}

\item{do.retro}{T/F Do retrospective analysis? Default = TRUE}

\item{do.osa}{T/F Calculate one-step-ahead (OSA) residuals? Default = TRUE}

\item{seed}{Seed used for generating data}

\item{save.sdrep}{T/F Save every assessment model (memory intensive)
\itemize{
  \item \code{"TRUE"} Save every assessment model  
  \item \code{"FALSE"} Save the last assessment model (default)
}}

\item{save.last.em}{T/F Save the last assessment model (Default = FALSE) (use if separate.em = FALSE)}
}
\value{
a list of model output
}
\description{
A wrapper function to step through the feedback period, update the operating model,
refit the estimation model, and generate catch advice.
}
\seealso{
\code{\link{make_em_input}}, \code{\link{update_om_fn}}, \code{\link{advice_fn}}
}
