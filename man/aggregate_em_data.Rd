% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aggregate_em_data.R
\name{aggregate_em_data}
\alias{aggregate_em_data}
\title{Aggregate Catch and Index Data for Estimation Model}
\usage{
aggregate_em_data(
  data,
  info,
  aggregate_catch_info,
  aggregate_index_info,
  ind_em,
  n_fleets,
  n_indices
)
}
\arguments{
\item{data}{List. The operating model output, including:
\itemize{
  \item `$agg_catch` - Aggregated catch data.
  \item `$catch_paa` - Proportions-at-age for catch.
  \item `$agg_indices` - Aggregated index data.
  \item `$index_paa` - Proportions-at-age for indices.
}}

\item{aggregate_catch_info}{List. Contains fleet aggregation details, including:
\itemize{
  \item `$fleet_pointer` - Vector mapping original fleets to aggregated fleets (0 = exclude).
}}

\item{aggregate_index_info}{List. Contains index aggregation details, including:
\itemize{
  \item `$index_pointer` - Vector mapping original indices to aggregated indices (0 = exclude).
}}

\item{ind_em}{Vector. Indexing the estimation model years within the dataset.}

\item{n_fleets}{Integer. The number of fleets in the estimation model.}

\item{n_indices}{Integer. The number of indices in the estimation model.}

\item{em_info}{List. The estimation model inputs, including `par_inputs` and fleet/index settings.}
}
\value{
List. A modified `info` list containing aggregated `catch_info` and `index_info`:
  \itemize{
    \item `$catch_info$agg_catch` - Aggregated catch data for the estimation model.
    \item `$catch_info$catch_paa` - Aggregated catch-at-age data.
    \item `$index_info$agg_indices` - Aggregated index data.
    \item `$index_info$index_paa` - Aggregated index-at-age data.
  }
}
\description{
This function processes and aggregates catch and index data from the operating model
for use in the estimation model (`em_info`). It accounts for multiple fleets and indices
and adjusts aggregation based on fleet and index pointers.
}
