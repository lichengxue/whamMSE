% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_F_from_Catch.R
\name{get_F_from_Catch}
\alias{get_F_from_Catch}
\title{Estimate Fleet-Specific or Global Fishing Mortality (F) from Catch}
\usage{
get_F_from_Catch(om, Catch, year, method = "nlminb", by_fleet = FALSE)
}
\arguments{
\item{om}{A fitted operating model containing stock and fleet information.}

\item{Catch}{Numeric vector of catch values for each fleet.}

\item{year}{Integer. The year for which F is estimated.}

\item{method}{Character. The optimization method used for solving F. Options include:
\itemize{
  \item \code{"nlminb"} (default): Uses `nlminb` optimization.
  \item \code{"BFGS"}: Uses `optim` with the BFGS method.
}}

\item{by_fleet}{Logical. If TRUE, estimates F separately for each fleet. If FALSE, estimates a single global F (default = FALSE).}
}
\value{
A numeric vector of estimated fleet-specific F values.
}
\description{
This function estimates fleet-specific or global fishing mortality (F) using Newton's method.
}
