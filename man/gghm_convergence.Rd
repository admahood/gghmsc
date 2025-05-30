% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/functions.R
\name{gghm_convergence}
\alias{gghm_convergence}
\title{Plot model convergence indexes}
\usage{
gghm_convergence(
  Hm,
  beta = TRUE,
  V = FALSE,
  gamma = FALSE,
  omega = FALSE,
  title = "Model Convergence"
)
}
\arguments{
\item{beta}{if TRUE, plots the beta (env. filters) parameters}

\item{V}{if TRUE, plots the V parameters}

\item{gamma}{if TRUE, plots the gamma (traits) paramters}

\item{omega}{if TRUE, plots the omega (spp associations) parameters}

\item{title}{character string to customize}
}
\description{
Plot model convergence indexes
}
