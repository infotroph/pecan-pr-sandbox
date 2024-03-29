% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pda.bayesian.tools.R
\name{pda.bayesian.tools}
\alias{pda.bayesian.tools}
\title{Paramater Data Assimilation using BayesianTools}
\usage{
pda.bayesian.tools(
  settings,
  external.data = NULL,
  external.priors = NULL,
  external.formats = NULL,
  ensemble.id = NULL,
  params.id = NULL,
  param.names = NULL,
  prior.id = NULL,
  chain = NULL,
  iter = NULL,
  adapt = NULL,
  adj.min = NULL,
  ar.target = NULL,
  jvar = NULL,
  remote = FALSE,
  ...
)
}
\arguments{
\item{settings}{= a pecan settings list}

\item{external.data}{list of external inputs}

\item{external.priors}{list of external priors}

\item{external.formats}{bety formats used when function is used without a DB connection, e.g. remote}

\item{ensemble.id}{ensemble IDs}

\item{params.id}{id of pars}

\item{param.names}{names of pars}

\item{prior.id}{ids of priors}

\item{chain}{how many chains}

\item{iter}{how many iterations}

\item{adapt}{adaptation intervals}

\item{adj.min}{to be used in adjustment}

\item{ar.target}{acceptance rate target}

\item{jvar}{jump variance}

\item{remote}{logical, if TRUE no DB connection is established}

\item{...}{additional arguments}
}
\value{
settings
}
\description{
Paramater Data Assimilation using BayesianTools R Package
}
\author{
Istem Fer
}
