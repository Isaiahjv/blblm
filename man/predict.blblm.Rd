% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/blblm.R
\name{predict.blblm}
\alias{predict.blblm}
\title{Returns predictions from fitted boostrap model on test data (data unseen by the fitted model). Parallel computing os optional.}
\usage{
\method{predict}{blblm}(
  object,
  new_data,
  confidence = FALSE,
  level = 0.95,
  para = FALSE,
  core = 2,
  ...
)
}
\arguments{
\item{object}{linear model}

\item{new_data}{dataset used for prediction}

\item{confidence}{logical indicating confidence interval usage}

\item{level}{numeric indicating level of confidence. see previous.}

\item{para}{logical indicating parallel cimouting usage}

\item{core}{integer indicating number of core in parallel computing}

\item{...}{so on}
}
\description{
Returns predictions from fitted boostrap model on test data (data unseen by the fitted model). Parallel computing os optional.
}
