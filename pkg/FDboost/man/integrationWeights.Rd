\name{integrationWeights}
\alias{integrationWeights}
\alias{integrationWeightsLeft}
\title{Functions to compute integration weights}
\usage{
  integrationWeights(X1, xind, id = NULL)
}
\arguments{
  \item{X1}{matrix of functional variable}

  \item{xind}{index of functional variable}

  \item{id}{defaults to NULL if X1 is a matrix. identity
  variable if X1 is in long format.}
}
\description{
  Computes trapezoidal integration weights for a functional
  variable X1 on grid xind.
}
\details{
  The function \code{integrationWeights()} computes
  trapezoidal integration weights, that are symmetric. The
  function \code{integrationWeightsLeft()} computes
  weights, that take into account only the distance to the
  prior observation point.  Those weights are adequate for
  historical effects.
}

