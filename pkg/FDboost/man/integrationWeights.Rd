\name{integrationWeights}
\alias{integrationWeights}
\alias{integrationWeightsLeft}
\title{Functions to compute integration weights}
\usage{
  integrationWeights(X1, xind)
}
\arguments{
  \item{X1}{matrix of functional variable}

  \item{xind}{index of functional variable}
}
\description{
  Computes trapezoidal integration weights for a functional
  variable X1 on grid xind.
}
\details{
  The function \code{integrationWeights()} computes
  trapezoidal integraion weights, that are symmetric. The
  function \code{integrationWeightsLeft()} computes
  weights, that take into account soley the distance to the
  prior observation point.  Those weights are adequate for
  historical effects.
}

