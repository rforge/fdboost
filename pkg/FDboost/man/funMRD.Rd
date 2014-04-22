\name{funMRD}
\alias{funMRD}
\title{Functional MRD}
\usage{
  funMRD(object, overTime = TRUE, breaks = object$yind,
    global = FALSE, ...)
}
\arguments{
  \item{object}{fitted FDboost-object with regular
  response}

  \item{overTime}{per default the functional MRD is
  calculated over time if \code{overTime=FALSE}, the MRD is
  calculated per curve}

  \item{breaks}{an optional vector or number giving the
  time-points at which the model is evaluated. Can be
  specified as number of equidistant time-points or as
  vector of time-points. Defaults to the index of the
  response in the model.}

  \item{global}{logical. defaults to \code{FALSE}, if TRUE
  the global MRD like in a normal linear model is
  calculated}

  \item{...}{currently not used}
}
\value{
  Returns a vector with the calculated MRD and some extra
  information in attributes.
}
\description{
  Calculates the functional MRD for a fitted FDboost-object
}
\details{
  Formula to calculate MRD over time, \code{overTime=TRUE}:
  \cr \eqn{ MRD(t) = n^{-1} \sum_i |(Y_i(t) -
  \hat{Y}_i(t))^2|/|Y_i(t)| }

  Formula to calculate MRD over subjects,
  \code{overTime=FALSE}: \cr \eqn{ MRD_i = \int |(Y_i(t) -
  \hat{Y}_i(t))^2|/|Y_i(t)| dt \approx G^{-1} \sum_g
  |(Y_i(t_g) - \hat{Y}_i(t_g))^2| / |Y_i(t)|}
}
\note{
  \code{breaks} cannot be changed in the case the
  \code{bsignal()} is used over the same domain as the
  response! In that case you would have to rename the index
  of the response or that of the covariates.
}

