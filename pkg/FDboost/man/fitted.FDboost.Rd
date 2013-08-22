\name{fitted.FDboost}
\alias{fitted.FDboost}
\title{Fitted values of a boosted functional regression model}
\usage{
  \method{fitted}{FDboost} (object, ...)
}
\arguments{
  \item{object}{a fitted \code{FDboost}-object}

  \item{...}{additional arguments passed on to
  \code{\link{predict.FDboost}}}
}
\value{
  matrix of fitted values
}
\description{
  Takes a fitted \code{FDboost}-object and computes the
  fitted values.
}
\seealso{
  \code{\link{FDboost}} for the model fit.
}

