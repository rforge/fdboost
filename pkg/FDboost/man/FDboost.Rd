\name{FDboost}
\alias{FDboost}
\title{Model-based Gradient Boosting for Functional Response}
\usage{
  FDboost(formula, timeformula, numInt = "equal", data,
    weights = NULL, offset = NULL,
    offset_control = o_control(), ...)
}
\arguments{
  \item{formula}{a symbolic description of the model to be
  fit.}

  \item{timeformula}{formula for the expansion over the
  index of the response. For a functional response
  \eqn{Y_i(t)} typically ~bbs(t) to obtain a smooth
  expansion of the effects along t.}

  \item{numInt}{integration scheme for the integration of
  the loss function. One of \code{c("equal", "Riemann")}
  meaning equal weights of 1 or trapezoidal Riemann
  weights. Alternativley a vector of length
  \code{nrow(response)} containing arbitray positive
  weights can be specified.}

  \item{data}{a data frame or list containing the variables
  in the model.}

  \item{weights}{(1) a numeric vector of weights for
  observational units, i.e. length(weights) has to be
  nrow(response), (2) alternatively weights can be
  specified for single observations then length(weights)
  has to be \code{nrow(response)}*\code{ncol(response)} per
  default weights is constantly 1.}

  \item{offset_control}{parameters for the calculation of
  the offset, defaults to \code{offset_control(k_min=20,
  silent=TRUE)}.}

  \item{offset}{a numeric vector to be used as offset over
  the index of the response (optional). If no offset is
  specified, per default a smooth time-specific offset is
  calculated and used within the model fit. If you do not
  want to use an offset you can set offset=0.}

  \item{...}{additional arguments passed to
  \code{\link[mboost]{mboost}}, including \code{offset},
  \code{family} and \code{control}.}
}
\value{
  on object of class \code{FDboost} that inherits from
  \code{mboost}. Special \code{\link{predict.FDboost}},
  \code{\link{coef.FDboost}} and \code{\link{plot.FDboost}}
  methods are available. The methods of
  \code{\link[mboost]{mboost}} ara available as well, e.g.
  \code{\link[mboost]{extract}}
}
\description{
  Gradient boosting for optimizing arbitrary loss
  functions, where component-wise models are utilized as
  base-learners in the case of functional response. This
  function is a wrapper for \code{mboost}'s
  \code{\link[mboost]{mboost}} and its siblings to fit
  models of the general form \cr \eqn{E(Y_i(t)) = g(\mu(t)
  + \int X_i(s)\beta(s,t)ds + f(z_{1i}, t) + f(z_{2i}) +
  z_{3i} \beta_3(t) + \dots }\cr with a functional (but not
  necessarily continuous) response \eqn{Y(t)}, response
  function \eqn{g}, (optional) smooth intercept
  \eqn{\mu(t)}, (multiple) functional covariates \eqn{X(t)}
  and scalar covariates \eqn{z_1}, \eqn{z_2}, etc.
}
\details{
  The functional response and functional covariates have to
  be supplied as n by <no. of evaluations> matrices, i.e.
  each row is one functional observation. For the model fit
  the matrix of the functional response evaluations
  \eqn{Y_i(t)} are stacked into one long vector.
}
\examples{
###############################################################################
print("to do")
###############################################################################
}
\author{
  Sarah Brockhaus, Torsten Hothorn
}
\references{
  Ivanescu, A., Staicu, A.-M., Scheipl, F. and Greven, S.
  (2012).  Penalized function-on-function regression.
  (submitted)
  \url{http://biostats.bepress.com/jhubiostat/paper240/}
}
\seealso{
  \code{\link[mboost]{mboost}} for the help of the wrapped
  function in package mboost.
}
\keyword{models,}
\keyword{nonlinear}

