\name{FDboost}
\alias{FDboost}
\title{Model-based Gradient Boosting for Functional Response}
\usage{
  FDboost(formula, timeformula, id = NULL,
    numInt = "equal", data, weights = NULL, offset = NULL,
    offset_control = o_control(), ...)
}
\arguments{
  \item{formula}{a symbolic description of the model to be
  fit.}

  \item{timeformula}{formula for the expansion over the
  index of the response. For a functional response
  \eqn{Y_i(t)} typically ~bbs(t) to obtain a smooth
  expansion of the effects along t. In the limiting case
  that \eqn{Y_i} is a scalar response use \code{~bols(1)},
  which sets up a base-learner for the scalar 1.}

  \item{id}{defaults to NULL which means that the response
  is a matrix with a regular time. If the response is given
  in long format for irregular observations, \code{id}
  contains the information which observations belong
  together. id should contain numbers 1, 2, 3, ...}

  \item{numInt}{integration scheme for the integration of
  the loss function. One of \code{c("equal", "Riemann")}
  meaning equal weights of 1 or trapezoidal Riemann
  weights. Alternatively a vector of length
  \code{nrow(response)} containing arbitrary positive
  weights can be specified.}

  \item{data}{a data frame or list containing the variables
  in the model.}

  \item{weights}{(1) a numeric vector of weights for
  observational units, i.e. \code{length(weights)} has to
  be \code{nrow(response)}, (2) alternatively weights can
  be specified for single observations then
  \code{length(weights)} has to be
  \code{nrow(response)}*\code{ncol(response)} per default
  weights is constantly 1.}

  \item{offset_control}{parameters for the calculation of
  the offset, defaults to \code{offset_control(k_min=20,
  silent=TRUE)}.}

  \item{offset}{a numeric vector to be used as offset over
  the index of the response (optional). If no offset is
  specified, per default a smooth time-specific offset is
  calculated and used within the model fit. If you do not
  want to use an offset you can set \code{offset=0}.}

  \item{...}{additional arguments passed to
  \code{\link[mboost]{mboost}}, including \code{offset},
  \code{family} and \code{control}.}
}
\value{
  on object of class \code{FDboost} that inherits from
  \code{mboost}. Special \code{\link{predict.FDboost}},
  \code{\link{coef.FDboost}} and \code{\link{plot.FDboost}}
  methods are available. The methods of
  \code{\link[mboost]{mboost}} are available as well, e.g.
  \code{\link[mboost]{extract}}.

  The \code{FDboost}-object is a named list containing:
  \item{...}{all elements of an
  \code{\link[mboost]{mboost}-object}} \item{yname}{the
  name of the response} \item{yind}{the observation points
  of the response, with its name as attribute}
  \item{data}{the data that was used for the model fit}
  \item{id}{NULL for a response over a regular grid,
  otherwise the id variable of the response}
  \item{predictOffset}{the function to predict the smooth
  offset} \item{offsetVec}{the offset for one trajectory
  for regular response and otherwise the offset for all
  trajectories} \item{callEval}{the evaluated function
  call} \item{timeformula}{the time-formula}
  \item{formulaFDboost}{the formula with which
  \code{FDboost} was called} \item{formulaMboost}{the
  formula with which \code{mboost} was called within
  \code{FDboost}}
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
  \eqn{Y_i(t)} are stacked into one long vector. If the
  functional response is supplied as a vector in long
  format, the argument \code{id} has to be specified.
}
\examples{
data("viscosity", package = "FDboost")
## set time-interval that should be modeled
interval <- "101"

## model time until "interval" and take log() of viscosity
end <- which(viscosity$timeAll==as.numeric(interval))
viscosity$vis <- log(viscosity$visAll[,1:end])
viscosity$time <- viscosity$timeAll[1:end]
# with(viscosity, funplot(time, vis, pch=16, cex=0.2))

## fit median regression model with 200 boosting iterations,
## step-length 0.2 and smooth time-specific offset
mod <- FDboost(vis ~ 1 + bols(T_C) + bols(T_A),
               timeformula=~bbs(time, lambda=100),
               numInt="Riemann", family=QuantReg(),
               offset=NULL, offset_control = o_control(k_min = 9),
               data=viscosity, control=boost_control(mstop = 100, nu = 0.4))
summary(mod)
}
\author{
  Sarah Brockhaus, Torsten Hothorn
}
\seealso{
  \code{\link[mboost]{mboost}} for the help of the wrapped
  function in package mboost. See
  \code{\link[FDboost]{bsignal}} and
  \code{\link[FDboost]{bbsc}} for possible base-learners
}
\keyword{models,}
\keyword{nonlinear}

