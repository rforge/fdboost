\name{bsignal}
\alias{bconcurrent}
\alias{bsignal}
\title{Base-learners for Functional Covarites}
\usage{
  bsignal(..., knots = 10, boundary.knots = NULL,
    degree = 3, differences = 2, df = 4, cyclic = FALSE)

  bconcurrent(..., knots = 10, boundary.knots = NULL,
    degree = 3, differences = 2, df = 4)
}
\arguments{
  \item{...}{matrix of functional data and the vector of
  observation points. The functional covariates have to be
  supplied as n by <no. of evaluations> matrices, i.e. each
  row is one functional observation. If no vector of
  observation points is supplied an equidistant grid over
  \eqn{ [0,1]} is assumed.}

  \item{knots}{either the number of knots or a vector of
  the positions of the interior knots (for more details see
  \code{\link[mboost]{bbs})}.}

  \item{boundary.knots}{boundary points at which to anchor
  the B-spline basis (default the range of the data). A
  vector (of length 2) for the lower and the upper boundary
  knot can be specified.}

  \item{degree}{degree of the regression spline.}

  \item{differences}{a non-negative integer, typically 1, 2
  or 3. If \code{differences} = \emph{k}, \emph{k}-th-order
  differences are used as a penalty (\emph{0}-th order
  differences specify a ridge penalty).}

  \item{df}{trace of the hat matrix for the base-learner
  defining the base-learner complexity. Low values of
  \code{df} correspond to a large amount of smoothing and
  thus to "weaker" base-learners.}

  \item{cyclic}{not yet implemented}
}
\value{
  Equally to the base-learners of the package mboost:

  An object of class \code{blg} (base-learner generator)
  with a \code{dpp} function.

  The call of \code{dpp} returns an object of class
  \code{bl} (base-learner) with a \code{fit} function. The
  call to \code{fit} finally returns an object of class
  \code{bm} (base-model).
}
\description{
  Base-learners that fit effects of functional covariates
}
\details{
  \code{bsignal} implements a base-learner for functional
  covariates to estimate an effect of the form
  \eqn{X_i(s)\beta(t,s)ds}. Defaults to a cubic B-spline
  basis with second difference penalties for
  \eqn{\beta(t,s)} in the direction of s and numerical
  integration over the entire range by using trapezoidal
  Riemann weights.

  \code{bconcurrent} implements a concurrent effect for a
  functional covariate on a functional response, i.e. an
  effect of the form \eqn{X_i(t)\beta(t)} for a functional
  response \eqn{Y_i(t)}.

  It is recommended to use centered functional covariates
  with \eqn{\sum_i X_i(s) = 0} for all \eqn{s} in
  \code{bconcurrent}- and \code{bsignal}-terms so that the
  global functional intercept can be interpreted as the
  global mean function.

  Cannot deal with any missing values in the covariates.
}
\examples{
###############################################################################
print("to do")
###############################################################################
}
\references{
  Ivanescu, A., Staicu, A.-M., Scheipl, F. and Greven, S.
  (2012).  Penalized function-on-function regression.
  (submitted)
  \url{http://biostats.bepress.com/jhubiostat/paper240/}
}
\seealso{
  \code{\link{FDboost}} for the model fit.
  \code{\link[mboost]{bbs}} and
  \code{\link[mboost]{brandom}} for the corresponding
  base-learners in mboost.
}
\keyword{models}

