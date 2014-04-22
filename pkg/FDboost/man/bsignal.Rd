\name{bsignal}
\alias{bconcurrent}
\alias{bhist}
\alias{bsignal}
\title{Base-learners for Functional Covarites}
\usage{
  bsignal(..., index = NULL, knots = 10,
    boundary.knots = NULL, degree = 3, differences = 2,
    df = 4, lambda = NULL, cyclic = FALSE, Z = NULL)

  bconcurrent(..., knots = 10, boundary.knots = NULL,
    degree = 3, differences = 2, df = 4)

  bhist(..., index = NULL, knots = 10,
    boundary.knots = NULL, degree = 3, differences = 2,
    df = 4, limits = "s<=t")
}
\arguments{
  \item{...}{matrix of functional data and the vector of
  observation points. The functional covariates have to be
  supplied as n by <no. of evaluations> matrices, i.e. each
  row is one functional observation. The base-learner
  \code{bhist} expects three arguments: functional
  covariate, index of functional covariate, index of
  functional response \eqn{[0,1]} is assumed.}

  \item{index}{a vector of integers for expanding the
  signal variable in \code{....} For example, bsignal(X, s,
  index = index) is equal to bsignal(X[index,], s), where
  index is an integer of length greater or equal to
  length(x).}

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

  \item{lambda}{smoothing penalty}

  \item{cyclic}{if cyclic = TRUE the fitted coefficient
  function coincides at the boundaries (useful for cyclic
  covariates such as day time etc.).}

  \item{Z}{a transformation matrix for the design-matrix
  over the index of the covariate. Z can be calculated as
  the transformation matrix for a sum-to-zero constraint in
  the case that all trajectories have the same mean (then a
  shift in the coefficient function is not identifiable).}

  \item{limits}{defaults to "s<=t" for an historical effect
  with s<=t, otherwise specifies the integration limits
  s_{hi, i}, s_{lo, i}: either one of "s<t" or "s<=t" for
  (s_{hi, i}, s_{lo, i}) = (0, t) or a function that takes
  s as the first and t as the second argument and returns
  TRUE for combinations of values (s,t) if s falls into the
  integration range for the given t. This is an
  experimental feature and not well tested yet; use at your
  own risk.}
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
### example for scalar response and two functional covariates
data(fuel)
modFuel <- FDboost(heatan ~ bsignal(UVVIS, uvvis.lambda, knots=40, df=4)
           + bsignal(NIR, nir.lambda, knots=40, df=4),
           timeformula=~bols(1), data=fuel)
summary(modFuel)
## plot(modFuel, rug=FALSE)
}
\references{
  Scheipl, F., Staicu, A.-M., and Greven, S. (2014),
  Functional Additive Mixed Models, Journal of
  Computational and Graphical Statistics, in press, DOI
  10.1080/10618600.2014.901914.
  \url{http://arxiv.org/abs/1207.5947}
}
\seealso{
  \code{\link{FDboost}} for the model fit.
}
\keyword{models}

