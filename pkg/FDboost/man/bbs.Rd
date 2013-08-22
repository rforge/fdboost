\name{bbs}
\alias{bbs}
\title{Base-learner for Scalar Covariate}
\usage{
  bbs(..., by = NULL, index = NULL, knots = 10,
    boundary.knots = NULL, outer.knots = TRUE, degree = 3,
    differences = 2, df = 4, lambda = NULL, center = FALSE,
    cyclic = FALSE)
}
\arguments{
  \item{...}{one or more predictor variables or one matrix
  or data frame of predictor variables.}

  \item{by}{an optional variable defining varying
  coefficients, either a factor or numeric variable.}

  \item{index}{a vector of integers for expanding the
  variables in \code{...}.}

  \item{knots}{either the number of knots or a vector of
  the positions of the interior knots (for more details see
  \code{\link[mboost]{bbs}}).}

  \item{boundary.knots}{boundary points at which to anchor
  the B-spline basis (default the range of the data). A
  vector (of length 2) for the lower and the upper boundary
  knot can be specified.}

  \item{outer.knots}{logical, the outer knots are outside
  of the range (default) or just repetitions the minimum
  and maximum value}

  \item{degree}{degree of the regression spline.}

  \item{differences}{a non-negative integer, typically 1, 2
  or 3. If \code{differences} = \emph{k}, \emph{k}-th-order
  differences are used as a penalty (\emph{0}-th order
  differences specify a ridge penalty).}

  \item{df}{trace of the hat matrix for the base-learner
  defining the base-learner complexity. Low values of
  \code{df} correspond to a large amount of smoothing and
  thus to "weaker" base-learners.}

  \item{lambda}{smoothing penalty, computed from \code{df}
  when \code{df} is specified.}

  \item{center}{not implemented yet}

  \item{cyclic}{if \code{cyclic = TRUE} the fitted values
  coincide at the boundaries (useful for cyclic covariates
  such as day time etc.).}
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
  The base-learner \code{bbs} of the package FDboost is per
  default the same as the base-learner \code{bbs} of
  package mboost but it has the extra argument
  \code{outer.knots}.
}
\seealso{
  \code{\link{FDboost}} for the model fit.
  \code{\link[mboost]{bbs}} for the corresponding
  base-learners in \code{mboost}.
}
\keyword{models}

