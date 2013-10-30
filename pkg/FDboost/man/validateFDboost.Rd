\name{validateFDboost}
\alias{validateFDboost}
\title{Validate model by resampling over curves}
\usage{
  validateFDboost(object, response = NULL,
    type = c("curves", "kfold"), B = 5, folds = NULL,
    grid = 1:mstop(object), getCoefCV = FALSE,
    mrdDelete = 0, ...)
}
\arguments{
  \item{object}{fitted FDboost-object}

  \item{response}{you can specify a response vector to
  calculate predictions errors. Defaults fo NULL which
  means that the response of the fitted model is used.}

  \item{type}{character argument for specifying the
  cross-validation method. Currently cross-validation over
  curves (\code{curves}) and random cross-validation
  (\code{kfold}) are impelemted.}

  \item{B}{number of folds for random cross-validation.}

  \item{folds}{a weight matrix with number of rows equal to
  the number of observations. The number of columns
  corresponds to the number of cross-validation runs. Can
  be computed using function \code{cvMa} or \code{cv}.}

  \item{grid}{the grid over which the optimal mstop is
  searched}

  \item{getCoefCV}{logical, defaults to FALSE. Should the
  coefficient-functions/surfaces be computed for all the
  models on the sampled data?}

  \item{mrdDelete}{Delete values that are mrdDelete percent
  smaller then the mean of the response. Defaults to 0
  which means that only response values beeing 0 are not
  used in the calculaiton of the MRD (= mean relative
  deviation)}

  \item{...}{further arguments passed to mclapply if
  parallel=TRUE, otherwise ignored}
}
\description{
  Validate the model fit by refitting the model $n$ times.
  Each time leaving out one observation.
}
\details{
  Calculates honest estimates of prediction errors as the
  curve that should be predicted is not part of the model
  fit.
}
\note{
  Use argument \code{mc.cores = 1L} to set the numbers of
  cores that is used in parallel computation.
}

