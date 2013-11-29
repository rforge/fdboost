\name{validateFDboost}
\alias{cvMa}
\alias{validateFDboost}
\title{Cross-Validation over Curves}
\usage{
  validateFDboost(object, response = NULL,
    weights = model.weights(object),
    folds = cvMa(ydim = object$ydim, weights = model.weights(object), type = "bootstrap"),
    grid = 1:mstop(object), getCoefCV = TRUE,
    mrdDelete = 0, ...)
}
\arguments{
  \item{object}{fitted FDboost-object}

  \item{response}{you can specify a response vector to
  calculate predictions errors. Defaults fo NULL which
  means that the response of the fitted model is used.}

  \item{folds}{a weight matrix with number of rows equal to
  the number of observations.}

  \item{grid}{the grid over which the optimal mstop is
  searched}

  \item{getCoefCV}{logical, defaults to TRUE Should the
  coefficients and predictions be computed for all the
  models on the sampled data?}

  \item{mrdDelete}{Delete values that are mrdDelete percent
  smaller then the mean of the response. Defaults to 0
  which means that only response values beeing 0 are not
  used in the calculaiton of the MRD (= mean relative
  deviation)}

  \item{ydim}{dimensions of reponse-matrix}

  \item{weights}{a numeric vector of weights for the model
  to be cross-validated. if weights=NULL all weights are
  taken to be 1.}

  \item{type}{character argument for specifying the
  cross-validation method. Currently (stratified)
  bootstrap, k-fold cross-validation and subsampling are
  implemented. The argument curves implies that a
  cross-validation leaving out one curve at a time is
  performed.}

  \item{B}{number of folds, per default 25 for
  \code{bootstrap} and \code{subsampling} and 10 for
  \code{kfold}.}

  \item{prob}{percentage of observations to be included in
  the learning samples for subsampling.}

  \item{strata}{a factor of the same length as
  \code{weights} for stratification.}

  \item{id}{id-variable to sample upon ids instead of
  sampling single observations. (Only interesting in the
  case that several observations per individual were
  made.)}

  \item{...}{further arguments passed to mclapply if
  parallel=TRUE, otherwise ignored}
}
\value{
  \code{cvMa} retruns a matrix of weights to be used in
  \code{cvrisk} or \code{validateFDboost}.
}
\description{
  Cross-Validation over Curves.
}
\details{
  The funciton \code{validateFDboost} calculates honest
  estimates of prediction errors as the curve/observations
  that should be predicted is not part of the model fit.
  The function \code{cvMa} can be used to build an
  appropriate weight matrix to be used with \code{cvrisk}
  or \code{validateFDboost}.  In \code{cvMa} whole
  trajectories are sampled. The probability for each
  trajectory to enter a fold is equal over all
  trajectories.  If \code{strata} is defined sampling is
  performed in each stratum separately thus preserving the
  distribution of the \code{strata} variable in each fold.
  If \code{id} is defined sampling is performed on the
  level of \code{id} thus sampling individuals.
}
\note{
  Use argument \code{mc.cores = 1L} to set the numbers of
  cores that is used in parallel computation.
}
\examples{
Ytest <- matrix(rnorm(15), ncol=3) # 5 trajectories, each with 3 observations
cvMa(ydim=c(5,3), type="bootstrap", B=4)
}
\seealso{
  \code{\link{cvrisk}} to perform cross-validation.
}

