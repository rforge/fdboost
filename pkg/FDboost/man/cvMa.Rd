\name{cvMa}
\alias{cvMa}
\title{Function to set up folds for a response-matrix}
\usage{
  cvMa(ydim = NULL, Y = NULL, weights = NULL,
    type = c("bootstrap", "kfold", "subsampling"),
    B = ifelse(type == "kfold", 10, 25), prob = 0.5,
    strata = NULL, id = NULL)
}
\arguments{
  \item{ydim}{dimensions of reponse-matrix}

  \item{Y}{response-matrix}

  \item{weights}{a numeric vector of weights for the model
  to be cross-validated.}

  \item{type}{character argument for specifying the
  cross-validation method. Currently (stratified)
  bootstrap, k-fold cross-validation and subsampling are
  implemented.}

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
}
\value{
  \code{cvMa} retruns a matrix of weights to be used in
  \code{cvrisk}.
}
\description{
  Function to set up folds for a response-matrix that are
  used in cross-validation
}
\details{
  It is sufficient to specify one of \code{ydim} and
  \code{Y}.

  The function \code{cvMa} can be used to build an
  appropriate weight matrix to be used with \code{cvrisk}.
  In \code{cvMa()} whole trajectories are sampled.  If
  \code{strata} is defined sampling is performed in each
  stratum separately thus preserving the distribution of
  the \code{strata} variable in each fold.  If \code{id} is
  defined sampling is performed on the level of \code{id}
  thus sampling whole individuals.
}
\examples{
Ytest <- matrix(rnorm(15), ncol=3) # 5 trajectories, each with 3 observations
cvMa(Y=Ytest, type="bootstrap", B=4) # 4 columns with bootstrap-weights
cvMa(ydim=c(5,3), type="bootstrap", B=4) # the same
}
\seealso{
  \code{\link{cvrisk}} to perform cross-validation.
}
\keyword{models,}
\keyword{regression}

