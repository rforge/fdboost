\name{mstop.validateFDboost}
\alias{mstop.validateFDboost}
\alias{plot.validateFDboost}
\title{Methods for objects of calss validateFDboost}
\usage{
  \method{mstop}{validateFDboost} (object, ...)

  \method{plot}{validateFDboost} (x,
    risk = c("median", "mean"), modObject = NULL,
    predictNA = FALSE, names.arg = NULL, ask = TRUE, ...)
}
\arguments{
  \item{x}{object of class validateFDboost}

  \item{object}{object of class validateFDboost}

  \item{risk}{which risk is minimized to obtain the optimal
  stopping iteration? defaults to the median}

  \item{modObject}{if the original model object of class
  \code{FDboost} is given predicted values of the whole
  model can be compared to the predictions of the
  cross-validated models}

  \item{predictNA}{should missing values in the response be
  predicted? Defaults to FALSE.}

  \item{names.arg}{names of the observed curves}

  \item{ask}{par(ask=ask)}

  \item{...}{additional arguments passed to callies.}
}
\description{
  Methods for objects that are fitted to determine the
  optimal mstop and the prediciton error of a model fitted
  by FDboost.
}
\details{
  \code{plot.validateFDboost} plots cross-validated risk,
  RMSE, MRD, measured and predicted values and residuals as
  determined by validateFDboost.
  \code{mstop.validateFDboost} extracts the optimal mstop
  by minimizing the median risk.
}

