\name{plot.FDboost}
\alias{plot.FDboost}
\alias{plotPredicted}
\title{Plot the fit or the coefficients of a boosted functional regression model}
\usage{
  \method{plot}{FDboost} (x, raw = FALSE, rug = TRUE,
    which = NULL, includeOffset = TRUE, ask = TRUE,
    n1 = 100, n2 = 40, n3 = 20, n4 = 11,
    onlySelected = TRUE, perspPlot = FALSE,
    commonRange = FALSE, ...)

  plotPredicted(x, subset = 1:x$ydim[1],
    posLegend = "topleft", ...)
}
\arguments{
  \item{x}{a fitted \code{FDboost}-object}

  \item{raw}{logical defaults to FALSE. If raw=FALSE for
  each effect the estimated function/surface is calcuated
  If raw=TRUE the coefficients of the model are returned.}

  \item{rug}{when TRUE (default) then the covariate to
  which the plot applies is displayed as a rug plot at the
  foot of each plot of a 1-d smooth, and the locations of
  the covariates are plotted as points on the contour plot
  representing a 2-d smooth.}

  \item{which}{a subset of base-learners to take into
  account for plotting. a list is returned.}

  \item{includeOffset}{logical, defaults to TRUE. Should
  the offset be included in the plot of the intercept
  (default) or should it be plotted separately.}

  \item{ask}{logical, defaults to TRUE, if several effects
  are plotted the user has to hit Return to see next plot.}

  \item{n1}{see below}

  \item{n2}{see below}

  \item{n3}{n1, n2, n3 give the number of gridpoints for
  1-/2-/3-dimensional smooth terms used in the marginal
  equidistant grids over the range of the covariates at
  which the estimated effects are evaluated.}

  \item{n4}{gives the number of points for the third
  dimension in a 3-dimensional smooth term}

  \item{onlySelected,}{logical, defaults to TRUE. Only plot
  effects that were selected in at least one boosting
  iteration}

  \item{perspPlot}{logical, defaults to FALSE, If TRUE,
  perspective plots (\code{\link[graphics]{persp}}) for 2-
  and 3-dimensional effects are drawn. If FALSE,
  image/contour-plots (\code{\link[graphics]{image}},
  \code{\link[graphics]{contour}}) are drawn for 2- and
  3-dimensional effects.}

  \item{commonRange}{logical, defaults to FALSE, if TRUE
  the range over all effects is the same (so far not
  implemented)}

  \item{subset}{subset of the observed response curves and
  their predictions that is plotted. Per default all
  observations are plotted.}

  \item{posLegend}{location of the legend, if a legend is
  drawn automatically (only in plotPredicted). The default
  is "topleft".}

  \item{...}{other arguments, not used.}
}
\description{
  Takes a fitted \code{FDboost}-object produced by
  \code{\link{FDboost}()} and plots the fitted effects or
  the coefficient-functions/surfaces.
}
\seealso{
  \code{\link{FDboost}} for the model fit and
  \code{\link{coef.FDboost}} for the calculation of the
  coefficient functions.
}

