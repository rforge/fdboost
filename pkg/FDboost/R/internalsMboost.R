
####################################################
#### internal utility functions of mboost 2.4-2  ###
####################################################

#######from helpers.R
### rows without missings in Matrices, matrices and data.frames
Complete.cases <- function(x) {
  if (isMATRIX(x)) return(rowSums(is.na(x)) == 0)
  complete.cases(x)
}


####### from helpers.R
### try to find duplicated entries / observations with missing values
### for more efficient memory handling
get_index <- function(x) {
  
  if (isMATRIX(x)) {
    ### handle missing values only
    cc <- Complete.cases(x)
    nd <- which(cc)
    index <- match(1:nrow(x), nd)
  } else {
    ### handle single variables (factors / numerics) factors
    if (length(x) == 1) {
      x <- x[[1]]
      nd <- which(!duplicated(x))
      nd <- nd[complete.cases(x[nd])]
      index <- match(x, x[nd])
      ### go for data.frames with >= 2 variables
    } else {
      tmp <- do.call("paste", x)
      nd <- which(!duplicated(tmp))
      nd <- nd[complete.cases(x[nd,])]
      index <- match(tmp, tmp[nd])
    }
  }
  return(list(nd, index))
}

####### from helpers.R
### check for classical or Matrix matrices
isMATRIX <- function(x)
  is.matrix(x) || is(x, "Matrix")



####### from bl.R
### hyper parameters for P-splines baselearner (including tensor product P-splines)
hyper_bbs <- function(mf, vary, knots = 20, boundary.knots = NULL, degree = 3,
                      differences = 2, df = 4, lambda = NULL, center = FALSE,
                      cyclic = FALSE, constraint = "none", deriv = 0L) {
  
  knotf <- function(x, knots, boundary.knots) {
    if (is.null(boundary.knots))
      boundary.knots <- range(x, na.rm = TRUE)
    ## <fixme> At the moment only NULL or 2 boundary knots can be specified.
    ## Knot expansion is done automatically on an equidistand grid.</fixme>
    if ((length(boundary.knots) != 2) || !boundary.knots[1] < boundary.knots[2])
      stop("boundary.knots must be a vector (or a list of vectors) ",
           "of length 2 in increasing order")
    if (length(knots) == 1) {
      knots <- seq(from = boundary.knots[1],
                   to = boundary.knots[2], length = knots + 2)
      knots <- knots[2:(length(knots) - 1)]
    }
    list(knots = knots, boundary.knots = boundary.knots)
  }
  nm <- colnames(mf)[colnames(mf) != vary]
  if (is.list(knots)) if(!all(names(knots) %in% nm))
    stop("variable names and knot names must be the same")
  if (is.list(boundary.knots)) if(!all(names(boundary.knots) %in% nm))
    stop("variable names and boundary.knot names must be the same")
  if (!identical(center, FALSE) && cyclic)
    stop("centering of cyclic covariates not yet implemented")
  ret <- vector(mode = "list", length = length(nm))
  names(ret) <- nm
  for (n in nm)
    ret[[n]] <- knotf(mf[[n]], if (is.list(knots)) knots[[n]] else knots,
                      if (is.list(boundary.knots)) boundary.knots[[n]]
                      else boundary.knots)
  if (cyclic & constraint != "none")
    stop("constraints not implemented for cyclic B-splines")
  stopifnot(is.numeric(deriv) & length(deriv) == 1)
  list(knots = ret, degree = degree, differences = differences,
       df = df, lambda = lambda, center = center, cyclic = cyclic,
       Ts_constraint = constraint, deriv = deriv)
}

### cyclic B-splines
### adapted version of mgcv:cSplineDes from S.N. Wood
cbs <- function (x, knots, boundary.knots, degree = 3, deriv = 0L) {
  nx <- names(x)
  x <- as.vector(x)
  ## handling of NAs
  nax <- is.na(x)
  if (nas <- any(nax))
    x <- x[!nax]
  
  knots <- c(boundary.knots[1], knots, boundary.knots[2])
  nKnots <- length(knots)
  ord <- degree + 1
  xc <- knots[nKnots - ord + 1]
  knots <- c(boundary.knots[1] -
               (boundary.knots[2] - knots[(nKnots - ord + 1):(nKnots - 1)]),
             knots)
  ind <- x > xc
  X <- splineDesign(knots, x, ord, derivs = rep(deriv, length(x)), outer.ok = TRUE)
  x[ind] <- x[ind] - boundary.knots[2] + boundary.knots[1]
  if (sum(ind)) {
    Xtmp <- splineDesign(knots, x[ind], ord, derivs = rep(deriv, length(x[ind])),
                         outer.ok = TRUE)
    X[ind, ] <- X[ind, ] + Xtmp
  }
  ## handling of NAs
  if (nas) {
    tmp <- matrix(NA, length(nax), ncol(X))
    tmp[!nax, ] <- X
    X <- tmp
  }
  ## add attributes
  attr(X, "degree") <- degree
  attr(X,"knots") <- knots
  attr(X,"boundary.knots") <- boundary.knots
  if (deriv != 0)
    attr(X, "deriv") <- deriv
  dimnames(X) <- list(nx, 1L:ncol(X))
  return(X)
}

bsplines <- function(x, knots, boundary.knots, degree,
                     Ts_constraint = "none", deriv = 0L){
  nx <- names(x)
  x <- as.vector(x)
  ## handling of NAs
  nax <- is.na(x)
  if (nas <- any(nax))
    x <- x[!nax]
  ## use equidistant boundary knots
  dx <- diff(boundary.knots)/(length(knots) + 1)
  bk_lower <- seq(boundary.knots[1] - degree * dx, boundary.knots[1],
                  length = degree + 1)
  bk_upper <- seq(boundary.knots[2], boundary.knots[2] + degree * dx,
                  length = degree + 1)
  ## complete knot mesh
  k <- c(bk_lower, knots, bk_upper)
  ## construct design matrix
  X <- splineDesign(k, x, degree + 1, derivs = rep(deriv, length(x)),
                    outer.ok = TRUE)
  ## handling of NAs
  if (nas) {
    tmp <- matrix(NA, length(nax), ncol(X))
    tmp[!nax, ] <- X
    X <- tmp
  }
  ### constraints; experimental
  D <- diag(ncol(X))
  D[lower.tri(D)] <- 1
  X <- switch(Ts_constraint, "none" = X,
              "increasing" = X %*% D,
              "decreasing" = -X %*% D)
  ## add attributes
  attr(X, "degree") <- degree
  attr(X, "knots") <- knots
  attr(X, "boundary.knots") <- list(lower = bk_lower, upper = bk_upper)
  if (Ts_constraint != "none")
    attr(X, "Ts_constraint") <- Ts_constraint
  if (Ts_constraint != "none")
    attr(X, "D") <- D
  if (deriv != 0)
    attr(X, "deriv") <- deriv
  dimnames(X) <- list(nx, 1L:ncol(X))
  return(X)
}


##################################################################################

# further utility functions of library mboost, bl.R
# necessary to copy them into FDboost?

### extract variables names from base-learner
names.blg <- function(x)
  x$get_names()

### extract data from base-learner
model.frame.blg <- function(formula, ...)
  formula$model.frame(...)

# ### extract coefficients
# coef.bm_lin <- function(object, ...) {
#   ret <- as.vector(object$model)
#   names(ret) <- object$Xnames
#   ret
# }
# 
# ### extract fitted values
# fitted.bm <- function(object)
#   object$fitted()
# 
# ### extract hatmatrix
# hatvalues.bl_lin <- function(model)
#   model$hatvalues()
# 
# ### data preprocessing (plug in weights)
# dpp <- function(object, weights)
#   UseMethod("dpp", object)
# 
# dpp.blg <- function(object, weights)
#   object$dpp(weights)
# 
# ### actually fit a base-learner to response y
# fit <- function(object, y)
#   UseMethod("fit", object)
# 
# fit.bl <- function(object, y)
#   object$fit(y)

