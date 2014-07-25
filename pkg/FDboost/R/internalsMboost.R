
####################################################
#### internal utility functions of mboost 2.3-0  ###
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


####### from helpers.R
### check for classical or Matrix matrices
isMATRIX <- function(x)
  is.matrix(x) || is(x, "Matrix")
