#' Functions to compute integration weights
#' 
#' Computes trapezoidal integration weights for a functional variable X1 on grid xind.
#' @param X1 matrix of functional variable
#' @param xind index of functional variable
#' @param id defaults to NULL if X1 is a matrix. identity variable if X1 is in long format.
#' 
#' @aliases integrationWeightsLeft
#' 
#'  @details The function \code{integrationWeights()} computes trapezoidal integration weights, 
#'  that are symmetric. The function \code{integrationWeightsLeft()} computes weights,
#'  that take into account only the distance to the prior observation point. 
#'  Those weights are adequate for historical effects.
#'   
#' 
#' @export
################################# 
# Trapezoidal integration weights for a functional variable X1 on grid xind
# corresponds to mean of left and right Riemann integration sum
integrationWeights <- function(X1, xind, id=NULL){
  
  if(is.null(id)) if(ncol(X1)!=length(xind) ) stop("Dimension of xind and X1 do not match")
  
  # compute integraion weights for irregular data in long format
  if(!is.null(id)){
    Lneu <- tapply(xind, id, FUN = function(x) colMeans(rbind(c(0,diff(x)), c(diff(x), 0))) )
    Lneu <- unlist(Lneu)
    names(Lneu) <- NULL    
    return(Lneu) 
  }  
  
  #Li <- c(diff(xind)[1]/2, diff(xind)[2:(length(xind)-1)], diff(xind)[length(xind)-1]/2)
  
  # Riemann integration weights: mean weights between left and right sum
  # \int^b_a f(t) dt = sum_i (t_i-t_{i-1})*(f(t_i)) (left sum)
  # Li <- c(0,diff(xind))
  
  # trapezoidal
  # \int^b_a f(t) dt = .5* sum_i (t_i - t_{i-1}) f(t_i) + f(t_{i-1}) = 
  #  (t_2 - t_1)/2 * f(a=t_1) + sum^{nx-1}_{i=2} ((t_i - t_{i-1})/2 + (t_{i+1} - t_i)/2) * f(t_i) 
  #			+ (t_{nx} - t_{nx-1})/2 * f(b=t_{nx})
  Li <- colMeans(rbind(c(0,diff(xind)), c(diff(xind), 0))) 
  
  # alternative calculation of trapezoidal weights
  #diffs <- diff(xind)
  #nxgrid <- length(xind)
  #Li <- 0.5 * c(diffs[1],  filter(diffs, filter=c(1,1))[-(nxgrid-1)], diffs[(nxgrid-1)] )
    
  L <- matrix(Li, nrow=nrow(X1), ncol=ncol(X1), byrow=TRUE)
  
  # taking into account missing values
  if(any(is.na(X1))){
    Lneu <- sapply(1:nrow(X1), function(i){
      x <- X1[i,]
      
      if(!any(is.na(x))){
        l <- L[i, ] # no missing values in curve i
      }else{
        xindL <- xind # lower
        xindL[is.na(x)] <- NA
        xindU <- xindL # upper
        
        if(is.na(xindL[1])){ # first observation is missing
          xindL[1] <- xind[1] - diff(c(xind[1], xind[2])) 
        } 
        if(is.na(xindU[length(xind)])){ # last observation is missing
          xindU[length(xind)] <- xind[length(xind)] + diff(c(xind[length(xind)-1], xind[length(xind)])) 
        }
        
        xindL <- na.locf(xindL, na.rm=FALSE) # index for lower sum
        xindU <- na.locf(xindU, fromLast=TRUE, na.rm=FALSE) # index for upper sum
        
        l <- colMeans(rbind(c(0,diff(xindL)), c(diff(xindU), 0))) # weight is 0 for missing values
      }
      return(l)
    }
           )
    
    return(t(Lneu))  
    
  }else{ 
    return(L)
  }
}

# # test integrationWeights()
# xind <- seq(0,1,l=5)
# xind <- c(0, 0.2, 0.4, 0.5, 1)
# colMeans(rbind(c(0,diff(xind)), c(diff(xind), 0)))
# X1 <- matrix(xind^2, ncol=5)
# intW <- integrationWeights(X1, xind)[1,]
# plot(X1[1,]~xind, type="b")
# points(rep(0,5)~cumsum(intW), col=2)
# 
# X1 <- matrix(c(1:5, 1:5+1, 1:5+2, 1:5+3), ncol=5, byrow=TRUE)
# X1[1,1] <- NA
# X1[1,2] <- NA
# X1[2,2] <- NA
# X1[3,5] <- NA
# xind <- c(2,4,6,8,10)
# 
# intW <- integrationWeights(X1, xind)
# rowSums(intW*X1, na.rm=TRUE) -c(0, 5, 10, 15)
# matplot(xind, t(X1), type="b")



#### Computes Riemann-weights that only take into account the distance to the previous 
# observation point
# important for bhist()
#' @export
integrationWeightsLeft <- function(X1, xind){
  
  if( ncol(X1)!=length(xind) ) stop("Dimension of xind and X1 do not match")
  
  # use lower Riemann sum
  Li <- diff(xind)
  Li <- c(Li[1], Li)
  L <- matrix(Li, nrow=nrow(X1), ncol=ncol(X1), byrow=TRUE)
  
  return(L)
}


################################################################################
### syntax for base learners is modified code of the package mboost, see bl.R



################################################################################
################################################################################
# Base-learners for functional covariates 

### model.matrix for P-splines base-learner of signal matrix mf
### with index/time as attribute
X_bsignal <- function(mf, vary, args) {
  
  stopifnot(is.data.frame(mf))
  xname <- attr(mf[[1]], "xname")
  X1 <- as.matrix(mf[[1]])
  xind <- attr(mf[[1]], "signalIndex")
  
  #   stopifnot(is.list(mf))
  #   xname <- names(mf)[1]
  #   X1 <- mf[[1]]
  #   xind <- mf[[2]]
  
  if(ncol(X1)!=length(xind)) stop("Dimension of signal matrix and its index do not match.")
    
  ### Construct spline basis over index xind of X1 
  if(is.null(args$boundary.knots))  args$boundary.knots <- range(xind, na.rm = TRUE)
  knots <- seq(from = args$boundary.knots[1], to = args$boundary.knots[2], length = args$knots + 2)
  knots <- knots[2:(length(knots) - 1)]
  
  # B-spline basis of specified degree  
  #X <- bs(xind, knots=knots, degree=args$degree, intercept=TRUE) # old version   
  X <- mboost:::bsplines(xind, knots=knots, boundary.knots=args$boundary.knots, 
                         degree=args$degree)
  
  # use cyclic splines
  if (args$cyclic) {
    X <- mboost:::cbs(xind,
             knots = knots,
             boundary.knots = args$boundary.knots,
             degree = args$degree)
  }
    
  colnames(X) <- paste(xname, 1:ncol(X), sep="")
      
  ### use the transformation matrix Z if necessary
  useZ <- FALSE
  ### Check whether integral over trajectories is different then centering is advisable
  if(is.null(args$Z) && all(rowMeans(X1)-mean(rowMeans(X1)) < .Machine$double.eps *10^10)){
    message(paste("All trajectories in ", xname, " have the same mean. Coefficient function is centered.", sep=""))
    useZ <- TRUE
  }
  
  ### in case that a Z-matrix is given in bsignal()
  if(!is.null(args$Z)) useZ <- TRUE  
  #browser()
  
  if(useZ  | !is.null(args$Z)){
    
    #----------------------------------
    ### <SB> Calculate constraints
    
    # If the argument Z is not NULL use the given Z (important for prediction!)
    if(is.null(args$Z)){
      C <- t(X) %*% rep(1, nrow(X))
      Q <- qr.Q(qr(C), complete=TRUE) # orthonormal matrix of QR decomposition
      Z <- Q[  , 2:ncol(Q)] # only keep last columns    
    }else Z <- args$Z
    
    ### Transform design and penalty matrix
    X <- X %*% Z
    #K <- t(Z) %*% K %*% Z
    
    attr(X, "Z") <- Z # remember Z
    ### myZ <<- Z
    #---------------------------------- 
  }else{
    Z <- NULL
  }
  
  #print("X_bsignal")
  #print(Z[1:3,1:3])
  
  ### Weighting with matrix of functional covariate
  L <- integrationWeights(X1=X1, xind=xind)
  # Design matrix is product of weighted X1 and basis expansion over xind 
  X <- (L*X1) %*% X
  attr(X, "Z") <- Z # remember Z
  
  ### Penalty matrix: product differences matrix
  if (args$differences > 0){
    if (!args$cyclic) {
      K <- diff(diag(ncol(X)), differences = args$differences)
    } else {
      ## cyclic P-splines
      differences <- args$differences
      K <- diff(diag(ncol(X) + differences),
                differences = differences)
      tmp <- K[,(1:differences)]   # save first "differences" columns
      K <- K[,-(1:differences)]    # drop first "differences" columns
      indx <- (ncol(X) - differences + 1):(ncol(X))
      K[,indx] <- K[,indx] + tmp   # add first "differences" columns
    }
  } else {
    if (args$differences != 0)
      stop(sQuote("differences"), " must be an non-negative integer")
    K <- diag(ncol(X))
  }
  
  ### penalty matrix is squared difference matrix
  K <- crossprod(K)
 
  ## compare specified degrees of freedom to dimension of null space
  if (!is.null(args$df)){
    rns <- ncol(K) - qr(as.matrix(K))$rank # compute rank of null space
    if (rns == args$df)
      warning( sQuote("df"), " equal to rank of null space ",
               "(unpenalized part of P-spline);\n  ",
               "Consider larger value for ", sQuote("df"),
               " or set ", sQuote("center = TRUE"), ".", immediate.=TRUE)
    if (rns > args$df)
      stop("not possible to specify ", sQuote("df"),
           " smaller than the rank of the null space\n  ",
           "(unpenalized part of P-spline). Use larger value for ",
           sQuote("df"), " or set ", sQuote("center = TRUE"), ".")
  }
  return(list(X = X, K = K, Z = Z))
}

###############################################################################

#' Base-learners for Functional Covarites
#' 
#' Base-learners that fit effects of functional covariates 
#' 
#' @param ... matrix of functional data and the vector of observation points.
#' The functional covariates have to be supplied as n by <no. of evaluations> 
#' matrices, i.e. each row is one functional observation.
#' The base-learner \code{bhist} expects three arguments: functional covariate,
#' index of functional covariate, index of functional response 
#' \eqn{[0,1]} is assumed.
#' @param index a vector of integers for expanding the signal variable in \code{....} 
#' For example, bsignal(X, s, index = index) is equal to bsignal(X[index,], s), 
#' where index is an integer of length greater or equal to length(x).
#' @param knots either the number of knots or a vector of the positions 
#' of the interior knots (for more details see \code{\link[mboost]{bbs})}.
#' @param boundary.knots boundary points at which to anchor the B-spline basis 
#' (default the range of the data). A vector (of length 2) 
#' for the lower and the upper boundary knot can be specified.
#' @param degree degree of the regression spline.
#' @param differences a non-negative integer, typically 1, 2 or 3. 
#' If \code{differences} = \emph{k}, \emph{k}-th-order differences are used as 
#' a penalty (\emph{0}-th order differences specify a ridge penalty).
#' @param df trace of the hat matrix for the base-learner defining the 
#' base-learner complexity. Low values of \code{df} correspond to a 
#' large amount of smoothing and thus to "weaker" base-learners.
#' @param lambda smoothing penalty
#' @param cyclic if cyclic = TRUE the fitted coefficient function coincides at the boundaries 
#' (useful for cyclic covariates such as day time etc.).
#' @param Z a transformation matrix for the design-matrix over the index of the covariate.
#' Z can be calculated as the transformation matrix for a sum-to-zero constraint in the case
#' that all trajectories have the same mean 
#' (then a shift in the coefficient function is not identifiable).
#' @param limits defaults to "s<=t" for an historical effect with s<=t, 
#' otherwise specifies the integration limits s_{hi, i}, s_{lo, i}: 
#' either one of "s<t" or "s<=t" for (s_{hi, i}, s_{lo, i}) = (0, t) or a 
#' function that takes s as the first and t as the second argument and returns 
#' TRUE for combinations of values (s,t) if s falls into the integration range for the given t. 
#' This is an experimental feature and not well tested yet; use at your own risk. 
#' 
#' @aliases bconcurrent 
#' 
#' @details \code{bsignal} implements a base-learner for functional covariates to  
#' estimate an effect of the form \eqn{X_i(s)\beta(t,s)ds}. Defaults to a cubic  
#' B-spline basis with second difference penalties for \eqn{\beta(t,s)} in the direction 
#' of s and numerical integration over the entire range by using trapezoidal 
#' Riemann weights. 
#' 
#' \code{bconcurrent} implements a concurrent effect for a functional covariate
#' on a functional response, i.e. an effect of the form \eqn{X_i(t)\beta(t)} for
#' a functional response \eqn{Y_i(t)}. 
#' 
#' It is recommended to use centered functional covariates with 
#' \eqn{\sum_i X_i(s) = 0} for all \eqn{s} in \code{bconcurrent}- and 
#' \code{bsignal}-terms so that the global functional intercept 
#' can be interpreted as the global mean function. 
#' 
#' Cannot deal with any missing values in the covariates.
#' 
#' @return Equally to the base-learners of the package mboost: 
#' 
#' An object of class \code{blg} (base-learner generator) with a 
#' \code{dpp} function. 
#' 
#' The call of \code{dpp} returns an object of class 
#' \code{bl} (base-learner) with a \code{fit} function. The call to 
#' \code{fit} finally returns an object of class \code{bm} (base-model).
#' 
#' @seealso \code{\link{FDboost}} for the model fit. 
#' @keywords models
#' @references Scheipl, F., Staicu, A.-M., and Greven, S. (2014), 
#' Functional Additive Mixed Models, Journal of Computational and Graphical Statistics, 
#' in press, DOI 10.1080/10618600.2014.901914.
#' \url{http://arxiv.org/abs/1207.5947} 
#' @examples 
#' ### example for scalar response and two functional covariates 
#' data(fuel)
#' modFuel <- FDboost(heatan ~ bsignal(UVVIS, uvvis.lambda, knots=40, df=4) 
#'            + bsignal(NIR, nir.lambda, knots=40, df=4), 
#'            timeformula=~bols(1), data=fuel)
#' summary(modFuel)
#' ## plot(modFuel, rug=FALSE)
#' @export
### P-spline base-learner for signal matrix with index vector
bsignal <- function(...,  index = NULL, #by = NULL,
                    knots = 10, boundary.knots = NULL, degree = 3, differences = 2, df = 4, 
                    lambda = NULL, #center = FALSE, 
                    cyclic = FALSE, Z=NULL
){
  
  if (!is.null(lambda)) df <- NULL
  
  cll <- match.call()
  cll[[1]] <- as.name("bsignal")
  #print(cll)
  
  mfL <- list(...)
  if(length(mfL)>2) stop("bsignal has too many arguments")
  if(!is.matrix(mfL[[1]])) stop("signal has to be a matrix")
  
  varnames <- all.vars(cll)
  if(length(mfL)==1){ 
    mfL[[2]] <- 1:ncol(mfL[[1]]); cll[[3]] <- "xind" 
    varnames <- c(all.vars(cll), "xindDefault")
  }
  
  # Reshape mfL so that it is the dataframe of the signal with the index as attribute
  mf <- mfL[[1]]
  colnames(mf) <- paste(cll[[2]], 1:ncol(mf), sep="_")
  attr(mf, "signalIndex") <- mfL[[2]]
  xname <- varnames[1]
  indname <- varnames[2]
  attr(mf, "xname") <- xname
  attr(mf, "indname") <- indname
  
  mf <- data.frame("z"=I(mf))
  names(mf) <- as.character(cll[[2]])
  
#   if(all(round(colSums(mf), 4)!=0)){
#     warning(xname, " is not centered. 
#     Functional covariates should be mean-centered in each measurement point.")
#   }
  
  #   mf <- mfL
  #   names(mf) <- varnames
  
  vary <- ""
  
  CC <- all(mboost:::Complete.cases(mf))
  #  CC <- all(mboost:::Complete.cases(mf[1]))
  if (!CC)
    warning("base-learner contains missing values;\n",
            "missing values are excluded per base-learner, ",
            "i.e., base-learners may depend on different",
            " numbers of observations.")
  
  #index <- NULL
  
  ret <- list(model.frame = function() 
    if (is.null(index)) return(mf) else{
      mftemp <- mf
      mf <- mftemp[index,,drop = FALSE] # this is necessary to pass the attributes
      attributes(mftemp[,xname])
      attr(mf[,xname], "signalIndex") <- attr(mftemp[,xname], "signalIndex")
      attr(mf[,xname], "xname") <- attr(mftemp[,xname], "xname")
      attr(mf[,xname], "indname") <- attr(mftemp[,xname], "indname")
      return(mf)
    } ,
              get_call = function(){
                cll <- deparse(cll, width.cutoff=500L)
                if (length(cll) > 1)
                  cll <- paste(cll, collapse="")
                cll
              },
              get_data = function() mf,
              get_index = function() index,
              get_vary = function() vary,
              get_names = function(){
                attr(xname, "indname") <- indname 
                xname 
              }, #colnames(mf),
              set_names = function(value) {
                #if(length(value) != length(colnames(mf)))
                if(length(value) != names(mf[1]))
                  stop(sQuote("value"), " must have same length as ",
                       sQuote("names(mf[1])"))
                for (i in 1:length(value)){
                  cll[[i+1]] <<- as.name(value[i])
                }
                attr(mf, "names") <<- value
              })
  class(ret) <- "blg"
  
  #browser()
  #print("bsignal")
  #print(Z)
  
  ret$dpp <- mboost:::bl_lin(ret, Xfun = X_bsignal,
                    args = list(mf, vary, knots = knots, boundary.knots =
                      boundary.knots, degree = degree, differences = differences,
                                df = df, lambda = lambda, center = FALSE, cyclic = cyclic,
                                Z=Z))
  
  return(ret)
}

# ### hyper parameters for bsignal base-learner
# # add the parameter Z
# hyper_bsignal <- function(df = NULL, lambda = 0, intercept = TRUE,
#                        contrasts.arg = "contr.treatment", Z=NULL)
#   list(df = df, lambda = lambda,
#        intercept = intercept,
#        contrasts.arg = contrasts.arg,
#        Z=NULL)

# testX <- I(matrix(rnorm(40), ncol=5))
# s <- seq(0,1,l=5)
# test <- bsignal(testX, s)
# test <- bsignal(testX)
# test$get_names()
# test$get_data()
# #test$dpp(rep(1,nrow(testX)))

########################


#################################
# Base-learner for concurrent effect of functional covariate

### model.matrix for P-splines base-learner of signal matrix mf
X_conc <- function(mf, vary, args) {
  
  stopifnot(is.data.frame(mf))
  xname <- attr(mf[[1]], "xname")
  X1 <- as.matrix(mf[[1]])
  class(X1) <- "matrix"
  xind <- attr(mf[[1]], "signalIndex")
  
  #   stopifnot(is.list(mf))
  #   xname <- names(mf)[1]
  #   X1 <- mf[[1]]
  #   xind <- mf[[2]]
  
  if(ncol(X1)!=length(xind)) stop("Dimension of signal matrix and its index do not match.")
  
  ### Construct spline basis over index xind of X1 
  if(is.null(args$boundary.knots))  args$boundary.knots <- range(xind, na.rm = TRUE)
  knots <- seq(from = args$boundary.knots[1], to = args$boundary.knots[2], length = args$knots + 2)
  knots <- knots[2:(length(knots) - 1)]
  
  # B-spline basis of specified degree  
  #X <- bs(xind, knots=knots, degree=args$degree, intercept=TRUE) # old version   
  X <- mboost:::bsplines(xind, knots=knots, boundary.knots=args$boundary.knots, 
                         degree=args$degree) 
  
  #   # to do: extra feature: cyclic splines
  #   if (args$cyclic) {
  #     X <- mboost:::cbs(xind,
  #              knots = knots,
  #              boundary.knots = args$boundary.knots,
  #              degree = args$degree)
  #   }
  
  colnames(X) <- paste(xname, 1:ncol(X), sep="")
  
  # set up design matrix for concurrent model
  listCol <- list()
  for(i in 1:ncol(X1)){
    listCol[[i]] <- X1[,i]
  }
  X1des <- as.matrix(bdiag(listCol)) 
  
  # Design matrix is product of expanded X1 and basis expansion over xind 
  X <- (X1des) %*% X
  
  ### Penalty matrix: product differences matrix
  differenceMatrix <- diff(diag(ncol(X)), differences = args$differences)
  K <- crossprod(differenceMatrix)
  
  ## compare specified degrees of freedom to dimension of null space
  if (!is.null(args$df)){
    rns <- ncol(K) - qr(as.matrix(K))$rank # compute rank of null space
    if (rns == args$df)
      warning( sQuote("df"), " equal to rank of null space ",
               "(unpenalized part of P-spline);\n  ",
               "Consider larger value for ", sQuote("df"),
               " or set ", sQuote("center = TRUE"), ".", immediate.=TRUE)
    if (rns > args$df)
      stop("not possible to specify ", sQuote("df"),
           " smaller than the rank of the null space\n  ",
           "(unpenalized part of P-spline). Use larger value for ",
           sQuote("df"), " or set ", sQuote("center = TRUE"), ".")
  }
  return(list(X = X, K = K))
}

#' @rdname bsignal
#' @export
### P-spline base learner for signal matrix with index vector
bconcurrent <- function(..., #by = NULL, index = NULL, 
                        knots = 10, boundary.knots = NULL, degree = 3, differences = 2, df = 4 
                        #lambda = NULL, center = FALSE, cyclic = FALSE
){
  
  #  if (!is.null(lambda)) df <- NULL
  
  cll <- match.call()
  cll[[1]] <- as.name("bconcurrent")
  #print(cll)
  
  mfL <- list(...)
  if(length(mfL)>2) stop("bconcurrent has too many arguments")
  if(!is.matrix(mfL[[1]])) stop("signal has to be a matrix")
  
  varnames <- all.vars(cll)
  if(length(mfL)==1){ 
    mfL[[2]] <- 1:ncol(mfL[[1]]); cll[[3]] <- "xind" 
    varnames <- c(all.vars(cll), "xindDefault")
  }
  
  # Reshape mfL so that it is the dataframe of the signal with the index as attribute
  mf <- mfL[[1]]
  colnames(mf) <- paste(cll[[2]], 1:ncol(mf), sep="_")
  attr(mf, "signalIndex") <- mfL[[2]]
  xname <- varnames[1]
  indname <- varnames[2]
  attr(mf, "xname") <- xname
  attr(mf, "indname") <- indname
  
  mf <- data.frame("z"=I(mf))
  names(mf) <- as.character(cll[[2]]) 
  
#   if(all(round(colSums(mf, na.rm = TRUE), 4)!=0)){
#     warning(xname, " is not centered. 
#     Functional covariates should be mean-centered in each measurement point.")
#   }
  
  #   mf <- mfL
  #   names(mf) <- varnames
  
  vary <- ""
  
  CC <- all(mboost:::Complete.cases(mf))
  #  CC <- all(mboost:::Complete.cases(mf[1]))
  if (!CC)
    warning("base-learner contains missing values;\n",
            "missing values are excluded per base-learner, ",
            "i.e., base-learners may depend on different",
            " numbers of observations.")
  
  index <- NULL  
  
  ret <- list(model.frame = function() 
    if (is.null(index)) return(mf) else return(mf[index,,drop = FALSE]),
              get_call = function(){
                cll <- deparse(cll, width.cutoff=500L)
                if (length(cll) > 1)
                  cll <- paste(cll, collapse="")
                cll
              },
              get_data = function() mf,
              get_index = function() index,
              get_vary = function() vary,
              get_names = function(){
                attr(xname, "indname") <- indname 
                xname 
              }, #colnames(mf),
              set_names = function(value) {
                #if(length(value) != length(colnames(mf)))
                if(length(value) != names(mf[1]))
                  stop(sQuote("value"), " must have same length as ",
                       sQuote("names(mf[1])"))
                for (i in 1:length(value)){
                  cll[[i+1]] <<- as.name(value[i])
                }
                attr(mf, "names") <<- value
              })
  class(ret) <- "blg"
  
  ret$dpp <- mboost:::bl_lin(ret, Xfun = X_conc,
                             args = list(mf, vary, knots = knots, boundary.knots = boundary.knots, 
                                         degree = degree, differences = differences,
                                         df = df, lambda = NULL, center = FALSE, cyclic = FALSE))
  return(ret)
}

# testX <- I(matrix(rnorm(40), ncol=5))
# s <- seq(0,1,l=5)
# test <- bconcurrent(testX, s)
# #test <- bconcurrent(testX)
# test$get_names()
# test$get_data()
# #test$dpp(rep(1,nrow(testX)))






#################################
#### Base-learner for historic effect of functional covariate
### with integral over s<=t

### model.matrix for P-splines base-learner of signal matrix mf
X_hist <- function(mf, vary, args) {
  
  stopifnot(is.data.frame(mf))
  xname <- attr(mf[[1]], "xname")
  X1 <- as.matrix(mf[[1]])
  class(X1) <- "matrix"
  xind <- attr(mf[[1]], "signalIndex")
  yind <- attr(mf[[1]], "indexY")
  nobs <- nrow(X1)
  
  # get id-variable
  id <- attr(mf[,xname], "id") # for data in long format
  
  #   stopifnot(is.list(mf))
  #   xname <- names(mf)[1]
  #   X1 <- mf[[1]]
  #   xind <- mf[[2]]
  
  if(ncol(X1)!=length(xind)) stop("Dimension of signal matrix and its index do not match.")
  
  ### Construct spline basis over index xind of X1 
  if(is.null(args$boundary.knots))  args$boundary.knots <- range(xind, na.rm = TRUE)
  knots <- seq(from = args$boundary.knots[1], to = args$boundary.knots[2], length = args$knots + 2)
  knots <- knots[2:(length(knots) - 1)]
  
  # B-spline basis of specified degree     
  B.s <- mboost:::bsplines(xind, knots=knots, boundary.knots=args$boundary.knots, 
                           degree=args$degree) 
  
  colnames(B.s) <- paste(xname, 1:ncol(B.s), sep="")
  
  # Weighting with matrix of functional covariate
  L <- integrationWeightsLeft(X1=X1, xind=xind)
  X1L <- L*X1
  
#   # set up design matrix for historical model and s<=t with s and t equal to xind
#   # expand matrix of original observations to lower triangular matrix 
#   X1des0 <- matrix(0, ncol=ncol(X1), nrow=ncol(X1)*nrow(X1))
#   for(i in 1:ncol(X1des0)){
#     #print(nrow(X1)*(i-1)+1)
#     X1des0[(nrow(X1)*(i-1)+1):nrow(X1des0) ,i] <- X1L[,i] # use fun. variable * integration weights
#   }
  
  ## set up design matrix for historical model according to limit()
  # use the argument limits (Code taken of function ff(), package refund)
  limits <- args$limits
  if (!is.null(limits)) {
    if (!is.function(limits)) {
      if (!(limits %in% c("s<t", "s<=t"))) {
        stop("supplied <limits> argument unknown")
      }
      if (limits == "s<t") {
        limits <- function(s, t) {
          s < t
        }
      }
      else {
        if (limits == "s<=t") {
          limits <- function(s, t) {
            (s < t) | (s == t)
          }
        }
      }
    }
  }
  
  ### expand the design matrix for all observations (yind is equal for all observations!)
  ### the response is a vector (y1(t1), y2(t1), ... , yn(t1), yn(tG))
  if(is.null(id)){
    X1des <- X1L[rep(1:nobs, times=length(yind)), ]
  } else{ # yind is over all observations in long format
    X1des <- X1L[id, ] 
  }
  ### use function limits to set up design matrix according to function limits 
  ### by setting 0 at the time-points that should not be used
  if (!is.null(limits)) {
    ## at the moment: xind and yind are the same for all observations
    ## expand yind by replication to the yind of all observations together
    if(is.null(id)){
      ind0 <- !t(outer( xind, rep(yind, each=nobs), limits) )
    } else{ # yind is over all observations in long format
      ind0 <- !t(outer( xind, yind, limits) )
    } 
    X1des[ind0] <- 0
  }
  
  # Design matrix is product of expanded X1 and basis expansion over xind 
  X1des <- X1des %*% B.s
  
  # if xind and yind are equal B.t and B.s are equal as well
  if(length(xind)==length(yind) && all(xind==yind) ){
    # design matrix over index of response for one response
    B.t <- B.s
  } else{
    # design matrix over index of response for one response
    B.t <- mboost:::bsplines(yind, knots=knots, boundary.knots=args$boundary.knots, 
                          degree=args$degree) 
  }
  # stack design-matrix of response n times
  if(is.null(id)){
    B.t <- B.t[rep(1:length(yind), each=nobs), ]
  }  
  # in long format B.t is already over all time-points in the response
    
  # calculate row-tensor
  # X <- (X1 %x% t(rep(1, ncol(X2))) ) * ( t(rep(1, ncol(X1))) %x% X2  )
  X <- X1des[,rep(1:ncol(X1des), each=ncol(B.t))] * B.t[,rep(1:ncol(B.t), times=ncol(X1des))] 
  
  ### Penalty matrix: product differences matrix
  K1 <- diff(diag(ncol(X1des)), differences = args$differences)
  K1 <- crossprod(K1)
  K2 <- diff(diag(ncol(B.t)), differences = args$differences)
  K2 <- crossprod(K2)  
  K <- kronecker(K2, diag(ncol(X1des))) +
    kronecker(diag(ncol(B.t)), K1)
  
  ### <FIXME> necessary if K1 and K2 are not the same anyway?
  if(!is.null(id)){
    K <- kronecker(K1, diag(ncol(B.t))) +
      kronecker(diag(ncol(X1des)), K2)
  }
    
  
  ## compare specified degrees of freedom to dimension of null space
  if (!is.null(args$df)){
    rns <- ncol(K) - qr(as.matrix(K))$rank # compute rank of null space
    if (rns == args$df)
      warning( sQuote("df"), " equal to rank of null space ",
               "(unpenalized part of P-spline);\n  ",
               "Consider larger value for ", sQuote("df"),
               " or set ", sQuote("center = TRUE"), ".", immediate.=TRUE)
    if (rns > args$df)
      stop("not possible to specify ", sQuote("df"),
           " smaller than the rank of the null space\n  ",
           "(unpenalized part of P-spline). Use larger value for ",
           sQuote("df"), " or set ", sQuote("center = TRUE"), ".")
  }
  return(list(X = X, K = K))
}

### 
### model.matrix for P-splines base-learner of signal matrix mf
### for irregular response
X_hist2 <- function(mf, vary, args) {
  
  stopifnot(is.data.frame(mf))
  xname <- attr(mf[[1]], "xname")
  X1 <- as.matrix(mf[[1]])
  class(X1) <- "matrix"
  xind <- attr(mf[[1]], "signalIndex")
  yind <- attr(mf[[1]], "indexY")
  nobs <- nrow(X1)
  
  # get id-variable
  id <- attr(mf[,xname], "id") # for data in long format
    
  ###### EXTRA LINE in comparison to X_hist
  ## important for prediction, otherwise id=NULL and yind is multiplied accordingly
  if(is.null(id)) id <- 1:nrow(X1)
  
  ## check yind 
  stopifnot(length(yind)==length(id))
  
  #   stopifnot(is.list(mf))
  #   xname <- names(mf)[1]
  #   X1 <- mf[[1]]
  #   xind <- mf[[2]]
  
  if(ncol(X1)!=length(xind)) stop("Dimension of signal matrix and its index do not match.")
  
  ### Construct spline basis over index xind of X1 
  if(is.null(args$boundary.knots))  args$boundary.knots <- range(xind, na.rm = TRUE)
  knots <- seq(from = args$boundary.knots[1], to = args$boundary.knots[2], length = args$knots + 2)
  knots <- knots[2:(length(knots) - 1)]
  
  # B-spline basis of specified degree     
  B.s <- mboost:::bsplines(xind, knots=knots, boundary.knots=args$boundary.knots, 
                           degree=args$degree) 
  
  colnames(B.s) <- paste(xname, 1:ncol(B.s), sep="")
  
  # Weighting with matrix of functional covariate
  L <- integrationWeightsLeft(X1=X1, xind=xind)
  X1L <- L*X1
  
  #   # set up design matrix for historical model and s<=t with s and t equal to xind
  #   # expand matrix of original observations to lower triangular matrix 
  #   X1des0 <- matrix(0, ncol=ncol(X1), nrow=ncol(X1)*nrow(X1))
  #   for(i in 1:ncol(X1des0)){
  #     #print(nrow(X1)*(i-1)+1)
  #     X1des0[(nrow(X1)*(i-1)+1):nrow(X1des0) ,i] <- X1L[,i] # use fun. variable * integration weights
  #   }
  
  ## set up design matrix for historical model according to limit()
  # use the argument limits (Code taken of function ff(), package refund)
  limits <- args$limits
  if (!is.null(limits)) {
    if (!is.function(limits)) {
      if (!(limits %in% c("s<t", "s<=t"))) {
        stop("supplied <limits> argument unknown")
      }
      if (limits == "s<t") {
        limits <- function(s, t) {
          s < t
        }
      }
      else {
        if (limits == "s<=t") {
          limits <- function(s, t) {
            (s < t) | (s == t)
          }
        }
      }
    }
  }
  
  ### expand the design matrix for all observations (yind is equal for all observations!)
  ### the response is a vector (y1(t1), y2(t1), ... , yn(t1), yn(tG))
  if(is.null(id)){
    X1des <- X1L[rep(1:nobs, times=length(yind)), ]
  } else{ # yind is over all observations in long format
    X1des <- X1L[id, ] 
  }
  ### use function limits to set up design matrix according to function limits 
  ### by setting 0 at the time-points that should not be used
  if (!is.null(limits)) {
    if(is.null(id)){ # yind is the same for all observations
      ind0 <- !t(outer( xind, rep(yind, each=nobs), limits) )
    } else{ # yind is over all observations in long format
      ind0 <- !t(outer( xind, yind, limits) )
    } 
    X1des[ind0] <- 0
  }
  
  # Design matrix is product of expanded X1 and basis expansion over xind 
  X1des <- X1des %*% B.s
  
  # if xind and yind are equal B.t and B.s are equal as well
  if(length(xind)==length(yind) && all(xind==yind) ){
    # design matrix over index of response for one response
    B.t <- B.s
  } else{
    # design matrix over index of response for one response
    B.t <- mboost:::bsplines(yind, knots=knots, boundary.knots=args$boundary.knots, 
                             degree=args$degree) 
  }
  # stack design-matrix of response n times
  if(is.null(id)){
    B.t <- B.t[rep(1:length(yind), each=nobs), ]
  }  
  # in long format B.t is already over all time-points in the response

  # calculate row-tensor
  # X <- (X1 %x% t(rep(1, ncol(X2))) ) * ( t(rep(1, ncol(X1))) %x% X2  )
  X <- X1des[,rep(1:ncol(X1des), each=ncol(B.t))] * B.t[,rep(1:ncol(B.t), times=ncol(X1des))] 
  
  ### Penalty matrix: product differences matrix
  K1 <- diff(diag(ncol(X1des)), differences = args$differences)
  K1 <- crossprod(K1) 
  K2 <- diff(diag(ncol(B.t)), differences = args$differences)
  K2 <- crossprod(K2) 
  K <- kronecker(K2, diag(ncol(X1des))) +
    kronecker(diag(ncol(B.t)), K1)
  
  ### <FIXME> necessary if K1 and K2 are not the same anyway?
  if(!is.null(id)){
    K <- kronecker(K1, diag(ncol(B.t))) +
      kronecker(diag(ncol(X1des)), K2)
  }
  
  
  ## compare specified degrees of freedom to dimension of null space
  if (!is.null(args$df)){
    rns <- ncol(K) - qr(as.matrix(K))$rank # compute rank of null space
    if (rns == args$df)
      warning( sQuote("df"), " equal to rank of null space ",
               "(unpenalized part of P-spline);\n  ",
               "Consider larger value for ", sQuote("df"),
               " or set ", sQuote("center = TRUE"), ".", immediate.=TRUE)
    if (rns > args$df)
      stop("not possible to specify ", sQuote("df"),
           " smaller than the rank of the null space\n  ",
           "(unpenalized part of P-spline). Use larger value for ",
           sQuote("df"), " or set ", sQuote("center = TRUE"), ".")
  }
  return(list(X = X, K = K))
}


### P-spline base learner for signal matrix with index vector
### for historical model according to function limit, defaults to s<=t
#' @rdname bsignal
#' @export
bhist <- function(..., index = NULL, #by = NULL, 
                  knots = 10, boundary.knots = NULL, degree = 3, differences = 2, df = 4,
                  #lambda = NULL, center = FALSE, cyclic = FALSE
                  limits="s<=t"
){
  
  #  if (!is.null(lambda)) df <- NULL
    
  cll <- match.call()
  cll[[1]] <- as.name("bhist")
  #print(cll)
  
  mfL <- list(...)
  if(length(mfL)>3) stop("bhist has too many arguments")
  if(length(mfL)<3) stop("bhist has too few arguments")
  if(!is.matrix(mfL[[1]])) stop("signal has to be a matrix")
  
  varnames <- all.vars(cll)

  ### not necessary as user has to specify three arguments
#   # if no index is specified use an equidistant index in 0,1
#   if(length(mfL)==1){ 
#     mfL[[2]] <- 1:ncol(mfL[[1]]); cll[[3]] <- "xind" 
#     varnames <- c(all.vars(cll), "xindDefault")
#   }
#   # if no index for the response variable is specified use the index of the signal
#   if(length(mfL)==2){ 
#     mfL[[3]] <- mfL[[2]]; cll[[4]] <- cll[[3]]
#     varnames <- c(all.vars(cll), cll[[3]])
#   }
  
  if(!is.atomic(mfL[[2]])) stop("index of signal has to be a vector")
  if(!is.atomic(mfL[[3]])) stop("index of response has to be a vector")
  
  # compare range of index signal and index response
  # minimal value of the signal-index has to be smaller than the response-index
  if(min(mfL[[2]]) < min(mfL[[3]]) ) stop("Index of response has values before index of signal.")
  
  # Reshape mfL so that it is the dataframe of the signal with 
  # the index of the signal and the index of the response as attributes
  mf <- mfL[[1]]
  colnames(mf) <- paste(cll[[2]], 1:ncol(mf), sep="_")
  xname <- varnames[1]
  indname <- varnames[2]
  indnameY <- varnames[3]
  attr(mf, "xname") <- xname
  attr(mf, "indname") <- indname 
  attr(mf, "signalIndex") <- mfL[[2]]
  attr(mf, "indnameY") <- indnameY
  attr(mf, "indexY") <- mfL[[3]]
  attr(mf, "id") <- index
  
  mf <- data.frame("z"=I(mf))
  names(mf) <- as.character(cll[[2]]) 
    
  #   if(all(round(colSums(mf, na.rm = TRUE), 4)!=0)){
  #     warning(xname, " is not centered. 
  #     Functional covariates should be mean-centered in each measurement point.")
  #   }
  
  #   mf <- mfL
  #   names(mf) <- varnames
  
  vary <- ""
  
  CC <- all(mboost:::Complete.cases(mf))
  #  CC <- all(mboost:::Complete.cases(mf[1]))
  if (!CC)
    warning("base-learner contains missing values;\n",
            "missing values are excluded per base-learner, ",
            "i.e., base-learners may depend on different",
            " numbers of observations.")
  
  #index <- NULL 
  #browser()
    
  ret <- list(model.frame = function() 
    if (is.null(index)) return(mf) else{
      mftemp <- mf
      mf <- mftemp[index,,drop = FALSE] # this is necessary to pass the attributes
      attributes(mftemp[,xname])
      attr(mf[,xname], "signalIndex") <- attr(mftemp[,xname], "signalIndex")
      attr(mf[,xname], "xname") <- attr(mftemp[,xname], "xname")
      attr(mf[,xname], "indname") <- attr(mftemp[,xname], "indname")
      attr(mf[,xname], "id") <- attr(mftemp[,xname], "id")
      return(mf)
    } ,
              get_call = function(){
                cll <- deparse(cll, width.cutoff=500L)
                if (length(cll) > 1)
                  cll <- paste(cll, collapse="")
                cll
              },
              get_data = function() mf,
              ## set index to NULL, as the index is treated within X_hist()
              ##get_index = function() index, 
              get_index = function() NULL,
              get_vary = function() vary,
              get_names = function(){
                attr(xname, "indname") <- indname 
                attr(xname, "indnameY") <- indnameY
                xname 
              }, #colnames(mf),
              set_names = function(value) {
                #if(length(value) != length(colnames(mf)))
                if(length(value) != names(mf[1]))
                  stop(sQuote("value"), " must have same length as ",
                       sQuote("names(mf[1])"))
                for (i in 1:length(value)){
                  cll[[i+1]] <<- as.name(value[i])
                }
                attr(mf, "names") <<- value
              })
  class(ret) <- "blg"
  
  if(is.null(index)){
    ### X-hist is for data in wide format with regular response
    ret$dpp <- mboost:::bl_lin(ret, Xfun = X_hist,
                               args = list(mf, vary, knots = knots, boundary.knots = boundary.knots, 
                                           degree = degree, differences = differences,
                                           df = df, lambda = NULL, center = FALSE, cyclic = FALSE,
                                           limits = limits)) # extra feature limit!
  }else{
    ### X_hist2 is for data in long format: sets id 1:n
    ret$dpp <- mboost:::bl_lin(ret, Xfun = X_hist2,
                               args = list(mf, vary, knots = knots, boundary.knots = boundary.knots, 
                                           degree = degree, differences = differences,
                                           df = df, lambda = NULL, center = FALSE, cyclic = FALSE,
                                           limits = limits)) # extra feature limit! 
  }
  #browser()
  return(ret)
}

# testX <- I(matrix(rnorm(40), ncol=5))
# s <- seq(0,1,l=5)
# test <- bhist(testX, s, knots=5)
# #test <- bhist(testX)
# test$get_names()
# test$get_data()
# #test$dpp(rep(1,nrow(testX)))
# extract(test)








#################################
# Base-learner with constraints for smooth varying scalar covariate

### almost equal to X_bbs() of package mboost
### difference: implements sum-to-zero-constraint over index of response
X_bbsc <- function(mf, vary, args) {
  
  stopifnot(is.data.frame(mf))
  mm <- lapply(which(colnames(mf) != vary), function(i) {
    X <- mboost:::bsplines(mf[[i]],
                  knots = args$knots[[i]]$knots,
                  boundary.knots = args$knots[[i]]$boundary.knots,
                  degree = args$degree)
    if (args$cyclic) {
      X <- mboost:::cbs(mf[[i]],
               knots = args$knots[[i]]$knots,
               boundary.knots = args$knots[[i]]$boundary.knots,
               degree = args$degree)
    }
    class(X) <- "matrix"
    return(X)
  })  ### options
  MATRIX <- any(sapply(mm, dim) > c(500, 50)) || (length(mm) > 1)
  MATRIX <- MATRIX && options("mboost_useMatrix")$mboost_useMatrix
  if (MATRIX) {
    diag <- Diagonal
    cbind <- cBind
    for (i in 1:length(mm)){
      tmp <- attributes(mm[[i]])[c("degree", "knots", "Boundary.knots")]
      mm[[i]] <- Matrix(mm[[i]])
      attributes(mm[[i]])[c("degree", "knots", "Boundary.knots")] <- tmp
    }
  }
  
  if (length(mm) == 1) {
    X <- mm[[1]]
    if (vary != "") {
      by <- model.matrix(as.formula(paste("~", vary, collapse = "")),
                         data = mf)[ , -1, drop = FALSE] # drop intercept
      DM <- lapply(1:ncol(by), function(i) {
        ret <- X * by[, i]
        colnames(ret) <- paste(colnames(ret), colnames(by)[i], sep = ":")
        ret
      })
      if (is(X, "Matrix")) {
        X <- do.call("cBind", DM)
      } else {
        X <- do.call("cbind", DM)
      }
    }
    if (args$differences > 0){
      if (!args$cyclic) {
        K <- diff(diag(ncol(mm[[1]])), differences = args$differences)
      } else {
        ## cyclic P-splines
        differences <- args$differences
        K <- diff(diag(ncol(mm[[1]]) + differences),
                  differences = differences)
        tmp <- K[,(1:differences)]   # save first "differences" columns
        K <- K[,-(1:differences)]    # drop first "differences" columns
        indx <- (ncol(mm[[1]]) - differences + 1):(ncol(mm[[1]]))
        K[,indx] <- K[,indx] + tmp   # add first "differences" columns
      }
    } else {
      if (args$differences != 0)
        stop(sQuote("differences"), " must be an non-neative integer")
      K <- diag(ncol(mm[[1]]))
    }
    
    if (vary != "" && ncol(by) > 1){       # build block diagonal penalty
      K <- kronecker(diag(ncol(by)), K)
    }
    if (args$center) {
      tmp <- attributes(X)[c("degree", "knots", "Boundary.knots")]
      X <- tcrossprod(X, K) %*% solve(tcrossprod(K))
      attributes(X)[c("degree", "knots", "Boundary.knots")] <- tmp
      K <- diag(ncol(X))
    } else {
      K <- crossprod(K)
    }
  }
  
  if (length(mm) == 2) {
    suppressMessages(
      X <- kronecker(mm[[1]], matrix(1, ncol = ncol(mm[[2]]))) *
        kronecker(matrix(1, ncol = ncol(mm[[1]])), mm[[2]])
    )
    if (vary != "") {
      by <- model.matrix(as.formula(paste("~", vary, collapse = "")),
                         data = mf)[ , -1, drop = FALSE] # drop intercept
      DM <- X * by[,1]
      if (ncol(by) > 1){
        for (i in 2:ncol(by))
          DM <- cbind(DM, (X * by[,i]))
      }
      X <- DM
      ### <FIXME> Names of X if by is given
    }
    if (args$differences > 0){
      if (!args$cyclic) {
        Kx <- diff(diag(ncol(mm[[1]])), differences = args$differences)
        Ky <- diff(diag(ncol(mm[[2]])), differences = args$differences)
      } else {
        ## cyclic P-splines
        differences <- args$differences
        Kx <- diff(diag(ncol(mm[[1]]) + differences),
                   differences = differences)
        Ky <- diff(diag(ncol(mm[[2]]) + differences),
                   differences = differences)
        
        tmp <- Kx[,(1:differences)]   # save first "differences" columns
        Kx <- Kx[,-(1:differences)]    # drop first "differences" columns
        indx <- (ncol(mm[[1]]) - differences + 1):(ncol(mm[[1]]))
        Kx[,indx] <- Kx[,indx] + tmp   # add first "differences" columns
        
        tmp <- Ky[,(1:differences)]   # save first "differences" columns
        Ky <- Ky[,-(1:differences)]    # drop first "differences" columns
        indx <- (ncol(mm[[2]]) - differences + 1):(ncol(mm[[2]]))
        Ky[,indx] <- Ky[,indx] + tmp   # add first "differences" columns
      }
    } else {
      if (args$differences != 0)
        stop(sQuote("differences"), " must be an non-negative integer")
      Kx <- diag(ncol(mm[[1]]))
      Ky <- diag(ncol(mm[[2]]))
    }
    
    Kx <- crossprod(Kx)
    Ky <- crossprod(Ky)
    suppressMessages(
      K <- kronecker(Kx, diag(ncol(mm[[2]]))) +
        kronecker(diag(ncol(mm[[1]])), Ky)
    )
    if (vary != "" && ncol(by) > 1){       # build block diagonal penalty
      suppressMessages(K <- kronecker(diag(ncol(by)), K))
    }
    if (!identical(args$center, FALSE)) {
      ### L = \Gamma \Omega^1/2 in Section 2.3. of Fahrmeir et al.
      ### (2004, Stat Sinica), always
      L <- eigen(K, symmetric = TRUE, EISPACK = FALSE)
      L$vectors <- L$vectors[,1:(ncol(X) - args$differences^2), drop = FALSE]
      L$values <- sqrt(L$values[1:(ncol(X) - args$differences^2), drop = FALSE])
      L <- L$vectors %*% (diag(length(L$values)) * (1/L$values))
      X <- as(X %*% L, "matrix")
      K <- as(diag(ncol(X)), "matrix")
    }
  }
  
  if (length(mm) > 2)
    stop("not possible to specify more than two variables in ",
         sQuote("..."), " argument of smooth base-learners")
    
  #----------------------------------
  ### <SB> Calculate constraints

  # If the argument Z is not NULL use the given Z (important for prediction!)
  if(is.null(args$Z)){
    C <- t(X) %*% rep(1, nrow(X))
    Q <- qr.Q(qr(C), complete=TRUE) # orthonormal matrix of QR decomposition
    Z <- Q[  , 2:ncol(Q)] # only keep last columns    
  }else Z <- args$Z
  
  ### Transform design and penalty matrix
  X <- X %*% Z
  K <- t(Z) %*% K %*% Z
  
  attr(X, "Z") <- Z # attribute is not used at the moment
  #----------------------------------
    
  ## compare specified degrees of freedom to dimension of null space
  if (!is.null(args$df)){
    rns <- ncol(K) - qr(as.matrix(K))$rank # compute rank of null space
    if (rns == args$df)
      warning( sQuote("df"), " equal to rank of null space ",
               "(unpenalized part of P-spline);\n  ",
               "Consider larger value for ", sQuote("df"),
               " or set ", sQuote("center = TRUE"), ".", immediate.=TRUE)
    if (rns > args$df)
      stop("not possible to specify ", sQuote("df"),
           " smaller than the rank of the null space\n  ",
           "(unpenalized part of P-spline). Use larger value for ",
           sQuote("df"), " or set ", sQuote("center = TRUE"), ".")
  }
  return(list(X = X, K = K, Z = Z))
}


###############################################################################

#' Constrained Base-learners for Scalar Covariates
#' 
#' Constrained base-learners for fitting effects of scalar covariates in models 
#' with functional response
#' 
#' @param ... one or more predictor variables or one matrix or data 
#' frame of predictor variables. 
#' @param by an optional variable defining varying coefficients, 
#' either a factor or numeric variable.
#' @param index a vector of integers for expanding the variables in \code{...}.
#' @param knots either the number of knots or a vector of the positions 
#' of the interior knots (for more details see \code{\link[mboost]{bbs}}).
#' @param boundary.knots boundary points at which to anchor the B-spline basis 
#' (default the range of the data). A vector (of length 2) 
#' for the lower and the upper boundary knot can be specified.
#' @param degree degree of the regression spline.
#' @param differences a non-negative integer, typically 1, 2 or 3. 
#' If \code{differences} = \emph{k}, \emph{k}-th-order differences are used as 
#' a penalty (\emph{0}-th order differences specify a ridge penalty).
#' @param df trace of the hat matrix for the base-learner defining the 
#' base-learner complexity. Low values of \code{df} correspond to a 
#' large amount of smoothing and thus to "weaker" base-learners.
#' @param lambda smoothing penalty, computed from \code{df} when 
#' \code{df} is specified.
#' @param K in \code{bolsc} it is possible to specify the penalty matrix K
#' @param center not implemented yet 
#' @param cyclic  if \code{cyclic = TRUE} the fitted values coincide at 
#' the boundaries (useful for cyclic covariates such as day time etc.).
#' @param contrasts.arg Note that a special \code{contrasts.arg} exists in 
#' package \code{mboost}, namely "contr.dummy". This contrast is used per default 
#' in \code{brandomc}. It leads to a 
#' dummy coding as returned by \code{model.matrix(~ x - 1)} were the 
#' intercept is implicitly included but each factor level gets a 
#' separate effect estimate (for more details see \code{\link[mboost]{brandom}}).
#' @param intercept if intercept = TRUE an intercept is added to the design matrix 
#' of a linear base-learner. 
#' 
#' @details The base-learners \code{bbsc}, \code{bolsc} and \code{brandomc} are 
#' basically the base-learners \code{\link[mboost]{bbs}}, \code{\link[mboost]{bols}} and 
#' \code{\link[mboost]{brandom}} with additional identifiability constraints. 
#' Instead of the default identifiability constraints 
#' (\eqn{\sum_{i,t} \hat f(x_i, t) = 0}) 
#' implemented in \code{mboost} for tensor product terms whose 
#' marginal terms include the index of the functional 
#' response \eqn{t} constraints that enforce 
#' \eqn{\sum_i \hat f(z_i, x_i, t) = 0} for all \eqn{t} are used, so that 
#' effects varying over \eqn{t} can be interpreted as deviations 
#' from the global functional intercept. 
#' 
#' Cannot deal with any missing values in the covariates.
#' 
#' @return Equally to the base-learners of the package mboost: 
#' 
#' An object of class \code{blg} (base-learner generator) with a 
#' \code{dpp} function. 
#' 
#' The call of \code{dpp} returns an object of class 
#' \code{bl} (base-learner) with a \code{fit} function. The call to 
#' \code{fit} finally returns an object of class \code{bm} (base-model).
#' 
#' @seealso \code{\link{FDboost}} for the model fit. 
#' \code{\link[mboost]{bbs}}, \code{\link[mboost]{bols}} and \code{\link[mboost]{brandom}} for the 
#' corresponding base-learners in mboost.
#' @references Scheipl, F., Staicu, A.-M., and Greven, S. (2014), 
#' Functional Additive Mixed Models, Journal of Computational and Graphical Statistics, 
#' in press, DOI 10.1080/10618600.2014.901914.
#' \url{http://arxiv.org/abs/1207.5947}
#' @keywords models
#' @aliases brandomc
#' @export
bbsc <- function(..., by = NULL, index = NULL, knots = 10, boundary.knots = NULL,
                degree = 3, differences = 2, df = 4, lambda = NULL, center = FALSE,
                cyclic = FALSE) {
  
  if (!is.null(lambda)) df <- NULL
  
  cll <- match.call()
  cll[[1]] <- as.name("bbsc")
  
  mf <- list(...)
  if (length(mf) == 1 && ((is.matrix(mf[[1]]) || is.data.frame(mf[[1]])) &&
    ncol(mf[[1]]) > 1 )) {
    mf <- as.data.frame(mf[[1]])
  } else {
    mf <- as.data.frame(mf)
    cl <- as.list(match.call(expand.dots = FALSE))[2][[1]]
    colnames(mf) <- sapply(cl, function(x) deparse(x))
  }
  stopifnot(is.data.frame(mf))
  if(!(all(sapply(mf, is.numeric)))) {
    if (ncol(mf) == 1){
      warning("cannot compute ", sQuote("bbsc"),
              " for non-numeric variables; used ",
              sQuote("bols"), " instead.")
      return(bols(mf, by = by, index = index))
    }
    stop("cannot compute bbsc for non-numeric variables")
  }
  vary <- ""
  if (!is.null(by)){
    mf <- cbind(mf, by)
    colnames(mf)[ncol(mf)] <- vary <- deparse(substitute(by))
  }
  
  CC <- all(mboost:::Complete.cases(mf))
  if (!CC)
    warning("base-learner contains missing values;\n",
            "missing values are excluded per base-learner, ",
            "i.e., base-learners may depend on different",
            " numbers of observations.")
  ### option
  DOINDEX <- (nrow(mf) > options("mboost_indexmin")[[1]])
  if (is.null(index)) {
    if (!CC || DOINDEX) {
      index <- mboost:::get_index(mf)
      mf <- mf[index[[1]],,drop = FALSE]
      index <- index[[2]]
    }
  }
  
  ret <- list(model.frame = function()
    if (is.null(index)) return(mf) else return(mf[index,,drop = FALSE]),
              get_call = function(){
                cll <- deparse(cll, width.cutoff=500L)
                if (length(cll) > 1)
                  cll <- paste(cll, collapse="")
                cll
              },
              get_data = function() mf,
              get_index = function() index,
              get_vary = function() vary,
              get_names = function() colnames(mf),
              set_names = function(value) {
                if(length(value) != length(colnames(mf)))
                  stop(sQuote("value"), " must have same length as ",
                       sQuote("colnames(mf)"))
                for (i in 1:length(value)){
                  cll[[i+1]] <<- as.name(value[i])
                }
                attr(mf, "names") <<- value
              })
  class(ret) <- "blg"
  
  ret$dpp <- mboost:::bl_lin(ret, Xfun = X_bbsc,
                    args = mboost:::hyper_bbs(mf, vary, knots = knots, boundary.knots =
                      boundary.knots, degree = degree, differences = differences,
                      df = df, lambda = lambda, center = center, cyclic = cyclic))
  return(ret)
}

# z2 <- rnorm(17)
# blz <- bbsc(z=z2, df=3)
# blz$get_call()
# blz$get_names()
# str(blz$get_data())
# str(blz$dpp(weights=rep(1,17)))

#################
### model.matrix for constrained ols base-learner with penalty matrix K
X_olsc <- function(mf, vary, args) {
  
  if (mboost:::isMATRIX(mf)) {
    X <- mf
    contr <- NULL
  } else {
    ### set up model matrix
    fm <- paste("~ ", paste(colnames(mf)[colnames(mf) != vary],
                            collapse = "+"), sep = "")
    fac <- sapply(mf[colnames(mf) != vary], is.factor)
    if (any(fac)){
      if (!is.list(args$contrasts.arg)){
        ## first part needed to prevent warnings from calls such as
        ## contrasts.arg = contr.treatment(4, base = 1):
        if (is.character(args$contrasts.arg) &&
              args$contrasts.arg == "contr.dummy"){
          if (!args$intercept)
            stop('"contr.dummy" can only be used with ',
                 sQuote("intercept = TRUE"))
          fm <- paste(fm, "-1")
        } else {
          txt <- paste("list(", paste(colnames(mf)[colnames(mf) != vary][fac],
                                      "= args$contrasts.arg", collapse = ", "),")")
          args$contrasts.arg <- eval(parse(text=txt))
        }
      } else {
        ## if contrasts are given as list check if "contr.dummy" is specified
        if (any(args$contrasts.arg == "contr.dummy"))
          stop('"contr.dummy"',
               " can only be used for all factors at the same time.\n",
               "Use ", sQuote('contrasts.arg = "contr.dummy"'),
               " to achieve this.")
      }
    } else {
      args$contrasts.arg <- NULL
    }
    X <- model.matrix(as.formula(fm), data = mf, contrasts.arg = args$contrasts.arg)
    if (!is.null(args$contrasts.arg) && args$contrasts.arg == "contr.dummy")
      attr(X, "contrasts") <- lapply(attr(X, "contrasts"),
                                     function(x) x <- "contr.dummy")
    contr <- attr(X, "contrasts")
    if (!args$intercept)
      X <- X[ , -1, drop = FALSE]
    MATRIX <- any(dim(X) > c(500, 50)) && any(fac)
    MATRIX <- MATRIX && options("mboost_useMatrix")$mboost_useMatrix
    if (MATRIX) {
      diag <- Diagonal
      cbind <- cBind
      if (!is(X, "Matrix"))
        X <- Matrix(X)
    }
    if (vary != "") {
      by <- model.matrix(as.formula(paste("~", vary, collapse = "")),
                         data = mf)[ , -1, drop = FALSE] # drop intercept
      DM <- lapply(1:ncol(by), function(i) {
        ret <- X * by[, i]
        colnames(ret) <- paste(colnames(ret), colnames(by)[i], sep = ":")
        ret
      })
      if (is(X, "Matrix")) {
        X <- do.call("cBind", DM)
      } else {
        X <- do.call("cbind", DM)
      }
    }
  }
  
  #----------------------------------
  ## <SB> use given penalty-matrix K
  if(is.null(args$K)){
    ### <FIXME> penalize intercepts???
    ### set up penalty matrix
    ANOVA <- (!is.null(contr) && (length(contr) == 1)) && (ncol(mf) == 1)
    K <- diag(ncol(X))
    ### for ordered factors use difference penalty
    if (ANOVA && any(sapply(mf[, names(contr), drop = FALSE], is.ordered))) {
      K <- diff(diag(ncol(X) + 1), differences = 1)[, -1, drop = FALSE]
      if (vary != "" && ncol(by) > 1){       # build block diagonal penalty
        suppressMessages(K <- kronecker(diag(ncol(by)), K))
      }
      K <- crossprod(K)
    }
  }else{
    ## <SB> check dimensions of given K
    stopifnot(dim(args$K)==rep(ncol(X), 2))
    K <- args$K
  }
  #----------------------------------
  
  #----------------------------------
  ### <SB> Calculate constraints
  
  # If the argument Z is not NULL use the given Z (important for prediction!)
  if(is.null(args$Z)){
    C <- t(X) %*% rep(1, nrow(X))
    Q <- qr.Q(qr(C), complete=TRUE) # orthonormal matrix of QR decompositon
    Z <- Q[  , 2:ncol(Q)] # only keep last columns    
  }else Z <- args$Z
  
  ### Transform design and penalty matrix
  X <- X %*% Z
  K <- t(Z) %*% K %*% Z
  
  attr(X, "Z") <- Z # attr not used at the moment
  #----------------------------------
  
  ### </FIXME>
  if (is(X, "Matrix") && !is(K, "Matrix"))
    K <- Matrix(K)
  
  ### <SB> return the transformation matrix Z as well
  list(X = X, K = K, Z = Z)
}


#' @rdname bbsc
#' @export
### Linear base-learner, potentially Ridge-penalized (but not by default)
### one can specify the penalty matrix K
### with sum-to-zero constraint over index of response
bolsc <- function(..., by = NULL, index = NULL, intercept = TRUE, df = NULL,
                 lambda = 0, K=NULL, contrasts.arg = "contr.treatment") {
  
  if (!is.null(df)) lambda <- NULL
  
  cll <- match.call()
  cll[[1]] <- as.name("bolsc")
  
  mf <- list(...)
  if (length(mf) == 1 && ((mboost:::isMATRIX(mf[[1]]) || is.data.frame(mf[[1]])) &&
                            ncol(mf[[1]]) > 1 )) {
    mf <- mf[[1]]
    ### spline bases should be matrices
    if (mboost:::isMATRIX(mf) && !is(mf, "Matrix"))
      class(mf) <- "matrix"
  } else {
    mf <- as.data.frame(mf)
    cl <- as.list(match.call(expand.dots = FALSE))[2][[1]]
    colnames(mf) <- sapply(cl, function(x) as.character(x))
  }
  if(!intercept && !any(sapply(mf, is.factor)) &&
       !any(sapply(mf, function(x){uni <- unique(x);
                                   length(uni[!is.na(uni)])}) == 1)){
    ## if no intercept is used and no covariate is a factor
    ## and if no intercept is specified (i.e. mf[[i]] is constant)
    if (any(sapply(mf, function(x) abs(mean(x, na.rm=TRUE) / sd(x,na.rm=TRUE))) > 0.1))
      ## if covariate mean is not near zero
      warning("covariates should be (mean-) centered if ",
              sQuote("intercept = FALSE"))
  }
  vary <- ""
  if (!is.null(by)){
    stopifnot(is.data.frame(mf))
    mf <- cbind(mf, by)
    colnames(mf)[ncol(mf)] <- vary <- deparse(substitute(by))
  }
  
  CC <- all(mboost:::Complete.cases(mf))
  if (!CC)
    warning("base-learner contains missing values;\n",
            "missing values are excluded per base-learner, ",
            "i.e., base-learners may depend on different",
            " numbers of observations.")
  ### option
  DOINDEX <- is.data.frame(mf) &&
    (nrow(mf) > options("mboost_indexmin")[[1]] || is.factor(mf[[1]]))
  if (is.null(index)) {
    ### try to remove duplicated observations or
    ### observations with missings
    if (!CC || DOINDEX) {
      index <- mboost:::get_index(mf)
      mf <- mf[index[[1]],,drop = FALSE]
      index <- index[[2]]
    }
  }
  
  ret <- list(model.frame = function()
    if (is.null(index)) return(mf) else return(mf[index,,drop = FALSE]),
              get_call = function(){
                cll <- deparse(cll, width.cutoff=500L)
                if (length(cll) > 1)
                  cll <- paste(cll, collapse="")
                cll
              },
              get_data = function() mf,
              get_index = function() index,
              get_names = function() colnames(mf),
              get_vary = function() vary,
              set_names = function(value) {
                if(length(value) != length(colnames(mf)))
                  stop(sQuote("value"), " must have same length as ",
                       sQuote("colnames(mf)"))
                for (i in 1:length(value)){
                  cll[[i+1]] <<- as.name(value[i])
                }
                attr(mf, "names") <<- value
              })
  class(ret) <- "blg"
  
  ret$dpp <- mboost:::bl_lin(ret, Xfun = X_olsc, args = hyper_olsc(
    df = df, lambda = lambda, K = K, # use penalty matrix as argument
    intercept = intercept, contrasts.arg = contrasts.arg,
    Z=NULL)) # Z in args not used at the moment
  return(ret)
}

### hyper parameters for olsc base-learner
# add the parameters Z and K
hyper_olsc <- function(df = NULL, lambda = 0, K=NULL, intercept = TRUE,
                      contrasts.arg = "contr.treatment", Z=NULL)
  list(df = df, lambda = lambda, K=K,
       intercept = intercept,
       contrasts.arg = contrasts.arg,
       Z=NULL)


#' @rdname bbsc
#' @export
# random-effects (Ridge-penalized ANOVA) base-learner
# almost equal to brandom, but with sum-to-zero-constraint over index of t
brandomc <- function (..., contrasts.arg = "contr.dummy", df = 4) {
  cl <- cltmp <- match.call()
  if (is.null(cl$df))
    cl$df <- df
  if (is.null(cl$contrasts.arg))
    cl$contrasts.arg <- contrasts.arg
  cl[[1L]] <- as.name("bolsc")
  ret <- eval(cl, parent.frame())
  cltmp[[1]] <- as.name("brandomc")
  assign("cll", cltmp, envir = environment(ret$get_call))
  ret
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



