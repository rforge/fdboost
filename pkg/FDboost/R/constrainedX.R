#' Constrained row tensor product 
#' 
#' EXPERIMENTAL! Combining single base-learners to form new, more complex base-learners, with
#' an identifiability constraint suitable for functional response. 
#' @param bl1 base-learner 1, e.g. \code{bols(x1)}
#' @param bl2 base-learner 2, e.g. \code{bols(x2)}
#' 
#' @details Similar to \code{\%X\%} in package mboost, see \code{\link[mboost]{\%X\%}}, a tensor product of two 
#' or more linear base-learners is returned by \code{\%Xc\%}. 
#' But \code{\%Xc\%} applies a sum-to-zero constraint to the design matrix suitable for
#' functional response if an interaction of two scalar covariates is specified. Before the constraint 
#' is applied an intercept-column is added to the design matrix, thus the two base-learners that are 
#' connected by \code{\%Xc\%} should both not contain an intercept. 
#' 
#' @examples  
#' ######## Example for function-on-scalar-regression with interaction effect of two scalar covariates 
#' data("viscosity", package = "FDboost") 
#' ## set time-interval that should be modeled
#' interval <- "101"
#' 
#' ## model time until "interval" and take log() of viscosity
#' end <- which(viscosity$timeAll == as.numeric(interval))
#' viscosity$vis <- log(viscosity$visAll[,1:end])
#' viscosity$time <- viscosity$timeAll[1:end]
#' # with(viscosity, funplot(time, vis, pch = 16, cex = 0.2))
#' 
#' ## fit median regression model with 100 boosting iterations,
#' ## step-length 0.4 and smooth time-specific offset
#' ## the factors are coded such that the intercept is the global median 
#' ## no integration weights are used!
#' mod1 <- FDboost(vis ~ 1 + bolsc(T_C, df=2) + bolsc(T_A, df=2) + 
#'                 bolsc(T_C, df=2) %Xc% bolsc(T_A, df=1),
#'                 timeformula = ~bbs(time, df=3),
#'                 numInt = "equal", family = QuantReg(),
#'                 offset = NULL, offset_control = o_control(k_min = 9),
#'                 data = viscosity, control=boost_control(mstop = 100, nu = 0.4))
#' 
#' @export
"%Xc%" <- function(bl1, bl2) {
  
  if (is.list(bl1) && !inherits(bl1, "blg"))
    return(lapply(bl1, "%Xc%", bl2 = bl2))
  
  if (is.list(bl2) && !inherits(bl2, "blg"))
    return(lapply(bl2, "%Xc%", bl1 = bl1))
  
  cll <- paste(bl1$get_call(), "%Xc%",
               bl2$get_call(), collapse = "")
  stopifnot(inherits(bl1, "blg"))
  stopifnot(inherits(bl2, "blg"))
  
  stopifnot(!any(colnames(model.frame(bl1)) %in%
                   colnames(model.frame(bl2))))
  mf <- cbind(model.frame(bl1), model.frame(bl2))
  index1 <- bl1$get_index()
  index2 <- bl2$get_index()
  if (is.null(index1)) index1 <- 1:nrow(mf)
  if (is.null(index2)) index2 <- 1:nrow(mf)
  
  mfindex <- cbind(index1, index2)
  index <- NULL
  
  CC <- all(Complete.cases(mf))
  if (!CC)
    warning("base-learner contains missing values;\n",
            "missing values are excluded per base-learner, ",
            "i.e., base-learners may depend on different",
            " numbers of observations.")
  ### option
  DOINDEX <- (nrow(mf) > options("mboost_indexmin")[[1]])
  if (is.null(index)) {
    if (!CC || DOINDEX) {
      index <- mboost:::get_index(mfindex)
      mf <- mf[index[[1]],,drop = FALSE]
      index <- index[[2]]
    }
  }
  
  vary <- ""
  
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
    ## <FIXME> Is this all we want to change if we set names here?
    set_names = function(value) attr(mf, "names") <<- value)
  ## </FIXME>
  class(ret) <- "blg"
  
  args1 <- environment(bl1$dpp)$args
  args2 <- environment(bl2$dpp)$args
  l1 <- args1$lambda
  l2 <- args2$lambda
  if (!is.null(l1) && !is.null(l2)) {
    args <- list(lambda = 1, df = NULL)
  } else {
    args <- list(lambda = NULL,
                 df = ifelse(is.null(args1$df), 1, args1$df) *
                   ifelse(is.null(args2$df), 1, args2$df))
  }
  
  Xfun <- function(mf, vary, args) {
    
    newX1 <- environment(bl1$dpp)$newX
    newX2 <- environment(bl2$dpp)$newX
    
    X1 <- newX1(mf[, bl1$get_names(), drop = FALSE],
                prediction = args$prediction)
    K1 <- X1$K
    X1 <- X1$X
    if (!is.null(l1)) K1 <- l1 * K1
    MATRIX <- options("mboost_useMatrix")$mboost_useMatrix
    if (MATRIX & !is(X1, "Matrix"))
      X1 <- Matrix(X1)
    if (MATRIX & !is(K1, "Matrix"))
      K1 <- Matrix(K1)
    
    X2 <- newX2(mf[, bl2$get_names(), drop = FALSE],
                prediction = args$prediction)
    K2 <- X2$K
    X2 <- X2$X
    if (!is.null(l2)) K2 <- l2 * K2
    if (MATRIX & !is(X2, "Matrix"))
      X2 <- Matrix(X2)
    if (MATRIX & !is(K2, "Matrix"))
      K2 <- Matrix(K2)
    suppressMessages(
      X <- kronecker(X1, Matrix(1, ncol = ncol(X2),
                                dimnames = list("", colnames(X2))),
                     make.dimnames = TRUE) *
        kronecker(Matrix(1, ncol = ncol(X1),
                         dimnames = list("", colnames(X1))),
                  X2, make.dimnames = TRUE)
    )
    suppressMessages(
      K <- kronecker(K1, diag(ncol(X2))) +
        kronecker(diag(ncol(X1)), K2)
    )
    
    
    ## add an intercept column to the design matrix
    X <- cbind(1, X)
    ## fixme what to to with the penalty??
    K <- bdiag(K[1,1], K)
    
    #----------------------------------
    ### <SB> Calculate constraints
    if(is.null(args$prediction)) args$prediction <- FALSE
    
    ## If model is fitted -> compute Z; but if model is predicted use the Z from the model fit
    ## if(!args$prediction){
    ## compute QR-decompotition only once
    if(is.null(args$Z)){
      C <- t(X) %*% rep(1, nrow(X))
      Q <- qr.Q(qr(C), complete=TRUE) # orthonormal matrix of QR decomposition
      args$Z <- Q[  , 2:ncol(Q)] # only keep last columns    
    } 
    
    ### Transform design and penalty matrix
    X <- X %*% args$Z
    K <- t(args$Z) %*% K %*% args$Z
    #print(args$Z)
    #----------------------------------
    
    list(X = X, K = K, args = args)  
  }
  
  ## compute the transformation matrix Z
  temp <- Xfun(mf = mf, vary = vary, args = args)
  args$Z <- temp$args$Z
  rm(temp)
  
  ret$dpp <- bl_lin(ret, Xfun = Xfun, args = args)
  
  return(ret)
}



#' Kronecker product of two base-learners with anisotropic penalty 
#' 
#' EXPERIMENTAL! Ignores weights that are specified in the model-call!
#' @param bl1 base-learner 1, e.g. \code{bbs(x1)}
#' @param bl2 base-learner 2, e.g. \code{bbs(x2)}
#' 
#' @details EXPERIMENTAL! 
#' 
#' When \code{\%O\%} is called with a specification of \code{df} in both base-learners, 
#' e.g. \code{bbs(x1, df = df1) \%O\% bbs(t, df = df2)}, the global \code{df} for the 
#' Kroneckered base-learner is computed as \code{df = df1 * df2}. 
#' And thus the penalty has only one smoothness parameter lambda resulting in an isotropic penalty, 
#' \eqn{P = lambda*[( P1 o I ) +  (I o P2)]}, with overall penalty \eqn{P}, Kronecker product \eqn{o}, 
#' marginal penalty matrices \eqn{P1, P2} and identity matrices \eqn{I}.  
#' (Currie et al. (2006) introduced the generalized linear array model, which has a design matrix that 
#' is composed of the Kronecker product of two marginal design matrices and 
#' see Brockhaus et al. (2015) for the application of array models to functional data.)  
#' 
#' A Kronecker product with anisotropic penalty is obtained by \code{\%A\%}, which allows for 
#' a different amount of smoothness in the two directions. 
#' For example \code{bbs(x1, df = df1) \%A\% bbs(t, df = df2)} results in computing, two
#' different values for lambda for the two marginal design matrices and a global value of 
#' lambda to adjust for the global \code{df}, i.e. 
#' \eqn{P = lambda*[( lambda1*P1 o I ) +  (I o lambda2*P2) ]}, with Kronecker product \eqn{o}, 
#' where \eqn{lambda1} is computed for \eqn{df1}, \eqn{lambda2} is computed for \eqn{df2}, 
#' and \eqn{lambda} is computed such that the global \eqn{df} hold \eqn{df = df1 * df2}.   
#' 
#' If the formula in \code{FDboost} contains base-learners connected by \code{\%O\%} or \code{\%A\%}, 
#' those effects are not expanded with \code{timeformula}, allowing for model specifications 
#' with different effects in time-direction.  
#' 
#' @references 
#' Brockhaus, S., Scheipl, F., Hothorn, T. and Greven, S. (2015): 
#' The functional linear array model. Statistical Modelling, 15(3), 279-300. 
#' 
#' Currie, I.D., Durban, M. and Eilers P.H.C. (2006):  
#' Generalized linear array models with applications to multidimensional smoothing. 
#' Journal of the Royal Statistical Society, Series B-Statistical Methodology, 68(2), 259-280.
#' 
#' @export
"%A%" <- function(bl1, bl2) {
  
  #   if (is.list(bl1) && !inherits(bl1, "blg"))
  #     return(lapply(bl1, "%X%", bl2 = bl2))
  #   
  #   if (is.list(bl2) && !inherits(bl2, "blg"))
  #     return(lapply(bl2, "%X%", bl1 = bl1))
  
  cll <- paste(bl1$get_call(), "%A%",
               bl2$get_call(), collapse = "")
  stopifnot(inherits(bl1, "blg"))
  stopifnot(inherits(bl2, "blg"))
  
  mf1 <- model.frame(bl1)
  mf2 <- model.frame(bl2)
  stopifnot(!any(colnames(mf1) %in%
                   colnames(mf2)))
  mf <- c(mf1, mf2)
  stopifnot(all(complete.cases(mf[[1]])))
  stopifnot(all(complete.cases(mf[[2]])))
  
  index <- NULL
  
  vary <- ""
  
  ret <- list(model.frame = function()
    return(mf),
    get_call = function(){
      cll <- deparse(cll, width.cutoff=500L)
      if (length(cll) > 1)
        cll <- paste(cll, collapse="")
      cll
    },
    get_data = function() mf,
    get_index = function() index,
    get_vary = function() vary,
    get_names = function() names(mf),
    ## <FIXME> Is this all we want to change if we set names here?
    set_names = function(value) attr(mf, "names") <<- value)
  ## </FIXME>
  class(ret) <- "blg"
  
  args1 <- environment(bl1$dpp)$args
  args2 <- environment(bl2$dpp)$args
  l1 <- args1$lambda
  l2 <- args2$lambda
  if (xor(is.null(l1), is.null(l2)))
    stop("you cannot mix lambda and df in ",
         sQuote("%A%"))
  if (!is.null(l1) && !is.null(l2)) {
    ### there is no common lambda!
    args <- list(lambda = NA, df = NA)
  } else {
    
    args <- list(lambda = NULL,
                 df = ifelse(is.null(args1$df), 1, args1$df) *
                   ifelse(is.null(args2$df), 1, args2$df))
    
    ### <SB> anisotropic penalty matrix 
    df1 <- args1$df
    df2 <- args2$df

    ## case that df equals nr columns of design matrix -> no penalty -> lambda = 0
    if( abs(ncol(environment(bl1$dpp)$X) - df1) < .Machine$double.eps*10^10){ 
      args$lambda1 <- 0 
    }else{
      ## <FIXME> specify weights!! important if model fit is only on part of data
      al1 <- df2lambda(X = environment(bl1$dpp)$X,
                       df = df1, lambda = NULL,
                       dmat = environment(bl1$dpp)$K, weights = 1,  XtX = NULL)
      args$lambda1 <- al1["lambda"]
    }

    
    ## case that df equals nr columns of design matrix
    if( abs(ncol(environment(bl2$dpp)$X) - df2) < .Machine$double.eps*10^10){ 
      args$lambda2 <- 0 
    }else{
      ## <FIXME> specify weights!! important if model fit is only on part of data
      al2 <- df2lambda(X = environment(bl2$dpp)$X,
                       df = df2, lambda = NULL,
                       dmat = environment(bl2$dpp)$K, weights = 1,  XtX = NULL)
      args$lambda2 <- al2["lambda"]
    }
    
    l1 <- args$lambda1
    l2 <- args$lambda2 

  }
  
  Xfun <- function(mf, vary, args) {

    newX1 <- environment(bl1$dpp)$newX
    newX2 <- environment(bl2$dpp)$newX
    
    X1 <- newX1(as.data.frame(mf[bl1$get_names()]),
                prediction = args$prediction)
    K1 <- X1$K
    X1 <- X1$X
    if (!is.null(l1)) K1 <- l1 * K1
    MATRIX <- options("mboost_useMatrix")$mboost_useMatrix
    if (MATRIX & !is(X1, "Matrix"))
      X1 <- Matrix(X1)
    if (MATRIX & !is(K1, "Matrix"))
      K1 <- Matrix(K1)
    
    X2 <- newX2(as.data.frame(mf[bl2$get_names()]),
                prediction = args$prediction)
    K2 <- X2$K
    X2 <- X2$X
    if (!is.null(l2)) K2 <- l2 * K2
    if (MATRIX & !is(X2, "Matrix"))
      X2 <- Matrix(X2)
    if (MATRIX & !is(K2, "Matrix"))
      K2 <- Matrix(K2)
    suppressMessages(
      K <- kronecker(K2, diag(ncol(X1))) +
        kronecker(diag(ncol(X2)), K1)
    )
    list(X = list(X1 = X1, X2 = X2), K = K)
  }
  
  ret$dpp <- mboost:::bl_lin_matrix(ret, Xfun = Xfun, args = args)
  
  return(ret)
}



