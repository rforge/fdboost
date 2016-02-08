#' Constrained row tensor product 
#' 
#' Combining single base-learners to form new, more complex base-learners, with
#' an identifiability constraint suitable for functional response. 
#' @param bl1 base-learner 1, e.g. \code{bols(x1)}
#' @param bl2 base-learner 2, e.g. \code{bols(x2)}
#' 
#' @details Similar to \code{\%X\%} in package mboost, see ?"\%X\%", a tensor product of two 
#' or more linear base-learners is returned by \code{\%Xc\%}. 
#' But \code{\%Xc\%} applies a sum-to-zero constraint to the design matrix suitable for
#' functional response if an interaction of two scalar covariates is specified.
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
#'                 bols(T_C, df=2) %Xc% bols(T_A, df=1),
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
      index <- get_index(mfindex)
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

