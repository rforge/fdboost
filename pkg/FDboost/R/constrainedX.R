#' Constrained row tensor product 
#' 
#' EXPERIMENTAL! 
#' 
#' Combining single base-learners to form new, more complex base-learners, with
#' an identifiability constraint to center the interaction around the intercept and
#' around the two main effects. Suitable for functional response. 
#' @param bl1 base-learner 1, e.g. \code{bols(x1)}
#' @param bl2 base-learner 2, e.g. \code{bols(x2)}
#' 
#' @details Similar to \code{\%X\%} in package mboost, see \code{\link[mboost]{\%X\%}}, 
#' a row tensor product of linear base-learners is returned by \code{\%Xc\%}. 
#' \code{\%Xc\%} applies a sum-to-zero constraint to the design matrix suitable for
#' functional response if an interaction of two scalar covariates is specified 
#' in the case that the model contains a global intercept and both main effects, 
#' as the interaction is centerd around the intercept and centered around the two main effects. 
#' See Web Appendix A of Brockhaus et al. (2015) for details on how to enforce the constraint 
#' for the functional intercept.   
#' Use e.g. in a model call to \code{FDboost}, following the scheme, 
#' \code{y ~ 1 + bolsc(x1) + bolsc(x2) + bols(x1) \%Xc\% bols(x2)}, 
#' where \code{1} induces a global intercept and \code{x1}, \code{x2} are factor variables.  
#' 
#' @references 
#' Brockhaus, S., Scheipl, F., Hothorn, T. and Greven, S. (2015): 
#' The functional linear array model. Statistical Modelling, 15(3), 279-300.
#' 
#' @author Sarah Brockhaus, David Ruegamer
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
#' ## fit model with interaction that is centered around the intercept 
#' ## and the two main effects 
#' mod1 <- FDboost(vis ~ 1 + bolsc(T_C, df=1) + bolsc(T_A, df=1) + 
#'                 bols(T_C, df=2) %Xc% bols(T_A, df=1),
#'                 timeformula = ~bbs(time, df=3),
#'                 numInt = "equal", family = QuantReg(),
#'                 offset = NULL, offset_control = o_control(k_min = 9),
#'                 data = viscosity, control=boost_control(mstop = 100, nu = 0.4))
#'                 
#' ## check centering around intercept
#' colMeans(predict(mod1, which = 4))
#' 
#' ## check centering around main effects
#' colMeans(predict(mod1, which = 4)[viscosity$T_A == "low", ])
#' colMeans(predict(mod1, which = 4)[viscosity$T_A == "high", ])
#' colMeans(predict(mod1, which = 4)[viscosity$T_C == "low", ])
#' colMeans(predict(mod1, which = 4)[viscosity$T_C == "low", ])
#'
#' ## find optimal mstop using cvrsik() or validateFDboost()
#' ## ... 
#' 
#' ## look at interaction effect in one plot
#' # funplot(mod1$yind, predict(mod1, which=4))
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
  
  ## Check that the used base-learners contain an intercept
  used_bl <- c( deparse(match.call()$bl1[[1]]), 
                deparse(match.call()$bl2[[1]]) )
  if(any(used_bl == "bolsc")) stop("Use bols instead of bolsc with %Xc%.")
  if(any(used_bl == "brandomc")) stop("Use brandom instead of brandomc with %Xc%.")
  if( (!is.null(match.call()$bl1$intercept) &&  match.call()$bl1$intercept != TRUE) |
      (!is.null(match.call()$bl2$intercept) &&  match.call()$bl2$intercept != TRUE) ){
    stop("Set intercept = TRUE in base-learners used with %Xc%.")
  }
  
  if(any(!used_bl %in% c("bols", "brandom")) ){
    warning("%Xc% is only tested for base-learners bols and brandom with factor variables!")
  }
  
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
    
    ## <SB> set prediciton to FALSE, if it is NULL
    if(is.null(args$prediction)) args$prediction <- FALSE
    
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
    
    ## If model is fitted -> compute Z; but if model is predicted use the Z from the model fit
    ## if(!args$prediction){
    ## compute QR-decompotition only once
    if(is.null(args$Z)){
      ## put all effects of the two main effects + intercept into the constraints 
      C <- t(X) %*% cbind(rep(1, nrow(X)), X1[ , -1], X2[ , -1])
      qr_C <- qr(C)
      if( any(class(qr_C) == "sparseQR") ){
        rank_C <- qr_C@Dim[2]
      }else{
        rank_C <- qr_C$rank 
      } 
      Q <- qr.Q(qr_C, complete=TRUE) # orthonormal matrix of QR decomposition
      args$Z <- Q[  , (rank_C + 1) : ncol(Q) ] # only keep last columns    
    } 
    
    ### Transform design and penalty matrix
    X <- X %*% args$Z
    K <- t(args$Z) %*% K %*% args$Z
    ## print(args$Z)
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
#' EXPERIMENTAL! 
#' 
#' Ignores weights that are specified in the model-call for the computation of
#' the anisotropic penalty!
#' @param bl1 base-learner 1, e.g. \code{bbs(x1)}
#' @param bl2 base-learner 2, e.g. \code{bbs(x2)}
#' 
#' @details EXPERIMENTAL! 
#' 
#' When \code{\%O\%} is called with a specification of \code{df} in both base-learners, 
#' e.g. \code{bbs(x1, df = df1) \%O\% bbs(t, df = df2)}, the global \code{df} for the 
#' Kroneckered base-learner is computed as \code{df = df1 * df2}. 
#' And thus the penalty has only one smoothness parameter lambda resulting in an isotropic penalty, 
#' \deqn{P = lambda * [(P1 o I) + (I o P2)],} 
#' with overall penalty \eqn{P}, Kronecker product \eqn{o}, 
#' marginal penalty matrices \eqn{P1, P2} and identity matrices \eqn{I}.  
#' (Currie et al. (2006) introduced the generalized linear array model, which has a design matrix that 
#' is composed of the Kronecker product of two marginal design matrices, which was implemented in mboost 
#' as \code{\%O\%}.  
#' See Brockhaus et al. (2015) for the application of array models to functional data.)  
#' 
#' In contrast, a Kronecker product with anisotropic penalty is obtained by \code{\%A\%}, 
#' which allows for a different amount of smoothness in the two directions. 
#' For example \code{bbs(x1, df = df1) \%A\% bbs(t, df = df2)} results in computing two
#' different values for lambda for the two marginal design matrices and a global value of 
#' lambda to adjust for the global \code{df}, i.e. 
#' \deqn{P = lambda * [(lambda1 * P1 o I) +  (I o lambda2 * P2)],} 
#' with Kronecker product \eqn{o}, 
#' where \eqn{lambda1} is computed individually for \eqn{df1} and \eqn{P1}, 
#' \eqn{lambda2} is computed individually for \eqn{df2}  and \eqn{P2}, 
#' and \eqn{lambda} is computed such that the global \eqn{df} hold \eqn{df = df1 * df2}. 
#' For the computation of \eqn{lambda1} and \eqn{lambda2} weights specified in the model 
#' call are not used, this implies that \eqn{lambda1} and \eqn{lambda2} are equal over 
#' folds of \code{cvrisk}. The computation of the global \eqn{lambda} considers the 
#' specified \code{weights}, such the global \eqn{df} are correct.        
#' 
#' If the \code{formula} in \code{FDboost} contains base-learners connected by 
#' \code{\%O\%} or \code{\%A\%}, 
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
#' @examples  
#' ######## Example for anisotropic penalty  
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
#' ## isotropic penalty, as timeformula is kroneckered to each effect using %O%...
#' mod1 <- FDboost(vis ~ 1 + 
#'                 bolsc(T_C, df=1) + 
#'                 bolsc(T_A, df=1) + 
#'                 bols(T_C, df=1) %Xc% bols(T_A, df=1),
#'                 timeformula = ~ bbs(time, df=3),
#'                 numInt = "equal", family = QuantReg(),
#'                 offset = NULL, offset_control = o_control(k_min = 9),
#'                 data = viscosity, control=boost_control(mstop = 200, nu = 0.4))
#' ## cf. the formula that is passed to mboost
#' mod1$formulaMboost
#' 
#' ## anisotropic effects using %A%
#' mod1a <- FDboost(vis ~ 1 + 
#'                 bolsc(T_C, df=1) %A% bbs(time, df=3) + 
#'                 bolsc(T_A, df=1) %A% bbs(time, df=3) + 
#'                 bols(T_C, df=1) %Xc% bols(T_A, df=1) %A% bbs(time, df=3),
#'                 timeformula = ~ bbs(time, df=3),
#'                 numInt = "equal", family = QuantReg(),
#'                 offset = NULL, offset_control = o_control(k_min = 9),
#'                 data = viscosity, control=boost_control(mstop = 200, nu = 0.4)) 
#'                 
#' ## optimize mstop for mod1 and mod1a
#' ## ...
#'                 
#' ## compare estimated coefficients
#' \dontrun{
#' par(mfrow=c(4, 2))
#' plot(mod1, which = 1)
#' plot(mod1a, which = 1)
#' plot(mod1, which = 2)
#' plot(mod1a, which = 2)
#' plot(mod1, which = 3)
#' plot(mod1a, which = 3)
#' funplot(mod1$yind, predict(mod1, which=4))
#' funplot(mod1$yind, predict(mod1a, which=4))
#' }
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
    
    ### <SB> anisotropic penalty matrix 
    df1 <- args1$df
    df2 <- args2$df
    
    args <- list(lambda = NULL,
                 df = ifelse(is.null(df1), 1, df1) *
                   ifelse(is.null(df2), 1, df2))

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



