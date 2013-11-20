#################################################################################
#' Model-based Gradient Boosting for Functional Response
#' 
#' Gradient boosting for optimizing arbitrary loss functions, where component-wise models 
#' are utilized as base-learners in the case of functional response.
#' This function is a wrapper for \code{mboost}'s \code{\link[mboost]{mboost}} and its 
#' siblings to fit models of the general form 
#' \cr \eqn{E(Y_i(t)) = g(\mu(t) + \int X_i(s)\beta(s,t)ds + f(z_{1i}, t) + 
#' f(z_{2i}) + z_{3i} \beta_3(t) + \dots }\cr
#' with a functional (but not necessarily continuous) response \eqn{Y(t)}, response function \eqn{g},
#' (optional) smooth intercept \eqn{\mu(t)}, (multiple) functional covariates \eqn{X(t)} 
#' and scalar covariates \eqn{z_1}, \eqn{z_2}, etc. 
#' 
#' @details The functional response and functional covariates have to be
#' supplied as n by <no. of evaluations> matrices, i.e. each row is one
#' functional observation. For the model fit the matrix of the functional
#' response evaluations \eqn{Y_i(t)} are stacked into one long vector.
#' 
#' @param formula a symbolic description of the model to be fit.
#' @param timeformula formula for the expansion over the index of the response. 
#' For a functional response \eqn{Y_i(t)} typically ~bbs(t) to obtain a smooth 
#' expansion of the effects along t.
#' @param numInt integration scheme for the integration of the loss function.
#' One of \code{c("equal", "Riemann")} meaning equal weights of 1 or 
#' trapezoidal Riemann weights.
#' Alternativley a vector of length \code{nrow(response)} containing arbitray 
#' positive weights can be specified.
#' @param data a data frame or list containing the variables in the model.
#' @param weights (1) a numeric vector of weights for observational units, 
#' i.e. length(weights) has to be nrow(response),
#' (2) alternatively weights can be specified for single observations then
#' length(weights) has to be \code{nrow(response)}*\code{ncol(response)}
#' per default weights is constantly 1. 
#' @param offset_control parameters for the calculation of the offset, 
#' defaults to \code{offset_control(k_min=20, silent=TRUE)}.  
#' @param offset a numeric vector to be used as offset over the index of the response (optional).
#' If no offset is specified, per default a smooth time-specific offset is calculated and used 
#' within the model fit. If you do not want to use an offset you can set offset=0.
#' @param ... additional arguments passed to \code{\link[mboost]{mboost}}, 
#' including \code{offset}, \code{family} and \code{control}.
#' @return on object of class \code{FDboost} that inherits from \code{mboost}.
#' Special \code{\link{predict.FDboost}}, \code{\link{coef.FDboost}} and 
#' \code{\link{plot.FDboost}} methods are available. 
#' The methods of \code{\link[mboost]{mboost}} ara available as well, 
#' e.g. \code{\link[mboost]{extract}}
#' 
#' @author Sarah Brockhaus, Torsten Hothorn
#' @seealso \code{\link[mboost]{mboost}} for the help of the wrapped function in package mboost.
#' @references Ivanescu, A., Staicu, A.-M., Scheipl, F. and Greven, S. (2012). 
#'  Penalized function-on-function regression. (submitted) 
#' \url{http://biostats.bepress.com/jhubiostat/paper240/}
#' @keywords models, nonlinear 
#' @examples  
#' ###############################################################################
#' print("to do")
#' ###############################################################################
#' @export
#' @import mboost Matrix 
#' @importFrom splines bs splineDesign
#' @importFrom mgcv gam s
#' @importFrom zoo na.locf
#' @importFrom nnls nnls
FDboost <- function(formula,          ### response ~ xvars
                    timeformula,      ### time
                    numInt="equal",   ### option for approximation of integral over loss
                    data,             ### list of response, time, xvars
                    weights = NULL,   ### optional
                    offset = NULL,    ### optional
                    offset_control = o_control(),
                    #offset_control = list(k_min=20, silent=TRUE),
                    ...)              ### goes directly to mboost
{
  dots <- list(...)

  ### save formula of FDboost
  formulaFDboost <- formula
  
  stopifnot(class(formula)=="formula")
  stopifnot(class(timeformula)=="formula")
   
  ### extract response; a numeric matrix (for the moment)
  yname <- all.vars(formula)[1]
  response <- data[[yname]]
  data[[yname]] <- NULL
  
  ### extract time;
  # <FixMe> Only keep first variable, 
  # otherwise it is impossible to have a vector for knots
  yind <- all.vars(timeformula)[[1]]
  stopifnot(length(yind) == 1)
  nameyind <- yind
  assign(yind, data[[yind]])
  time <- data[[yind]]
  stopifnot(is.numeric(time))
  data[[yind]] <- NULL
  attr(time, "nameyind") <- nameyind
  
  ### extract covariates
  # data <- as.data.frame(data)
  if(length(all.vars(formula)) > 1){
    data <- data[all.vars(formula)[-1]]    
  }else data <- list(NULL)  # <SB> intercept-model without covariates
      
  ### get covariates that are modelled constant over time
  # code of function pffr() of package refund
  tf <- terms.formula(formula, specials=c("c"))
  trmstrings <- attr(tf, "term.labels")
  terms <- sapply(trmstrings, function(trm) as.call(parse(text=trm))[[1]], simplify=FALSE) 
  #ugly, but getTerms(formula)[-1] does not work for terms like I(x1:x2) 
  frmlenv <- environment(formula)
  where.c <- attr(tf, "specials")$c - 1    # indices of scalar offset terms
  
  #transform: c(foo) --> foo
  if(length(where.c)){ 
    trmstrings[where.c] <- sapply(trmstrings[where.c], function(x){
      sub("\\)$", "", sub("^c\\(", "", x)) #c(BLA) --> BLA
    })
    blconstant <- "bols(ONEtime, intercept = FALSE)"
    assign("ONEtime", rep(1.0, length(time)))
  }
    
  ### check dimensions
  ### response has trajectories as rows
  stopifnot(is.matrix(response))
  # <SB> dataframe is list and can contain time-points of functional covariables of arbitrary length
  # if (nrow(data) > 0) stopifnot(nrow(response) == nrow(data))
  nr <- nrow(response)
  stopifnot(ncol(response) == length(time))
  nc <- ncol(response)
  dresponse <- as.vector(response) # column-wise stacking of response
  
  ### save original dimensions of response
  ydim <- dim(response)
  
  ### intercept
  assign("ONEx", rep(1.0, nr))
  
#   ### FIXME: plausibility check: 
#   # for constrained effects in the formula the model should include an intercept
#   grep("bbsc", trmstrings) + grep("brandomc", trmstrings)
  
  ### compose mboost formula
  cfm <- paste(deparse(formula), collapse = "") 
  cfm <- strsplit(cfm, "~")[[1]]
  #xfm <- strsplit(cfm[2], "\\+")[[1]]
  xfm <- trmstrings
  ### replace "1" with intercept base learner
  if ( any( gsub(" ", "", strsplit(cfm[2], "\\+")[[1]]) ==  "1")){
    if( grepl("lambda", deparse(timeformula)) || 
         ( grepl("bols", deparse(timeformula)) &  !grepl("df", deparse(timeformula))) ){
      xfm <- c("bols(ONEx, intercept = FALSE)", xfm)
    } else{
      xfm <- c("bols(ONEx, intercept = FALSE, df=1)", xfm)
    }
    where.c <- where.c + 1
  }
    
  yfm <- strsplit(cfm[1], "\\+")[[1]]
  tfm <- paste(deparse(timeformula), collapse = "")
  tfm <- strsplit(tfm, "~")[[1]]
  tfm <- strsplit(tfm[2], "\\+")[[1]]
  
  ## set up formula for effects constant in time
  if(length(where.c) > 0){
    ### set c_df to the df/lambda in timeformula
    ##<FIXME> Does this make sense for bols() base-learner?
    if( grepl("lambda", tfm) || 
          ( grepl("bols", tfm) &  !grepl("df", tfm)) ){
      c_lambda <- eval(parse(text=paste(tfm, "$dpp(rep(1.0,", length(time), "))$df()", sep="")))["lambda"]
      cfm <- paste("bols(ONEtime, intercept = FALSE, lambda = ", c_lambda ,")")
    } else{
      c_df <- eval(parse(text=paste(tfm, "$dpp(rep(1.0,", length(time), "))$df()", sep="")))["df"]
      cfm <- paste("bols(ONEtime, intercept = FALSE, df = ", c_df ,")")
    }
  }

  tmp <- outer(xfm, tfm, function(x, y) paste(x, y, sep = "%O%"))
  # do not expand an effect bconcurrent() by a Kronecker product
  if( length(c(grep("bconcurrent", tmp), grep("bhis", tmp)) ) > 0 ) 
    tmp[c(grep("bconcurrent", tmp), grep("bhis", tmp))] <- xfm[c(grep("bconcurrent", tmp), grep("bhis", tmp))]   
  if(length(where.c) > 0) 
    tmp[where.c] <- outer(xfm[where.c], cfm, function(x, y) paste(x, y, sep = "%O%"))
  xpart <- paste(as.vector(tmp), collapse = " + ")
  fm <- as.formula(paste("dresponse ~ ", xpart))
  
  ### expand weights for observations
  if (is.null(weights)) weights <- rep(1, nr)
  w <- weights
  if (length(w) == nr) w <- rep(w, nc) # expand weights if they are only on the columns
  if(length(w)!=nc*nr) stop("Dimensions of weights do not match the dimensions of the response.") # check dimensions of w 
  
  ### multiply integration weights numInt to weights and w
  if(is.numeric(numInt)){
    if(length(numInt)!=length(time)) stop("Length of integration weights and time vector are not equal.")
    weights <- weights*numInt
    w <- rep(weights, each=nr)
  }else{
    if(!numInt %in% c("equal", "Riemann")) warning("argument numInt is ignored as it is neither numeric nor one of (\"equal\", \"Riemann\")")
    if(numInt=="Riemann"){ 
      w <- w*as.vector(integrationWeights(X1=response, time))
    }
  }
  
  ### set weights of missing values to 0
  if(sum(is.na(dresponse))>0){
    #warning(paste("The response contains", sum(is.na(dresponse)) ,"missing values. The corresponding weights are set to 0."))
    w[which(is.na(dresponse))] <- 0
  }
  
    
  ## per default add smooth time-specific offset 
  ## idea: allow to use an offset linear in time?
  if(is.null(offset)){
    message("Use a smooth offset.") 
    ### <FixMe> is the use of family@offset correct?
    #meanY <- colMeans(response, na.rm=TRUE)
    if(! "family" %in% names(dots) ){
      offsetFun <- Gaussian()@offset
    } else offsetFun <- dots$family@offset
    meanY <- c()
    # do a linear interpolation of the response to prevent bias because of missing values
    # only use responses with less than 90% missings for the calculation of the offset
    # only use response curves whose weights are not completely 0 (important for resampling methods)
    meanNA <- apply(response, 1, function(x) mean(is.na(x)))
    responseInter <- t(apply(response[meanNA < 0.9 & rowSums(matrix(w, ncol=nc))!=0 , ], 1, function(x) 
      approx(time, x, rule=offset_control$rule, xout=time)$y))
    # check wether first ore last columns of the response contain soley NA
    # then use the values of the next column
    if(any(apply(responseInter, 2, function(x) all(is.na(x)) ) )){
      warning("Column of interpolated response contains nothing but NA.")
      allNA <- apply(responseInter, 2, function(x) all(is.na(x)) ) 
      allNAlower <- allNAupper <- allNA
      allNAupper[1:round(ncol(responseInter)/2)] <- FALSE # missing columns with low index 
      allNAlower[round(ncol(responseInter)/2):ncol(responseInter)] <- FALSE # missing columns with high index 
      responseInter[, allNAlower] <- responseInter[, max(which(allNAlower))+1]
      responseInter[, allNAupper] <- responseInter[, min(which(allNAlower))-1]
    }
  
    for(i in 1:nc){
      try(meanY[i] <- offsetFun(responseInter[,i], 1*!is.na(responseInter[,i]) ), silent=TRUE)
    }
    if( is.null(meanY) ||  any(is.na(meanY)) ){
      warning("Mean offset cannot be computed by family@offset(). Use a weighted mean instead.")
      meanY <- c()
      for(i in 1:nc){
        meanY[i] <- Gaussian()@offset(responseInter[,i], 1*!is.na(responseInter[,i]) )
      }
    }   
    ### <FixMe> is the computation of k ok? 
    modOffset <- try( gam(meanY ~ s(time, bs="ad", 
                                    k = min(offset_control$k_min, round(length(time)/2))  )), 
                      silent=offset_control$silent )
    if(any(class(modOffset)=="try-error")){
      warning("Could not fit the smooth offset by adaptive splines (default), 
              use a simple spline expansion with 5 df instead.")
      if(round(length(time)/2) < 8) warning("Most likely because of too few time-points.")
      modOffset <- lm(meanY ~ bs(time, df=5))
    } 
    offsetVec <- modOffset$fitted.values 
    predictOffset <- function(time){
      ret <- predict(modOffset, newdata=data.frame(time=time))
      names(ret) <- NULL
      ret      
    } 
    offset <- as.vector(matrix(offsetVec, ncol=ncol(response), nrow = nrow(response), byrow=TRUE))    
  }else{
    if(length(offset)==1 && offset==0){
      offsetVec <- NULL
      offset <- NULL # use one constant offset in mboost()
      predictOffset <- function(time) 0      
    }else{ # expand the offset to the long vector like dresponse
      if(length(offset)!=nc) stop("Dimensions of offset and response do not match.")
      offsetVec <- offset
      offset <- as.vector(matrix(offset, ncol=ncol(response), nrow = nrow(response), byrow=TRUE))
      ### <FixMe> use a more sophisticated model to estimate the time-specific offset? 
      modOffset <- lm(offsetVec ~ bs(time, df=length(offsetVec)-2))
      predictOffset <- function(time){
        ret <- predict(modOffset, newdata=data.frame(time=time))
        names(ret) <- NULL
        ret      
      }      
    }   
  }
      
  if (length(data) > 0) {
    ### mboost isn't happy with nrow(data) == 0
    ret <- mboost(fm, data = data, weights = w, offset=offset, ...) 
  } else {
    ret <- mboost(fm, weights = w, offset=offset, ...)
  }
  
  ### assign new class (e.g. for specialized predictions)
  class(ret) <- c("FDboost", class(ret))
  
  ### reset weights for cvrisk etc., expanding works OK in bl_lin_matrix!
  # ret$"(weights)" <- weights
  # <SB> do not reset weights as than the integration weights are lost
  
  ret$yname <- yname
  ret$ydim <- ydim  
  ret$yind <- time
  ret$data <- data
    
  # if the offset is just an integer the prediction gives back this integer
  ret$predictOffset <- predictOffset
  if(is.null(offsetVec)) ret$predictOffset <- function(time) ret$offset
  
  # offsetVec is an integer if no smooth offset was calculated
  ret$offsetVec <- offsetVec
  if(is.null(offsetVec)) ret$offsetVec <- ret$offset 
      
  # save the call
  ret$call <- match.call()
  
  # save the evaluated call
  ret$callEval <- ret$call
  ret$callEval[-1] <- lapply(ret$call[-1], function(x){  
    eval(x, parent.frame(3)) # use the environment of the call to FDboost()
  })
  ret$callEval$data <- NULL # do not save data and weights in callEval to save memory
  ret$callEval$weights <- NULL
  
  # save formulas as character strings to save memory
  ret$timeformula <- paste(deparse(timeformula), collapse="")
  ret$formulaFDboost <- paste(deparse(formulaFDboost), collapse="")
  ret$formulaMboost <- paste(deparse(fm), collapse="")
  
#   ret$timeformula <- timeformula
#   ret$formulaFDboost <- formulaFDboost
#   ret$formulaMboost <- fm
  
  ret
}

