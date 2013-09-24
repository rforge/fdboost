
#' Function to truncate time in functional data 
#' 
#' @param funVar names of functional variables that should be truncated
#' @param time name of time variable
#' @param newtime new time vector that should be used. Must be part of the old timeline.
#' @param data list containig all the data
#' @note All variables that are not part if \code{funVar}, or \code{time}
#' are simply copied into the new data list
#' @return A list with the data containing all variables of the original dataset
#' with the variables of \code{funVar} truncated according to \code{newtime}.
#' @export 
truncateTime <- function(funVar, time, newtime, data){
  
  stopifnot(c(funVar, time) %in% names(data))
  stopifnot(newtime %in% data[[time]])
  
  ret <- data
  ret[[time]] <- newtime
  
  for(i in 1:length(funVar)){
    ret[[funVar[i]]] <- ret[[funVar[i]]][ , data[[time]] %in% newtime] 
  }  
  rm(data)
  return(ret)  
}
# dat <- fda::growth
# dat$hgtm <- t(dat$hgtm[,1:10])
# dat$hgtf <- t(dat$hgtf[,1:10])
# datTr <- truncateTime(funVar=c("hgtm","hgtf"), time="age", newtime=1:16, data=dat)
# 
# par(mfrow=c(1,2))
# with(dat, funplot(age, hgtm, main="Original data"))
# with(datTr, funplot(age, hgtm, main="Yearly data"))
# par(mfrow=c(1,1))

#' Plot functional data with linear interpolation of missing values 
#' 
#' @param x optional, time-vector for plotting 
#' @param y matrix of functional data with functions in rows and measured times in columns
#' @param rug logical. Should rugs be plotted? Defaults to TRUE.
#' @param ... further arguments passed to \code{\link[graphics]{matplot}}.
#' 
#' @details All observations are marked by a small cross (\code{pch=3}).
#' Missing values are imputed by linear interpolation. Parts that are
#' interpolated are plotted by dotted lines, parts with nonmissing values as solid lines.
#' @export
# with(fda::growth, funplot(age, t(hgtm)))
funplot <- function(x, y, rug=TRUE, ...){
  
  ### Get further arguments passed to the matplot-functions
  dots <- list(...)
  
  getArguments <- function(x, dots=dots){
    if(any(names(dots) %in% names(x))){
      dots[names(dots) %in% names(x)]
    }else list()
  }
  
  plotWithArgs <- function(plotFun, args, myargs){        
    args <- c(myargs[!names(myargs) %in% names(args)], args)        
    do.call(plotFun, args)            
  }
  
  argsMatplot  <- getArguments(x=c(formals(graphics::matplot), par(), main="", sub=""), dots=dots)
      
  # Deal with missing values: interpolate data
  if (missing(x)) {
    if (missing(y)) 
      stop("must specify at least one of 'x' and 'y'")
    else{
      x <- seq_len(NCOL(y))
      xlabel <- "index"
    } 
  } else xlabel <- deparse(substitute(x))
  
  ylabel <- if (!missing(y)) 
    deparse(substitute(y))
  
  time <- x
  
  stopifnot(length(x) == ncol(y))

  # Checke weather there are at least two values per row for the interpolation
  atLeast2values <- apply(y, 1, function(x) sum(is.na(x)) < length(x)-1 )
  if(any(!atLeast2values)) warning(sum(!atLeast2values), " rows contain less than 2 non-missing values.")  
  
  # Linear interpolation per row
  yint <- matrix(NA, ncol=ncol(y), nrow=nrow(y))
  yint[atLeast2values, ] <- t(apply(y[atLeast2values, ], 1, function(x) approx(time, x, xout=time)$y))
   
  # Plot the observed points
  plotWithArgs(matplot, args=argsMatplot, 
               myargs=list(x=time, y=t(y), xlab=xlabel, ylab=ylabel, type="p", pch=3) )
  
  # Plot solid lines for parts of the function without missing values
  plotWithArgs(matplot, args=argsMatplot, 
               myargs=list(x=time, y=t(y), type="l", lty=1, add=TRUE) )
 
  # Plot dotted lines for parts of the function without missing values
  plotWithArgs(matplot, args=argsMatplot, 
               myargs=list(x=time, y=t(yint), type="l", lty=3, add=TRUE) )
  
  if(rug) rug(time, 0.01)
  
#   matplot(time, t(y), xlab=xlabel, ylab="", type="p", pch=3)
#   matplot(time, t(y), type="l", pch=1, lty=1, add=TRUE)
#   matplot(time, t(yint), type="l", pch=1, lty=3, add=TRUE)  
}



#####################################################################################

#' @rdname plot.FDboost
#' @export
#' 
### function to plot the observed response and the predicted values of a model
plotPredicted <- function(x, subset=1:x$ydim[1], posLegend="topleft", ...){
  
  stopifnot("FDboost" %in% class(x))
  
  response <- matrix(x$response, nrow=x$ydim[1], ncol=x$ydim[2])[subset, ] 
  pred <- fitted(x)[subset, ]
  pred[is.na(response)] <- NA
  
  ylim <- range(response, pred, na.rm = TRUE)
  
  # Observed values
  funplot(x$yind, response, lwd=1, pch=1, ylim=ylim, lty=3, 
          ylab=x$yname, xlab=attr(x$yind, "nameyind"), ...)
  funplot(x$yind, pred, lwd=1, pch=2, add=TRUE, ...)
  # predicted values
  legend(posLegend, legend=c("observed","predicted"), col=1, pch=1:2)  
}


#####################################################################################

#' @rdname plot.FDboost
#' @export
#' 
### function to plot the residuals
plotResiduals <- function(x, subset=1:x$ydim[1], posLegend="topleft", ...){
  
  stopifnot("FDboost" %in% class(x))
  
  response <- matrix(x$response, nrow=x$ydim[1], ncol=x$ydim[2])[subset, ] 
  pred <- fitted(x)[subset, ]
  pred[is.na(response)] <- NA
  
  # Observed - predicted values
  funplot(x$yind, response-pred, ylab=x$yname, xlab=attr(x$yind, "nameyind"), ...) 
}


#####################################################################################
### Goodness of fit

# function to get y, yhat and time
getYYhatTime <- function(object, breaks=object$yind){
  
  y <- matrix(object$response, nrow=object$ydim[1], ncol=object$ydim[2]) 
  time <- object$yind
  
  ### use the original time variable
  if(all(breaks==object$yind)){ 
    yhat <- matrix(object$fitted(), nrow=object$ydim[1], ncol=object$ydim[2]) 
  }else{ ### use a time variables according to breaks
    if(length(breaks)==1){ # length of equidistant time-points
      time <- seq( min(object$yind), max(object$yind), l=breaks)      
    }else{ # time-points ot be evaluated
      time <- breaks
    }
    # Interpolate observed values
    yInter <- t(apply(y, 1, function(x) approx(object$yind, x, xout=time)$y))
    # Get dataframe to predict values at time
    newdata <- list()
    for(j in 1:length(object$baselearner)){
      datVarj <- object$baselearner[[j]]$get_data()
      if(grepl("bconcurrent", names(object$baselearner)[j])){
        datVarj <- t(apply(datVarj[[1]], 1, function(x) approx(object$yind, x, xout=time)$y))
        datVarj <- list(datVarj)
      } 
      names(datVarj) <- names(object$baselearner[[j]]$get_data())
      newdata <- c(newdata, datVarj)
    }
    newdata[[attr(object$yind, "nameyind")]] <- time 
    yhatInter <- predict(object, newdata=newdata)
    
    y <- yInter
    yhat <- yhatInter
  }
  
  return(list(y=y, yhat=yhat, time=time))
}



#' Functional R-squared
#' 
#' Claculates the functional R-squared for a fitted FDboost-object
#' 
#' @param object fitted FDboost-object
#' @param overTime per default the functional R-squared is calcualted over time
#' if \code{overTime=FALSE}, the R-squared is calculated per curve
#' @param breaks an optional vector or number giving the timepoints at which the model is evaluated.
#' Can be specified as number of equidistant time-points or as vector of time-points.
#' Defaults to the index of the response in the model.
#' @param global logical. defaults to \code{FALSE}, 
#' if TRUE the global R-squared like in a normal linear model is calculated
#' @param ... currently not used
#' 
#' @note \code{breaks} cannot be changed in the case the \code{bsignal()} 
#' is used over the same domain
#' as the response! In that case you would have to rename the index of the response or that 
#' of the covariates.
#' 
#' @details \code{breaks} should be set to some grid, if there are many
#' missing values or time-points with very few observations in the dataset.
#' Otherwiese at these points of t the variance will be almost 0 
#' (or even 0 if there is only one observation at a timepoint),
#' and then the prediction by the local means \eqn{\mu(t)} is locally very good.
#' The observations are interpolated linearily if necessary.
#' 
#' Formula to calculate R-squared over time, \code{overTime=TRUE}: \cr
#' \eqn{R^2(t) = 1 - \sum_{i}( Y_i(t) - \hat{Y}_i(t))^2 /  \sum_{i}( Y_i(t) - \bar{Y}(t) )^2 } 
#' 
#' Formula to calculate R-squared over subjects, \code{overTime=FALSE}: \cr
#' \eqn{R^2_i = 1 - \int (Y_i(t) - \hat{Y}_i(t))^2 dt /  \int (Y_i(t) - \bar{Y}(t))^2 dt }
#' 
#' @references Ramsay, J., Silverman, B. (2006). Functional data analysis. 
#' Wiley Online Library. chapter 16.3
#' 
#' @return Returns a vector with the calculated R-squared and some extra information in attributes.
#' 
#' @export
funRsquared <- function(object, overTime=TRUE, breaks=object$yind, global=FALSE, ...){
  
  # Get y, yhat and time of the model fit
  temp <- getYYhatTime(object=object, breaks=breaks)
  y <- temp$y
  yhat <- temp$yhat
  time <- temp$time
  
  stopifnot(dim(y)==dim(yhat))
  stopifnot(dim(y)[2]==length(time))
  
  if(global){
    ret <- 1 - ( sum((y-yhat)^2, na.rm=TRUE)  / sum( (y-mean(y, na.rm=TRUE))^2, na.rm=TRUE) )
    attr(ret, "name") <- "global R-squared"
    return(ret)
  }
  
  # Mean function over time (matrix containing the mean in each t in the whole column)
  mut <- matrix(colMeans(y, na.rm=TRUE), nrow=nrow(y), ncol=ncol(y), byrow=TRUE)
  
  ### for each time-point t 
  if(overTime){ 
    # numerator cannot be 0
    num <- colSums((y - mut)^2, na.rm=TRUE)
    num[round(num, 2)==0] <- NA
    ret <- 1 - ( colSums((y - yhat)^2, na.rm=TRUE) / num ) # over t 
    
    # Set values of R-squared to NA if too many values are missing
    ret[apply(y, 2, function(x) sum(is.na(x))>0.5*length(x) )] <- NA
    
    attr(ret, "name") <- "R-squared over time"
    attr(ret, "time") <- time
    attr(ret, "missings") <- apply(y, 2, function(x) sum(is.na(x))/length(x) )
    
  }else{ ### for each subject i
    # numerator cannot be 0
    num <- rowSums((y - mut)^2, na.rm=TRUE)
    num[round(num, 2)==0] <- NA    
    ret <- 1 - ( rowSums((y - yhat)^2, na.rm=TRUE) /  num ) # over subjects
    
    attr(ret, "name") <- "R-squared over subjects"
    attr(ret, "missings") <- apply(y, 1, function(x) sum(is.na(x))/length(x))
  }
  return(ret)
}


# R2t <- funRsqrt(mod)
# plot(attr(R2t, "time"), R2t, ylim=c(0,1), type="b")
# points(attr(R2t, "time"), attr(R2t, "missings"), col=2, type="b")
# abline(h=c(0, 1))
# 
# R2i <- funRsqrt(mod, overTime=FALSE)
# plot(R2i, type="b", ylim=c(0,1))
# points(attr(R2i, "missings"), col=2, type="b")
# abline(h=c(0, 1))

#' Functional MSE
#' 
#' Claculates the functional MSE for a fitted FDboost-object
#' 
#' @param object fitted FDboost-object
#' @param overTime per default the functional R-squared is calcualted over time
#' if \code{overTime=FALSE}, the R-squared is calculated per curve
#' @param breaks an optional vector or number giving the timepoints at which the model is evaluated.
#' Can be specified as number of equidistant time-points or as vector of time-points.
#' Defaults to the index of the response in the model.
#' @param global logical. defaults to \code{FALSE}, 
#' if TRUE the global R-squared like in a normal linear model is calculated
#' @param relative logical. defaults to \code{FALSE}. If \code{TRUE} the MSE is standardized
#' by the global variance of the response \cr
#' \eqn{ n^{-1} \int  \sum_i (Y_i(t) - \bar{Y})^2 dt \approx  G^{-1} n^{-1} \sum_g \sum_i (Y_i(t_g) - \bar{Y})^2 } 
#' @param root take the square root of the MSE
#' @param ... currently not used
#' 
#' @note \code{breaks} cannot be changed in the case the \code{bsignal()} 
#' is used over the same domain
#' as the response! In that case you would have to rename the index of the response or that 
#' of the covariates.
#' 
#' @details 
#' Formula to calculate MSE over time, \code{overTime=TRUE}: \cr
#' \eqn{ MSE(t) = n^{-1} \sum_i (Y_i(t) - \hat{Y}_i(t))^2 } 
#' 
#' Formula to calculate MSE over subjects, \code{overTime=FALSE}: \cr
#' \eqn{ MSE_i = \int (Y_i(t) - \hat{Y}_i(t))^2 dt  \approx G^{-1} \sum_g (Y_i(t_g) - \hat{Y}_i(t_g))^2}
#' 
#' @return Returns a vector with the calculated MSE and some extra information in attributes.
#' 
#' @export
funMSE <- function(object, overTime=TRUE, breaks=object$yind, global=FALSE, 
                   relative=FALSE, root=FALSE, ...){
  
  # Get y, yhat and time of the model fit
  temp <- getYYhatTime(object=object, breaks=breaks)
  y <- temp$y
  yhat <- temp$yhat
  time <- temp$time
  
  if(global){
    ret <- mean((y-yhat)^2, na.rm=TRUE)
    attr(ret, "name") <- "global MSE"
  }else{
    ### for each time-point t 
    if(overTime){     
      ret <- colMeans((y - yhat)^2, na.rm=TRUE)   
      attr(ret, "name") <- "MSE over time"
      attr(ret, "time") <- time
      attr(ret, "missings") <- apply(y, 2, function(x) sum(is.na(x))/length(x))     
    }else{ 
      ### for each subject i
      ret <- rowMeans((y - yhat)^2, na.rm=TRUE)
      attr(ret, "name") <- "MSE over subjects"      
      attr(ret, "missings") <- apply(y, 1, function(x) sum(is.na(x))/length(x))
    }    
  }
  
  if(relative){
    variY <- mean( (y - mean(y, na.rm=TRUE))^2, na.rm=TRUE ) # global variance of y
    ret <- ret / variY
    attr(ret, "name") <- paste("relative", attr(ret, "name"))
  }
  
  if(root){
    ret <- sqrt(ret)
    attr(ret, "name") <- paste("root", attr(ret, "name"))
  }
  
  return(ret)
}


#' Functional MRD
#' 
#' Claculates the functional MRD for a fitted FDboost-object
#' 
#' @param object fitted FDboost-object
#' @param overTime per default the functional MRD is calcualted over time
#' if \code{overTime=FALSE}, the MRD is calculated per curve
#' @param breaks an optional vector or number giving the timepoints at which the model is evaluated.
#' Can be specified as number of equidistant time-points or as vector of time-points.
#' Defaults to the index of the response in the model.
#' @param global logical. defaults to \code{FALSE}, 
#' if TRUE the global MRD like in a normal linear model is calculated
#' @param ... currently not used
#' 
#' @note \code{breaks} cannot be changed in the case the \code{bsignal()} 
#' is used over the same domain
#' as the response! In that case you would have to rename the index of the response or that 
#' of the covariates.
#' 
#' @details 
#' Formula to calculate MRD over time, \code{overTime=TRUE}: \cr
#' \eqn{ MRD(t) = n^{-1} \sum_i |(Y_i(t) - \hat{Y}_i(t))^2|/|Y_i(t)| } 
#' 
#' Formula to calculate MRD over subjects, \code{overTime=FALSE}: \cr
#' \eqn{ MRD_i = \int |(Y_i(t) - \hat{Y}_i(t))^2|/|Y_i(t)| dt  \approx G^{-1} \sum_g |(Y_i(t_g) - \hat{Y}_i(t_g))^2| / |Y_i(t)|}
#' 
#' @return Returns a vector with the calculated MRD and some extra information in attributes.
#' 
#' @export
funMRD <- function(object, overTime=TRUE, breaks=object$yind, global=FALSE,  ...){
  
  # Get y, yhat and time of the model fit
  temp <- getYYhatTime(object=object, breaks=breaks)
  y <- temp$y
  yhat <- temp$yhat
  time <- temp$time
  
  # You cannot use observations that are 0, so set them to NA
  y1 <- y
  y1[ round(y1, 1) == 0 ] <- NA
  
  if(global){
    ret <- mean( abs((y1 - yhat) / y1), na.rm=TRUE  )
    attr(ret, "name") <- "global MRD"
  }else{
    ### for each time-point t 
    if(overTime){     
      ret <- colMeans( abs((y1 - yhat) / y1), na.rm=TRUE  )   
      attr(ret, "name") <- "MRD over time"
      attr(ret, "time") <- time
      attr(ret, "missings") <- apply(y, 2, function(x) sum(is.na(x))/length(x))     
    }else{ 
      ### for each subject i
      ret <- rowMeans( abs((y1 - yhat) / y1), na.rm=TRUE  )
      attr(ret, "name") <- "MRD over subjects"      
      attr(ret, "missings") <- apply(y, 1, function(x) sum(is.na(x))/length(x))
    }    
  }
  
  return(ret)
}


