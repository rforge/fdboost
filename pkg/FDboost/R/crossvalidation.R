#' Validate model by resampling over curves
#' 
#' Validate the model fit by refitting the model $n$ times. Each time leaving out one observation.
#' 
#' @param object fitted FDboost-object
#' @param response you can specify a response vector to calculate predictions errors. 
#' Defaults fo NULL which means that the response of the fitted model is used.
#' @param type character argument for specifying the cross-validation method.
#' Currently cross-validation over curves (\code{curves}) 
#' and random cross-validation (\code{kfold}) are impelemted.
#' @param B number of folds for random cross-validation.
#' @param folds a weight matrix with number of rows equal to the number of observations. 
#' The number of columns corresponds to the number of cross-validation runs. 
#' Can be computed using function \code{cvMa} or \code{cv}.
#' @param grid the grid over which the optimal mstop is searched 
#' @param getCoefCV logical, defaults to FALSE. Should the coefficient-functions/surfaces
#'  be computed for all the models on the sampled data?
#'  @param mrdDelete Delete values that are mrdDelete percent smaller then the mean
#'  of the response. Defaults to 0 which means that only response values beeing 0
#'  are not used in the calculaiton of the MRD (= mean relative deviation)
#' @param ... further arguments passed to mclapply if parallel=TRUE, otherwise ignored
#' 
#' @details Calculates honest estimates of prediction errors as the curve
#' that should be predicted is not part of the model fit. 
#' 
#' @note Use argument \code{mc.cores = 1L} to set the numbers of cores that is used in 
#' parallel computation.
#' 
#' @export
validateFDboost <- function(object, response=NULL,
                            type = c("curves", "kfold"), B = 5, folds=NULL,
                            grid=1:mstop(object), getCoefCV=FALSE, mrdDelete=0, ...){
  
  type <- match.arg(type)
  
  nObs <- object$ydim[1] # number of curves
  Gy <- object$ydim[2] # number of time-points per curve
  
  # weights of mboost-model 
  weights <- model.weights(object)
  # matrix(weights, nrow=nObs, ncol=Gy)
  
  if(is.null(response)) response <- object$response # response as vector!
  #plot(response)
  #points(response, col=3)
  
  index <- rep(1:object$ydim[1], times=object$ydim[2])
  
  if(is.null(folds)){
    if(type=="curves"){
      # set up folds so that always one curve is left out for the estimation
      folds1 <- kronecker( rep(1, l=Gy), -diag(nObs)+1) # folds are row-wise!
      folds <- folds1*weights
      # matrix(folds[,1], nrow=nObs, ncol=Gy)
      attr(folds, "type") <- "over curves"
      B <- nObs
    }
    if(type=="kfold"){
      # set up folds so that cross-validation is done at random over all curves
      folds1 <- cv(weights=rep(1, length(weights)), type="kfold", B=B)
      folds <- folds1*weights
    } 
  }
  
  # out-of-bag-weights: i.e. the left out curve/ the left out observations
  OOBweights <- matrix(rep(weights, ncol(folds)), ncol = ncol(folds))
  OOBweights[folds > 0] <- 0 
  # matrix(OOBweights[,1], nrow=nObs, ncol=Gy)
  
  #   # matrix of measured responses
  #   resp <-  matrix(object$response, nrow=nObs, ncol=Gy)
  #   
  #   # Matrix of predictions using all terms in the model
  #   predFun <- matrix(nrow=nObs, ncol=Gy)
  #   
  #   # List of predictions using only one term at a time
  #   predFunCoef <- replicate(length(object$baselearner), predFun, simplify=FALSE)
  #   coefCV <- list()
  
  # Function to suppress the warning of missings in the response
  h <- function(w){
    if( any( grepl( "response contains missing values;", w) ) ) 
      invokeRestart( "muffleWarning" )  
  }
  
  ###### Function to fit the model
  # function working with FDboost, thus the smooth offset is recalculated in each model
  dummyfct <- function(weights, oobweights) {
    
    # create data frame for model fit and use the weights vector for the CV
    dathelp <- object$data
    dathelp[[attr(object$yind, "nameyind")]] <- object$yind
    dathelp[[object$yname]] <- matrix(object$response, ncol=object$ydim[2])
    
    call <- object$callEval
    call$data <- dathelp
    
    # use weights of training data
    call$weights <- weights
    
    # as weights of original model fit are used they should not be transformed again
    call$numInt <- "equal" 
    
    #call$control <- boost_control(risk="oobag")
    
    #     # <FIXME> is there a more elegant way to make sure that formula, timeformula are known in
    #     # the environment of the model fit?
    #     call$formula <- get("formulaFDboost", environment(object$predictOffset))
    #     call$timeformula <- get("timeformula", environment(object$predictOffset))
    
    mod <- withCallingHandlers(suppressMessages(eval(call)), warning = h) # suppress the warning of missing responses    
    #test <- FDboost(object$formulaFDboost, timeformula=object$timeformula, data=dathelp)
    
    mod[max(grid)]
    
    # oob risk
    ## risk <- mod$risk()[grid] # this risk is out-of-bag if the argument in control is set to "oob"
    # calculate the risk out of bag
    
    # get risk function of the family
    riskfct <- get("family", environment(mod$update))@risk
    
    # oobweights über riskfct() wie aus mboost
    risk <- sapply(grid, function(g){riskfct( mod$response, mod[g]$fitted(), w=oobweights)})
    
    #     # do not include oobweights in ()^2! gives same results as calculation im mboost, for family=Gaussian
    #     risk <- sapply(grid, function(g){
    #       sum( oobweights*(((mod$response -  mod[g]$fitted())))^2 ) } )
    
    # mod[grid[which.min(risk)]]
    
    # prediction for all observations, not only oob! 
    # -> oob-predictions have to be merged out of predGrid
    predGrid <- predict(mod, aggregate="cumsum", unlist=FALSE)
    predGrid <- sapply(predGrid, function(x) as.vector(x) )[,grid] # save vectors of predictions in matrix
    
    # estimations of the coefficients
    # use the coef.FDboost() function to calculate the coefficients 
    
    # Compute coefficients for all mstops in grid -> takes some time!
    coefCV <- list()
    
    if(getCoefCV){
      # <ToDo> deal with j=1 in more intelligent way.
      for(j in 1:length(object$baselearner)){
        coefCV[[j]] <- coef(mod, which=j, n1 = 20, n2 = 20, n3 = 15, n4 = 10)$smterms[[1]]
        coefCV[[j]]$value <- lapply(grid, function(g){
          ret <- coef(mod[g], which=j, n1 = 50, n2 = 20, n3 = 15, n4 = 10)$smterms[[1]]$value
          if(j>1) attr(ret, "offset") <- NULL # as offset is the same within one model
          return(ret)
        })
        names(coefCV[[j]]$value) <- grid
      }
    }
    
    return(list(risk=risk, predGrid=predGrid, coefCV=coefCV, mod=mod))  # return model to be able to compute coefficients
  }
  
  #   ###### Function to fit the model
  #   # This code is ok, if offset=0, i.e. offset=NULL in mboost,
  #   # otherwise the offset of the whole model is used in the optimization!
  #   fitfct <- object$update  
  #   dummyfct <- function(weights, oobweights) {
  #     
  #     # fit the model as mboost-model!
  #     mod <- withCallingHandlers(fitfct(weights = weights, oobweights = oobweights), 
  #                                warning = h) # suppress the warning of missing responses
  #     mod[max(grid)]
  #     
  #     # oob risk
  #     risk <- mod$risk()[grid] 
  #     
  #     # prediction for all observations, not only oob! 
  #     # -> oob-predictions have to be merged out of predGrid
  #     predGrid <- predict(mod, aggregate="cumsum")[,grid] 
  #     
  #     # estimations of the coefficients
  #     # use the coef.FDboost() function to calculate the coefficients 
  #     mod$yind <- object$yind
  #     mod$predictOffset <- function(time) mod$offset
  #     mod$offsetVec <- mod$offset
  #     class(mod)
  #     class(mod) <- c("FDboost", class(mod))
  #     
  #     # Compute coefficients all mstops in grid -> takes some time!
  #     coefCV <- list()
  #     
  #     if(getCoefCV){
  #       # <ToDo> deal with j=1 in more intelligent way.
  #       for(j in 1:length(object$baselearner)){
  #         coefCV[[j]] <- coef(mod, which=j, n1 = 20, n2 = 20, n3 = 15, n4 = 10)$smterms[[1]]
  #         coefCV[[j]]$value <- lapply(grid, function(g){
  #           ret <- coef(mod[g], which=j, n1 = 50, n2 = 20, n3 = 15, n4 = 10)$smterms[[1]]$value
  #           if(j>1) attr(ret, "offset") <- NULL # as offset is the same within one model
  #           return(ret)
  #         })
  #         names(coefCV[[j]]$value) <- grid
  #       } 
  #       #     # plot example coefficients one for a small number on the grid, the other for a large one
  #       #     matplot(t(coefCV[[2]]$value[[10]]), type="l", lty=1)
  #       #     matplot(t(coefCV[[2]]$value[[3]]), type="l", add=TRUE, lty=3)  
  #       #     matplot(t(coefCV[[13]]$value[[10]]), type="l", lty=1)
  #       #     matplot(t(coefCV[[13]]$value[[3]]), type="l", add=TRUE, lty=3) 
  #     } 
  #     rm(mod)
  #     
  #     return(list(risk=risk, predGrid=predGrid, coefCV=coefCV)) #mod=mod
  #   }
  #       
  
  ### computation of models on partitions of data
  if(Sys.info()["sysname"]=="Linux"){
    modRisk <- mclapply(1:ncol(folds),
                        function(i) dummyfct(weights = folds[, i],
                                             oobweights = OOBweights[, i]), ...)
  }else{
    modRisk <- mclapply(1:ncol(folds),
                        function(i) dummyfct(weights = folds[, i], 
                                             oobweights = OOBweights[, i]), mc.cores=1)
  }
  
  # str(modRisk, max.level=2)
  # str(modRisk, max.level=5)
  
  # function call without parallelization
  if(FALSE){
    modRisk <- list()
    for(i in 1:ncol(folds)){
      modRisk[[i]] <- dummyfct(weights = folds[, i], oobweights = OOBweights[, i]) 
    }
  }
  
  # check whether model fit worked in all iterations
  modFitted <- sapply(modRisk, function(x) class(x)=="list")
  if(any(!modFitted)){
    
    # stop() or warning()?
    if(sum(!modFitted) > sum(modFitted)) warning("More than half of the models could not be fitted.")
    
    warning("Model fit did not work in fold ", paste(which(!modFitted), collapse=", "))
    modRisk <- modRisk[modFitted]
    OOBweights <- OOBweights[,modFitted]  
    folds1 <- folds1[,modFitted] 
    folds <- folds[,modFitted]
  } 
  
  # matrix(model.weights(modRisk[[1]]$mod), ncol=Gy) 
  
  ####### restructure the results
  
  ## get the out-of-bag risk
  oobrisk <- t(sapply(modRisk, function(x) x$risk))
  ### <SB> transformation like in mboost - "mean(oobrisk)"
  oobrisk <- oobrisk / colSums(OOBweights) 
  colnames(oobrisk) <- grid
  rownames(oobrisk) <- which(modFitted)
  
  ## check for curves with extreme risk-values at the global median
  riskOptimal <- oobrisk[ , which.min(apply(oobrisk, 2, median))]
  bound <- median(riskOptimal) + 1.5*(quantile(riskOptimal, 0.75) - quantile(riskOptimal, 0.25))
  
  if(any(riskOptimal>bound)){
    message("Curves with high values in oobrisk:")
    message(paste("Curve ", which(riskOptimal>bound), ": " ,
                  round(riskOptimal[which(riskOptimal>bound)], 2), collapse=",  ", sep="" )  )  
  }

  
  ## predict response for all mstops in grid out of bag
  # predictions for each response are in a vector!
  oobpreds0 <- lapply(modRisk, function(x) x$predGrid)
  oobpreds <- matrix(nrow=nrow(oobpreds0[[1]]) , ncol=ncol(oobpreds0[[1]]))
  for(j in 1:length(oobpreds0)){
    oobpreds[folds1[,j]==0] <- oobpreds0[[j]][folds1[,j]==0] 
  }
  colnames(oobpreds) <- grid
  rm(oobpreds0)
  
  #   # use this predictions to calculate the MSE for each mstop in grid
  #   # results are slightly different due to handlig of weights and missings
  #   round(colMeans((oobpreds - object$response)^2, na.rm=TRUE)); round(colMeans(oobrisk))
  #   plot(colMeans((oobpreds - object$response)^2, na.rm=TRUE), colMeans(oobrisk))
  
  # Almost equal: look at 2nd mstop in grid
  # mean(oobrisk[, 2]); mean((oobpreds[,2] - response)^2, na.rm=TRUE)
  
  ############################
  # calculate mean squared error and root mean squared error for each mstop!
  mse <- colMeans((oobpreds - response)^2, na.rm=TRUE)
  
  # Calculate MSE for each curve
  # apply over the columns of oobpreds, i.e. over grid
  # tapply over the individual curves denoted by index
  mseCurves <- apply( ((oobpreds - response)^2), 2, function(xvec){
    tapply(xvec, index, function(x) mean(x, na.rm=TRUE))
  } )
  
  # Check the calculation of mseCurves and mse
  # nPerCurve <- tapply(response, index, function(x) sum(!is.na(x)))
  # weighted.mean(mseCurves[,1], nPerCurve)
  
  # compare: colMeans(oobrisk); mse
  # they are slightly different 
  
  rmse <- sqrt(mse)
  rmseCurves <- sqrt(mseCurves)
  
  # mean relative deviation
  # do not use values that are 0
  response1 <- response
  sum(round(response1, 1) == 0, na.rm=TRUE)
  response1[ round(response1, 1) == 0 ] <- NA
  
  if(mrdDelete>0){
    # do not use points that are less than 20% of the overall mean without 0 observations
    respNotused <- 0.2*mean(abs(response1), na.rm=TRUE) > abs(response1)
    #if(sum(respNotused, na.rm=TRUE)>0) print(sum(respNotused, na.rm=TRUE))
    response1[ respNotused ] <- NA
  }
  
  mrd <- colMeans( abs((oobpreds - response1) / response1), na.rm=TRUE ) 
  
  if(mrdDelete>0){
    attr(mrd, "notUsed") <- sum(respNotused, na.rm=TRUE)
  }
  
  mrdCurves <- apply( (abs((oobpreds - response1) / response1)), 2, function(xvec){
    tapply(xvec, index, function(x) mean(x, na.rm=TRUE))
  } )
  
  # Check the calculation of mrdCurves and mrd
  # nPerCurve <- tapply(response, index, function(x) sum(!is.na(x)))
  # weighted.mean(mrdCurves[,1], nPerCurve)
  
  # get the coefficient estimates
  if(getCoefCV){
    coefCV <- lapply(modRisk, function(x) x$coefCV)
  }else{
    coefCV <- NULL
  } 
  
  # calculate coefficients for the median mstop
  coefCV <- list()
  predCV <- list()
  
  if(getCoefCV){
    
    optimalMstop <- grid[which.min(apply(oobrisk, 2, median))]

    # estiamtes of coefficients
    for(l in 1:length(modRisk[[1]]$mod$baselearner)){
      # estimate the coefficients for the model of the first fold
      coefCV[[l]] <- coef(modRisk[[1]]$mod[optimalMstop], 
                          which=l, n1 = 20, n2 = 20, n3 = 15, n4 = 10)$smterms[[1]]
      coefCV[[l]]$value <- lapply(1:length(modRisk), function(g){
        ret <- coef(modRisk[[g]]$mod[optimalMstop], 
                    which=l, n1 = 50, n2 = 20, n3 = 15, n4 = 10)$smterms[[1]]$value
        return(ret)
      })
    }
    
    niceNames <- c("offset", lapply(coefCV, function(x) x$main))
    
    # predictions of terms based on the coefficients for each model
    for(l in 1:(length(modRisk[[1]]$mod$baselearner)+1)){
      predCV[[l]] <- t(sapply(1:length(modRisk), function(g){
        if(l==1){
          ret <- attr(predict(modRisk[[g]]$mod[optimalMstop], which=l), "offset")          
          if( !is.null(dim(ret)) ){
            ret <- ret[1,]
          }else{
            ret <- rep(ret, modRisk[[1]]$mod$ydim[2])
          }
        }else{
          ret <- predict(modRisk[[g]]$mod[optimalMstop], which=l-1)
          if( !is.null(dim(ret)) ){
            ret <- ret[g,]
          }else{
            ret <- rep(0, modRisk[[1]]$mod$ydim[2])
          }
        }
        return(ret)
      }))
      names(predCV)[l] <- niceNames[l]
#       matplot(modRisk[[1]]$mod$yind, t(predCV[[l]]), type="l", 
#               main=names(predCV)[l], xlab=attr(modRisk[[1]]$mod$yind, "nameyind"), ylab="coef")
    }
    
  }
  
  rm(modRisk)
  
  ret <- list(response=response, yind=object$yind,
              folds=folds, grid=grid,
              coefCV=coefCV,
              predCV=predCV,
              oobpreds=oobpreds, 
              oobrisk=oobrisk, 
              oobriskMean=colMeans(oobrisk),
              mseCurves=mseCurves, rmseCurves=rmseCurves, mrdCurves=mrdCurves, 
              mse=mse, rmse=rmse, mrd=mrd)
  
  class(ret) <- "validateFDboost"
  
  return(ret) 
}



#' @rdname plot.validateFDboost
#' @method mstop validateFDboost
#' @export
#' 
# Function to extract the optimal stopping iteration
mstop.validateFDboost <- function(object, ...){
  
  dots <- list(...)
  
  if(is.null(dots$risk)){
    risk <- "median"
  }else{
    risk <- dots$risk
  }
  
  if(risk=="median"){
    riskMedian <- apply(object$oobrisk, 2, median) 
    mstop <- object$grid[which.min(riskMedian)]
    attr(mstop, "risk") <- "minimize median risk"
  }else{
    riskMean <- colMeans(object$oobrisk)
    mstop <- object$grid[which.min(riskMean)]
    attr(mstop, "risk") <- "minimize mean risk"
  }  
  return(mstop) 
}


#' Methods for objects of class validateFDboost
#' 
#' Methods for objects that are fitted to determine the optimal mstop and the 
#' prediciton error of a model fitted by FDboost.
#' 
#' @param x object of class validateFDboost
#' @param object object of class validateFDboost
#' @param risk which risk is minimized to obtain the optimal stopping iteration?
#' defaults to the median
#' @param modObject if the original model object of class \code{FDboost} is given 
#' predicted values of the whole model can be compared to the predictions of the cross-validated models
#' @param predictNA should missing values in the response be predicted? Defaults to FALSE.
#' @param names.arg names of the observed curves
#' @param ask par(ask=ask)
#' @param commonRange, plot predicted coefficients on a common range, defaults to TRUE
#' @param showNumbers show number of curve in plot of predicted coefficients, defaults to FALSE
#' @param terms, logical, defaults to TRUE plot the added terms (default) or the coefficients?
#' @param probs vector of qunatiles to be used in the plotting of 2-dimensional coefficients surfaces,
#' defaults to \code{probs=c(0.25, 0.5, 0.75)}
#' @param ... additional arguments passed to callies.
#' 
#' @details \code{plot.validateFDboost} plots cross-validated risk, RMSE, MRD, measured and predicted values 
#' and residuals as determined by validateFDboost.
#' \code{mstop.validateFDboost} extracts the optimal mstop by minimizing the median risk.
#' 
#' @aliases mstop.validateFDboost
#' 
#' @method plot validateFDboost
#' 
#' @export
plot.validateFDboost <- function(x, risk=c("median","mean"),
                                 modObject=NULL, predictNA=FALSE, 
                                 names.arg=NULL, ask=TRUE, ...){
  
  # get the optimal mstop
  risk <- match.arg(risk)
  mopt <- mstop(x, risk=risk)
  # get the position in the grid of mopt
  mpos <- which(x$grid==mopt)
  
  par(ask=ask)
  
  # Plot the cross validated risk
  ylim <- c(min(x$oobrisk), quantile(x$oobrisk, 0.97))
  matplot(colnames(x$oobrisk), t(x$oobrisk), type="l", col="grey", 
          ylim=ylim, xlab="Number of boosting iterations", ylab="Squared Error",
          main=attr(x$folds, "type"), 
          sub=attr(mopt, "risk"))
  
  riskMean <- colMeans(x$oobrisk)
  lines(colnames(x$oobrisk), riskMean, lty=1)
  mOptMean <- x$grid[which.min(riskMean)]
  lines(c(mOptMean, mOptMean), 
        c(min(c(0, ylim[1] * ifelse(ylim[1] < 0, 2, 0.5))), 
          riskMean[paste(mOptMean)]), lty = 1)
  
  riskMedian <- apply(x$oobrisk, 2, median) 
  lines(colnames(x$oobrisk), riskMedian, lty=2)
  mOptMedian <- x$grid[which.min(riskMedian)]
  lines(c(mOptMedian, mOptMedian), 
        c(min(c(0, ylim[1] * ifelse(ylim[1] < 0, 2, 0.5))), 
          riskMedian[paste(mOptMedian)]), lty = 2)
  
  legend("topright", legend=paste(c(mOptMean, mOptMedian), c("(mean)","(median)")),
         lty=c(1,2), col=c("black","black"))
  
  #mOptOverModels <- apply(x$oobrisk, 1, which.min)
  #abline(v=mOptOverModels, lty=3)
  
  # Plot RMSE and MRD for optimal mstop
  if(is.null(names.arg)){
    names.arg = 1:length(x$rmseCurves[,mpos])
  } 
  barplot(x$rmseCurves[,mpos], main="RMSE", names.arg = names.arg, las=2) # 
  abline(h=x$rmse[mpos], lty=2)
  barplot(x$mrdCurves[,mpos], main="MRD", names.arg = names.arg, las=2)
  abline(h=x$mrd[mpos], lty=2)
  
  # Plot the predictions for the optimal mstop
  response <- x$response
  pred <- x$oobpreds[,mpos]
  if(!predictNA){
    pred[is.na(response)] <- NA
  }
  
  predMat <- matrix(pred, ncol=length(x$yind))
  responseMat <- matrix(response, ncol=length(x$yind))
  
  ylim <- range(response, pred, na.rm = TRUE)
  funplot(x$yind, responseMat, lwd = 1, pch = 1, ylim = ylim,  
          ylab = "", xlab = attr(x$yind, "nameyind"), main="Measured and Predicted Values", ...)
  funplot(x$yind, predMat, lwd = 1.5, pch = 2, add = TRUE, ...)
  
  # Plot residuals for the optimal mstop
  funplot(x$yind, responseMat-predMat, ylab = "", xlab = attr(x$yind, "nameyind"), 
          main="Residuals", ...)
  abline(h=0, lty=2, col="grey")
  
  # Plot coefficients
  
  #   # example: plot coeficients of 5th effect for folds 1-4 each for the optimal mstop
  #   plot(modObject, which=5, n1 = 20, n2 = 20, n3 = 15, n4 = 10, levels=seq(-.4, 1.4, by=0.4))
  #   contour(x$coefCV[[1]][[5]]$x,
  #           x$coefCV[[1]][[5]]$y,
  #           t(x$coefCV[[1]][[5]]$value[[mpos]]), lty=2, add=TRUE, col=2, levels=seq(-.4, 1.4, by=0.4))
  #   contour(x$coefCV[[2]][[5]]$x,
  #           x$coefCV[[2]][[5]]$y,
  #           t(x$coefCV[[2]][[5]]$value[[mpos]]), lty=2, add=TRUE, col=2, levels=seq(-.4, 1.4, by=0.4))
  #   
  #   image(x$coefCV[[1]][[5]]$value[[mpos]], type="l", lty=1)
  
  #   matplot(x$coefCV[[1]][[5]]$value[[mpos]], type="l", lty=1)
  #   matplot(x$coefCV[[2]][[5]]$value[[mpos]], type="l", lty=2, add=TRUE)
  #   matplot(x$coefCV[[3]][[5]]$value[[mpos]], type="l", lty=3, add=TRUE)
  #   matplot(x$coefCV[[4]][[5]]$value[[mpos]], type="l", lty=4, add=TRUE)
  
  # old possibility to plot the coefficients
  #   if(!is.null(modObject)){  
  #     for(j in 1:length(x$predFunCoef)){      
  #       predObject <- predict(modObject, which=j)
  #       if(j==1) predObject <- predObject + attr(predObject, "offset")
  #       funplot(x$yind, x$predFunCoef[[j]], lty=3, xlab="t", ylab="",
  #               ylim=range(x$predFunCoef[[j]], predObject, na.rm=TRUE))
  #       funplot(x$yind, predObject, pch=1, add=TRUE)
  #     } 
  #   }
  par(ask=FALSE)  
}


#' @rdname plot.validateFDboost
#' @export
#' 
plotPredCoef <- function(x, commonRange=TRUE, showNumbers=FALSE, ask=TRUE, 
                         terms=TRUE,         
                         probs=c(0.25, 0.5, 0.75), # quantiles of variables to use for plotting
                         ...){
  
  stopifnot(any(class(x)=="validateFDboost"))
  
  par(ask=ask)
  
  ylim <- NULL
  
  if(terms){
    if(commonRange){
      ylim <- range(x$predCV[-1]) # exclude offset in calculation of common range
    }
    
    for(l in 1:length(x$predCV)){
      
      if(l==1){ # do not use common range for offset
        matplot(x$yind, t(x$predCV[[l]]), type="l", 
                main=names(x$predCV)[l], xlab=attr(x$yind, "nameyind"), ylab="coef", ylim=NULL, ...)
      }else{
        matplot(x$yind, t(x$predCV[[l]]), type="l", 
                main=names(x$predCV)[l], xlab=attr(x$yind, "nameyind"), ylab="coef", ylim=ylim, ...)
      }
      
      if(showNumbers){
        matplot(x$yind, t(x$predCV[[l]]), add=TRUE )
      }
    }
  }else{ # plot coefficients
    
    if(commonRange){
      ylim <- range(lapply(x$coefCV[-1], function(x) range(x$value))) # exclude offset in calculation of common range
    }
    
    for(l in 2:length(x$coefCV)){ # loop over effects
      
      # coef() of a certain term
      temp <- x$coefCV[l][[1]]
      
      if(temp$dim==2){
        quantx <- quantile(temp$x, probs=probs, type=1)
        quanty <- quantile(temp$y, probs=probs, type=1)
        
        for(j in 1:length(probs)){ 
          
          # impute matrix of 0 if effect was never chosen
          temp$value[sapply(temp$value, function(x) is.null(dim(x)))] <- list(matrix(0, ncol=20, nrow=20))
          
          myCol <- sapply(temp$value, function(x) x[, quanty[j]==temp$y]) # first column
          matplot(temp$x, myCol, type="l", xlab=temp$xlab, ylim=ylim,
                  main=paste(temp$main, " at ", probs[j]*100, "% of ", temp$ylab, sep=""), ylab="coef")
          
          if(showNumbers){
            matplot(temp$x, myCol, add=TRUE )
          }
        }
        
        for(j in 1:length(probs)){  
          myRow <- sapply(temp$value, function(x) x[quantx[j]==temp$x, ]) # first column
          matplot(temp$x, myRow, type="l", xlab=temp$ylab, ylim=ylim,
                  main=paste(temp$main, " at ", probs[j]*100, "% of ", temp$xlab, sep=""), ylab="coef")
          
          if(showNumbers){
            matplot(temp$x, myRow, add=TRUE )
          }
        }
      }
      
    } # end loop over effects
  }
  par(ask=FALSE)
}




#' Function to set up folds for a response-matrix  
#' 
#' Function to set up folds for a response-matrix that are used in cross-validation
#' 
#' @param ydim dimensions of reponse-matrix
#' @param Y response-matrix
#' @param weights a numeric vector of weights for the model to be cross-validated.
#' @param type character argument for specifying the cross-validation 
#' method. Currently (stratified) bootstrap, k-fold cross-validation 
#' and subsampling are implemented.
#' @param B number of folds, per default 25 for \code{bootstrap} and
#' \code{subsampling} and 10 for \code{kfold}.
#' @param prob percentage of observations to be included in the learning samples 
#' for subsampling.
#' @param strata a factor of the same length as \code{weights} for stratification.
#' @param id id-variable to sample upon ids instead of sampling single observations. 
#' (Only interesting in the case that several observations per individual were made.)
#' @seealso \code{\link{cvrisk}} to perform cross-validation.
#' @return \code{cvMa} retruns a matrix of weights to be used in \code{cvrisk}.
#' @keywords models, regression
#' 
#' @details
#'   It is sufficient to specify one of \code{ydim} and \code{Y}. 
#'   
#'   The function \code{cvMa} can be used to build an appropriate 
#'   weight matrix to be used with \code{cvrisk}. 
#'   In \code{cvMa()} whole trajectories are sampled.
#'   If \code{strata} is defined 
#'   sampling is performed in each stratum separately thus preserving 
#'   the distribution of the \code{strata} variable in each fold. 
#'   If \code{id} is defined sampling is performed on the level of \code{id} 
#'   thus sampling whole individuals.
#' 
#' @examples
#' Ytest <- matrix(rnorm(15), ncol=3) # 5 trajectories, each with 3 observations 
#' cvMa(Y=Ytest, type="bootstrap", B=4) # 4 columns with bootstrap-weights
#' cvMa(ydim=c(5,3), type="bootstrap", B=4) # the same
#' @export
# Modified version of funciton cv() taken from package mboost
# add option id to sample on the level of id if there are reapeated measures
cvMa <- function(ydim=NULL, Y=NULL, weights=NULL, 
                 type = c("bootstrap", "kfold", "subsampling"), 
                 B = ifelse(type == "kfold", 10, 25), prob = 0.5, strata = NULL, id=NULL){
  
  ####################### functions taken from mboost
  # Weights for cross-validation and boosting
  # functions for cross-validation and bootstrap weights, 
  # taken from the package mboost
  cvkfold <- function (n, k){
    if (k > n/2) 
      stop("k > n/2")
    fl <- floor(n/k)
    folds <- c(rep(c(rep(0, fl), rep(1, n)), k - 1), rep(0, n * k - (k - 1) * (fl + n)))
    matrix(folds, nrow = n)[sample(1:n), , drop = FALSE]
  }
  cvboot <- function(n, B, weights) rmultinom(B, n, weights/sum(weights))
  
  cvsub <- function(n, prob, B) {
    k <- floor(n * prob)
    indx <- rep(c(0, 1), c(n - k, k))
    replicate(B, sample(indx))[sample(1:n),, drop = FALSE]
  }
  #######################
  
  if(is.null(ydim) & is.null(Y)) stop("You have to specify Y or ydim")
  
  ncolY <- if(is.null(ydim)) ncol(Y) else ydim[2]
  nrowY <- if(is.null(ydim)) nrow(Y) else ydim[1]  
  if(is.null(weights)) weights <- rep(1, nrowY)
  
  type <- match.arg(type)
  n <- length(weights)
  
  if (nrowY!=n) stop("Weights and Y/ydim have to match")
  if (is.null(strata)) strata <- gl(1, n)
  if (!is.factor(strata)) stop(sQuote("strata"), " must be a factor")
  
  
  # sample individual trajectories
  if(is.null(id)){
    folds <- matrix(0, nrow = n, ncol = B)
    for (s in levels(strata)) {
      indx <- which(strata == s)
      folds[indx, ] <- switch(type, 
                              bootstrap = cvboot(length(indx), B = B, weights[indx]), 
                              kfold = cvkfold(length(indx),k = B) * weights[indx], 
                              subsampling = cvsub(length(indx), prob = prob, 
                                                          B = B)*weights[indx])
    } 
    # sample on level of id
  }else{
    folds <- matrix(0, nrow = length(unique(id)), ncol = B)
    for (s in levels(strata)) {
      indx <- which(strata == s)
      ids <- id[which(strata == s)] # id in strata
      ids1 <- unique(ids)
      weightids <- as.numeric(table(ids))
      indx <- 1:length(ids1)
      folds[indx, ] <- switch(type, 
                              bootstrap = cvboot(length(indx), B = B, weightids), 
                              kfold = cvkfold(length(indx), k = B) * weightids, 
                              subsampling = cvsub(length(indx), prob = prob, 
                                                          B = B) * weightids)
    }
    # expand folds over ids
    if(length(unique(weightids))==1){ # equal design
      folds <- apply(folds, 2, function(x) rep(x, each=weightids[1]))
    }else{
      foldsLong <- c()
      for(i in 1:length(unique(id))){
        foldsLong <- rbind(foldsLong, apply(folds[i,,drop=FALSE], 2, 
                                            function(x) rep(x, each=weightids[i])))
      }
      folds <- foldsLong
    }        
  }
  
  # expand folds over the functional measures of the response
  # foldsMa <- apply(folds, 2, function(x) rep(x, each=ncolY))
  foldsMa <- apply(folds, 2, function(x) rep(x, times=ncolY))
  attr(foldsMa, "type") <- paste(B, "-fold ", type, sep = "")
  
  return(foldsMa)
}


