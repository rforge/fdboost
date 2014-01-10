
#' @rdname validateFDboost
#' @export
# wrapper for function cv() of mboost, additional type "curves"
# add option id to sample on the level of id if there are reapeated measures
cvMa <- function(ydim, weights=NULL, 
                 type = c("bootstrap", "kfold", "subsampling", "curves"), 
                 B = ifelse(type == "kfold", 10, 25), prob = 0.5, strata = NULL, id=NULL){
  
  ncolY <- ydim[2]
  nrowY <- ydim[1]  
  if(is.null(weights)) weights <- rep(1, nrowY*ncolY)
  
  type <- match.arg(type)
  n <- length(weights)
  
  if ( (nrowY*ncolY) != n) stop("Weights and ydim do not match.")


  if(type=="curves"){
    if(!is.null(id)) warning("Sampling is done over curves not on the level of id!")
    # set up folds so that always one curve is left out for the estimation
    foldsMa <- kronecker( rep(1, l=ncolY), -diag(nrowY)+1)*weights # folds are row-wise!
    # matrix(folds[,1], nrow=nrowY, ncol=ncolY)
    B <- nrowY
  }else{
    # expand folds over the functional measures of the response
    if(is.null(id)){
      folds <- cv(weights=rep(1, nrowY), type = type, B = B, prob = prob, strata = strata)
      #foldsMa <- apply(folds, 2, function(x) rep(x, times=ncolY))*weights #the same
      foldsMa <- folds[rep(1:nrowY, times=ncolY), ]*weights
    }else{
      stopifnot(length(id)==nrowY)
      folds <- cv(weights=rep(1, length(unique(id))), type = type, B = B, prob = prob, strata = strata)
      folds <- folds[id, ] # reapeat folds according to id
      foldsMa <- folds[rep(1:nrowY, times=ncolY), ]*weights
    }
    
  }
  attr(foldsMa, "type") <- paste(B, "-fold ", type, sep = "")  
  return(foldsMa)
}


#' Cross-Validation over Curves
#' 
#' Cross-Validation over Curves. 
#' 
#' @param object fitted FDboost-object
#' @param response you can specify a response vector to calculate predictions errors. 
#' Defaults fo NULL which means that the response of the fitted model is used.
#' @param weights a numeric vector of weights for the model to be cross-validated.
#' if weights=NULL all weights are taken to be 1.
#' @param folds a weight matrix with number of rows equal to the number of observations. 
#' @param grid the grid over which the optimal mstop is searched 
#' @param getCoefCV logical, defaults to TRUE Should the coefficients and predictions
#'  be computed for all the models on the sampled data?
#' @param mrdDelete Delete values that are mrdDelete percent smaller then the mean
#'  of the response. Defaults to 0 which means that only response values beeing 0
#'  are not used in the calculaiton of the MRD (= mean relative deviation)
#' @param ydim dimensions of reponse-matrix
#' @param type character argument for specifying the cross-validation 
#' method. Currently (stratified) bootstrap, k-fold cross-validation 
#' and subsampling are implemented.
#' The argument curves implies that a cross-validation leaving out one curve at 
#' a time is performed.
#' @param B number of folds, per default 25 for \code{bootstrap} and
#' \code{subsampling} and 10 for \code{kfold}.
#' @param prob percentage of observations to be included in the learning samples 
#' for subsampling.
#' @param strata a factor of the same length as \code{weights} for stratification.
#' @param id id-variable to sample upon ids instead of sampling single observations. 
#' (Only interesting in the case that several observations per individual were made.)
#' @param ... further arguments passed to mclapply if parallel=TRUE, otherwise ignored
#' 
#' @details The funciton \code{validateFDboost} calculates honest estimates 
#' of prediction errors as the curve/observations
#' that should be predicted is not part of the model fit.
#'   The function \code{cvMa} can be used to build an appropriate 
#'   weight matrix to be used with \code{cvrisk} or \code{validateFDboost}. 
#'   In \code{cvMa} whole trajectories are sampled. The probability for each 
#'   trajectory to enter a fold is equal over all trajectories. 
#'   If \code{strata} is defined 
#'   sampling is performed in each stratum separately thus preserving 
#'   the distribution of the \code{strata} variable in each fold. 
#'   If \code{id} is defined sampling is performed on the level of \code{id} 
#'   thus sampling individuals. 
#' 
#' @note Use argument \code{mc.cores = 1L} to set the numbers of cores that is used in 
#' parallel computation.
#' 
#' @seealso \code{\link{cvrisk}} to perform cross-validation.
#' @return \code{cvMa} retruns a matrix of weights to be used in \code{cvrisk} or \code{validateFDboost}.
#' @examples
#' Ytest <- matrix(rnorm(15), ncol=3) # 5 trajectories, each with 3 observations 
#' cvMa(ydim=c(5,3), type="bootstrap", B=4) 
#' 
#' @aliases cvMa
#' 
#' @export
validateFDboost <- function(object, response=NULL, weights=model.weights(object), 
                            folds=cvMa(ydim=object$ydim, weights=model.weights(object), type="bootstrap"),
                            grid=1:mstop(object), getCoefCV=TRUE, mrdDelete=0, ...){
  
  if(is.null(folds)){
    warning("is.null(folds), per default folds=cvMa(ydim=object$ydim, weights=model.weights(object), type=\"bootstrap\")")
    folds <- cvMa(ydim=object$ydim, weights=model.weights(object), type="bootstrap")
  }
  
  type <- attr(folds, "type")
  call <- match.call()
  
  nObs <- object$ydim[1] # number of curves
  Gy <- object$ydim[2] # number of time-points per curve
  
  # if weights=NULL set all weights to 1 
  if(is.null(weights)){
    weights <- rep(1.0, nObs*Gy)
  }
  # matrix(weights, nrow=nObs, ncol=Gy)
  
  if(is.null(response)) response <- object$response # response as vector!
  #plot(response)
  #points(response, col=3)
  
  # index of observations that belong to the same trajectory
  index <- rep(1:object$ydim[1], times=object$ydim[2])
  
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
    
    # <FIXME> check whether the risk is calculated oob?
    #call$control <- boost_control(risk="oobag")
    
    mod <- withCallingHandlers(suppressMessages(eval(call)), warning = h) # suppress the warning of missing responses    
    #test <- FDboost(formula(object$formulaFDboost), timeformula=formula(object$timeformula), data=dathelp)
    
    mod[max(grid)]
    
    # oob risk
    ## risk <- mod$risk()[grid] # this risk is out-of-bag if the argument in control is set to "oob"
    # calculate the risk out of bag
    
    # get risk function of the family
    riskfct <- get("family", environment(mod$update))@risk
    
    # oobweights using riskfct() like in mboost
    risk <- sapply(grid, function(g){riskfct( mod$response, mod[g]$fitted(), w=oobweights)})
    
    #     # do not include oobweights in ()^2! gives same results as calculation im mboost, for family=Gaussian
    #     risk <- sapply(grid, function(g){
    #       sum( oobweights*(((mod$response -  mod[g]$fitted())))^2 ) } )
    
    # mod[grid[which.min(risk)]]
    
    # matrix(oobweights, nObs, Gy)
    
    relMSE <- simplify2array(mclapply(grid, function(g){
      sum(((mod$response - mod[g]$fitted())[oobweights==1])^2) /
        sum( ( (mod$response - mean(mod$response[oobweights]))[oobweights==1] )^2 )
    }, mc.cores=1) )
    
    # prediction for all observations, not only oob! 
    # -> oob-predictions have to be merged out of predGrid
    predGrid <- predict(mod, aggregate="cumsum", unlist=FALSE)
    predGrid <- sapply(predGrid, function(x) as.vector(x) )[,grid] # save vectors of predictions in matrix
    
    # predict oob and save responses that were predicted
    predOOB <- predict(mod, aggregate="cumsum", unlist=FALSE)
    keepOOB <- matrix(oobweights, nObs, Gy)
    keepRow <- apply(keepOOB, 1, function(x) any(x!=0))
    if(any(keepRow)){
      predOOB <- sapply(predOOB, function(x) as.vector(x[keepRow,]) )[,grid]
      respOOB <- as.vector(matrix(mod$response, nObs, Gy)[keepRow,] )
      attr(respOOB, "curves") <- which(keepRow)
    }else{
      predOOB <- NULL
      respOOB <- NULL
    }
    
    return(list(risk=risk, predGrid=predGrid, predOOB=predOOB, respOOB=respOOB, relMSE=relMSE, mod=mod))  
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
  #         coefCV[[j]] <- coef(mod, which=j, n1 = 50, n2 = 20, n3 = 15, n4 = 10)$smterms[[1]]
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
  
  ## check for folds with extreme risk-values at the global median
  riskOptimal <- oobrisk[ , which.min(apply(oobrisk, 2, median))]
  bound <- median(riskOptimal) + 1.5*(quantile(riskOptimal, 0.75) - quantile(riskOptimal, 0.25))
  
  # fold equals curve if type="curves"
  if(any(riskOptimal>bound)){
    message("Fold with high values in oobrisk (median is ", round(median(riskOptimal),2), "):")
    message(paste("In fold ", which(riskOptimal>bound), ": " ,
                  round(riskOptimal[which(riskOptimal>bound)], 2), collapse=",  ", sep="" )  )  
  }
  
  ## only makes sense for type="curves" with leaving-out one curve per fold!!
  if(grepl( "curves", type)){
  ## predict response for all mstops in grid out of bag
  # predictions for each response are in a vector!
  # if a curve was several times not in the training data the last prediction is taken
  oobpreds0 <- lapply(modRisk, function(x) x$predGrid)
  oobpreds <- matrix(nrow=nrow(oobpreds0[[1]]) , ncol=ncol(oobpreds0[[1]]))
  for(j in 1:length(oobpreds0)){
    oobpreds[folds[,j]==0] <- oobpreds0[[j]][folds[,j]==0] 
  }
  colnames(oobpreds) <- grid
  rm(oobpreds0)
  
  #   # use this predictions to calculate the MSE for each mstop in grid
  #   # results are slightly different due to handlig of weights and missings
  #   round(colMeans((oobpreds - object$response)^2, na.rm=TRUE)); round(colMeans(oobrisk))
  #   plot(colMeans((oobpreds - object$response)^2, na.rm=TRUE), colMeans(oobrisk))
  
  # Almost equal: look at 2nd mstop in grid
  # mean(oobrisk[, 2]); mean((oobpreds[,2] - response)^2, na.rm=TRUE)
  }else{
    oobpreds <- NULL
  }
    
  # alternative OOB-prediction: works for general folds not only oob
  predOOB <- lapply(modRisk, function(x) x$predOOB)
  predOOB <- do.call('rbind', predOOB)
  colnames(predOOB) <- grid
  respOOB <- lapply(modRisk, function(x) x$respOOB)
  respOOB <- do.call('c', respOOB)
  indexOOB <- lapply(modRisk, function(x) attr(x$respOOB, "curves"))
  indexOOB <- lapply(indexOOB, function(x) rep(x, times=Gy) )
  indexOOB <- unlist(indexOOB)
  attr(respOOB, "index") <- indexOOB
  
  ############################
  # calculate mean squared error and root mean squared error for each mstop!
#  mseOld <- colMeans((oobpreds - response)^2, na.rm=TRUE)
  mse <- colMeans((predOOB - respOOB)^2, na.rm=TRUE)
  
  # Calculate MSE for each curve
#   # apply over the columns of oobpreds, i.e. over grid
#   # tapply over the individual curves denoted by index
#   mseCurvesOld <- apply( ((oobpreds - response)^2), 2, function(xvec){
#     tapply(xvec, index, function(x) mean(x, na.rm=TRUE))
#   } )
  mseCurves <- apply( ((predOOB - respOOB)^2), 2, function(xvec){
    tapply(xvec, indexOOB, function(x) mean(x, na.rm=TRUE))
  } )
  
  # Check the calculation of mseCurves and mse
  # nPerCurve <- tapply(response, index, function(x) sum(!is.na(x)))
  # weighted.mean(mseCurves[,1], nPerCurve); mse[1]
  
  # compare: cbind(colMeans(oobrisk), mse, mseOld)
  # they are slightly different 
  
  rmse <- sqrt(mse)
  rmseCurves <- sqrt(mseCurves)
  
  ## get the relMSE
  relMSE <- t(sapply(modRisk, function(x) x$relMSE)) # no division with b_i as 1/b_i reduces in relMSE
  colnames(relMSE) <- grid
  rownames(relMSE) <- which(modFitted)
  
  ### mean relative deviation
  # do not use values that are 0
  respOOB1 <- respOOB
  sum(round(respOOB, 1) == 0, na.rm=TRUE)
  respOOB1[ round(respOOB1, 1) == 0 ] <- NA
  
  if(mrdDelete>0){
    # do not use points that are less than 20% of the overall mean without 0 observations
    respNotused <- mrdDelete*mean(abs(respOOB1), na.rm=TRUE) > abs(respOOB1)
    #if(sum(respNotused, na.rm=TRUE)>0) print(sum(respNotused, na.rm=TRUE))
    respOOB1[ respNotused ] <- NA
  }
  
  mrd <- colMeans( abs((predOOB - respOOB1) / respOOB1), na.rm=TRUE ) 
  
  if(mrdDelete>0){
    attr(mrd, "notUsed") <- sum(respNotused, na.rm=TRUE)
  }
  
  mrdCurves <- apply( abs((predOOB - respOOB1) / respOOB1), 2, function(xvec){
    tapply(xvec, indexOOB, function(x) mean(x, na.rm=TRUE))
  } )
  
  # Check the calculation of mrdCurves and mrd
  # nPerCurve <- tapply(response, index, function(x) sum(!is.na(x)))
  # weighted.mean(mrdCurves[,1], nPerCurve); mrd[1]

  ##### OLD
#   # get the coefficient estimates
#   if(getCoefCV){
#     coefCV <- lapply(modRisk, function(x) x$coefCV)
#   }else{
#     coefCV <- NULL
#   } 
  
  # calculate coefficients for the median mstop
  coefCV <- list()
  predCV <- list()
  
  if(getCoefCV){
    
    # use median of oobrisk!
    optimalMstop <- grid[which.min(apply(oobrisk, 2, median))]

    ### estimates of coefficients
    timeHelp <- seq(min(modRisk[[1]]$mod$yind), max(modRisk[[1]]$mod$yind), l=50)
    for(l in 1:length(modRisk[[1]]$mod$baselearner)){
      # estimate the coefficients for the model of the first fold
      coefCV[[l]] <- coef(modRisk[[1]]$mod[optimalMstop], 
                          which=l, n1 = 50, n2 = 20, n3 = 15, n4 = 10)$smterms[[1]]
      if(l==1){
        coefCV[[l]]$offset <- matrix(ncol=50, nrow=length(modRisk))
        coefCV[[l]]$offset[1,] <- modRisk[[1]]$mod$predictOffset(time=timeHelp)
      }
      attr(coefCV[[l]]$value, "offset") <- NULL # as offset is the same within one model

      # add estimates for the models of the other folds
      coefCV[[l]]$value <- lapply(1:length(modRisk), function(g){
        ret <- coef(modRisk[[g]]$mod[optimalMstop], 
                    which=l, n1 = 50, n2 = 20, n3 = 15, n4 = 10)$smterms[[1]]$value
        if(l==1){
          coefCV[[l]]$offset[g,] <- modRisk[[g]]$mod$predictOffset(time=timeHelp)
        }
        attr(ret, "offset") <- NULL # as offset is the same within one model
        return(ret)
      })
    }
    
    niceNames <- c("offset", lapply(coefCV, function(x) x$main))
    
    ### predictions of terms based on the coefficients for each model
    # only makes sense for type="curves" with leaving-out one curve per fold!!
    if(grepl( "curves", type)){
      for(l in 1:(length(modRisk[[1]]$mod$baselearner)+1)){
        predCV[[l]] <- t(sapply(1:length(modRisk), function(g){
          if(l==1){ # save offset of model
            ret <- attr(predict(modRisk[[g]]$mod[optimalMstop], which=l), "offset")          
            if( !is.null(dim(ret)) ){
              ret <- ret[1,]
            }else{
              ret <- rep(ret, modRisk[[1]]$mod$ydim[2])
            }
          }else{ # other effects
            ret <- predict(modRisk[[g]]$mod[optimalMstop], which=l-1) # model g
            if( !is.null(dim(ret)) ){
              ret <- ret[g,] # save g-th row
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
    
  } # end of if(getCoefCV)
  
  rm(modRisk)
  
  ret <- list(response=response, yind=object$yind,
              folds=folds, grid=grid,
              coefCV=coefCV,
              predCV=predCV,
              oobpreds=oobpreds, 
              predOOB=predOOB, respOOB=respOOB,
              oobrisk=oobrisk, 
              oobriskMean=colMeans(oobrisk),
              mseCurves=mseCurves, rmseCurves=rmseCurves, mrdCurves=mrdCurves, 
              mse=mse, rmse=rmse, mrd=mrd, relMSE=relMSE)
  
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
#' @param which a subset of base-learners to take into account for plotting
#' @param commonRange, plot predicted coefficients on a common range, defaults to TRUE
#' @param showQuantiles plot the 0.05 and the 0.95 Quantile of coefficients in 1-dim effects
#' @param showNumbers show number of curve in plot of predicted coefficients, defaults to FALSE
#' @param terms, logical, defaults to TRUE plot the added terms (default) or the coefficients?
#' @param probs vector of qunatiles to be used in the plotting of 2-dimensional coefficients surfaces,
#' defaults to \code{probs=c(0.25, 0.5, 0.75)}
#' @param ylim values for limits of y-axis
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
  if(!is.null(x$oobpreds[,mpos])){
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
  }
  
#   }else{
#     response <- x$respOOB
#     pred <- x$predOOB[,mpos]
#     if(!predictNA){
#       pred[is.na(response)] <- NA
#     }
#     indexOOB <- attr(x$respOOB, "index")
#     predMat <- matrix(pred, ncol=length(x$yind))
#     responseMat <- matrix(response, ncol=length(x$yind))
#   }


  
  # Plot coefficients
  
  #   # example: plot coeficients of 5th effect for folds 1-4 each for the optimal mstop
  #   plot(modObject, which=5, n1 = 50, n2 = 20, n3 = 15, n4 = 10, levels=seq(-.4, 1.4, by=0.4))
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
plotPredCoef <- function(x, which=NULL, 
                         commonRange=TRUE, showNumbers=FALSE, showQuantiles=TRUE,
                         ask=TRUE, 
                         terms=TRUE,         
                         probs=c(0.25, 0.5, 0.75), # quantiles of variables to use for plotting
                         ylim=NULL, ...){
  
  stopifnot(any(class(x)=="validateFDboost"))
  
  if(is.null(which)) which <- 1:length(x$coefCV)
  
  if(length(which)>1) par(ask=ask)
  
  if(terms){
    
    if( length(x$predCV)==0 ){
      warning("terms=TRUE, but predCV is empty.")
      return(NULL)
    }
    
    if(commonRange & is.null(ylim)){ 
      if(length(x$yind)>1){
        ylim <- range(x$predCV[-1]) # exclude offset in calculation of common range
      }else{
        ylim <- range(x$predCV[which])
      } 
    }
    
    for(l in which){
      
      if(l==1 && length(x$yind)>1){ # do not use common range for offset for functional response
        matplot(x$yind, t(x$predCV[[l]]), type="l", 
                main=names(x$predCV)[l], xlab=attr(x$yind, "nameyind"), ylab="coef", ylim=NULL, ...)
      }else{
        matplot(x$yind, t(x$predCV[[l]]), type="l", 
                main=names(x$predCV)[l], xlab=attr(x$yind, "nameyind"), ylab="coef", ylim=ylim, ...)
      }
      
      if(showNumbers){
        matplot(x$yind, t(x$predCV[[l]]), add=TRUE, ...)
      }
    }
  }else{ # plot coefficients
    
    if(commonRange & is.null(ylim)){
      if(length(x$yind)>1){
        ylim <- range(lapply(x$coefCV[-1], function(x) x$value)) # exclude offset in calculation of common range
      }else{
        ylim <- range(lapply(x$coefCV[which], function(x) x$value))
      } 
    }
    
    for(l in which){ # loop over effects
      
      # coef() of a certain term
      temp <- x$coefCV[l][[1]]
      
      # set the range for each effect individually
      if(FALSE) ylim <- range(temp$value)
      
      # plot the estimated offset
      if(l==1 && length(x$yind)>1){
        #temp <- x$coefCV[[1]]
        myMat <- temp$offset
        
        matplot(seq(min(x$yind), max(x$yind), length=ncol(myMat)), t(myMat), type="l", 
                xlab=attr(x$yind,"nameyind"),
                main="offset", ylab="coef", ...)
        
        if(showNumbers){
          matplot(seq(min(x$yind), max(x$yind), length=ncol(myMat)), t(myMat), add=TRUE, ...)
        }
        
      }
      
      if(temp$dim==1){
        # impute vector of 0 if effect was never chosen
        temp$value[sapply(temp$value, function(x) length(x)==1)] <- list(rep(0, 50))
        myMat <- sapply(temp$value, function(x) x) # as one matrix
        
        if(showQuantiles){
          
          matplot(temp$x, myMat, type="l", xlab=temp$xlab,
                  main=temp$main, ylab="coef", ylim=ylim, 
                  col=rgb(0.6,0.6,0.6, alpha=0.5), ...)
          
          if(showNumbers){
            matplot(temp$x, myMat, add=TRUE, col=rgb(0.6,0.6,0.6, alpha=0.5), ...)
          }
          
          lines(temp$x, rowMeans(myMat), col=1, lwd=2)
          lines(temp$x, apply(myMat, 1, quantile, 0.95), col=2, lwd=2, lty=2)
          lines(temp$x, apply(myMat, 1, quantile, 0.05), col=2, lwd=2, lty=2)
        }else{
          
          matplot(temp$x, myMat, type="l", xlab=temp$xlab,
                  main=temp$main, ylab="coef", ylim=ylim, ...)
          
          if(showNumbers){
            matplot(temp$x, myMat, add=TRUE, ...)
          }
          
        }
      }
      
      if(temp$dim==2){
        
        if(!is.factor(temp$x)){
          quantx <- quantile(temp$x, probs=probs, type=1)
        } else quantx <- temp$x
        
        quanty <- quantile(temp$y, probs=probs, type=1)
        
        if(!is.factor(temp$x)){ # temp$x is metric
          
          for(j in 1:length(quanty)){ 
            
            # impute matrix of 0 if effect was never chosen
            temp$value[sapply(temp$value, function(x) is.null(dim(x)))] <- list(matrix(0, ncol=20, nrow=20))
            
            # set lower triangular matrix to NA for historic effect
            if(grepl("bhist", temp$main)){
              for(k in 1:length(temp$value)){
                temp$value[[k]][outer(1:nrow(temp$value[[k]]), 1:ncol(temp$value[[k]]), "<=")==FALSE] <- NA
              }
            }
            
            myCol <- sapply(temp$value, function(x) x[, quanty[j]==temp$y]) # first column
            matplot(temp$x, myCol, type="l", xlab=temp$xlab, ylim=ylim,
                    main=paste(temp$main, " at ", probs[j]*100, "% of ", temp$ylab, sep=""), ylab="coef", ...)
            
            if(showNumbers){
              matplot(temp$x, myCol, add=TRUE, ...)
            }
          }
          
          for(j in 1:length(quantx)){  
            myRow <- sapply(temp$value, function(x) x[quantx[j]==temp$x, ]) # first column
            matplot(temp$x, myRow, type="l", xlab=temp$ylab, ylim=ylim,
                    main=paste(temp$main, " at ", probs[j]*100, "% of ", temp$xlab, sep=""), ylab="coef", ...)
            
            if(showNumbers){
              matplot(temp$x, myRow, add=TRUE, ...)
            }
          }
        }else{ # temp$x is factor
          for(j in 1:length(quantx)){ 
            
            # impute matrix of 0 if effect was never chosen
            temp$value[sapply(temp$value, function(x) is.null(dim(x)))] <- list(matrix(0, ncol=20, nrow=length(quantx)))
            
            if(is.null(temp$z)){
              myRow <- sapply(temp$value, function(x) x[quantx[j]==temp$x, ]) # first column
              matplot(temp$y, myRow, type="l", xlab=temp$ylab, ylim=ylim,
                      main=paste(temp$main, " at ", temp$xlab,"=" ,quantx[j], sep=""), ylab="coef", ...)
            }else{
              quantz <- temp$z
              myRow <- sapply(temp$value, function(x) x[quantx[j]==temp$x & quantz[j]==temp$z, ]) # first column
              matplot(temp$y, myRow, type="l", xlab=temp$ylab, ylim=ylim,
                      main=paste(temp$main, " at ", temp$xlab,"=" , quantx[j], ", " , 
                                 temp$zlab,"=" ,quantz[j], sep=""), ylab="coef", ...)
            }
            if(showNumbers){
              matplot(temp$y, myRow, add=TRUE )
            }
          } 
        }
        
      } # end if(temp$dim==2)
      
    } # end loop over effects
  }
  par(ask=FALSE)
}


