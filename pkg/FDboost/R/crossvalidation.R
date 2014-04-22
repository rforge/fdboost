
#' @rdname validateFDboost
#' @export
# wrapper for function cv() of mboost, additional type "curves"
# add option id to sample on the level of id if there are repeated measures
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

#' @rdname validateFDboost
#' @export
# wrapper for function cv() of mboost, additional type "curves"
# create folds for data in long format
cvLong <- function(id, weights=rep(1, l=length(id)), 
                 type = c("bootstrap", "kfold", "subsampling", "curves"), 
                 B = ifelse(type == "kfold", 10, 25), prob = 0.5, strata = NULL){
  type <- match.arg(type)
  n <- length(weights)
  
  if(type=="curves"){
    # set up folds so that always one curve is left out for the estimation
    folds <- -diag(length(unique(id)))+1
    foldsLong <- folds[id,]*weights    
    B <- length(unique(id))
  }else{
    # expand folds over the functional measures of the response
    folds <- cv(weights=rep(1, length(unique(id))), type = type, B = B, prob = prob, strata = strata)
    foldsLong <- folds[id,]*weights
  }
  attr(foldsLong, "type") <- paste(B, "-fold ", type, sep = "")  
  return(foldsLong)
}



#' Cross-Validation and Bootstrapping over Curves
#' 
#' Cross-Validation and bootstrapping over curves to compute the empirical risk for 
#' hyper-parameter selection and to compute resampled coefficients and predictions for
#' the models. 
#' 
#' @param object fitted FDboost-object
#' @param response you can specify a response vector to calculate predictions errors. 
#' Defaults to NULL which means that the response of the fitted model is used.
#' @param folds a weight matrix with number of rows equal to the number of observed trajectories.  
#' @param grid the grid over which the optimal mstop is searched 
#' @param getCoefCV logical, defaults to TRUE Should the coefficients and predictions
#' be computed for all the models on the sampled data?
#' @param riskopt how is the optimal stopping iteration determined. Defaults to the mean
#' but median is possible as well. 
#' @param mrdDelete Delete values that are mrdDelete percent smaller then the mean
#' of the response. Defaults to 0 which means that only response values being 0 
#' are not used in the calculation of the MRD (= mean relative deviation) 
#' @param refitSmoothOffset logical, should the offset be refitted in each learning sample? 
#' Defaults to TRUE. 
#' In \code{\link[mboost]{cvrisk}} the offset of the original model \code{object} is used.
#' @param ydim dimensions of response-matrix
#' @param weights a numeric vector of weights for the model to be cross-validated.
#' if weights=NULL all weights are taken to be 1.
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
#' @param ... further arguments passed to mclapply
#' 
#' @details The function \code{validateFDboost} calculates honest estimates 
#' of prediction errors as the curves/observations
#' that should be predicted are not part of the model fit.
#' The functions \code{cvMa} and \code{cvLong} can be used to build an appropriate 
#' weight matrix to be used with \code{cvrisk}. The function \code{cvMa} is appropriate for
#' responses observed over a regular grid. The function \code{cvLong} for responsese that
#' are observed over a irregular grid.  
#' Note, that the function \code{validateFDboost} expects folds that give weights
#' per trajectory and that \code{cvrisk} expectes folds that give weights for single observations.
#' In \code{cvMa} and \code{cvLong} trajectories are sampled. The probability for each 
#' trajectory to enter a fold is equal over all trajectories. 
#' If \code{strata} is defined 
#' sampling is performed in each stratum separately thus preserving 
#' the distribution of the \code{strata} variable in each fold. 
#' If \code{id} is defined sampling is performed on the level of \code{id} 
#' thus sampling individuals. 
#' 
#' @note Use argument \code{mc.cores = 1L} to set the numbers of cores that is used in 
#' parallel computation. On Windows only 1 core is possible, \code{mc.cores = 1}, which is the default.
#' 
#' @seealso \code{\link{cvrisk}} to perform cross-validation.
#' @return \code{cvMa} and \code{cvLong} return a matrix of weights to be used in \code{cvrisk}. 
#' \code{validateFDboost} returns an validateFDboost-object, which is a named list containing: 
#' \item{response}{the response}
#' \item{yind}{the observation points of the response}
#' \item{id}{the id variable of the response}
#' \item{folds}{folds that were used}
#' \item{grid}{grid of possible numbers of boosting iterations}
#' \item{coefCV}{if \code{getCoefCV} is TRUE the estimated coefficient functions in the folds}
#' \item{predCV}{if \code{getCoefCV} is TRUE the out-of-bag predicted values of the response}
#' \item{oobpreds}{if the type of folds is curves the out-of-bag predictions for each trajectory}
#' \item{oobrisk}{the out-of-bag risk}
#' \item{oobriskMean}{the out-of-bag risk at the minimal mean risk}
#' \item{oobmse}{the out-of-bag mean squared error (MSE)}
#' \item{oobrelMSE}{the out-of-bag relative mean squared error (relMSE)}
#' \item{oobmrd}{the out-of-bag mean relative deviation (MRD)}
#' \item{oobrisk0}{the out-of-bag risk without consideration of integration weights}
#' \item{oobmse0}{the out-of-bag mean squared error (MSE) without consideration of integration weights}
#' 
#' @examples
#' Ytest <- matrix(rnorm(15), ncol=3) # 5 trajectories, each with 3 observations 
#' cvMa(ydim=c(5,3), type="bootstrap", B=4)  
#' 
#' data("viscosity", package = "FDboost") 
#' ## set time-interval that should be modeled
#' interval <- "101"
#' 
#' ## model time until "interval" and take log() of viscosity
#' end <- which(viscosity$timeAll==as.numeric(interval))
#' viscosity$vis <- log(viscosity$visAll[,1:end])
#' viscosity$time <- viscosity$timeAll[1:end]
#' # with(viscosity, funplot(time, vis, pch=16, cex=0.2))
#' 
#' ## fit median regression model with 200 boosting iterations,
#' ## step-length 0.2 and
#' ## smooth time-specific offset
#' mod <- FDboost(vis ~ 1 + bols(T_C) + bols(T_A),
#'                timeformula=~bbs(time, lambda=100),
#'                numInt="Riemann", family=QuantReg(),
#'                offset=NULL, offset_control = o_control(k_min = 9),
#'                data=viscosity, control=boost_control(mstop = 100, nu = 0.4))
#' 
#' \dontrun{
#' ## for the example B is set to a small value so that bootstrap is faster               
#' val1 <- validateFDboost(mod, folds=cv(rep(1, 64), B=3) )
#' # plot(val1)
#' mstop(val1)
#' }
#' 
#' @aliases cvMa cvLong
#' 
#' @export
validateFDboost <- function(object, response=NULL,  
                            #folds=cvMa(ydim=object$ydim, weights=model.weights(object), type="bootstrap"),
                            folds=cv(rep(1, object$ydim[1]), type="bootstrap"),
                            grid=1:mstop(object), getCoefCV=TRUE, riskopt=c("mean","median"), 
                            mrdDelete=0, refitSmoothOffset=TRUE, ...){
  
#   if(is.null(folds)){
#     warning("is.null(folds), per default folds=cvMa(ydim=object$ydim, weights=model.weights(object), type=\"bootstrap\")")
#     folds <- cvMa(ydim=object$ydim, weights=model.weights(object), type="bootstrap")
#   }
  
  type <- attr(folds, "type")
  if(is.null(type)) type <- "unknown"
  call <- match.call()
  riskopt <- match.arg(riskopt)
  
  if(is.null(object$id)){
    nObs <- object$ydim[1] # number of curves
    Gy <- object$ydim[2] # number of time-points per curve
  }else{
    nObs <- length(unique(object$id)) # number of curves
    Gy <- NULL # number of time-points per curve
  }
  
  if(is.null(response)) response <- object$response # response as vector!
  #plot(response)
  #points(response, col=3)
  
  # id of observations that belong to the same trajectory
  if(is.null(object$id)){
    id <- rep(1:object$ydim[1], times=object$ydim[2])
  }else{
    id <- object$id
  }
  
  # save integration weights of original model
  intWeights <- model.weights(object)
  
  # out-of-bag-weights: i.e. the left out curve/ the left out observations
  OOBweights <- matrix(1, ncol = ncol(folds), nrow=nrow(folds))
  OOBweights[folds > 0] <- 0 
  
  #   # matrix of measured responses for regular observations
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
  
  ### get yind in long format
  yindLong <- object$yind
  if(is.null(object$id)){
    yindLong <- rep(object$yind, each=object$ydim[1])
  }
  ### compute ("length of each trajectory")^-1 in the response
  if(length(object$yind)>1){
    lengthTi1 <- 1/tapply(yindLong[!is.na(response)], id[!is.na(response)], function(x) max(x) - min(x))
    if(any(is.infinite(lengthTi1))) lengthTi1[is.infinite(lengthTi1)] <- max(lengthTi1[is.infinite(lengthTi1)])
  }else{
    lengthTi1 <- rep(1, l=length(response))
  }
  
  ###### Function to fit the model
  # function working with FDboost, thus the smooth offset is recalculated in each model
  dummyfct <- function(weights, oobweights) {
    
    # create data frame for model fit and use the weights vector for the CV
    dathelp <- object$data
    dathelp[[attr(object$yind, "nameyind")]] <- object$yind
    
    if(is.null(object$id)){
      dathelp[[object$yname]] <- matrix(object$response, ncol=object$ydim[2])
    }else{
      dathelp[[object$yname]] <- object$response
    }
    
    call <- object$callEval
    call$data <- dathelp
    
    # use weights of training data expanded by id to suitable length
    call$weights <- weights[id]
    
    # use call$numInt of original model fit, as weights contains only resampling weights
 
    # Using the offset of object with the following settings 
    # call$control <- boost_control(risk="oobag")
    # call$oobweights <- oobweights[id]
    if(refitSmoothOffset==FALSE && is.null(call$offset) ){
      if(is.null(object$id)){
        call$offset <- matrix(object$offset, ncol=object$ydim[2])[1,]
      }else{
        call$offset <- object$offset
      }
    } 
    # the model is the same for  
    # mod <- object$update(weights = weights, oobweights = oobweights) # (cvrisk) 
    # and
    # mod <- withCallingHandlers(suppressMessages(eval(call)), warning = h)
    # and then the risk can be computed by
    # risk <- mod$risk()[grid] 
    
    ## compute the model by FDboost() - the offset is computed on learning sample
    mod <- withCallingHandlers(suppressMessages(eval(call)), warning = h) # suppress the warning of missing responses    
    mod <- mod[max(grid)]
        
    # compute weights for stanardizing risk, mse, ...
    oobwstand <- lengthTi1[id]*oobweights[id]*intWeights*(1/sum(oobweights))
    
    ############# compute risk, mse, relMSE and mrd
    # get risk function of the family
    riskfct <- get("family", environment(mod$update))@risk
    
    ####################
    ### compute risk and mse without integration weights, like in cvrisk
    risk0 <- sapply(grid, function(g){riskfct( response, mod[g]$fitted(), 
                                               w=oobweights[id])}) / sum(oobweights[id])
    
    mse0 <- simplify2array(mclapply(grid, function(g){
      sum( ((response - mod[g]$fitted())^2*oobweights[id]), na.rm=TRUE )
    }, mc.cores=1) ) /sum(oobweights[id])
    ####################
    
    # oobweights using riskfct() like in mboost, but with different weights!
    risk <- sapply(grid, function(g){riskfct( response, mod[g]$fitted(), w=oobwstand)})
    
    ### mse (mean squared error) equals risk in the case of familiy=Gaussian()
    mse <- simplify2array(mclapply(grid, function(g){
      sum( ((response - mod[g]$fitted())^2*oobwstand), na.rm=TRUE )
    }, mc.cores=1) )
    
#     ### mse2 equals mse in the case of equal grids without missings at the ends
#     mse2 <- simplify2array(mclapply(grid, function(g){
#       sum( ((response - mod[g]$fitted())^2*oobweights[id]*intWeights), na.rm=TRUE )
#     }, mc.cores=1) ) / (sum(oobweights)* (max(mod$yind)-min(mod$yind) ) )

    ### compute overall mean of response in learning sample
    meanResp <- sum(response*intWeights*lengthTi1[id]*weights[id], na.rm=TRUE) / sum(weights)
    
    # # compute overall mean of response in whole sample
    # meanResp <- sum(response*intWeights*lengthTi1[id], na.rm=TRUE) / nObs
    
    ### compute relative mse
    relMSE <- simplify2array(mclapply(grid, function(g){
      sum( ((response - mod[g]$fitted())^2*oobwstand), na.rm=TRUE ) /
        sum( ((response - meanResp)^2*oobwstand), na.rm=TRUE )
    }, mc.cores=1) )
    
    ### mean relative deviation
    resp0 <- response
    resp0[abs(resp0) <= mrdDelete | round(resp0, 1) == 0] <- NA
    mrd <- simplify2array(mclapply(grid, function(g){
      sum( abs(resp0 - mod[g]$fitted())/abs(resp0)*oobwstand, na.rm=TRUE )
    }, mc.cores=1) )
    
    rm(resp0, meanResp)
    
    
    ####### prediction for all observations, not only oob! 
    # -> oob-predictions have to be merged out of predGrid
    predGrid <- predict(mod, aggregate="cumsum", unlist=FALSE)
    
    if(is.null(object$id)){
      predGrid <- sapply(predGrid, function(x) as.vector(x) )[,grid] # save vectors of predictions in matrix
    }else{
      predGrid <- predGrid[,grid] # save vectors of predictions for grid in matrix
    }
    
    ## <FIXME> are those calculations necessary?
    if(is.null(object$id)){
      # predict oob and save responses that were predicted
      predOOB <- predict(mod, aggregate="cumsum", unlist=FALSE)
      keepOOB <- matrix(oobweights, nObs, Gy)
      keepRow <- apply(keepOOB, 1, function(x) any(x!=0))
      if(any(keepRow)){
        if( length(object$yind) > 1){ # functional response
          predOOB <- sapply(predOOB, function(x) as.vector(x[keepRow, ]) )[,grid]
        }else{ # scalar response
          predOOB <- sapply(predOOB, function(x) as.vector(x[keepRow, ]) )[,grid]
        }
        respOOB <- as.vector(matrix(mod$response, nObs, Gy)[keepRow,] )
        attr(respOOB, "curves") <- which(keepRow)
      }else{
        predOOB <- NULL
        respOOB <- NULL
      }
      ### <FIXME> compute those variables for response in long-format?
    }else{
      # predict oob and save responses that were predicted
      predOOB <- predict(mod, aggregate="cumsum", unlist=FALSE)[,grid]
      keepidOOB <- (oobweights!=0)[mod$id]
      if(any(keepidOOB)){
        predOOB <- predOOB[keepidOOB,]
        respOOB <- mod$response[keepidOOB]
        attr(respOOB, "curves") <- unique(mod$id[keepidOOB])
      }else{
        predOOB <- NULL
        respOOB <- NULL
      } 
    }
    
    return(list(risk=risk, predGrid=predGrid, predOOB=predOOB, respOOB=respOOB, 
                mse=mse, relMSE=relMSE, mrd=mrd, risk0=risk0, mse0=mse0, mod=mod))  
  } 
  
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
  ## get out-of-bag mse
  oobmse <- t(sapply(modRisk, function(x) x$mse))  
  ## get out-of-bag relMSE
  oobrelMSE <- t(sapply(modRisk, function(x) x$relMSE))  
  ## get out-of-bag mrd
  oobmrd <- t(sapply(modRisk, function(x) x$mrd))
  
  colnames(oobrisk) <- colnames(oobmse) <- colnames(oobrelMSE) <- colnames(oobmrd) <- grid
  rownames(oobrisk) <- rownames(oobmse) <- rownames(oobrelMSE) <- rownames(oobmrd) <- which(modFitted)

  ## get out-of-bag risk without integration weights
  oobrisk0 <- t(sapply(modRisk, function(x) x$risk0))
  ## get out-of-bag mse without integration weights
  oobmse0 <- t(sapply(modRisk, function(x) x$mse0))
  
  colnames(oobrisk0) <- colnames(oobmse0) <- grid
  rownames(oobrisk0) <- rownames(oobmse0) <- which(modFitted)
  
  ############# check for folds with extreme risk-values at the global median
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
    # <FIXME> only works for regular response! implement depending on id!
    oobpreds0 <- lapply(modRisk, function(x) x$predGrid)
    oobpreds <- matrix(nrow=nrow(oobpreds0[[1]]) , ncol=ncol(oobpreds0[[1]]))
    for(j in 1:length(oobpreds0)){
      oobpreds[folds[,j]==0] <- oobpreds0[[j]][folds[,j]==0] 
    }
    colnames(oobpreds) <- grid
    rm(oobpreds0)
    
  }else{
    oobpreds <- NULL
  }
    
#   # alternative OOB-prediction: works for general folds not only oob
#   predOOB <- lapply(modRisk, function(x) x$predOOB)
#   predOOB <- do.call('rbind', predOOB)
#   colnames(predOOB) <- grid
#   respOOB <- lapply(modRisk, function(x) x$respOOB)
#   respOOB <- do.call('c', respOOB)
#   indexOOB <- lapply(modRisk, function(x) attr(x$respOOB, "curves"))
#   if(is.null(object$id)){
#     indexOOB <- lapply(indexOOB, function(x) rep(x, times=Gy) )
#     indexOOB <- unlist(indexOOB)
#   }else{
#     indexOOB <- names(unlist(indexOOB))[unlist(indexOOB)]
#   }
#   attr(respOOB, "index") <- indexOOB

  coefCV <- list()
  predCV <- list()

  if(getCoefCV){
    
    if(riskopt=="median"){
      print("median")
      optimalMstop <- grid[which.min(apply(oobrisk, 2, median))]
    }else{
      print("mean")
      optimalMstop <- grid[which.min(apply(oobrisk, 2, mean))]
    }  
    
    attr(coefCV, "risk") <- paste("minimize", riskopt, "risk")
    
    ### estimates of coefficients
    timeHelp <- seq(min(modRisk[[1]]$mod$yind), max(modRisk[[1]]$mod$yind), l=40)
    for(l in 1:length(modRisk[[1]]$mod$baselearner)){
      # estimate the coefficients for the model of the first fold
      coefCV[[l]] <- coef(modRisk[[1]]$mod[optimalMstop], 
                          which=l, n1 = 40, n2 = 20, n3 = 15, n4 = 10)$smterms[[1]]
#       if(l==1){
#         coefCV[[l]]$offset <- matrix(ncol=40, nrow=length(modRisk))
#         coefCV[[l]]$offset[1,] <- modRisk[[1]]$mod$predictOffset(time=timeHelp)
#       }
      attr(coefCV[[l]]$value, "offset") <- NULL # as offset is the same within one model

      # add estimates for the models of the other folds
      coefCV[[l]]$value <- lapply(1:length(modRisk), function(g){
        ret <- coef(modRisk[[g]]$mod[optimalMstop], 
                    which=l, n1 = 40, n2 = 20, n3 = 15, n4 = 10)$smterms[[1]]$value
#         if(l==1){
#           coefCV[[l]]$offset[g,] <- modRisk[[g]]$mod$predictOffset(time=timeHelp)
#         }
        attr(ret, "offset") <- NULL # as offset is the same within one model
        return(ret)
      })
    }
    
    ## predict offset
    offset <- sapply(1:length(modRisk), function(g){
      # offset is vector of length yind or numeric of length 1 for constant offset
      ret <- modRisk[[g]]$mod$predictOffset(time=timeHelp)
      if( length(ret)==1 & length(object$yind) > 1 ) ret <- rep(ret, length(timeHelp))
      return(ret)
    })
    
    attr(coefCV, "offset") <- offset
    
    niceNames <- c("offset", lapply(coefCV, function(x) x$main))
    
    ### predictions of terms based on the coefficients for each model
    # only makes sense for type="curves" with leaving-out one curve per fold!!
    if(grepl( "curves", type)){
      for(l in 1:(length(modRisk[[1]]$mod$baselearner)+1)){
        predCV[[l]] <- t(sapply(1:length(modRisk), function(g){
          if(l==1){ # save offset of model
            # offset is vector of length yind or numeric of length 1 for constant offset
            ret <- modRisk[[g]]$mod[optimalMstop]$predictOffset(object$yind) 
            # regular data or scalar response
            if(is.null(object$id)){
              if( length(ret)==1 ) ret <- rep(ret, modRisk[[1]]$mod$ydim[2])
            # irregular data
            }else{
              if( length(ret)==1 ){ ret <- rep(ret, sum(object$id==g)) }else{ ret <- ret[object$id==g] }
            }
          }else{ # other effects
            ret <- predict(modRisk[[g]]$mod[optimalMstop], which=l-1) # model g
            if(!(l-1) %in% selected(modRisk[[g]]$mod[optimalMstop]) ){ # effect was never chosen
              if(is.null(object$id)){
                ret <- matrix(0, ncol=modRisk[[1]]$mod$ydim[2], nrow=modRisk[[1]]$mod$ydim[1])
              }else{
                ret <- matrix(0, nrow=length(object$id), ncol=1)
              } 
            }
            if(is.null(object$id)){
              ret <- ret[g,] # save g-th row = preds for g-th observations
            }else{
              ret <- ret[object$id==g,] # save preds of g-th observations
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
  
  ret <- list(response=response, yind=object$yind, id=object$id,
              folds=folds, grid=grid,
              coefCV=coefCV,
              predCV=predCV,
              oobpreds=oobpreds, 
              #predOOB=predOOB, respOOB=respOOB,
              oobrisk=oobrisk, 
              oobriskMean=colMeans(oobrisk),
              oobmse=oobmse,
              oobrelMSE=oobrelMSE,
              oobmrd=oobmrd,
              oobrisk0=oobrisk0, 
              oobmse0=oobmse0)
  
  class(ret) <- "validateFDboost"
  
  return(ret) 
}



#' @rdname plot.validateFDboost
#' @method mstop validateFDboost
#' @export
#' 
# Function to extract the optimal stopping iteration
mstop.validateFDboost <- function(object, riskopt=c("mean", "median"), ...){
  
  dots <- list(...)
  riskopt <- match.arg(riskopt)
  
  if(riskopt=="median"){
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
#' prediction error of a model fitted by FDboost.
#' 
#' @param x object of class validateFDboost
#' @param object object of class validateFDboost
#' @param riskopt how the risk is minimized to obtain the optimal stopping iteration; 
#' defaults to the mean, can be changed to the median.
#' @param modObject if the original model object of class \code{FDboost} is given 
#' predicted values of the whole model can be compared to the predictions of the cross-validated models
#' @param predictNA should missing values in the response be predicted? Defaults to FALSE.
#' @param names.arg names of the observed curves
#' @param ask par(ask=ask)
#' @param which a subset of base-learners to take into account for plotting. 
#' In the case of \code{plot.validateFDboost} the diagnostic plots that are given. 
#' @param pers plot coefficient surfaces as persp-plots? Defaults to TRUE.
#' @param commonRange, plot predicted coefficients on a common range, defaults to TRUE
#' @param showQuantiles plot the 0.05 and the 0.95 Quantile of coefficients in 1-dim effects
#' @param showNumbers show number of curve in plot of predicted coefficients, defaults to FALSE
#' @param terms, logical, defaults to TRUE plot the added terms (default) or the coefficients?
#' @param probs vector of quantiles to be used in the plotting of 2-dimensional coefficients surfaces,
#' defaults to \code{probs=c(0.25, 0.5, 0.75)}
#' @param ylim values for limits of y-axis
#' @param ... additional arguments passed to callies.
#' 
#' @details \code{plot.validateFDboost} plots cross-validated risk, RMSE, MRD, measured and predicted values 
#' and residuals as determined by \code{validateFDboost}.
#' \code{mstop.validateFDboost} extracts the optimal mstop by minimizing the median risk.
#' 
#' @aliases mstop.validateFDboost
#' 
#' @method plot validateFDboost
#' 
#' @export
plot.validateFDboost <- function(x, riskopt=c("mean", "median"), 
                                 which=1, 
                                 modObject=NULL, predictNA=FALSE, 
                                 names.arg=NULL, ask=TRUE, ...){
  
  # get the optimal mstop
  riskopt <- match.arg(riskopt)
  mopt <- mstop(x, riskopt=riskopt)
  # get the position in the grid of mopt
  mpos <- which(x$grid==mopt)
  
  if(length(which)>1) par(ask=ask)
  
  if(1 %in% which){
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
  }
  
  if(any(c(2,3) %in% which)){
    # Plot RMSE and MRD for optimal mstop
    if(is.null(names.arg)){
      names.arg <- seq(along=x$oobrisk[,mpos])
    }
    stopifnot(length(names.arg)==length(x$oobrisk[,mpos]))
  }
  
  # plot risk
  if(2 %in% which){
    barplot(x$oobrisk[,mpos], main="risk", names.arg = names.arg, las=2) # 
    #abline(h=x$rmse[mpos], lty=2)
  }
  
  # plot MRD
  if(3 %in% which){
    barplot(x$oobmrd[,mpos], main="MRD", names.arg = names.arg, las=2)
    #abline(h=x$mrd[mpos], lty=2)
  }
    
  # Plot the predictions for the optimal mstop
  if(!is.null(x$oobpreds[,mpos])){
    response <- x$response
    pred <- x$oobpreds[,mpos]
    if(!predictNA){
      pred[is.na(response)] <- NA
    }
    
    predMat <- pred
    responseMat <- response
    
    if(is.null(x$id)) predMat <- matrix(pred, ncol=length(x$yind))
    if(is.null(x$id)) responseMat <- matrix(response, ncol=length(x$yind))
    
    if(4 %in% which){
      ylim <- range(response, pred, na.rm = TRUE)
      funplot(x$yind, responseMat, id=x$id, lwd = 1, pch = 1, ylim = ylim,  
              ylab = "", xlab = attr(x$yind, "nameyind"), 
              main="Measured and Predicted Values", ...)
      funplot(x$yind, predMat, id=x$id, lwd = 2, pch = 2, add = TRUE, ...)
      posLegend <- "topleft"
      legend(posLegend, legend=c("observed","predicted"), col=1, pch=1:2)  
    }
    
    if(5 %in% which){
      # Plot residuals for the optimal mstop
      funplot(x$yind, responseMat-predMat, id=x$id, ylab = "", xlab = attr(x$yind, "nameyind"), 
              main="Residuals", ...)
      abline(h=0, lty=2, col="grey")
    }
  }
  
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
plotPredCoef <- function(x, which=NULL, pers=TRUE,
                         commonRange=TRUE, showNumbers=FALSE, showQuantiles=TRUE,
                         ask=TRUE, 
                         terms=TRUE,         
                         probs=c(0.05, 0.5, 0.95), # quantiles of variables to use for plotting
                         ylim=NULL, ...){
  
  stopifnot(any(class(x)=="validateFDboost"))
  
  ### Get further arguments passed to the different plot-functions
  dots <- list(...)
  
  getArguments <- function(x, dots=dots){
    if(any(names(dots) %in% names(x))){
      dots[names(dots) %in% names(x)]
    }else list()
  }
  
  #argsPlot <- getArguments(x=formals(graphics::plot.default), dots=dots)
  argsPlot <- getArguments(x=c(formals(graphics::plot.default), par()), dots=dots)
  argsMatplot  <- getArguments(x=c(formals(graphics::matplot), par()), dots=dots)
  argsFunplot  <- getArguments(x=c(formals(funplot), par()), dots=dots)

  argsPersp <- getArguments(x=formals(getS3method("persp", "default")), dots=dots)
  
  plotWithArgs <- function(plotFun, args, myargs){        
    args <- c(myargs[!names(myargs) %in% names(args)], args)        
    do.call(plotFun, args)            
  }    
  
  if(is.null(which)) which <- 1:length(x$coefCV)
  
  if(length(which)>1) par(ask=ask)
  
  if(terms){
    
    if(all(which==1:length(x$coefCV))){
      which <- 1:(length(x$coefCV)+1)
    }else{
      which <- which+1 
    }
    
    if( length(x$predCV)==0 ){
      warning("terms=TRUE, but predCV is empty.")
      return(NULL)
    }
    
    if(commonRange & is.null(ylim)){ 
      ylim <- range(x$predCV[which])
    }
    # <FIXME> use extra ylim for functional offset?
    
    ### loop over base-learners
    for(l in which){
      
      if( is.null(x$id) ){
        
        if(length(x$yind)>1){
          if(l==1 && length(x$yind)>1){ 
            matplot(x$yind, t(x$predCV[[l]]), type="l", 
                    main=names(x$predCV)[l], xlab=attr(x$yind, "nameyind"), ylab="coef", ylim=ylim, ...)
          }else{
            matplot(x$yind, t(x$predCV[[l]]), type="l", 
                    main=names(x$predCV)[l], xlab=attr(x$yind, "nameyind"), ylab="coef", ylim=ylim, ...)
            
            if(showNumbers){
              matplot(x$yind, t(x$predCV[[l]]), add=TRUE, ...)
            }
          }
        }else{
          plot(1:length(x$predCV[[l]]), x$predCV[[l]], main=names(x$predCV)[l], xlab="index", ylab="coef")
        }
        
      }else{
        funplot(x$yind, unlist(x$predCV[[l]]), 
                id=x$id, type="l", 
                main=names(x$predCV)[l], xlab=attr(x$yind, "nameyind"), ylab="coef", ylim=ylim, ...)
        
        if(showNumbers){
          points(x$yind, unlist(x$predCV[[l]]), type="p", pch=paste(x$id))
        }
      } 
    } # end for-loop 
    
  }else{ # plot coefficients
    
    if(commonRange & is.null(ylim)){
      if(length(x$yind)>1){
        ylim <- range(lapply(x$coefCV[which], function(x) x$value)) 
      }else{
        ylim <- range(lapply(x$coefCV[which], function(x) x$value))
      } 
    }
    
    for(l in which){ # loop over effects
      
      # coef() of a certain term
      temp <- x$coefCV[l][[1]]
      
      # set the range for each effect individually
      if(FALSE) ylim <- range(temp$value)
      
      ## plot the estimated offset
      if(l==1 && length(x$yind)>1){
        myMat <- attr(x$coefCV, "offset")
        
        if(showQuantiles){ 
          
          timeHelp <- seq(min(x$yind), max(x$yind), length=nrow(myMat))
          plotWithArgs(matplot, args=argsMatplot, 
                       myargs=list(x=timeHelp, y=myMat, type="l", xlab=temp$xlab,
                                   main=temp$main, ylab="coef", ylim=NULL, 
                                   col=rgb(0.6,0.6,0.6, alpha=0.5)))
          
          if(showNumbers){
            matplot(timeHelp, myMat, add=TRUE, col=rgb(0.6,0.6,0.6, alpha=0.5), ...)
          }
          lines(timeHelp, rowMeans(myMat), col=1, lwd=2)
          lines(timeHelp, apply(myMat, 1, quantile, 0.95), col=2, lwd=2, lty=2)
          lines(timeHelp, apply(myMat, 1, quantile, 0.05), col=2, lwd=2, lty=2)
          
        }else{
          
          plotWithArgs(matplot, args=argsMatplot, 
                       myargs=list(x=timeHelp, y=myMat, type="l", xlab=temp$xlab,
                                   main=temp$main, ylab="coef", ylim=NULL))
          
          if(showNumbers){
            matplot(timeHelp, myMat, add=TRUE, ...)
          }
        }
      }
      
      if(temp$dim==1){
        # impute vector of 0 if effect was never chosen
        temp$value[sapply(temp$value, function(x) length(x)==1)] <- list(rep(0, 40))
        myMat <- sapply(temp$value, function(x) x) # as one matrix
        
        if(showQuantiles){
          
          plotWithArgs(matplot, args=argsMatplot, 
                       myargs=list(x=temp$x, y=myMat, type="l", xlab=temp$xlab,
                                   main=temp$main, ylab="coef", ylim=ylim, 
                                   col=rgb(0.6,0.6,0.6, alpha=0.5)))
          
          if(showNumbers){
            matplot(temp$x, myMat, add=TRUE, col=rgb(0.6,0.6,0.6, alpha=0.5), ...)
          }
          
          lines(temp$x, rowMeans(myMat), col=1, lwd=2)
          lines(temp$x, apply(myMat, 1, quantile, 0.95), col=2, lwd=2, lty=2)
          lines(temp$x, apply(myMat, 1, quantile, 0.05), col=2, lwd=2, lty=2)
          
        }else{
          
          plotWithArgs(matplot, args=argsMatplot, 
                       myargs=list(x=temp$x, y=myMat, type="l", xlab=temp$xlab,
                                   main=temp$main, ylab="coef", ylim=ylim))
          
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
        
        # set lower triangular matrix to NA for historic effect
        if(grepl("bhist", temp$main)){
          for(k in 1:length(temp$value)){
            temp$value[[k]][outer(1:nrow(temp$value[[k]]), 1:ncol(temp$value[[k]]), "<=")==FALSE] <- NA
          }
        }
        
        if(!is.factor(temp$x)){ # temp$x is metric
          
          # impute matrix of 0 if effect was never chosen
          temp$value[sapply( temp$value, function(x){ is.null(dim(x)) || dim(x)[2]==1 })] <- list(matrix(0, ncol=20, nrow=20))
          
          # plot coefficient surfaces at different pointwise quantiles
          if(pers){ 
            matvec <- sapply(temp$value, c)
            for(k in 1:length(probs)){
              tempZ <- matrix(apply(matvec, 1, quantile, probs=probs[k], na.rm=TRUE), ncol=length(temp$x))

              plotWithArgs(persp, args=argsPersp, 
                           myargs=list(x=temp$x, y=temp$y, z=tempZ,
                                       ticktype="detailed", theta=30, phi=30,
                                       xlab=temp$xlab, ylab=temp$ylab, zlab="coef", 
                                       zlim=if(is.null(ylim)) range(matvec, na.rm=TRUE) else ylim,  
                                       main=paste(temp$main, " at ", probs[k]*100, "%-quantile", sep=""), 
                                       col=getColPersp(tempZ)))
            }  
          }else{ # do 2-dim plots
            
            for(j in 1:length(quanty)){ 
              
              myCol <- sapply(temp$value, function(x) x[, quanty[j]==temp$y]) # first column
              
              if(showQuantiles){ 
                #browser()
                plotWithArgs(matplot, args=argsMatplot, 
                             myargs=list(x=temp$y, y=myCol, type="l", xlab=temp$ylab,
                                         main=paste(temp$main, " at ", probs[j]*100, "% of ", temp$xlab, sep=""), 
                                         ylab="coef", ylim=ylim, 
                                         col=rgb(0.6,0.6,0.6, alpha=0.5)))

                if(showNumbers){
                  matplot(temp$y, myCol, add=TRUE, col=rgb(0.6,0.6,0.6, alpha=0.5), ...)
                }
                lines(temp$y, rowMeans(myCol), col=1, lwd=2)
                lines(temp$y, apply(myCol, 1, quantile, 0.95, na.rm = TRUE), col=2, lwd=2, lty=2)
                lines(temp$y, apply(myCol, 1, quantile, 0.05, na.rm = TRUE), col=2, lwd=2, lty=2)
                
              }else{
                plotWithArgs(matplot, args=argsMatplot, 
                             myargs=list(x=temp$y, y=myCol, type="l", xlab=temp$ylab,
                                         main=paste(temp$main, " at ", probs[j]*100, "% of ", temp$xlab, sep=""), 
                                         ylab="coef", ylim=ylim))
                if(showNumbers){
                  matplot(temp$y, myCol, add=TRUE, ...)
                }
              }
              
            } # end loop over quanty
            
            for(j in 1:length(quantx)){  
              myRow <- sapply(temp$value, function(x) x[quantx[j]==temp$x, ]) # first column
              
              if(showQuantiles){
                
                plotWithArgs(matplot, args=argsMatplot, 
                             myargs=list(x=temp$x, y=myRow, type="l", xlab=temp$xlab,
                                         main=paste(temp$main, " at ", probs[j]*100, "% of ", temp$ylab, sep=""), 
                                         ylab="coef", ylim=ylim, 
                                         col=rgb(0.6,0.6,0.6, alpha=0.5)))
                if(showNumbers){
                  matplot(temp$x, myRow, add=TRUE, col=rgb(0.6,0.6,0.6, alpha=0.5), ...)
                }
                lines(temp$x, rowMeans(myRow), col=1, lwd=2)
                lines(temp$x, apply(myRow, 1, quantile, 0.95, na.rm = TRUE), col=2, lwd=2, lty=2)
                lines(temp$x, apply(myRow, 1, quantile, 0.05, na.rm = TRUE), col=2, lwd=2, lty=2)
                
              }else{
                plotWithArgs(matplot, args=argsMatplot, 
                             myargs=list(x=temp$x, y=myRow, type="l", xlab=temp$xlab,
                                         main=paste(temp$main, " at ", probs[j]*100, "% of ", temp$ylab, sep=""), 
                                         ylab="coef", ylim=ylim))                
                if(showNumbers){
                  matplot(temp$x, myRow, add=TRUE, ...)
                }
              }
            }
            
          } # end else
          
        }else{ # temp$x is factor
          for(j in 1:length(quantx)){ 
            
            # impute matrix of 0 if effect was never chosen
            temp$value[sapply(temp$value, function(x) is.null(dim(x)))] <- list(matrix(0, ncol=20, nrow=length(quantx)))
            
            if(is.null(temp$z)){
              
              myRow <- sapply(temp$value, function(x) x[quantx[j]==temp$x, ]) # first column
              
              if(showQuantiles){
                
                plotWithArgs(matplot, args=argsMatplot, 
                             myargs=list(x=temp$y, y=myRow, type="l", xlab=temp$ylab,
                                         main=paste(temp$main, " at ", temp$xlab,"=" ,quantx[j], sep=""), 
                                         ylab="coef", ylim=ylim, 
                                         col=rgb(0.6,0.6,0.6, alpha=0.5)))

                if(showNumbers){
                  matplot(temp$y, myRow, add=TRUE, col=rgb(0.6,0.6,0.6, alpha=0.5), ...)
                }
                lines(temp$y, rowMeans(myRow), col=1, lwd=2)
                lines(temp$y, apply(myRow, 1, quantile, 0.95, na.rm = TRUE), col=2, lwd=2, lty=2)
                lines(temp$y, apply(myRow, 1, quantile, 0.05, na.rm = TRUE), col=2, lwd=2, lty=2)
                
              }else{
                matplot(temp$y, myRow, type="l", xlab=temp$ylab, ylim=ylim,
                        main=paste(temp$main, " at ", temp$xlab,"=" ,quantx[j], sep=""), ylab="coef", ...)
              } 
              
            }else{
              quantz <- temp$z
              myRow <- sapply(temp$value, function(x) x[quantx[j]==temp$x & quantz[j]==temp$z, ]) # first column
              
              if(showQuantiles){
                
                plotWithArgs(matplot, args=argsMatplot, 
                             myargs=list(x=temp$y, y=myRow, type="l", xlab=temp$ylab,
                                         main=paste(temp$main, " at ", temp$xlab,"=" , quantx[j], ", " , 
                                                    temp$zlab,"=" ,quantz[j], sep=""), 
                                         ylab="coef", ylim=ylim, 
                                         col=rgb(0.6,0.6,0.6, alpha=0.5)))
                
#                 matplot(temp$y, myRow, type="l", xlab=temp$xlab,
#                         main=paste(temp$main, " at ", temp$xlab,"=" , quantx[j], ", " , 
#                                    temp$zlab,"=" ,quantz[j], sep=""), 
#                         ylab="coef", ylim=ylim, 
#                         col=rgb(0.6,0.6,0.6, alpha=0.5), ...) 
                if(showNumbers){
                  matplot(temp$y, myRow, add=TRUE, col=rgb(0.6,0.6,0.6, alpha=0.5), ...)
                }
                lines(temp$y, rowMeans(myRow), col=1, lwd=2)
                lines(temp$y, apply(myRow, 1, quantile, 0.95, na.rm = TRUE), col=2, lwd=2, lty=2)
                lines(temp$y, apply(myRow, 1, quantile, 0.05, na.rm = TRUE), col=2, lwd=2, lty=2)
                
              }else{
                plotWithArgs(matplot, args=argsMatplot, 
                             myargs=list(x=temp$y, y=myRow, type="l", xlab=temp$ylab,
                                         main=paste(temp$main, " at ", temp$xlab,"=" , quantx[j], ", " , 
                                                    temp$zlab,"=" ,quantz[j], sep=""), 
                                         ylab="coef", ylim=ylim))
              }
              
            }
            if(showNumbers){
              matplot(temp$y, myRow, add=TRUE, ...)
            }
          } 
        }
        
      } # end if(temp$dim==2)
      
    } # end loop over effects
  }
  par(ask=FALSE)
}


