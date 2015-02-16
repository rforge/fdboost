
#' Summary of a boosted functional regression model 
#' 
#'  Takes a fitted \code{FDboost}-object and produces a summary.
#' 
#' @param object a fitted \code{FDboost}-object
#' @param ... currently not used
#' @seealso \code{\link{FDboost}} for the model fit.
#' @return a list with summary information 
#' @method summary FDboost
#' @export
### similar to summary.mboost()
summary.FDboost <- function(object, ...) {
  
  ret <- list(object = object, selprob = NULL)
  xs <- selected(object)
  nm <- variable.names(object)
  selprob <- tabulate(xs, nbins = length(nm)) / length(xs)
  names(selprob) <- names(nm)
  selprob <- sort(selprob, decreasing = TRUE)
  ret$selprob <- selprob[selprob > 0]
  class(ret) <- "summary.FDboost"
  
  ### only show one unique offset value
  #ret$object$offset <- unique(ret$object$offset)
  
  return(ret)
}

#' Print a boosted functional regression model 
#' 
#'  Takes a fitted \code{FDboost}-object and produces a print on the console.
#' 
#' @param x a fitted \code{FDboost}-object
#' @param ... currently not used
#' @seealso \code{\link{FDboost}} for the model fit.
#' @return a list with information on the model 
#' @method print FDboost
#' @export
### similar to print.mboost()
print.FDboost <- function(x, ...) {
  
  cat("\n")
  cat("\t Model-based Boosting with Functional Response\n")
  cat("\n")
  if (!is.null(x$call))
    cat("Call:\n", deparse(x$call), "\n\n", sep = "")
  show(x$family)
  cat("\n")
  cat("Number of boosting iterations: mstop =", mstop(x), "\n")
  cat("Step size: ", x$control$nu, "\n")
  
  if(length(unique(x$offset))<10){
    cat("Offset: ", round(unique(x$offset), 3), "\n")
  }else{
    cat("Offset: ", round(unique(x$offset), 3)[1:3], "..." ,
        round(unique(x$offset), 3)[length(unique(x$offset))-3+1:3], "\n")
  }

  cat("Number of baselearners: ", length(variable.names(x)), "\n")
  cat("\n")
  invisible(x)
  
}


#' Prediction for boosted functional regression model 
#' 
#'  Takes a fitted \code{FDboost}-object produced by \code{\link{FDboost}()} and produces 
#'  predictions given a new set of values for the model covariates or the original 
#'  values used for the model fit. This is a wrapper
#'  function for \code{\link[mboost]{predict.mboost}()}
#' 
#' @param object a fitted \code{FDboost}-object
#' @param newdata  A named list or a data frame containing the values of the model 
#' covariates at which predictions are required.
#' If this is not provided then predictions corresponding to the original data are returned. 
#' If \code{newdata} is provided then it should contain all the variables needed for 
#' prediction, in the format supplied to \code{FDboost}, i.e., 
#' functional predictors must be supplied as matrices with each row corresponding to 
#' one observed function.
#' @param which a subset of base-learners to take into account for computing predictions 
#' or coefficients. If which is given (as an integer vector corresponding to base-learners) 
#' a list is returned. 
#' @param unlist logical, defaults to TRUE. Should predictions be returned in matrix form 
#' (default) or as a list 
#' @param ...  additional arguments passed on to \code{\link[mboost]{predict.mboost}()}.
#' 
#' @seealso \code{\link{FDboost}} for the model fit 
#' and \code{\link{plotPredicted}} for a plot of the observed values and their predictions.
#' @return a matrix or list of predictions depending on values of unlist and which 
#' @method predict FDboost
#' @export
# predict function: wrapper for predict.mboost()
## <TODO> check which
## <TODO> check unlist
predict.FDboost <- function(object, newdata = NULL, which=NULL, unlist=TRUE, ...) {
  
  stopifnot(any(class(object)=="FDboost")) 
  # print("Prediction FDboost") 
  
  classObject <- class(object)
  class(object) <- "mboost"
  
  ### Prepare data so that the function predict.mboost() can be used
  if(!is.null(newdata)){
    
    # check class of elements in newdata: objects of class "matrix" have to be saved "AsIs"
    if(class(newdata) != "data.frame"){
      newdata <- lapply(newdata, function(x) if(class(x)=="matrix") return(I(x)) else return(x)  )      
    } 
    
    # FIXME: allow for more general which in combination with newdata
    if(length(which)>1) stop("If newdata not equals NULL only which=NULL or which of length 1 can be specified.")
    
    # dummy variable to fit the intercept
    n <- NROW(newdata[[1]])
    newdata[["ONEx"]] <- rep(1.0, n)
    
    # save time-variable (index over response)
    nameyind <- attr(object$yind, "nameyind")
    lengthYind <- length(newdata[[nameyind]])
    if(lengthYind==0) stop("Index of response, ", nameyind, ", must be specified and have length >0.")
    assign(nameyind, newdata[[nameyind]])
    newdata$ONEtime <- rep(1.0, lengthYind)
    
    # In the case of bsignal(), bconcurrent() and bhist() it is necessary 
    # to add the index of the signal-matrix as attribute
    posBsignal <- grep("bsignal", names(object$baselearner))
    posBconc <- grep("bconcurrent", names(object$baselearner))
    posBhist <- grep("bhist", names(object$baselearner))
    whichHelp <- which
    if(is.null(which)) whichHelp <- 1:length(object$baselearner)
    posBsignal <- whichHelp[whichHelp %in% posBsignal]
    posBconc <- whichHelp[whichHelp %in% posBconc]
    posBhist <- whichHelp[whichHelp %in% posBhist]
    newdataConc <- list() # data to predict concurrent effect
    newdataHist <- list() # data to predict historic effect 
    indname_all <- c() # save the names of all indices
    if(length(c(posBsignal, posBconc, posBhist))>0){
      #if(!is.list(newdata)) newdata <- list(newdata)
      for(i in c(posBsignal, posBconc, posBhist)){
        form <- strsplit(object$baselearner[[i]]$get_call(), "%O%")[[1]][1]
        form <- gsub("\"", "", form, fixed=TRUE) # delete all quatation marks
        form <- gsub("\\", "", form, fixed=TRUE) # delete all backslashes
        formula_help <- formula(paste("resHelp ~", form))
        xname <- all.vars(formula_help)[2]
        indname <- if(length(all.vars(formula_help))>=3) all.vars(formula_help)[3] else "xindDefault" 
        indname_all <- c(indname_all, indname)
        if(i %in% c(posBhist, posBconc)){
          indnameY <- attr(object$yind, "nameyind")
        } else{
          indnameY <- NULL
        }
        
        #if(length(newdata[[indname]])!=ncol(newdata[[xname]])){
        #  stop(paste("Dimensions of index", indname, "and signal", xname, "do not match."))
        #} 
        
        attr(newdata[[xname]], "indname") <- indname
        attr(newdata[[xname]], "xname") <- xname
        attr(newdata[[xname]], "signalIndex") <-  if(indname!="xindDefault") newdata[[indname]] else seq(0,1,l=ncol(newdata[[xname]]))
        
        if(i %in% c(posBhist, posBconc)){
          attr(newdata[[indnameY]], "indnameY") <-  indnameY
          attr(newdata[[xname]], "indexY") <-  if(indnameY!="xindDefault") newdata[[indnameY]] else seq(0,1,l=ncol(newdata[[xname]]))
          if(any(classObject=="FDboostLong")){
            id <- newdata[[ attr(object$id, "nameid") ]]
            #if(is.null(id)) warning("There is no id-variable called,", attr(object$id, "nameid"),  "in newdata.")
            
            attr(newdata[[xname]], "id") <-  id
          } 
        } 
        
        # save data of concurrent effects
        if(i %in% posBconc) newdataConc[[xname]] <- newdata[[xname]]
        
        # save data of historic effects
        if(i %in% posBhist) newdataHist[[xname]] <- newdata[[xname]]
        
        ## delete signalIndex if newdata is a list 
        # - not a good idea if you have the same index for several functional covariates
        #if(!is.data.frame(newdata) & grepl("signal", form) ){
        #  newdata <- newdata[-which(names(newdata)==indname)]
        #} 
      }
    }
    
    ## Check dimensions of newdata
    #if(length(unique(sapply(newdata, NROW))) != 1){
    #  stop(paste("Dimensions of newdata do not match:", 
    #             paste(names(newdata), sapply(newdata, NROW), sep=": ", collapse=", ")))
    #}
    
    ### <ToDo> check compatibility of type = c("link", "response", "class")
    ### and aggregate = c("sum", "cumsum", "none") when they are not in the default!
    
    # Predict concurrent effect extra
    predMboostConc <- 0
    if(length(posBconc) > 0){ 
      if(length(unique(sapply(newdataConc, NROW))) != 1) stop(paste("Dimensions of newdata do not match:", 
                                                                    paste(names(newdataConc), sapply(newdataConc, NROW), sep=": ", collapse=", ")))
      predMboostConc <- rowSums(predict(object=object, newdata=as.data.frame(newdataConc), 
                                        which=posBconc, ...))
      whichHelp <- whichHelp[-which(whichHelp %in% posBconc)]
    }
    
    # Predict historic effect extra
    predMboostHist <- 0
    if(length(posBhist) > 0){ 
      if(length(unique(sapply(newdataHist, NROW))) != 1) stop(paste("Dimensions of newdata do not match:", 
                                                                    paste(names(newdataHist), sapply(newdataHist, NROW), sep=": ", collapse=", ")))
      predMboostHist <- predict(object=object, newdata=as.data.frame(newdataHist), 
                                        which=posBhist, ...)
      predMboostHist <- rowSums(predMboostHist)
      whichHelp <- whichHelp[-which(whichHelp %in% posBhist)]
    }
    
    # Predict effect of offset
    predOffset <- object$offsetVec # offset is just an integer 
    #if(1 %in% whichHelp && grepl("ONEx", names(object$baselearner)[[1]])){
    if(length(object$offsetVec)>1){ # offset is a smooth function
      if(is.null(object$id)){
        predOffset <- rep(object$predictOffset(newdata[[nameyind]]), each=n)
      }else{
        predOffset <- object$predictOffset(newdata[[nameyind]])
      }  
      names(predOffset) <- NULL
    }
            
    # Prediction using the function predict.mboost() without concurrent and historic effects
    predMboost0 <- 0
    
    if(length(whichHelp) > 0){
      if(!is.null(object$ydim)){ # regular grid
        # for for models containing %O%, mboost$predict allows to give newdata as a list
        predMboost0 <- predict(object=object, newdata=newdata, which=whichHelp, ...)
        predMboost0 <- rowSums(predMboost0)
      }else{ # long format
        # for models WITHOUT %O%, mboost$predict expects a data.frame;
        # a list is NOT possible
        ## only keep the necessary variables in the dataframe
        vars <- all.vars(formula(object$formulaMboost)[[3]])
        if(is.null(object$id)){ # get rid of id and index of Y in the case of regular data 
          vars <- vars[!vars %in% c(attr(object$id, "nameid"), indname_all)]
        }
        vars <- vars[vars %in% names(newdata)]
        
        ## check whether t has the same length as the variables
        if(length(vars) > 2 & "ONEx" %in% vars){ newdata[["ONEx"]] <- rep(1L, NROW(newdata[[vars[2]]])) }
        if(length(unique(lapply(newdata[vars], NROW)))!=1){
          stop("Only can predict with newdata for irregular response for equal NROW of all variables and the index of the response.")  
        } 
        newdata0 <- as.data.frame(newdata[vars])
        # add a variable t, for y(t) to avoid error that the variable is not in dataframe
        # newdata0[,attr(object$yind, "nameyind")] <- seq(min(object$yind), max(object$yind), 
        #                                                l=nrow(newdata0))
        
        # predict(object=object, newdata=data.frame(newdata[vars]), which=whichHelp)
        predMboost0 <- predict(object=object, newdata=newdata0, which=whichHelp, ...)
        predMboost0 <- rowSums(predMboost0)
      }
    }
    
    if(length(predMboostConc)!=1 & length(predMboost0)!=1) stopifnot(dim(predMboostConc)==dim(predMboost0))
    if(length(predMboostHist)!=1 & length(predMboost0)!=1) stopifnot(dim(predMboostHist)==dim(predMboost0))
    
    # Sum up the prediction of the offset, concurrent effects, historic effects and other effects
    # if the whole model is predicted
    if(is.null(which)){
      predMboost <-  predMboost0 + predMboostConc + predMboostHist + predOffset
    }else{
      predMboost <-  predMboost0 + predMboostConc + predMboostHist
      attr(predMboost, "offset") <- predOffset
      # matrix(predOffset, ncol=ncol(predMboost)) 
    }  
    
#     # add the scalar offset in the case that the intercept is part of the prediction
#     if(length(predOffset) == 1 && 
#          ( 1 %in% whichHelp && grepl("ONEx", names(object$baselearner)[[1]])  ){
#       predMboost <- predMboost + object$offset
#     } 
            
  }else{ # is.null(newdata)
    n <- object$ydim[1]
    lengthYind <- object$ydim[2] 

    # predict.mboost() does not like a 0 in which
    if(length(which)>1 && 0 %in% which){
      which <- which[which!=0]
      predMboost <- predict(object=object, newdata=NULL, which=which, ...) #
    }
    
    # Prediction using the function predict.mboost() 
    predMboost <- predict(object=object, newdata=NULL, which=which, ...)
    ## <SB> as newdata=NULL, the offset of the original model is correct for the prediction
    
#     # add the scalar offset in the case that the intercept is part of the prediction
#     if(length(predOffset) == 1 && 
#          ( 1 %in% whichHelp && grepl("ONEx", names(object$baselearner)[[1]])  ){
#            predMboost <- predMboost + object$offset
#          }     
  }
  
  # save the offset as matrix if it is a vector
  offsetTemp <- attr(predMboost, "offset")
  
  # model-estimation in long format, by specifying id
  if(is.null(object$ydim)) return(predMboost)
    
  ### Reshape prediction of mboost for FDboost
  if(!is.list(predMboost)){
    predMboost <- list(predMboost)
    listByHand <- TRUE
  }else listByHand <- FALSE
  
  # Reshape prediciton in vector to the prediction as matrix
  reshapeToMatrix <- function(pred){
    if(length(pred)==1) return(pred)
    if(is.matrix(pred)){
      fit <- list()
      # reshape the fit into a list of matrices
      for(terms in 1:ncol(pred)){
        fit[[terms]] <- matrix(pred[,terms], 
                               nrow=n, ncol=lengthYind, byrow=FALSE)
      }
      names(fit) <- colnames(pred)
    }else{
      fit <- list(matrix(pred, nrow=n, ncol=lengthYind, byrow=FALSE))
    }                   
    return(fit)    
  }
  
  ret <- vector("list", length=length(predMboost))
  for(i in 1:length(predMboost)){
    ret[[i]] <- reshapeToMatrix(pred=predMboost[[i]])
  } 
  if(listByHand)  ret <- unlist(ret, recursive=FALSE)
  if(unlist & length(ret)==1) ret <- ret[[1]]
  if(is.null(names(ret)) & length(ret)==length(predMboost)) names(ret) <- names(predMboost)
  if(!is.null(which)) attr(ret, "offset") <- reshapeToMatrix(offsetTemp)[[1]]
  return(ret) 
}


#' Fitted values of a boosted functional regression model 
#' 
#'  Takes a fitted \code{FDboost}-object and computes the fitted values.
#' 
#' @param object a fitted \code{FDboost}-object
#' @param ... additional arguments passed on to \code{\link{predict.FDboost}}
#' 
#' @seealso \code{\link{FDboost}} for the model fit.
#' @return matrix of fitted values
#' @method fitted FDboost
#' @export
### similar to fitted.mboost() but returns the fitted values as matrix
fitted.FDboost <- function(object, ...) {
  args <- list(...)
  if (length(args) == 0) {
    if(is.null(object$id)){
      ret <- matrix(object$fitted(), nrow=object$ydim[1])
    }else{
      ret <- object$fitted()
    }
  } else {
    ret <- predict(object, newdata=NULL, ...)
  }
  ret
}
 
#' Residual values of a boosted functional regression model 
#' 
#' Takes a fitted \code{FDboost}-object and computes the residuals.
#' 
#' @param object a fitted \code{FDboost}-object
#' @param ... not used
#' 
#' @details The residual is missing if the corresponding value of the response was missing.
#' @seealso \code{\link{FDboost}} for the model fit.
#' @return matrix of residual values
#' @method residuals FDboost
#' @export
### residuals (the current negative gradient)
residuals.FDboost <- function(object, ...){
  
  if(is.null(object$id)){
    resid <- matrix(object$resid(), nrow=object$ydim[1])
    resid[is.na(object$response)] <- NA 
  }else{
    resid <- object$resid()
    resid[is.na(object$response)] <- NA 
  }
  resid
}


#' Coefficients of boosted functional regression model 
#' 
#' Takes a fitted \code{FDboost}-object produced by \code{\link{FDboost}()} and 
#' returns estimated coefficient functions/surfaces \eqn{\beta(t), \beta(s,t)} and 
#' estimated smooth effects \eqn{f(z), f(x,z)} or \eqn{f(x, z, t)}. 
#' Not implemented for smooths in more than 3 dimensions.
#' 
#' @param object a fitted \code{FDboost}-object
#' @param raw  logical defaults to FALSE.
#' If raw=FALSE for each effect the estimated function/surface is calculated
#' If raw=TRUE the coefficients of the model are returned. 
#' @param which a subset of base-learners for which the coefficients
#' should be computed (numeric vector), 
#' defaults to NULL which is the same as \code{which=1:length(object$baselearner)}
#' In the special case of \code{which=0}, only the coefficients of the offset are returned.
#' @param computeCoef defaults to TRUE, if FALSE only the names of the terms are returned
#' @param n1 see below
#' @param n2 see below
#' @param n3 n1, n2, n3 give the number of grid-points for 1-/2-/3-dimensional 
#' smooth terms used in the marginal equidistant grids over the range of the 
#' covariates at which the estimated effects are evaluated.
#' @param n4 gives the number of points for the third dimension in a 3-dimensional smooth term
#' @param ... other arguments, not used.
#' 
#' @return If \code{raw==FALSE}, a list containing 
#' \itemize{
#'  \item \code{pterms} a matrix containing the parametric / non-functional coefficients 
#'  \item \code{smterms} a named list with one entry for each smooth term in the model. 
#'  Each entry contains
#'     \itemize{
#'          \item \code{x, y, z} the unique grid-points used to evaluate the smooth/coefficient function/coefficient surface
#'          \item \code{xlim, ylim, zlim} the extent of the x/y/z-axes
#'          \item \code{xlab, ylab, zlab} the names of the covariates for the x/y/z-axes
#'          \item \code{value} a vector/matrix/list of matrices containing the coefficient values 
#'          \item \code{dim} the dimensionality of the effect
#'          \item \code{main} the label of the smooth term (a short label)
#' }} 
#' @method coef FDboost
#' @export
### similar to coef.pffr() by Fabian Scheipl in package refund 
coef.FDboost <- function(object, raw=FALSE, which=NULL, computeCoef=TRUE, 
                         n1=40, n2=40, n3=20, n4=10,...){
  
  if(raw){
    return(object$coefficients)  
  } else {
    
    # delete an extra 0 in which as the offset is always returned
    if( length(which) > 1 && 0 %in% which) which <- which[which!=0]
    
    # List to be returned
    ret <- list()
    
    ## <FIXME> all linear terms should be part of pterms?
    # ret$pterms <- NULL 
    
    ## offset as first element
    ret$offset$x <- seq( min(object$yind), max(object$yind), l=n1)
    ret$offset$xlab <- attr(object$yind, "nameyind")
    ret$offset$xlim <- range(object$yind)
    ret$offset$value <- object$predictOffset(ret$offset$x)
    if(length(ret$offset$value)==1) ret$offset$value <- rep(ret$offset$value, n1)
    ret$offset$dim <- 1
    ret$offset$main <- "offset"
    
    # For the special case of which=0, only return the coefficients of the offset
    if(!is.null(which) & length(which)==1 && which==0){
      if(computeCoef){
        return(ret)
      }else{
        return("offset")
      }  
    } 
    
    getCoefs <- function(i){
      ## this constructs a grid over the range of the covariates
      ## and returns estimated values on this grid, with 
      ## by-variables set to 1
      ## cf. mgcv:::plots.R (plot.mgcv.smooth etc..) for original code
      
      safeRange <- function(x){
        if(is.factor(x)) return(c(NA, NA))
        return(range(x, na.rm=TRUE))
      }
      
      makeDataGrid <- function(trm){        
        ### <TODO>  delete attr(d, "xm") <- xg and change code accordingly
        varnms <- trm$get_names()
        yListPlace <- NULL
        zListPlace <- NULL
        
        #generate grid of values in range of original data
        if(trm$dim==1){
          ng <- n1
          varnms <- varnms[!varnms %in% c("ONEx", "ONEtime")] 
          # Extra setup of dataframe in the case of a functional covariate
          if(grepl("bsignal", trm$get_call()) | grepl("bconcurrent", trm$get_call())){
            x <- attr(trm$model.frame()[[1]], "signalIndex")
            xg <- seq(min(x), max(x),length=ng) 
            varnms[1] <- attr(trm$model.frame()[[1]], "indname")            
          }else{
            x <- trm$model.frame()[[varnms]]
            xg <- if(is.factor(x)) {
              sort(unique(x))
            } else seq(min(x), max(x), length=ng)
          }
          d <- list(xg)  # data.fame
          names(d) <- varnms
          attr(d, "xm") <- xg
          attr(d, "xname") <- varnms
          # For effect constant over index of response: add dummy-index so that length in clear
          if(attr(object$yind, "nameyind") != varnms){
            d[[attr(object$yind, "nameyind")]] <- seq(min(object$yind), max(object$yind),length=ng) 
          }           
        }        
        
        if(trm$dim > 1){          
          ng <- ifelse(trm$dim==2, n2, n3)           
          
          ### get variables for x, y and eventually z direction
          
          # Extra setup of dataframe in the case of a functional covariate
          if(grepl("bsignal", trm$get_call()) | grepl("bhist", trm$get_call())){
            x <- attr(trm$model.frame()[[1]], "signalIndex")
            xg <- seq(min(x), max(x),length=ng) 
            varnms[1] <- attr(trm$model.frame()[[1]], "indname")             
          }else{ # scalar covariate
            x <- trm$model.frame()[[1]]
            xg <- if(is.factor(x)) {
              sort(unique(x))
            } else seq(min(x), max(x),length=ng)            
          }
          yListPlace <- ifelse(grepl("by", trm$get_call()), 3, 2)
          if(!grepl("bhist", trm$get_call())){
            y <- trm$model.frame()[[yListPlace]] # index of response or second scalar covariate
          }else{
            y <- object$yind
            varnms[2] <- attr(object$yind, "nameyind")
            if(varnms[1]==varnms[2]) varnms[1] <- paste(varnms[2], "_cov", sep="")
          }
          yg <- if(is.factor(y)) {
            sort(unique(y))
          } else seq(min(y), max(y),length=ng)
          if(length(varnms)==2){
            d <- list(xg, yg)  # data.frame
            attr(d, "xm") <- xg
            attr(d, "ym") <- yg    
          } else {
            zListPlace <- ifelse(grepl("by", trm$get_call()), 2, 3)
            z <- trm$model.frame()[[zListPlace]]
            zg <- if(is.factor(z)) {
              sort(unique(z))
            }else{
              if(grepl("by", trm$get_call())){ 1 }else{ seq(min(z), max(z), length=n4) }
            } 
            d <- list(xg, yg, zg)  # data.frame
            ## special case of factor by-variable 
            #if(grepl("by", trm$get_call()) && grepl("bols", trm$get_call()) && is.factor(z)){
            #  d <- list(rep(xg, length(yg)), rep(yg, each=length(xg)), zg)
            #}
            # special case of factor by-variable 
            if(grepl("by", trm$get_call()) && grepl("bols", trm$get_call()) && is.factor(z)){
              d <- list(rep(xg, length(zg)), yg, rep(zg, each=length(zg)))
            }
            #<FIXME> what to do if %X% was used???
            if(grepl("%X%", trm$get_call())){
              #d <- list( rep(xg), rep(yg, each=length(xg)), zg)
              d <- list(xg, yg=1, zg)
            } 
            attr(d, "xm") <- d[[1]]
            attr(d, "ym") <- d[[2]]
            attr(d, "zm") <- d[[3]]
          }
          names(d) <- varnms[c(1, yListPlace, zListPlace)] # colnames(d) <- varnms
          attr(d, "varnms") <- varnms[c(1, yListPlace, zListPlace)] 
        }
        
        
        ## add dummy signal to data for bsignal()
        if(grepl("bsignal", trm$get_call()) ){
          d[[ trm$get_names()[1] ]] <- I(diag(ng)/integrationWeights(diag(ng), d[[varnms[1]]] ))
        }
        
        ## add dummy signal to data for bhist()
        ### <FIXME> is dummy-signal for bhist() correct?
        ## standardisation weights depending on t must be multiplied to the final \beta(s,t)
        ## as they cannot be included into the variable x(s)
        ## use intFun() to compute the integration weights
        # ls(environment(trm$dpp))
        if(grepl("bhist", trm$get_call()) ){
          ## temp <- I(diag(ng)/integrationWeightsLeft(diag(ng), d[[varnms[1]]]))
          ## use intFun() of the bl to compute the integration weights
          temp <- environment(trm$dpp)$args$intFun(diag(ng), d[[varnms[1]]])
          ##if(environment(trm$dpp)$args$stand=="transform"){
          ##  xindStand <- (d[[varnms[1]]] - min(d[[varnms[1]]])) / (max(d[[varnms[1]]]) - min(d[[varnms[1]]]))
          ##  temp <- environment(trm$dpp)$args$intFun(diag(ng), xindStand)
          ##}
          d[[attr(trm$model.frame()[[1]], "xname")]] <- I(diag(ng)/temp)
        }
        
        ## add dummy signal to data for bconcurrent()
        if(grepl("bconcurrent", trm$get_call())){
          d[[ trm$get_names()[1] ]] <- I(matrix(rep(1.0, ng^2), ncol=ng))
        }
        
        if(trm$get_vary() != ""){
          d$by <- 1
          colnames(d) <- c(head(colnames(d),-1), trm$get_vary())
        } 
        
        # set time-variable to 1, if response is a scalar
        if(length(object$yind)==1){
          if( all(attr(d, "zm") == d[[attr(object$yind, "nameyind")]]) ){
            attr(d, "zm") <- 1
          }
          d[[attr(object$yind, "nameyind")]] <- 1 
        }
        
        return(d)
      }
      
      getP <- function(trm, d){
        #return an object similar to what plot.mgcv.smooth etc. returns 
        if(trm$dim==1){
          predHelp <- predict(object, which=i, newdata=d)
          if(!is.matrix(predHelp)){ X <- predHelp
          }else{
            X <- if(any(trm$get_names() %in% c("ONEtime"))){ # effect constant in t 
              predHelp[,1]
            }else{ predHelp[1,] # smooth intercept/ concurrent effect                
            }            
          } 
          P <- list(x=attr(d, "xm"), xlab=attr(d, "xname"), xlim=safeRange(attr(d, "xm")))
        ## trm$dim > 1
        }else{
          varnms <- attr(d, "varnms")
          if(trm$dim==2){
            X <- predict(object, newdata=d, which=i)
            ## for bhist(), multiply with standardisation weights if necessary
            ## you need the args$vecStand from the prediction of X, constructed here
            if(grepl("bhist", trm$get_call()) && 
                 environment(trm$dpp)$args$stand %in% c("length","time")){  
              Lnew <- environment(trm$dpp)$args$intFun(d[[3]], d[[varnms[1]]])
              ## Standardize with exact length of integration interval
              ##  (1/t-t0) \int_{t0}^t f(s) ds
              if(environment(trm$dpp)$args$stand == "length"){
                ind0 <- !t(outer( d[[varnms[1]]], d[[varnms[2]]], environment(trm$dpp)$args$limits) )
                Lnew[ind0] <- 0
                ## integration weights in s-direction always sum exactly to 1, 
                vecStand <- rowSums(Lnew)
              }           
              ## use time of current observation for standardization
              ##  (1/t) \int_{t0}^t f(s) ds
              if(environment(trm$dpp)$args$stand=="time"){
                yindHelp <- d[[varnms[1]]]
                yindHelp[yindHelp==0] <- Lnew[1,1] 
                vecStand <- yindHelp
              }
              X <- t(t(X)*vecStand) 
            }
            P <- list(x=attr(d, "xm"), y=attr(d, "ym"), xlab=varnms[1], ylab=varnms[2],
                      ylim=safeRange(attr(d, "ym")), xlim=safeRange(attr(d, "xm")),
                      z=attr(d, "zm"), zlab=varnms[3])
          }else{
            if(trm$dim==3){
              values3 <- seq(min(d[[varnms[3]]]), max(d[[varnms[3]]]), l=n4)
              xygrid <- expand.grid(d[[varnms[1]]], d[[varnms[2]]] )
              X <- lapply(values3, function(x){
                d1 <- list(xygrid[,1], xygrid[,2], x)
                names(d1) <- varnms
                matrix(predict(object, newdata=d1, which=i), ncol=length(d[[varnms[1]]]))
              }) 
              P <- list(x=attr(d, "xm"), y=attr(d, "ym"), z=values3, 
                   xlab=varnms[1], ylab=varnms[2], zlab=varnms[3],
                   ylim=safeRange(attr(d, "ym")), xlim=safeRange(attr(d, "xm")), zlim=safeRange(attr(d, "zm")))
            }
          }
        }
        if(!is.null(object$ydim)){
          P$value <- X
        }else{
          #### <FIXME> is dimension, byrow correct??
          P$value <- matrix(X, nrow=length(attr(d, "xm")) )
        }
        #P$coef <- cbind(d, "value"=P$value)    
        P$dim <- trm$dim
        return(P)
      }
      
      trm <- object$baselearner[[i]]
      trm$dim <- length(trm$get_names())
      if(any(grepl("ONEx", trm$get_names()), 
             grepl("ONEtime", trm$get_names()))) trm$dim <- trm$dim - 1 
            
      if( grepl("bhist", trm$get_call()) ){
        trm$dim <- 2
      }
      
      # If a by-variable was specified, reduce number of dimensions
      # as smooth linear effect in several groups can be plotted in one plot 
      if( grepl("by =", trm$get_call()) && grepl("bols", trm$get_call()) || 
            grepl("by =", trm$get_call()) && grepl("bbs", trm$get_call()) ) trm$dim <- trm$dim - 1
      
      # <FIXME> what to do with bbs(..., by=factor)?
      
      if(trm$dim > 3){
        warning("can't deal with smooths with more than 3 dimensions, returning NULL for ", 
                shrtlbls[i])
        return(NULL)
      }
      
      d <- makeDataGrid(trm) 
      ### <FIXME> better solution for %X% in base-learner!!!
      if(!is.null(object$ydim) & any(grepl("%X%", trm$get_call()) ) ) trm$dim <- trm$dim - 1
      
      ## it is necessary to expand the dataframe!
      if(is.null(object$ydim) && !grepl("bconcurrent", trm$get_call())){
        #browser()
        #print(attr(d, "varnms"))
        vari <- names(d)[1]
        if(is.factor(d[[vari]])){
          d[[vari]] <- d[[vari]][ rep(1:NROW(d[[vari]]), times=length(d[[attr(object$yind ,"nameyind")]]) ) ]
          if(trm$dim>1) d[[attr(object$yind ,"nameyind")]] <- rep(d[[attr(object$yind ,"nameyind")]], 
                                                                  each=length(unique(d[[vari]]))  )
          }else{
          # expand signal variable
          if( grepl("bhist", trm$get_call()) | grepl("bsignal", trm$get_call()) ){
            vari <- names(d)[!names(d) %in% attr(d, "varnms")]
            d[[vari]] <- d[[vari]][ rep(1:NROW(d[[vari]]), times=NROW(d[[vari]])), ]
            
          }else{ # expand scalar variable
            vari <- names(d)[1]
            if(vari!=attr(object$yind ,"nameyind")) d[[vari]] <- d[[vari]][ rep(1:NROW(d[[vari]]), times=NROW(d[[vari]])) ]
          } 
          # expand yind 
          if(trm$dim>1) d[[attr(object$yind ,"nameyind")]] <- rep(d[[attr(object$yind ,"nameyind")]], 
                                                                  each=length(d[[attr(object$yind ,"nameyind")]])) 
        }
      }

      P <- getP(trm, d)  

      # get proper labeling
      P$main <- shrtlbls[i]
      
      return(P)
    } # end of function getCoefs()
      

    ## Function to obtain nice short names for labeling of plots
    shortnames <- function(x){
      if(substr(x,1,1)=="\"") x <- substr(x, 2, nchar(x)-1)
      
      sign <- "%O%"
      if( !grepl(sign, x) ) sign <- "%X%"
      
      xpart <- unlist(strsplit(x, sign)) 
      
      for(i in 1:length(xpart)){
        xpart[i] <- gsub("\\\"", "'", xpart[i], fixed=TRUE)
        nvar <- length(all.vars(formula(paste("Y~", xpart[i])))[-1])
        commaSep <- unlist(strsplit(xpart[i], ","))  
        
        # shorten the name to first variable and delte x= if present
        if(grepl("=", commaSep[1])){
          temp <- unlist(strsplit(commaSep[1], "="))
          temp[1] <- unlist(strsplit(temp[1], "(", fixed=TRUE))[1]
          if(substr(temp[2], 1, 1)==" ") temp[2] <- substr(temp[2], 2, nchar(temp[2]))
          xpart[i] <- paste(temp[1], "(", temp[2], ")", sep="")
        }else{
          xpart[i] <- paste(commaSep[1], ")", sep="")
        }
        #xpart[i] <- if(length(commaSep)==1){
        #  paste(paste(commaSep[1:nvar], collapse=","), sep="")
        #}else paste(paste(commaSep[1:nvar], collapse=","), ")", sep="") 
        #if(substr(xpart[i], nchar(xpart[i])-2, nchar(xpart[i])-1)==") "){
        #  xpart[i] <- substr(xpart[i], 1, nchar(xpart[i])-2)
        #}
      }
      ret <- xpart
      if(length(xpart)>1) ret <- paste(xpart, collapse=paste("", sign))
      ret
    }
    
    if(is.null(which)) which <- 1:length(object$baselearner)
    
    ## short names for the terms, if shortnames() does not work, use the original names
    shrtlbls <- try(unlist(lapply(names(object$baselearner), shortnames)))
    if(class(shrtlbls)=="try-error") shrtlbls <- names(object$baselearner) 
    
    if(computeCoef){
      ### smooth terms
      ret$smterms <- lapply(which, getCoefs)         
      names(ret$smterms) <- sapply(seq_along(ret$smterms), function(i){
        ret$smterms[[i]]$main
      })
      return(ret)
    }else{
      return(shrtlbls[which])
    }  
  }
}

# help function to color perspective plots - col1 positive values, col2 negative values
getColPersp <- function(z, col1="tomato", col2="lightblue"){
  nrz <- nrow(z)
  ncz <- ncol(z)
  
  # Compute the z-value at the facet centres
  zfacet <- z[-1, -1] + z[-1, -ncz] + z[-nrz, -1] + z[-nrz, -ncz]
  
  # use the colors col1 and col2 for negative and positive values
  colfacet <- matrix(nrow=nrow(zfacet), ncol=ncol(zfacet))
  colfacet[zfacet < 0] <- col2
  colfacet[zfacet > 0] <- col1
  colfacet[zfacet == 0] <- "white"
  
  return(colfacet) 
}


#' Plot the fit or the coefficients of a boosted functional regression model 
#' 
#' Takes a fitted \code{FDboost}-object produced by \code{\link{FDboost}()} and 
#' plots the fitted effects or the coefficient-functions/surfaces.
#' 
#' @param x a fitted \code{FDboost}-object
#' @param raw  logical defaults to FALSE.
#' If raw=FALSE for each effect the estimated function/surface is calculated
#' If raw=TRUE the coefficients of the model are returned. 
#' @param rug when TRUE (default) then the covariate to which the plot applies is 
#' displayed as a rug plot at the foot of each plot of a 1-d smooth, 
#' and the locations of the covariates are plotted as points on the contour plot 
#' representing a 2-d smooth.
#' @param which a subset of base-learners to take into account for plotting. 
#' a list is returned. 
#' @param includeOffset logical, defaults to TRUE. Should the offset be included in the plot of the intercept (default)
#' or should it be plotted separately.
#' @param ask logical, defaults to TRUE, if several effects are plotted the user
#' has to hit Return to see next plot.
#' @param n1 see below
#' @param n2 see below
#' @param n3 n1, n2, n3 give the number of grid-points for 1-/2-/3-dimensional 
#' smooth terms used in the marginal equidistant grids over the range of the 
#' covariates at which the estimated effects are evaluated.
#' @param n4 gives the number of points for the third dimension in a 3-dimensional smooth term
#' @param onlySelected, logical, defaults to TRUE. Only plot effects that were 
#' selected in at least one boosting iteration
#' @param pers logical, defaults to FALSE, 
#' If TRUE, perspective plots (\code{\link[graphics]{persp}}) for 
#' 2- and 3-dimensional effects are drawn.
#' If FALSE, image/contour-plots (\code{\link[graphics]{image}}, 
#' \code{\link[graphics]{contour}}) are drawn for 2- and 3-dimensional effects. 
#' @param commonRange logical, defaults to FALSE, 
#' if TRUE the range over all effects is the same (so far not implemented)
#' 
#' @param subset subset of the observed response curves and their predictions that is plotted. 
#' Per default all observations are plotted.
#' @param posLegend location of the legend, if a legend is drawn automatically 
#' (only used in plotPredicted). The default is "topleft".
#' @param lwdObs lwd of observed curves (only used in plotPredicted)
#' @param lwdPred lwd of predicted curves (only used in plotPredicted)
#' @param ... other arguments, passed to \code{funplot} (only used in plotPredicted)
#' 
#' @aliases plotPredicted plotResiduals
#' 
#' @seealso \code{\link{FDboost}} for the model fit and 
#' \code{\link{coef.FDboost}} for the calculation of the coefficient functions. 
#' 
#' @method plot FDboost
#' @export
### function to plot raw values or coefficient-functions/surfaces of a model 
plot.FDboost <- function(x, raw=FALSE, rug=TRUE, which=NULL, 
                         includeOffset=TRUE, ask=TRUE,
                         n1=40, n2=40, n3=20, n4=11,
                         onlySelected=TRUE, pers=FALSE, commonRange=FALSE, ...){
  
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

  argsImage <- getArguments(x=formals(graphics::image.default), dots=dots)
  argsContour <- getArguments(x=formals(graphics::contour.default), dots=dots)
  argsPersp <- getArguments(x=formals(getS3method("persp", "default")), dots=dots)
  
  plotWithArgs <- function(plotFun, args, myargs){        
    args <- c(myargs[!names(myargs) %in% names(args)], args)        
    do.call(plotFun, args)            
  }
  
  #     #### function by Fabian Scheipl
  #     # colors in rgb
  #     alpha <- function(x, alpha=25){
  #       tmp <- sapply(x, col2rgb)
  #       tmp <- rbind(tmp, rep(alpha, length(x)))/255
  #       return(apply(tmp, 2, function(x) do.call(rgb, as.list(x))))
  #     }      
  #     clrs <- alpha( rainbow(x$ydim[1]), 125) 
  
  ### get the effects to be plotted
  whichSpecified <- which
  if(is.null(which)) which <- 1:length(x$baselearner) 
  
  if(onlySelected){
#     sel <- selected(x)
#     if( !1 %in% sel  && length(x$offsetVec) > 1 && grepl("ONEx", names(x$baselearner)[[1]])){
#       sel <- c(1, sel) # plot the offset as first effect
#     } 
    which <- intersect(which, c(0, selected(x)))
  }
  
  # In the case that intercept and offset should be plotted and the intercept was never selected
  # plot the offset
  if( (1 %in% whichSpecified | is.null(whichSpecified))  & !1 %in% which & length(x$yind)>1) which <- c(0, which)
  
  if(length(which)==0){
    warning("Nothing selected for plotting.")
    return(NULL)
  } 

  ### plot coefficients of model (smooth curves and surfaces)
  if(!raw){ 
        
    # compute the coefficients of the smooth terms that should be plotted
    coefMod <- coef(x, which=which, n1=n1, n2=n2, n3=n3, n4=n4)
    terms <- coefMod$smterms
    offsetTerms <- coefMod$offset
    bl_data <- lapply(x$baselearner[which], function(x) x[["get_data"]]()) 
    
    # plot nothing but the offset
    if(length(which)==1 && which==0){
      terms <- list(offsetTerms)
      bl_data <- c(offset = list( list(x$yind) ), bl_data)
      names(bl_data[[1]]) <- attr(x$yind, "nameyind")
    }
    
    # include the offset in the plot of the intercept
    if( includeOffset && 1 %in% which && grepl("ONEx", names(terms)[1]) ){
      terms[[1]]$value <- terms[[1]]$value + matrix(offsetTerms$value, ncol=1, nrow=n1)
      terms[[1]]$main <- paste("offset", "+", terms[[1]]$main)
    }
    
    # plot the offset as extra effect
    # case 1: the offset should be included as extra plot
    # case 2: the whole model is plotted, but the intercept-base-learner was never selected
    if( (!includeOffset | (includeOffset & !1 %in% which)) & 
         is.null(whichSpecified)){
      terms <- c(offset = list(offsetTerms), terms)
      bl_data <- c(offset = list( list(x$yind) ), bl_data)
      names(bl_data[[1]]) <- attr(x$yind, "nameyind")     
    } 
   
    if((length(terms)>1 | terms[[1]]$dim==3) & ask) par(ask=TRUE)
    
#     ### <TODO> implement common range
#     if(commonRange){
#       #range <- range(fit)
#       #range[1] <- range[1]-0.02*diff(range)
#     }else range <- NULL
    
    for(i in 1:length(terms)){
      trm <- terms[[i]] 
      
      if(grepl("bhist", trm$main)){
        ## set 0 to NA so that beta only has values in its domain
        # get the limits- function
        limits <- get("args", (environment(x$baselearner[[which[i]]]$dpp)))$limits
        ## do not use it if stand  == "transform"
        ## if(get("args", (environment(x$baselearner[[which[i]]]$dpp)))$stand != "transform"){
        trm$value[!outer(trm$x, trm$y, limits)] <- NA
        # }        
      }
      
      # plot for 1-dim effects
      if(trm$dim==1){
        if(length(trm$value)==1) trm$value <- rep(trm$value, l=length(trm$x)) 
        
        if(!"add" %in% names(dots)){
          plotWithArgs(plot, args=argsPlot, 
                       myargs=list(x=trm$x, y=trm$value, xlab=trm$xlab, main=trm$main, 
                                   ylab="coef", type="l"))
        }else{
          plotWithArgs(lines, args=argsPlot, 
                       myargs=list(x=trm$x, y=trm$value, xlab=trm$xlab, main=trm$main, 
                                   ylab="coef", type="l"))          
        }
        
          if(rug & !is.factor(x=trm$x)){
            if(grepl("bconcurrent", trm$main) | grepl("bsignal", trm$main)){
              rug(attr(bl_data[[i]][[1]], "signalIndex"), ticksize = 0.02)
            }else ifelse(length(unique(bl_data[[i]][[1]]))!=1,
                         rug(bl_data[[i]][[1]], ticksize = 0.02),
                         rug(bl_data[[i]][[2]], ticksize = 0.02))
          } 
      } 
      
      # plot with factor variable
      if(trm$dim==2 & (is.factor(trm$x) | is.factor(trm$y)) | is.factor(trm$z) ){
        if(is.factor(trm$y)){ # effect with by-variable (by-variable is factor)
          plotWithArgs(matplot, args=argsMatplot, 
                       myargs=list( x=trm$z, y=t(trm$value), xlab=trm$ylab, main=trm$main, 
                                    ylab="coef", type="l", col=as.numeric(trm$y) ) )
          if(rug){
            #rug(bl_data[[i]][[3]], ticksize = 0.02) 
            rug(x$yind, ticksize = 0.02)
          }
        }else{ # effect of factor variable
          plotWithArgs(matplot, args=argsMatplot, 
                       myargs=list(x=trm$y, y=t(trm$value), xlab=trm$ylab, main=trm$main, 
                                   ylab="coef", type="l"))
          if(rug){
            #rug(bl_data[[i]][[2]], ticksize = 0.02) 
            rug(x$yind, ticksize = 0.02)
          }
        }
        
      }else{
        # persp-plot for 2-dim effects
        if(trm$dim==2 & pers){
          if(length(unique(as.vector(trm$value)))==1){
            # persp() gives error if only a flat plane should be drawn
            plot(y=trm$value[1,], x=trm$x, main=trm$main, type="l", xlab=trm$ylab, 
                 ylab="coef")
          }else{  
            range <- range(trm$value, na.rm = TRUE)
            if(range[1]==range[2]) range <- range(0, range)
            zlim <- c(range[1] - 0.05*(range[2] - range[1]), 
                      range[2] + 0.05*(range[2] - range[1]))
            plotWithArgs(persp, args=argsPersp,
                         myargs=list(x=trm$x, y=trm$y, z=trm$value, xlab=paste("\n", trm$xlab), 
                                     ylab=paste("\n", trm$ylab), zlab=paste("\n", "coef"), 
                                     main=trm$main, theta=30, zlim=zlim,
                                     phi=30, ticktype="detailed",
                                     col=getColPersp(trm$value)))
          } 
        }
        # image for 2-dim effects
        if(trm$dim==2 & !pers){        
          plotWithArgs(image, args=argsImage,
                       myargs=list(x=trm$y, y=trm$x, z=t(trm$value), xlab=trm$ylab, ylab=trm$xlab, 
                                   main=trm$main, col = terrain.colors(length(trm$x)^2)))          
          plotWithArgs(contour, args=argsContour,
                       myargs=list(trm$y, trm$x, z=t(trm$value), add = TRUE))
          
          if(rug){
            ##points(expand.grid(bl_data[[i]][[1]], bl_data[[i]][[2]]))
            if(grepl("bhist", trm$main)){
              rug(x$yind, ticksize = 0.02)
            }else{
              ifelse(grepl("by", trm$main),
                     rug(bl_data[[i]][[3]], ticksize = 0.02),
                     rug(bl_data[[i]][[2]], ticksize = 0.02))
            }
            ifelse(grepl("bsignal", trm$main) | grepl("bhist", trm$main),
              rug(attr(bl_data[[i]][[1]], "signalIndex"), ticksize = 0.02, side=2),
              rug(bl_data[[i]][[1]], ticksize = 0.02, side=2))
          }
        }        
      }
      ### 3 dim plots
      # persp-plot for 3-dim effects
      if(trm$dim==3 & pers){
        for(j in 1:length(trm$z)){
          plotWithArgs(persp, args=argsPersp,
                       myargs=list(x=trm$x, y=trm$y, z=trm$value[[j]], xlab=paste("\n", trm$xlab), 
                                   ylab=paste("\n", trm$ylab), zlab=paste("\n", "coef"), 
                                   theta=30, phi=30, ticktype="detailed", 
                                   zlim=range(trm$value), col=getColPersp(trm$value[[j]]), 
                                   main= paste(trm$zlab ,"=", round(trm$z[j],2), ": ", trm$main, sep=""))
            )         
        }
      }
      # image for 3-dim effects
      if(trm$dim==3 & !pers){
        for(j in 1:length(trm$z)){
          plotWithArgs(image, args=argsImage,
            myargs=list(x=trm$x, y=trm$y, z=trm$value[[j]], xlab=trm$xlab, ylab=trm$ylab,
                        col = terrain.colors(length(trm$x)^2), zlim=range(trm$value),
                        main= paste(trm$zlab ,"=", round(trm$z[j],2), ": ", trm$main, sep="")))
          plotWithArgs(contour, args=argsContour,
                       myargs=list(trm$x, trm$y, trm$value[[j]], xlab=trm$xlab, add = TRUE))
          if(rug){
            points(bl_data[[i]][[1]], bl_data[[i]][[2]]) 
          }
#           plotWithArgs(filled.contour, args=list(),
#                        myargs=list(x=trm$x, y=trm$y, z=trm$value[[j]], xlab=trm$xlab, ylab=trm$ylab,
#                                    zlim=range(trm$value), color.palette=terrain.colors,
#                                    main= paste(trm$zlab ,"=", round(trm$z[j],2), ": ", trm$main, sep="")))
        }        
      }
    } # end for-loop
    if(length(terms)>1 & ask) par(ask=FALSE) 
    
  ### plot smooth effects as they are estimated for the original data
  }else{
        
    ################################          
    # predict the effects using the original data
    terms <- predict(x, which=which)
    
    # convert matrix into a list, each list entry for one effect
    if(is.null(x$ydim) & !is.null(dim(terms))){
      temp <- list()
      for(i in 1:ncol(terms)){
        temp[[i]] <- terms[,i]
      }
      names(temp) <- colnames(terms)
      terms <- temp
      rm(temp)
    }
    
    if(length(which)==1 && which==0) terms <- attr(terms, "offset")
    
    if(length(which)==1 & is.null(x$id)) terms <- list(terms)
    if(class(terms)!="list") terms <- list(terms) 
    
    shrtlbls <- try(coef(x, which=which, computeCoef=FALSE))# get short names
    if(class(shrtlbls)=="try-error"){
      shrtlbls <- names(x$baselearner)[which[which!=0]]
      if(0 %in% which) shrtlbls <- c("offset", which)
    }
    if(is.null(shrtlbls)) shrtlbls <- "offset" 
    time <- x$yind
    
    # include the offset in the plot of the intercept
    if( includeOffset && 1 %in% which && grepl("ONEx", shrtlbls[1]) ){
      terms[[1]] <- terms[[1]] + x$offset
      shrtlbls[1] <- paste("offset", "+", shrtlbls[1])
    }   
    if(length(which)>1 & ask) par(ask=TRUE)
    if(commonRange){
      range <- range(terms)
      range[1] <- range[1]-0.02*diff(range)
    }else range <- NULL
    for(i in 1:length(terms)){
      # set values of predicted effect to missing if response is missing
      if(sum(is.na(x$response))>0) terms[[i]][is.na(x$response)] <- NA
      if(is.null(dots$ylim)){
        range <- range(terms[[i]], na.rm=TRUE)
      }else{
        range <- dots$ylim
      }
      
      if(length(time)>1){
        
        plotWithArgs(funplot, args=argsFunplot, 
                     myargs=list(x=time, y=terms[[i]], id=x$id, type="l", ylab="effect", lty=1, rug=FALSE,
                                 xlab=attr(time, "nameyind"), ylim=range, main=shrtlbls[i]))
        if(rug) rug(time)
      }else{
        plotWithArgs(plot, args=argsPlot, 
                     myargs=list(x=x$response-x$offset, y=terms[[i]], type="p", ylab="effect", 
                                 xlab="response-offset", ylim=range, main=shrtlbls[i])) 
      }
    }
    if(length(which)>1 & ask) par(ask=FALSE) 
    }       
}



########################################################################################

#' Extract information of a base-learner
#' 
#' Takes a base-learner and extracts information.
#' 
#' @param object a base-learner
#' @param what  a character specifying the quantities to extract.
#' This can be a subset of "design" (default; design matrix), 
#' "penalty" (penalty matrix) and "index" (index of ties used to expand 
#' the design matrix)
#' @param asmatrix a logical indicating whether the the returned matrix should be 
#' coerced to a matrix (default) or if the returned object stays as it is 
#' (i.e., potentially a sparse matrix). This option is only applicable if extract 
#' returns matrices, i.e., what = "design" or what = "penalty".
#' @param expand a logical indicating whether the design matrix should be expanded 
#' (default: FALSE). This is useful if ties where taken into account either manually 
#' (via argument index in a base-learner) or automatically for data sets with many 
#' observations. expand = TRUE is equivalent to extract(B)[extract(B, what = "index"),] 
#' for a base-learner B.
#' @param ... currently not used
#' @method extract blg
#' @seealso \code{\link[mboost]{extract}} for the extract functions of the package mboost
extract.blg <- function(object, what = c("design", "penalty", "index"),
                        asmatrix = FALSE, expand = FALSE, ...){
  what <- match.arg(what)
  
  if(grepl("%O%", object$get_call())){
    object <- object$dpp( rep(1, NROW(object$model.frame()[[1]])) )    
  }else{
    object <- object$dpp(rep(1, nrow(object$model.frame())))
  }  
  return(extract(object, what = what,
                 asmatrix = asmatrix, expand = expand))
}
