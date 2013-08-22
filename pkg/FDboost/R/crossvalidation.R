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
