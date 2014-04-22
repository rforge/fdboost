
### models for paper:

## full
## (d)spec h2o
## (d)spec
## spec
## NIR

### use original observations without preprocessing!
library(FDboost)

###### Load data and compute the first derivative
data(fuel)
str(fuel)

# # normalize the wavelength to 0-1
# fuel$nir.lambda0 <- (fuel$nir.lambda - min(fuel$nir.lambda)) / 
#   (max(fuel$nir.lambda) - min(fuel$nir.lambda)) 
# fuel$uvvis.lambda0 <- (fuel$uvvis.lambda - min(fuel$uvvis.lambda)) / 
#   (max(fuel$uvvis.lambda) - min(fuel$uvvis.lambda))

# compute first derivatives as first order differences
fuel$dUVVIS <- t(apply(fuel$UVVIS, 1, diff))
fuel$dNIR <- t(apply(fuel$NIR, 1, diff)) 

# get the wavelength for the derivatives
fuel$duvvis.lambda <- fuel$uvvis.lambda[-1]
fuel$dnir.lambda <- fuel$nir.lambda[-1]
# fuel$duvvis.lambda0 <- fuel$uvvis.lambda0[-1]
# fuel$dnir.lambda0 <- fuel$nir.lambda0[-1]


################################################################################################
######### Model-fit for heatan with fun cov and derivatives


################################################################################################

if(FALSE){
  
  modH2O <- FDboost(h2o ~ bsignal(UVVIS, uvvis.lambda, knots=40, df=4) 
                    + bsignal(NIR, nir.lambda, knots=40, df=4)
                    + bsignal(dUVVIS, duvvis.lambda, knots=40, df=4) 
                    + bsignal(dNIR, dnir.lambda, knots=40, df=4), 
                    timeformula=~bols(1), data=fuel)
  
  set.seed(212)
  cvmH2O <- suppressWarnings(cvrisk(modH2O, grid=seq(100, 5000, by=100), 
                                    folds=cv( model.weights(modH2O), 
                                              type = "bootstrap", B = 10), mc.cores=10))
  
  par(mfrow=c(1,2))
  plot(cvmH2O)
  
  modH2O[mstop(cvmH2O)]
  #modH2O[2400]
  
  #### create new variable of predicted h2o
  h2o.fit <- modH2O$fitted()
  
  plot(fuel$h2o, h2o.fit)
  abline(0,1)  
}

################################################################################################
# fit all models with different effects
allTypes <- c("full", "dspech2o", "dspec", "spec", "NIR")

relmse <- list()
mse <- list()

for(i in 1:length(allTypes)){
  
  # i <- 1
  print("#########################################")
  print(allTypes[i])
  
  ## set up formula depending on allTypes[i]
  formula <- switch(allTypes[i],
                    full = heatan ~ bsignal(UVVIS, uvvis.lambda, knots=40, df=4.41) 
                    + bsignal(NIR, nir.lambda, knots=40, df=4.41)
                    + bsignal(dUVVIS, duvvis.lambda, knots=40, df=4.41) 
                    + bsignal(dNIR, dnir.lambda, knots=40, df=4.41)
                    + bbs(h2o.fit, knots=10, df=4.41)
                    + bsignal(UVVIS, uvvis.lambda, df=2.1) %X% bbs(h2o.fit, knots=10, df=2.1)
                    + bsignal(NIR, nir.lambda, df=2.1) %X% bbs(h2o.fit, knots=10, df=2.1)
                    + bsignal(dUVVIS, duvvis.lambda, df=2.1) %X% bbs(h2o.fit, knots=10, df=2.1)
                    + bsignal(dNIR, dnir.lambda, df=2.1) %X% bbs(h2o.fit, knots=10, df=2.1),
                    dspech2o = heatan ~ bsignal(UVVIS, uvvis.lambda, knots=40, df=4.41) 
                    + bsignal(NIR, nir.lambda, knots=40, df=4.41)
                    + bsignal(dUVVIS, duvvis.lambda, knots=40, df=4.41) 
                    + bsignal(dNIR, dnir.lambda, knots=40, df=4.41)
                    + bbs(h2o.fit, knots=10, df=4.41),
                    dspec = heatan ~ bsignal(UVVIS, uvvis.lambda, knots=40, df=4.41) 
                    + bsignal(NIR, nir.lambda, knots=40, df=4.41)
                    + bsignal(dUVVIS, duvvis.lambda, knots=40, df=4.41) 
                    + bsignal(dNIR, dnir.lambda, knots=40, df=4.41),
                    spec = heatan ~ bsignal(UVVIS, uvvis.lambda, knots=40, df=4.41) 
                    + bsignal(NIR, nir.lambda, knots=40, df=4.41),
                    NIR = heatan ~ bsignal(NIR, nir.lambda, knots=40, df=4.41)
  )
    
  ####################################################################################
  
  #### do a model fit with the variables selected in stabsel:
  mod <- FDboost(formula,
                 timeformula=~bols(1),
                 data=fuel) 
  # plot(mod, rug=FALSE)
  
  ####### cross-validation to determine variation in coefficients and predictions
  set.seed(2703)
  val <- validateFDboost(mod, folds=cv(model.weights(mod), type = "bootstrap", B = 50),
                         grid = 10:500, mc.cores=10)
  
  relmse[[i]] <- val$oobrelMSE
  mse[[i]] <- val$oobrisk # = val$oobmse # (as family=Gaussian)
  
  print(mstop(val))
  
#   save(val, file=paste(allTypes[i], "_val.RData", sep="") )
#   load("val.RData")
#  
#   mopt <- val$grid[which.min(colMeans(val$oobrisk))]
#   print(mopt)
#   mod[mopt]
#   print(summary(mod))
#   
#   colMeans(val$oobrisk)[paste(mopt)]
#   colMeans(val$relMSE)[paste(mopt)]
#   
  
  ##################
  
  ### plot predictions on whole dataset into graphic
  pdf(paste(allTypes[i], "_valCoef.pdf", sep=""))
  par(mar=c(5, 3, 1, 1), cex.axis=1.5, cex.lab=1.5)
  
  plotPredCoef(val, terms=FALSE, commonRange=TRUE, which=1, rug=FALSE)
  plot(mod, which=1, add=TRUE, lwd=2, col=4, lty=5, onlySelected=FALSE, rug=FALSE)
  
  if(length(mod$baselearner)>1){
    plotPredCoef(val, terms=FALSE, commonRange=TRUE, which=2, rug=FALSE)
    plot(mod, which=2, add=TRUE, lwd=2, col=4, lty=5, onlySelected=FALSE, rug=FALSE)
  }
  
  if(length(mod$baselearner)>2){
    plotPredCoef(val, terms=FALSE, commonRange=TRUE, which=3, rug=FALSE)
    plot(mod, which=3, add=TRUE, lwd=2, col=4, lty=5, onlySelected=FALSE, rug=FALSE)
  }
  
  if(length(mod$baselearner)>3){ 
    plotPredCoef(val, terms=FALSE, commonRange=TRUE, which=4, rug=FALSE)
    plot(mod, which=4, add=TRUE, lwd=2, col=4, lty=5, onlySelected=FALSE, rug=FALSE)
  }
  
  if(length(mod$baselearner)>4){
    plotPredCoef(val, terms=FALSE, commonRange=TRUE, which=5, rug=FALSE)
    plot(mod, which=5, add=TRUE, lwd=2, col=4, lty=5, onlySelected=FALSE, rug=FALSE)
  }  
  dev.off()  

  ###############
  if(allTypes[i]=="spec"){
    pdf("spec_valCoef2.pdf")
    par(mar=c(5, 3, 1, 1), cex.axis=1.5, cex.lab=1.5)
    plot(mod, which=1, lwd=2, col="white", main="", lty=5, rug=FALSE,
         ylab="", xlab="wavelength [nm]", ylim=range(val$coefCV[[1]]$value) )
    plotPredCoef(val, terms=FALSE, commonRange=TRUE, which=1, add=TRUE)
    plot(mod, which=1, lwd=2, col=4, main="", lty=5, rug=FALSE, add=TRUE)
    
    legend("topright", legend=c("data", "BS", "mean BS", "5, 95% BS"), 
           lty=c(5,1,1,2), col=c(4,8,1,2), cex=1.5,
           lwd=c(2,1,2,2))
    
    plot(mod, which=2, lwd=2, col="white", main="", lty=5, rug=FALSE, 
         ylab="", xlab="wavelength [nm]", ylim=range(val$coefCV[[2]]$value) )
    plotPredCoef(val, terms=FALSE, commonRange=TRUE, which=2, add=TRUE)
    plot(mod, which=2, lwd=2, col=4, main="", lty=5, rug=FALSE, add=TRUE)
    dev.off() 
  }
}

names(relmse) <- names(mse) <- allTypes

### get relMSE and MSE in optimal mstop
mopt <- lapply(mse, function(x) names(which.min(colMeans(x))) )

relmseopt <- lapply(1:length(allTypes), function(i) relmse[[i]][ ,mopt[[i]] ] )
mseopt <- lapply(1:length(allTypes), function(i) mse[[i]][ ,mopt[[i]] ] )

# give nice names
names(relmseopt) <- names(mseopt) <- c("full","(d)spec H2O", "(d)spec", "spec", "NIR")

pdf("relMSE.pdf", width=11, height=7)
par(mar=c(3, 5, 1, 1), cex.axis=1.8, cex.lab=1.8)
boxplot(relmseopt, ylab="relMSE")
boxplot(mseopt, ylab="MSE")
dev.off()





