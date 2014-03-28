###############################################################################
# run the simulations for boosting of functional data
# code based on code by Fabian Scheipl 
# author: Sarah Brockhaus
###############################################################################

print(R.Version()$version.string)

library(refund)
library(FDboost)
pathResults <- NULL
pathModels <- NULL

source("simfuns.R")



###################################### M=50

set.seed(18102012)

settings <- makeSettings(
  M=c(50),
  ni=c(1),
  Gy=c(30),
  Gx=c(30),
  snrEps=c(1,2),
  snrE=c(0),
  snrB=c(2),
  scenario=c(3),
  balanced=c(TRUE),
  nuisance=c(0),
  rep=1)

length(settings)

usecores <- 1
options(cores=usecores)
M50N1G30 <- try(doSim(settings=settings, cores=usecores))

save(M50N1G30, file=paste(pathResults, "NEWplotModelsM50N1G30.Rdata", sep=""))

# 
# 
# set.seed(18102012)
# 
# settings <- makeSettings(
#   M=c(50),
#   ni=c(1),
#   Gy=c(100),
#   Gx=c(100),
#   snrEps=c(1,5),
#   snrE=c(0),
#   snrB=c(2),
#   scenario=c(1:5),
#   balanced=c(TRUE),
#   nuisance=c(0,4,16),
#   rep=1)
# 
# length(settings)
# 
# usecores <- 1
# options(cores=usecores)
# M50N1G100 <- try(doSim(settings=settings, cores=usecores))
# 
# save(M50N1G100, file=paste(pathResults, "NEWplotModelsM50N1G100.Rdata", sep=""))
# 
# 
# 
# # ###################################### M=100
# 
# set.seed(18102012)
# 
# settings <- makeSettings(
#   M=c(100),
#   ni=c(1),
#   Gy=c(30),
#   Gx=c(30),
#   snrEps=c(1,5),
#   snrE=c(0),
#   snrB=c(2),
#   scenario=c(1:4),
#   balanced=c(TRUE),
#   nuisance=c(0,4,16),
#   rep=1)
# 
# length(settings)
# 
# usecores <- 1
# options(cores=usecores)
# M100N1G30 <- try(doSim(settings=settings, cores=usecores))
# 
# save(M100N1G30, file=paste(pathResults, "plotModelsM100N1G30.Rdata", sep=""))
# 
# 
# 
# set.seed(18102012)
# 
# settings <- makeSettings(
#   M=c(100),
#   ni=c(1),
#   Gy=c(100),
#   Gx=c(100),
#   snrEps=c(1,5),
#   snrE=c(0),
#   snrB=c(2),
#   scenario=c(1:4),
#   balanced=c(TRUE),
#   nuisance=c(0,4,16),
#   rep=1)
# 
# length(settings)
# 
# usecores <- 1
# options(cores=usecores)
# M100N1G100 <- try(doSim(settings=settings, cores=usecores))
# 
# save(M100N1G100, file=paste(pathResults, "plotModelsM100N1G100.Rdata", sep=""))
# 
# 


# # ###################################### M=200
# 
# set.seed(18102012)
# 
# settings <- makeSettings(
#   M=c(200),
#   ni=c(1),
#   Gy=c(30),
#   Gx=c(30),
#   snrEps=c(1,5),
#   snrE=c(0),
#   snrB=c(2),
#   scenario=c(1:4),
#   balanced=c(TRUE),
#   nuisance=c(0,4,16),
#   rep=1)
# 
# length(settings)
# 
# usecores <- 1
# options(cores=usecores)
# M200N1G30 <- try(doSim(settings=settings, cores=usecores))
# 
# save(M200N1G30, file=paste(pathResults, "plotModelsM200N1G30.Rdata", sep=""))
# 
# 
# 
# set.seed(18102012)
# 
# settings <- makeSettings(
#   M=c(200),
#   ni=c(1),
#   Gy=c(100),
#   Gx=c(100),
#   snrEps=c(1,5),
#   snrE=c(0),
#   snrB=c(2),
#   scenario=c(1:4),
#   balanced=c(TRUE),
#   nuisance=c(0,4,16),
#   rep=1)
# 
# length(settings)
# 
# usecores <- 1
# options(cores=usecores)
# M200N1G100 <- try(doSim(settings=settings, cores=usecores))
# 
# save(M200N1G100, file=paste(pathResults, "plotModelsM200N1G100.Rdata", sep=""))
