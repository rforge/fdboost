###############################################################################
# changed code of Fabian Scheipl lfpr3_analysis.R
# Author: Sarah Brockhaus
###############################################################################

# M=50
load("M50N1G30pffr.Rdata")
M50N1G30 <- M50N1G30[,-1]

load("M50N1G100pffr.Rdata") 
M50N1G100 <- M50N1G100[,-1]

# M=100
load("M100N1G30pffr.Rdata")
M100N1G30 <- M100N1G30[,-1]

load("M100N1G100pffr.Rdata")
M100N1G100 <- M100N1G100[,-1]

# M=200
load("M200N1G30pffr.Rdata")
M200N1G30 <- M200N1G30[,-1]

load("M200N1G100pffr.Rdata")
M200N1G100 <- M200N1G100[,-1]

# M=50
load("M50N1G30FDboost.Rdata")
M50N1G30FDboost <- M50N1G30FDboost[,1:ncol(M50N1G30)]

load("M50N1G100FDboost.Rdata")
M50N1G100FDboost <- M50N1G100FDboost[,1:ncol(M50N1G30)]

# M=100
load("M100N1G30FDboost.Rdata")
M100N1G30FDboost <- M100N1G30FDboost[,1:ncol(M50N1G30)]

load("M100N1G100FDboost.Rdata")
M100N1G100FDboost <- M100N1G100FDboost[,1:ncol(M50N1G30)]

# M=200
load("M200N1G30FDboost.Rdata")
M200N1G30FDboost <- M200N1G30FDboost[,1:ncol(M50N1G30)]

load("M200N1G100FDboost.Rdata")
M200N1G100FDboost <- M200N1G100FDboost[,1:ncol(M50N1G30)]

    
  resNames <- ls()[grep("^M", ls())]
#  resNames <- ls()[grep("^S", ls())]
      
  res <- get(resNames[1])
  rStart <- 2
  if(any(class(res)=="try-error")){
    res <- get(resNames[2])
    rStart <- 3
  }
    
  for(r in rStart:length(resNames)){    
    if(any(class(get(resNames[r]))=="try-error")){
        get(resNames[r]) <- NULL
        print(paste(resNames[r], "has class try-error."))
    }else{
        res <- rbind(res, get(resNames[r])) 
      }     
  }
  
  #summary(res)
  dim(res)
  
  res$Mf <- factor(res$M)
  res$nif <- factor(res$ni)
  res$Gyf <- factor(res$Gy)
  res$snrEpsf <- factor(res$snrEps)
  res$snrEf <- factor(res$snrE)
  res$snrBf <- factor(res$snrB)
  res$scenario <- factor(res$scenario)
  res$balancedf <- factor(res$balanced)
  
  res$modelf <- factor(res$model)
  res$model0 <- as.numeric(res$model)-1   
  
  res$Mnif <- res$Mf:res$nif

  # total number of trajectories
  res$N <- res$M*res$ni
  res$Nf <- factor(res$N)

res$mod <- factor(res$model, labels = c("PFFR", "FDboost") )
res$nuisanceV <- factor(res$nuisance)
  
  save(res, file="boosting.Rdata")
#}

#load("boosting.Rdata")

# library(plyr)
# maxErrors <- ddply(res, ~ M + ni + scenario + snrEps, function(d) {
#             return(d[which.max(d$imsey),])
#         })
# minErrors <- ddply(res, ~ M + ni + scenario + snrEps, function(d) {
#             return(d[which.min(d$imsey),])
#         })
# medianErrors <- ddply(res, ~ M + ni + scenario + snrEps, function(d) {
#             return(d[which(rank(d$imsey)==floor(nrow(d)/2)),])
#         })

# for each scenario: errors on M/n vs snrE/snrB and snrE

library(mgcv)
library(ggplot2)


res$"M" <- factor(res$M) # , labels = paste("M:", names(table(res$M)))
#res$"n_i" <- factor(res$ni)
res$"G" <- factor(res$Gy) #, labels = paste("G:", names(table(res$Gy)))
#res$`$\\operatorname{SNR}_\\varepsilon$` <- factor(res$snrEps)
#res$`$\\operatorname{SNR}_B$` <- factor(res$snrB)
res$"snrEps" <- factor(res$snrEps)
res$"snrB" <- factor(res$snrB)  
res$"NG" <- factor(res$Nf:res$G) #, 

################## plot results


pdf("ComputationTime.pdf", width=9, height=7)
ggplot(subset(res, !is.na(time.elapsed)), aes(y=time.elapsed, x=snrEps, fill=mod)) + 
  geom_boxplot(outlier.size=.6) + 
  facet_grid(N~G , labeller="label_both") +
  #theme_clear(base_size=14)   + 
  theme(legend.position="top", legend.direction="horizontal") +
  scale_fill_manual(name = "", values=c("white", "grey80")) +                    
  scale_colour_manual(name = "", values=c("grey40", "black")) +
  labs(x="snrEps", y="time") +
  scale_y_continuous(breaks=c(1, 2, 5, 10, 20, 60, 120, 300, 600, 1200, 2700, 
                              5400, 10800, 21600, 43200), 
                     trans="log10",
                     labels=c("1s", "2s", "5s", "10s", "20s", "1 min", "2 min", "5 min", 
                              "10 min", "20 min", "45 min", "90 min", "3h", "6h", "12h")) +
  #labs(title="Computation time") + 
  xlab(bquote(SNR[epsilon]))  +
  theme(text=element_text(size = 30))
dev.off()



pdf("reliMSEy.pdf", width=9, height=7)
ggplot(subset(res, !is.na(relmsey)), aes(y=relmsey, fill=mod, colour=mod, x=snrEps)) +
  geom_boxplot(outlier.size=.6) +
  facet_grid( N  ~  G, labeller="label_both") + #G
  scale_y_log10() + 
  theme(legend.position="top", legend.direction="horizontal") +
  scale_fill_manual(name = "", values=c("white", "grey80")) +              
  scale_colour_manual(name = "", values=c("grey40", "black")) + 
  #labs(title="riMSEy")
  #labs(title="reliMSE(Y(t))") + 
  ylab("reliMSE(Y(t))") + xlab(bquote(SNR[epsilon])) +
  theme(text=element_text(size = 30))
dev.off()


pdf("reliMSEg0.pdf", width=9, height=7)
ggplot(subset(res, !is.na(relmseg0)), aes(y=relmseg0, fill=mod, colour=mod, x=snrEps)) +
  geom_boxplot(outlier.size=.6) +
  facet_grid( N  ~  G, labeller="label_both") + #G
  scale_y_log10() + 
  theme(legend.position="top", legend.direction="horizontal") +
  scale_fill_manual(name = "", values=c("white", "grey80")) +              
  scale_colour_manual(name = "", values=c("grey40", "black")) + 
  ylab(bquote(reliMSE(beta[0](t)))) + xlab(bquote(SNR[epsilon])) +
  theme(text=element_text(size = 30))
dev.off()

pdf("reliMSEfx1.pdf", width=9, height=7)
ggplot(subset(res, !is.na(relmsefx1)), aes(y=relmsefx1, fill=mod, colour=mod, x=snrEps)) +
  geom_boxplot(outlier.size=.6) +
  facet_grid( N  ~  G, labeller="label_both") + #G
  scale_y_log10() + 
  theme(legend.position="top", legend.direction="horizontal") +
  scale_fill_manual(name = "", values=c("white", "grey80")) +              
  scale_colour_manual(name = "", values=c("grey40", "black")) + 
  #labs(title=bquote(reliMSE(beta[1](s,t)))) +
  ylab(bquote(reliMSE(beta[1](s,t)))) + xlab(bquote(SNR[epsilon])) +
  theme(text=element_text(size = 30))
dev.off()

pdf("reliMSEfx2.pdf", width=9, height=7)
ggplot(subset(res, !is.na(relmsefx2)), aes(y=relmsefx2, fill=mod, colour=mod, x=snrEps)) +
  geom_boxplot(outlier.size=.6) +
  facet_grid( N  ~  G, labeller="label_both") + #G
  scale_y_log10() + 
  theme(legend.position="top", legend.direction="horizontal") +
  scale_fill_manual(name = "", values=c("white", "grey80")) +              
  scale_colour_manual(name = "", values=c("grey40", "black")) + 
  #labs(title=bquote(reliMSE(beta[2](s,t)))) +
  ylab(bquote(reliMSE(beta[2](s,t)))) + xlab(bquote(SNR[epsilon])) +
  theme(text=element_text(size = 30))
dev.off()


