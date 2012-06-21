#######################
## Model 5 1916-2008 ##
#######################

model <- 5

## Use robust methods
library(robust)
lmRobControl <- lmRob.control(mxr=100,mxf=100,trace=F)

## Set start and end years of full model
start <- 1916
end <- 2008

## Set start and end years of partial model
## Roughly 70-30 split

p1b.start <- 1953 # Brest training data start
p1b.end <- 2008 # Brest training data end
p2b.start <- 1916 # Brest prediction start
p2b.end <- 1943 # Brest prediction end

p1n.start <- p1b.start # Newlyn training data start
p1n.end <- p1b.end # Newlyn training data end
p2n.start <- p2b.start # Newlyn prediction start
p2n.end <- p2b.end # Newlyn prediction end

source("~/Dropbox/BrestNewlyn/brestNewlynFunctionsACRE.R")


# Variance reduction
#> (tg.lmRob.newlyn$r.sq-tg.lmRob.newlynTotalP$r.sq)
#[1] 0.07981396 cf 0.1263007
#> tg.lmRob.newlyn$r.sq
#[1] 0.7305994 cf 0.7766394
#> (tg.lmRob.brest$r.sq-tg.lmRob.brestTotalP$r.sq)
#[1] 0.1131965 cf 0.1295884
#> tg.lmRob.brest$r.sq
#[1] 0.6617371 cf 0.6778706
#> (tg.lmRob.brest_22$r.sq-tg.lmRob.brestTotalP_22$r.sq)
#[1] 0.06603757 cf 0.09613512
#> tg.lmRob.brest_22$r.sq
#[1] 0.692955 cf 0.7228666

#> c(coef(tg.lmRob.brest)[2]/10.045, coef(tg.lmRob.newlyn)[2]/10.045)
#       t        t 
#1.634816 1.721597 cf 1.462220 1.828716
#> c(sqrt(diag(tg.lmRob.brest$cov))[2]/10.045,  sqrt(diag(tg.lmRob.newlyn$cov))[2]/10.045)
#        t         t 
#0.1830476 0.1311725 cf 0.08604919 0.09671664

#> c(coef(tg.lmRob.brest_22)[2]/10.045, coef(tg.lmRob.newlyn)[2]/10.045)
#       t        t 
#1.791970 1.721597 cf 1.794447 1.828716
#> c(sqrt(diag(tg.lmRob.brest_22$cov))[2]/10.045,  sqrt(diag(tg.lmRob.newlyn$cov))[2]/10.045)
#        t         t 
#0.1696494 0.1311725 cf 0.10862234 0.09671664


x11()
#png(file="model5Newlyn.png")
mNewlynTotalP <- mean(newlynTotalP, na.rm=T)
plot(time,newlynTotalP-mNewlynTotalP,col='blue',type='l', ylim=c(-1000,1000), ann=F, lwd=2)
lines(tg.lmRob.newlyn$x[,2],tg.lmRob.newlyn$fitted-mNewlynTotalP,col='red')
lines(tg.lmRob.newlyn.partial$x[,2],tg.lmRob.newlyn.partial$fit-mNewlynTotalP,col='cyan', lwd=2)
lines(time[p2p.start:p2p.end],tg.lmRob.newlyn.pred$fit-mNewlynTotalP,col='orange', lwd=2)
grid(lwd=2)
title(main="Newlyn: projection to 1916", ylab="Pressure Anomaly [Pa]", xlab="Year")
legend(1916,3000, c("Total pressure", "Regressed pressure", "Training set", "Prediction"),
       fill=NA, col=c("blue", "red","cyan","orange"), lty=1, lwd=2, bty='n', border=NA)
#dev.off()


##*******************************************************************************************************
##
## Multi-plot of Brest
##
x11()
#png(file="model5Brest.png")
## first image
split.screen(c(2,1))
screen(1)
mBrestTotalP <- mean(brestTotalP, na.rm=T)
plot(time,brestTotalP-mBrestTotalP,col='blue',type='l', ylim=c(-1700,1700), ann=F, lwd=2)
lines(tg.lmRob.brest$x[,2],tg.lmRob.brest$fitted-mBrestTotalP,col='red')
lines(time[p1p.start:p1p.end],tg.lmRob.brest.partial$fit-mBrestTotalP,col='cyan', lwd=2)
lines(time[p2p.start:p2p.end],tg.lmRob.brest.pred$fit-mBrestTotalP,col='orange', lwd=2)
grid(lwd=2)
title(main="Brest: projection to 1916", ylab="Pressure Anomaly [Pa]")
legend(1910,1700, c("Total pressure", "Regressed pressure", "Training set", "Prediction"),
       fill=NA, col=c("blue", "red","cyan","orange"), lty=1, lwd=2, bty='n', border=NA)
## second image
screen(2)
mBrestTotalP_22<- mean(brestTotalP_22, na.rm=T)
plot(time,brestTotalP_22-mBrestTotalP_22,col='blue',type='l', ylim=c(-1700,1700), ann=F, lwd=2)
lines(tg.lmRob.brest_22$x[,2],tg.lmRob.brest_22$fitted-mBrestTotalP_22,col='red')
lines(time[p1p.start:p1p.end],tg.lmRob.brest.partial_22$fit-mBrestTotalP_22,col='cyan', lwd=2)
lines(time[p2p.start:p2p.end],tg.lmRob.brest.pred_22$fit-mBrestTotalP_22,col='orange', lwd=2)
grid(lwd=2)
title(main="Brest22: projection to 1916 with 22mm 'correction' pre-1943.", ylab="Pressure Anomaly [Pa]", xlab="Year")
legend(1910,1700, c("Total pressure", "Regressed pressure", "Training set", "Prediction"),
       fill=NA, col=c("blue", "red","cyan","orange"), lty=1, lwd=2, bty='n', border=NA)
close.screen( all=TRUE)
#dev.off()
