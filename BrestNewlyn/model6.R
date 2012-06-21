#######################
## Model 6 1850-2008 ##
#######################

model <- 6

## No Newlyn......

## Use robust methods
library(robust)
lmRobControl <- lmRob.control(mxr=100,mxf=100,trace=F)

## Set start and end years of full model
start <- 1850
end <- 2008

## Set start and end years of partial model
## Roughly 70-30 split

p1b.start <- 1953 # Brest training data start
p1b.end <- 2008 # Brest training data end
p2b.start <- 1850 # Brest prediction start
p2b.end <- 1943 # Brest prediction end

source("~/Dropbox/BrestNewlyn/brestNewlynFunctions.R")


# Variance reduction
#> (tg.lmRob.brest$r.sq-tg.lmRob.brestTotalP$r.sq)
#[1] 0.1295884
#> tg.lmRob.brest$r.sq
#[1] 0.6778706
#> (tg.lmRob.brest_22$r.sq-tg.lmRob.brestTotalP_22$r.sq)
#[1] 0.09613512
#> tg.lmRob.brest_22$r.sq
#[1] 0.7228666

#> c(coef(tg.lmRob.brest)[2]/10.045, coef(tg.lmRob.newlyn)[2]/10.045)
#       t        t 
#1.462220 1.828716
#> c(sqrt(diag(tg.lmRob.brest$cov))[2]/10.045,  sqrt(diag(tg.lmRob.newlyn$cov))[2]/10.045)
#        t         t 
#0.08604919 0.09671664

#> c(coef(tg.lmRob.brest_22)[2]/10.045, coef(tg.lmRob.newlyn)[2]/10.045)
#       t        t 
#1.794447 1.828716
#> c(sqrt(diag(tg.lmRob.brest_22$cov))[2]/10.045,  sqrt(diag(tg.lmRob.newlyn$cov))[2]/10.045)
#        t         t 
#0.10862234 0.09671664

##*******************************************************************************************************
##
## Multi-plot of Brest
##
split.screen(c(2,1))
## first image
screen(1)
plot(time,brestTotalP,col='blue',type='l', ylim=c(-1000,3000), ann=F)
lines(tg.lmRob.brest$x[,2],tg.lmRob.brest$fitted,col='red')
lines(time[p1p.start:p1p.end],tg.lmRob.brest.partial$fit,col='cyan', lwd=2)
lines(time[p2p.start:p2p.end],tg.lmRob.brest.pred$fit,col='orange', lwd=2)
grid(lwd=2)
title(main="Brest: projection to 1850", ylab="Pressure [Pa]", xlab="Year")
legend(1850,3000, c("Total pressure", "Regressed pressure", "Training set", "Prediction"),
       fill=NA, col=c("blue", "red","cyan","orange"), lty=1, lwd=2, bty='n', border=NA)
## second image
screen(2)
plot(time,brestTotalP_22,col='blue',type='l', ylim=c(-1000,3000), ann=F)
lines(tg.lmRob.brest_22$x[,2],tg.lmRob.brest_22$fitted,col='red')
lines(time[p1p.start:p1p.end],tg.lmRob.brest.partial_22$fit,col='cyan', lwd=2)
lines(time[p2p.start:p2p.end],tg.lmRob.brest.pred_22$fit,col='orange', lwd=2)
grid(lwd=2)
title(main="Brest22: projection to 1850 with 22mm 'correction' pre-1943", ylab="Pressure [Pa]", xlab="Year")
legend(1850,3000, c("Total pressure", "Regressed pressure", "Training set", "Prediction"),
       fill=NA, col=c("blue", "red","cyan","orange"), lty=1, lwd=2, bty='n', border=NA)

close.screen( all=TRUE)
