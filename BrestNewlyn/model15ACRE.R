########################
## Model 15 1916-2008 ##
########################

model <- 15

## Uses N Pacific EOFs

## Use robust methods
library(robust)
lmRobControl <- lmRob.control(mxr=100,mxf=100,trace=F)

## Set start and end years of full model
start <- 1916
end <- 2008

## Set start and end years of partial model
## Roughly 70-30 split

p1b.start <- 1953 # San Francisco training data start
p1b.end <- 2008 # San Francisco training data end
p2b.start <- 1916 # San Francisco prediction start
p2b.end <- 1952 # San Francisco prediction end

p1n.start <- p1b.start # Key West training data start
p1n.end <- p1b.end # Key West training data end
p2n.start <- p2b.start # Key West prediction start
p2n.end <- p2b.end # Key West prediction end

nstns <- 2 # Only consider 1 EOF which contain 91.2% of the data

source("~/Dropbox/BrestNewlyn/keyWestFunctionsACRE.R")

# Variance reduction
#> (tg.lmRob.keywest$r.sq-tg.lmRob.keywestTotalP$r.sq)
#[1] 0.3700237 cf 0.273999
#> tg.lmRob.keywest$r.sq
#[1] 0.7254665 cf 0.7001414
#> c(coef(tg.lmRob.sanfran)[2]/10.045, coef(tg.lmRob.keywest)[2]/10.045)
#       t        t 
#1.718610 1.559207 cf 1.822623 1.798708 
#> c(sqrt(diag(tg.lmRob.sanfran$cov))[2]/10.045,  sqrt(diag(tg.lmRob.keywest$cov))[2]/10.045)
#        t         t 
#0.6081806 0.7192874 cf 0.2615218 0.2927652

x11()
#png(file="model15KeyWest.png")
mKeywestTotalP <- mean(keywestTotalP, na.rm=T)
plot(time,keywestTotalP-mKeywestTotalP,col='blue',type='l', ylim=c(-1700,1700), ann=F, lwd=2)
lines(tg.lmRob.keywest$x[,2],tg.lmRob.keywest$fitted-mKeywestTotalP,col='red')
lines(tg.lmRob.keywest.partial$x[,2],tg.lmRob.keywest.partial$fit-mKeywestTotalP,col='cyan', lwd=2)
lines(time[p2p.start:p2p.end],tg.lmRob.keywest.pred$fit-mKeywestTotalP,col='orange', lwd=2)
title(main="Model 15: Key West", xlab="Year", ylab="Pressure Anomaly [Pa]")
grid(lwd=2)
#dev.off()
x11()
#png(file="model15SanFran.png")
mSanfranTotalP <- mean(sanfranTotalP, na.rm=T)
plot(time,sanfranTotalP-mSanfranTotalP,col='blue',type='l', ylim=c(-1700,1700), ann=F, lwd=2)
lines(tg.lmRob.sanfran$x[,2],tg.lmRob.sanfran$fitted-mSanfranTotalP,col='red')
lines(time[p1p.start:p1p.end],tg.lmRob.sanfran.partial$fit-mSanfranTotalP,col='cyan', lwd=2)
lines(time[p2p.start:p2p.end],tg.lmRob.sanfran.pred$fit-mSanfranTotalP,col='orange', lwd=2)
title(main="Model 15: San Francisco", xlab="Year", ylab="Pressure Anomaly [Pa]")
grid(lwd=2)
#dev.off()
