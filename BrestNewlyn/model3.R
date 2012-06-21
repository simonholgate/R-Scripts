#######################
## Model 3 1953-2008 ##
#######################

model <- 3

## Use robust methods
library(robust)
lmRobControl <- lmRob.control(mxr=100,mxf=100,trace=F)

## Set start and end years of full model
start <- 1953
end <- 2008

## Set start and end years of partial model
## Roughly 70-30 split

p1b.start <- 1953 # Brest training data start
p1b.end <- 1989 # Brest training data end
p2b.start <- 1990 # Brest prediction start
p2b.end <- 2008 # Brest prediction end

p1n.start <- p1b.start # Newlyn training data start
p1n.end <- p1b.end # Newlyn training data end
p2n.start <- p2b.start # Newlyn prediction start
p2n.end <- p2b.end # Newlyn prediction end

source("~/Dropbox/BrestNewlyn/brestNewlynFunctions.R")

## Calculate correlation array. What is the relationship between each of the pressure series and the sea level time-series?
## Columns are Brest, Newlyn then Pa
corr.array <- array(NA,dim=c(met$hadnstns,met$hadnstns))
corr.data <- cbind(brestTotalP,newlynTotalP,Pa)
for(i in 1:met$hadnstns){
  for(j in 1:met$hadnstns){
    corr.array[i,j] <- cor.test(corr.data[,i], corr.data[,j], method="p", alternative="t", na.action=na.fail)$estimate
  }
}

### How does a model based on the most correlated 9 components, compare?
#tg.cor.sorted.brest <- sort(abs(corr.array[1,3:37]), index=T, decreasing=T)$ix[1:9]
#data.brest.partial <- data.frame(msl=brestTotalP[1:37], t=seq(from=start,to=1989), Pa=Pa[1:37,tg.cor.sorted.brest])
#data.brest.new <- data.frame(t=seq(from=1990,to=end), Pa=Pa[38:56,3:37])
#tg.cor.lmRob.brest.partial <- lmRob(msl ~ ., data=data.brest.partial)
#tg.cor.lmRob.brest.pred <- predict.lmRob(tg.cor.lmRob.brest.partial, se.fit=T, newdata=data.brest.new, interval='confidence')

## The predictions from the most correlated pressures are worse than from those with the R^2 technique


# Variance reduction
#> (tg.lmRob.newlyn$r.sq-tg.lmRob.newlynTotalP$r.sq)
#[1] 0.273999
#> tg.lmRob.newlyn$r.sq
#[1] 0.7001414
#> c(coef(tg.lmRob.brest)[2]/10.045, coef(tg.lmRob.newlyn)[2]/10.045)
#       t        t 
#1.822623 1.798708 
#> c(sqrt(diag(tg.lmRob.brest$cov))[2]/10.045,  sqrt(diag(tg.lmRob.newlyn$cov))[2]/10.045)
#        t         t 
#0.2615218 0.2927652

x11()
plot(time,newlynTotalP,col='blue',type='l', ylim=c(0,3000))
lines(time[1:55],tg.lmRob.newlyn$fitted,col='red')
lines(time[p1p.start:p1p.end],tg.lmRob.newlyn.partial$fit,col='cyan', lwd=2)
lines(time[p2p.start:p2p.end],tg.lmRob.newlyn.pred$fit,col='orange', lwd=2)
grid(lwd=2)
x11()
plot(time,brestTotalP,col='blue',type='l', ylim=c(0,3000))
lines(time,tg.lmRob.brest$fitted,col='red')
lines(time[p1p.start:p1p.end],tg.lmRob.brest.partial$fit,col='cyan', lwd=2)
lines(time[p2p.start:p2p.end],tg.lmRob.brest.pred$fit,col='orange', lwd=2)
grid(lwd=2)
