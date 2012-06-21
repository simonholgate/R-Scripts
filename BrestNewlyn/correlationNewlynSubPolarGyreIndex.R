library(fields)
library(robust)
source("~/Dropbox/BrestNewlyn/matrixMethods.R")
load("~/Dropbox/brestNewlynData/analysis/paper/correlationACRE/brestNewlyn.tot.ps.RData")

newlyn.full.p.start <- 1916
newlyn.full.p.end <- 2008
newlyn.full.p.yrs <- c(newlyn.full.p.start:newlyn.full.p.end)

hatun.period <- c(which(newlyn.full.p.yrs>=1960):which(newlyn.full.p.yrs==2004))
hatun.yrs <- newlyn.full.p.yrs[hatun.period]
hatun.tot.p <- newlyn.tot.p.full[hatun.period]

filt.hatun.tot.p <- filter(hatun.tot.p, c(0.25, 0.5, 0.25), method="conv", sides=2)

plot(hatun.yrs, -hatun.tot.p, type='l', col='red', lwd=2)
lines(hatun.yrs, -filt.hatun.tot.p, col='blue', lwd=2, lty=2)
grid()

hakkinen.period <- c(which(newlyn.full.p.yrs==1993):which(newlyn.full.p.yrs==2001))
hakkinen.tot.p <- newlyn.tot.p.full[hakkinen.period]

hakkinen <- read.table(file="hakkinenData.txt", sep=",", colClasses=c("numeric", "numeric"),
                       col.names=c("year", "pc1"))
##x11()
##plot(hakkinen$year, hakkinen$pc1, type='l', col='blue')

intsct <- intersect(which(hakkinen$year>=1993), which(hakkinen$year<=2002))
hakkinen.yrs <- 1993:2001
hakkinen.annual <- hakkinen$pc1[intsct]
dim(hakkinen.annual) <- c(12,length(hakkinen.yrs))
hakkinen.annual <- colMeans(hakkinen.annual)
## Re-scaling Hakkinen data to fit Newlyn
##range(hakkinen.annual) -> -0.1181667  0.1391667
##which(hatun.yrs==1993) -> 34
##which(hatun.yrs==2002) -> 43
##diff(range(-hatun.tot.p[34:43])) -> 1164.665
##mean(-hatun.tot.p[34:43]) -> -102490.1
lines(hakkinen.yrs, hakkinen.annual/0.139/2*1164-102490, col='magenta', lwd=2)

hatun <- read.table(file="hatunData.txt", sep=",", colClasses=c("numeric", "numeric"),
                       col.names=c("year", "pc1"))

##plot(hatun$year, hatun$pc1, type='l', col='blue', lty=3)
## Re-scaling Hatun data to fit Newlyn
##range(hatun$pc1) -> -8.1 10.3
##diff(range(-hatun.tot.p)) -> 1660.483
##mean(-hatun.tot.p[34:43]) -> -102232.8
lines(hatun$year, hatun$pc1/10.3/2*1660-102232, col='cyan', lwd=2, lty=3)

## Normalised plot
plot(hatun.yrs, (-hatun.tot.p - mean(-hatun.tot.p))/1076.7611, type='l', col='red', lwd=2, ylim=c(-1,1), ann=F)
##lines(hatun.yrs, (-filt.hatun.tot.p - mean(-filt.hatun.tot.p, na.rm=T))/764.2831, col='blue', lwd=2, lty=2)
grid(lwd=2)
title(main="Newlyn sea level and gyre index", ylab="Normalised sea level/gyre index", xlab="Years")
lines(hakkinen.yrs, (hakkinen.annual - mean(hakkinen.annual))/0.1314722, col='magenta', lwd=2)
lines(hatun$year, (hatun$pc1 - mean(hatun$pc1))/9.218182, col='cyan', lwd=2)

intsct <- intersect(which(hatun$year>=1993), which(hatun$year<2002))
hakkinen.hatun <- hatun$pc1[intsct]
  
## Normalised plot over Hakkinen yrs
plot(hakkinen.yrs, (-hakkinen.tot.p - mean(-hakkinen.tot.p))/254.2495, type='l', col='red', lwd=2, ylim=c(-1,1), ann=F)
##lines(hatun.yrs, (-filt.hatun.tot.p - mean(-filt.hatun.tot.p, na.rm=T))/764.2831, col='blue', lwd=2, lty=2)
grid(lwd=2)
title(main="Newlyn sea level and gyre index", ylab="Normalised sea level/gyre index", xlab="Years")
lines(hakkinen.yrs, (hakkinen.annual - mean(hakkinen.annual))/0.1314722, col='magenta', lwd=2)
lines(hakkinen.yrs, (hakkinen.hatun - mean(hakkinen.hatun))/9.666667, col='cyan', lwd=2)

## Normalised plot over Hakkinen yrs shown over Hatun yrs
plot(hatun.yrs, (-hatun.tot.p - mean(-hakkinen.tot.p))/254.2495, type='l', col='red', lwd=2, ylim=c(-1,1), ann=F)
##lines(hatun.yrs, (-filt.hatun.tot.p - mean(-filt.hatun.tot.p, na.rm=T))/764.2831, col='blue', lwd=2, lty=2)
grid(lwd=2)
title(main="Newlyn sea level and gyre index", ylab="Normalised sea level/gyre index", xlab="Years")
lines(hakkinen.yrs, (hakkinen.annual - mean(hakkinen.annual))/0.1314722, col='magenta', lwd=2)
lines(hatun$year, (hatun$pc1 - mean(hatun$pc1))/9.666667, col='cyan', lwd=2)
