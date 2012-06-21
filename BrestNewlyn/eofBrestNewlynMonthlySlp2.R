# Method based on Bjornsson and Venegas

# Monthly data
eof <- new.env()
load("/home/simonh/Dropbox/brestNewlynData/pressure/brestNewlynHadSLP2Paper2.RData",envir=eof)

# slpNewlynStns is a matrix of 1920 rows and 37 columns, with each column a time-series for a given station
# Remove mean from each column
eof$dmSlpNewlynStns <- eof$slpNewlynStns[,3:37]
eof$cmSlpNewlynStns <- colMeans(eof$slpNewlynStns)
eof$numStns <- length(eof$cmSlpNewlynStns)
for(i in 1:eof$numStns){
eof$dmSlpNewlynStns[,i] <- eof$dmSlpNewlynStns[,i] - eof$cmSlpNewlynStns[i]
}

# Form the co-variance matrix
eof$covDmSlpNewlynStns <- t(eof$dmSlpNewlynStns) %*% eof$dmSlpNewlynStns

# Calculate eigenvalues and eigenvectors
eof$eigCovDmSlpNewlynStns <- eigen(eof$covDmSlpNewlynStns)

# Normalise the EOFs so that they are dimensionless and the maximum value is 100
for (i in 1:eof$numStns){
  eof$eigCovDmSlpNewlynStns$vectors[,i] <- 100*eof$eigCovDmSlpNewlynStns$vectors[,i]/max(eof$eigCovDmSlpNewlynStns$vectors[,i])
}

# EOF1
eof$EOF1 <- eof$eigCovDmSlpNewlynStns$vectors[,1]
dim(eof$EOF1) <- c(7,5)

# EOF2
eof$EOF2 <- eof$eigCovDmSlpNewlynStns$vectors[,2]
dim(eof$EOF2) <- c(7,5)

# Find expansion coefficients
#eof$C <- diag(x=eof$eigCovDmSlpNewlynStns$values,nrow=eof$numStns, ncol=eof$numStns)
eof$C <- diag(x=eof$eigCovDmSlpNewlynStns$values)
#for (i in 1:numStns){
eof$PC <- eof$dmSlpNewlynStns %*% eof$eigCovDmSlpNewlynStns$vectors[,1]

#Pattern of pressures for a single month
eof$slpNewlynStnsArray <- eof$slpNewlynStns[1,3:37]
dim(eof$slpNewlynStnsArray) <- c(7,5)

#quartz()
x11()
contour(seq(from=-20,to=10,by=5),seq(from=40,to=60,by=5), eof$slpNewlynStnsArray)
# Map
library(maps)
library(mapdata)

#quartz()
x11()
map("worldHires", xlim=c(-25,15), ylim=c(35,65), interior=F, fill=F, col="grey50", resolution=0)
map.axes()
#EOF1
contour(seq(from=-20,to=10,by=5),seq(from=40,to=60,by=5), eof$EOF1, lwd=2, add=T)
#EOF2
contour(seq(from=-20,to=10,by=5),seq(from=40,to=60,by=5), eof$EOF2, lty='dotdash', add=T)

# Time series pattern
eof$tpcaSlpNewlynStns <- prcomp(eof$slpNewlynStns, scale=TRUE)
summary(eof$tpcaSlpNewlynStns)
eof$pTpcaSlpNewlynStns1 <- predict(eof$tpcaSlpNewlynStns)[,1]
eof$pTpcaSlpNewlynStns2 <- predict(eof$tpcaSlpNewlynStns)[,2]

#quartz()
x11()
plot(1:1920,eof$pTpcaSlpNewlynStns1, type='l', col='red')
lines(1:1920,eof$pTpcaSlpNewlynStns2, col='blue')

# Annual data
ann <- new.env()
tmp <- new.env()
# HadSLP2r data is 1850-2009 = 160 years = 1920 months. 37 cols - Newlyn, Brest and 35 grid points
load("~/Dropbox/brestNewlynData/pressure/brestNewlynHadSLP2Paper2.RData", envir=tmp)
tmp$hadnstns <- 37
ann$hadSLP2rAnnualPressureArray <- array(NA, dim=c(160,tmp$hadnstns))
for(i in 1:tmp$hadnstns){
  tmp$hadSLP2rMonthlyPressureArray <- tmp$slpNewlynStns[,i]
  dim(tmp$hadSLP2rMonthlyPressureArray) <- c(12,160)
  ann$hadSLP2rAnnualPressureArray[,i] <- colMeans(tmp$hadSLP2rMonthlyPressureArray, na.rm=T)
}
ann$hadSLP2rAnnualTime <- c(1850:2009)
rm(tmp)

# Spatial loading pattern
ann$spcaSlpNewlynStns <- prcomp(t(ann$hadSLP2rAnnualPressureArray), scale=TRUE)
summary(ann$spcaSlpNewlynStns)
#Pattern of pressures for a single year
ann$slpNewlynStnsArray <- ann$hadSLP2rAnnualPressureArray[1,3:37]
dim(ann$slpNewlynStnsArray) <- c(7,5)
#quartz()
x11()
contour(seq(from=-20,to=10,by=5),seq(from=40,to=60,by=5), ann$slpNewlynStnsArray)
# Map
library(maps)
library(mapdata)
#quartz()
x11()
map("worldHires", xlim=c(-25,15), ylim=c(35,65), interior=F, fill=F, col="grey50", resolution=0)
map.axes()
#EOF1
ann$pSpcaSlpNewlynStns1 <- predict(ann$spcaSlpNewlynStns)[,1]
ann$pSpcaSlpNewlynStns1Array <- ann$pSpcaSlpNewlynStns1[3:37]
dim(ann$pSpcaSlpNewlynStns1Array) <- c(7,5)
contour(seq(from=-20,to=10,by=5),seq(from=40,to=60,by=5), ann$pSpcaSlpNewlynStns1Array, lwd=2, add=T)
#EOF2
ann$pSpcaSlpNewlynStns2 <- predict(ann$spcaSlpNewlynStns)[,2]
ann$pSpcaSlpNewlynStns2Array <- ann$pSpcaSlpNewlynStns2[3:37]
dim(ann$pSpcaSlpNewlynStns2Array) <- c(7,5)
contour(seq(from=-20,to=10,by=5),seq(from=40,to=60,by=5), ann$pSpcaSlpNewlynStns2Array, lty='dotdash', add=T)
#EOF3
ann$pSpcaSlpNewlynStns3 <- predict(ann$spcaSlpNewlynStns)[,3]
ann$pSpcaSlpNewlynStns3Array <- ann$pSpcaSlpNewlynStns3[3:37]
dim(ann$pSpcaSlpNewlynStns3Array) <- c(7,5)
contour(seq(from=-20,to=10,by=5),seq(from=40,to=60,by=5), ann$pSpcaSlpNewlynStns3Array, lty='dotted', add=T)

# Time series pattern
ann$tpcaSlpNewlynStns <- prcomp(ann$hadSLP2rAnnualPressureArray, scale=TRUE)
summary(ann$tpcaSlpNewlynStns)
ann$pTpcaSlpNewlynStns1 <- predict(ann$tpcaSlpNewlynStns)[,1]
ann$pTpcaSlpNewlynStns2 <- predict(ann$tpcaSlpNewlynStns)[,2]
ann$pTpcaSlpNewlynStns3 <- predict(ann$tpcaSlpNewlynStns)[,3]

#quartz()
x11()
plot(ann$hadSLP2rAnnualTime,ann$pTpcaSlpNewlynStns1, type='l', col='red')
lines(ann$hadSLP2rAnnualTime,ann$pTpcaSlpNewlynStns2, col='blue')
lines(ann$hadSLP2rAnnualTime,ann$pTpcaSlpNewlynStns3, col='magenta')
