# Monthly data
pca <- new.env()
load("/home/simonh/Dropbox/brestNewlynData/pressure/brestNewlynHadSLP2Paper2.RData",envir=pca)

# Spatial loading pattern
pca$pcaSlpNewlynStns <- prcomp(t(pca$slpNewlynStns), scale=TRUE)
summary(pca$pcaSlpNewlynStns)
#Pattern of pressures for a single month
pca$slpNewlynStnsArray <- pca$slpNewlynStns[1,3:37]
dim(pca$slpNewlynStnsArray) <- c(7,5)

# Map
library(maps)
library(mapdata)

#quartz()
x11()
map("worldHires", xlim=c(-25,15), ylim=c(35,65), interior=F, fill=F, col="grey50", resolution=0)
map.axes()
contour(seq(from=-20,to=10,by=5),seq(from=40,to=60,by=5), pca$slpNewlynStnsArray, add=T)

#quartz()
x11()
map("worldHires", xlim=c(-25,15), ylim=c(35,65), interior=F, fill=F, col="grey50", resolution=0)
map.axes()
#EOF1
pca$pPcaSlpNewlynStns1 <- predict(pca$pcaSlpNewlynStns)[,1]
pca$pPcaSlpNewlynStns1Array <- pca$pPcaSlpNewlynStns1[3:37]
dim(pca$pPcaSlpNewlynStns1Array) <- c(7,5)
contour(seq(from=-20,to=10,by=5),seq(from=40,to=60,by=5), pca$pPcaSlpNewlynStns1Array, lwd=2, add=T)
#EOF2
pca$pPcaSlpNewlynStns2 <- predict(pca$pcaSlpNewlynStns)[,2]
pca$pPcaSlpNewlynStns2Array <- pca$pPcaSlpNewlynStns2[3:37]
dim(pca$pPcaSlpNewlynStns2Array) <- c(7,5)
contour(seq(from=-20,to=10,by=5),seq(from=40,to=60,by=5), pca$pPcaSlpNewlynStns2Array, lty='dotdash', add=T)

# Time series pattern
pca$tpcaSlpNewlynStns <- prcomp(pca$slpNewlynStns, scale=TRUE)
summary(pca$tpcaSlpNewlynStns)
pca$pTpcaSlpNewlynStns1 <- predict(pca$tpcaSlpNewlynStns)[,1]
pca$pTpcaSlpNewlynStns2 <- predict(pca$tpcaSlpNewlynStns)[,2]

#quartz()
x11()
plot(1:1920,pca$pTpcaSlpNewlynStns1, type='l', col='red')
lines(1:1920,pca$pTpcaSlpNewlynStns2, col='blue')

# Annual data
pcaann <- new.env()
tmp <- new.env()
# HadSLP2r data is 1850-2009 = 160 years = 1920 months. 37 cols - Newlyn, Brest and 35 grid points
load("~/Dropbox/brestNewlynData/pressure/brestNewlynHadSLP2Paper2.RData", envir=tmp)
tmp$hadnstns <- 37
pcaann$hadSLP2rAnnualPressureArray <- array(NA, dim=c(160,tmp$hadnstns))
for(i in 1:tmp$hadnstns){
  tmp$hadSLP2rMonthlyPressureArray <- tmp$slpNewlynStns[,i]
  dim(tmp$hadSLP2rMonthlyPressureArray) <- c(12,160)
  pcaann$hadSLP2rAnnualPressureArray[,i] <- colMeans(tmp$hadSLP2rMonthlyPressureArray, na.rm=T)
}
pcaann$hadSLP2rAnnualTime <- c(1850:2009)
rm(tmp)

# Spatial loading pattern
pcaann$spcaSlpNewlynStns <- prcomp(t(pcaann$hadSLP2rAnnualPressureArray), scale=TRUE)
summary(pcaann$spcaSlpNewlynStns)
#Pattern of pressures for a single year
pcaann$slpNewlynStnsArray <- pcaann$hadSLP2rAnnualPressureArray[1,3:37]
dim(pcaann$slpNewlynStnsArray) <- c(7,5)

#quartz()
x11()
map("worldHires", xlim=c(-25,15), ylim=c(35,65), interior=F, fill=F, col="grey50", resolution=0)
map.axes()
contour(seq(from=-20,to=10,by=5),seq(from=40,to=60,by=5), pcaann$slpNewlynStnsArray, add=T)

#quartz()
x11()
map("worldHires", xlim=c(-25,15), ylim=c(35,65), interior=F, fill=F, col="grey50", resolution=0)
map.axes()
#EOF1
pcaann$pSpcaSlpNewlynStns1 <- predict(pcaann$spcaSlpNewlynStns)[,1]
pcaann$pSpcaSlpNewlynStns1Array <- pcaann$pSpcaSlpNewlynStns1[3:37]
dim(pcaann$pSpcaSlpNewlynStns1Array) <- c(7,5)
contour(seq(from=-20,to=10,by=5),seq(from=40,to=60,by=5), pcaann$pSpcaSlpNewlynStns1Array, lwd=2, add=T)
#EOF2
pcaann$pSpcaSlpNewlynStns2 <- predict(pcaann$spcaSlpNewlynStns)[,2]
pcaann$pSpcaSlpNewlynStns2Array <- pcaann$pSpcaSlpNewlynStns2[3:37]
dim(pcaann$pSpcaSlpNewlynStns2Array) <- c(7,5)
contour(seq(from=-20,to=10,by=5),seq(from=40,to=60,by=5), pcaann$pSpcaSlpNewlynStns2Array, lty='dotdash', add=T)
#EOF3
pcaann$pSpcaSlpNewlynStns3 <- predict(pcaann$spcaSlpNewlynStns)[,3]
pcaann$pSpcaSlpNewlynStns3Array <- pcaann$pSpcaSlpNewlynStns3[3:37]
dim(pcaann$pSpcaSlpNewlynStns3Array) <- c(7,5)
contour(seq(from=-20,to=10,by=5),seq(from=40,to=60,by=5), pcaann$pSpcaSlpNewlynStns3Array, lty='dotted', add=T)

# Time series pattern
pcaann$tpcaSlpNewlynStns <- prcomp(pcaann$hadSLP2rAnnualPressureArray, scale=TRUE)
summary(pcaann$tpcaSlpNewlynStns)
pcaann$pTpcaSlpNewlynStns1 <- predict(pcaann$tpcaSlpNewlynStns)[,1]
pcaann$pTpcaSlpNewlynStns2 <- predict(pcaann$tpcaSlpNewlynStns)[,2]
pcaann$pTpcaSlpNewlynStns3 <- predict(pcaann$tpcaSlpNewlynStns)[,3]

#quartz()
x11()
plot(pcaann$hadSLP2rAnnualTime,pcaann$pTpcaSlpNewlynStns1, type='l', col='red')
lines(pcaann$hadSLP2rAnnualTime,pcaann$pTpcaSlpNewlynStns2, col='blue')
lines(pcaann$hadSLP2rAnnualTime,pcaann$pTpcaSlpNewlynStns3, col='magenta')
