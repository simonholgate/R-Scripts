# Using R method

# Monthly data
eof <- new.env()
load("/home/simonh/Dropbox/brestNewlynData/pressure/brestNewlynHadSLP2Paper2.RData",envir=eof)

# slpNewlynStns is a matrix of 1920 rows and 37 columns, with each column a time-series for a given station
# We only want columns 3:37 as 1 and 2 are Newlyn and Brest

# Remove mean from each column
eof$dmSlpNewlynStns <- eof$slpNewlynStns[,3:37]
eof$cmSlpNewlynStns <- colMeans(eof$dmSlpNewlynStns)
eof$numStns <- length(eof$cmSlpNewlynStns)
eof$numRows <- length(eof$slpNewlynStns[,1])

for(i in 1:eof$numStns){
eof$dmSlpNewlynStns[,i] <- eof$dmSlpNewlynStns[,i] - eof$cmSlpNewlynStns[i]
}

# Calculate EOFs
eof$pcaSlpMonthly <- princomp(eof$dmSlpNewlynStns, scale=TRUE, scores=TRUE)

# % variance explained
# Get the variance explaned from summary(eof$pcaSlpMonthly)

# EOF1
eof$PCA1 <- eof$pcaSlpMonthly$loadings[,1]
dim(eof$PCA1) <- c(7,5)

# EOF2
eof$PCA2 <- eof$pcaSlpMonthly$loadings[,2]
dim(eof$PCA2) <- c(7,5)

# EOF3
eof$PCA3 <- eof$pcaSlpMonthly$loadings[,3]
dim(eof$PCA3) <- c(7,5)

# Map
library(maps)
library(mapdata)

#quartz()
x11()
map("worldHires", xlim=c(-25,15), ylim=c(35,65), interior=F, fill=F, col="grey50", resolution=0)
map.axes()
#EOF1
contour(seq(from=-20,to=10,by=5),seq(from=40,to=60,by=5), eof$PCA1, lwd=2, add=T, col='red')
#EOF2
contour(seq(from=-20,to=10,by=5),seq(from=40,to=60,by=5), eof$PCA2, lty='dotdash', add=T, col='blue')
#EOF2
contour(seq(from=-20,to=10,by=5),seq(from=40,to=60,by=5), eof$PCA3, lty='dotted', add=T, col='green')

#quartz()
x11()
plot(1:eof$numRows,(eof$pcaSlpMonthly$scores[,1]), type='l', col='red')
lines(1:eof$numRows,(eof$pcaSlpMonthly$scores[,2]), col='blue')
#plot(1:eof$numRows,eof$PC[,1], type='l', col='blue')
#lines(1:eof$numRows,eof$PC[,2], col='green')

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
ann$pcaSlpNewlynStns <- princomp(ann$hadSLP2rAnnualPressureArray[,3:37], scale=TRUE)

summary(ann$pcaSlpNewlynStns)

# EOF1
ann$PCA1 <- ann$pcaSlpNewlynStns$loadings[,1]
dim(ann$PCA1) <- c(7,5)

# EOF2
ann$PCA2 <- ann$pcaSlpNewlynStns$loadings[,2]
dim(ann$PCA2) <- c(7,5)

# EOF3
ann$PCA3 <- ann$pcaSlpNewlynStns$loadings[,3]
dim(ann$PCA3) <- c(7,5)

# Map
#quartz()
x11()
map("worldHires", xlim=c(-25,15), ylim=c(35,65), interior=F, fill=F, col="grey50", resolution=0)
map.axes()
#EOF1
contour(seq(from=-20,to=10,by=5),seq(from=40,to=60,by=5), ann$PCA1, lwd=2, add=T)
#EOF2
contour(seq(from=-20,to=10,by=5),seq(from=40,to=60,by=5), ann$PCA2, lty='dotdash', add=T)
#EOF3
contour(seq(from=-20,to=10,by=5),seq(from=40,to=60,by=5), ann$PCA3, lty='dotted', add=T)

# Time series pattern
#quartz()
x11()
plot(ann$hadSLP2rAnnualTime,ann$pcaSlpNewlynStns$scores[,1], type='l', col='red')
lines(ann$hadSLP2rAnnualTime,ann$pcaSlpNewlynStns$scores[,2], col='blue')
lines(ann$hadSLP2rAnnualTime,ann$pcaSlpNewlynStns$scores[,3], col='magenta')
