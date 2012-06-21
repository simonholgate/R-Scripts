#######################
# ENACT dimensions:   #
#                     #
# Xmin:   1/3         #
# Xres:   1/3         #
# Xn:    1080         #
#                     #
# Ymin:  -82.0        #
# Yres:  varies w lat #
# Yn:     915         #
#                     #
# Timeslices: 639     #
# Time span:          #
# 14.10.92 - 5.1.05   #
# Time interval: 7 d  #
#######################

# Load matrixMethods for use later
source("~/workspace/RScripts/matrixMethods.R")
# Lerwick 60 09 N  01 08 W => 60.15 -1.1333
# intersect(which(lon>=(360-1.3)),which(lon<=(360-1.0)))
# intersect(which(lat>=60.0),which(lat<=60.2))
# Lerwick is at x=1077, y=685

# Wick 58 26 N  03 05 W => 58.4333 -3.0833
# intersect(which(lon>=(360-3.3)),which(lon<=(360-3.0)))
# intersect(which(lat>=58.3),which(lat<=58.5))
# Wick is at x=1071, y=675

# Invergordon 57 41 N  04 10 W => 57.6833 -4.1666
# intersect(which(lon>=(360-4.3)),which(lon<=(360-4.0)))
# intersect(which(lat>=57.5),which(lat<=57.7))
# 1068, 671 is NA so use
# Invergordon is at x=1066, y=671

# Stornoway 58 12 N  06 23 W => 58.2 -6.3833
# intersect(which(lon>=(360-6.5)),which(lon<=(360-6.2)))
# intersect(which(lat>=58.1),which(lat<=58.3))
# Stornoway is at x=1061, y=674

# Aberdeen I 57 09 N  02 05 W => 57.15 -2.0833
# intersect(which(lon>=(360-2.2)),which(lon<=(360-2.0)))
# intersect(which(lat>=57.1),which(lat<=57.2))
# Aberdeen is at x=1074, y=668

# Torshavn 62 01 N  06 46 W => 60.02 -6.77
# intersect(which(lon>=(360-6.9)),which(lon<=(360-6.5)))
# intersect(which(lat>=61.9),which(lat<=62.1))
# Torshavn is at x=1060, y=697
lerwickStnsLatLon <- c(60.15, -1.1333, 58.4333, -3.0833, 57.6833, 
  -4.1666, 58.2, -6.3833, 57.15, -2.0833, 60.02, -6.77)
dim(lerwickStnsLatLon) <- c(2,6)
# High res grid
lon<-seq(from=1/3, length=1080, by=1/3)
junk<-read.table('~/diskx/altimetry/lats.lis')
lat<-junk$V1
# jd (Julian Day) in this context appears to be days since 1950 
daytable<-read.table('~/diskx/altimetry/daytable.lis', 
  col.names=c("seq","jd","year","month","day","yearday"))

xres<-1/3
l<-1080
m<-915

# 1 deg grid
#xres<-1
#lon <- seq(from=1, to=360)
#lat <- seq(from=-90, to=90)

#l<-360
#m<-180

lm<-l*m

#n<-639
#n<-513
# Make n a factor of 4 so that we can reduce the array into pseudo months
#n<-512
n<-624
time<-seq(as.Date("1992/10/14"), by="7 days", length.out=n)

# Just worry about n time slices to begin with
lonR <- c(1030:1080,1:30)
latR <- 580:700
l1 <- length(lonR)
m1 <- length(latR)
lm1 <- l1*m1
lons <- seq(from=(lon[lonR[1]]-360),to=lon[lonR[l1]], by=xres)
dataArray<-array(NA,dim=c(l1,m1,n))

# Read the filtered time series that CWH has produced and see which grid
# points have only null values (=9999.0). These can then be flagged as
# land/ice points and removed from the time series.
# The format of hstack_resid (which has annual and semi annual signals removed
# - see the email Altimetry.txt) has time as the first axis, longitude as the
# second and latitude as the third.

# We'll use hmaps.dat in place of hstack.dat for EOFs
inConn <- file("~/diskx/altimetry/hmaps.dat", "rb") # little endian
#inConn <- file("~/diskx/altimetry/hmaps_sm5_1deg.dat", "rb") # big endian
for (i in 1:n) {
  data<-readBin(inConn, what="numeric", n = lm, size = 4, endian='little')
#  data<-readBin(inConn, what="numeric", n = lm, size = 4, endian='big')
  dim(data)<-c(l,m)
  dataArray[,,i] <- data[lonR,latR]
}
close(inConn)

library(fields)

Z <- data[lonR,latR]
Z[which(Z==9999)] <- NA

#postscript(file="nAtl.ps")
par(family="HersheySans")
image.plot(lons, lat[latR], Z)
world(add=T)
#dev.off()
#Z[34,53] <- -50
#Z[37,46] <- -50

#stop("Imported data")

##################################
# Extract TS at Lerwick, Wick,   #
# Invergordon and Stornoway      #
# Lerwick is at x=1077, y=685    #
# Wick is at x=1071, y=675       #
# Invergordon is at x=1066, y=671#
# Stornoway is at x=1061, y=674  #
# Aberdeen is at x=1074, y=668   #
# Torshavn is at x=1060, y=697   #
# Both need to be converted to   #
# the new grid co-ords           #
##################################

lerwickAltim <- dataArray[which(lonR==1077), which(latR==685),]
wickAltim <- dataArray[which(lonR==1071), which(latR==675),]
invergordonAltim <- dataArray[which(lonR==1066), which(latR==671),]
stornowayAltim <- dataArray[which(lonR==1061), which(latR==674),]
aberdeenAltim <- dataArray[which(lonR==1074), which(latR==668),]
torshavnAltim <- dataArray[which(lonR==1060), which(latR==697),]

x11()
par(family="HersheySans")
plot(time,stornowayAltim, col='magenta', type='l')
lines(time,lerwickAltim, col='red')
lines(time,wickAltim, col='blue')
lines(time,invergordonAltim, col='orange')
lines(time,aberdeenAltim, col='green')
lines(time,torshavnAltim, col='cyan')

# Calculate linear trends from the data
seqTime <- c(1:length(time))
altimRateArray <- vector(length=6, mode='numeric')
junk <- lm(lerwickAltim ~ seqTime)
altimRateArray[1] <- junk$coef[2]
junk <- lm(wickAltim ~ seqTime)
altimRateArray[2] <- junk$coef[2]
junk <- lm(invergordonAltim ~ seqTime)
altimRateArray[3] <- junk$coef[2]
junk <- lm(stornowayAltim ~ seqTime)
altimRateArray[4] <- junk$coef[2]
junk <- lm(aberdeenAltim ~ seqTime)
altimRateArray[5] <- junk$coef[2]
junk <- lm(torshavnAltim ~ seqTime)
altimRateArray[6] <- junk$coef[2]
#> time[2]-time[1]
#Time difference of 7 days
#> altimRateArray*52
#[1] 0.2208210 0.2098690 0.2011339 0.2152284 0.1314888 0.2215287
# Filter the altimetry with 53 week triangular filter
trian53ptFilter<-c(1:52,seq(from=51,to=1,by=-1))
trian53ptFilter<-trian53ptFilter/sum(trian53ptFilter)
filtLerwickAltim<-filter(lerwickAltim,trian53ptFilter,method="c",sides=2)
filtWickAltim<-filter(wickAltim,trian53ptFilter,method="c",sides=2)
filtInvergordonAltim<-filter(invergordonAltim,trian53ptFilter,method="c",sides=2)
filtStornowayAltim<-filter(stornowayAltim,trian53ptFilter,method="c",sides=2)
filtAberdeenAltim<-filter(aberdeenAltim,trian53ptFilter,method="c",sides=2)
filtTorshavnAltim<-filter(torshavnAltim,trian53ptFilter,method="c",sides=2)

x11()
par(family="HersheySans")
plot(time,filtStornowayAltim, col='orange', type='l', lwd=2, ylim=c(-3,5), ann=F)
title(ylab="Sea Level [mm]", xlab="Year", main="Altimetry at TGs Around N Scotland")	
lines(time,filtWickAltim, col='blue', lwd=2)
lines(time,filtInvergordonAltim, col='green', lwd=2)
lines(time,filtAberdeenAltim, col='magenta', lwd=2)
lines(time,filtTorshavnAltim, col='cyan', lwd=2)
lines(time,filtLerwickAltim, col='red', lwd=2)

legend(x=as.Date("1992/7/1"), y=5, 
  legend=c('Lerwick','Aberdeen I','Invergordon','Stornoway','Wick','Torshavn'),
  col=c('red','magenta','green', 'orange', 'blue', 'cyan'), lwd=2, pch=21)
  
stop("End of altimetry extraction")

##############################
# EOF part - uses SVD method #
##############################

# Flatten 3D array to 2D
dim(dataArray) <- c(lm1,n)
dataArray <- rotate270.matrix(dataArray)
# Replace missing value of 9999 with NA
dataArray[which(dataArray==9999)]<-NA
# For now, skip using monthly means
#monDataArr <- dataArray
# Data is every 7 days. Form pseudo months by averaging
# 4 weeks together
monDataArr <- array(NA, dim=c((n/4),lm1))
for (i in 1:lm1) {
  junk <- dataArray[,i]
# Remove mean seasonal cycle
  dim(junk) <- c(52,(n/52))
  junk <- junk - rowMeans(junk)
  dim(junk) <- c(n)

  dim(junk) <- c(4,(n/4))
  monDataArr[,i] <- colMeans(junk)
} 

# Remove mean from each column
meanDA <- colMeans(monDataArr)

# land is where the colMeans are NA
sea <- which(is.finite(meanDA))
lensea <- length(sea)
for (i in 1:lm1) {
  monDataArr[,i] <- monDataArr[,i]-meanDA[i]
}
latestMonDataArr <- monDataArr[(n/4),]
#latestMonDataArr <- rotate90.matrix(latestMonDataArr)
dim(latestMonDataArr) <- c(l1,m1)
latestMonDataArr<-flip.matrix(latestMonDataArr)
latestMonDataArr<-mirror.matrix(latestMonDataArr)
latestMonDataArr[34,53] <- -50
latestMonDataArr[37,46] <- -50

# Reduce the array size by removing land points
monDataArr<-monDataArr[,sea]

# Calculate the SVD
#S <- svd(monDataArr, nu=52, nv=4000)
S <- svd(monDataArr)
# S is a list of the componets formed from the SVD X = UDV'. 
# The components are named u,d,v

# First ncomp eof components
ncomp<-20
PC <- array(NA, dim=c((n/4),ncomp))
eof <- array(NA, dim=c((n/4),lensea,ncomp))
for (i in 1:ncomp){
  PC[,i] <- monDataArr %*% S$v[,i]
  eof[,,i] <- PC[,i] %*% t(S$v[,i])
}

# Put back the sea
fullEOFarr <- array(NA, dim=c((n/4),lm1, ncomp))
fullEOFarr[,sea,] <- eof

eof1 <- fullEOFarr[,,1]
eof1 <- rotate90.matrix(eof1)
dim(eof1) <- c(l1,m1,(n/4))

x11()
par(family="HersheySans")
image.plot(lons, lat[latR], eof1[,,1])
world(add=T)
# Reconstruction of latest map
eofLatest <- fullEOFarr[(n/4),,]
eofLatest <- rotate270.matrix(eofLatest)
eofLatest <- colSums(eofLatest)
dim(eofLatest) <- c(l1,m1)

x11()
par(family="HersheySans")
image.plot(lons, lat[latR], eofLatest, zlim=c(-22,42))
world(add=T)

x11()
par(family="HersheySans")
image.plot(lons, lat[latR], latestMonDataArr, zlim=c(-22,42))
world(add=T)

# Amount of variance explained
require('stats')
eigVals <- S$d^2/lensea
ame <- eigVals*100/sum(eigVals)

stop("EOFs completed")
#
# Cost function, following Church et al, J Clim, 2004, and Kaplan et al, J
# Clim, 2000
#
# Equation to be minimised is:
# S(alpha) = t(K*Ur*alpha - Ho)*M(-1)*(K*Ur*alpha - Ho) + t(alpha)*Lambda*alpha
#
# where:
# Ho is the tide gauge observations, K is the sampling vector defined to be
# equal to 1 where the TG data are available and 0 otherwise, Lambda is the
# diagonal matrix of the eigenvalues of the covariance matrix.
# M is the error covariance matrix, given by:
# M = R + K*U'*Lambda'*t(U')*t(K)
# Here, R is the variace of the instrumental error, and the 2nd term of the
# RHS repressents errors in omission introduced by deleting higher order EOFs.
# The primes denote the truncated matrices
# Following Church et al, it is assumed that the instrumental variance is 4mm

# Select synthetic TG sites as the N & S perimeter of the region selected
#lonR <- 145:205
#latR <- 90:140

# K has the same size as monDataArr (Ho in the above)
K <- array(0,dim(monDataArr))
for (i in 1:l1){
	K[,i]<-1
}
for (i in (lm1-l1):lm1){
	K[,i]<-1
}

TG <- K*monDataArr

# P denotes the truncated or "prime" matrix
Lambda <- array(0, dim=c(n/4,n/4))
LambdaP <- array(0, dim=c(n/4,n/4))

for (i in 1:ncomp){
  Lambda[i,i] <- eigVals[i]
}
for (i in (ncomp+1):(n/4)){
  LambdaP[i,i] <- eigVals[i]
}

V <- array(0, dim(S$v))
VP <- array(0, dim(S$v))

V[,1:ncomp] <- S$v[,1:ncomp]
VP[,(ncomp+1):(n/4)] <- S$v[,(ncomp+1):(n/4)]

R <- diag(4^2/(n/4),nrow=lensea)

M <- R + t(K)*VP%*%LambdaP%*%t(VP)*K

# S(alpha) = t(K*Ur*alpha - Ho)*M(-1)*(K*Ur*alpha - Ho) + t(alpha)*Lambda*alpha
Salp <- array(NA,dim=dim(PC))
for (i in 1:(n/4)){
  Salp[i,] <- t(K%*%V%*%PC[i,] - TG)%*%solve(M)%*%(K%*%V%*%PC[i,] - TG) 
  	+ t(PC[i,])%*%solve(Lambda)%*%PC[i,]
}
