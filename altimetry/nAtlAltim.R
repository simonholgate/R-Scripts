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
#                     #
#######################

# Load matrixMethods for use later
source("~/workspace/RScripts/matrixMethods.R")
# Newlyn 50 06 N  05 33 W => 50.1 -5.55
# Newlyn is at x=86, y=90
#> (-5.55-xmin)/xres
#[1] 85.69998
#> (50.1-ymin)/yres
#[1] 89.90001

# Use jet colors for images
#source('~/bin/RScripts/jet.colors.R')

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
n<-468

# Just worry about n time slices to begin with
# Focus on the N Atlantic covering 100W-15E, 5S-80N : [(781:1080,1:45),443:877]
#lonR <- c(260:360,1:15)
#latR <- 85:170
lonR <- c(781:1080,1:45)
latR <- 443:877
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
#dev.off()
#stop("Imported data")

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
image.plot(lon[lonR], lat[latR], eof1[,,1])

# Reconstruction of latest map
eofLatest <- fullEOFarr[(n/4),,]
eofLatest <- rotate270.matrix(eofLatest)
eofLatest <- colSums(eofLatest)
dim(eofLatest) <- c(l1,m1)

x11()
par(family="HersheySans")
image.plot(lon[lonR], lat[latR], eofLatest)
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
