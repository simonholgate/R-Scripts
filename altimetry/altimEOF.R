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

#l<-1080
#m<-915

# 1 deg grid
lon <- seq(from=1, to=360)
lat <- seq(from=-90, to=90)

l<-360
m<-180

lm<-l*m

#n<-639
#n<-513
# Make n a factor of 4 so that we can reduce the array into pseudo months
#n<-512
n<-52

# Just worry about n time slices to begin with
# Focus on the Pacific covering [105:255,90:155]
l1 <- 80
m1 <- 50
lm1 <- l1*m1
dataArray<-array(NA,dim=c(l1,m1,n))

# Read the filtered time series that CWH has produced and see which grid
# points have only null values (=9999.0). These can then be flagged as
# land/ice points and removed from the time series.
# The format of hstack_resid (which has annual and semi annual signals removed
# - see the email Altimetry.txt) has time as the first axis, longitude as the
# second and latitude as the third.

# We'll use hmaps.dat in place of hstack.dat for EOFs
#inConn <- file("~/diskx/altimetry/hmaps.dat", "rb")
inConn <- file("~/diskx/altimetry/hmaps_sm5_1deg.dat", "rb")
for (i in 1:n) {
  data<-readBin(inConn, what="numeric", n = lm, size = 4, endian='big')
  dim(data)<-c(l,m)
  dataArray[,,i] <- data[145:224,90:139]
}
close(inConn)

library(fields)

Z <- data[145:224,90:139]
Z[which(Z==9999)] <- NA
par(family="HersheySans")
image.plot(lon[145:224], lat[90:139], Z)

############
# EOF part #
############

# Flatten 3D array to 2D
dim(dataArray) <- c(lm1,n)
dataArray <- rotate270.matrix(dataArray)
monDataArr <- dataArray
# Data is every 7 days. Form pseudo months by averaging
# 4 weeks together
#monDataArr <- array(NA, dim=c((n/4),lm1))
#for (i in 1:lm1) {
#  junk <- dataArray[,i]
#  dim(junk) <- c(4,(n/4))
#  monDataArr[,i] <- colMeans(junk)
#} 

# Remove mean from each column
meanDA <- colMeans(monDataArr)

# land is where the rowMeans are 9999
sea <- which(meanDA!=9999)
lensea <- length(sea)
for (i in 1:lm1) {
  monDataArr[,i] <- monDataArr[,i]-meanDA[i]
}

# Reduce the array size by removing land points
monDataArr<-monDataArr[,sea]

# Form covariance matrix
covMat <- t(monDataArr) %*% monDataArr

# Calculate eigenvalues and eigenvectors of the covariance matrix
ev <- eigen(covMat, symmetric=TRUE)

# First ncomp eof components
ncomp<-10
PC <- array(NA, dim=c(n,lensea))
eof <- array(NA, dim=c(n,lensea,ncomp))
for (i in 1:ncomp){
  PC[,i] <- monDataArr %*% ev$vectors[,i]
  eof[,,i] <- PC[,i] %*% t(ev$vectors[,i])
}

# Put back the sea
fullEOFarr <- array(NA, dim=c(n,lm1, ncomp))
fullEOFarr[,sea,] <- eof

eof1 <- fullEOFarr[,,1]
eof1 <- rotate90.matrix(eof1)
dim(eof1) <- c(l1,m1,n)

eofNow <- fullEOFarr[1,,]
eofNow <- rotate90.matrix(eofNow)
eofNow <- colSums(eofNow)
dim(eofNow) <- c(l1,m1,n)

x11()
par(family="HersheySans")
image.plot(lon[145:224], lat[90:139], eof1[,,1])
# Amount of variance explained
require('stats')
ame <- diag(ev$values)/sum(diag(ev$values))

