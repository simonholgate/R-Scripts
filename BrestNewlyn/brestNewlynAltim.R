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
# intersect(which(lon>=(360-5.7)),which(lon<=(360-5.3)))
# intersect(which(lat>=50.0),which(lat<=50.2))
# Newlyn is at x=1063, y=632

# Brest 48 23 N  04 30 W => 48.3833 -4.5
# intersect(which(lon>=(360-4.7)),which(lon<=(360-4.3)))
# intersect(which(lat>=48.2),which(lat<=48.5))
# Brest is at x=1066, y=625

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
# Focus on the Bay of Biscay area
lonR <- c(1030:1080,1:30)
latR <- 580:670
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
#Znew <- Z
#Znew[34,53] <- -50
#Znew[37,46] <- -50
#image.plot(lons, lat[latR], Znew, zlim=c(-25,55))

#stop("Imported data")

##################################
# Extract TS at Brest and Newlyn #
# Newlyn is at x=1063, y=632     #
# Brest is at x=1066, y=625      #
# Both need to be converted to   #
# the new grid co-ords           #
##################################
newlynAltim <- dataArray[which(lonR==1063), which(latR==632),]
brestAltim <- dataArray[which(lonR==1066), which(latR==625),]
x11()
par(family="HersheySans")
plot(1:624,brestAltim, type='l', col='red')
lines(1:624,newlynAltim, col='blue')

# Calculate rates
altimRates <- vector(length=2, mode='numeric')
junk <- lm(newlynAltim ~ c(1:624))
altimRates[1] <- junk$coef[2]*52
junk <- lm(brestAltim ~ c(1:624))
altimRates[2] <- junk$coef[2]*52

# Filter the altimetry with 53 week triangular filter
trian53ptFilter<-c(1:27,seq(from=26,to=1,by=-1))
trian53ptFilter<-trian53ptFilter/sum(trian53ptFilter)
filtNewlynAltim<-filter(newlynAltim,trian53ptFilter,method="c",sides=2)
filtBrestAltim<-filter(brestAltim,trian53ptFilter,method="c",sides=2)

filtNewlynAltim<-filtNewlynAltim-mean(filtNewlynAltim, na.rm=T)
filtBrestAltim<-filtBrestAltim-mean(filtBrestAltim, na.rm=T)

x11()
par(family="HersheySans")
plot(time,filtNewlynAltim, col='orange', type='l', lwd=2, ylim=c(-5.5,7), ann=F)
title(ylab="Sea Level [mm]", xlab="Year", main="Altimetry at Brest and Newlyn")	
lines(time,filtBrestAltim, col='blue', lwd=2)
legend(x=as.Date("1992/7/1"), y=7, 
  legend=c('Newlyn','Brest'),
  col=c('orange', 'blue'), lwd=2, pch=21)

stop("Extracted data")

# Tide gauges
oracle<-new.env()
load("newlynBrestTG.RData", envir=oracle)

# Calculate rates
tgRates <- vector(length=2, mode='numeric')
junk <- lm(oracle$newlynMonthly ~ c(1:1116))
tgRates[1] <- junk$coef[2]*12
junk <- lm(oracle$brestMonthly ~ c(1:1116))
tgRates[2] <- junk$coef[2]*12
#> tgRates
#[1] 1.724267 1.235732

oracle$time <- seq.Date(from=as.Date("1914/1/15"), to=as.Date("2006/12/15"), by="1 month")

x11()
par(family="HersheySans")
plot(oracle$time, oracle$newlynMonthly, type='l', col='orange', ylab="[mm]")
lines(oracle$time, oracle$brestMonthly, col='blue')

oracle$post1992 <- which(oracle$time >= as.Date("1992/10/1"))
oracle$seqPost1992 <- c(1:length(oracle$post1992))

# Calculate rates
tgRates1992 <- vector(length=2, mode='numeric')
lm.newlyn <- lm(oracle$newlynMonthly[oracle$post1992] ~ oracle$seqPost1992)
tgRates1992[1] <- lm.newlyn$coef[2]*12
lm.brest <- lm(oracle$brestMonthly[oracle$post1992] ~ oracle$seqPost1992)
tgRates1992[2] <- lm.brest$coef[2]*12
#> tgRates1992
#[1] 4.2104058 0.7185255
x11()
par(family="HersheySans")
# abline doesn't work here as we're using a Date object
plot(oracle$time[oracle$post1992], oracle$newlynMonthly[oracle$post1992], type='l', col='orange', ylab="[mm]")
lines(oracle$time[oracle$post1992[!is.na(oracle$newlynMonthly[oracle$post1992])]], fitted(lm.newlyn), lty=2, col='orange')
lines(oracle$time[oracle$post1992], oracle$brestMonthly[oracle$post1992], col='blue')
lines(oracle$time[oracle$post1992], fitted(lm.brest), lty=2, col='blue')

# Explore rates post October 1952 when Brest restarted
oracle$post1952 <- which(oracle$time >= as.Date("1952/10/1"))
oracle$seqPost1952 <- c(1:length(oracle$post1952))

# Calculate rates
tgRates1952 <- vector(length=2, mode='numeric')
lm.newlyn1952 <- lm(oracle$newlynMonthly[oracle$post1952] ~ oracle$seqPost1952)
tgRates1952[1] <- lm.newlyn1952$coef[2]*12
lm.brest1952 <- lm(oracle$brestMonthly[oracle$post1952] ~ oracle$seqPost1952)
tgRates1952[2] <- lm.brest1952$coef[2]*12

#> tgRates1952
#[1] 1.716404 1.671995

x11()
par(family="HersheySans")
# abline doesn't work here as we're using a Date object
plot(oracle$time[oracle$post1952], oracle$newlynMonthly[oracle$post1952], type='l', col='orange', ylab="[mm]")
lines(oracle$time[oracle$post1952[!is.na(oracle$newlynMonthly[oracle$post1952])]], fitted(lm.newlyn1952), lty=2, col='orange')
lines(oracle$time[oracle$post1952], oracle$brestMonthly[oracle$post1952], col='blue')
lines(oracle$time[oracle$post1952[!is.na(oracle$brestMonthly[oracle$post1952])]], fitted(lm.brest1952), lty=2, col='blue')

# Explore rates pre January 1940
oracle$pre1940 <- intersect(which(oracle$time >= as.Date("1915/5/15")), 
  which(oracle$time < as.Date("1940/1/1")))
oracle$seqpre1940 <- c(1:length(oracle$pre1940))

# Calculate rates
tgRates1940 <- vector(length=2, mode='numeric')
lm.newlyn1940 <- lm(oracle$newlynMonthly[oracle$pre1940] ~ oracle$seqpre1940)
tgRates1940[1] <- lm.newlyn1940$coef[2]*12
lm.brest1940 <- lm(oracle$brestMonthly[oracle$pre1940] ~ oracle$seqpre1940)
tgRates1940[2] <- lm.brest1940$coef[2]*12

#> tgRates1940
#[1] 2.5869099 0.8822277

x11()
par(family="HersheySans")
# abline doesn't work here as we're using a Date object
plot(oracle$time[oracle$pre1940], oracle$newlynMonthly[oracle$pre1940], type='l', col='orange', ylab="[mm]")
lines(oracle$time[oracle$pre1940[!is.na(oracle$newlynMonthly[oracle$pre1940])]], fitted(lm.newlyn1940), lty=2, col='orange')
lines(oracle$time[oracle$pre1940], oracle$brestMonthly[oracle$pre1940], col='blue')
lines(oracle$time[oracle$pre1940[!is.na(oracle$brestMonthly[oracle$pre1940])]], fitted(lm.brest1940), lty=2, col='blue')

# For comparison, 15/1/1920-15/12/1939 gives rates of:
#[1] 3.101257 2.214049

# Explore rates for comparison with Marcos et al, 2007
oracle$post1993 <- intersect(which(oracle$time >= as.Date("1993/1/1")),which(oracle$time <= as.Date("2002/9/1"))) 
oracle$seqPost1993 <- c(1:length(oracle$post1993))

# Calculate rates
tgRates1993 <- vector(length=2, mode='numeric')
lm.newlyn1993 <- lm(oracle$newlynMonthly[oracle$post1993] ~ oracle$seqPost1993)
tgRates1993[1] <- lm.newlyn1993$coef[2]*12
lm.brest1993 <- lm(oracle$brestMonthly[oracle$post1993] ~ oracle$seqPost1993)
tgRates1993[2] <- lm.brest1993$coef[2]*12

#> tgRates1993
#[1] 3.166487 2.873809

x11()
par(family="HersheySans")
# abline doesn't work here as we're using a Date object
plot(oracle$time[oracle$post1993], oracle$newlynMonthly[oracle$post1993], type='l', col='orange', ylab="[mm]")
lines(oracle$time[oracle$post1993[!is.na(oracle$newlynMonthly[oracle$post1993])]], fitted(lm.newlyn1993), lty=2, col='orange')
lines(oracle$time[oracle$post1993], oracle$brestMonthly[oracle$post1993], col='blue')
lines(oracle$time[oracle$post1993[!is.na(oracle$brestMonthly[oracle$post1993])]], fitted(lm.brest1993), lty=2, col='blue')

#####
# Plots of difference in rates between Newlyn and Brest
#####

# Read Delfzijl data
delfzijl.mm <- read.table("150011Delzijl.mm.txt", col.names=c("Year", "Height"), 
  colClasses=c("numeric","numeric"))

dDays <- round((delfzijl.mm$Year %% 1)*365)
dYrs <- floor(delfzijl.mm$Year)
dDate <- strptime(paste(dYrs, dDays), "%Y %j")

# Explore rates post October 1952 when Brest restarted
dPost1952 <- which(as.Date(dDate) >= as.Date("1952/10/1"))
dSeqPost1952 <- c(1:length(dPost1952))

# Calculate rates
lm.delfzijl1952 <- lm(delfzijl.mm$Height[dPost1952] ~ dSeqPost1952)
tgRates1952[3] <- lm.delfzijl1952$coef[2]*12

# Explore rates pre January 1940
dPre1940 <- intersect(which(as.Date(dDate) >= as.Date("1915/5/15")), 
  which(as.Date(dDate) < as.Date("1940/1/1")))
dSeqpre1940 <- c(1:length(dPre1940))

# Calculate rates
lm.delfzijl1940 <- lm(delfzijl$Height[dPre1940] ~ dSeqpre1940)
tgRates1940[3] <- lm.delfzijl1940$coef[2]*12

x11()
par(family="HersheySans")
plot(oracle$time[oracle$pre1940], (delfzijl$Height[dPre1940]+7000), col='magenta', type='l', ylab="[mm]")
lines(oracle$time[oracle$pre1940], oracle$newlynMonthly[oracle$pre1940], type='l', col='orange')
lines(oracle$time[oracle$pre1940[!is.na(oracle$newlynMonthly[oracle$pre1940])]], fitted(lm.newlyn1940), lty=2, col='orange')
lines(oracle$time[oracle$pre1940], oracle$brestMonthly[oracle$pre1940], col='blue')
lines(oracle$time[oracle$pre1940[!is.na(oracle$brestMonthly[oracle$pre1940])]], fitted(lm.brest1940), lty=2, col='blue')
lines(oracle$time[oracle$pre1940[!is.na(delfzijl$Height[dPre1940])]], fitted(lm.delfzijl1940)+7000, lty=2, col='magenta')

# Produce annual means from the data
junk <- oracle$newlynMonthly
dim(junk) <- c(12,93)
newlynAnnual <- colMeans(junk, na.rm=T)

junk <- oracle$brestMonthly
dim(junk) <- c(12,93)
brestAnnual <- colMeans(junk, na.rm=T)

junk <- delfzijl.mm$Height[which(as.Date(dDate) >= as.Date("1914/1/1"))]
dim(junk) <- c(12,93)
delfzijlAnnual <- colMeans(junk, na.rm=T)

#
annualYrs = seq.Date(from=as.Date("1914/1/1"),to=as.Date("2006/12/1"), by="1 year")
x11()
par(family="HersheySans")
plot(annualYrs, delfzijlAnnual-mean(delfzijlAnnual, na.rm=T), col='magenta', type='l', ylab="[mm]")
lines(annualYrs, newlynAnnual-mean(newlynAnnual, na.rm=T), col='orange')
lines(annualYrs, brestAnnual-mean(brestAnnual, na.rm=T), col='blue')

# SLP correction
load("newlynStnsSlpRates.RData")
tgRates1996 <- tgRates1996 + stnsSlpRates1996
tgRates <- tgRates + stnsSlpRates
#> tgRates
#[1] 1.725052 1.236479
#> tgRates1996
#[1] 4.2277234 0.7303179
#> altimRates
#[1] 0.1740214 0.3086302

# CGPS rates (from Simon Williams)
cgpsRates <- vector(length=2, mode='numeric')
# Newlyn -0.21mm/yr +/-0.4
cgpsRates[1] <- -0.21
# Brest -0.54mm/yr error > +/-0.4
cgpsRates[2] <- -0.54
#
#> tgRates+cgpsRates
#[1] 1.5150523 0.6964787

# Model data
# Filter the altimetry with 53 week triangular filter
trian13ptFilter<-c(1:7,seq(from=6,to=1,by=-1))
trian13ptFilter<-trian13ptFilter/sum(trian13ptFilter)
filtNewlynModel<-filter(newlynMonthlyMean,trian13ptFilter,method="c",sides=2)
filtNewlynModel<-filtNewlynModel-mean(filtNewlynModel, na.rm=T)
filtBrestModel<-filter(brestMonthlyMean,trian13ptFilter,method="c",sides=2)
filtBrestModel<-filtBrestModel-mean(filtBrestModel, na.rm=T)

x11()
par(family="HersheySans")
plot(monthsArray,filtNewlynModel, col='orange', type='l', lwd=2, ylim=c(-0.1,0.1), ann=F)
title(ylab="Sea Level [mm]", xlab="Year", main="Model Data from Brest and Newlyn")	
lines(monthsArray,filtBrestModel, col='blue', lwd=2)
legend(x=as.Date("1960/7/1"), y=0.1, 
  legend=c('Newlyn','Brest'),
  col=c('orange', 'blue'), lwd=2, pch=21)
stop("Finished")

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
