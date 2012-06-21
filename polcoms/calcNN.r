#rm(list=ls())
load("globalMeans.RData")

plot(ts(globalSLPMonthlyMean,start=c(1960,1),freq=12))

# SST anomalies
SSTMat <- matrix(globalSSTMonthlyMean,nrow=12)
SSTClim <- apply(SSTMat,1,mean)
aSSTMat <- SSTMat - SSTClim
aSST <- matrix(aSSTMat)

# SSH anomalies
SSHMat <- matrix(globalSSHMonthlyMean,nrow=12)
SSHClim <- apply(SSHMat,1,mean)
aSSHMat <- SSHMat - SSHClim
aSSH <- matrix(aSSHMat)

# SLP anomalies
SLPMat <- matrix(globalSLPMonthlyMean,nrow=12)
SLPClim <- apply(SLPMat,1,mean)
aSLPMat <- SLPMat - SLPClim
aSLP <- matrix(aSLPMat)

# WindE anomalies
WindEMat <- matrix(globalWindEMonthlyMean,nrow=12)
WindEClim <- apply(WindEMat,1,mean)
aWindEMat <- WindEMat - WindEClim
aWindE <- matrix(aWindEMat)

# WindN anomalies
WindNMat <- matrix(globalWindNMonthlyMean,nrow=12)
WindNClim <- apply(WindNMat,1,mean)
aWindNMat <- WindNMat - WindNClim
aWindN <- matrix(aWindNMat)

# time series plot
ts.obj <- data.frame('aSST'=aSST,'aSSH'=aSSH,
			'aSLP'=aSLP,'aWindE'=aWindE,
			'aWindN'=aWindN)

plot(ts(ts.obj, start=c(1960,1), freq=12),
			main="Time series plot")

mycors <- cor(ts.obj)

# Use NN transfer function
# Predictors
myX <- data.frame('aSST'=aSST,
			'aSLP'=aSLP,'aWindE'=aWindE,
			'aWindN'=aWindN)
# Predictands
myY <- data.frame('aSSH'=aSSH)

# These are NN parameters
nNodes <- 3
nIter <- 600 # Number iterations
step<-10 # reporting steps
errc<-rep(0,nIter/step) # Calibration error
errv<-errc # Validation error
nSamp <- dim(myY)[1] # total number of samples
nCal<-nSamp*4/5 # Size of calibration set
vIter <- seq(0,nIter,by=step)

# standardization

yMin<-apply(myY,2,min)
yMax<-apply(myY,2,max)
xMin<-apply(myX,2,min)
xMax<-apply(myX,2,max)
Y<-scale(myY,center=yMin,scale=(yMax-yMin))
X<-scale(myX,center=xMin,scale=(xMax-xMin))
#Z<-scale(myZ,center=xMin,scale=(xMax-xMin)) # For historical predictors

# validation and calibration set for NN
samp <- sample(1:nSamp,nCal)
Xcal=X[samp,]
Xval=X[-samp,]
Ycal=Y[samp,]
Yval=Y[-samp,]

# initialise NN
dd.nn <- nnet(Xcal, Ycal, size = nNodes, decay = 0.001, maxit = 2, 
		trace=F,linout=T) 

op <- par(mfrow=c(2,2))

for(i in 1:(nIter/step) )
{
# iterate network
		dd.nn <- nnet(Xcal, Ycal, size = nNodes, decay = 0.00001, 
				maxit = step, Wts=dd.nn$wts,trace=FALSE,linout=T)

# predict on calibration set
		Yec<-predict(dd.nn,Xcal,type="raw")
		plot(c(as.matrix(Yec)),c(as.matrix(Ycal)), type='p',xlab='pred', ylab='obs', main='Calibration')
# predict on validation set
		Yev<-predict(dd.nn,Xval,type="raw")
		plot(c(as.matrix(Yev)),c(as.matrix(Yval)), type='p',xlab='pred', ylab='obs', main='Verification')

# Get mean squared errors
		errc[i]<-mean((Yec-Ycal)^2)
		errv[i]<-mean((Yev-Yval)^2)
		plot(vIter[1:i],errc[1:i],type='l', xlab='Iteration', ylab='error',main='Calibration',ylog=T)
		plot(vIter[1:i],errv[1:i],type='l', xlab='Iteration', ylab='error',main='Verification',ylog=T)	
}

# find error minima
cat("errc min at iter =",which.min(errc)*step,"  value = ",min(errc),"\n")
cat("errv min at iter =",which.min(errv)*step,"  value = ",min(errv),"\n")
cat("errc+errv min =",which.min(errc+errv)*step,"  value = ",min(errc+errv),"\n")
iterstop<-which.min(errc+errv)*step

# rerun nnet to minimum obtained in previous step
dd.nn <- nnet(X, Y, size = nNodes, decay = 0.00001, maxit = iterstop, trace=T,linout=T)
yp<-predict(dd.nn,X,type="raw")
summary(lm(yp ~ Y))

par(op)
