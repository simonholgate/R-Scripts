# Look at change in rates over the period
meanRates177 <- array(NA,dim=c(4,177))
meanRates177[1,] <- 
  as.vector(mean(as.data.frame(s177PgrCorrectedRates),na.rm=TRUE))
# Calculate rates of first and second halves of record fixing the intercept
# point between them, using matrices
A <- matrix(data=0, nrow=55, ncol=3)
# Fill matrix with ones in the first column (the intercept). The second and
# third columns represent the first and second half rates, respectively
A[,1] <- 1
A[1:27,2] <- -27:-1
A[29:55,3] <- 1:27
for ( i in 1:177 ){
# We must delete any rows where data is missing
  j <- which(is.finite(stns177Annual[,i]))
  Atemp <- A[j,]
  junk <- solve(t(Atemp)%*%Atemp)%*%t(Atemp)%*%stns177Annual[j,i]
  meanRates177[2,i] <- junk[2,1]
  meanRates177[3,i] <- junk[3,1]
}
meanRates177[4,] <- meanRates177[3,] - meanRates177[2,]
# Remove La Libertad II (137) and Ko Taphao Noi (81) which have broken records
# Also remove FORT PHRACHULA CHOMKLAO and KUSHIRO which have broken records
meanRates177[,c(81,137,83,91)] <- NA
# Average into latitude bands
#zonalLat <- seq(from=-75, to=75, by=15)
zonalLat <- seq(from=-70, to=70, by=20)
lenLat <- length(zonalLat)
zonalMeanRates <- vector(mode="numeric", length=lenLat)
zonalMeanAccel <- vector(mode="numeric", length=lenLat)
zonalMeanRatesSE <- vector(mode="numeric", length=lenLat)
zonalMeanAccelSE <- vector(mode="numeric", length=lenLat)
for ( i in 1:lenLat ){
  j <- intersect(
    which(stns177Lat>(zonalLat[i]-10)), which(stns177Lat<=(zonalLat[i]+10)))
#    which(stns177Lat>(zonalLat[i]-7.5)), which(stns177Lat<=(zonalLat[i]+7.5)))
  zonalMeanRates[i] <- mean(meanRates177[1,j], na.rm=TRUE)
  if(length(j)>0){
    zonalMeanRatesSE[i] <- sd(meanRates177[1,j], na.rm=TRUE)/sqrt(length(j))
    zonalMeanAccelSE[i] <- sd(meanRates177[4,j], na.rm=TRUE)/sqrt(length(j))
  }
  zonalMeanAccel[i] <- mean(meanRates177[4,j], na.rm=TRUE)
}
#
plot(zonalLat, zonalMeanRates, ylim=c(0,4), xlim=c(-50,50), ann=FALSE)
title(main='Zonal mean rates of sea level change from 175 stations', 
  ylab='Rate of sea level change [mm/yr]', xlab='Latitude [deg]')
arrows(zonalLat, zonalMeanRates-zonalMeanRatesSE, 
  zonalLat, zonalMeanRates+zonalMeanRatesSE, col="blue", lwd=2, 
  code=3, angle=90, length=.05)
grid()

x11()
plot(zonalLat, zonalMeanAccel, ylim=c(-3,3), xlim=c(-50,50), ann=FALSE)
title(main='Zonal mean rates of sea level acceleration from 175 stations', 
  ylab='Rate of sea level acceleration [mm/yr]', xlab='Latitude [deg]')
arrows(zonalLat, zonalMeanAccel-zonalMeanAccelSE, 
  zonalLat, zonalMeanAccel+zonalMeanAccelSE, col="blue", lwd=2, 
  code=3, angle=90, length=.05)
grid()
#x11()
#plot(stns177Lat, meanRates177[1,],xlim=c(-90,90))
#x11()
#plot(stns177Lat, meanRates177[4,],xlim=c(-90,90))
