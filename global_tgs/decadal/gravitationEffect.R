# Look at change in rates over the period
meanRates177 <- array(NA,dim=c(4,177))
meanRates177[1,] <- 
  as.vector(mean(as.data.frame(s177PgrCorrectedRates),na.rm=TRUE))
meanRates177[2,] <-
  as.vector(mean(as.data.frame(s177PgrCorrectedRates[1:22,]),na.rm=TRUE))
meanRates177[3,] <-
  as.vector(mean(as.data.frame(s177PgrCorrectedRates[23:44,]),na.rm=TRUE))
meanRates177[4,] <- meanRates177[3,] - meanRates177[2,]
# Remove La Libertad II (137) and Ko Taphao Noi (81) which have broken records
meanRates177[,c(81,137)] <- NA
# Average into latitude bands
#zonalLat <- seq(from=-75, to=75, by=15)
zonalLat <- seq(from=-80, to=80, by=20)
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
plot(zonalLat, zonalMeanRates, ylim=c(0,4), ann=FALSE)
title(main='Zonal mean rates of sea level change from 175 stations', 
  ylab='Rate of sea level change [mm/yr]', xlab='Latitude [deg]')
arrows(zonalLat, zonalMeanRates-zonalMeanRatesSE, 
  zonalLat, zonalMeanRates+zonalMeanRatesSE, col="blue", lwd=2, 
  code=3, angle=90, length=.05)

x11()
plot(zonalLat, zonalMeanAccel, ylim=c(-3,3), ann=FALSE)
title(main='Zonal mean rates of sea level acceleration from 175 stations', 
  ylab='Rate of sea level acceleration [mm/yr]', xlab='Latitude [deg]')
arrows(zonalLat, zonalMeanAccel-zonalMeanAccelSE, 
  zonalLat, zonalMeanAccel+zonalMeanAccelSE, col="blue", lwd=2, 
  code=3, angle=90, length=.05)

