# Script to plot the correlation between two points against their great circle
# distance from each other. This is to identify teleconnections i.e. high
# correlations which are far apart from each other.
gcd <- array(NA,dim=c(9,9))
for (i in 1:9) {
  for (j in 1:9) {
    gcd[i,j] <- greatcircle(nineStationsLatLon[i,1],nineStationsLatLon[i,2],
      nineStationsLatLon[j,1], nineStationsLatLon[j,2])
  }
}
plot(abs(as.vector(gcd)/1000), as.vector(pgrCorrelationMatrix[1:9,1:9]), ann=F)
title(main="Correlation coefficient versus Great Circle Distance", 
  xlab="Distance [km]", ylab="Correlation Coefficient between pairs")

grid()
