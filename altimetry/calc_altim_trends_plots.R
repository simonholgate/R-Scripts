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

op <- par(no.readonly = TRUE)
par(lab=c(8,6,7))
plot(zonalLat, zonalMeanLat, lwd=2, xlim=c(-50,50), ann=FALSE)
#title(main='Zonal mean rates of sea level change from altimetry',
title(ylab='Rate of sea level change [mm/yr]', xlab='Latitude [deg]')
grid(col="black")
dev2bitmap(file='zonalRatesSLCAltimetry.png')

x11()
par(lab=c(5,6,7))
#image.plot(lon,lat,trendArray,zlim=c(-15,15))
library(fields)
image.plot(lon,lat,trendArray,zlim=c(-50,50))
dev2bitmap(file='ratesSLCAltimetry.png')

par(lab=c(8,6,7))

x11()
plot(zonalLat, zonalMeanLatMass, lwd=2, xlim=c(-50,50), ann=FALSE)
#title(main='Zonal mean rates of sea level change due to mass (from altimetry)',
title(ylab=expression(
  paste('Rate of sea level change due to mass [mm ',yr^{-1},']')), 
  xlab=expression(paste('Latitude [',degree,']')))
grid(col="black")
dev2bitmap(file='zonalRatesSLMassCAltimetry.png')

x11()
plot(zonalLat, zonalPacMeanLat, lwd=2, xlim=c(-50,50), ann=FALSE, pch=2)
points(zonalLat, zonalAtlMeanLat, lwd=2, ann=FALSE, pch=1)
title(main='Zonal rates of sea level change (Atlantic and Pacific) (from altimetry)')
title(
  ylab=expression(paste('Rate of sea level change [mm ',yr^{-1},']')), 
  xlab=expression(paste('Latitude [',degree,']')))
grid(col="black")
legend(-45,0.5,c("Atlantic","Pacific"),pt.lwd=2,pch=c(1,2))
dev2bitmap(file='zonalRatesSLCPacAtlAltimetry.png')

x11()
plot(zonalLat, zonalAtlMeanLat-zonalPacMeanLat, lwd=2, 
  xlim=c(-50,50), ann=FALSE)
title(main='Difference in zonal SLC rates (Atl - Pac) (from altimetry)')
title(
  ylab=expression(paste('Rate of sea level change [mm ',yr^{-1},']')), 
  xlab=expression(paste('Latitude [',degree,']')))
grid(col="black")
dev2bitmap(file='zonalRatesSLCDiffPacAtlAltimetry.png')

x11()
par(lab=c(5,6,7))
library(fields)
# Rates calculated are rates per 7 days so multiply by 52 weeks to get mm/yr
image.plot(lon,lat,trendArray1*52,zlim=c(-50,50))
dev2bitmap(file='ratesSLCAltimetry1.png')

x11()
par(lab=c(5,6,7))
image.plot(lon,lat,trendArray2*52,zlim=c(-50,50))
dev2bitmap(file='ratesSLCAltimetry2.png')

x11()
par(lab=c(5,6,7))
image.plot(lon,lat,(trendArray2-trendArray1)*52,zlim=c(-50,50))
dev2bitmap(file='diffRatesSLCAltimetry12.png')

par(op)
