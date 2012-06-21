postscript(file="eAtlRates.ps")
op <- par()
par(oma=c( 0,0,3,5), family="HersheySans")
image(lats_e3043, 1955:1995, masked_e_atl_rates, ann=FALSE,
  col=tim.colors(10),
  zlim=c(-20,20))
#axis(3,at=c(200,400,600,800),
#  labels=c("-3.98","9.8","36.79","51.33"))
#  labels=c("50.59", "35.27", "8.18", "-5.62"))
title(xlab="Latitude", ylab="Year")
par(oma=c( 0,0,0,1))
image.plot(legend.only=TRUE, zlim=c(-20,20), legend.lab="mm/yr")
par <- op
title(main="Model Decadal Rates: East side of Atlantic")
dev.off()

postscript(file="wAtlRates.ps")
par(oma=c( 0,0,3,5), family="HersheySans")
image(1:length(mask_w[,1]), 1955:1995, masked_w_atl_rates, xlab="Station",
  ylab="Year", col=tim.colors(10),
  zlim=c(-20,20))
axis(3,at=c(200,400,600,800,1000,1200,1400),
  labels=c("-9.80", "8.64", "15.74", "27.06", "39.00", "48.30", "59.18"))
 # labels=c("51.03","44.42","29.54","20.63","10.72","-0.47","-23.46"))
title(ylab="Latitude")
par(oma=c( 0,0,0,1))
image.plot(legend.only=TRUE, zlim=c(-20,20), legend.lab="mm/yr")
par <- op
title(main="Model Decadal Rates: West side of Atlantic")
dev.off()

postscript(file="cAtlRates.ps")
par(oma=c( 0,0,3,5), family="HersheySans")
image(1:988, 1955:1995, masked_c_atl_rates, xlab="Station",
  ylab="Year", col=tim.colors(10),
  zlim=c(-20,20))
axis(3,at=c(200,400,600,800,1000,1200,1400),
  labels=c("-9.80", "8.64", "15.74", "27.06", "39.00", "48.30", "59.18"))
 # labels=c("51.03","44.42","29.54","20.63","10.72","-0.47","-23.46"))
title(ylab="Latitude")
par(oma=c( 0,0,0,1))
image.plot(legend.only=TRUE, zlim=c(-20,20), legend.lab="mm/yr")
par <- op
title(main="Model Decadal Rates: Central Atlantic")
dev.off()
#
# Using yr as the x axis and limiting latitudinal range to 30-43deg N to fit
# with Vassil

postscript(file="eAtlRates30_43.ps")
op <- par()
par(oma=c( 0,0,3,5))
image(1955:1995, e3043, 
  rotate270.matrix(masked_e_atl_rates[e3043,]), 
  col=tim.colors(10), ann=FALSE, axes=FALSE,
  zlim=c(-20,20))
axis(2,at=c(560,580,600,620,640),
  labels=c("32.35", "34.11", "36.79", "38.27", "40.98"))
axis(1,at=c(1960,1970,1980,1990),
  labels=c("1960", "1970", "1980", "1990"))
title(ylab="Latitude", xlab="Year")
par(oma=c( 0,0,0,1))
image.plot(legend.only=TRUE, zlim=c(-20,20), legend.lab="mm/yr")
par <- op
title(main="Model Decadal Rates: East side of Atlantic 30-43N")
dev.off()

postscript(file="wAtlRates30_43.ps")
par(oma=c( 0,0,3,5))
image(1955:1995, w3043, 
  rotate270.matrix(masked_w_atl_rates[w3043,]), 
  col=tim.colors(10), ann=FALSE, axes=FALSE,
  zlim=c(-20,20))
axis(2,at=c(920,940,960,980,1000,1020,1040),
  labels=c("30.15", "32.55", "34.31", "36.41", "39.00", "40.80", "41.86"))
axis(1,at=c(1960,1970,1980,1990),
  labels=c("1960", "1970", "1980", "1990"))
title(ylab="Latitude", xlab="Year")
par(oma=c( 0,0,0,1))
image.plot(legend.only=TRUE, zlim=c(-20,20), legend.lab="mm/yr")
par <- op
title(main="Model Decadal Rates: West side of Atlantic 30-43N")
dev.off()

postscript(file="cAtlRates30_43.ps")
par(oma=c( 0,0,3,5), family="HersheySans")
image(1955:1995, c3043, 
  rotate270.matrix(masked_c_atl_rates[c3043,]), 
  col=tim.colors(10), ann=FALSE, axes=FALSE,
  zlim=c(-20,20))
axis(2,at=c(920,940,960,980,1000,1020,1040),
  labels=c("30.15", "32.55", "34.31", "36.41", "39.00", "40.80", "41.86"))
axis(1,at=c(1960,1970,1980,1990),
  labels=c("1960", "1970", "1980", "1990"))
title(ylab="Latitude", xlab="Year")
par(oma=c( 0,0,0,1))
image.plot(legend.only=TRUE, zlim=c(-20,20), legend.lab="mm/yr")
par <- op
title(main="Model Decadal Rates: Central Atlantic 30-43N")
dev.off()

# Remove 1960-1990 mean from data to show the anomaly

postscript(file="eAtlAnom30_43.ps")
op <- par()
par(oma=c( 0,0,3,5))
image(1950:2000, e3043, 
  rotate270.matrix(anom_e_atl[e3043,]), 
  col=tim.colors(20), ann=FALSE, axes=FALSE,
  zlim=c(-100,100))
axis(2,at=c(560,580,600,620,640),
  labels=c("32.35", "34.11", "36.79", "38.27", "40.98"))
axis(1,at=c(1950,1960,1970,1980,1990,2000),
  labels=c("1950", "1960", "1970", "1980", "1990", "2000"))
title(ylab="Latitude", xlab="Year")
par(oma=c( 0,0,0,1))
image.plot(legend.only=TRUE, zlim=c(-100,100), legend.lab="mm/yr")
par <- op
title(main="Model Anomaly: East side of Atlantic 30-43N")
dev.off()

postscript(file="wAtlAnom30_43.ps")
par(oma=c( 0,0,3,5))
image(1950:2000, w3043, 
  rotate270.matrix(anom_w_atl[w3043,]), 
  col=tim.colors(20), ann=FALSE, axes=FALSE,
  zlim=c(-100,100))
axis(2,at=c(920,940,960,980,1000,1020,1040),
  labels=c("30.15", "32.55", "34.31", "36.41", "39.00", "40.80", "41.86"))
axis(1,at=c(1950, 1960,1970,1980,1990,2000),
  labels=c("1950","1960", "1970", "1980", "1990", "2000"))
title(ylab="Latitude", xlab="Year")
par(oma=c( 0,0,0,1))
image.plot(legend.only=TRUE, zlim=c(-100,100), legend.lab="mm/yr")
par <- op
title(main="Model Anomaly: West side of Atlantic 30-43N")
dev.off()
