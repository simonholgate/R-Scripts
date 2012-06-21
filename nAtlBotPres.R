# R script to read bottom pressure from Vassil's N. Atlantic model
# Data is Fortran sequential mode, is real*4, runs from 1950 to 2000 and 
# is monthly means. Missing value = 1e-30
library(fields)
# Domain is Mercator projection:
# xn=503 LINEAR  from -98.5 at spacing 0.234375 deg (max lats = 19.15625)
# yn=533 and non-linear (-35.08 to 65.37)
xn <- 503
yn <- 533
ndata <- xn*yn
nyr <- 51
years <- c(1950:2000)

lons <- seq(from=-98.5, by=0.234375, length=503)
lats <- c(  
     -35.08,-34.89,-34.69,-34.50,-34.31,-34.11,-33.92,-33.72,-33.53,
     -33.33,-33.14,-32.94,-32.74,-32.55,-32.35,-32.15,-31.95,-31.75,
     -31.55,-31.35,-31.15,-30.95,-30.75,-30.55,-30.35,-30.15,-29.94,
     -29.74,-29.54,-29.33,-29.13,-28.92,-28.72,-28.51,-28.30,-28.10,
     -27.89,-27.68,-27.48,-27.27,-27.06,-26.85,-26.64,-26.43,-26.22,
     -26.01,-25.80,-25.59,-25.38,-25.17,-24.95,-24.74,-24.53,-24.31,
     -24.10,-23.89,-23.67,-23.46,-23.24,-23.03,-22.81,-22.59,-22.38,
     -22.16,-21.94,-21.73,-21.51,-21.29,-21.07,-20.85,-20.63,-20.41,
     -20.19,-19.97,-19.75,-19.53,-19.31,-19.09,-18.87,-18.65,-18.42,
     -18.20,-17.98,-17.76,-17.53,-17.31,-17.08,-16.86,-16.64,-16.41,
     -16.19,-15.96,-15.74,-15.51,-15.28,-15.06,-14.83,-14.60,-14.38,
     -14.15,-13.92,-13.70,-13.47,-13.24,-13.01,-12.78,-12.55,-12.33,
     -12.10,-11.87,-11.64,-11.41,-11.18,-10.95,-10.72,-10.49,-10.26,
     -10.03, -9.80, -9.56, -9.33, -9.10, -8.87, -8.64, -8.41, -8.18,
      -7.94, -7.71, -7.48, -7.25, -7.01, -6.78, -6.55, -6.32, -6.08,
      -5.85, -5.62, -5.38, -5.15, -4.92, -4.68, -4.45, -4.21, -3.98,
      -3.75, -3.51, -3.28, -3.05, -2.81, -2.58, -2.34, -2.11, -1.87,
      -1.64, -1.41, -1.17, -0.94, -0.70, -0.47, -0.23,  0.00,  0.23,
       0.47,  0.70,  0.94,  1.17,  1.41,  1.64,  1.87,  2.11,  2.34,
       2.58,  2.81,  3.05,  3.28,  3.51,  3.75,  3.98,  4.21,  4.45,
       4.68,  4.92,  5.15,  5.38,  5.62,  5.85,  6.08,  6.32,  6.55,
       6.78,  7.01,  7.25,  7.48,  7.71,  7.94,  8.18,  8.41,  8.64,
       8.87,  9.10,  9.33,  9.56,  9.80, 10.03, 10.26, 10.49, 10.72,
      10.95, 11.18, 11.41, 11.64, 11.87, 12.10, 12.33, 12.55, 12.78,
      13.01, 13.24, 13.47, 13.70, 13.92, 14.15, 14.38, 14.60, 14.83,
      15.06, 15.28, 15.51, 15.74, 15.96, 16.19, 16.41, 16.64, 16.86,
      17.08, 17.31, 17.53, 17.76, 17.98, 18.20, 18.42, 18.65, 18.87,
      19.09, 19.31, 19.53, 19.75, 19.97, 20.19, 20.41, 20.63, 20.85,
      21.07, 21.29, 21.51, 21.73, 21.94, 22.16, 22.38, 22.59, 22.81,
      23.03, 23.24, 23.46, 23.67, 23.89, 24.10, 24.31, 24.53, 24.74,
      24.95, 25.17, 25.38, 25.59, 25.80, 26.01, 26.22, 26.43, 26.64,
      26.85, 27.06, 27.27, 27.48, 27.68, 27.89, 28.10, 28.30, 28.51,
      28.72, 28.92, 29.13, 29.33, 29.54, 29.74, 29.94, 30.15, 30.35,
      30.55, 30.75, 30.95, 31.15, 31.35, 31.55, 31.75, 31.95, 32.15,
      32.35, 32.55, 32.74, 32.94, 33.14, 33.33, 33.53, 33.72, 33.92,
      34.11, 34.31, 34.50, 34.69, 34.89, 35.08, 35.27, 35.46, 35.65,
      35.84, 36.03, 36.22, 36.41, 36.60, 36.79, 36.97, 37.16, 37.35,
      37.53, 37.72, 37.90, 38.09, 38.27, 38.46, 38.64, 38.82, 39.00,
      39.19, 39.37, 39.55, 39.73, 39.91, 40.09, 40.27, 40.45, 40.63,
      40.80, 40.98, 41.16, 41.33, 41.51, 41.68, 41.86, 42.03, 42.21,
      42.38, 42.55, 42.73, 42.90, 43.07, 43.24, 43.41, 43.58, 43.75,
      43.92, 44.09, 44.26, 44.42, 44.59, 44.76, 44.92, 45.09, 45.25,
      45.42, 45.58, 45.75, 45.91, 46.07, 46.24, 46.40, 46.56, 46.72,
      46.88, 47.04, 47.20, 47.36, 47.52, 47.68, 47.83, 47.99, 48.15,
      48.30, 48.46, 48.61, 48.77, 48.92, 49.08, 49.23, 49.38, 49.53,
      49.69, 49.84, 49.99, 50.14, 50.29, 50.44, 50.59, 50.74, 50.88,
      51.03, 51.18, 51.33, 51.47, 51.62, 51.76, 51.91, 52.05, 52.20,
      52.34, 52.48, 52.63, 52.77, 52.91, 53.05, 53.19, 53.33, 53.47,
      53.61, 53.75, 53.89, 54.02, 54.16, 54.30, 54.44, 54.57, 54.71,
      54.84, 54.98, 55.11, 55.25, 55.38, 55.51, 55.64, 55.78, 55.91,
      56.04, 56.17, 56.30, 56.43, 56.56, 56.69, 56.82, 56.94, 57.07,
      57.20, 57.33, 57.45, 57.58, 57.70, 57.83, 57.95, 58.08, 58.20,
      58.32, 58.45, 58.57, 58.69, 58.81, 58.93, 59.06, 59.18, 59.30,
      59.42, 59.53, 59.65, 59.77, 59.89, 60.01, 60.12, 60.24, 60.36,
      60.47, 60.59, 60.70, 60.82, 60.93, 61.04, 61.16, 61.27, 61.38,
      61.49, 61.61, 61.72, 61.83, 61.94, 62.05, 62.16, 62.27, 62.38,
      62.49, 62.59, 62.70, 62.81, 62.92, 63.02, 63.13, 63.23, 63.34,
      63.44, 63.55, 63.65, 63.76, 63.86, 63.96, 64.07, 64.17, 64.27,
      64.37, 64.47, 64.57, 64.67, 64.77, 64.87, 64.97, 65.07, 65.17,
       65.27, 65.37)

fileName <- "pb_50-00av.dat"
# Open binary file for reading
inConn <- file(fileName, "rb")

# Set up main sea surface height array
bp <- array(NA,dim=c(xn,yn,nyr))

for (i in 1:nyr){
# Read junk first integer that Fortran writes
  junk1 <- readBin(inConn, what="integer", n = 1, size = NA, endian = "little")
# Read year of data
  junk <- readBin(inConn, what="numeric", n = ndata, size = 4, endian = "little")
# Read junk last integer that Fortran writes
  junk1 <- readBin(inConn, what="integer", n = 1, size = NA, endian = "little")
  dim(junk) <- c(xn,yn)
  bp[,,i] <- junk
}

# End of data reading
close(inConn)

# Interpolate lats onto regular grid
library(fields)
#bpinterp <- array(NA,dim=c(xn,yn,nyr))
#latlin <- seq(from=-35.08, to=65.37, length=yn)
#grid.list <- make.surface.grid(list(lons, latlin))

#for (k in 1:nyr){
#  for (i in 1:xn){
#    junk <- approx(x=lats, y=bp[i,,k], xout=latlin)
#    bpinterp[i,,k] <- junk$y
#  } 
#}

# Find edge of continent for later
#contours <- contourLines(lons,lats,bp[,,1],levels=c(-1e29))

# Plot
bp[which(bp<=-1e29)]<-NA
image.plot(lons,lats,bp[,,1])
#lines(contours[[1]]$x,contours[[1]]$y, col='black', type='l')

## Export for GMT
#library(gmt)
#dataArray <- cbind(grid.list,as.vector(bpinterp[,,1]))
#bpinterp[which(bpinterp<=-1e2)]<-NaN
#r2gmt(dataArray,"data.gmt.xyz")

#x11()
#bpinterp[which(bpinterp<=-1e2)]<-NA
#image.plot(lons,latlin,bpinterp[,,1])

#nearest <- function (x,y) {
## function to return the index of the nearest value of vector x to value y
#  sorted <- sort(c(x,y))
## remove any duplicates in case y is identical with a member of x
#  usorted <- unique(sorted)
#  index <- which(usorted == y)
#  closest <- which.min(c(abs(y-usorted[index+1]), abs(y-usorted[index-1])))
#  if (closest == 1) {
#    nindex <- index+1
#  } else {
#    nindex <- index-1
#  }
#  xindex <- which(x == usorted[nindex])
#  return(xindex)
#}

## Identify East and West coast sea level
#lenContour <- length(contours[[1]]$x)
#coastalLons <- vector(length=lenContour, mode="integer")
#coastalLats <- vector(length=lenContour, mode="integer")
#coastalSsh <- array(NA,dim=c(lenContour,nyr))
#for (i in 1:length(contours[[1]]$x)){
#  coastalLons[i] <- nearest(lons,contours[[1]]$x[i])
#  coastalLats[i] <- nearest(lats,contours[[1]]$y[i])
## Extract sea levels at the coastal points
#  coastalSsh[i,] <- bp[coastalLons[i],coastalLons[i],]
#}
#image.plot(c(1:lenContour),years,coastalSsh)

# Masking data
#source("~/bin/RScripts/writeXPM.R")
#writeXPM("bp",bp[,,1])

source("~/bin/RScripts/readXPM.R")
source("~/bin/RScripts/matrixMethods.R")

mask_ew <- readXPM("mask_ew_atl")
mask_ew[which(mask_ew==1)] <- NA
mask_ew[which(mask_ew==2)] <- 1
mask_ew <- rotate270.matrix(mask_ew)

temp <- mask_ew
temp[291:503,] <- NA
mask_w <- which(is.finite(temp), arr.ind=TRUE)
temp <- mask_ew
temp[1:290,] <- NA
mask_e <- which(is.finite(temp), arr.ind=TRUE)
rm(temp)

masked_w_atl <- array(NA, dim=c(length(mask_w[,1]),nyr))
masked_e_atl <- array(NA, dim=c(length(mask_e[,1]),nyr))

for (i in 1:nyr){
  for (j in 1:length(mask_e[,1])){
    masked_e_atl[j,i] <- bp[mask_e[j,1],mask_e[j,2],i]
  }
  for (j in 1:length(mask_w[,1])){
    masked_w_atl[j,i] <- bp[mask_w[j,1],mask_w[j,2],i]
  }
}

lats_e <- lats[mask_e[,2]]
lats_w <- lats[mask_w[,2]]

postscript(file="eAtl.ps")
image.plot(1:length(mask_e[,1]), 1950:2000, masked_e_atl, xlab="station",
  ylab="Year", main="East side of Atlantic")
dev.off()
postscript(file="wAtl.ps")
image.plot(1:length(mask_w[,1]), 1950:2000, masked_w_atl, xlab="station",
  ylab="Year", main="West side of Atlantic")
dev.off()
postscript(file="atlMask.ps")
image.plot(lons,lats,mask_ew*bp[,,1], xlab="Lon", ylab="Lat", 
  main="E-W Mask")
dev.off()

# Decadal rates
X <- 1:10
masked_w_atl_rates <- array(NA, dim=c(length(mask_w[,1]), 41))
masked_e_atl_rates <- array(NA, dim=c(length(mask_e[,1]), 41))
for (i in 1:length(mask_w[,1])){
  for (j in 1:41) {
    fit <- lm(masked_w_atl[i,j:(j+9)] ~ X, x=TRUE)
    masked_w_atl_rates[i,j] <- fit$coeff[2]
  }
}
for (i in 1:length(mask_e[,1])){
  for (j in 1:41) {
    fit <- lm(masked_e_atl[i,j:(j+9)] ~ X, x=TRUE)
    masked_e_atl_rates[i,j] <- fit$coeff[2]
  }
}

postscript(file="eAtlRates.ps")
image.plot(1:length(mask_e[,1]), 1955:1995, masked_e_atl_rates, xlab="station",
  ylab="Year", main="East side of Atlantic")
dev.off()
postscript(file="wAtlRates.ps")
image.plot(1:length(mask_w[,1]), 1955:1995, masked_w_atl_rates, xlab="station",
  ylab="Year", main="West side of Atlantic")
dev.off()

op <- par()
par(oma=c( 0,0,3,5))
image(1:length(mask_e[,1]), 1955:1995, masked_e_atl_rates, xlab="Station",
  ylab="Year", col=tim.colors(10),
  zlim=c(-0.04,0.04))
axis(3,at=c(200,400,600,800),
  labels=c("50.59", "35.27", "8.18", "-5.62"))
#  labels=c("-3.98","9.8","36.79","51.33"))
title(ylab="Latitude")
par(oma=c( 0,0,0,1))
image.plot(legend.only=TRUE, zlim=c(-0.04,0.04), legend.lab="mm/yr")
par <- op
title(main="East side of Atlantic")

par(oma=c( 0,0,3,5))
image(1:length(mask_w[,1]), 1955:1995, masked_w_atl_rates, xlab="Station",
  ylab="Year", col=tim.colors(10),
  zlim=c(-0.04,0.04))
axis(3,at=c(200,400,600,800,1000,1200,1400),
  labels=c("51.03","44.42","29.54","20.63","10.72","-0.47","-23.46"))
#  labels=c("-9.80", "8.64", "15.74", "27.06", "39.00", "48.30", "59.18"))
title(ylab="Latitude")
par(oma=c( 0,0,0,1))
image.plot(legend.only=TRUE, zlim=c(-0.04,0.04), legend.lab="mm/yr")
par <- op
title(main="West side of Atlantic")
