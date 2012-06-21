###############################################################################
# brestNewlynWindStressACRE.R
# 
# Reads the fields of annual mean wind stress (calculated from the daily winds 
# from the 20th Century Reanalysis Project) and extract the time series at Brest
# and Newlyn of the N and E components. 
#
# Author: simonh
###############################################################################

# Newlyn 50 06 N  05 33 W => 50.1 -5.55
# Brest 48 23 N  04 30 W

library(fields)
source("~/Dropbox/BrestNewlyn/matrixMethods.R")

l <- 192
m <- 94
lm <- l*m
numYears <- 138

lat <- c(88.542, 86.6531, 84.7532, 82.8508, 80.9473, 79.0435, 77.1394, 75.2351, 
    73.3307, 71.4262, 69.5217, 67.6171, 65.7125, 63.8079, 61.9033, 59.9986, 
    58.0939, 56.1893, 54.2846, 52.3799, 50.4752, 48.5705, 46.6658, 44.7611, 
    42.8564, 40.9517, 39.047, 37.1422, 35.2375, 33.3328, 31.4281, 29.5234, 
    27.6186, 25.7139, 23.8092, 21.9044, 19.9997, 18.095, 16.1902, 14.2855, 
    12.3808, 10.47604, 8.57131, 6.66657, 4.76184, 2.8571, 0.952368, 
    -0.952368, -2.8571, -4.76184, -6.66657, -8.57131, -10.47604, -12.3808, 
    -14.2855, -16.1902, -18.095, -19.9997, -21.9044, -23.8092, -25.7139, 
    -27.6186, -29.5234, -31.4281, -33.3328, -35.2375, -37.1422, -39.047, 
    -40.9517, -42.8564, -44.7611, -46.6658, -48.5705, -50.4752, -52.3799, 
    -54.2846, -56.1893, -58.0939, -59.9986, -61.9033, -63.8079, -65.7125, 
    -67.6171, -69.5217, -71.4262, -73.3307, -75.2351, -77.1394, -79.0435, 
    -80.9473, -82.8508, -84.7532, -86.6531, -88.542)

lon <- c(0, 1.875, 3.75, 5.625, 7.5, 9.375, 11.25, 13.125, 15, 16.875, 18.75, 
    20.625, 22.5, 24.375, 26.25, 28.125, 30, 31.875, 33.75, 35.625, 37.5, 
    39.375, 41.25, 43.125, 45, 46.875, 48.75, 50.625, 52.5, 54.375, 56.25, 
    58.125, 60, 61.875, 63.75, 65.625, 67.5, 69.375, 71.25, 73.125, 75, 
    76.875, 78.75, 80.625, 82.5, 84.375, 86.25, 88.125, 90, 91.875, 93.75, 
    95.625, 97.5, 99.375, 101.25, 103.125, 105, 106.875, 108.75, 110.625, 
    112.5, 114.375, 116.25, 118.125, 120, 121.875, 123.75, 125.625, 127.5, 
    129.375, 131.25, 133.125, 135, 136.875, 138.75, 140.625, 142.5, 144.375, 
    146.25, 148.125, 150, 151.875, 153.75, 155.625, 157.5, 159.375, 161.25, 
    163.125, 165, 166.875, 168.75, 170.625, 172.5, 174.375, 176.25, 178.125, 
    180, 181.875, 183.75, 185.625, 187.5, 189.375, 191.25, 193.125, 195, 
    196.875, 198.75, 200.625, 202.5, 204.375, 206.25, 208.125, 210, 211.875, 
    213.75, 215.625, 217.5, 219.375, 221.25, 223.125, 225, 226.875, 228.75, 
    230.625, 232.5, 234.375, 236.25, 238.125, 240, 241.875, 243.75, 245.625, 
    247.5, 249.375, 251.25, 253.125, 255, 256.875, 258.75, 260.625, 262.5, 
    264.375, 266.25, 268.125, 270, 271.875, 273.75, 275.625, 277.5, 279.375, 
    281.25, 283.125, 285, 286.875, 288.75, 290.625, 292.5, 294.375, 296.25, 
    298.125, 300, 301.875, 303.75, 305.625, 307.5, 309.375, 311.25, 313.125, 
    315, 316.875, 318.75, 320.625, 322.5, 324.375, 326.25, 328.125, 330, 
    331.875, 333.75, 335.625, 337.5, 339.375, 341.25, 343.125, 345, 346.875, 
    348.75, 350.625, 352.5, 354.375, 356.25, 358.125)

inStressE1871 <- file("tau_x_1871.bin", "rb")
inStressN1871 <- file("tau_y_1871.bin", "rb")
inStressE1921 <- file("tau_x_1921.bin", "rb")
inStressN1921 <- file("tau_y_1921.bin", "rb")

stressE1871 <- array(NA, dim=c(l, m))
stressN1871 <- array(NA, dim=c(l, m))
stressE1921 <- array(NA, dim=c(l, m))
stressN1921 <- array(NA, dim=c(l, m))


  junk <- readBin(inStressE1871, "integer", size=4, n = 1)
  junkw <- readBin(inStressE1871, "numeric", size=8, n = lm)
  dim(junkw) <- c(l,m)
  stressE1871 <- mirror.matrix(junkw)
  junk <- readBin(inStressE1871, "integer", size=4, n = 1)

  junk <- readBin(inStressN1871, "integer", size=4, n = 1)
  junkw <- readBin(inStressN1871, "numeric", size=8, n = lm)
  dim(junkw) <- c(l,m)
  stressN1871 <- mirror.matrix(junkw)
  junk <- readBin(inStressN1871, "integer", size=4, n = 1)

  junk <- readBin(inStressE1921, "integer", size=4, n = 1)
  junkw <- readBin(inStressE1921, "numeric", size=8, n = lm)
  dim(junkw) <- c(l,m)
  stressE1921 <- mirror.matrix(junkw)
  junk <- readBin(inStressE1921, "integer", size=4, n = 1)

  junk <- readBin(inStressN1921, "integer", size=4, n = 1)
  junkw <- readBin(inStressN1921, "numeric", size=8, n = lm)
  dim(junkw) <- c(l,m)
  stressN1921 <- mirror.matrix(junkw)
  junk <- readBin(inStressN1921, "integer", size=4, n = 1)

  
  close(inStressE1871)
  close(inStressN1871)
  close(inStressE1921)
  close(inStressN1921)
    
#x11()
#negs <- which(stressE1871[,,1] <= 0)
#negE <- mirror.matrix(stressE1871[,,1])
#poss <- which(stressE1871[,,1] > 0)
#posE <- mirror.matrix(stressE1871[,,1])
#posE[negs] <- NA
#negE[poss] <- NA
#contour(lon,(-1*lat),posE, col='red')
#contour(lon,(-1*lat),negE, col='blue', add=T)
#world(shift=T, add=T)

#x11()
#negs <- which(stressN1871[,,1] <= 0)
#negN <- mirror.matrix(stressN1871[,,1])
#poss <- which(stressN1871[,,1] > 0)
#posN <- mirror.matrix(stressN1871[,,1])
#posN[negs] <- NA
#negN[poss] <- NA
#contour(lon,(-1*lat),posN, col='red')
#contour(lon,(-1*lat),negN, col='blue', add=T)
#world(shift=T, add=T)

x11()
image.plot(lon,(-1*lat),stressE1871, zlim=c(-50, 50))
world(shift=T, add=T)

x11()
image.plot(lon,(-1*lat),stressN1871, zlim=c(-50, 50))
world(shift=T, add=T)

x11()
image.plot(lon,(-1*lat),stressE1921, zlim=c(-50, 50))
world(shift=T, add=T)

x11()
image.plot(lon,(-1*lat),stressN1921, zlim=c(-50, 50))
world(shift=T, add=T)

## x11()
## wslon <- seq(from=-180, to = 178.125, by = 1.875)
## wsNAnnualMean <- stressN1871
## wsNAnnualMean[1:96,,] <- stressN1871[97:192,,]
## wsNAnnualMean[97:192,,] <- stressN1871[1:96,,]
## filled.contour(wslon,(-1*lat), wsNAnnualMean[,,1], color=tim.colors, nlevels=20,
## plot.axes = { world(add=T) }, zlim=c(-0.12, 0.15))

## x11()
## wsEAnnualMean <- stressE1871
## wsEAnnualMean[1:96,,] <- stressE1871[97:192,,]
## wsEAnnualMean[97:192,,] <- stressE1871[1:96,,]
## filled.contour(wslon,(-1*lat),wsEAnnualMean[,,1],
##                color=tim.colors, nlevels=20,
## plot.axes = { world(add=T) }, zlim=c(-0.12, 0.15))

## x11()
## lonLat <- expand.grid(wslon[seq(from=2, to=192, by=4)],
##                       (-1*lat[seq(from=2, to=94, by=4)]))
## plot(lonLat[,1], lonLat[,2], type='p', pch=".", xlab="Longitude",
##      ylab="Latitude")
## arrow.plot(lonLat[,1], lonLat[,2],
##            u=as.vector(wsEAnnualMean[seq(from=2, to=192, by=4),
##              seq(from=2, to=94, by=4),138]),
##            v=as.vector(wsNAnnualMean[seq(from=2, to=192, by=4),
##              seq(from=2, to=94, by=4),138]),
##            arrow.ex=.2, col='magenta', length=.05,
##            true.angle=T, plot.axes = { world(add=T) })
## arrow.plot(0,100, u=0.1, v=0, arrow.ex=.2, col='magenta',
##            length=.05, true.angle=T)
## text(10, 100, "0.1 m/s", pos=3)
## points(-20,27, pch=19, col="magenta")


##save(file="~/data/ACRE/wsACRE.RData", wsEAnnualMean, stressE1871, wsNAnnualMean, stressN1871, wslon, lon, lat)
