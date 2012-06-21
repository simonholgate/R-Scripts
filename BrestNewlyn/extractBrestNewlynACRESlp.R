# Find the nearest ACRE grid box for a number of sites and interpolate the monthly SLP to a point 
library(fields)
source("/home/simonh/workspace/RScripts/matrixMethods.R")

ntotstns <- (2+35)
nlon<-180
nlat<-91
nyr<-138
nmon<-12
xlon<-seq(from=0,by=2,length=nlon)
ylat<-seq(from=90,by=-2,length=nlat)
ylatr <- seq(from=-90,by=2,length=nlat)

# Newlyn 50 06 N  05 33 W => 50.1 -5.55
# Brest 48 23 N  04 30 W => 48.38 -4.5
# Add extra stations from around the model domain here similar to Thompson (1980)
# The 35 stations are on a regular 5 deg. grid covering 40N-60N, 20W-10E
# as I can do that with this dataset. 

location<-new.env()

location$lat <- c(48.38, 50.1,
                  40, 40, 40, 40, 40, 40, 40,
                  45, 45, 45, 45, 45, 45, 45,
                  50, 50, 50, 50, 50, 50, 50, 
                  55, 55, 55, 55, 55, 55, 55, 
                  60, 60, 60, 60, 60, 60, 60)
location$lon <- c(-4.5, -5.55,
                  -20, -15, -10, -5, 0, 5, 10,
                  -20, -15, -10, -5, 0, 5, 10,
                  -20, -15, -10, -5, 0, 5, 10,
                  -20, -15, -10, -5, 0, 5, 10,
                  -20, -15, -10, -5, 0, 5, 10)

# interp.surface needs values within the same range as the grid i.e. 0-360
location$lon360 <- c(360-4.5, 360-5.55,
                      360-20, 360-15, 360-10, 360-5, 0, 5, 10,
                     360-20, 360-15, 360-10, 360-5, 0, 5, 10,
                     360-20, 360-15, 360-10, 360-5, 0, 5, 10,
                     360-20, 360-15, 360-10, 360-5, 0, 5, 10,
                     360-20, 360-15, 360-10, 360-5, 0, 5, 10)

location$lonlat <- cbind(location$lon360, location$lat)

library(ncdf)
source("/home/simonh/workspace/RScripts/matrixMethods.R")

nc <- open.ncdf("~/data/ACRE/prmsl.mon.mean.nc")
slp.id <- nc$var[[2]]
slp <- get.var.ncdf(nc, slp.id)
lengthT <- nc$dim[[4]]$len
t <- seq.Date(from=as.Date("1871-01-15"),to=as.Date("2008-12-15"),length.out=lengthT)
close.ncdf(nc)

# For each of the lengthT months, we need to interpolate the grid to the locations we've chosen
slpNewlynStns <- array(NA,dim=c(lengthT,ntotstns))
for (i in 1:lengthT) {
  obj <- list(x=xlon, y=ylat, z=mirror.matrix(slp[,,i]))
  slpNewlynStns[i,] <- interp.surface(obj,location$lonlat)
}

# Save the data
save(slpNewlynStns, file="brestNewlynACREslp.RData")
