# Find the nearest ACRE grid box for Key West and San Francisco and interpolate the monthly SLP to a point 
library(fields)
source("/home/simonh/workspace/RScripts/matrixMethods.R")

ntotstns <- 2
nlon<-180
nlat<-91
nyr<-138
nmon<-12
xlon<-seq(from=0,by=2,length=nlon)
ylat<-seq(from=90,by=-2,length=nlat)
ylatr <- seq(from=-90,by=2,length=nlat)
sanfran.lon <- -122.467
sanfran.lat <- 37.8
keywest.lon <- -81.8
keywest.lat <- 24.55

location<-new.env()

location$lat <- c(keywest.lat,sanfran.lat)

location$lon <- c(keywest.lon,sanfran.lon)


# interp.surface needs values within the same range as the grid i.e. 0-360
location$lon360 <- c(360+keywest.lon, 360+sanfran.lon)

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
slpKeyWestStns <- array(NA,dim=c(lengthT,ntotstns))
for (i in 1:lengthT) {
  obj <- list(x=xlon, y=ylat, z=mirror.matrix(slp[,,i]))
  slpKeyWestStns[i,] <- interp.surface(obj,location$lonlat)
}

# Save the data
save(slpKeyWestStns, file="keyWestACREslp.RData")
