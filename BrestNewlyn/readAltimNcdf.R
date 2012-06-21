# Script to read Laurent Testut's along track altimetry data
library(ncdf)
nc <- open.ncdf("~/diskx/altimetry/netcdf_brest_newlin/track-raw.TPJ1.070.ref.sla.nc")

# Get variable information

#     print(paste("File",nc$filename,"contains",nc$nvars,"variables"))
#     for( i in 1:nc$nvars ) {
#             v <- nc$var[[i]]
#             print(paste("Here is information on variable number",i))
#             print(paste("   Name: ",v$name))
#             print(paste("   Units:",v$units))
#             print(paste("   Missing value:",v$missval))
#             print(paste("   # dimensions :",v$ndims))
#             print(paste("   Variable size:",v$varsize))
#             }

#[1] "File ~/diskx/altimetry/netcdf_brest_newlin/track-raw.TPJ1.070.ref.sla.nc contains 8 variables"

#[1] "Here is information on variable number 1"
#[1] "   Name:  lon"
#[1] "   Units: deg"
#[1] "   Missing value: 1e+30"
#[1] "   # dimensions : 1"
#[1] "   Variable size: 313"
#[1] "Here is information on variable number 2"
#[1] "   Name:  lat"
#[1] "   Units: deg"
#[1] "   Missing value: 1e+30"
#[1] "   # dimensions : 1"
#[1] "   Variable size: 313"
#[1] "Here is information on variable number 3"
#[1] "   Name:  mssh"
#[1] "   Units: m"
#[1] "   Missing value: 1e+30"
#[1] "   # dimensions : 1"
#[1] "   Variable size: 313"
#[1] "Here is information on variable number 4"
#[1] "   Name:  cycle"
#[1] "   Units: "
#[1] "   Missing value: NA"
#[1] "   # dimensions : 1"
#[1] "   Variable size: 526"
#[1] "Here is information on variable number 5"
#[1] "   Name:  time"
#[1] "   Units: Julian days since 1950-01-01"
#[1] "   Missing value: 1e+30"
#[1] "   # dimensions : 2"
#[1] "   Variable size: 526" "   Variable size: 313"
#[1] "Here is information on variable number 6"
#[1] "   Name:  sla"
#[1] "   Units: m"
#[1] "   Missing value: 1e+30"
#[1] "   # dimensions : 2"
#[1] "   Variable size: 526" "   Variable size: 313"
#[1] "Here is information on variable number 7"
#[1] "   Name:  tide"
#[1] "   Units: m"
#[1] "   Missing value: 1e+30"
#[1] "   # dimensions : 2"
#[1] "   Variable size: 526" "   Variable size: 313"
#[1] "Here is information on variable number 8"
#[1] "   Name:  mog2d"
#[1] "   Units: m"
#[1] "   Missing value: 1e+30"
#[1] "   # dimensions : 2"
#[1] "   Variable size: 526" "   Variable size: 313"

lon70 <- get.var.ncdf(nc, nc$var[[1]])
lat70 <- get.var.ncdf(nc, nc$var[[2]])
mssh70 <- get.var.ncdf(nc, nc$var[[3]])
time70 <- get.var.ncdf(nc, nc$var[[5]])
sla70 <- get.var.ncdf(nc, nc$var[[6]])

sla_na70 <- sla70[1,]

close.ncdf(nc)

nc <- open.ncdf("~/diskx/altimetry/netcdf_brest_newlin/track-raw.TPJ1.239.ref.sla.nc")
lon239 <- get.var.ncdf(nc, nc$var[[1]])
lat239 <- get.var.ncdf(nc, nc$var[[2]])
mssh239 <- get.var.ncdf(nc, nc$var[[3]])
time239 <- get.var.ncdf(nc, nc$var[[5]])
sla239 <- get.var.ncdf(nc, nc$var[[6]])
sla_na239 <- sla239[1,]
close.ncdf(nc)

nc <- open.ncdf("~/diskx/altimetry/netcdf_brest_newlin/track-raw.TPJ1.248.ref.sla.nc")
lon248 <- get.var.ncdf(nc, nc$var[[1]])
lat248 <- get.var.ncdf(nc, nc$var[[2]])
mssh248 <- get.var.ncdf(nc, nc$var[[3]])
time248 <- get.var.ncdf(nc, nc$var[[5]])
sla248 <- get.var.ncdf(nc, nc$var[[6]])
sla_na248 <- sla248[1,]
close.ncdf(nc)

nc <- open.ncdf("~/diskx/altimetry/netcdf_brest_newlin/track-raw.TPJ1.061.ref.sla.nc")
lon61 <- get.var.ncdf(nc, nc$var[[1]])
lat61 <- get.var.ncdf(nc, nc$var[[2]])
mssh61 <- get.var.ncdf(nc, nc$var[[3]])
time61 <- get.var.ncdf(nc, nc$var[[5]])
sla61 <- get.var.ncdf(nc, nc$var[[6]])
sla_na61 <- sla61[1,]
close.ncdf(nc)

sla_na70[which(sla_na70>99.9)] <- NA
nn70 <- which(mssh70>99.9)
lon_na70<-lon70
lon_na70[nn70] <- NA
lat_na70<-lat70
lat_na70[nn70] <- NA

sla_na239[which(sla_na239>99.9)] <- NA
nn239 <- which(mssh239>99.9)
lon_na239<-lon239
lon_na239[nn239] <- NA
lat_na239<-lat239
lat_na239[nn239] <- NA

sla_na248[which(sla_na248>99.9)] <- NA
nn248 <- which(mssh248>99.9)
lon_na248<-lon248
lon_na248[nn248] <- NA
lat_na248<-lat248
lat_na248[nn248] <- NA

sla_na61[which(sla_na61>99.9)] <- NA
nn61 <- which(mssh61>99.9)
lon_na61<-lon61
lon_na61[nn61] <- NA
lat_na61<-lat61
lat_na61[nn61] <- NA

mssh_na61<-mssh61
mssh_na61[nn61] <- NA
mssh_na70<-mssh70
mssh_na70[nn70] <- NA
mssh_na239<-mssh239
mssh_na239[nn239] <- NA
mssh_na248<-mssh248
mssh_na248[nn248] <- NA

# Write out file for use in GMT
fins <- which(is.finite(sla_na70))
sla_vec70 <- sla_na70[fins]
lon_vec70 <- lon_na70[fins]
lat_vec70 <- lat_na70[fins]

fins <- which(is.finite(sla_na239))
sla_vec239 <- sla_na239[fins]
lon_vec239 <- lon_na239[fins]
lat_vec239 <- lat_na239[fins]

fins <- which(is.finite(sla_na248))
sla_vec248 <- sla_na248[fins]
lon_vec248 <- lon_na248[fins]
lat_vec248 <- lat_na248[fins]

fins <- which(is.finite(sla_na61))
sla_vec61 <- sla_na61[fins]
lon_vec61 <- lon_na61[fins]
lat_vec61 <- lat_na61[fins]

# Convert sla to mm
sla_vec70 <-sla_vec70*1000
sla_vec239 <-sla_vec239*1000
sla_vec248 <-sla_vec248*1000
sla_vec61 <-sla_vec61*1000

lon_vec <- c(lon_vec70, lon_vec239, lon_vec248, lon_vec61)
lat_vec <- c(lat_vec70, lat_vec239, lat_vec248, lat_vec61)
sla_vec <- c(sla_vec70, sla_vec239, sla_vec248, sla_vec61)

xyz <- cbind(lon_vec,lat_vec,sla_vec)
write.table(xyz, "sla.xyz", col.names=FALSE, row.names=FALSE)

# For fields
library(fields)

# Need to work out how to convert sla_na into colours to plot at points
# along the line
split.screen( rbind(c(0, .8,0,1), c(.8,1,0,1)))
par(family="HersheySans")

screen(1)
plot(lon_vec,lat_vec, type='n', ann=F) 
title(xlab='Lon', ylab='Lat', main='Tracks 70, 61, 239 & 248')
ribbon.plot(lon_vec70, lat_vec70, sla_vec70, lwd=10, zlim=c(-200,200))
ribbon.plot(lon_vec239, lat_vec239, sla_vec239, lwd=10, zlim=c(-200,200))
ribbon.plot(lon_vec248, lat_vec248, sla_vec248, lwd=10, zlim=c(-200,200))
ribbon.plot(lon_vec61, lat_vec61, sla_vec61, lwd=10, zlim=c(-200,200))
world(add=T)

screen(2)
image.plot(zlim=c(-200,200),legend.only=TRUE, smallplot=c(.1,.2, .3,.7),
  col=tim.colors(), legend.lab='SLA [mm]')

close.screen( all=TRUE)

x11()
split.screen( rbind(c(0, .8,0,1), c(.8,1,0,1)))
par(family="HersheySans")

screen(1)
plot(lon_vec,lat_vec, type='n', ann=F) 
title(xlab='Lon', ylab='Lat', main='Tracks 70, 61, 239 & 248')
ribbon.plot(lon_na70, lat_na70, mssh_na70, lwd=10, zlim=c(48,62))
ribbon.plot(lon_na239, lat_na239, mssh_na239, lwd=10, zlim=c(48,62))
ribbon.plot(lon_na248, lat_na248, mssh_na248, lwd=10, zlim=c(48,62))
ribbon.plot(lon_na61, lat_na61, mssh_na61, lwd=10, zlim=c(48,62))
world(add=T)

screen(2)
image.plot(zlim=c(48,62),legend.only=TRUE, smallplot=c(.1,.2, .3,.7),
  col=tim.colors(), legend.lab='MSSH [mm]')

close.screen( all=TRUE)


# Calculate trends at each point along track

trends70 <- vector(mode="numeric", length=length(mssh70))
X <- 1:dim(sla70)[1]
for (i in 1:length(mssh70)){
  sla_nn <- which(sla70[,i] <= 99)
  if (length(sla_nn)>0){
    junk <- lm(sla70[sla_nn,i] ~ X[sla_nn])
    trends70[i] <- junk$coeff[2]*1000*52
  }
}
trends70[nn70] <- NA

trends61 <- vector(mode="numeric", length=length(mssh61))
X <- 1:dim(sla61)[1]
for (i in 1:length(mssh61)){
  sla_nn <- which(sla61[,i] <= 99)
  if (length(sla_nn)>0){
    junk <- lm(sla61[sla_nn,i] ~ X[sla_nn])
    trends61[i] <- junk$coeff[2]*1000*52
  }
}
trends61[nn61] <- NA

trends239 <- vector(mode="numeric", length=length(mssh239))
X <- 1:dim(sla239)[1]
for (i in 1:length(mssh239)){
  sla_nn <- which(sla239[,i] <= 99)
  if (length(sla_nn)>0){
    junk <- lm(sla239[sla_nn,i] ~ X[sla_nn])
    trends239[i] <- junk$coeff[2]*1000*52
  }
}
trends239[nn239] <- NA

trends248 <- vector(mode="numeric", length=length(mssh248))
X <- 1:dim(sla248)[1]
for (i in 1:length(mssh248)){
  sla_nn <- which(sla248[,i] <= 99)
  if (length(sla_nn)>0){
    junk <- lm(sla248[sla_nn,i] ~ X[sla_nn])
    trends248[i] <- junk$coeff[2]*1000*52
  }
}
trends248[nn248] <- NA

x11()
split.screen( rbind(c(0, .8,0,1), c(.8,1,0,1)))
par(family="HersheySans")

screen(1)
plot(lon_vec,lat_vec, type='n', ann=F) 
title(xlab='Lon', ylab='Lat', main='Tracks 70, 61, 239 & 248')
ribbon.plot(lon_na70, lat_na70, trends70, lwd=10, zlim=c(-15,15))
ribbon.plot(lon_na239, lat_na239, trends239, lwd=10, zlim=c(-15,15))
ribbon.plot(lon_na248, lat_na248, trends248, lwd=10, zlim=c(-15,15))
ribbon.plot(lon_na61, lat_na61, trends61, lwd=10, zlim=c(-15,15))
world(add=T)

screen(2)
image.plot(zlim=c(-15,15),legend.only=TRUE, smallplot=c(.1,.2, .3,.7),
  col=tim.colors(), legend.lab='Trends [mm/yr]')

close.screen( all=TRUE)

# Find cross-over points.
#library(spam)
source("~/workspace/RScripts/greatcircle.R")

x70 <- cbind(lon_vec70, lat_vec70)
y61 <- cbind(lon_vec61, lat_vec61)

y239 <- cbind(lon_vec239, lat_vec239)
x248 <- cbind(lon_vec248, lat_vec248)

near7061 <- array(NA, dim=c(length(lon_vec70),length(lon_vec61)))
for (i in 1:length(lon_vec70)){
  for (j in 1:length(lon_vec61)){
    near7061[i,j] <- greatcircle(x70[i,2], x70[i,1], y61[j,2], y61[j,1])
  }
}
nearest7061 <- which(near7061==min(near7061),arr.ind=T)
#nearest7061 <- array(NA, dim=c(2,2))
#nearest7061[1,] <- x70[nearest[1],]
#nearest7061[2,] <- y61[nearest[2],]

near70239 <- array(NA, dim=c(length(lon_vec70),length(lon_vec239)))
for (i in 1:length(lon_vec70)){
  for (j in 1:length(lon_vec239)){
    near70239[i,j] <- greatcircle(x70[i,2], x70[i,1], y239[j,2], y239[j,1])
  }
}
nearest70239 <- which(near70239==min(near70239),arr.ind=T)
#nearest70239 <- array(NA, dim=c(2,2))
#nearest70239[1,] <- x70[nearest[1],]
#nearest70239[2,] <- y239[nearest[2],]

near24861 <- array(NA, dim=c(length(lon_vec248),length(lon_vec61)))
for (i in 1:length(lon_vec248)){
  for (j in 1:length(lon_vec61)){
    near24861[i,j] <- greatcircle(x248[i,2], x248[i,1], y61[j,2], y61[j,1])
  }
}
nearest24861 <- which(near24861==min(near24861),arr.ind=T)
#nearest24861 <- array(NA, dim=c(2,2))
#nearest24861[1,] <- x248[nearest[1],]
#nearest24861[2,] <- y61[nearest[2],]

near248239 <- array(NA, dim=c(length(lon_vec248),length(lon_vec239)))
for (i in 1:length(lon_vec248)){
  for (j in 1:length(lon_vec239)){
    near248239[i,j] <- greatcircle(x248[i,2], x248[i,1], y239[j,2], y239[j,1])
  }
}
nearest248239 <- which(near248239==min(near248239),arr.ind=T)


# x-points
sla_x7061 <- sla70[,nearest7061[1]]
fins <- which(sla_x7061 > 99.9)
sla_x7061[fins] <- NA
time7061 <- time70[,nearest7061[1]]
fins <- which(time7061<100)
time7061[fins] <- NA

sla_x6170 <- sla61[,nearest7061[2]]
fins <- which(sla_x6170 > 99.9)
sla_x6170[fins] <- NA
time6170 <- time61[,nearest7061[2]]
fins <- which(time6170<100)
time6170[fins] <- NA

sla_x70239 <- sla70[,nearest70239[1]]
fins <- which(sla_x70239 > 99.9)
sla_x70239[fins] <- NA
time70239 <- time70[,nearest70239[1]]
fins <- which(time70239<100)
time70239[fins] <- NA

sla_x23970 <- sla239[,nearest70239[2]]
fins <- which(sla_x23970 > 99.9)
sla_x23970[fins] <- NA
time23970 <- time239[,nearest70239[2]]
fins <- which(time23970<100)
time23970[fins] <- NA

sla_x24861 <- sla248[,nearest24861[1]]
fins <- which(sla_x24861 > 99.9)
sla_x24861[fins] <- NA
time24861 <- time248[,nearest24861[1]]
fins <- which(time24861<100)
time24861[fins] <- NA

sla_x61248 <- sla61[,nearest24861[2]]
fins <- which(sla_x61248 > 99.9)
sla_x61248[fins] <- NA
time61248 <- time61[,nearest24861[2]]
fins <- which(time61248<100)
time61248[fins] <- NA

sla_x248239 <- sla248[,nearest248239[1]]
fins <- which(sla_x248239 > 99.9)
sla_x248239[fins] <- NA
time248239 <- time239[,nearest248239[1]]
fins <- which(time248239<100)
time248239[fins] <- NA

sla_x239248 <- sla239[,nearest248239[2]]
fins <- which(sla_x239248 > 99.9)
sla_x239248[fins] <- NA
time239248 <- time239[,nearest248239[2]]
fins <- which(time239248<100)
time239248[fins] <- NA

xpoints <- cbind(time7061, sla_x7061, time6170, sla_x6170, 
  time70239, sla_x70239, time23970, sla_x23970, time24861, sla_x24861, 
  time61248, sla_x61248, time248239, sla_x248239, time239248, sla_x239248)
  
write.table(xpoints, file="xpoints.txt", na="NaN", col.names=F, row.names=F)