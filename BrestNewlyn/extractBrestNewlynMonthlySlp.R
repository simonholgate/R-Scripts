#
library(fields)

nstns<-8
nextrastns<-17
ntotstns <- nstns+nextrastns
nlon<-72
nlat<-37
nyr<-158
nmon<-12
xlon<-seq(from=-180,by=5,length=nlon)
xlon2 <- seq(from=180,by=5,length=nlon)
ylat<-seq(from=90,by=-5,length=nlat)
x<-seq(from=1,length=nlon)
y<-seq(from=1,length=nlat)

# Add extra stations from around the model domain here similar to Thompson (1980)
# St Mawgan - Latitude: 50.454243N  Longitude: 4.99915W
#> (-4.99915-xmin)/xres -> 89
#> (50.454243-ymin)/yres -> 93
# Thorney Island - latitude 50.8166667 longitude -0.9166667
#> (-0.9166667-xmin)/xres -> 114
#> (50.8166667-ymin)/yres -> 96
# Shoeburyness - Latitude: 51.55; Longitude: 0.833
#> (0.833-xmin)/xres -> 124
#> (51.55-ymin)/yres -> 103
# Gorleston - Latitude, 52.5833, Longitude, 1.7167
#> (1.7167-xmin)/xres -> 129
#> (52.5833-ymin)/yres -> 112
# Kilnsea - Latitude, 53.6167, Longitude, 0.1333
#> (0.1333-xmin)/xres -> 120
#> (53.6167-ymin)/yres -> 122
# Eskdalemuir - Latitude: 55.317; Longitude: -3.2
#> (-3.2-xmin)/xres -> 100
#> (55.317-ymin)/yres -> 137
# Kinloss - Latitude: 57.65; Longitude: -3.567
#> (-3.567-xmin)/xres -> 98
#> (57.65-ymin)/yres -> 158
# Ronaldsway - Latitude: 54.07620. Longitude: -4.62333
#> (-4.62333-xmin)/xres -> 91
#> (54.07620-ymin)/yres -> 126
# Mumbles - Latitude: 51.567; Longitude: -3.983
#> (-3.983-xmin)/xres -> 95
#> (51.567-ymin)/yres -> 103
# The next eight stations are slightly aribtrary as I don't have Thompson's
# actually lats and lons for these. I'm also adding a site and making them
# more regular as I can with this dataset.
# Sites 10-18 are then:
# 45N 20W (45-ymin)/yres (-20-xmin)/xres -> 44 1,
# 45N 10W (45-ymin)/yres (-10-xmin)/xres -> 44 59,
# 45N 0W (45-ymin)/yres (0-xmin)/xres -> 44 119,
# 55N 10E (55-ymin)/yres (10-xmin)/xres -> 134 179,
# 65N 10E (65-ymin)/yres (10-xmin)/xres -> 224 179,
# 65N 10W (65-ymin)/yres (-10-xmin)/xres -> 224 59,
# 65N 20W (65-ymin)/yres (-20-xmin)/xres -> 224 1,
# 55N 20W (55-ymin)/yres (-20-xmin)/xres -> 134 1

polcoms<-new.env()
# extraMetStns in lon, lat order
polcoms$extraMetStns <- c(89,93,114,96,124,103,129,112,120,122,100,137,98,158,
  91,126,95,103,1,44,59,44,119,44,179,134,179,224,59,224,1,224,1,134)
#dim(polcoms$extraMetStns) <- c(2,17)


polcoms$x <- 1:198
polcoms$y <- 1:224
polcoms$lon <- seq(from=-19.83333, by=(1/6), length=198)
polcoms$lat <- seq(from=40.11111, by=(1/9), length=224)
polcoms$indices <- c(86,90,92,74,2,50,2,100,2,150, 50, 50, 50,100, 50,150, polcoms$extraMetStns)
dim(polcoms$indices)<-c(2,ntotstns)
polcoms$indices <- t(polcoms$indices)
polcoms$lonlat <- cbind(polcoms$lon[polcoms$indices[,1]], polcoms$lat[polcoms$indices[,2]])

# Need to find which point in HadSLP2 grid is closest to POLCOMS points
lonLat <- expand.grid(xlon,ylat)
xy <- expand.grid(x,y)

# This works with spam 0.15-5 but is a bit of a palaver.....
polcoms$near <- nearest.dist(polcoms$lonlat, lonLat, method="greatcircle", delta=100, upper=T)
gc()

polcoms$lenlonlat <- length(polcoms$lonlat[,1])
polcoms$nearest <- vector(length=polcoms$lenlonlat, mode="integer")
for (i in 1:polcoms$lenlonlat){
  polcoms$gtzero <- which(polcoms$near[i,]>0)
# There may be some points outside of the HadSLP2 region. If so these will
# return all zeros
  if (length(polcoms$gtzero>0)){
    polcoms$nearest[i] <- polcoms$gtzero[which.min(polcoms$near[i,polcoms$gtzero])]
  } else {
    polcoms$nearest[i] <- NA
  }
}

# Now extract these grid points from the pressure dataset
hadLonLat <- lonLat[polcoms$nearest,]
hadXY <- xy[polcoms$nearest,]

slpMonthlyArray<-array(NA,dim=c(nlon,nlat,nmon,nyr))

con<-file("~/diskx/HadSLP2r/hadSLP2_kij_1850-2007.bin","rb")

for (i in 1:nlat){
  for (j in 1:nlon){
    for (k in 1:nyr){
      slp<-readBin(con,what="numeric", size=4, n=nmon, endian='little')
      slpMonthlyArray[j,i,,k]<-slp
    }
  }
}

close(con)

dim(slpMonthlyArray)<-c(nlon,nlat,nmon*nyr)

# We now need to flip each grid so that it is the "right" way
# round when we come to extract the lats and lons

# Create time series of pressures for the nstns stations
slpNewlynStns <- array(0,dim=c((nmon*nyr),ntotstns))


for (count in 1:ntotstns){
    slpNewlynStns[,count] <- slpMonthlyArray[hadXY[count,1],hadXY[count,2],]
}

# Save the data
save(slpNewlynStns, file="brestNewlynHadSLP2.RData")
