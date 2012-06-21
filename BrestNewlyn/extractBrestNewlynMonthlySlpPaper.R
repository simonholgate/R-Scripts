#
library(fields)

nstns<-2
nextrastns<-30
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

# Newlyn is at (86,90) and Brest is at (92,74)

# Add extra stations from around the model domain here similar to Thompson (1980)
# The 35 stations are on a regular 5 deg. grid covering 40N-55N, 10W-5E
# as I can do that with this dataset.
# Sites 3-37 are then:
# 42.5N 17.5W c((42.5-ymin)/yres, (-17.5-xmin)/xres) -> 22 14
# 42.5N 12.5W c((42.5-ymin)/yres, (-12.5-xmin)/xres) -> 22 44
# 42.5N 7.5W c((42.5-ymin)/yres, (-7.5-xmin)/xres) -> 22 74
# 42.5N 2.5W c((42.5-ymin)/yres, (-2.5-xmin)/xres) -> 22 104
# 42.5N 2.5E c((42.5-ymin)/yres, (2.5-xmin)/xres) -> 22 134
# 42.5N 7.5E c((42.5-ymin)/yres, (7.5-xmin)/xres) -> 22 164
# 42.5N 12.5E c((42.5-ymin)/yres, (12.5-xmin)/xres) -> 22 194
# 47.5N 17.5W c((47.5-ymin)/yres, (-17.5-xmin)/xres) -> 67 14
# 47.5N 12.5W c((47.5-ymin)/yres, (-12.5-xmin)/xres) -> 67 44
# 47.5N 7.5W c((47.5-ymin)/yres, (-7.5-xmin)/xres) -> 67 74
# 47.5N 2.5W c((47.5-ymin)/yres, (-2.5-xmin)/xres) -> 67 104
# 47.5N 2.5E c((47.5-ymin)/yres, (2.5-xmin)/xres) -> 67 134
# 47.5N 7.5E c((47.5-ymin)/yres, (7.5-xmin)/xres) -> 67 164
# 47.5N 12.5E c((47.5-ymin)/yres, (12.5-xmin)/xres) -> 67 194
# 52.5N 17.5W c((52.5-ymin)/yres, (-17.5-xmin)/xres) -> 112 14
# 52.5N 12.5W c((52.5-ymin)/yres, (-12.5-xmin)/xres) -> 112 44
# 52.5N 7.5W c((52.5-ymin)/yres, (-7.5-xmin)/xres) -> 112 74
# 52.5N 2.5W c((52.5-ymin)/yres, (-2.5-xmin)/xres) -> 112 104
# 52.5N 2.5E c((52.5-ymin)/yres, (2.5-xmin)/xres) -> 112 134
# 52.5N 7.5E c((52.5-ymin)/yres, (7.5-xmin)/xres) -> 112 164
# 52.5N 12.5E c((52.5-ymin)/yres, (12.5-xmin)/xres) -> 112 194
# 57.5N 17.5W c((57.5-ymin)/yres, (-17.5-xmin)/xres) -> 157 14
# 57.5N 12.5W c((57.5-ymin)/yres, (-12.5-xmin)/xres) -> 157 44
# 57.5N 7.5W c((57.5-ymin)/yres, (-7.5-xmin)/xres) -> 157 74
# 57.5N 2.5W c((57.5-ymin)/yres, (-2.5-xmin)/xres) -> 157 104
# 57.5N 2.5E c((57.5-ymin)/yres, (2.5-xmin)/xres) -> 157 134
# 57.5N 7.5E c((57.5-ymin)/yres, (7.5-xmin)/xres) -> 157 164
# 57.5N 12.5E c((57.5-ymin)/yres, (12.5-xmin)/xres) -> 157 194
# 62.5N 17.5W c((62.5-ymin)/yres, (-17.5-xmin)/xres) -> 202 14
# 62.5N 12.5W c((52.5-ymin)/yres, (-12.5-xmin)/xres) -> 202 44
# 62.5N 7.5W c((62.5-ymin)/yres, (-7.5-xmin)/xres) -> 202 74
# 62.5N 2.5W c((62.5-ymin)/yres, (-2.5-xmin)/xres) -> 202 104
# 62.5N 2.5E c((62.5-ymin)/yres, (2.5-xmin)/xres) -> 202 134
# 62.5N 7.5E c((62.5-ymin)/yres, (7.5-xmin)/xres) -> 202 164
# 62.5N 12.5E c((62.5-ymin)/yres, (12.5-xmin)/xres) -> 202 194
polcoms<-new.env()

# extraMetStns in lon, lat order
polcoms$extraMetStns <- c(90,86,74,92,
                          44,22,74,22,104,22,134,22,164,22,194,22,
                          44,67,74,67,104,67,134,67,164,67,194,67,
                          44,112,74,112,104,112,134,112,164,112,194,112,
                          44,157,74,157,104,157,134,157,164,157,194,157,
                          44,202,74,202,104,202,134,202,164,202,194,202)
dim(polcoms$extraMetStns) <- c(2,ntotstns)

polcoms$x <- 1:198
polcoms$y <- 1:224
polcoms$lon <- seq(from=-19.83333, by=(1/6), length=198)
polcoms$lat <- seq(from=40.11111, by=(1/9), length=224)
polcoms$indices <- polcoms$extraMetStns
#dim(polcoms$indices)<-c(2,ntotstns)
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
save(slpNewlynStns, file="brestNewlynHadSLP2Paper.RData")
