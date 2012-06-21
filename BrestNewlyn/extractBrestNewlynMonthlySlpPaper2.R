# Find the nearest HadSLP2r grid box for a number of sites
library(fields)

ntotstns <- (2+35)
nlon<-72
nlat<-37
nyr<-160
nmon<-12
xlon<-seq(from=-180,by=5,length=nlon)
xlon2 <- seq(from=180,by=5,length=nlon)
ylat<-seq(from=90,by=-5,length=nlat)
x<-seq(from=1,length=nlon)
y<-seq(from=1,length=nlat)

# Newlyn 50 06 N  05 33 W => 50.1 -5.55
# Brest 48 23 N  04 30 W => 48.38 -4.5
# Add extra stations from around the model domain here similar to Thompson (1980)
# The 35 stations are on a regular 5 deg. grid covering 40N-60N, 20W-10E
# as I can do that with this dataset. 

location<-new.env()

location$lat <- c(50, 50,
                  40, 40, 40, 40, 40, 40, 40,
                  45, 45, 45, 45, 45, 45, 45,
                  50, 50, 50, 50, 50, 50, 50, 
                  55, 55, 55, 55, 55, 55, 55, 
                  60, 60, 60, 60, 60, 60, 60)
location$lon <- c(-5, -5,
                  -20, -15, -10, -5, 0, 5, 10,
                  -20, -15, -10, -5, 0, 5, 10,
                  -20, -15, -10, -5, 0, 5, 10,
                  -20, -15, -10, -5, 0, 5, 10,
                  -20, -15, -10, -5, 0, 5, 10) 

location$lonlat <- cbind(location$lon, location$lat)

location$nearest <- vector(mode="integer", length=ntotstns)
# Not working properly nearest.dist so....
for (i in 1:ntotstns){
  location$nearest[i] <- intersect(which(lonLat[,1]==location$lonlat[i,1]), which(lonLat[,2]==location$lonlat[i,2]))
}

# Now extract these grid points from the pressure dataset
hadLonLat <- lonLat[location$nearest,]
hadXY <- xy[location$nearest,]

slpMonthlyArray<-array(NA,dim=c(nlon,nlat,nmon,nyr))

con<-file("~/diskx/HadSLP2r/hadSLP2_kij_1850-2009.bin","rb")

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
save(slpNewlynStns, file="brestNewlynHadSLP2Paper2.RData")
