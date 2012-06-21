# Extracts the altimetry time series over the NW European shelf for use in Altim
# and comparison with POLCOMMS

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

library(fields)
# Use jet colors for images
source('~/bin/RScripts/jet.colors.R')
lon<-seq(from=1/3, length=1080, by=1/3)
junk<-read.table('~/scratch/data/altimetry/lats.lis')
lat<-junk$V1
# jd (Julian Day) in this context appears to be days since 1950 
daytable<-read.table('~/scratch/data/altimetry/daytable.lis', 
  col.names=c("seq","jd","year","month","day","yearday"))

l<-1080
m<-915
lm<-l*m
n<-639

lonLatArray<-array(NA,dim=c(lm,3))
lonLatArray[,1]<-rep(lon,m)
lonLatArray[,2]<-as.vector(t(array(rep(lat,l),dim=c(m,l))))

nwLonLatArray<-array(NA,dim=c(99*127,3))
nwDataArray<-array(NA,dim=c(639,99*127))
# Read the filtered time series that CWH has produced and see which grid
# points have only null values (=9999.0). These can then be flagged as
# land/ice points and removed from the time series.
# The format of hstack_resid (which has annual and semi annual signals removed
# - see the email Altimetry.txt) has time as the first axis, longitude as the
# second and latitude as the third.
#inConn <- file("~/diskx/altimetry/hstack_resid.dat", "rb")
inConn <- file("~/diskx/altimetry/hstack_sm5_resid.dat", "rb")
j<-1
for (i in 1:lm) {
  if ((lonLatArray[i,1]>=340.1666)||(lonLatArray[i,1]<=13.0)){
    if ((lonLatArray[i,2]>=40.1111)&&(lonLatArray[i,2]<=64.88889)){
      nwLonLatArray[j,1]<-lonLatArray[i,1]
      nwLonLatArray[j,2]<-lonLatArray[i,2]
      data<-readBin(inConn, what="numeric", n = n, size = 4)
# Check for a "sea" point. If all the data is 9999 then which(data!=9999)
# will have a length of 0 and the if statement evaluates to FALSE
      if (length(which(data!=9999))) {
        nwLonLatArray[j,3]<-1
        nwDataArray[,j]<-data
      } else {
        nwLonLatArray[j,3]<-0
      }
      j<-j+1
    } else {
# Not reading data so need to skip forward by 639*4 bytes
#      seek(inConn,origin="current",where=(639*4))
      data<-readBin(inConn, what="numeric", n = n, size = 4)
    }
  } else {
# Not reading data so need to skip forward by 639*4 bytes
#    seek(inConn,origin="current",where=(639*4))
    data<-readBin(inConn, what="numeric", n = n, size = 4)
  }
}
close(inConn)

# The number of points which are not permanent land or ice
numNWSeaPoints<-length(which(nwLonLatArray[,3]==1))
# The lons and lats of the sea points
nwSeaPointArray<-nwLonLatArray[which(nwLonLatArray[,3]==1),1:2]
# Convert Longitudes to degrees W instead of 360
j <- which(nwSeaPointArray[,1]>300)
nwSeaPointArray[j,1] <- nwSeaPointArray[j,1]-360

#nwAltimDataArray<-nwDataArray[,which(nwLonLatArray[,3]==1)]
nwAltimDataArray<-nwDataArray
nwAltimDataArray[which(nwAltimDataArray==9999.0)]<-0

nwZerosArray<-array(NA,dim=c(639,1))
for (i in 1:639){
  j<-which(nwAltimDataArray[i,]==0)
  if (length(j)){
    nwZerosArray[i]<-0
  } else {
    nwZerosArray[i]<-1
  }
}

# Let's take the 4 weekly mean of the data to get rid of some blanks
# Put NAs in place of zeros
nwAltimDataArray[which(nwAltimDataArray==0)]<-NA
#smNWAltimDataArray<-array(NA,dim=c(158,length(which(nwLonLatArray[,3]==1))))
smNWAltimDataArray<-array(NA,dim=c(158,length(nwLonLatArray[,3])))
for (i in 1:158){
  for (j in 1:length(nwLonLatArray[,3])){
    smNWAltimDataArray[i,j]<-mean(nwAltimDataArray[(i*4-3):(i*4),j],na.rm=TRUE)
  }
}
smNWAltimDataArray[which(is.na(smNWAltimDataArray))]<-0

nwEurope <- array(NA,dim=c(99,127,639))

for (ii in 1:639){

# Take month 1, January, and create a 2D array for visualising
  jan <- nwDataArray[1,]

  dim(jan) <- c(99,127)
  ice<-which(jan==9999)
  jan[ice]<-NA
  jann <- jan[c(40:99,1:39),]

  nwEurope[,,ii] <- jann
}

nwlat <- seq(from=40.1111, length=127, by=1/3)
nwlon <- seq(from=(340.1666-360), length=99, by=1/3)
# Visualise some data...
image.plot(nwlon,nwlat,nwEurope[,,1])
