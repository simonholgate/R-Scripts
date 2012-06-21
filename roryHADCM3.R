#######################
# HADCM3 dimensions:  #
#                     #
# Xmin:   -1.25       #
# Xres:   1.25        #
# Xn:     288         #
#                     #
# Ymin:  -90.625      #
# Yres:   1.25        #
# Yn:     144         #
#                     #
# Timeslices: 1200    #
#                     #
#######################

# Newlyn 50 06 N  05 33 W => 50.1 -5.55
# Newlyn is at x=86, y=90
#> (-5.55-xmin)/xres
#[1] 85.69998
#> (50.1-ymin)/yres
#[1] 89.90001

library(fields)
# Use jet colors for images
source('~/bin/RScripts/jet.colors.R')


lon<-seq(from=-1.25, length=288, by=1.25)
lat<-seq(from=-90.625, length=144, by=1.25)

l<-288
m<-144
lm<-l*m
n<-1200

lonLatArray<-array(NA,dim=c(lm,3))
lonLatArray[,1]<-rep(lon,m)
lonLatArray[,2]<-as.vector(t(array(rep(lat,l),dim=c(m,l))))

# Read the filtered time series that Rory has produced and see which grid
# points have only null values (=-1.9e19). These can then be flagged as
# land points and removed from the time series.
# The format of shma_str.dat (which has annual and semi annual signals
# removed) has longitude as the first axis, latitude as the
# second and time as the third.
inConn <- file("/diskx/users/rjbi/share/simon/data/sshma_str.dat", "rb")
data<-readBin(inConn, what="numeric", n = lm, size = 4)
close(inConn)
## A land point = -1.9e19
j<-which(data < -1.0e19)
lonLatArray[j,3]<-1
## Not a land point
j<-which(data > -1.0e19)
lonLatArray[j,3]<-0

#for (i in 1:lm) {
## A land point = -1.9e19
#  if (identical(all.equal(data[i],-1.9e19),TRUE)) {
#    lonLatArray[i,3]<-1
## Not a land point
#  } else {
#    lonLatArray[i,3]<-0
#  }
#}

j<-which(lonLatArray[,3]==1)
# The number of points which are not permanent land or ice
numSeaPoints<-length(j)
# The lons and lats of the sea points
seaPointArray<-lonLatArray[j,1:2]

data[j]<-NA
map<-array(data,dim=c(288,144))
image.plot(lon,lat,map,zlim=c(-40,40),col=jet.colors(100))

map2<-array(NA,dim=c(288,144))
map2[1:144,]<-map[145:288,]
map2[145:288,]<-map[1:144,]
lon2<-lon
lon2[1:144]<-lon[145:288]-360
lon2[145:288]<-lon[1:144]
image.plot(lon2,lat,map2,col=jet.colors(100),zlim=c(-40,40))

NAtlanticMask<-array(NA,dim=c(288,144))
i<-intersect(which(lon2>=-100),which(lon2<=15))
j<-intersect(which(lat>=0),which(lat<=75))
NAtlanticMask[i,j]<-map2[i,j]
image.plot(lon2[i],lat[j],NAtlanticMask[i,j],col=jet.colors(100),zlim=c(-40,40))

k<-which(is.finite(NAtlanticMask))
NAtlanticMask[k]<-1


