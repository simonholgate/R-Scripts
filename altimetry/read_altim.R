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

# Newlyn 50 06 N  05 33 W => 50.1 -5.55
# Newlyn is at x=86, y=90
#> (-5.55-xmin)/xres
#[1] 85.69998
#> (50.1-ymin)/yres
#[1] 89.90001

library(fields)
# Use jet colors for images
source('~/bin/RScripts/jet.colors.R')
lon<-seq(from=1/3, length=1080, by=1/3)
junk<-read.table('~/diskx/altimetry/lats.lis')
lat<-junk$V1
# jd (Julian Day) in this context appears to be days since 1950 
daytable<-read.table('~/diskx/altimetry/daytable.lis', 
  col.names=c("seq","jd","year","month","day","yearday"))

l<-1080
m<-915
lm<-l*m
n<-639

lonLatArray<-array(NA,dim=c(lm,3))
lonLatArray[,1]<-rep(lon,m)
lonLatArray[,2]<-as.vector(t(array(rep(lat,l),dim=c(m,l))))

# Read the filtered time series that CWH has produced and see which grid
# points have only null values (=9999.0). These can then be flagged as
# land/ice points and removed from the time series.
# The format of hstack_resid (which has annual and semi annual signals removed
# - see the email Altimetry.txt) has time as the first axis, longitude as the
# second and latitude as the third.
inConn <- file("/diskx/users/simonh/altimetry/hstack_sm5_resid.dat", "rb")
for (i in 1:lm) {
  data<-readBin(inConn, what="numeric", n = n, size = 4)
  if (length(which(data==9999))) {
    lonLatArray[i,3]<-1
  } else {
    lonLatArray[i,3]<-0
  }
}
close(inConn)

plot(daytable$seq[1:639],data,type='l',lwd=2)

# The number of points which are not permanent land or ice
numSeaPoints<-length(which(lonLatArray[,3]==1))
# The lons and lats of the sea points
seaPointArray<-lonLatArray[which(lonLatArray[,3]==1),1:2]

