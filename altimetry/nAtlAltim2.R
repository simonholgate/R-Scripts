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

# Load matrixMethods for use later
source("~/workspace/RScripts/matrixMethods.R")
# Newlyn 50 06 N  05 33 W => 50.1 -5.55
# Newlyn is at x=86, y=90
#> (-5.55-xmin)/xres
#[1] 85.69998
#> (50.1-ymin)/yres
#[1] 89.90001

# Use jet colors for images
#source('~/bin/RScripts/jet.colors.R')

# High res grid
lon<-seq(from=1/3, length=1080, by=1/3)
junk<-read.table('~/diskx/altimetry/lats.lis')
lat<-junk$V1
# jd (Julian Day) in this context appears to be days since 1950 
daytable<-read.table('~/diskx/altimetry/aviso/ref_h_anom/daytable.txt', 
  col.names=c("seq","jd","year","month","day","yearday"))

xres<-1/3
l<-1080
m<-915

# 1 deg grid
#xres<-1
#lon <- seq(from=1, to=360)
#lat <- seq(from=-90, to=90)

#l<-360
#m<-180

lm<-l*m

#n<-639
#n<-513
# Make n a factor of 4 so that we can reduce the array into pseudo months
#n<-512
n<-818

# Just worry about n time slices to begin with
# Focus on the N Atlantic covering 100W-15E, 5S-80N : [(781:1080,1:45),443:877]
#lonR <- c(260:360,1:15)
#latR <- 85:170
lonR <- c(781:1080,1:45)
latR <- 443:877
l1 <- length(lonR)
m1 <- length(latR)
lm1 <- l1*m1
lons <- seq(from=(lon[lonR[1]]-360),to=lon[lonR[l1]], by=xres)
dataArray<-array(NA,dim=c(l1,m1,n))

# Read the filtered time series that CWH has produced and see which grid
# points have only null values (=9999.0). These can then be flagged as
# land/ice points and removed from the time series.
# The format of hstack_resid (which has annual and semi annual signals removed
# - see the email Altimetry.txt) has time as the first axis, longitude as the
# second and latitude as the third.

# We'll use hmaps.dat in place of hstack.dat for EOFs
inConn <- file("~/diskx/altimetry/aviso/ref_h_anom/ref_h_anom_stackmaps.dat", "rb") # little endian
#inConn <- file("~/diskx/altimetry/hmaps_sm5_1deg.dat", "rb") # big endian
for (i in 1:n) {
  data<-readBin(inConn, what="numeric", n = lm, size = 4, endian='little')
#  data<-readBin(inConn, what="numeric", n = lm, size = 4, endian='big')
  dim(data)<-c(l,m)
  dataArray[,,i] <- data[lonR,latR]
}
close(inConn)

library(fields)

Z <- data[lonR,latR]
Z[which(Z==9999)] <- NA
#postscript(file="nAtl.ps")

par(family="HersheySans")
image.plot(lons, lat[latR], Z)
#dev.off()
stop("Imported data")

