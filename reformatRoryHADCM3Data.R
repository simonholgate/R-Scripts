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

library(fields)
# Use jet colors for images
source('~/bin/RScripts/jet.colors.R')

lon<-seq(from=-1.25, length=288, by=1.25)
lat<-seq(from=-90.625, length=144, by=1.25)

l<-288
m<-144
lm<-l*m
n<-1200

# Read the filtered time series that Rory has produced and see which grid
# points have only null values (=-1.9e19). These can then be flagged as
# land points and removed from the time series.
# The format of shma_str.dat (which has annual and semi annual signals
# removed) has longitude as the first axis, latitude as the
# second and time as the third.
inConn <- file("/diskx/users/rjbi/share/simon/data/sshma_str.dat", "rb")
data<-readBin(inConn, what="numeric", n = lm*n, size = 4)
close(inConn)
#
print("Read data")

# Reformat into a 3d array
timeseries<-array(data,dim=c(lm,n))
#
print("Reformatted data")

# Free some memory
rm(data)
#
print("Freed data")

# Transpose 
timeseries<-t(timeseries)
#
print("Transposed data")

# Write out a timeseries array
dim(timeseries)<-c(lm*n,1)
timeseries<-as.vector(timeseries)
outConn <- file("/diskx/users/simonh/roryHADCM3/ssh_stack.dat", "wb")
for (i in 1:n) {
  writeBin(timeseries[((i*lm)-lm+1):(i*lm)],outConn,size=4)
}
close(outConn)
#
print("Wrote data")
