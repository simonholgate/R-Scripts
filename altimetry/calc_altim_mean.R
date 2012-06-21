###########################################################
# Script to calculate annual RMS values of sea level from #
# altimetry                                               #
###########################################################

# 1. We're going to read in the altimetry file in time series order rather than
# map order
# 2. We'll calculate annual averages from the time series
# 3. From the vector of annual means, we'll remove the mean and any linear trend
# 4. The root-mean-square of the residuals will be calculated and stored in a
# matrix to be displayed as an image

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
# Timeslices: 818     #
#                     #
#######################

l<-1080
#l<-360
m<-915
#m<-180
lm<-l*m
n<-818
#n<-513
nstack <- 1000

lon<-seq(from=1/3, length=l, by=1/3)
junk<-read.table('~/diskx/altimetry/lats.lis')
lat<-junk$V1
# jd (Julian Day) in this context appears to be days since 1950 
daytable<-read.table('~/diskx/altimetry/aviso/ref_h_anom/daytable.txt',
  col.names=c("seq","jd","year","month","day","yearday"), 
  colClasses=c("integer"))

time <- seq.Date(from=as.Date(
  paste(daytable$year[1],daytable$month[1],daytable$day[1], sep="/")),
                 to=as.Date(
  paste(daytable$year[n],daytable$month[n],daytable$day[n], sep="/")),
                 length.out=n)

# We're only interested in complete years
yearEnds <- array(NA, dim=c(15,3))
for (i in 1:15){
  yearEnds[i,1] <- max(which(time<=as.Date(paste((1992+i),"/1/1",sep=""))))+1
  yearEnds[i,2] <- max(which(time<=as.Date(paste((1992+i),"/12/31",sep=""))))
  yearEnds[i,3] <- length(yearEnds[i,1]:yearEnds[i,2])
}
dataLength <- length(yearEnds[1,1]:yearEnds[15,2])

meanArray <- array(NA, dim=c(l,m))

# Read the time series that CWH has produced and see which grid
# points have only null values (=9999.0). These can then be flagged as
# land/ice points and removed from the time series.
# The format of hmaps (see the email Altimetry.txt)
# has longitude as the first axis, 
# latitude as the second and time as the third.
# In contrast, stacktime has time as the first dimension
inConn <- file("~/diskx/altimetry/aviso/ref_h_anom/ref_h_anom_stacktime.dat",
               "rb")

for (i in 1:m) {
  
  message(paste(i))
  
  for (j in 1:l) {
    data <- readBin(inConn, what="numeric", n = nstack, size = 4)
    
# Just select the time slice we're interested in
    slice <- data[yearEnds[1,1]:yearEnds[15,2]]
    
# Set a limit of 10% of missing data in total
    if (length(which(slice==9999))>=(dataLength*.1)) {
# Do nothing
      next
    } else {
      
# We must set any rows where data is missing to NA 
      missing <- which(data==9999)
      data[missing] <- NA
      
# Convert from cm to mm
      data <- data*10
      
      meanArray[j,i] <- mean(data, na.rm=T)
    }
  }
}

close(inConn)

x11()
library(fields)
image.plot(lon, lat, meanArray)
