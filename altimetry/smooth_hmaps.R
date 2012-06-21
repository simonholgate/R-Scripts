# Smooths the ENACT altimetry maps using FFT

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

# The format of hmaps.dat (see the email Altimetry.txt) has time as the last 
# axis, with longitude first and latitude second.
inConn <- file("/diskx/users/cwh/msla2/hmaps.dat", "rb")
j<-1
#for (i in 1:n){
  data<-readBin(inConn, what="numeric", n = lm, size = 4)
  slArray<-array(data,dim=c(1080,915))
# Linearly interpolate missing values zonally. If all values are 9999 then
# ignore as that won't be smoothed anyway
  for (lats in 1:m){
    missDataIndices<-which(slArray[,lats]==9999)
    if (length(missDataIndices)){# There is some missing data, otherwise there 
				 # is no interpolation to do
      dataIndices<-which(slArray[,lats]!=9999)
      if (length(dataIndices)){# There is at least some data, otherwise there 
			       # is no interpolation to do either 
			       # - leave as 9999
# If there is only one data value then the whole zone takes on that value
# Otherwise wee interpolate between the values taking care to wrap around the
# globe where necessary
        if ((length(j)==1){
          slArray[,lats]<-slArray[j,lats]
        } else {
          for (gap in 1:(length(j)-1)) {
            x<-c(slArray[dataIndices[gap],lats],
                 slArray[dataIndices[gap+1],lats])
            approx(x,n=(dataIndices[gap+1]-dataIndices[gap]-1))

  image.plot(lon,lat,slArray)
#}
close(inConn)
