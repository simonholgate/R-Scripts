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


lon<-seq(from=1/3, length=1080, by=1/3)
junk<-read.table('~/diskx/altimetry/lats.lis')
lat<-junk$V1
# jd (Julian Day) in this context appears to be days since 1950 
daytable<-read.table('~/diskx/altimetry/aviso/ref_h_anom/daytable.txt',
  col.names=c("seq","jd","year","month","day","yearday"), 
  colClasses=c("integer"))

l<-1080
#l<-360
m<-915
#m<-180
lm<-l*m
n<-818
#n<-513

meanArray <- vector(mode="numeric",length=n)
# Weight by area mask
weightVector <- cos(lat*pi/180)
weightMask <- array(weightVector, dim=c(m,l))
weightMask <- t(weightMask)
weightMask <- as.vector(weightMask)

# Read the time series that CWH has produced and see which grid
# points have only null values (=9999.0). These can then be flagged as
# land/ice points and removed from the time series.
# The format of hmaps (see the email Altimetry.txt) has longitude as the first axis, 
# latitude as the second and time as the third.
#inConn <- file("/diskx/users/simonh/altimetry/hstack_sm5_resid.dat", "rb")
#inConn <- file("~/diskx/altimetry/hmaps.dat", "rb")
inConn <- file("~/diskx/altimetry/aviso/ref_h_anom/ref_h_anom_stackmaps.dat", "rb")
for (i in 1:n) {
  data<-readBin(inConn, what="numeric", n = lm, size = 4)
#  data<-readBin(inConn, what="numeric", n = n, size = 4, endian="big")
  if (all(data==9999)) {
# Do nothing
  } else {
    message(as.character(i))
# We must delete any rows where data is missing 
    j <- which(data!=9999)
# Convert from cm to mm
    data <- data*10
# Weight by area
    meanArray[i] <- mean(data[j]*weightMask[j], na.rm = T)
  }
}
close(inConn)

x11()
plot(daytable$jd[1:n], meanArray, type='l', col='blue', lwd=2)

# Remove annual and semi annual signals
w1 <- 2*pi/(365)
w2 <- 2*pi/(365/2)
time <- daytable$jd[1:n]
mdata <- data.frame(meanArray=meanArray, t=time,   
  w1cos=cos(w1*time), w1sin=sin(w1*time), 
  w2cos=cos(w2*time), w2sin=sin(w2*time))

model <- lm(meanArray ~  w1cos + w1sin + w2cos + w2sin, data=mdata)
resid <- model$resid

# Calculate annual means from the residuals for 1993 to 2007
nyears<-14
annMeans <- vector(mode="numeric", length=nyears)
annSDs <- vector(mode="numeric", length=nyears)
annJD <- vector(mode="numeric", length=nyears)
for(i in 1:nyears){
  junk <- which(daytable$year==(i+1992))
  annMeans[i] <- mean(resid[junk], na.rm=T)
  annSDs[i] <- sd(resid[junk], na.rm=T)
  annJD[i] <- mean(daytable$jd[junk], na.rm=T)
}

# Calculate 5 year running mean trends
trendMidPoints <- c(1996:2006)
trends <- vector(mode='numeric', length=11)
for(i in 1:11){
  junk <- lm(annMeans[i:(i+4)] ~ c(1:5))
  trends[i] <- junk$coef[2]
}
x11()
plot(trendMidPoints, trends, type='l', lwd=2)
