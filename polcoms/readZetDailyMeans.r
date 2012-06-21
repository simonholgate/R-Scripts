###############################################################################
# Variation of readzet.r to read the 2D daily mean sea level calculated from
# POLCOMMS model output. 
# 
# The arrays were written from R with the makeZetDailyMeans.r script and are
# little endian real*4 numbers. Each month contains (day in month - 2) daily
# means.
#
# Simon Holgate, August 2005
#
###############################################################################

#######################
# _*S2 dimensions:    #
#                     #
# *_Xmin:   -19.83333 #
# Xres:   1/6         #
# Xn:     198         #
#                     #
# Ymin:   40.11111    #
# Yres:   1/9         #
# Yn:     224         #
#######################

# Newlyn 50 06 N  05 33 W => 50.1 -5.55
# Newlyn is at x=86, y=90
#> (-5.55-xmin)/xres
#[1] 85.69998
#> (50.1-ymin)/yres
#[1] 89.90001

# Lerwick 60 09 N 01 08 W => 60.15 -1.13
# Lerwick is at x=113, y=181
#> (-1.13-xmin)/xres
#[1] 112.2200
#> (60.15-ymin)/yres
#[1] 180.35

library(fields)
# Use jet colors for images
source('~/bin/RScripts/jet.colors.R')
lat<-seq(from=40.11111,by=1/9,length=224)
lon<-seq(from=-19.83333,by=1/6,length=198)

l<-198
m<-224
lm<-l*m
n<-20
#
zet<-array(0,dim=c(l,m,39))
# Weights for Doodson filter
weightArray<-array(
 c(1,0,1,0,0,1,0,1,1,0,2,0,1,1,0,2,1,1,2,0,2,1,1,2,0,1,1,0,2,0,1,1,0,1,0,0,1,0,1),dim=c(39,1))
# We're not interested in u or v so just use one array
#uvb<-array(0,dim=c(l,m))

# Number of months to iterate over this will be determined by the number of
# file names in inname
inname<-c( 
"zet_UBVB.S12run401.jan60",
"zet_UBVB.S12run401.feb60",
"zet_UBVB.S12run401.mar60",
"zet_UBVB.S12run401.apr60",
"zet_UBVB.S12run401.may60",
"zet_UBVB.S12run401.jun60",
"zet_UBVB.S12run401.jul60",
"zet_UBVB.S12run401.aug60",
"zet_UBVB.S12run401.sep60",
"zet_UBVB.S12run401.oct60",
"zet_UBVB.S12run401.nov60",
"zet_UBVB.S12run401.dec60" 
)

leapYrDaysArray<-array(data=
  c(31,29,31,30,31,30,31,31,30,31,30,31),
  dim=c(12,1))
daysArray<-array(data=
  c(31,28,31,30,31,30,31,31,30,31,30,31),
  dim=c(12,1))
yearsArray<-array(data=
  c(1960,1961,1962,1963,1964,1965,1966,1967,1968,1969),
  dim=c(10,1))
#monthArray<-array(data=
#  c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),
#  dim=c(12,1))
monthArray<-array(data=
  c("01","02","03","04","05","06","07","08","09","10","11","12"),
  dim=c(12,1))
newlynMonthlyMean<-array(NA,dim=c(12,1))
lerwickMonthlyMean<-array(NA,dim=c(12,1))

# Test for leap year
if ((yearsArray[1] %% 4)==0) {
  if ((yearsArray[1] %% 100)==0) {
    if ((yearsArray[1] %% 400)==0) {
      leapFlag<-1
    } else {
      leapFlag<-0
    }
  } else {
    leapFlag<-1
  }
  leapFlag<-1
} else {
  leapFlag<-0
}

nmons<-length(inname)

zettMonthlyMean<-array(NA,dim=c(l,m,12))
zettDailyMeanArray<-array(NA,dim=c(l,m,31))

outMonthlyMean <- file(paste("zett.mm.1960.dat",sep=""), "wb")

for (ii in 1:nmons) {

  inConn <- file(inname[ii], "rb")
  outDailyMean <- file(
    paste("dailyMeans/zett.dm.",monthArray[ii],".1960.dat",sep=""), "wb")

  if ( leapFlag == 1) {
  # Leap year
    hours<-24*leapYrDaysArray[ii]
    days<-leapYrDaysArray[ii]
  } else {
  # Regular year
    hours<-24*daysArray[ii]
    days<-daysArray[ii]
  }
  # Number of time slices to read. Output is once per hour.
  ntimes<-hours

  newlyn<-array(NA,dim=c(ntimes,1))
  lerwick<-array(NA,dim=c(ntimes,1))

####################
# Read header info #
####################
# Read junk first integer that Fortran writes
  junk1<-readBin(inConn, what="integer", n = 1, size = NA, endian = "big")
# Read data
  junk<-readBin(inConn, what="integer", n = 4, size = NA, endian = "big")
# Read junk last integer that Fortran writes
  junk1<-readBin(inConn, what="integer", n = 1, size = NA, endian = "big")

  ll<-junk[1]
  mm<-junk[2]
  nn<-junk[3]
  npsea<-junk[4]
  
  isea<-array(NA,dim=npsea)
  jsea<-array(NA,dim=npsea)

# Read junk first integer that Fortran writes
  junk1<-readBin(inConn, what="integer", n = 1, size = NA, endian = "big")
# Read data
  isea<-readBin(inConn, what="integer", n = npsea, size = NA, endian = "big")
# Read junk last integer that Fortran writes
  junk1<-readBin(inConn, what="integer", n = 1, size = NA, endian = "big")

# Read junk first integer that Fortran writes
  junk1<-readBin(inConn, what="integer", n = 1, size = NA, endian = "big")
# Read data
  jsea<-readBin(inConn, what="integer", n = npsea, size = NA, endian = "big")
# Read junk last integer that Fortran writes
  junk1<-readBin(inConn, what="integer", n = 1, size = NA, endian = "big")

# Read junk first integer that Fortran writes
  junk1<-readBin(inConn, what="integer", n = 1, size = NA, endian = "big")
# Read data
  junk<-readBin(inConn, what="integer", n = 4, size = NA, endian = "big")
# Read junk last integer that Fortran writes
  junk1<-readBin(inConn, what="integer", n = 1, size = NA, endian = "big")

  ll<-junk[1]
  mm<-junk[2]
  nn<-junk[3]
  npusea<-junk[4]

  iusea<-array(NA,dim=npusea)
  jusea<-array(NA,dim=npusea)

# Read junk first integer that Fortran writes
  junk1<-readBin(inConn, what="integer", n = 1, size = NA, endian = "big")
# Read data
  iusea<-readBin(inConn, what="integer", n = npusea, size = NA, endian = "big")
# Read junk last integer that Fortran writes
  junk1<-readBin(inConn, what="integer", n = 1, size = NA, endian = "big")

# Read junk first integer that Fortran writes
  junk1<-readBin(inConn, what="integer", n = 1, size = NA, endian = "big")
# Read data
  jusea<-readBin(inConn, what="integer", n = npusea, size = NA, endian = "big")
# Read junk last integer that Fortran writes
  junk1<-readBin(inConn, what="integer", n = 1, size = NA, endian = "big")

  zett<-array(NA,dim=npsea)
# We're not interested in u or v so just use one array
  uvbt<-array(NA,dim=npusea)

####################
# Read data arrays #
####################
  for (jk in 1:ntimes) {
# Read junk first integer that Fortran writes
    junk1<-readBin(inConn, what="integer", n = 1, size = NA, endian = "big")
# Read data
    itimt<-readBin(inConn, what="integer", n = 1, size = NA, endian = "big")
# Read junk last integer that Fortran writes
    junk1<-readBin(inConn, what="integer", n = 1, size = NA, endian = "big")

# Read junk first integer that Fortran writes
    junk1<-readBin(inConn, what="integer", n = 1, size = NA, endian = "big")
# Read sea level data
    zett<-readBin(inConn, what="numeric", n = npsea, size = 4, endian = "big")
# Read junk last integer that Fortran writes
    junk1<-readBin(inConn, what="integer", n = 1, size = NA, endian = "big")

# Read junk first integer that Fortran writes
    junk1<-readBin(inConn, what="integer", n = 1, size = NA, endian = "big")
# Read data
    uvbt<-readBin(inConn, what="numeric", n = npusea, size = 4, endian = "big")
# Read junk last integer that Fortran writes
    junk1<-readBin(inConn, what="integer", n = 1, size = NA, endian = "big")

# Read junk first integer that Fortran writes
    junk1<-readBin(inConn, what="integer", n = 1, size = NA, endian = "big")
# Read data
    uvbt<-readBin(inConn, what="numeric", n = npusea, size = 4, endian = "big")
# Read junk last integer that Fortran writes
    junk1<-readBin(inConn, what="integer", n = 1, size = NA, endian = "big")

# Check for whether the first few arrays can be skipped
    if ((jk <= 17 ) | (jk >= (ntimes-17))) {
      timeLev <- 1
    } else {
      if (timeLev < 39) {
        timeLev <- timeLev+1
      } else {
# Filter the tides from the data
        zetMean<-array(0,dim=c(l,m))
        for (ip in 1:npsea) {
          if (!is.nan(zet[isea[ip],jsea[ip],1])){
            zetMean[isea[ip],jsea[ip]] <- 
              sum(zet[isea[ip],jsea[ip],]*weightArray)/30
          } else {
            zetMean[isea[ip],jsea[ip]] <- NA
          }
          writeBin(zetMean[isea[ip],jsea[ip]], outDailyMean)
        }
# Write out mean to file
# Place mean in daily mean array
        zettDailyMeanArray[,,ii] <- zetMean
# Shift the usable "old" data to the beginning of the array
        old <- array(0,dim=c(l,m,39))
        old[,,c(1,15)]<-zet[,,c(25,39)]
        zet <- old
        timeLev <- 16
      }
    }
    for (ip in 1:npsea) {
      zet[isea[ip],jsea[ip],timeLev]<-zett[ip]
    }

# Produce time series at Newlyn
    newlyn[jk]<-zet[86,90,timeLev]
    lerwick[jk]<-zet[113,180,timeLev]
#    print(paste(jk,zet[86,90],newlyn[jk]))

# 
    zet[which(zet==0)]<-NA

#    if (nchar(as.character(jk))<2){
#      jkstr<-paste("00",jk,sep="")
#    } else {
#      if (nchar(as.character(jk))<3){
#      jkstr<-paste("0",jk,sep="")
#      } else {
#      jkstr<-paste(jk,sep="")
#      }
#    }
#    postscript(
#      file=paste("/local/simonhdata/data/polcomms/ps/zett",monthArray[ii],"1960.",jkstr,".ps",
#      sep=""))
#      image.plot(lon,lat,zet[c(1:198),c(1:224)],zlim=c(-4,4),
#        col=jet.colors(100))
#    dev.off()

  }
#
# Doodson filter on hourly values
# See http://www.pol.ac.uk/psmslh/gloup/doodson_X0.html for details
# Running weighted mean of 19 values either side of a central value
# Weights are (1010010110201102112 0 2112011020110100101)/30
  filtNewlyn<-array(NA,dim=c(hours,1))
  filtLerwick<-array(NA,dim=c(hours,1))
  meanArray<-array(NA,dim=c(39,1))
  weightArray<-array(
  c(1,0,1,0,0,1,0,1,1,0,2,0,1,1,0,2,1,1,2,0,2,1,1,2,0,1,1,0,2,0,1,1,0,1,0,0,1,0,1),dim=c(39,1))
  for (i in 1:(hours-38)) {
    newlynMeanArray<-(newlyn[i:(i+38),1]*weightArray)
    lerwickMeanArray<-(lerwick[i:(i+38),1]*weightArray)
    filtNewlyn[i+19]<-sum(newlynMeanArray)/30
    filtLerwick[i+19]<-sum(lerwickMeanArray)/30
  }
  dayMeanNewlyn<-array(NA,c((days-2),1))
  dayMeanLerwick<-array(NA,c((days-2),1))
  middays<-seq(36,hours,24)
  for (i in 1:(days-2)) {
    dayMeanNewlyn[i]<-filtNewlyn[middays[i]]
    dayMeanLerwick[i]<-filtLerwick[middays[i]]
  }

  newlynMonthlyMean[ii] <- sum(dayMeanNewlyn)/(days-2)
  lerwickMonthlyMean[ii] <- sum(dayMeanLerwick)/(days-2)
  zettMonthlyMean[ii] <- sum(zettDailyMeanArray, na.rm=TRUE)/(days-2)

  close(inConn)
  close(outDailyMean)

}
# End of data reading

# Write monthly means out
for (ip in 1:npsea) {
  writeBin(zettMonthlyMean[isea[ip],jsea[ip]],outMonthlyMean)
}
close(outMonthlyMean)

#
zet[which(zet==0)]<-NA
image.plot(lon,lat,zet[c(1:198),c(1:224)],zlim=c(-4,4),
  col=jet.colors(100))
