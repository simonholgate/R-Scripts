###############################################################################
# Variation of readzet.r which calculates the 2D daily mean sea level from the
# POLCOMMS model output. This requires reading 39 2D arrays (19 either side of
# the midday of the day for which the mean id being calculated).
#
# With the first hour of model output being 0100 hours, we must skip the 1st
# 16 2D sea level arrays before reading arrays 17-55, filtering those with the
# Doodson filter before reading another 24 arrays in and calculating the next
# daily mean array
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
# intersect(which(isea==86),which(jsea==90)) => ip=17701
#> (-5.55-xmin)/xres
#[1] 85.69998
#> (50.1-ymin)/yres
#[1] 89.90001

#library(fields)
# Use jet colors for images
source('~/bin/RScripts/jet.colors.R')
lat<-seq(from=40.11111,by=1/9,length=224)
lon<-seq(from=-19.83333,by=1/6,length=198)

l<-198
m<-224
lm<-l*m
n<-20
#
# Weights for Doodson filter
weightArray<-array(
 c(1,0,1,0,0,1,0,1,1,0,2,0,1,1,0,2,1,1,2,0,2,1,1,2,0,1,1,0,2,0,1,1,0,1,0,0,1,0,1),dim=c(39,1))

# Number of months to iterate over this will be determined by the number of
# file names in inname
inname<-c( 
"zet_UBVB.S12run405.jan",
"zet_UBVB.S12run405.feb",
"zet_UBVB.S12run405.mar",
"zet_UBVB.S12run405.apr",
"zet_UBVB.S12run405.may",
"zet_UBVB.S12run405.jun",
"zet_UBVB.S12run405.jul",
"zet_UBVB.S12run405.aug",
"zet_UBVB.S12run405.sep",
"zet_UBVB.S12run405.oct",
"zet_UBVB.S12run405.nov",
"zet_UBVB.S12run405.dec" 
)

leapYrDaysArray<-array(data=
  c(31,29,31,30,31,30,31,31,30,31,30,31),
  dim=c(12,1))
daysArray<-array(data=
  c(31,28,31,30,31,30,31,31,30,31,30,31),
  dim=c(12,1))
yearsArray<-array(data=
  c(1990), dim=c(1,1))
#  c(1990,1991,1992,1993,1994,1995,1996,1997,1998,1999), dim=c(10,1))
#  c(1960,1961,1962,1963,1964,1965,1966,1967,1968,1969,
#  1970,1971,1972,1973,1974,1975,1976,1977,1978,1979,
#  1980,1981,1982,1983,1984,1985,1986,1987,1988,1989,
#  1990,1991,1992,1993,1994,1995,1996,1997,1998,1999), dim=c(40,1))

#monthArray<-array(data=
#  c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),
#  dim=c(12,1))
monthArray<-array(data=
  c("01","02","03","04","05","06","07","08","09","10","11","12"),
  dim=c(12,1))

numYears<-length(yearsArray)
newlynMonthlyMean<-array(NA,dim=c(12,numYears))

for (jj in 1:numYears) {
# Test for leap year
  if ((yearsArray[jj] %% 4)==0) {
    if ((yearsArray[jj] %% 100)==0) {
      if ((yearsArray[jj] %% 400)==0) {
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

  zettMonthlyMean<-array(NA,dim=c(28420,12))
  zettDailyMean<-array(NA,dim=c(28420,29))

  outMonthlyMean <- file(paste("zett.mm.",yearsArray[jj],".dat",sep=""), "wb")

  for (ii in 1:nmons) {

    fname <- paste("S12run405ZetFiles/",inname[ii], substr(as.character(yearsArray[jj]),3,4),sep="")
    inConn <- file(fname, "rb")
    outDailyMean <- file(paste("dailyMeans/zett.dm.",monthArray[ii],".",   
      yearsArray[jj],".dat",sep=""), "wb")

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

    meanCount<-1
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

    zet<-array(0,dim=c(npsea,39))

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
      if ((jk <= 17 ) | (jk > (ntimes-16))) {
        timeLev <- 1
      } else {
        if (timeLev < 39) {
          timeLev <- timeLev+1
        } else {
# Filter the tides from the data
          zetMean<-array(0,dim=npsea)
          for (ip in 1:npsea) {
            zetMean[ip] <- sum(zet[ip,]*weightArray)/30
          }
# Write out mean to file
#          writeBin(as.vector(zetMean), outDailyMean, size=4)
          print(paste("Writing daily mean: ",meanCount," jk: ",jk
            ," Month: ",ii," Year: ",yearsArray[jj]))
# Place mean in daily mean array
          zettDailyMean[,meanCount] <- zetMean
          meanCount <- meanCount+1
# Shift the usable "old" data to the beginning of the array
          old <- array(0,dim=c(npsea,39))
          old[,c(1:15)]<-zet[,c(25:39)]
          zet <- old
          timeLev <- 16
        }
      }
      zet[,timeLev]<-zett

# Produce time series at Newlyn
      newlyn[jk]<-zet[17701,timeLev]

    }
#
# Doodson filter on hourly values
# See http://www.pol.ac.uk/psmslh/gloup/doodson_X0.html for details
# Running weighted mean of 19 values either side of a central value
# Weights are (1010010110201102112 0 2112011020110100101)/30
    filtNewlyn<-array(NA,dim=c(hours,1))
    meanArray<-array(NA,dim=c(39,1))
#    weightArray<-array(
#    c(1,0,1,0,0,1,0,1,1,0,2,0,1,1,0,2,1,1,2,0,2,1,1,2,0,1,1,0,2,0,1,1,0,1,0,0,1,0,1),dim=c(39,1))
    for (i in 1:(hours-38)) {
      meanArray<-(newlyn[i:(i+38),1]*weightArray)
      filtNewlyn[i+19]<-sum(meanArray)/30
    }
    dayMeanNewlyn<-array(NA,c((days-2),1))
    middays<-seq(36,hours,24)
    for (i in 1:(days-2)) {
      dayMeanNewlyn[i]<-filtNewlyn[middays[i]]
    }

    newlynMonthlyMean[ii,jj] <- sum(dayMeanNewlyn)/(days-2)
    for (jk in 1:npsea) {
      zettMonthlyMean[jk,ii] <- sum(zettDailyMean[jk,], na.rm=TRUE)/(days-2)
    }

    close(inConn)
    close(outDailyMean)

  }
# End of data reading

# Write monthly means out
  zetMeanOut<-array(0,dim=npsea)
  for (ii in 1:nmons) {
    writeBin(as.vector(zettMonthlyMean[,ii]),outMonthlyMean, size=4)
  }
  close(outMonthlyMean)
###########################
## Print out monthly maps #
###########################
#
#  for (ii in 1:12) {
#    if (nchar(as.character(ii))<2){
#      iistr<-paste("0",ii,sep="")
#    } else {
#      iistr<-paste(ii,sep="")
#    }
#    postscript(file=paste("/local/simonhdata/data/polcomms/ps/zett",
#      iistr,yearsArray[jj],".ps",sep=""))
#    zetMonthlyArray<-array(NA,dim=c(l,m))
#    for (ip in 1:npsea) {
#      zetMonthlyArray[isea[ip],jsea[ip]]<-zettMonthlyMean[ip,ii]
#    }
#    image.plot(lon,lat,zetMonthlyArray[c(1:198),c(1:224)],zlim=c(-1,1),
#      col=jet.colors(100))
#    dev.off()
#  }

# End of years loop
}
#
zetArray<-array(NA,dim=c(l,m))
for (ip in 1:npsea) {
  zetArray[isea[ip],jsea[ip]]<-zet[ip,1]
}
#image.plot(lon,lat,zetArray[c(1:198),c(1:224)],zlim=c(-4,4),
#  col=jet.colors(100))
#
#x11()
zetMonthlyArray<-array(NA,dim=c(l,m))
for (ip in 1:npsea) {
  zetMonthlyArray[isea[ip],jsea[ip]]<-zettMonthlyMean[ip,12]
}
#image.plot(lon,lat,zetMonthlyArray[c(1:198),c(1:224)],zlim=c(-1,1),
#  col=jet.colors(100))
#
#x11()
zetDailyArray<-array(NA,dim=c(l,m))
for (ip in 1:npsea) {
  zetDailyArray[isea[ip],jsea[ip]]<-zettDailyMean[ip,29]
}
#image.plot(lon,lat,zetDailyArray[c(1:198),c(1:224)],zlim=c(-1.5,1.5),
#  col=jet.colors(100))
#
#
# Create mean annual cycle for Newlyn
newlynMeanAnnualCycle <- array(NA,c(12,1))
for (ii in 1:12) {
  newlynMeanAnnualCycle[ii] <- mean(newlynMonthlyMean[ii,],na.rm=TRUE)
}
#plot(c(1:12),newlynMeanAnnualCycle,type="b")
write.table(newlynMeanAnnualCycle,"newlynAnnCycleFromPOLCOMS1960To1989.txt",quote=F,row.names=F,col.names=F)
