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
# Brest 48 23 N  04 30 W
# Brest is at x=92, y=74
#> (-4.5-min)/xres
#[1] 91.99998
#> (48.38-ymin)/yres
#[1] 74.42001

library(fields)
# Use jet colors for images
source('./jet.colors.R')
lat<-seq(from=40.11111,by=1/9,length=224)
lon<-seq(from=-19.83333,by=1/6,length=198)

l<-198
m<-224
lm<-l*m
n<-20
#
zet<-array(0,dim=c(l,m))
# We're not interested in u or v so just use one array
#uvb<-array(0,dim=c(l,m))

# Number of months to iterate over this will be determined by the number of
# file names in inname
inname<-c( 
"zet_UBVB.S12run401.jan60"
#"zet_UBVB.S12run401.feb60",
#"zet_UBVB.S12run401.mar60",
#"zet_UBVB.S12run401.apr60",
#"zet_UBVB.S12run401.may60",
#"zet_UBVB.S12run401.jun60",
#"zet_UBVB.S12run401.jul60",
#"zet_UBVB.S12run401.aug60",
#"zet_UBVB.S12run401.sep60",
#"zet_UBVB.S12run401.oct60",
#"zet_UBVB.S12run401.nov60",
#"zet_UBVB.S12run401.dec60" 
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
monthArray<-array(data=
  c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"),
  dim=c(12,1))
newlynMonthlyMean<-array(NA,dim=c(12,1))
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
}

nmons<-length(inname)

for (ii in 1:nmons) {

  print('Initializing...')
  inConn <- file(inname[ii], "rb")

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

# Array for som analysis
# Might take some time this
  matts<-array(NA,dim=c(ntimes,npsea))

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

  for (jk in 1:ntimes) {
  print(paste('Reading',jk,'of',ntimes))
# Read junk first integer that Fortran writes
    junk1<-readBin(inConn, what="integer", n = 1, size = NA, endian = "big")
# Read data
    itimt<-readBin(inConn, what="integer", n = 1, size = NA, endian = "big")
# Read junk last integer that Fortran writes
    junk1<-readBin(inConn, what="integer", n = 1, size = NA, endian = "big")

# Read junk first integer that Fortran writes
    junk1<-readBin(inConn, what="integer", n = 1, size = NA, endian = "big")
# Read data
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

    for (ip in 1:npsea) {
      zet[isea[ip],jsea[ip]]<-zett[ip]
    }

# Add zett to the ts matrix matts
    matts[jk,] <- zett

# Produce time series at Newlyn
    newlyn[jk]<-zet[86,90]
#    print(paste(jk,zet[86,90],newlyn[jk]))

# 
    zet[which(zet==0)]<-NA

    if (nchar(as.character(jk))<2){
      jkstr<-paste("00",jk,sep="")
    } else {
      if (nchar(as.character(jk))<3){
      jkstr<-paste("0",jk,sep="")
      } else {
      jkstr<-paste(jk,sep="")
      }
    }
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
  meanArray<-array(NA,dim=c(39,1))
  weightArray<-array(
  c(1,0,1,0,0,1,0,1,1,0,2,0,1,1,0,2,1,1,2,0,2,1,1,2,0,1,1,0,2,0,1,1,0,1,0,0,1,0,1),dim=c(39,1))
  for (i in 1:(hours-38)) {
    meanArray<-(newlyn[i:(i+38),1]*weightArray)
    filtNewlyn[i+19]<-sum(meanArray)/30
  }
  dayMeanNewlyn<-array(NA,c((days-2),1))
  middays<-seq(36,hours,24)
  for (i in 1:(days-2)) {
    dayMeanNewlyn[i]<-filtNewlyn[middays[i]]
  }

  newlynMonthlyMean[ii] <- sum(dayMeanNewlyn)/(days-2)

  close(inConn)

}
# End of data reading
zet[which(zet==0)]<-NA
image.plot(lon,lat,zet[c(1:198),c(1:224)],zlim=c(-4,4),
  col=jet.colors(100))

