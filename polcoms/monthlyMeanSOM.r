###############################################################################
# Reads the monthly mean sea level from the POLCOMMS model output and
# produces self organising maps (SOM) from which to explore spatial patterns
# in the data. 
# This requires reading 120 arrays (10 years x 12 months) 
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
# Number of months to iterate over this will be determined by the number of
# file names in inname
inname<-c( 
"zet_UBVB.S12run401.jan",
"zet_UBVB.S12run401.feb",
"zet_UBVB.S12run401.mar",
"zet_UBVB.S12run401.apr",
"zet_UBVB.S12run401.may",
"zet_UBVB.S12run401.jun",
"zet_UBVB.S12run401.jul",
"zet_UBVB.S12run401.aug",
"zet_UBVB.S12run401.sep",
"zet_UBVB.S12run401.oct",
"zet_UBVB.S12run401.nov",
"zet_UBVB.S12run401.dec" 
)

yearsArray<-array(data=
  c(1960:1989), dim=c(30,1))

numYears<-length(yearsArray)

####################
# Read header info #
####################
isea<-array(NA,dim=npsea)
jsea<-array(NA,dim=npsea)
inNpseaIseaJsea <- file("npseaiseajsea.dat", "rb")
npsea<-readBin(inNpseaIseaJsea, what="integer", n = 1, size = NA)
isea<-readBin(inNpseaIseaJsea, what="integer", n = npsea, size = NA)
jsea<-readBin(inNpseaIseaJsea, what="integer", n = npsea, size = NA)
close(inNpseaIseaJsea)

nmons <- length(inname)*numYears

# Note that the orientation of the array is different for the SOM!!
zettMonthlyMeans <- array(NA,dim=c(nmons,npsea))
meanCount <- 1
newlynMonthlyMean<-array(NA,dim=c(12,numYears))

for (jj in 1:numYears) {

  inMonthlyMean <- file(paste("zett.mm.",yearsArray[jj],".dat",sep=""), "rb")

  for (ii in 1:12) {
    zettMonthlyMeans[meanCount,] <- 
      readBin(inMonthlyMean, what="numeric", n = npsea, size=4)
    newlynMonthlyMean[ii,jj] <- zettMonthlyMeans[meanCount,17701]

    meanCount <- meanCount + 1
  }
  close(inMonthlyMean)

# End of years loop
}
#
# Load in the 200m bathymetry mask
load("bathy/bathMask.R")
# Create vectorised mask
bathyMaskVector <- array(NA,dim=c(npsea,1))
for (ip in 1:npsea) {
  bathyMaskVector[ip] <- bathyMask[isea[ip],jsea[ip]]
}
# Mask the monthly array data and remove mean annual cycle from data
longTermMonthlyMeanArray<-array(NA,dim=c(12,npsea))
zettMonthlyMeanArray<-array(t(zettMonthlyMeans),dim=c(npsea,12,30))
for (ip in 1:npsea){
  for (ii in 1:12) {
    longTermMonthlyMeanArray[ii,ip] <- mean(zettMonthlyMeanArray[ip,ii,])
  }
}

maskedMonthlyArray<-array(NA,dim=c(nmons,npsea))
DSMonthlyArray<-array(NA,dim=c(nmons,npsea))
DSMaskedMonthlyArray<-array(NA,dim=c(nmons,npsea))
monthlyIndex <- rep(1:12,30)

for (ii in 1:nmons) {
  DSMonthlyArray[ii,]<-zettMonthlyMeans[ii,]-  
    longTermMonthlyMeanArray[monthlyIndex[ii],]
  maskedMonthlyArray[ii,]<-zettMonthlyMeans[ii,]*bathyMaskVector
  DSMaskedMonthlyArray[ii,]<-DSMonthlyArray[ii,]*bathyMaskVector
}
#
zetMonthlyArray<-array(NA,dim=c(l,m))
for (ip in 1:npsea) {
  zetMonthlyArray[isea[ip],jsea[ip]]<-zettMonthlyMeans[360,ip]
}
image.plot(lon,lat,zetMonthlyArray[c(1:198),c(1:224)],zlim=c(-1,1),
  col=jet.colors(100))
#
# Create mean annual cycle for Newlyn
newlynMeanAnnualCycle <- array(NA,c(12,1))
for (ii in 1:12) {
  newlynMeanAnnualCycle[ii] <- mean(newlynMonthlyMean[ii,],na.rm=TRUE)
}
plot(c(1:12),newlynMeanAnnualCycle,type="b")
#
# Delete masked area for SOM
j<-which(is.finite(maskedMonthlyArray[1,]))
maskedMonthlyArray<-maskedMonthlyArray[,j]
DSMaskedMonthlyArray<-DSMaskedMonthlyArray[,j]
