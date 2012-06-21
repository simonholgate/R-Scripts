# Script to read in the POLCOMMS S12 (12km) bathymetry and produce a mask of
# the <200m depth regions
junk<-read.table("S12.bathy",na.strings="0.00")
bathy<-array(NA,dim=c(198,224))
x<-1
y<-224
for (ii in 1:2464) {
  for (jj in 1:18) {
    bathy[x,y]<-junk[ii,jj]
    x<-x+1
    if (x>198) {
      x<-1
      y<-y-1
    }
  }
}
deep<-which(bathy>200)
bathyShallow<-bathy
bathyShallow[deep]<-NA
bathyMask<-bathyShallow
mask<-which(!is.na(bathyShallow))
bathyMask[mask]<-1
library(fields)
image(bathyMask)

