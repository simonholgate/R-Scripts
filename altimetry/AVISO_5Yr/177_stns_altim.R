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

library(spam)

# Load matrixMethods for use later
source("~/workspace/RScripts/matrixMethods.R")

# High res grid
lon<-seq(from=1/3, length=1080, by=1/3)
junk<-read.table('~/diskx/altimetry/lats.lis')
lat<-junk$V1

# Get lons and lats for all the 177 stations
stntable <-
  read.table('~/diskx/177StationsUpdate2009/5YrMeans/includedTideGauges.mod.txt',
  col.names=c("ccode","scode","sname","lon","lat","nyrs"), 
  colClasses=c("character","character","character","numeric","numeric","integer"),
  sep=";")

message('Read stntable')

nstns <- length(stntable$lon)
x <- cbind(stntable$lon,stntable$lat)

loc <- array(NA, dim=c(nstns,2))

for (i in 1:nstns){
# We need to define an area around each station before  doing the expansion
# otherwise the grid is too big. We'll make this +/- 2 degs
  lonrange <- intersect(which(lon>=(x[i,1]-2)),which(lon<=(x[i,1]+2)))
  latrange <- intersect(which(lat>=(x[i,2]-2)),which(lat<=(x[i,2]+2)))

  y <- expand.grid(lon[lonrange],lat[latrange])
  near <- nearest.dist(rbind(x[i,],x[i,]),y)
  
  gtzero <- which(near[1,]>0)
# There are some station points outside of the altimetry region so these will
# return all zeros
  if (length(gtzero>0)){
    nearest <- gtzero[which.min(near[1,gtzero])]
  } else {
    nearest <- NA
  }
# Now we need to refer this point back to the larger grid  
# - how many of these will end up on land?
  stnloc <- y[nearest,]
  loc[i,1] <- which(lon==as.numeric(stnloc[1]))
  loc[i,2] <- which(lat==as.numeric(stnloc[2]))
}
message('Calculated nearest')

# jd (Julian Day) in this context appears to be days since 1950 
daytable<-read.table('~/diskx/altimetry/aviso/ref_h_anom/daytable.txt',
  col.names=c("seq","jd","year","month","day","yearday"), 
  colClasses=c("integer"))

xres<-1/3
l<-1080
m<-915
lm<-l*m
n<-818

time<-seq(as.Date("1992/10/14"), by="7 days", length.out=n)

dataArray<-array(NA,dim=c(n,nstns))

# Read the filtered time series that CWH has produced and see which grid
# points have only null values (=9999.0). These can then be flagged as
# land/ice points and removed from the time series.

inConn <- 
  file("~/diskx/altimetry/aviso/ref_h_anom/ref_h_anom_stackmaps.dat", "rb")

for (i in 1:n) {
  message(paste('Reading data for slice: ', as.character(i)))
  
  data<-readBin(inConn, what="numeric", n = lm, size = 4, endian='little')
  dim(data)<-c(l,m)
# Extract stn data at locs
  for(j in 1:nstns){
    px <- loc[j,1]
    py <- loc[j,2]
    dataArray[i,j] <- data[px, py]
# On first timeslice only, if data isn't available at that point, 
# examine surrounding points
    if(i==1){
      if(dataArray[i,j]==9999.0) {

        surround <- c(px+1, py, px-1, py, px, py+1, px, py-1, 
          px+1, py+1, px-1, py+1, px+1, py-1, px-1, py-1)
        dim(surround) <- c(2,8)
        
        junk <- data[t(surround)]
        k <- min(which(junk!=9999.0))  # k is infinite if all are 9999 
          if(is.finite(k)){
            dataArray[i,j] <- data[surround[1,k], surround[2,k]]
# Update that location for next time
            loc[j,1] <- surround[1,k]
            loc[j,2] <- surround[2,k] 
          } else {
# There are 5 locations where the surrounding squares are all 9999 so expand the
# examined region by another square for these:
#> which(is.na(trends[1,]))
#[1]  12  66 139 140 156
#> stntable$sname[ which(is.na(trends[10,]))]
#[1] "OSLO"          "DUBLIN"        "BUENOS AIRES"  "PALERMO"      
#[5] "WASHINGTON DC"
            surround <- c(
              px+2, py, px+2, py+1, px+2, py-1, px+2, py+2, px+2, py-2, 
              px-2, py, px-2, py+1, px-2, py-1, px-2, py+2, px-2, py-2,
              px+1, py+2, px, py+2, px-1, py+2, 
              px+1, py-2, px, py-2, px-1, py-2
              )
            dim(surround) <- c(2,16)
            junk <- data[t(surround)]
            k <- min(which(junk!=9999.0))  # k is infinite if all are 9999 
            if(is.finite(k)){
              dataArray[i,j] <- data[surround[1,k], surround[2,k]]
# Update that location for next time
              loc[j,1] <- surround[1,k]
              loc[j,2] <- surround[2,k] 
            } else {
            # Do nothing, leave location unchanged
            }           

          }  
        
      }
    }
  }
}
close(inConn)
message('Closed inConn')
# Replace 9999.0 by NA
junk <- which(dataArray==9999.0)
dataArray[junk] <- NA

# Convert data from cm to mm
dataArray <- dataArray*10

#stop("Imported data")

# Remove annual and semi annual signals
w1 <- 2*pi/(365)
w2 <- 2*pi/(365/2)
jd <- daytable$jd[1:n]
residArray <- array(NA, dim=c(n,nstns))

for (i in 1:nstns){
  if(length(which(is.na(dataArray[,i])))==818){
# Do nothing and leave residArray as all NAs
  } else {  
    mdata <- data.frame(dataArray=dataArray[,i], t=jd,   
      w1cos=cos(w1*jd), w1sin=sin(w1*jd), 
      w2cos=cos(w2*jd), w2sin=sin(w2*jd))

    model <- lm(dataArray ~  w1cos + w1sin + w2cos + w2sin, data=mdata, x=TRUE)
    indx <- as.integer(rownames(model$x))
    residArray[indx,i] <- model$resid
  }
}
message('Removed Seasonal signals')
# Calculate annual means from the residuals for 1993 to 2007
nyears<-15
annMeans <- array(NA, dim=c(nyears,nstns))
annJD <- vector(mode="numeric", length=nyears)
for(i in 1:nyears){
  junk <- which(daytable$year==(i+1992))
  annMeans[i,] <- colMeans(residArray[junk,], na.rm=T)
  annJD[i] <- mean(daytable$jd[junk], na.rm=T)
}
message('Calculated means')
# Calculate 5 year running mean trends
trendMidPoints <- c(1996:2006)
trends <- array(NA,dim=c(11,nstns))
for(i in 1:11){
  for(j in 1:nstns){
    if(length(which(is.na(annMeans[,j])))==nyears){
# Do nothing and leave residArray as all NAs
    } else {
      junk <- lm(annMeans[i:(i+4),j] ~ c(1:5))
      trends[i,j] <- junk$coef[2]
    }
  }
}
message('Calculated trends')
# Have a look at tg locs - add an area surrounding it so that it's clearer
library(fields)
mask <- data
mask[which(mask<9999.0)] <- NA
for(j in 1:nstns){
  px <- loc[j,1]
  py <- loc[j,2]
  surround <- c(px+1, py, px-1, py, px, py+1, px, py-1, 
    px+1, py+1, px-1, py+1, px+1, py-1, px-1, py-1)
  dim(surround) <- c(2,8)
        
  mask[t(surround)] <- 9990
}
image.plot(lon,lat,mask,col=rainbow(5),nlevel=2)
message('Plotted locations')