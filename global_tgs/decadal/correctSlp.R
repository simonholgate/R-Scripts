library(fields)
source("~/Work/SeaLevelRise/RScripts/jet.colors.R")

nlon<-72
nlat<-37
nyr<-128
xlon<-seq(from=-180,by=5,length=nlon)
ylat<-seq(from=-90,by=5,length=nlat)

#load("slpAnnualArray.RData")
slpAnnualArray<-array(NA,dim=c(nlon,nlat,nyr))
con<-file("../pressure/hadSLP1_1871-1998.ann.bin","rb")

for (i in 1:nyr){
  slp<-readBin(con,what="integer", size=4, n=(nlat*nlon))

  slp<-array(slp,dim=c(nlon,nlat))
  slpAnnualArray[,,i]<-slp
}
close(con)
# Convert from Pa to mb
slpAnnualArray<-slpAnnualArray*0.01
# In converting to mm from mb remember that 1mb = 9.948mm of sea surface
# height
slpAnnualArray<-slpAnnualArray*9.948

stns177Lon<-stns177[,4]
stns177Lat<-stns177[,5]
j<-which(stns177Lon<180)

# Remember that the pressure longitudes go from 180 to 535.
# This means that 0 -> 360 and 180 -> 540

stns177Lon535 <- stns177Lon
stns177Lon535[j] <- stns177Lon535[j]+360

i1 <- as.integer(stns177Lon535/5 - 35)
i2 <- i1+1

j <- which(i2 == 73)
i2[j] <- 1

j1 <- as.integer(((90 - stns177Lat)/5) + 1)
j2 <- j1+1

# Error checking
i <- which(i1<=0)
j <- which(i1>72)
if ((length(i) > 0) | (length(j) > 0)){
  stop(paste("Error in i1:", as.character(i1[union(i,j)])))
}

i <- which(i2<=0)
j <- which(i2>72)
if ((length(i) > 0) | (length(j) > 0)){
  stop(paste("Error in i2:", as.character(i2[union(i,j)])))
}

i <- which(j1<=0)
j <- which(j1>37)
if ((length(i) > 0) | (length(j) > 0)){
  stop(paste("Error in j1:", as.character(j1[union(i,j)])))
}

i <- which(j2<=0)
j <- which(j2>37)
if ((length(i) > 0) | (length(j) > 0)){
  stop(paste("Error in j2:", as.character(j2[union(i,j)])))
}

# corners are top left, top right, bot right, bot left
# Now, because of the fact that we subtracted 35 instead of adding 1 above
# (due to the weird co-ordinate scheme of the Hadley data set) we need to add
# 35 instead of subtracting 1 here
x1 <- stns177Lon535 - as.numeric(i1+35)*5.
x2 <- 5. - x1
y1 <- 90.0 - as.numeric(j1-1)*5. - stns177Lat
y2 <- 5. - y1
#if(x1.lt.0.0.or.x2.lt.0.0.or.y1.lt.-90.0.or.y2.lt.-90.0) then
#  write(6,*) 'error x1,x2,y1,y2,blon,blat',x1,x2,y1,y2,blon,blat
#       stop
#endif
# Error checking
i <- which(x1<0)
if (length(i) > 0) {
  stop(paste("Error in x1:", as.character(x1[i])))
}

i <- which(x2<0)
if (length(i) > 0){
  stop(paste("Error in x2:", as.character(x2[i])))
}

i <- which(y1< -90)
if (length(i) > 0) {
  stop(paste("Error in y1:", as.character(y1[i])))
}

i <- which(y2< -90)
if (length(i) > 0) {
  stop(paste("Error in y2:", as.character(y2[i])))
}

wt <- array(NA, dim=c(4,177))
wt[1,] <- x2*y2
wt[2,] <- x1*y2
wt[3,] <- x1*y1
wt[4,] <- x2*y1

wttot <- vector(mode="numeric", length=177)
for (i in 1:177) {
  wttot[i] <- sum(wt[,i])
  for (j in 1:4) {
    wt[j,i] <- wt[j,i]/wttot[i]
  }
}

# Create time series of pressures for the 177 stations
ssh177Stns <- array(0,dim=c(nyr,177))
#for (i in 1:177) {
#  ssh177Stns[,i] <- slpAnnualArray[xcoord[i],ycoord[i],]
#}
ap <- array(data=NA, dim=c(nyr,177))

for (icorn in 1:4) {
  if(icorn==1) i<-i1
  if(icorn==1) j<-j1
  if(icorn==2) i<-i2
  if(icorn==2) j<-j1
  if(icorn==3) i<-i2
  if(icorn==3) j<-j2
  if(icorn==4) i<-i1
  if(icorn==4) j<-j2

  for (count in 1:177){
    ap[,count] <- slpAnnualArray[i[count],j[count],]
    ssh177Stns[,count] <- ssh177Stns[,count] + ap[,count]*wt[icorn,count]
  }
}


#j<-which(stns177Lon>180)

#stns177LonNeg <- stns177Lon
#stns177LonNeg[j]<- -1*(360-stns177Lon[j])

#xcoord <- ceiling(((stns177LonNeg%/%5+round((stns177LonNeg%%5)/5))*5 
#    - xlon[1])/5) + 1
#ycoord <- ceiling(((stns177Lat%/%5+round((stns177Lat%%5)/5))*5
#    - ylat[1])/5) + 1

# Calculate rates of ssh change
s177SlpAnnualYrs<-c(1871:1998)
s177SlpRates <- array(dim = c(46,177))
s177SlpDecadeMidPoints <- array(dim = c(46,1))
for (i in 1:177){
  for (j in 1948:1989) {
    interval <- findInterval(j:(j+9),s177SlpAnnualYrs)
    fit <- lsfit(s177SlpAnnualYrs[interval], ssh177Stns[interval, i])
    coeffs <- coef(fit)
    s177SlpRates[(j-1947),i] <- coeffs[2]
    s177SlpDecadeMidPoints[j-1947] <- j+4
  }
}

# Save the data
#save(s177SlpDecadeMidPoints,s177SlpRates,file="stns177SlpRates.RData")
