# Correction for inverted barometer based of Phil's Jabber pro-fortran code
library(fields)

nlon<-72
nlat<-37
nyr<-158
nmon<-nyr*12

xlon<-seq(from=-180,by=5,length=nlon)
ylat<-seq(from=-90,by=5,length=nlat)

#load("slpMonthlyArray.RData")
slpMonthlyArray<-array(NA,dim=c(nlon,nlat,nmon))
con<-file("~/diskx/HadSLP2r/hadSLP2r.monthly.bin","rb")

for (i in 1:nmon){
  slp<-readBin(con,what="integer", size=4, n=(nlat*nlon))

  slp<-array(slp,dim=c(nlon,nlat))
  slpMonthlyArray[,,i]<-slp
}
close(con)
# Convert from Pa to mb
slpMonthlyArray<-slpMonthlyArray*0.01
# In converting to mm from mb remember that 1mb = 9.948mm of sea surface
# height
slpMonthlyArray<-slpMonthlyArray*9.948

stnsLon<-stns[,4]
stnsLat<-stns[,5]
j<-which(stnsLon<180)

# Remember that the pressure longitudes go from 180 to 535.
# This means that 0 -> 360 and 180 -> 540

stnsLon535 <- stnsLon
stnsLon535[j] <- stnsLon535[j]+360

i1 <- as.integer(stnsLon535/5 - 35)
i2 <- i1+1

j <- which(i2 == 73)
i2[j] <- 1

j1 <- as.integer(((90 - stnsLat)/5) + 1)
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
x1 <- stnsLon535 - as.numeric(i1+35)*5.
x2 <- 5. - x1
y1 <- 90.0 - as.numeric(j1-1)*5. - stnsLat
y2 <- 5. - y1

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

wt <- array(NA, dim=c(4,numStns))
wt[1,] <- x2*y2
wt[2,] <- x1*y2
wt[3,] <- x1*y1
wt[4,] <- x2*y1

wttot <- vector(mode="numeric", length=numStns)
for (i in 1:numStns) {
  wttot[i] <- sum(wt[,i])
  for (j in 1:4) {
    wt[j,i] <- wt[j,i]/wttot[i]
  }
}

# Create time series of pressures for the stations
sshStns <- array(0,dim=c(nyr,numStns))

ap <- array(data=NA, dim=c(nyr,numStns))

for (icorn in 1:4) {
  if(icorn==1) i<-i1
  if(icorn==1) j<-j1
  if(icorn==2) i<-i2
  if(icorn==2) j<-j1
  if(icorn==3) i<-i2
  if(icorn==3) j<-j2
  if(icorn==4) i<-i1
  if(icorn==4) j<-j2

  for (count in 1:numStns){
    ap[,count] <- slpMonthlyArray[i[count],j[count],]
    sshStns[,count] <- sshStns[,count] + ap[,count]*wt[icorn,count]
  }
}


# Calculate rates of ssh change
slpMonthlyYrs<-seq(from=(1850+1/24),to=(2007+23/24), by=1/12)
slpRates <- array(dim = c(lenRatesP1,numStns))
slpMidPoints <- c(startRate:endRateP1)
for (i in 1:numStns){
  for (j in startYr:lastRate) {
    interval <- findInterval(j:(j+rateLenM1),slpAnnualYrs)
    fit <- lsfit(slpAnnualYrs[interval], sshStns[interval, i])
    coeffs <- coef(fit)
    slpRates[(j-startYrM1),i] <- coeffs[2]
  }
}

# Save the data
save(slpMidPoints,slpRates,file="stnsSlpRates.RData")
