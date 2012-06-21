library(fields)

nstns<-2
nlon<-72
nlat<-37
nyr<-158
nmon<-12
xlon<-seq(from=-180,by=5,length=nlon)
ylat<-seq(from=-90,by=5,length=nlat)

newlynStnsLatLon <- c(50.1, -5.55, 48.3833, 4.5)
dim(newlynStnsLatLon) <- c(2,2)

#load("slpAnnualArray.RData")
#slpMonthlyArray<-array(NA,dim=c(nlon,nlat,nmon,nyr))
slpAnnualArray<-array(NA,dim=c(nlon,nlat,nyr))

con<-file("~/diskx/HadSLP2r/hadSLP2_kij_1850-2007.bin","rb")

for (i in 1:nlat){
  for (j in 1:nlon){
    for (k in 1:nyr){
      slp<-readBin(con,what="numeric", size=4, n=nmon)
#      slpMonthlyArray[j,i,,k]<-slp
      slpAnnualArray[j,i,k]<-mean(slp, na.rm=TRUE)
    }
  }
}
#for (i in 1:nmon){
#  slp<-readBin(con,what="numeric", size=4, n=(nlat*nlon))
#
#  slp<-array(slp,dim=c(nlon,nlat))
#  slpAnnualArray[,,i]<-slp
#}
close(con)

#dim(slpMonthlyArray)<-c(nlon,nlat,nmon*nyr)

# Convert from Pa to mb
slpAnnualArray<-slpAnnualArray*0.01
# In converting to mm from mb remember that 1mb = 9.948mm of sea surface
# height
slpAnnualArray<-slpAnnualArray*9.948

newlynStnsLon<-newlynStnsLatLon[2,]
newlynStnsLat<-newlynStnsLatLon[1,]
j<-which(newlynStnsLon<180)

# Remember that the pressure longitudes go from 180 to 535.
# This means that 0 -> 360 and 180 -> 540

newlynStnsLon535 <- newlynStnsLon
newlynStnsLon535[j] <- newlynStnsLon535[j]+360

i1 <- as.integer(newlynStnsLon535/5 - 35)
i2 <- i1+1

j <- which(i2 == 73)
i2[j] <- 1

j1 <- as.integer(((90 - newlynStnsLat)/5) + 1)
j2 <- j1+1

# Error checking
i <- which(i1<=0)
j <- which(i1>nlon)
if ((length(i) > 0) | (length(j) > 0)){
  stop(paste("Error in i1:", as.character(i1[union(i,j)])))
}

i <- which(i2<=0)
j <- which(i2>nlon)
if ((length(i) > 0) | (length(j) > 0)){
  stop(paste("Error in i2:", as.character(i2[union(i,j)])))
}

i <- which(j1<=0)
j <- which(j1>nlat)
if ((length(i) > 0) | (length(j) > 0)){
  stop(paste("Error in j1:", as.character(j1[union(i,j)])))
}

i <- which(j2<=0)
j <- which(j2>nlat)
if ((length(i) > 0) | (length(j) > 0)){
  stop(paste("Error in j2:", as.character(j2[union(i,j)])))
}

# corners are top left, top right, bot right, bot left
# Now, because of the fact that we subtracted 35 instead of adding 1 above
# (due to the weird co-ordinate scheme of the Hadley data set) we need to add
# 35 instead of subtracting 1 here
x1 <- newlynStnsLon535 - as.numeric(i1+35)*5.
x2 <- 5. - x1
y1 <- 90.0 - as.numeric(j1-1)*5. - newlynStnsLat
y2 <- 5. - y1
#
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

wt <- array(NA, dim=c(4,nstns))
wt[1,] <- x2*y2
wt[2,] <- x1*y2
wt[3,] <- x1*y1
wt[4,] <- x2*y1

wttot <- vector(mode="numeric", length=nstns)
for (i in 1:nstns) {
  wttot[i] <- sum(wt[,i])
  for (j in 1:4) {
    wt[j,i] <- wt[j,i]/wttot[i]
  }
}

# Create time series of pressures for the nstns stations
sshNewlynStns <- array(0,dim=c(nyr,nstns))
#
ap <- array(data=NA, dim=c(nyr,nstns))

for (icorn in 1:4) {
  if(icorn==1) i<-i1
  if(icorn==1) j<-j1
  if(icorn==2) i<-i2
  if(icorn==2) j<-j1
  if(icorn==3) i<-i2
  if(icorn==3) j<-j2
  if(icorn==4) i<-i1
  if(icorn==4) j<-j2

  for (count in 1:nstns){
    ap[,count] <- slpAnnualArray[i[count],j[count],]
    sshNewlynStns[,count] <- sshNewlynStns[,count] + ap[,count]*wt[icorn,count]
  }
}

# Calculate rates of ssh change over the period of each tide gauge:
# Newlyn 1914-2006
# Brest 1914-2006

stnsAnnualYrs<-seq.Date(from=as.Date("1850/7/1"),to=as.Date("2007/7/1"), by="1 year")
#nstnyrs <- length(newlynStnsSlpAnnualYrs)
stnAnnualYrs <- seq.Date(from=as.Date("1914/7/1"),to=as.Date("2006/7/1"), by="1 year")
stnAnnualYrs1996 <- seq.Date(from=as.Date("1996/7/1"),to=as.Date("2006/7/1"), by="1 year")

stnsSlpRates <- vector(length=nstns, mode="numeric")
stnsSlpRates1996 <- vector(length=nstns, mode="numeric")

interval <- c(which(stnsAnnualYrs==stnAnnualYrs[1]):
	which(stnsAnnualYrs==stnAnnualYrs[length(stnAnnualYrs)]))
fit <- lm(sshNewlynStns[interval, 1] ~ seq(1:length(interval)))
coeffs <- coef(fit)
stnsSlpRates[1] <- coeffs[2]

fit <- lm(sshNewlynStns[interval, 2] ~ seq(1:length(interval)))
coeffs <- coef(fit)
stnsSlpRates[2] <- coeffs[2]

# Altimetry period

interval <- c(which(stnsAnnualYrs==stnAnnualYrs1996[1]):
	which(stnsAnnualYrs==stnAnnualYrs[length(stnAnnualYrs)]))
fit <- lm(sshNewlynStns[interval, 1] ~ seq(1:length(interval)))
coeffs <- coef(fit)
stnsSlpRates1996[1] <- coeffs[2]

fit <- lm(sshNewlynStns[interval, 2] ~ seq(1:length(interval)))
coeffs <- coef(fit)
stnsSlpRates1996[2] <- coeffs[2]
#newlynStnsSlpDecadeMidPoints <- seq.Date(as.Date("1960/1/1"),as.Date("2005/1/1"), by="5 years")
# Save the data
save(stnsSlpRates, stnsSlpRates1996, file="newlynStnsSlpRates.RData")
