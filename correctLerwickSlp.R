library(fields)

nstns<-6
nlon<-72
nlat<-37
nyr<-158
nmon<-12
xlon<-seq(from=-180,by=5,length=nlon)
ylat<-seq(from=-90,by=5,length=nlat)

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

lerwickStnsLon<-lerwickStnsLatLon[2,]
lerwickStnsLat<-lerwickStnsLatLon[1,]
j<-which(lerwickStnsLon<180)

# Remember that the pressure longitudes go from 180 to 535.
# This means that 0 -> 360 and 180 -> 540

lerwickStnsLon535 <- lerwickStnsLon
lerwickStnsLon535[j] <- lerwickStnsLon535[j]+360

i1 <- as.integer(lerwickStnsLon535/5 - 35)
i2 <- i1+1

j <- which(i2 == 73)
i2[j] <- 1

j1 <- as.integer(((90 - lerwickStnsLat)/5) + 1)
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
x1 <- lerwickStnsLon535 - as.numeric(i1+35)*5.
x2 <- 5. - x1
y1 <- 90.0 - as.numeric(j1-1)*5. - lerwickStnsLat
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
sshLerwickStns <- array(0,dim=c(nyr,nstns))
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
    sshLerwickStns[,count] <- sshLerwickStns[,count] + ap[,count]*wt[icorn,count]
  }
}

# Calculate rates of ssh change over the period of each tide gauge:
# Lerwick 1957-2007
# Wick 1965-2007
# Invergordon 1960-1971
# Stornoway 1977-2007
# Aberdeen 1957-2007
# Torshavn 1957-2006

stnsAnnualYrs<-seq.Date(from=as.Date("1850/7/1"),to=as.Date("2007/7/1"), by="1 year")
#nstnyrs <- length(lerwickStnsSlpAnnualYrs)
lerwickStnAnnualYrs <- seq.Date(from=as.Date("1957/7/1"),to=as.Date("2007/7/1"), by="1 year")
wickStnAnnualYrs <- seq.Date(from=as.Date("1965/7/1"),to=as.Date("2007/7/1"), by="1 year")
invergordonStnAnnualYrs <- seq.Date(from=as.Date("1960/7/1"),to=as.Date("1971/7/1"), by="1 year")
stornowayStnAnnualYrs <- seq.Date(from=as.Date("1977/7/1"),to=as.Date("2007/7/1"), by="1 year")
aberdeenStnAnnualYrs <- seq.Date(from=as.Date("1957/7/1"),to=as.Date("2007/7/1"), by="1 year")
torshavnStnAnnualYrs <- seq.Date(from=as.Date("1957/7/1"),to=as.Date("2006/7/1"), by="1 year")

stnsSlpRates <- vector(length=nstns, mode="numeric")

interval <- c(which(stnsAnnualYrs==lerwickStnAnnualYrs[1]):
	which(stnsAnnualYrs==lerwickStnAnnualYrs[length(lerwickStnAnnualYrs)]))
fit <- lm(sshLerwickStns[interval, 1] ~ seq(1:length(interval)))
coeffs <- coef(fit)
stnsSlpRates[1] <- coeffs[2]

interval <- c(which(stnsAnnualYrs==wickStnAnnualYrs[1]):
	lwhich(stnsAnnualYrs==wickStnAnnualYrs[ength(wickStnAnnualYrs)]))
fit <- lm(sshLerwickStns[interval, 2] ~ seq(1:length(interval)))
coeffs <- coef(fit)
stnsSlpRates[2] <- coeffs[2]

interval <- c(which(stnsAnnualYrs==invergordonStnAnnualYrs[1]):
	which(stnsAnnualYrs==invergordonStnAnnualYrs[length(invergordonStnAnnualYrs)]))
fit <- lm(sshLerwickStns[interval, 3] ~ seq(1:length(interval)))
coeffs <- coef(fit)
stnsSlpRates[3] <- coeffs[2]

interval <- c(which(stnsAnnualYrs==stornowayStnAnnualYrs[1]):
	which(stnsAnnualYrs==stornowayStnAnnualYrs[length(stornowayStnAnnualYrs)]))
fit <- lm(sshLerwickStns[interval, 4] ~ seq(1:length(interval)))
coeffs <- coef(fit)
stnsSlpRates[4] <- coeffs[2]

interval <- c(which(stnsAnnualYrs==aberdeenStnAnnualYrs[1]):
	which(stnsAnnualYrs==aberdeenStnAnnualYrs[length(aberdeenStnAnnualYrs)]))
fit <- lm(sshLerwickStns[interval, 5] ~ seq(1:length(interval)))
coeffs <- coef(fit)
stnsSlpRates[5] <- coeffs[2]

interval <- c(which(stnsAnnualYrs==torshavnStnAnnualYrs[1]):
	which(stnsAnnualYrs==torshavnStnAnnualYrs[length(torshavnStnAnnualYrs)]))
fit <- lm(sshLerwickStns[interval, 6] ~ seq(1:length(interval)))
coeffs <- coef(fit)
stnsSlpRates[6] <- coeffs[2]

#lerwickStnsSlpDecadeMidPoints <- seq.Date(as.Date("1960/1/1"),as.Date("2005/1/1"), by="5 years")
# Save the data
save(stnsSlpRates,file="lerwickStnsSlpRates.RData")
