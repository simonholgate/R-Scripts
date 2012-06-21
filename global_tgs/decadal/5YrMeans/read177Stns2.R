# Script to read in the 177 stations included in the Holgate and Woodworth
# (2004) analysis.
# V2 is updated to use the HadSLP2 pressure data
#
stns177 <- read.table(file="includedTideGauges.mod.txt", sep=";",
col.names=c("ccode", "scode", "sname", "slondec", "slatdec", "numYrs"),
colClasses=c("character","character","character","numeric","numeric","integer"))
#
dataList177 <- list()
dataArray <- array(NA,dim=c(length(stns177[,1]), 7))
dataArray[,1] <- as.numeric(stns177$ccode)+(as.numeric(stns177$scode)/1000)
dataArray[,2] <- as.numeric(stns177$slondec)
dataArray[,3] <- as.numeric(stns177$slatdec)

# Set number of years for rates to be calculated over, start and end dates of
# the series etc
rateLen <- 5
rateLenM1 <- rateLen-1
startYr<-1948
startYrM1<-startYr-1
endYr<-2007
lenYrs<-length(startYr:endYr)
lenYrsM1<-lenYrs-1

startRate<-(startYr+0.5)+rateLen/2
endRate<-(endYr+0.5)-rateLen/2
endRateP1<-endRate+1
lastRate<-endYr-rateLenM1
lenRates<-length(startRate:endRate)
lenRatesP1<-lenRates+1

# 70 % completeness
#pc70<-round(0.7*rateLen)
pc70<-3

# Now get the data from the database
source("get177StationsAnnual.R")

# Calculate the rates 
source("long177Recs.R")

# Correct for pressure changes
source("correctSlp2.R")
source("pressureCorrectedRates177_2.R")

# PGR correct these rates
source("pgrCorrectRates2.R")

# Group into the 13 regions of the original study and calculate a global mean
source("group13Regions3.R")

# Compare with rates calculated previously
global <- new.env()
load("nineStationsAndGlobalRates.RData", envir = global)

x11()
plot(s177DecadeMidPoints,s177GlobalMean,type='l', lwd=2, ylim=c(-10,11))
lines(global$slpDecadeMidPoints,global$pressureCorrectedMean, lwd=2, col='blue')
lines(global$globalRates55Yrs[,1],global$globalRates55Yrs[,2], lwd=2, col='red')
grid(col="grey30")

# Look at principal components of the array
#pc.pressureCorrectedRates <- princomp(pressureCorrectedRates)
#x11()
#plot(pc.pressureCorrectedRates)
#summary(pc.pressureCorrectedRates)

save(file="s177Rates2.RData", s177DecadeMidPoints, s177GlobalMean,
  s177GlobalSE)
#
# Check completeness of data
#
completeRegions <- vector(mode="integer", length=56)
completeStations <- vector(mode="integer", length=56)
for (i in 1:lenRatesP1) {
  completeRegions[i] <- length(which(is.finite(s177GlobalArray[i,])))
  completeStations[i] <- length(which(is.finite(s177PgrCorrectedRates[i,])))
}
x11()
plot(s177DecadeMidPoints, completeRegions, type="s")
x11()
plot(s177DecadeMidPoints, completeStations, type="s")
