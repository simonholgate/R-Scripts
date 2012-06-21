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
#
x11()
plot(s177DecadeMidPoints, s177GlobalMean, col='black', lwd=2, type='l', 
  lty='dashed', ylim=range(s177GlobalArray, na.rm=T))
lines(s177DecadeMidPoints, s177GlobalArray[,1], col='red', lwd=2)
lines(s177DecadeMidPoints, s177GlobalArray[,2], col='blue', lwd=2)
lines(s177DecadeMidPoints, s177GlobalArray[,3], col='magenta', lwd=2)
lines(s177DecadeMidPoints, s177GlobalArray[,4], col='cyan', lwd=2)
lines(s177DecadeMidPoints, s177GlobalArray[,5], col='green', lwd=2)
lines(s177DecadeMidPoints, s177GlobalArray[,6], col='orange', lwd=2)
lines(s177DecadeMidPoints, s177GlobalArray[,7], col='darkblue', lwd=2, lty='dotdash')
lines(s177DecadeMidPoints, s177GlobalArray[,8], col='grey', lwd=2, lty='dotdash')
lines(s177DecadeMidPoints, s177GlobalArray[,9], col='yellow', lwd=2, lty='dotdash') 
lines(s177DecadeMidPoints, s177GlobalArray[,10], col='violet', lwd=2, lty='dotdash')
lines(s177DecadeMidPoints, s177GlobalArray[,11], col='brown', lwd=2, lty='dotdash') 
lines(s177DecadeMidPoints, s177GlobalArray[,12], col='salmon', lwd=2, lty='dotdash')
lines(s177DecadeMidPoints, s177GlobalArray[,13], col='mintcream', lwd=2, lty='dotdash')
lines(s177DecadeMidPoints, s177GlobalArray[,13], col='turquoise', lwd=2, lty='dotdash')
#
x11()
plot(s177DecadeMidPoints, s177GlobalMean, type='l', lwd=2, lty='dotdash', xlab='Year', ylab='Rate [mm/yr]')
lines(altimetry$time_j1[84:162], altimetry$rates_j1, col='blue', lwd=2)
lines(aviso$trendMidPoints, aviso$trends, col='magenta', lwd=2)
lines(altimetry$time_tx[84:356],altimetry$rates_tx,col='red', lwd=2)
lines(aviso_stns$trendMidPoints, aviso_stns$s177GlobalMean, col='darkgreen', lwd=2)
legend(x=1950, y=-4, 
  legend=c("177 TGs","Leuliette Tx","Leuliette J1","AVISO Global", "AVISO TGs locs"),
  col=c("black",'red','blue','magenta', 'darkgreen'),
  lwd=2, lty=c('dotdash','solid','solid','solid','solid'))
#
x11()
plot(s177DecadeMidPoints, s177GlobalMean, type='l', lwd=2, lty='dotdash', 
  xlab='Year', ylab='Rate [mm/yr]', xlim=c(1995,2007))
lines(altimetry$time_j1[84:162], altimetry$rates_j1, col='blue', lwd=2)
lines(aviso$trendMidPoints, aviso$trends, col='magenta', lwd=2)
lines(altimetry$time_tx[84:356],altimetry$rates_tx,col='red', lwd=2)
lines(aviso_stns$trendMidPoints, aviso_stns$s177GlobalMean, col='darkgreen', lwd=2)
legend(x=1995, y=-2, 
  legend=c("177 TGs","Leuliette Tx","Leuliette J1","AVISO Global", "AVISO TGs locs"),
  col=c("black",'red','blue','magenta', 'darkgreen'),
  lwd=2, lty=c('dotdash','solid','solid','solid','solid'))
#
x11()
nine_stns <- new.env()
load("../../nineStationsUpdate2009/nineStations5YrRates.RData", envir=nine_stns)
plot(s177DecadeMidPoints, s177GlobalMean, type='l', lwd=2, lty='dotdash', 
  xlab='Year', ylab='Rate [mm/yr]', xlim=c(1914,2007), ylim=c(-11,11))
lines(nine_stns$pentadMidPoints, nine_stns$pgrCorrWeightedMean, col='violet', lwd=2)
lines(altimetry$time_j1[84:162], altimetry$rates_j1, col='blue', lwd=2)
lines(aviso$trendMidPoints, aviso$trends, col='magenta', lwd=2)
lines(altimetry$time_tx[84:356],altimetry$rates_tx,col='red', lwd=2)
lines(aviso_stns$trendMidPoints, aviso_stns$s177GlobalMean, col='darkgreen', lwd=2)
legend(x=1920, y=-5, 
  legend=c("177 TGs","9 TGs","Leuliette Tx","Leuliette J1","AVISO Global", "AVISO TGs locs"),
  col=c("black",'violet','red','blue','magenta', 'darkgreen'),
  lwd=2, lty=c('dotdash','solid','solid','solid','solid','solid'))
