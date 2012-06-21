# Script to read in the many stations that have at least 12 of the last 16 years
# of data over the period of the satellite altimeters.
# Uses the HadSLP2 pressure data
#
stns <- read.table(file="includedTideGaugesAltimPeriod.txt", sep=";",
col.names=c("ccode", "scode", "sname", "slondec", "slatdec", "numYrs"),
colClasses=c("character","character","character","numeric","numeric","integer"))
#
dataList <- list()
dataArray <- array(NA,dim=c(length(stns[,1]), 7))
dataArray[,1] <- as.numeric(stns$ccode)+(as.numeric(stns$scode)/1000)
dataArray[,2] <- as.numeric(stns$slondec)
dataArray[,3] <- as.numeric(stns$slatdec)

# Set number of years for rates to be calculated over, start and end dates of
# the series etc

numStns <- length(stns$ccode)

rateLen <- 5
rateLenM1 <- rateLen-1
startYr<-1992
startYrM1<-startYr-1
endYr<-2008
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
#source("getStationsAnnual.R")
load("altimeterPeriodStnsAnnual.RData")

# Calculate the rates 
source("~/workspace/RScripts/global_tgs/pentadal/AltimetryPeriod/longRecs.R")

# Correct for pressure changes
source("~/workspace/RScripts/global_tgs/pentadal/AltimetryPeriod/correctSlp.R")
source("~/workspace/RScripts/global_tgs/pentadal/AltimetryPeriod/pressureCorrectedRates.R")

# PGR correct these rates
source("~/workspace/RScripts/global_tgs/pentadal/AltimetryPeriod/pgrCorrectRates.R")

# Group into the 13 regions of the original study and calculate a global mean
source("~/workspace/RScripts/global_tgs/pentadal/AltimetryPeriod/group13Regions.R")

# Compare with rates calculated previously
global <- new.env()
load("~/diskx/177Stations/nineStationsAndGlobalRates.RData", envir = global)

x11()
plot(decadeMidPoints,globalMean,type='l', lwd=2, ylim=c(-10,11))
lines(global$slpDecadeMidPoints,global$pressureCorrectedMean, lwd=2, col='blue')
lines(global$globalRates55Yrs[,1],global$globalRates55Yrs[,2], lwd=2, col='red')
grid(col="grey30")

save(file="altimetryPeriodRates.RData", decadeMidPoints, globalMean,
  globalSE)
#
# Check completeness of data
#
completeRegions <- vector(mode="integer", length=lenRatesP1)
completeStations <- vector(mode="integer", length=lenRatesP1)
completeStationsList <- vector(mode="list", length=lenRatesP1)
for (i in 1:lenRatesP1) {
  completeRegions[i] <- length(which(is.finite(globalArray[i,])))
  completeStations[i] <- length(which(is.finite(pgrCorrectedRates[i,])))
  completeStationsList[[i]] <- which(is.finite(pgrCorrectedRates[i,]))
}
x11()
plot(decadeMidPoints, completeRegions, type="s")
x11()
plot(decadeMidPoints, completeStations, type="s")
#
x11()
plot(decadeMidPoints, globalMean, col='black', lwd=2, type='l', 
  lty='dashed', ylim=range(globalArray, na.rm=T))
lines(decadeMidPoints, globalArray[,1], col='red', lwd=2)
lines(decadeMidPoints, globalArray[,2], col='blue', lwd=2)
lines(decadeMidPoints, globalArray[,3], col='magenta', lwd=2)
lines(decadeMidPoints, globalArray[,4], col='cyan', lwd=2)
lines(decadeMidPoints, globalArray[,5], col='green', lwd=2)
lines(decadeMidPoints, globalArray[,6], col='orange', lwd=2)
lines(decadeMidPoints, globalArray[,7], col='darkblue', lwd=2, lty='dotdash')
lines(decadeMidPoints, globalArray[,8], col='grey', lwd=2, lty='dotdash')
lines(decadeMidPoints, globalArray[,9], col='yellow', lwd=2, lty='dotdash') 
lines(decadeMidPoints, globalArray[,10], col='violet', lwd=2, lty='dotdash')
lines(decadeMidPoints, globalArray[,11], col='brown', lwd=2, lty='dotdash') 
lines(decadeMidPoints, globalArray[,12], col='salmon', lwd=2, lty='dotdash')
lines(decadeMidPoints, globalArray[,13], col='mintcream', lwd=2, lty='dotdash')
lines(decadeMidPoints, globalArray[,13], col='turquoise', lwd=2, lty='dotdash')
#
aviso_stns<-new.env()
load("~/diskx/177StationsUpdate2009/5YrMeans/aviso_stns/aviso_stns_5yrs.RData",
  envir=aviso_stns)
altimetry<-new.env()
load("~/diskx/177StationsUpdate2009/5YrMeans/leulietteAltimetry.RData", 
  envir=altimetry)
aviso<-new.env()
load("~/diskx/altimetry/5yrRates/aviso_ref_trends.RData", envir=aviso)
x11()
plot(decadeMidPoints, globalMean, type='l', lwd=2, lty='dotdash', xlab='Year', ylab='Rate [mm/yr]')
lines(altimetry$time_j1[84:162], altimetry$rates_j1, col='blue', lwd=2)
lines(aviso$trendMidPoints, aviso$trends, col='magenta', lwd=2)
lines(altimetry$time_tx[84:356],altimetry$rates_tx,col='red', lwd=2)
lines(aviso_stns$trendMidPoints, aviso_stns$s177GlobalMean, col='darkgreen', lwd=2)
legend(x=1995, y=-1, 
  legend=c("Altim Period TGs","Leuliette Tx","Leuliette J1","AVISO Global", "AVISO TGs locs"),
  col=c("black",'red','blue','magenta', 'darkgreen'),
  lwd=2, lty=c('dotdash','solid','solid','solid','solid'))
title('Comparison of Altimetry Period Data')
#
x11()
plot(decadeMidPoints, globalMean, type='l', lwd=2, lty='dotdash', 
  xlab='Year', ylab='Rate [mm/yr]', xlim=c(1995,2007))
lines(altimetry$time_j1[84:162], altimetry$rates_j1, col='blue', lwd=2)
lines(aviso$trendMidPoints, aviso$trends, col='magenta', lwd=2)
lines(altimetry$time_tx[84:356],altimetry$rates_tx,col='red', lwd=2)
lines(aviso_stns$trendMidPoints, aviso_stns$s177GlobalMean, col='darkgreen', lwd=2)
legend(x=1995, y=-1, 
  legend=c("Altim Period TGs","Leuliette Tx","Leuliette J1","AVISO Global", "AVISO TGs locs"),
  col=c("black",'red','blue','magenta', 'darkgreen'),
  lwd=2, lty=c('dotdash','solid','solid','solid','solid'))
#
x11()
nine_stns <- new.env()
load("~/diskx/nineStationsUpdate2009/nineStations5YrRates.RData", envir=nine_stns)
plot(decadeMidPoints, globalMean, type='l', lwd=2, lty='dotdash', 
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
