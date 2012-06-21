# GIA correct the sea level rates
load("pgrRates.RData")
#
pgrCorrectedRates <- array(NA,dim = c(lenRatesP1,numStns))
pgrCorrection <- array(NA,dim = c(numStns))
for (i in 1:numStns){
  tmp <- as.numeric(stns$ccode[i])+as.numeric(stns$scode[i])/1000
  # Use Kalinigrad (80.181) pgr value for Klaipeda (80.161)
  if (tmp==80.161) tmp <- 80.181
  # Use Santander II pgr value for Santander III
  if (tmp==200.013) tmp <- 200.012
  # Use Cape Ferguson pgr value for Townsville I
  if (tmp==680.051) tmp <- 680.054
  # Use Rosslyn Bay pgr value for Bundaberg, Burnett Heads or Brisbane
  if ((tmp==680.073) || (tmp==680.078)) tmp <- 680.069
  # Use Hillarys pgr value for Carnarvon
  if (tmp==680.479) tmp <- 680.472
  # Use Broome pgr value for Port Hedland
  if (tmp==680.494) tmp <- 680.497
  
  
  j <- which(pgr$csCode == tmp)

  pgrCorrection[i] <- pgr$pgrRate[j]
  pgrCorrectedRates[,i] <- pressureCorrectedRates[,i]-pgr$pgrRate[j]
}
#
# Calculate the mean record
pgrCorrMean <- array(NA,dim = c(lenRatesP1))
pgrCorrSD <- array(NA,dim = c(lenRatesP1))
pgrCorrSE <- array(NA,dim = c(lenRatesP1))
for (i in 1:lenRatesP1){
  pgrCorrMean[i] <- mean(pgrCorrectedRates[i,],na.rm=TRUE)
  pgrCorrSD[i] <- sd(pgrCorrectedRates[i,],na.rm=TRUE)
  pgrCorrSE[i] <- pgrCorrSD[i]/
    sqrt(length(which(is.finite(pgrCorrectedRates[i,]))))
}
#
# Demean stations
demeanedStations <- array(NA,dim = c(lenYrs,numStns))
for (i in 1:numStns){
  demeanedStations[,i] <- 
    stnsAnnual[,i] - mean(stnsAnnual[,i],na.rm=TRUE)
}
#
stnsMean <- array(NA,dim = c(lenYrsM1))
stnsSD <- array(NA,dim = c(lenYrsM1))
stnsSE <- array(NA,dim = c(lenYrsM1))
for (i in 1:lenYrs){
  stnsMean[i] <- mean(demeanedStations[i,],na.rm=TRUE)
  stnsSD[i] <- sd(demeanedStations[i,],na.rm=TRUE)
  stnsSE[i] <- stnsSD[i]/
    sqrt(length(which(is.finite(demeanedStations[i,]))))
}
#
par(las=1)
#quartz(height=6,width=6)
#x11(height=6,width=6)
#postscript("runningMeanDecadalRateSLR.ps")
par(lab=c(9,10,7))
plot(decadeMidPoints[1:lenRatesP1], pgrCorrMean, type="l",ann=FALSE, lwd=2, 
  ylim=c(-6,10), yaxs="i",xaxs="i")
lines(decadeMidPoints[1:lenRatesP1], pgrCorrMean+pgrCorrSE, lty="dashed")
lines(decadeMidPoints[1:lenRatesP1], pgrCorrMean-pgrCorrSE, lty="dashed")
grid(col="black", lwd=1)
title(main="Running mean decadal rate of sea level rise from the altimetry period records.",
xlab="Year", ylab="Rate of SLR [mm/yr]")
#
mean(pgrCorrMean,na.rm=TRUE)
#* (decadeMidPoints[length(decadeMidPoints)]-decadeMidPoints[1])
#
#quartz(height=6,width=6)
x11(height=6,width=6)
#postscript("demeanedAnnualSeaLevel.ps")
plot(stnsAnnualYrs, stnsMean, type="l", ann=FALSE, lwd=2,
ylim=c(-40,40))
lines(stnsAnnualYrs, stnsMean+stnsSE, lty="dashed")
lines(stnsAnnualYrs, stnsMean-stnsSE, lty="dashed")
grid(col="black", lwd=1)
title(main="Demeaned average annual sea level from the altimetry period records.",
xlab="Year", ylab="Demeaned annual sea level [mm]")
#
fit <- lsfit(stnsAnnualYrs, stnsMean)
coeffs <- coef(fit)
meanRate <- coeffs[2] + mean(pgrCorrection)

meanRate * (stnsAnnualYrs[length(stnsAnnualYrs)]-stnsAnnualYrs[1])
#
#
#quartz(height=6,width=6)
x11(height=6,width=6)
#postscript("integralRunningMeanRate.ps")
plot(decadeMidPoints[1:lenRatesP1], 
  cumsum(pgrCorrMean), type="l",ann=FALSE, lwd=2)
lines(decadeMidPoints[1:lenRatesP1], 
  cumsum(pgrCorrMean)+pgrCorrSE, lty="dashed")
lines(decadeMidPoints[1:lenRatesP1], 
  cumsum(pgrCorrMean)-pgrCorrSE, lty="dashed")
grid(col="black", lwd=1)
title(main=
"Integral of the running mean decadal rate of SLR from the altimetry period records.",
xlab="Year", ylab="Sea level [mm]")
#
#quartz(height=6,width=6)
x11(height=6,width=6)
plot(decadeMidPoints[1:lenRatesP1], cumsum(pgrCorrMean),ann=FALSE,
lty="twodash", type="l", xlim=c(startYr,endYr),ylim=c(-5,40),xaxs='i', lwd=2)
grid(col="black")

