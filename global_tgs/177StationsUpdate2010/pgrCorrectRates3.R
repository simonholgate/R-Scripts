# V3 uses robust methods
# Sept-2010
#library(robust)

# correct the sea level rates
load("pgrRates.RData")
#
s177PgrCorrectedRates <- array(NA,dim = c(59,177))
s177PgrCorrection <- array(NA,dim = c(177))
for (i in 1:177){
  tmp <- as.numeric(stns177[i,1])+as.numeric(stns177[i,2])/1000
  j <- which(pgr$csCode == tmp)
  s177PgrCorrection[i] <- pgr$pgrRate[j]
  s177PgrCorrectedRates[,i] <- s177PressureCorrectedRates[,i]-pgr$pgrRate[j]
}
#
# Calculate the mean record
s177PgrCorrMean <- array(NA,dim = c(59))
s177PgrCorrSD <- array(NA,dim = c(59))
s177PgrCorrSE <- array(NA,dim = c(59))
# Robust scale and median
s177PgrCorrRS <- array(NA,dim = c(59,2))
for (i in 1:59){
	finiteStns<-which(is.finite(s177PgrCorrectedRates[i,]))
  s177PgrCorrMean[i] <- mean(s177PgrCorrectedRates[i,],na.rm=TRUE)
  s177PgrCorrSD[i] <- sd(s177PgrCorrectedRates[i,],na.rm=TRUE)
  s177PgrCorrSE[i] <- s177PgrCorrSD[i]/
    sqrt(length(which(is.finite(s177PgrCorrectedRates[i,]))))
  RS <- s_Sn(s177PgrCorrectedRates[i,finiteStns], mu.too = T)
  s177PgrCorrRS[i,]<-RS
}
#
# Demean stations
demeaned177Stations <- array(NA,dim = c(62,177))
for (i in 1:177){
  demeaned177Stations[,i] <- 
    stns177Annual[,i] - mean(stns177Annual[,i],na.rm=TRUE)
}
#
stns177Mean <- array(NA,dim = c(62))
stns177SD <- array(NA,dim = c(62))
stns177SE <- array(NA,dim = c(62))
# Robust scale and median
stns177RS <- array(NA,dim = c(62,2))
for (i in 1:62){
	finiteStns<-which(is.finite(demeaned177Stations[i,]))
  stns177Mean[i] <- mean(demeaned177Stations[i,],na.rm=TRUE)
  stns177SD[i] <- sd(demeaned177Stations[i,],na.rm=TRUE)
  stns177SE[i] <- stns177SD[i]/
    sqrt(length(finiteStns))
  RS <- s_Sn(demeaned177Stations[i,finiteStns], mu.too = T)
  stns177RS[i,]<-RS
}
#
par(las=1)
#quartz(height=6,width=6)
#x11(height=6,width=6)
#postscript("runningMeanDecadalRateSLR.ps")
par(lab=c(9,10,7))
plot(s177DecadeMidPoints[1:58], s177PgrCorrMean[1:58], type="l",ann=FALSE, lwd=2, 
  ylim=c(-14,14), yaxs="i",xaxs="i")
lines(s177DecadeMidPoints[1:58], s177PgrCorrMean[1:58]+s177PgrCorrSE[1:58], lty="dashed")
lines(s177DecadeMidPoints[1:58], s177PgrCorrMean[1:58]-s177PgrCorrSE[1:58], lty="dashed")
grid(col="black", lwd=1)
title(main="Running mean decadal rate of sea level rise from the 177 records.",
xlab="Year", ylab="Rate of SLR [mm/yr]")
#
mean(s177PgrCorrMean,na.rm=TRUE)
#* (decadeMidPoints[length(decadeMidPoints)]-decadeMidPoints[1])
#
#quartz(height=6,width=6)
x11(height=6,width=6)
#postscript("demeanedAnnualSeaLevel.ps")
plot(stns177AnnualYrs, stns177Mean, type="l", ann=FALSE, lwd=2, ylim=c(-100,150))
lines(stns177AnnualYrs, stns177Mean+stns177SE, lty="dashed")
lines(stns177AnnualYrs, stns177Mean-stns177SE, lty="dashed")
grid(col="black", lwd=1)
title(main="Demeaned average annual sea level from the 177 records.",
xlab="Year", ylab="Demeaned annual sea level [mm]")
#
#fit <- lsfit(stns177AnnualYrs, stns177Mean)
#coeffs <- coef(fit)
stns177AnnualDF <- data.frame(Year=stns177AnnualYrs, Mean=stns177Mean)
fit <- lmRob(Mean ~ Year, data=stns177AnnualDF, x=T)
s177MeanRate <- coef(fit)[2] + mean(s177PgrCorrection)

s177MeanRate * (stns177AnnualYrs[length(stns177AnnualYrs)]-stns177AnnualYrs[1])
#
#
#quartz(height=6,width=6)
x11(height=6,width=6)
#postscript("integralRunningMeanRate.ps")
plot(s177DecadeMidPoints[1:58], 
  cumsum(s177PgrCorrMean[1:58]), type="l",ann=FALSE, lwd=2)
lines(s177DecadeMidPoints[1:58], 
  cumsum(s177PgrCorrMean[1:58])+s177PgrCorrSE[1:58], lty="dashed")
lines(s177DecadeMidPoints[1:58], 
  cumsum(s177PgrCorrMean[1:58])-s177PgrCorrSE[1:58], lty="dashed")
lines(s177DecadeMidPoints[1:58], 
		cumsum(s177PgrCorrRS[1:58,1]), lwd=2, col='red')
lines(s177DecadeMidPoints[1:58], 
		cumsum(s177PgrCorrRS[1:58,1])+s177PgrCorrRS[1:58,2], 
		lty="dashed", col='red')
lines(s177DecadeMidPoints[1:58], 
		cumsum(s177PgrCorrRS[1:58,1])-s177PgrCorrRS[1:58,2], 
		lty="dashed", col='red')
grid(col="black", lwd=1)
title(main=
"Integral of the running mean decadal rate of SLR from the 177 records.",
xlab="Year", ylab="Sea level [mm]")
#
#quartz(height=6,width=6)
x11(height=6,width=6)
plot(s177DecadeMidPoints[1:58], cumsum(s177PgrCorrMean[1:58]),ann=FALSE,
lty="twodash", type="l", xlim=c(1948,2009),ylim=c(-5,90),xaxs='i', lwd=2)
grid(col="black")

