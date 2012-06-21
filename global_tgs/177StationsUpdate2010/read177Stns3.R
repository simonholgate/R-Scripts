# Script to read in the 177 stations included in the Holgate and Woodworth
# (2004) analysis.
# V2 is updated to use the HadSLP2 pressure data
# V3 (Sept-2010) uses downloaded data and station IDs from new database

stns177 <- read.table(file="includedTideGauges.mod.txt", sep=";",
col.names=c("ccode", "scode", "sname", "slondec", "slatdec", "numYrs"),
colClasses=c("character","character","character","numeric","numeric","integer"))
#
dataList177 <- list()
dataArray <- array(NA,dim=c(length(stns177[,1]), 7))
dataArray[,1] <- as.numeric(stns177$ccode)+(as.numeric(stns177$scode)/1000)
dataArray[,2] <- as.numeric(stns177$slondec)
dataArray[,3] <- as.numeric(stns177$slatdec)

# Convert coutry/station codes into ids
source("get177StationIDs.R")

# Now get the data
source("get177StationsAnnual2.R")

# Calculate the rates 
source("long177Recs3.R")

# Correct for pressure changes
source("correctSlp3.R")
source("pressureCorrectedRates177_2.R")

# PGR correct these rates
source("pgrCorrectRates3.R")

# Group into the 13 regions of the original study and calculate a global mean
#source("group13Regions2.R")
source("group13Regions4.R")

# Compare with rates calculated previously
global <- new.env()
load("nineStationsAndGlobalRates.RData", envir = global)

#s177DecadeMidPoints<-s177DecadeMidPoints[1:44]
x11()
plot(s177DecadeMidPoints[1:58],s177GlobalMean[1:58],type='l', lwd=2, ylim=c(-3,7), xlab="Year", ylab="Rate [mm/yr]")
lines(global$slpDecadeMidPoints,global$pressureCorrectedMean, lwd=2, col='blue')
lines(global$globalRates55Yrs[,1],global$globalRates55Yrs[,2], lwd=2, col='red')
points(2000.5,s177GlobalMean[59], col='red', pch=19)
grid(col="grey30")

x11()
plot(s177DecadeMidPoints[1:58],cumsum(s177GlobalMean[1:58]),type='l', lwd=2, 
		ylim=c(-5, 90), xlab="Year", ylab="Sea Level Change [mm]")
lines(s177DecadeMidPoints[1:58],cumsum(s177GlobalMean[1:58])+s177GlobalRS[1:58,2], 
		lwd=2, lty='dashed')
lines(s177DecadeMidPoints[1:58],cumsum(s177GlobalMean[1:58])-s177GlobalRS[1:58,2], 
		lwd=2, lty='dashed')

#lines(s177DecadeMidPoints[1:58],cumsum(s177GlobalRS[1:58,1]), lwd=2, col='red')
#lines(s177DecadeMidPoints[1:58],cumsum(s177GlobalRS[1:58,1])+s177GlobalRS[1:58,2], 
#		lwd=2, col='red', lty='dashed')
#ines(s177DecadeMidPoints[1:58],cumsum(s177GlobalRS[1:58,1])-s177GlobalRS[1:58,2], 
#		lwd=2, col='red', lty='dashed')
grid(col="grey30")

# Look at principal components of the array
#pc.pressureCorrectedRates <- princomp(s177PressureCorrectedRates)
#x11()
#plot(pc.pressureCorrectedRates)
#summary(pc.pressureCorrectedRates)
x11()
plot(s177DecadeMidPoints[1:58], 
		s177GlobalRS_MAD[1:58,1], col='blue', type='l', lwd=2, ylim=c(-10,15))
lines(s177DecadeMidPoints[1:58], 
		s177GlobalRS_MAD[1:58,1]+s177GlobalRS_MAD[1:58,2], 
		lty="dashed", col='blue')
lines(s177DecadeMidPoints[1:58], 
		s177GlobalRS_MAD[1:58,1]-s177GlobalRS_MAD[1:58,2], 
		lty="dashed", col='blue')
grid(col="black", lwd=1)
title(main= "Robust decadal rate of SLR from the 177 records.",
		xlab="Year", ylab="Sea level [mm]")

x11()
png("robustRatesComparison.png")
plot(s177DecadeMidPoints[1:58],s177GlobalMean[1:58],type='l', lwd=2, 
		ylim=c(-10, 15), xlab="Year", ylab="Sea Level Change [mm]")
lines(s177DecadeMidPoints[1:58],s177GlobalMean[1:58]+s177GlobalSD[1:58], 
		lwd=2, lty='dashed')
lines(s177DecadeMidPoints[1:58],s177GlobalMean[1:58]-s177GlobalSD[1:58], 
		lwd=2, lty='dashed')
lines(s177DecadeMidPoints[1:58], 
		s177GlobalRS_MAD[1:58,1], col='red', lwd=2)
lines(s177DecadeMidPoints[1:58], 
		s177GlobalRS_MAD[1:58,1]+s177GlobalRS_MAD[1:58,2], 
		lty="dashed", col='red')
lines(s177DecadeMidPoints[1:58], 
		s177GlobalRS_MAD[1:58,1]-s177GlobalRS_MAD[1:58,2], 
		lty="dashed", col='red')
grid(col="grey30")
legend(1952,14,c("Mean+/-SD","Median+/-MAD"), col=c("black","red"), lty=c(1,1),
		lwd=c(2,2))
dev.off()
#
save(file="s177Rates2.RData", s177DecadeMidPoints, s177GlobalMean,
  s177GlobalSE)
#
