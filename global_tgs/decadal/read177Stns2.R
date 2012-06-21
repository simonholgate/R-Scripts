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
#source("group13Regions2.R")
source("group13Regions3.R")

# Compare with rates calculated previously
global <- new.env()
load("nineStationsAndGlobalRates.RData", envir = global)

#s177DecadeMidPoints<-s177DecadeMidPoints[1:44]
x11()
plot(s177DecadeMidPoints,s177GlobalMean,type='l', lwd=2, ylim=c(-3,6))
lines(global$slpDecadeMidPoints,global$pressureCorrectedMean, lwd=2, col='blue')
lines(global$globalRates55Yrs[,1],global$globalRates55Yrs[,2], lwd=2, col='red')
grid(col="grey30")

# Look at principal components of the array
pc.pressureCorrectedRates <- princomp(pressureCorrectedRates)
x11()
plot(pc.pressureCorrectedRates)
summary(pc.pressureCorrectedRates)

save(file="s177Rates2.RData", s177DecadeMidPoints, s177GlobalMean,
  s177GlobalSE)
#
