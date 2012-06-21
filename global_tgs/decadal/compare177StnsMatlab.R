# Script to read in the 177 stations included in the Holgate and Woodworth
# (2004) analysis
stns177 <- read.table(file="includedTideGauges.mod.txt", sep=";",
col.names=c("ccode", "scode", "sname", "slondec", "slatdec", "numYrs"),
colClasses=c("character","character","character","numeric","numeric","integer"))

# Now get the data from the database
source("get177StationsAnnual.R")

# Calculate the rates 
source("long177Recs.R")
# Correct for pressure changes
source("correctSlp2.R")
source("pressureCorrectedRates177_2.R")
# pgr correct
source("pgrCorrectRates2.R")

# Group into the 13 regions of the original study and calculate a global mean
source("group13Regions3.R")

# Compare with rates calculated previously
load("nineStationsAndGlobalRates.RData")

plot(s177DecadeMidPoints[1:42],s177GlobalMean,type='l', lwd=2, ylim=c(-3,6))
lines(slpDecadeMidPoints,pressureCorrectedMean, lwd=2, col='blue')
lines(globalRates55Yrs[,1],globalRates55Yrs[,2], lwd=2, col='red')
grid(col="grey30")

# Look at principal components of the array
pc.pressureCorrectedRates <- princomp(pressureCorrectedRates)
x11()
plot(pc.pressureCorrectedRates)
summary(pc.pressureCorrectedRates)

