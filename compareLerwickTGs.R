# Script showing exploration of the rate of sea level at Lerwick and it's
# relation to other TGs
#source('~/workspace/RScripts/getLerwickData2.R')
#source('~/workspace/RScripts/lerwickAltim.R')
#source('~/workspace/RScripts/compareLerwickNAO.R')
# oracle <- new.env()
#load('lerwickData.RData', env=oracle)
lerwickTime <- seq.Date(from=as.Date("1992/10/15"), length.out=144, by="1 month")
x11()
par(family='HersheySans')
plot(lerwickTime,lerwickMonthly[431:574]-mean(lerwickMonthly[431:574], na.rm=TRUE), type='l')
lines(time,lerwickAltim, col='red')
