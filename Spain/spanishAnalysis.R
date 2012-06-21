###############################################################################
# spainAnalysis.R
# 
# R script to analyze the Spanish sea level data obtained from 
# getSpanishData.R
#
# Author: simonh
###############################################################################

load("~/diskx/polcoms/spain/spanishData.RData")
numYrs <- length(annualYrs)

# Plot annual data from all sites 
x11()
par(family="HersheySans")

for (i in 1:numStations){
  plot(annualYrs, spanishAnnual[,i], type='l', col='blue', ann=F)
  title(main=paste(stations$ccode[i],stations$scode[i],gsub('[[:space:]]+$', '', stations$sname[i])), 
    ylab='Height [mm]', xlab='Year')
  Sys.sleep(4)
}

# Get first and last year of data
firstLast <- array(NA, dim=c(numStations, 3))
for (i in 1:numStations){
# First year of data
  for (j in 1:numYrs){
    if(is.na(spanishAnnual[j,i])) {
      next
    } else {
      firstLast[i,1] <- j
      break
    }
  } 
# Last year of data
  for (j in numYrs:1){
    if(is.na(spanishAnnual[j,i])) {
      next
    } else {
      firstLast[i,2] <- j
      break
    }
  } 
# Length
  firstLast[i,3] <- (firstLast[i,2] - firstLast[i,1]) + 1
}
#> mean(firstLast[,3])
#[1] 26.81481

# Calculate trends at all 27 sites
trends <- vector(mode='numeric', length=numStations)
for (i in 1:numStations){
  junk <- lm(spanishAnnual[,i] ~ annualYrs)
  trends[i] <- junk$coef[2]
}

# Examine the 12 Atlantic sites more carefully
atlantic <- cbind(stations[1:12,],firstLast[1:12,],trends[1:12])
colnames(atlantic) <- c("ccode","scode","sname","first","last","length","trend")

# Calculate running 10 year means for the Atlantic stations
tenYrTrends <- array(NA,dim=c((numYrs-9),12))
midpointYrs <- seq.Date(from=as.Date("1949/1/1"), length=(numYrs-9), by="1 year")
rainbow_colours <- rainbow(12) 

for (i in 1:(numYrs-9)) {
  for ( j in 1:12) {
     if (length(which(is.finite(spanishAnnual[i:(i+9),j])))>=7){
       junk <- lm(spanishAnnual[i:(i+9),j] ~ annualYrs[1:10], na.action = na.exclude)
       tenYrTrends[i,j] <- junk$coef[2]
     } else {
       tenYrTrends[i,j] <- NA
     }     
  }
}

meanTenYrTrend <- rowMeans(tenYrTrends, na.rm=T)
meanTenYrTrendLongRecs <- rowMeans(tenYrTrends[,c(11,8,7,3)], na.rm=T)
#> range(meanTenYrTrend)
#[1] -6.004575 12.744931

x11()
par(family="HersheySans")
for (i in 1:12){
  if (i == 1){
    plot(midpointYrs, tenYrTrends[,i], type='l', col=rainbow_colours[i], ann=F, ylim=c(-25,25))
    title(main='Decadal Rates of Sea Level Change at Spanish Atlantic Stations', 
      ylab='Rate [mm/yr]', xlab='Year')
   } else {
    lines(midpointYrs, tenYrTrends[,i], col=rainbow_colours[i], ann=F)
   }
}
lines(midpointYrs, meanTenYrTrend, lwd=2)

global<-new.env()
load("~/diskx/177Stations/globalRates55Yrs.RData", envir=global)
global$midpointYrs <- seq.Date(from=as.Date("1953/1/1"), to=as.Date("1997/1/1"), by="1 year")
x11()
par(family="HersheySans")

plot(midpointYrs, tenYrTrends[,11]-15, type='l', col='blue', 
  ann=F, ylim=c(-40,30), lwd=2)
#title(main='Decadal Rates of Sea Level Change at Long Atlantic Stations', 
title(ylab='Rate [mm/yr]', xlab='Year')
lines(midpointYrs, tenYrTrends[,8]-5, col='red', lwd=2)
lines(midpointYrs, tenYrTrends[,7], col='magenta', lwd=2)
lines(midpointYrs, tenYrTrends[,3]+10, col='orange', lwd=2)
#lines(midpointYrs, meanTenYrTrendLongRecs, lwd=2, lty=4)
lines(global$midpointYrs, global$globalRates55Yrs[,2], lwd=2, lty=4)
legend(x=as.Date("1955/7/1"), y=-22, 
  legend=c('Santander I','Coruna I','Coruna II', 'Vigo','Global Mean'),
  col=c('blue', 'red', 'magenta', 'orange', 'black'), lwd=2,
  lty=c(1,1,1,1,4))
#dev2bitmap(file="decadalratesAtlSpain.png", res=150)

x11()
par(family="HersheySans")
plot(annualYrs, spanishAnnual[,11]-mean(spanishAnnual[,11], na.rm=T)-100, type='l', col='blue', ann=F, ylim=c(-300,200), lwd=2)
lines(annualYrs, spanishAnnual[,8]-mean(spanishAnnual[,8], na.rm=T)-50, col='red', lwd=2)
lines(annualYrs, spanishAnnual[,7]-mean(spanishAnnual[,7], na.rm=T), col='magenta', lwd=2)
lines(annualYrs, spanishAnnual[,3]-mean(spanishAnnual[,3], na.rm=T)+50, col='orange', lwd=2)
legend(x=1980, y=-150, 
  legend=c('Santander I','Coruna I','Coruna II', 'Vigo'),
  col=c('blue', 'red', 'magenta', 'orange'), lwd=2)
title(ylab='Height [mm]', xlab='Year')
#dev2bitmap(file="annualDataAtlSpain.png", res=150)

save(file="atlanticSpain.RData", atlantic, tenYrTrends, midpointYrs, meanTenYrTrend,meanTenYrTrendLongRecs)

# Extend analysis used for the British Isles to Atlantic Spain using POLCOMS
# run
