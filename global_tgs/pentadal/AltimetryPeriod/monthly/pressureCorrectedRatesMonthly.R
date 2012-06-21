# R script to pressure correct the sea level rates which have previously been
# PGR corrected

pressureCorrectedRates <- array(NA,dim=c(lenRatesP1,numStns))
pressureCorrectedMean <- array(0,dim=c(lenRatesP1,1))
for (i in 1:numStns){
# Remember that a positive rate of pressure increase is equivalent to a
# negative rate of sea level rise (i.e. falling sea level). High pressure
# pushes sea level down.
  pressureCorrectedRates[,i] <-
    stnsRates[1:lenRatesP1,i] + slpRates[,i]
}

for (i in 1:lenRatesP1){
  pressureCorrectedMean[i] <- mean(pressureCorrectedRates[i,],
    na.rm=TRUE)
  dataArrayTmp <- dataArray
  dataArrayTmp[,5] <- pressureCorrectedRates[i,]
  dataList[i] <- list(data=dataArrayTmp)
}

