# R script to pressure correct the sea level rates which have previously been
# PGR corrected

s177PressureCorrectedRates <- array(NA,dim=c(56,177))
s177PressureCorrectedMean <- array(0,dim=c(56,1))
for (i in 1:177){
# Remember that a positive rate of pressure increase is equivalent to a
# negative rate of sea level rise (i.e. falling sea level). High pressure
# pushes sea level down.
  s177PressureCorrectedRates[,i] <-
    stns177Rates[1:56,i] + s177SlpRates[,i]
}

for (i in 1:56){
  s177PressureCorrectedMean[i] <- mean(s177PressureCorrectedRates[i,],
    na.rm=TRUE)
  dataArrayTmp <- dataArray
  dataArrayTmp[,5] <- s177PressureCorrectedRates[i,]
  dataList177[i] <- list(data=dataArrayTmp)
}

