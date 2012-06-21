# R script to pressure correct the sea level rates which have previously been
# PGR corrected

# Pressure corrections are only available up to 1998 meaning the last decade
# is centred on 1993, 2 years less than the sea level data.
# Hence pressure corrected array is length 46, not 48
s177PressureCorrectedRates <- array(NA,dim=c(46,177))
s177PressureCorrectedMean <- array(0,dim=c(46,1))
for (i in 1:177){
# Remember that a positive rate of pressure increase is equivalent to a
# negative rate of sea level rise (i.e. falling sea level). High pressure
# pushes sea level down.
  s177PressureCorrectedRates[,i] <-
    stns177Rates[1:46,i] + s177SlpRates[,i]
}

for (i in 1:46){
  s177PressureCorrectedMean[i] <- mean(s177PressureCorrectedRates[i,],
    na.rm=TRUE)
}

