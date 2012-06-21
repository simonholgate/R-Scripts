# Using data from "getStationsAnnual.R"
require(stats)
#
stnsRates <- array(dim = c((lenRatesP1),numStns))
decadeMidPoints <- c(startRate:endRateP1)

for (i in 1:numStns){
  for (j in startYr:lastRate) {
    interval <- findInterval(j:(j+rateLenM1),stnsAnnualYrs)
# We only want to calculate rates where the period is 70% complete
    if (length(which(is.finite(stnsAnnual[interval, i])))>=pc70){
      fit <- lsfit(stnsAnnualYrs[interval], stnsAnnual[interval, i])
      coeffs <- coef(fit)
      stnsRates[(j-startYrM1),i] <- coeffs[2]
    } else {
      stnsRates[(j-startYrM1),i] <- NA
    }
  }
}
