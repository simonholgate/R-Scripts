# Using data from "getStationsMonthly.R"
require(stats)
#
stnsRates <- array(dim = c(lenRatesP1,numStns))
midPoints <- seq(from=startRate, to=endRateP1, by=1/12)

for (i in 1:numStns){
  for (j in 1:lenRates) {
    interval <- j:(j+rateLenM1)
# We only want to calculate rates where the period is 70% complete
    if (length(which(is.finite(stnsMonthly[interval, i])))>=pc70){
      fit <- lsfit(stnsMonthlyYrs[interval], stnsMonthly[interval, i])
      coeffs <- coef(fit)
      stnsRates[j,i] <- coeffs[2]
    } else {
      stnsRates[j,i] <- NA
    }
  }
}
