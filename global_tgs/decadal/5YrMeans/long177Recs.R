# Using data from "get177StationsAnnual.R"
require(stats)
#
stns177Rates <- array(dim = c((lenRatesP1),177))
s177DecadeMidPoints <- c(startRate:endRateP1)

for (i in 1:177){
  for (j in startYr:lastRate) {
    interval <- findInterval(j:(j+rateLenM1),stns177AnnualYrs)
# We only want to calculate rates where the period is 70% complete
    if (length(which(is.finite(stns177Annual[interval, i])))>=pc70){
      fit <- lsfit(stns177AnnualYrs[interval], stns177Annual[interval, i])
      coeffs <- coef(fit)
      stns177Rates[(j-startYrM1),i] <- coeffs[2]
    } else {
      stns177Rates[(j-startYrM1),i] <- NA
    }
  }
}
