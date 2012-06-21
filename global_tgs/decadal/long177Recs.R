# Using data from "get177StationsAnnual.R"
require(stats)
#
stns177Rates <- array(dim = c(56,177))
s177DecadeMidPoints <- array(dim = c(56,1))

for (i in 1:177){
  for (j in 1948:1998) {
    interval <- findInterval(j:(j+9),stns177AnnualYrs)
# We only want to calculate rates where the decade is 70% complete
    if (length(which(is.finite(stns177Annual[interval, i])))>=7){
      fit <- lsfit(stns177AnnualYrs[interval], stns177Annual[interval, i])
      coeffs <- coef(fit)
      stns177Rates[(j-1947),i] <- coeffs[2]
    } else {
      stns177Rates[(j-1947),i] <- NA
    }
    s177DecadeMidPoints[j-1947] <- j+4.5
  }
}
#
