# Using data from "get177StationsAnnual2.R"
# V3 uses robust statistics
require(stats)
require(robust)
#
# Add a last decade for the altimetry period 1992-2009
stns177Rates <- array(dim = c(59,177))
s177DecadeMidPoints <- array(dim = c(58,1))

for (i in 1:177){
  for (j in 1948:2000) {
    interval <- findInterval(j:(j+9),stns177AnnualYrs)
# We only want to calculate rates where the decade is 70% complete
    if (length(which(is.finite(stns177Annual[interval, i])))>=7){
      data <- data.frame(Year=stns177AnnualYrs[interval], Mean=stns177Annual[interval, i])
      fit <- lmRob(Mean ~ Year, data=data)
      stns177Rates[(j-1947),i] <- coef(fit)[2]
    } else {
      stns177Rates[(j-1947),i] <- NA
    }
    s177DecadeMidPoints[j-1947] <- j+4.5
  }

  ####################
  # Altimetry period #
  ####################
    interval <- findInterval(1992:2009,stns177AnnualYrs)
# We only want to calculate rates where the decade is 70% complete
    if (length(which(is.finite(stns177Annual[interval, i])))>=12){
      data <- data.frame(Year=stns177AnnualYrs[interval], Mean=stns177Annual[interval, i])
      fit <- lmRob(Mean ~ Year, data=data)
      stns177Rates[59,i] <- coef(fit)[2]
    } else {
      stns177Rates[59,i] <- NA
    }
}
#
