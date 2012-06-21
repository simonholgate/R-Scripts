# R script to take the decadal rates of SLR and calculate principal components for
# comparsion with "synthetic" time series from model output

# Annual data is in stns177Annual
# Rates are in stns177Rates

# Avoid missing value problems by only choosing stations which are complete
completeStns <- vector(mode="integer", length=0)
for (i in 1:177){
  if (length(which(is.na(stns177Rates[1:46,i]))) == 0){
    completeStns <- c(completeStns,i)
  }
}

completeStnsRates <- stns177Rates[1:46,completeStns]

numCompleteStns <- length(completeStns)

# Detrend each time series
dtCompleteStnsRates <- array(NA, dim=c(46,numCompleteStns))
for (i in 1:numCompleteStns){
  m <- mean(completeStnsRates[,i])
  message(as.character(m))
  dtCompleteStnsRates[,i] <- completeStnsRates[,i] - m
}

# Find the covariance matrix
R <- cov(dtCompleteStnsRates)

# Find eigenvalues and eigenvectors
eigenV <- eigen(R)

# Find the explained variances
PCvars <- diag(eigenV$vectors)/sum(diag(eigenV$vectors))

# Highest PCvar is EOF1
EOForder <- rank(PCvars)
