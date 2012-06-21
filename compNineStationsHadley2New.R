#% compNineStationsHadley2.R
#% R script to calculate the pgr corrected sea level over 10 year periods 
#% where there is at least 50 of the last 60 years of data for a selection of stations
#%
#
startYr <- 1904
endYr <- 1995
#
len <- length(startYr:endYr)
#
#% Required years of data in shorter section 
nReqYrShort <- 7
#
pgrFile = '~/data/cabanes/accel/pgr.rates.stns.out.sed'
#
# Replace longCSCode with codes of the 9 stations
#
longCSCode <-
c(960.121,940.071,823.081,840.011,760.031,210.021,170.161,270.061,690.002)
#c(960.121,940.071,823.081,840.011,760.031,210.021,170.161,690.002)
#
pgrDataFile <- read.table(pgrFile)
pgrCSCode <- pgrDataFile[,1]
#
dataList <- list()
nineStationsTrendArray<-array(NA,dim=c(len,11))
#
yearMean<- numeric(0)
yearCount<- numeric(0)
yearSD<- numeric(0)
stationSet<- numeric(0)
#
diffCount <- 0
#
for (i in startYr:endYr) {
  diffTrend<- numeric(0)
  shortTrendFile = paste('~/data/cabanes/accel/hadleyAPC2/temp27.',as.character(i),'-',
    as.character(i+9),'.hadley2.apc.sed', sep="")
#
#% Open file of stations to get data from
  shortDataFile <- read.table(shortTrendFile)
  shortYears <- shortDataFile[,1]
#
#% Ensure that there are at least nReqYrShort years in the record
  longEnough <- which(shortYears>=nReqYrShort)
  shortDataFile <- shortDataFile[longEnough,]
#  
  shortCSCode <- shortDataFile[,2]+shortDataFile[,3]/1000
  shortCCode <- shortDataFile[,2]
#
#% Only keep shortData stations which have entry in longData
  keepStations <- numeric(0)
  for (j in 1:length(shortCSCode)){
    foundStation <- which(longCSCode == shortCSCode[j])
    if (length(foundStation)>0){
        keepStations <- c(keepStations,j)
    }
  }
#
  shortDataFile <- shortDataFile[keepStations,]
#
  shortCCode <- shortDataFile[,2]
  shortSCode <- shortDataFile[,3]
  shortCSCode <- shortCCode+shortSCode/1000
  stationSet <- union(stationSet,shortCSCode)
#
  dataArray <- array(NA,dim=c(length(shortCSCode), 7))
  shortLon <- shortDataFile[,4]
  shortLat <- shortDataFile[,5]
  shortTrend <- as.numeric(as.vector(shortDataFile[,6]))
  dataArray[,1] <- shortCSCode
  dataArray[,2] <- shortLon
  dataArray[,3] <- shortLat
  dataArray[,4] <- shortDataFile[,1]
  dataArray[,5] <- shortTrend
#
#% Now check there is a corresponding record in the pgr file and delete the
#% others from the pgr file
  keepStations<-numeric(0)
  for (j in 1:length(shortCSCode)) {
    foundStation <- which(pgrCSCode == shortCSCode[j])
    if (length(foundStation)>0) {
        keepStations <- c(keepStations,foundStation)
    } else {
      if (shortCSCode[j]==690.012){# Auckland II
        foundStation <- which(pgrCSCode==690.011) # Auckland
        keepStations <- c(keepStations,foundStation)
      } else {
        stop(paste('shortCSCode ',as.character(shortCSCode[j]),
          ' not found in pgr file', sep=""))
      }
    }
  }
#        
  pgr <- pgrDataFile[keepStations,2]
  dataArray[,6] <- pgr
  diffTrend <- shortTrend - pgr
  dataArray[,7] <- diffTrend
#
  for (j in 1:length(shortCSCode)){
    junk <- which(longCSCode==shortCSCode[j])
    nineStationsTrendArray[(i-startYr+1),junk] <- diffTrend[j]
  }
#
  diffCount <- diffCount + length(pgr)
#
  yearMean <- c(yearMean,mean(diffTrend))
  yearCount <- c(yearCount,length(diffTrend))
  yearSD <- c(yearSD,sd(diffTrend))

# 
  dataList[(i-startYr+1)] <- list(data=dataArray)

} # for loop

nineStationsTrendArray[,10] <- yearMean
#
message(paste('Total pgr corrected records: ', as.character(diffCount),
  sep=""));

#% Calculate mean difference of selections
#message(paste('Mean trend of pgr corrected records is ',   
#  as.character(mean(diffTrend)), sep=""))
yearSE <- yearSD/sqrt(yearCount)
# 
# Compare regions with the mean: which region is the most highly correlated?
X <- 1909:2000
#data <- cbind(t(X),t(yearMean),t(c(mSc,mScS)),t(c(mnE,mnES)),t(c(mMe,mMeS)),
#        t(c(meA,meAS)), t(c(musG,musGS)),t(c(mwA,mwAS)),t(c(mswA,mswAS)),
#        t(c(mseA,mseAS)),t(c(mAus,mAusS)), t(c(mpI,mpIS)),t(c(meP,mePS)), 
#        t(c(mseP,msePS)), t(c(mnwA,mnwAS)))
message(paste('Mean of lower and upper halves: ',as.character(mean(yearMean[1:46])),' ',as.character(mean(yearMean[47:92]))))
message(paste('Mean total: ',as.character(mean(yearMean))))
#
cumYearMean<-cumsum(yearMean)
rateCumYearMean <- lm(cumYearMean ~ X, x=TRUE)
#
# Correlations
corrArray <- array(NA, dim=c(10,10))
confArray <- array(NA, dim=c(10,10,2))
for (i in 1:10) {
  for (j in 1:10) {
    junk <- cor.test(nineStationsTrendArray[,i], nineStationsTrendArray[,j])
    confArray[i,j,1] <- junk$conf.int[1]
    confArray[i,j,2] <- junk$conf.int[2]
    if((junk$estimate >= confArray[i,j,1]) && 
      (junk$estimate <= confArray[i,j,2])) {
        corrArray[i,j] <- junk$estimate
    }
  }
}
library(fields)
corrArray <- upper.tri(corrArray)*corrArray
j <- which(corrArray==0)
corrArray[j] <- NA
image.plot(1:10,1:10,corrArray)
#
# Mean excluding Trieste - column 8
for (i in 1:92) {
  nineStationsTrendArray[i,11] <- 
    mean(nineStationsTrendArray[i,c(1:7,9)], na.rm=TRUE)
}
#lm(cumsum(nineStationsTrendArray[,11]) ~ X)
for (i in 1:10) {
  message(paste(as.character(i),
    as.character(mean(nineStationsTrendArray[,i], na.rm=TRUE))))
}
#
# Spline fit Trieste
#
# Remove duplicate values from Trieste
j <- which(!duplicated(nineStationsTrendArray[, 8]))
splineJunk<-interpSpline(X[j], nineStationsTrendArray[j,8], na.action=na.omit)
#plot( predict( splineJunk, X), type = "l", col="red" )
#lines(X,nineStationsTrendArray[,8], col='blue')

# Thinking about serial correlation, reduce the number of degrees of freedom
# by the amount indicated by the ACF.
# Use longest section of unbroken data from each station for calculation
annualData <- array(NA,dim=c(101,9))
annualData[1:100,] <- nineStationsAnnual
len <- vector(mode="numeric", length=9)
unBroken <- array(NA,dim=c(9,2))
for (i in 1:9) {
  nas <- which(is.na(annualData[,i]))
  nas <- c(1,nas,length(annualData[,i]))
  lennas <- length(nas)
  lens <- vector(mode="numeric", length=(lennas-1))
# First length is from beginning - can be 0
#  lens[1] <- nas[1]-1
  for (j in 1:(lennas-1)) {
    lens[j] <- nas[j+1]-nas[j]
  }
# last length is from end to last NA - can be 0
#  lens[lennas] <- length(annualData[,i])-nas[j]
  maxLen <- which.max(lens)
  if(maxLen == 1) {
    unBroken[i,1] <- 1
  } else {
    unBroken[i,1] <- nas[maxLen]+1
  }
  unBroken[i,2] <- nas[maxLen+1]-1
  len[i] <- length(unBroken[i,1]:unBroken[i,2])
}
#
# Effective degrees of freedom
n_pr<-vector(mode="integer", length=9)
# Calculate acf
for (i in 1:9){
  junk <- acf(nineStationsAnnual[unBroken[i,1]:unBroken[i,2],i], plot=FALSE)
  r1<-junk$acf[2]
  n_pr[i] <- round(((unBroken[i,2]-unBroken[i,1])+1)*(1-r1)/(1+r1))
}
# For "global" curve - see Maul and Martin 1993
junk <- acf(cumYearMean, plot=FALSE)
r1<-junk$acf[2]
n_prCYM <- round(length(cumYearMean)*(1-r1)/(1+r1))
#
# Calculate linear means from annual data
nineStationsArray <- array(NA,dim=c(9,9))
confArrayTrend <- array(NA, dim=c(9,2))
for (i in 1:9){
  fit <- lm(nineStationsAnnual[,i] ~ annualYrs, x=TRUE)
  junk <- confint(fit)
  confArrayTrend[i,] <- junk[2,]
  nineStationsArray[i,1] <- fit$coeff[2]
  nineStationsArray[i,2] <- pgrCorrection[i]
  nineStationsArray[i,3] <- nineStationsArray[i,1] - nineStationsArray[i,2]
# RMS of each fit
  nineStationsArray[i,5] <- sqrt(sum(fit$residuals^2)/(length(fit$residuals)-2))
# SE of each fit 
  nineStationsArray[i,6] <- sd(fit$residuals, na.rm=TRUE)/
    sqrt(length(which(is.finite(nineStationsAnnual[,i]))))
  bot <- sum(fit$x[,2]^2) - sum(fit$x[,2])*sum(fit$x[,2])/length(fit$x[,2])
  nineStationsArray[i,7] <- nineStationsArray[i,5]/sqrt(bot)
# RMS of each fit with corrected degrees of freedom
  nineStationsArray[i,8] <- sqrt(sum(fit$residuals^2)/n_pr[i])
# SE of each fit with corrected degrees of freedom - see Topping p106
  nineStationsArray[i,9] <- nineStationsArray[i,8]/sqrt(bot)
}
# SE for global curve with corrected degrees of freedom
botCYM <- sum(rateCumYearMean$x[,2]^2) - sum(rateCumYearMean$x[,2])*sum(rateCumYearMean$x[,2])/length(rateCumYearMean$x[,2])
seCYM <- sqrt(sum(rateCumYearMean$residuals^2)/n_prCYM)/sqrt(botCYM)

nineStationsMeanTrend <- mean(nineStationsArray[,3])
for (i in 1:9){
  nineStationsArray[i,4] <- (nineStationsArray[i,3] - nineStationsMeanTrend)^2
}
RMS <- sqrt(sum(nineStationsArray[,4])/(length(nineStationsArray[,4])-2))
#
# Comparison of rates of first and second halves
#
# Look at change in rates over the period
meanRates <- array(NA,dim=c(3,10))
# Calculate rates of first and second halves of record fixing the intercept
# point between them, using matrices
A <- matrix(data=0, nrow=100, ncol=3)
# Fill matrix with ones in the first column (the intercept). The second and
# third columns represent the first and second half rates, respectively
A[,1] <- 1
A[1:50,2] <- -49.5:-0.5
A[51:100,3] <- 0.5:49.5
for ( i in 1:9 ){
# We must delete any rows where data is missing
  j <- which(is.finite(nineStationsAnnual[,i]))
  Atemp <- A[j,]
  junk <- solve(t(Atemp)%*%Atemp)%*%t(Atemp)%*%nineStationsAnnual[j,i]
  meanRates[1,i] <- junk[2,1]
  meanRates[2,i] <- junk[3,1]
}
# Now look at "global" mean
#cumYearMean <- cumsum(yearMean)
A <- matrix(data=0, nrow=92, ncol=3)
# Fill matrix with ones in the first column (the intercept). The second and
# third columns represent the first and second half rates, respectively
A[,1] <- 1
A[1:46,2] <- -45.5:-0.5
A[47:92,3] <- 0.5:45.5
j <- which(is.finite(cumYearMean))
Atemp <- A[j,]
junk <- solve(t(Atemp)%*%Atemp)%*%t(Atemp)%*%cumYearMean
intercept12 <- junk[2,1]
meanRates[1,10] <- junk[2,1]
meanRates[2,10] <- junk[3,1]
# Calculate difference in rates
meanRates[3,] <- meanRates[2,] - meanRates[1,]
#
# Look at inflation factor
tau<-1:20
r_0<-1
alpha<-0.1
r_tau[tau] <- r_0*exp(-(alpha*tau)^2)
F_I <- 1 + sqrt(pi)*(r_0/alpha)

# Calculate SE of two halves rates
# Account for serial correlation in deg of freedom
n_prCYM12 <- round(46*(1-r1)/(1+r1))
botCYM12 <- sum(rateCumYearMean$x[1:46,2]^2) - sum(rateCumYearMean$x[1:46,2])*sum(rateCumYearMean$x[1:46,2])/length(rateCumYearMean$x[1:46,2])
# Need to find the residuals of the fits to the first and second halves
X1<-X[1:46]
X2<-X[47:92]
# Calculate constant in slope equations
# The intercept is fixed at intercept12
# y = mx + c
# y(42) = intercept12
# y(42) = meanRates[1,10]*(X1[46]+X2[1])/2 + c1
c1 <- intercept12 - meanRates[1,10]*(X1[46]+X2[1])/2
cYM1 <- meanRates[1,10]*X1 + c1
# y(43) = intercept12
# y(43) = meanRates[2,10]*X2[1] + c2
c2 <- intercept12 - meanRates[2,10]*(X1[46]+X2[1])/2
cYM2 <- meanRates[2,10]*X2 + c2
# Now calclate the residuals
resid1 <- cumYearMean[1:46]-cYM1
resid2 <- cumYearMean[47:92]-cYM2
#
seCYM1 <- sqrt(sum(resid1^2)/n_prCYM12)/sqrt(botCYM12)
seCYM2 <- sqrt(sum(resid2^2)/n_prCYM12)/sqrt(botCYM12)
#
# Look at mean over 1955-98 to comapre with 177 stations
fit <- lm(cumYearMean[47:90] ~ X[47:90])
mean(yearMean[47:90])
