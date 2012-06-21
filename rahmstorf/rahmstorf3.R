# Script to reproduce the data in Rahmstorf (ScienceExpress, 2006) from Church
# and White along with GISS global temperature data

# Load GISS global surface temperature data
load("~/data/GISSGlobalTemps/globalAnnualMeanTemps.RData")

# Load Church and White global sea level reconstruction
cw <- read.csv("church_white_grl_gmsl.lis", header=FALSE,
col.names=c("years", "GMSL", "SE"))
# The file has three columns:
# Column 1: time in years (1870.0417 = Jan 1870 1870.1250 = Feb 1870 etc)
# Column 2: GMSL in millimetres (the mean has been removed from this series)
# Column 3: One-sigma error in millimetres. This corresponds to the half-width
# of the darker-grey band on the bottom trace of figure 2(a) of the above paper.
# That is, the distance from the black GMSL line to the edge of the darker-grey
# error band.

# Produce annual mean sea levels from Church and White
# Years runs from 1870 to 2001
cwMonthly <- cw$GMSL
dim(cwMonthly) <- c(12,132)
cwAnnual <- colMeans(cwMonthly)
cwYears <- 1870.5:2001.5

# Create data frame so that were dealing with the same data consistently
data.all <- data.frame(years=c(1880:2001), tanom=annualMeanTemp[1:122], 
  sl=cwAnnual[11:132])

# Filter the time series first
filtCWAnnual<-filter(data.all$sl,rep(1/15,15), method='convolution',sides=2)
filtAnnMeanTemp<-filter(data.all$tanom,rep(1/15,15), method='convolution',
  sides=2)

diffFiltCWAnnual <- diff(filtCWAnnual)

junk <- lm(diffFiltCWAnnual[1:61] ~ filtAnnMeanTemp[2:62])
pred1 <- coef(junk)[2]*(filtAnnMeanTemp[2:122]-(-coef(junk)[1]/coef(junk)[2]))
junk <- lm(diffFiltCWAnnual[62:121] ~ filtAnnMeanTemp[63:122])
pred2 <- coef(junk)[2]*(filtAnnMeanTemp[2:122]-(-coef(junk)[1]/coef(junk)[2]))
predR <- 3.4*(filtAnnMeanTemp[2:122]+0.5)
plot(data.all$years[2:122], diffFiltCWAnnual, ann=FALSE, lwd=2, type='l')
lines(data.all$years[2:122], pred1, col='red', lwd=2)
lines(data.all$years[2:122], pred2, col='blue', lwd=2)
lines(data.all$years[2:122], predR, col='magenta', lwd=2)

# rms errors
rmsR <- sqrt(mean((diffFiltCWAnnual - predR)^2,na.rm=TRUE))
rmsR1st <- sqrt(mean((diffFiltCWAnnual[1:61] - predR[1:61])^2,na.rm=TRUE))
rmsR2nd <- sqrt(mean((diffFiltCWAnnual[62:121] - predR[62:121])^2,na.rm=TRUE))
rms1 <- sqrt(mean((diffFiltCWAnnual - pred1)^2,na.rm=TRUE))
rms11st <- sqrt(mean((diffFiltCWAnnual[1:61] - pred1[1:61])^2,na.rm=TRUE))
rms12nd <- sqrt(mean((diffFiltCWAnnual[62:121] - pred1[62:121])^2,na.rm=TRUE))
rms2 <- sqrt(mean((diffFiltCWAnnual - pred2)^2,na.rm=TRUE))
rms21st <- sqrt(mean((diffFiltCWAnnual[1:61] - pred2[1:61])^2,na.rm=TRUE))
rms22nd <- sqrt(mean((diffFiltCWAnnual[62:121] - pred2[62:121])^2,na.rm=TRUE))
#> c(rmsR,rmsR1st,rmsR2nd)
#[1] 0.6338026 0.5868715 0.6782877
#> c(rms1,rms11st,rms12nd)
#[1] 0.6401727 0.5606655 0.7121080
#> c(rms2,rms21st,rms22nd)
#[1] 0.8068231 0.9398007 0.6436764

# Remove the trend from the surface temperatures
junk <- lm(filtAnnMeanTemp ~ data.all$years)
trendTemps <- coef(junk)
detrendTemps <- junk$resid
# Remove the trend from the sea levels
junk <- lm(filtCWAnnual ~ data.all$years)
trendCWAnnual <-  coef(junk)
detrendCWAnnual <- junk$resid

# Remove the trend from the 1st half of the surface temperatures
junk <- lm(filtAnnMeanTemp[1:61] ~ data.all$years[1:61])
trendTemps1st <- coef(junk)
detrendTemps1st <- junk$resid
# Remove the trend from the sea levels
junk <- lm(filtCWAnnual[1:61] ~ data.all$years[1:61])
trendCWAnnual1st <-  coef(junk)
detrendCWAnnual1st <- junk$resid

# Remove the trend from the 2nd half of the surface temperatures
junk <- lm(filtAnnMeanTemp[62:122] ~ data.all$years[62:122])
trendTemps2nd <- coef(junk)
detrendTemps2nd <- junk$resid
# Remove the trend from the sea levels
junk <- lm(filtCWAnnual[62:122] ~ data.all$years[62:122])
trendCWAnnual2nd <-  coef(junk)
detrendCWAnnual2nd <- junk$resid

# Calculate annual rates of sea level change
diffDtCWAnnual <- diff(detrendCWAnnual)
diffCWAnnual <- diff(cwAnnual)
diffDtCWAnnual1st <- diff(detrendCWAnnual1st)
diffDtCWAnnual2nd <- diff(detrendCWAnnual2nd)

# Do 5 year binning
junk5y <- diffDtCWAnnual[2:121]
dim(junk5y) <- c(5,24) 
diffDtCWAnn5Yr <- colMeans(junk5y)

junk5y <- detrendTemps[2:121]
dim(junk5y) <- c(5,24) 
detrendTemps5Yr <- colMeans(junk5y)

junk5y <- diffCWAnnual[2:121]
dim(junk5y) <- c(5,24) 
diffCWAnn5Yr <- colMeans(junk5y)

junk5y <- filtAnnMeanTemp[2:121]
dim(junk5y) <- c(5,24) 
filtAnnMeanTemp5Yr <- colMeans(junk5y)

junk5y <- diffDtCWAnnual1st[1:60]
dim(junk5y) <- c(5,12) 
diffDtCWAnn5Yr1st <- colMeans(junk5y)

junk5y <- detrendTemps1st[2:61]
dim(junk5y) <- c(5,12) 
detrendTemps5Yr1st <- colMeans(junk5y)

junk5y <- diffDtCWAnnual2nd[1:60]
dim(junk5y) <- c(5,12) 
diffDtCWAnn5Yr2nd <- colMeans(junk5y)

junk5y <- detrendTemps2nd[2:61]
dim(junk5y) <- c(5,12) 
detrendTemps5Yr2nd <- colMeans(junk5y)

#x11()
#plot(data.all$years,data.all$tanom, col='blue',type='l', lwd=2, ann=FALSE)
#lines(data.all$years,filtAnnMeanTemp,lwd=2, col='red')
#title("GISS global surface temperature from Hansen et al", xlab="Year", 
#  ylab="Temperature [C]")
#legend("topleft",c("Annual mean","Filtered with 15 yr running mean"), 
#  col=c("blue","red"), lty=1, lwd=2)

#x11()
#plot(data.all$years,data.all$sl, type='l', lwd=2, col='blue', ann=FALSE)
#lines(data.all$years,filtCWAnnual, lwd=2, col='red')
#title("Annual mean sea level from Church and White", xlab="Year", 
#  ylab="Global mean sea level [mm]")
#legend("topleft",c("Annual mean","Filtered with 15 yr running mean"), 
#  col=c("blue","red"), lty=1, lwd=2)

# Use half the existing data to predict the other half.
#> lm(diffDtCWAnn5Yr ~ detrendTemps5Yr)
#Coefficients:
#    (Intercept)  detrendTemps5Yr  
#       -0.01872          6.71877
#
#> lm(diffDtCWAnnual ~ detrendTemps[1:107])
#Coefficients:
#        (Intercept)  detrendTemps[1:107]
#           -0.01375              6.68729

#x11()
#plot(detrendTemps[1:107],diffDtCWAnnual)
#points(detrendTemps[61:114],diffDtCWAnnual[54:107], col='blue')
#points(detrendTemps[8:60],diffDtCWAnnual[1:53], col='red')
#lines(detrendTemps,detrendTemps*6.71877-0.01872,lwd=2,col='magenta')
#lines(detrendTemps,detrendTemps*6.68729-0.01375,lwd=2,col='orange')

#> lm(diffDtCWAnnual1st ~ detrendTemps1st[1:53])
#Coefficients:
#       (Intercept)  detrendTemps[1:53]
#         0.2129                 8.6130
#x11()
#plot(data.all$years[8:60],diffDtCWAnnual1st, type='l')
#lines(data.all$years[8:60],detrendTemps1st[1:53]*8.6130+0.2129,lwd=2,
# col='magenta')

#> lm(diffDtCWAnnual2nd ~ detrendTemps2nd[1:53])
#Coefficients:
#          (Intercept)  detrendTemps2nd[1:53]
#               0.1619                 5.7867
#x11()
#plot(data.all$years[61:113],diffDtCWAnnual2nd, type='l')
#lines(data.all$years[61:113],detrendTemps2nd[1:53]*5.7867+0.1619,lwd=2,
# col='magenta')

# Using:
# dH/dt = a (T - To)
predA <- 6.68729*detrendTemps-0.01375+trendCWAnnual[2]
rmsA <- sqrt(mean((diffFiltCWAnnual[8:115] - predA)^2,na.rm=TRUE))
rmsA1st <- sqrt(mean((diffFiltCWAnnual[8:61] - predA[1:54])^2,na.rm=TRUE))
rmsA2nd <- sqrt(mean((diffFiltCWAnnual[62:115] - predA[55:108])^2,na.rm=TRUE))

predF <- 8.6130*detrendTemps+0.2129+trendCWAnnual1st[2]
rmsF <- sqrt(mean((diffFiltCWAnnual[8:115] - predF)^2,na.rm=TRUE))
rmsF1st <- sqrt(mean((diffFiltCWAnnual[8:61] - predF[1:54])^2,na.rm=TRUE))
rmsF2nd <- sqrt(mean((diffFiltCWAnnual[62:115] - predF[55:108])^2,na.rm=TRUE))

predS <- 5.7867*detrendTemps+0.1619+trendCWAnnual2nd[2]
rmsS <- sqrt(mean((diffFiltCWAnnual[8:115] - predS)^2,na.rm=TRUE))
rmsS1st <- sqrt(mean((diffFiltCWAnnual[8:61] - predS[1:54])^2,na.rm=TRUE))
rmsS2nd <- sqrt(mean((diffFiltCWAnnual[62:115] - predS[55:108])^2,na.rm=TRUE))

#plot(data.all$years[2:122], diffFiltCWAnnual, ann=FALSE, lwd=2, type='l')
#lines(data.all$years[2:122], predR, col='magenta', lwd=2)
#lines(data.all$years[8:115], predF, col='red', lwd=2)
#lines(data.all$years[8:115], predS, col='orange', lwd=2)
#lines(data.all$years[8:115], predA, col='blue', lwd=2)
#title(main='Rate of sea level change and reconstructions from SST', 
#  xlab='Year', ylab='Rate of sea level change [mm/yr]')
#legend("topleft", c("Original data","Rahmstorf reconstruction",
#  "Detrended recon", "1st half recon", "2nd half recon"), lwd=2,
#  col=c('black','magenta','blue','red','orange'))

rmsMat <- matrix(c(rmsR, rmsR1st, rmsR2nd, rmsA, rmsA1st, rmsA2nd, rmsF,
  rmsF1st, rmsF2nd, rmsS, rmsS1st, rmsS2nd), ncol=3, byrow=TRUE)
#
#plot(data.all$years[7:114],cumsum(predR[7:114])+data.all$sl[8], lwd=2, 
#  type='l',  ann=FALSE, xlim=c(1880,2000), ylim=c(-130, 120), col='magenta')
#points(data.all$years,data.all$sl, pch=20)
#lines(data.all$years[8:115],cumsum(predF)-(sum(predF[1:27]) - 
#  mean(data.all$sl[8:61])), lwd=2, col='red')
#lines(data.all$years[8:115],cumsum(predS)-(sum(predS[1:81]) - 
#  mean(data.all$sl[62:115])), lwd=2, col='blue')
#lines(data.all$years[8:115],cumsum(predA)+data.all$sl[8], lwd=2, col='orange')
#
source('getR2_3.R')

# Range of values from bootstrapping first half of the data using Simon
# Williams approach

# First expand dH/dt = a(T-T0) into dHdt = aT - aT0 so that aT0 is the
# intercept and a is the slope. Use the bootstrap results.
intercept1st <- -myslopes1st*myintercepts1st
intercept2nd <- -myslopes2nd*myintercepts2nd

meanSlope1st <- mean(myslopes1st, na.rm=TRUE)
meanIntercept1st <- mean(intercept1st, na.rm=TRUE)
meanSlope2nd <- mean(myslopes2nd, na.rm=TRUE)
meanIntercept2nd <- mean(intercept2nd, na.rm=TRUE)
# Calculate intercept on x axis or T0
meanT01st <- -mean(myintercepts1st, na.rm=TRUE)/meanSlope1st
meanT02nd <- -mean(myintercepts2nd, na.rm=TRUE)/meanSlope2nd

#> length(detrendTemps)
#[1] 108

# Use matrices 
A <- cbind(rep(1, 108), detrendTemps)
m1st <- matrix(c(meanSlope1st,meanIntercept1st),nrow=2)
m2nd <- matrix(c(meanSlope2nd,meanIntercept2nd),nrow=2)

S1st <- A %*% m1st
S2nd <- A %*% m2nd

# We use Q for integration
Q <- matrix(0,108,108)
Q[lower.tri(Q, diag=TRUE)] <- 1

Cm1st <- cov(cbind(intercept1st,myslopes1st),use="pairwise.complete.obs")
Cs1st <- A %*% Cm1st %*% t(A)
Cm2nd <- cov(cbind(intercept2nd,myslopes2nd),use="pairwise.complete.obs")
Cs2nd <- A %*% Cm2nd %*% t(A)

# Finally do the integration
Ci1st <- Q%*%Cs1st%*%t(Q)
Ci2nd <- Q%*%Cs2nd%*%t(Q)

ConfInt1st <- sqrt(diag(Ci1st))*2.5
#ConfInt1st <- sqrt(diag(Ci1st))
ConfInt2nd <- sqrt(diag(Ci2nd))*2.5

# Range of values from bootstrapping first half of the data
predRangeMean1st <- meanSlope1st*(detrendTemps-meanT01st)
predRangeMean2nd <- meanSlope2nd*(detrendTemps-meanT02nd)

