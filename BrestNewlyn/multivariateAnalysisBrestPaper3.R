###############################################################################
# multivariateAnalysisBrestPaper.R
# 
# Implementation of Thompson's methods (Thompson, 1980; 1986) to the tide gauge
# time series at Brest
# These are the same methods implemented in multivariateAnalysisNewlyn.R
#
# Version "Paper" is changed so that we use total pressure throughout, as in
# Thompson '86, rather than having it as a forcing term on the RHS (as
# in Thompson '80). This avoids the problem of aliasing.
# We include HadSLP2r data for the first time.
#
# In comparison with multivariateAnalysisBrest3.R, the file is restructured
# to follow the logic of the paper and produce the necessary outputs
#
# Author: simonh
###############################################################################

##################
# Tide Gauge data
##################

# Annual data #
tg_annual<-new.env()
tg_annual$brest <- read.table("~/Dropbox/brestNewlynData/analysis/paper/brestAnnual.rlrdata", col.names=c("Year", "Height", "Flag", "QC"), colClasses=c("numeric","numeric", "character", "character"), sep=";", na.string="-99999")
tg_annual$newlyn <- read.table("~/Dropbox/brestNewlynData/analysis/paper/newlynAnnual.rlrdata", col.names=c("Year", "Height", "Flag", "QC"), colClasses=c("numeric","numeric", "character", "character"), sep=";", na.string="-99999")

# Delete the 1916 value from both time-series so that the two series
# are zero in 1916
tg_annual$newlyn$Height <- tg_annual$newlyn$Height - tg_annual$newlyn$Height[1]
tg_annual$brest$Height <- tg_annual$brest$Height - tg_annual$brest$Height[which(tg_annual$brest$Year==1916)]

tg_annual$brestTS <- ts(data=tg_annual$brest$Height, start=tg_annual$brest$Year[1], end=tg_annual$brest$Year[202], frequency=1)
tg_annual$newlynTS <- ts(data=tg_annual$newlyn$Height, start=tg_annual$newlyn$Year[1], end=tg_annual$newlyn$Year[93], frequency=1)

# Use robust methods
library(robust)
# Calculate rate over entire Newlyn period (1916-2008)
tg_annual$newlyn.lmRob <- lmRob(Height ~ Year, data=tg_annual$newlyn, x=T)
#> summary(tg_annual$newlyn.lmRob)
#
#Call: lmRob(formula = Height ~ Year, data = tg_annual$newlyn)
#
#Residuals:
#         Min           1Q       Median           3Q          Max 
#-55.08867890 -16.59755299  -0.06796805  19.47641170  71.66280011 
#
#Coefficients:
#            Value         Std. Error    t value       Pr(>|t|)     
#(Intercept) -3324.8167012   225.3975135   -14.7509023     0.0000000
#Year            1.7248521     0.1149240    15.0086345     0.0000000
#
#Residual standard error: 25.6063 on 90 degrees of freedom
#1 observation deleted due to missingness 
#Multiple R-Squared: 0.649363 
#
#Test for Bias:
#            statistic   p-value
#M-estimate   1.576701 0.4545941
#LS-estimate  1.032041 0.5968911

tg_annual$newlyn.lm <- lm(Height ~ Year, data=tg_annual$newlyn, x=T)

tg_annual$newlynGap <- tg_annual$newlyn
tg_annual$newlynGap[29:37,2] <- NA
tg_annual$newlynGap.lmRob <- lmRob(Height ~ Year, data=tg_annual$newlynGap, x=T)
#> summary(tg_annual$newlynGap.lmRob)
#
#Call: lmRob(formula = Height ~ Year, data = tg_annual$newlynGap, x = T)
#
#Residuals:
#        Min          1Q      Median          3Q         Max 
#-54.4382075 -17.2897926  -0.4830823  20.1620392  72.2069140 
#
#Coefficients:
#            Value         Std. Error    t value       Pr(>|t|)     
#(Intercept) -3346.6535896   218.9146371   -15.2874821     0.0000000
#Year            1.7354878     0.1115382    15.5595806     0.0000000
#
#Residual standard error: 25.8677 on 81 degrees of freedom
#10 observations deleted due to missingness 
#Multiple R-Squared: 0.670672 
#
#Test for Bias:
#            statistic   p-value
#M-estimate  0.9900018 0.6095704
#LS-estimate 0.1021959 0.9501856

tg_annual$newlyn1920 <- tg_annual$newlyn[4:93,]
tg_annual$newlyn1920.lmRob <- lmRob(Height ~ Year, data=tg_annual$newlyn1920, x=T)
#> summary(tg_annual$newlyn1920.lmRob)
#
#Call: lmRob(formula = Height ~ Year, data = tg_annual$newlyn1920, x = T)
#
#Residuals:
#       Min         1Q     Median         3Q        Max 
#-54.816821 -17.374839  -0.340244  19.069125  72.136333 
#
#Coefficients:
#            Value         Std. Error    t value       Pr(>|t|)     
#(Intercept) -3284.9150560   253.1929433   -12.9739598     0.0000000
#Year            1.7046847     0.1289866    13.2159812     0.0000000
#
#Residual standard error: 25.2971 on 87 degrees of freedom
#1 observation deleted due to missingness 
#Multiple R-Squared: 0.631512 
#
#Test for Bias:
#            statistic   p-value
#M-estimate   1.310681 0.5192652
#LS-estimate  2.017716 0.3646352

# Brest
tg_annual$brest1916 <- tg_annual$brest[110:202,]
tg_annual$brest1916.lmRob <- lmRob(Height ~ Year, data=tg_annual$brest1916, x=T)
#> summary(tg_annual$brest1916.lmRob)
#
#Call: lmRob(formula = Height ~ Year, data = tg_annual$brest1916, x = T)
#
#Residuals:
#       Min         1Q     Median         3Q        Max 
#-64.440392 -21.397379   2.413695  21.416783  78.587177 
#
#Coefficients:
#            Value         Std. Error    t value       Pr(>|t|)     
#(Intercept) -2.768637e+03  2.851473e+02 -9.709498e+00  2.664535e-15
#Year         1.403992e+00  1.452050e-01  9.669034e+00  3.330669e-15
#
#Residual standard error: 31.2865 on 82 degrees of freedom
#9 observations deleted due to missingness 
#Multiple R-Squared: 0.546912 
#
#Test for Bias:
#            statistic   p-value
#M-estimate   3.681634 0.1586877
#LS-estimate -1.919591 1.0000000

tg_annual$brest1916.lm <- lm(Height ~ Year, data=tg_annual$brest1916, x=T)

tg_annual$brest1920 <- tg_annual$brest[114:202,]
tg_annual$brest1920.lmRob <- lmRob(Height ~ Year, data=tg_annual$brest1920, x=T)
#> summary(tg_annual$brest1920.lmRob)
#
#Call: lmRob(formula = Height ~ Year, data = tg_annual$brest1920, x = T)
#
#Residuals:
#       Min         1Q     Median         3Q        Max 
#-64.502480 -21.444073   2.452766  21.559213  68.635069 
#
#Coefficients:
#            Value         Std. Error    t value       Pr(>|t|)     
#(Intercept) -2.751219e+03  3.098006e+02 -8.880612e+00  1.838529e-13
#Year         1.395196e+00  1.575692e-01  8.854495e+00  2.067235e-13
#
#Residual standard error: 30.2054 on 78 degrees of freedom
#9 observations deleted due to missingness 
#Multiple R-Squared: 0.536609 
#
#Test for Bias:
#            statistic    p-value
#M-estimate   5.445354 0.06569865
#LS-estimate -1.087141 1.00000000

# Split Newlyn into two, missing out the years for which Brest has no data
#(1944-1952)
tg_annual$newlyn1st <- tg_annual$newlyn[1:28,]
tg_annual$newlyn1st.lmRob <- lmRob(Height ~ Year, data=tg_annual$newlyn1st, x=T)
#> summary(tg_annual$newlyn1st.lmRob)
#
#Call: lmRob(formula = Height ~ Year, data = tg_annual$newlyn1st)
#
#Residuals:
#        Min          1Q      Median          3Q         Max 
#-42.3303229 -16.0758758   0.8891626  19.4679803  45.8628900 
#
#Coefficients:
#            Value         Std. Error    t value       Pr(>|t|)     
#(Intercept) -4.194660e+03  1.149430e+03 -3.649340e+00  1.158514e-03
#Year         2.175424e+00  5.957085e-01  3.651826e+00  1.151148e-03
#
#Residual standard error: 27.3275 on 26 degrees of freedom
#Multiple R-Squared: 0.35296 
#
#Test for Bias:
#            statistic   p-value
#M-estimate   2.026173 0.3630966
#LS-estimate -3.657567 1.0000000 

tg_annual$newlyn2nd <- tg_annual$newlyn[38:93,]
tg_annual$newlyn2nd.lmRob <- lmRob(Height ~ Year, data=tg_annual$newlyn2nd, x=T)
#> summary(tg_annual$newlyn2nd.lmRob)
#
#Call: lmRob(formula = Height ~ Year, data = tg_annual$newlyn2nd)
#
#Residuals:
#        Min          1Q      Median          3Q         Max 
#-54.3645729 -16.5053531  -0.7995256  20.6354271  71.8490205 
#
#Coefficients:
#            Value         Std. Error    t value       Pr(>|t|)     
#(Intercept) -3.432688e+03  5.246446e+02 -6.542882e+00  2.451315e-08
#Year         1.778641e+00  2.649396e-01  6.713383e+00  1.304331e-08
#
#Residual standard error: 25.7296 on 53 degrees of freedom
#1 observation deleted due to missingness 
#Multiple R-Squared: 0.424608 
#
#Test for Bias:
#            statistic   p-value
#M-estimate  0.5926822 0.7435338
#LS-estimate 1.3345035 0.5131168

# Do the same for Brest
tg_annual$brest1st <- tg_annual$brest[110:137,]
tg_annual$brest1st.lmRob <- lmRob(Height ~ Year, data=tg_annual$brest1st, x=T)
#>  summary(tg_annual$brest1st.lmRob) 
#
#Call: lmRob(formula = Height ~ Year, data = tg_annual$brest1st)
#
#Residuals:
#       Min         1Q     Median         3Q        Max 
#-51.376299 -20.342542   7.834497  25.852077  92.391281 
#
#Coefficients:
#            Value         Std. Error    t value       Pr(>|t|)     
#(Intercept) -5.229346e+03  1.714529e+03 -3.050019e+00  5.210676e-03
#Year         2.681083e+00  8.885611e-01  3.017331e+00  5.642724e-03
#
#Residual standard error: 34.1086 on 26 degrees of freedom
#Multiple R-Squared: 0.23705 
#
#Test for Bias:
#            statistic   p-value
#M-estimate   2.256287 0.3236336
#LS-estimate -1.033872 1.0000000

tg_annual$brest2nd <- tg_annual$brest[147:202,]
tg_annual$brest2nd.lmRob <- lmRob(Height ~ Year, data=tg_annual$brest2nd, x=T)
#> summary(tg_annual$brest2nd.lmRob)
#
#Call: lmRob(formula = Height ~ Year, data = tg_annual$brest2nd)
#
#Residuals:
#       Min         1Q     Median         3Q        Max 
#-57.290964 -12.631163   1.627075  17.286877  82.963155 
#
#Coefficients:
#            Value         Std. Error    t value       Pr(>|t|)     
#(Intercept) -3.722835e+03  5.894965e+02 -6.315279e+00  5.302714e-08
#Year         1.883608e+00  2.975106e-01  6.331229e+00  4.997448e-08
#
#Residual standard error: 24.903 on 54 degrees of freedom
#Multiple R-Squared: 0.473446 
#
#Test for Bias:
#            statistic    p-value
#M-estimate   2.519923 0.28366491
#LS-estimate  4.709127 0.09493492

# Summary
# Newlyn Robust: 1.7248521     0.1149240
# Newlyn: 1.7648     0.1017
# Newlyn1920 Robust: 1.7046847     0.1289866
# NewlynGap Robust: 1.7354878     0.1115382

# Brest1916 Robust: 1.403992e+00  1.452050e-01
# Brest1916: 1.3639     0.1193
# Brest1920 Robust: 1.395196e+00  1.575692e-01

# Newlyn1st Robust: 2.175424e+00  5.957085e-01
# Newlyn2nd Robust: 1.778641e+00  2.649396e-01
# Brest1st Robust: 2.681083e+00  8.885611e-01
# Brest2nd Robust: 1.883608e+00  2.975106e-01

# Newlyn GIA ICE5G VM2 = 0.34 mm/yr (sinking)
# Brest GIA ICE5G VM2 = 0.30

# GIA corrected estimates
(1.883608e+00-0.30)-(1.778641e+00-0.34)
#0.144967
(2.681083e+00-0.30)-(2.175424e+00-0.34)
#0.545659
(1.403992e+00-0.30)-(1.7248521-0.34)
#-0.2808601
(1.395196e+00-0.30)-(1.7046847-0.34)
#-0.2694887

# GPS rates
# Brest -0.57 mm/yr
# Newlyn -0.27 mm/yr

# GPS corrected estimates
(1.883608e+00+0.57)-(1.778641e+00+0.27)
#0.404967
(2.681083e+00+0.57)-(2.175424e+00+0.27)
#0.805659
(1.403992e+00+0.57)-(1.7248521+0.27)
#-0.0208601
(1.395196e+00+0.57)-(1.7046847+0.27)
#-0.0094887

# Monthly data #
tg_monthly<-new.env()
load("~/Dropbox/brestNewlynData/newlynBrestTG.RData", envir=tg_monthly)
tg_monthly$time <- seq.Date(from=as.Date("1914/1/15"), to=as.Date("2006/12/15"), by="1 month")

# Delfzijl
tg_monthly$delfzijl <- read.table("~/Dropbox/brestNewlynData/analysis/paper/delfzijlMonthly.metdata", col.names=c("Year", "Height", "Flag", "QC"), colClasses=c("numeric","numeric", "character", "character"), sep=";", na.string="-99999")

rm(junk, selectString, drv, con, i, indices, res, dataIndices, envir=tg_monthly)

# Create annual means
#> tmp<-new.env()
#> tmp$delfzijl <-  tg_monthly$delfzijl$Height
#> dim(tmp$delfzijl) <- c(12,143)
#> tmp$delfzijlAnnual <- colMeans(tmp$delfzijl)
#> tg_annual$delfzijl <- data.frame(c(1865:2007), tmp$delfzijlAnnual)
#> colnames(tg_annual$delfzijl) <- c("Year", "Height")
#> rm(tmp)

tg_annual$delfzijl1916 <- tg_annual$delfzijl[52:143,]
tg_annual$delfzijl1916.lmRob <- lmRob(Height ~ Year, data=tg_annual$delfzijl1916, x=T)
#> summary(tg_annual$delfzijl1916.lmRob)
#
#Call: lmRob(formula = Height ~ Year, data = tg_annual$delfzijl1916, 
#    x = T)
#
#Residuals:
#        Min          1Q      Median          3Q         Max 
#-115.618246  -25.306544   -2.337109   27.070529   81.291655 
#
#Coefficients:
#            Value         Std. Error    t value       Pr(>|t|)     
#(Intercept) -3847.6430893   353.7797367   -10.8758153     0.0000000
#Year            1.9710311     0.1803746    10.9274314     0.0000000
#
#Residual standard error: 39.3837 on 90 degrees of freedom
#Multiple R-Squared: 0.556452 
#
#Test for Bias:
#             statistic   p-value
#M-estimate   1.8498831 0.3965546
#LS-estimate -0.6510794 1.0000000

# Compare just over the period of the model: 15/1/1960-15/12/1999
tg_monthly$mp.brest <- tg_monthly$brestMonthly[intersect(which(tg_monthly$time>=as.Date("1960/1/15")), which(tg_monthly$time<=as.Date("1999/12/15")))]
tg_monthly$mp.time <- tg_monthly$time[intersect(which(tg_monthly$time>=as.Date("1960/1/15")), which(tg_monthly$time<=as.Date("1999/12/15")))]
tg_monthly$mp.brest.orig <- tg_monthly$mp.brest
tg_monthly$mp.brest[which(is.na(tg_monthly$mp.brest))] <- mean(tg_monthly$mp.brest, na.rm=T)

tg_monthly$mp.newlyn <- tg_monthly$newlynMonthly[intersect(which(tg_monthly$time>=as.Date("1960/1/15")), which(tg_monthly$time<=as.Date("1999/12/15")))]
tg_monthly$mp.newlyn.orig <- tg_monthly$mp.newlyn
# Fill missing Newlyn value with the mean
tg_monthly$mp.newlyn[which(is.na(tg_monthly$mp.newlyn))] <- mean(tg_monthly$mp.newlyn, na.rm=T)

# Remove linear trend and intercept from Newlyn and Brest over model period
tg_monthly$mp.brest.lmRob <- lmRob(tg_monthly$mp.brest ~ as.Date(tg_monthly$mp.time))
tg_monthly$mp.newlyn.lmRob <- lmRob(tg_monthly$mp.newlyn ~ as.Date(tg_monthly$mp.time))
tg_monthly$mp.resid.diff <- tg_monthly$mp.newlyn.lmRob$resid-tg_monthly$mp.brest.lmRob$resid

# Correlation of residuals
#> cor.test(tg_monthly$mp.newlyn.lmRob$resid, tg_monthly$mp.brest.lmRob$resid)
#
#	Pearson's product-moment correlation
#
#data:  tg_monthly$mp.newlyn.lm$resid and tg_monthly$mp.brest.lm$resid 
#t = 57.0852, df = 478, p-value < 2.2e-16
#alternative hypothesis: true correlation is not equal to 0 
#95 percent confidence interval:
# 0.9213594 0.9444182 
#sample estimates:
#      cor 
#0.9338524 

# SD of each residual time series and their difference
#> sd(tg_monthly$mp.newlyn.lmRob$resid)
#[1] 80.84692
#> sd(tg_monthly$mp.brest.lmRob$resid)
#[1] 85.5429
#> sd(tg_monthly$mp.resid.diff)
#[1] 30.61029

# Calculate the trend of the residual
tg_monthly$mp.resid.diff.lmRob <- lmRob(tg_monthly$mp.resid.diff ~ as.Date(tg_monthly$mp.time))
#> tg_monthly$mp.resid.diff.lmRob
#
#Call:
#lm(formula = tg_monthly$mp.resid.diff ~ as.Date(tg_monthly$mp.time))
#
#Coefficients:
#        (Intercept)  as.Date(tg_monthly$mp.time)  
#          9.014e-16           -2.691e-19  

# Calculate the trend of the residual
tg_monthly$mp.resid.diff.lmRob <- lmRob(tg_monthly$mp.resid.diff ~ as.Date(tg_monthly$mp.time))
#> tg_monthly$mp.resid.diff.lmRob
#
#Call:
#lm(formula = tg_monthly$mp.resid.diff ~ as.Date(tg_monthly$mp.time))
#
#Coefficients:
#        (Intercept)  as.Date(tg_monthly$mp.time)  
#          9.014e-16           -2.691e-19  

##############
## Met data ##
##############
met<-new.env()
load("~/Dropbox/brestNewlynData/met_data/brestNewlynMetMonthlyMean.RData", envir=met)
# Create annual mean pressures
tmp <- new.env()
tmp$brestMonthlyPressureArray <- met$brestMetMonthlyMean[,1,]
dim(tmp$brestMonthlyPressureArray) <- c(12, 40)
met$brestAnnualPressureInterp <- colMeans(tmp$brestMonthlyPressureArray, na.rm=T)
tmp$newlynMonthlyPressureArray <- met$newlynMetMonthlyMean[,1,]
dim(tmp$newlynMonthlyPressureArray) <- c(12, 40)
met$newlynAnnualPressureInterp <- colMeans(tmp$newlynMonthlyPressureArray, na.rm=T)
load("~/Dropbox/brestNewlynData/pressure/plymouthPressures.RData", envir=tmp)
met$plymouthMonthlyPressure <- tmp$pressure
met$plymouthMonthlyTime <- tmp$time
tmp$plymouthMonthlyPressureArray <- met$plymouthMonthlyPressure
dim(tmp$plymouthMonthlyPressureArray) <- c(12, 155)
met$plymouthAnnualPressure <- colMeans(tmp$plymouthMonthlyPressureArray, na.rm=T)
met$plymouthAnnualTime <- seq.Date(from=as.Date("1850/07/01"), length=155, by ="1 year")

load("~/Dropbox/brestNewlynData/pressure/brestPressures.RData", envir=tmp)
tmp$brestMonthlyPressureArray <- tmp$pressure
dim(tmp$brestMonthlyPressureArray) <- c(12, 52)
met$brestAnnualPressure <- colMeans(tmp$brestMonthlyPressureArray, na.rm=T)
met$brestAnnualTime <- seq.Date(from=as.Date("1951/07/01"), length=52, by ="1 year")

load("~/Dropbox/brestNewlynData/met_data/brestNewlynMetMonthlyMeanInterpPaper.RData", envir=tmp)
# Dims of met$extraMetMonthlyMeanInterp are: 12  5 40 16
# reflecting 5 parameters collected from 16 stations (regular grid) for 12 months over 40 years
# The 5 parameters are:  slp, windE, windN, stressE, stressN
met$gridAnnualPressureArray <- array(NA, dim=c(40,16))
for(i in 1:16){
  tmp$gridMonthlyPressureArray <- tmp$extraMetMonthlyMeanInterp[,1,,i]
  met$gridAnnualPressureArray[,i] <- colMeans(tmp$gridMonthlyPressureArray, na.rm=T)
}

tmp <- new.env()
# HadSLP2r data is 1850-2007 = 158 years = 1896 months. 37 cols - Newlyn, Brest and 35 grid points
load("~/Dropbox/brestNewlynData/pressure/brestNewlynHadSLP2Paper2.RData", envir=tmp)
met$hadnstns <- 37
met$hadSLP2rAnnualPressureArray <- array(NA, dim=c(160,met$hadnstns))
for(i in 1:met$hadnstns){
  tmp$hadSLP2rMonthlyPressureArray <- tmp$slpNewlynStns[,i]
  dim(tmp$hadSLP2rMonthlyPressureArray) <- c(12,160)
  met$hadSLP2rAnnualPressureArray[,i] <- colMeans(tmp$hadSLP2rMonthlyPressureArray, na.rm=T)
}
met$hadSLP2rAnnualTime <- c(1850:2009)
rm(tmp)

#############
## S12R408 ##
#############
# Add model data from S12run408 with realistic deep ocean boundary condition 
# which varies over the 45 year of the run (1960-2004) unlike the previous model
# S12run405
ls(met)S12R408<-new.env()
load("~/Dropbox/brestNewlynData/brestNewlynMonthlyS12R408.RData", envir=S12R408)
S12R408$time <- S12R408$monthsArray

#############
## Model 1 ##
#############
# Follow Thompson (1980) and see if we can improve
# Uses local pressures from Plymouth
model1<-new.env()

# Just use annual pressures over the periods of the model data (Pa) ending in 2000 (end of model pressure) 
model1$plymouthPa <-
  met$plymouthAnnualPressure[intersect(which(met$plymouthAnnualTime>=as.Date("1960/7/01")), which(met$plymouthAnnualTime<=as.Date("1999/7/01")))]
model1$brestPa <-
  met$brestAnnualPressure[intersect(which(met$brestAnnualTime>=as.Date("1960/7/01")), which(met$brestAnnualTime<=as.Date("1999/7/01")))]
model1$time <- c(1960:1999)
tg_annual$brest6099 <- tg_annual$brest[154:193,]
tg_annual$newlyn6099 <- tg_annual$newlyn[45:84,]
# Calculate the total pressure
#model1$brestTotalP <- 1025*9.8*tg_annual$brest6099[,2]/1000 + model1$brestPa
#model1$newlynTotalP <- 1025*9.8*tg_annual$newlyn6099[,2]/1000 + model1$plymouthPa
model1$brestTotalP <- 1025*9.8*tg_annual$brest6099[,2]/1000 + met$brestAnnualPressureInterp
model1$newlynTotalP <- 1025*9.8*tg_annual$newlyn6099[,2]/1000 + met$newlynAnnualPressureInterp
model1$tg.lmRob.brestTotalP <- lmRob(model1$brestTotalP ~ model1$time)
model1$tg.lmRob.newlynTotalP <- lmRob(model1$newlynTotalP ~ model1$time)

model1$Pa <- array(NA,dim=c(40,(16*3)))
# Zero lag
for(i in 1:16){
  model1$Pa[,i] <- as.vector(met$gridAnnualPressureArray[,i])
}
# Lag 1 year
for(i in 17:32){
  model1$Pa[2:40,i] <- as.vector(met$gridAnnualPressureArray[,(i-16)])[1:39]
  model1$Pa[1,i] <- mean(as.vector(met$gridAnnualPressureArray[,(i-16)])[1:39], na.rm=T)
}
# Lag 2 years....
for(i in 33:48){
  model1$Pa[3:40,i] <- as.vector(met$gridAnnualPressureArray[,(i-32)])[1:38]
  model1$Pa[1:2,i] <- mean(as.vector(met$gridAnnualPressureArray[,(i-32)])[1:38], na.rm=T)
}

#########
# Brest #
#########

# What we are going to do now is calculate the sd of the residuals as each of the
# 48 components is added in turn 
model1$tg.sd.brest <- vector(mode="numeric", length=48)

for(i in 1:48){
  model1$data.brest <- data.frame(msl=model1$brestTotalP, t=seq(from=1960,to=1999), Pa=model1$Pa[,i])
  model1$tg.lmRob.brest <- lmRob(msl ~ t + Pa, data=model1$data.brest)
  
#  model1$tg.resid.brest <- model1$tg.lmRob.brest$resid
  model1$tg.sd.brest[i] <- model1$tg.lmRob.brest$scale
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model1$tg.sd.sorted.brest <- sort(model1$tg.sd.brest, index=TRUE)

# Now construct models of the top 9 components
model1$data.brest <- data.frame(msl=model1$brestTotalP,
  t=seq(from=1960,to=1999),
  Pa1=model1$Pa[,model1$tg.sd.sorted.brest$ix[1]], 
  Pa2=model1$Pa[,model1$tg.sd.sorted.brest$ix[2]], Pa3=model1$Pa[,model1$tg.sd.sorted.brest$ix[3]],
  Pa4=model1$Pa[,model1$tg.sd.sorted.brest$ix[4]], Pa5=model1$Pa[,model1$tg.sd.sorted.brest$ix[5]],
  Pa6=model1$Pa[,model1$tg.sd.sorted.brest$ix[6]], Pa7=model1$Pa[,model1$tg.sd.sorted.brest$ix[7]],
  Pa8=model1$Pa[,model1$tg.sd.sorted.brest$ix[8]], Pa9=model1$Pa[,model1$tg.sd.sorted.brest$ix[9]])

model1$tg.lmRob.brest <- lmRob(msl ~ t + Pa1 + Pa2 + Pa3 + Pa4 + Pa5 + Pa6 + Pa7 + Pa8 + Pa9, x=T, y=T,
                               data=model1$data.brest)

# Variance reduction
#> (model1$tg.lmRob.brestTotalP$scale-model1$tg.lmRob.brest$scale)/model1$tg.lmRob.brestTotalP$scale
#[1] 0.254534
#> model1$tg.lmRob.brest$scale/model1$tg.lmRob.brestTotalP$scale
#[1] 0.745466


##########
# Newlyn #
##########

# What we are going to do now is calculate the sd of the residuals as each of the
# 48 components is added in turn 
model1$tg.sd.newlyn <- vector(mode="numeric", length=48)

for(i in 1:48){
  model1$data.newlyn <- data.frame(msl=model1$newlynTotalP, t=seq(from=1960, to=1999), Pa=model1$Pa[,i])
  model1$tg.lmRob.newlyn <- lmRob(msl ~ t + Pa, data=model1$data.newlyn)
  
#  model1$tg.resid.newlyn <- model1$tg.lmRob.newlyn$resid
  model1$tg.sd.newlyn[i] <- model1$tg.lmRob.newlyn$scale
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model1$tg.sd.sorted.newlyn <- sort(model1$tg.sd.newlyn, index=TRUE)

# Now construct models of the top 9 components
model1$data.newlyn <- data.frame(msl=model1$newlynTotalP,
  t=seq(from=1960,to=1999),
  Pa1=model1$Pa[,model1$tg.sd.sorted.newlyn$ix[1]], 
  Pa2=model1$Pa[,model1$tg.sd.sorted.newlyn$ix[2]], Pa3=model1$Pa[,model1$tg.sd.sorted.newlyn$ix[3]],
  Pa4=model1$Pa[,model1$tg.sd.sorted.newlyn$ix[4]], Pa5=model1$Pa[,model1$tg.sd.sorted.newlyn$ix[5]],
  Pa6=model1$Pa[,model1$tg.sd.sorted.newlyn$ix[6]], Pa7=model1$Pa[,model1$tg.sd.sorted.newlyn$ix[7]],
  Pa8=model1$Pa[,model1$tg.sd.sorted.newlyn$ix[8]], Pa9=model1$Pa[,model1$tg.sd.sorted.newlyn$ix[9]])

model1$tg.lmRob.newlyn <- lmRob(msl ~ t + Pa1 + Pa2 + Pa3 + Pa4 + Pa5 + Pa6 + Pa7 + Pa8 + Pa9, x=T, y=T,
                               data=model1$data.newlyn)

# Variance reduction
#> (model1$tg.lmRob.newlynTotalP$scale-model1$tg.lmRob.newlyn$scale)/model1$tg.lmRob.newlynTotalP$scale
#[1] 0.3708184
#> model1$tg.lmRob.newlyn$scale/model1$tg.lmRob.newlynTotalP$scale
#[1] 0.6291816

#############
## Model 2 ##
#############
lmRobControl <- lmRob.control(mxr=100,mxf=100,trace=F)

# Follow Thompson (1980) and see if we can improve
# Uses local pressures from Plymouth
model2<-new.env()

# Just use annual pressures over the periods of the model data (Pa) ending in 2000 (end of model pressure) 
model2$newlynPa <-
  met$hadSLP2rAnnualPressureArray[intersect(which(met$hadSLP2rAnnualTime>=1960), which(met$hadSLP2rAnnualTime<=1999)),1]
model2$brestPa <-
  met$hadSLP2rAnnualPressureArray[intersect(which(met$hadSLP2rAnnualTime>=1960), which(met$hadSLP2rAnnualTime<=1999)),2]
model2$time <- c(1960:1999)

# Calculate the total pressure
model2$brestTotalP <- 1025*9.8*tg_annual$brest6099[,2]/1000 + model2$brestPa
model2$newlynTotalP <- 1025*9.8*tg_annual$newlyn6099[,2]/1000 + model2$newlynPa

model2$tg.lmRob.brestTotalP <- lmRob(model2$brestTotalP ~ model2$time, control=lmRobControl)
model2$tg.lmRob.newlynTotalP <- lmRob(model2$newlynTotalP ~ model2$time, control=lmRobControl)

# Some of the HadSLP2r pressure measurements are duplicated due to the large (5 deg) grid size
met$hadSLP2rUniq <- which(!duplicated(met$hadSLP2rAnnualPressureArray[1,]))
# First two points are Newlyn and Brest, which are the same
met$hadNUniq <- length(met$hadSLP2rUniq)
met$hadSLP2rUniq <- met$hadSLP2rUniq[2:met$hadNUniq]
met$hadNUniq <- length(met$hadSLP2rUniq)

model2$Pa <- array(NA,dim=c(40,(met$hadNUniq*3)))
# Zero lag
for(i in 1:met$hadNUniq){
  model2$Pa[,i] <- as.vector(met$hadSLP2rAnnualPressureArray[111:150,met$hadSLP2rUniq[i]])
}
# Lag 1 year
for(i in (met$hadNUniq+1):(met$hadNUniq*2)){
  model2$Pa[,i] <- as.vector(met$hadSLP2rAnnualPressureArray[110:149,met$hadSLP2rUniq[(i-met$hadNUniq)]])
}
# Lag 2 years....
for(i in ((met$hadNUniq*2)+1):(met$hadNUniq*3)){
  model2$Pa[,i] <- as.vector(met$hadSLP2rAnnualPressureArray[109:148,met$hadSLP2rUniq[(i-(met$hadNUniq*2))]])
}

#########
# Brest #
#########

# What we are going to do now is calculate the sd of the residuals as each of the
# 102 components is added in turn 
model2$tg.sd.brest <- vector(mode="numeric", length=(met$hadNUniq*3))

for(i in 1:(met$hadNUniq*3)){
  model2$data.brest <- data.frame(msl=model2$brestTotalP, t=seq(from=1960,to=1999), Pa=model2$Pa[,i])
  model2$tg.lmRob.brest <- lmRob(msl ~ t + Pa, data=model2$data.brest, control=lmRobControl)
  
#  model2$tg.resid.brest <- model2$tg.lmRob.brest$resid
  model2$tg.sd.brest[i] <- model2$tg.lmRob.brest$scale
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model2$tg.sd.sorted.brest <- sort(model2$tg.sd.brest, index=TRUE)

# Now construct models of the top 9 components
model2$data.brest <- data.frame(msl=model2$brestTotalP,
  t=seq(from=1960,to=1999),
  Pa1=model2$Pa[,model2$tg.sd.sorted.brest$ix[1]], 
  Pa2=model2$Pa[,model2$tg.sd.sorted.brest$ix[2]], Pa3=model2$Pa[,model2$tg.sd.sorted.brest$ix[3]],
  Pa4=model2$Pa[,model2$tg.sd.sorted.brest$ix[4]], Pa5=model2$Pa[,model2$tg.sd.sorted.brest$ix[5]],
  Pa6=model2$Pa[,model2$tg.sd.sorted.brest$ix[6]], Pa7=model2$Pa[,model2$tg.sd.sorted.brest$ix[7]],
  Pa8=model2$Pa[,model2$tg.sd.sorted.brest$ix[8]], Pa9=model2$Pa[,model2$tg.sd.sorted.brest$ix[9]])

model2$tg.lmRob.brest <- lmRob(msl ~ t + Pa1 + Pa2 + Pa3 + Pa4 + Pa5 + Pa6 + Pa7 + Pa8 + Pa9, x=T, y=T,
                               data=model2$data.brest, control=lmRobControl)

# Variance reduction
#> (model2$tg.lmRob.brest$r.sq-model2$tg.lmRob.brestTotalP$r.sq)
#[1] 0.3193968
#> model2$tg.lmRob.brest$r.sq
#[1] 0.606882

##########
# Newlyn #
##########

# What we are going to do now is calculate the sd of the residuals as each of the
# components is added in turn 
model2$tg.sd.newlyn <- vector(mode="numeric", length=(met$hadNUniq*3))

for(i in 1:(met$hadNUniq*3)){
  model2$data.newlyn <- data.frame(msl=model2$newlynTotalP, t=seq(from=1960, to=1999), Pa=model2$Pa[,i])
  model2$tg.lmRob.newlyn <- lmRob(msl ~ t + Pa, data=model2$data.newlyn, control=lmRobControl)
  
#  model2$tg.resid.newlyn <- model2$tg.lmRob.newlyn$resid
  model2$tg.sd.newlyn[i] <- model2$tg.lmRob.newlyn$scale
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model2$tg.sd.sorted.newlyn <- sort(model2$tg.sd.newlyn, index=TRUE)

# Now construct models of the top 9 components
model2$data.newlyn <- data.frame(msl=model2$newlynTotalP,
  t=seq(from=1960,to=1999),
  Pa1=model2$Pa[,model2$tg.sd.sorted.newlyn$ix[1]], 
  Pa2=model2$Pa[,model2$tg.sd.sorted.newlyn$ix[2]], Pa3=model2$Pa[,model2$tg.sd.sorted.newlyn$ix[3]],
  Pa4=model2$Pa[,model2$tg.sd.sorted.newlyn$ix[4]], Pa5=model2$Pa[,model2$tg.sd.sorted.newlyn$ix[5]],
  Pa6=model2$Pa[,model2$tg.sd.sorted.newlyn$ix[6]], Pa7=model2$Pa[,model2$tg.sd.sorted.newlyn$ix[7]],
  Pa8=model2$Pa[,model2$tg.sd.sorted.newlyn$ix[8]], Pa9=model2$Pa[,model2$tg.sd.sorted.newlyn$ix[9]])

model2$tg.lmRob.newlyn <- lmRob(msl ~ t + Pa1 + Pa2 + Pa3 + Pa4 + Pa5 + Pa6 + Pa7 + Pa8 + Pa9, x=T, y=T,
                               data=model2$data.newlyn, control=lmRobControl)

# Variance reduction
#> (model2$tg.lmRob.newlyn$r.sq-model2$tg.lmRob.newlynTotalP$r.sq)
#[1] 0.3942693
#> model2$tg.lmRob.newlyn$r.sq
#[1] 0.6150269

#############
## Model 3 ##
#############
# Follow Thompson (1980) and see if we can improve
# Uses local pressures from Plymouth and HadSLP2r
# Calculate over 1953-2008 period

model3<-new.env()

# Just use annual pressures over the periods of the model data (Pa) ending in 2000 (end of model pressure) 
model3$newlynPa <-
  met$hadSLP2rAnnualPressureArray[intersect(which(met$hadSLP2rAnnualTime>=1953), which(met$hadSLP2rAnnualTime<=2008)),1]
model3$brestPa <-
  met$hadSLP2rAnnualPressureArray[intersect(which(met$hadSLP2rAnnualTime>=1953), which(met$hadSLP2rAnnualTime<=2008)),2]

model3$time <- c(1953:2008)
tg_annual$brest5308 <- tg_annual$brest[147:202,]
tg_annual$newlyn5308 <- tg_annual$newlyn[38:93,]

# Calculate the total pressure
model3$brestTotalP <- 1025*9.8*tg_annual$brest5308[,2]/1000 + model3$brestPa
model3$newlynTotalP <- 1025*9.8*tg_annual$newlyn5308[,2]/1000 + model3$newlynPa

model3$tg.lmRob.brestTotalP <- lmRob(model3$brestTotalP ~ model3$time, control=lmRobControl)
model3$tg.lmRob.newlynTotalP <- lmRob(model3$newlynTotalP ~ model3$time, control=lmRobControl)

model3$Pa <- array(NA,dim=c(56,(met$hadNUniq*3)))
# Zero lag
for(i in 1:met$hadNUniq){
  model3$Pa[,i] <- as.vector(met$hadSLP2rAnnualPressureArray[104:159,met$hadSLP2rUniq[i]])
}
# Lag 1 year
for(i in (met$hadNUniq+1):(met$hadNUniq*2)){
  model3$Pa[,i] <- as.vector(met$hadSLP2rAnnualPressureArray[103:158,met$hadSLP2rUniq[(i-met$hadNUniq)]])
}
# Lag 2 years....
for(i in ((met$hadNUniq*2)+1):(met$hadNUniq*3)){
  model3$Pa[,i] <- as.vector(met$hadSLP2rAnnualPressureArray[102:157,met$hadSLP2rUniq[(i-(met$hadNUniq*2))]])
}

#########
# Brest #
#########

# What we are going to do now is calculate the sd of the residuals as each of the
# components is added in turn 
model3$tg.sd.brest <- vector(mode="numeric", length=(met$hadNUniq*3))

for(i in 1:(met$hadNUniq*3)){
  model3$data.brest <- data.frame(msl=model3$brestTotalP, t=seq(from=1953,to=2008), Pa=model3$Pa[,i])
  model3$tg.lmRob.brest <- lmRob(msl ~ t + Pa, data=model3$data.brest, control=lmRobControl)

  model3$tg.sd.brest[i] <- model3$tg.lmRob.brest$scale
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model3$tg.sd.sorted.brest <- sort(model3$tg.sd.brest, index=TRUE)

# Now construct models of the top 9 components
model3$data.brest <- data.frame(msl=model3$brestTotalP,
  t=seq(from=1953,to=2008),
  Pa1=model3$Pa[,model3$tg.sd.sorted.brest$ix[1]], 
  Pa2=model3$Pa[,model3$tg.sd.sorted.brest$ix[2]], Pa3=model3$Pa[,model3$tg.sd.sorted.brest$ix[3]],
  Pa4=model3$Pa[,model3$tg.sd.sorted.brest$ix[4]], Pa5=model3$Pa[,model3$tg.sd.sorted.brest$ix[5]],
  Pa6=model3$Pa[,model3$tg.sd.sorted.brest$ix[6]], Pa7=model3$Pa[,model3$tg.sd.sorted.brest$ix[7]],
  Pa8=model3$Pa[,model3$tg.sd.sorted.brest$ix[8]], Pa9=model3$Pa[,model3$tg.sd.sorted.brest$ix[9]])

model3$tg.lmRob.brest <- lmRob(msl ~ t + Pa1 + Pa2 + Pa3 + Pa4 + Pa5 + Pa6 + Pa7 + Pa8 + Pa9, x=T, y=T,
                               data=model3$data.brest, control=lmRobControl)

# Variance reduction
#> (model3$tg.lmRob.brest$r.sq-model3$tg.lmRob.brestTotalP$r.sq)
#[1] 0.2491007
#> model3$tg.lmRob.brest$r.sq
#[1] 0.7240837

##########
# Newlyn #
##########

# What we are going to do now is calculate the sd of the residuals as each of the
# components is added in turn 
model3$tg.sd.newlyn <- vector(mode="numeric", length=(met$hadNUniq*3))

for(i in 1:(met$hadNUniq*3)){
  model3$data.newlyn <- data.frame(msl=model3$newlynTotalP, t=seq(from=1953, to=2008), Pa=model3$Pa[,i])
  model3$tg.lmRob.newlyn <- lmRob(msl ~ t + Pa, data=model3$data.newlyn, control=lmRobControl)
  
  model3$tg.sd.newlyn[i] <- model3$tg.lmRob.newlyn$scale
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model3$tg.sd.sorted.newlyn <- sort(model3$tg.sd.newlyn, index=TRUE)

# Now construct models of the top 9 components
model3$data.newlyn <- data.frame(msl=model3$newlynTotalP,
  t=seq(from=1953,to=2008),
  Pa1=model3$Pa[,model3$tg.sd.sorted.newlyn$ix[1]], 
  Pa2=model3$Pa[,model3$tg.sd.sorted.newlyn$ix[2]], Pa3=model3$Pa[,model3$tg.sd.sorted.newlyn$ix[3]],
  Pa4=model3$Pa[,model3$tg.sd.sorted.newlyn$ix[4]], Pa5=model3$Pa[,model3$tg.sd.sorted.newlyn$ix[5]],
  Pa6=model3$Pa[,model3$tg.sd.sorted.newlyn$ix[6]], Pa7=model3$Pa[,model3$tg.sd.sorted.newlyn$ix[7]],
  Pa8=model3$Pa[,model3$tg.sd.sorted.newlyn$ix[8]], Pa9=model3$Pa[,model3$tg.sd.sorted.newlyn$ix[9]])

model3$tg.lmRob.newlyn <- lmRob(msl ~ t + Pa1 + Pa2 + Pa3 + Pa4 + Pa5 + Pa6 + Pa7 + Pa8 + Pa9, x=T, y=T,
                               data=model3$data.newlyn, control=lmRobControl)

# Variance reduction
#> (model3$tg.lmRob.newlyn$r.sq-model3$tg.lmRob.newlynTotalP$r.sq)
#[1] 0.273999
#> model3$tg.lmRob.newlyn$r.sq
#[1] 0.7001414
#> c(coef(model3$tg.lmRob.brest)[2]/10.045, coef(model3$tg.lmRob.newlyn)[2]/10.045)
#       t        t 
#1.822623 1.798708 
#> c(sqrt(diag(model3$tg.lmRob.brest$cov))[2]/10.045,  sqrt(diag(model3$tg.lmRob.newlyn$cov))[2]/10.045)
#        t         t 
#0.2615218 0.2927652

#############
## Model 4 ##
#############
# Follow Thompson (1980) and see if we can improve
# Uses pressures from HadSLP2r
# Calculate over 1916-1943 period

model4<-new.env()

# Just use annual pressures over the periods of the model data (Pa) ending in 2000 (end of model pressure) 
model4$newlynPa <-
  met$hadSLP2rAnnualPressureArray[intersect(which(met$hadSLP2rAnnualTime>=1916), which(met$hadSLP2rAnnualTime<=1943)),1]
model4$brestPa <-
  met$hadSLP2rAnnualPressureArray[intersect(which(met$hadSLP2rAnnualTime>=1916), which(met$hadSLP2rAnnualTime<=1943)),2]
model4$time <- c(1916:1943)
tg_annual$brest1643 <- tg_annual$brest[110:137,]
tg_annual$newlyn1643 <- tg_annual$newlyn[1:28,]

# Calculate the total pressure
model4$brestTotalP <- 1025*9.8*tg_annual$brest1643[,2]/1000 + model4$brestPa
model4$newlynTotalP <- 1025*9.8*tg_annual$newlyn1643[,2]/1000 + model4$newlynPa

model4$tg.lmRob.brestTotalP <- lmRob(model4$brestTotalP ~ model4$time, control=lmRobControl)
model4$tg.lmRob.newlynTotalP <- lmRob(model4$newlynTotalP ~ model4$time, x=T, control=lmRobControl)

model4$Pa <- array(NA,dim=c(28,(met$hadNUniq*3)))
# Zero lag
for(i in 1:met$hadNUniq){
  model4$Pa[,i] <- as.vector(met$hadSLP2rAnnualPressureArray[67:94,met$hadSLP2rUniq[i]])
}
# Lag 1 year
for(i in (met$hadNUniq+1):(met$hadNUniq*2)){
  model4$Pa[,i] <- as.vector(met$hadSLP2rAnnualPressureArray[66:93,met$hadSLP2rUniq[(i-met$hadNUniq)]])
}
# Lag 2 years....
for(i in ((met$hadNUniq*2)+1):(met$hadNUniq*3)){
  model4$Pa[,i] <- as.vector(met$hadSLP2rAnnualPressureArray[65:92,met$hadSLP2rUniq[(i-(met$hadNUniq*2))]])
}

#########
# Brest #
#########

# What we are going to do now is calculate the sd of the residuals as each of the
# components is added in turn 
model4$tg.sd.brest <- vector(mode="numeric", length=(met$hadNUniq*3))

for(i in 1:(met$hadNUniq*3)){
  model4$data.brest <- data.frame(msl=model4$brestTotalP, t=seq(from=1916,to=1943), Pa=model4$Pa[,i])
  model4$tg.lmRob.brest <- lmRob(msl ~ t + Pa, data=model4$data.brest, control=lmRobControl)
  
  model4$tg.sd.brest[i] <- model4$tg.lmRob.brest$scale
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model4$tg.sd.sorted.brest <- sort(model4$tg.sd.brest, index=TRUE)

# Now construct models of the top 9 components
model4$data.brest <- data.frame(msl=model4$brestTotalP,
  t=seq(from=1916,to=1943),
  Pa1=model4$Pa[,model4$tg.sd.sorted.brest$ix[1]], 
  Pa2=model4$Pa[,model4$tg.sd.sorted.brest$ix[2]], Pa3=model4$Pa[,model4$tg.sd.sorted.brest$ix[3]],
  Pa4=model4$Pa[,model4$tg.sd.sorted.brest$ix[4]], Pa5=model4$Pa[,model4$tg.sd.sorted.brest$ix[5]],
  Pa6=model4$Pa[,model4$tg.sd.sorted.brest$ix[6]], Pa7=model4$Pa[,model4$tg.sd.sorted.brest$ix[7]],
  Pa8=model4$Pa[,model4$tg.sd.sorted.brest$ix[8]], Pa9=model4$Pa[,model4$tg.sd.sorted.brest$ix[9]])

model4$tg.lmRob.brest <- lmRob(msl ~ t + Pa1 + Pa2 + Pa3 + Pa4 + Pa5 + Pa6 + Pa7 + Pa8 + Pa9, x=T, y=T,
                               data=model4$data.brest, control=lmRobControl)

# Variance reduction
#> (model4$tg.lmRob.brest$r.sq-model4$tg.lmRob.brestTotalP$r.sq)
#[1] 0.3982849
#> model4$tg.lmRob.brest$r.sq
#[1] 0.6360687

##########
# Newlyn #
##########

# What we are going to do now is calculate the sd of the residuals as each of the
# components is added in turn 
model4$tg.sd.newlyn <- vector(mode="numeric", length=(met$hadNUniq*3))

for(i in 1:(met$hadNUniq*3)){
  model4$data.newlyn <- data.frame(msl=model4$newlynTotalP, t=seq(from=1916, to=1943), Pa=model4$Pa[,i])
  model4$tg.lmRob.newlyn <- lmRob(msl ~ t + Pa, data=model4$data.newlyn, control=lmRobControl)
  
  model4$tg.sd.newlyn[i] <- model4$tg.lmRob.newlyn$scale
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model4$tg.sd.sorted.newlyn <- sort(model4$tg.sd.newlyn, index=TRUE)

# Now construct models of the top 9 components
model4$data.newlyn <- data.frame(msl=model4$newlynTotalP,
  t=seq(from=1916,to=1943),
  Pa1=model4$Pa[,model4$tg.sd.sorted.newlyn$ix[1]], 
  Pa2=model4$Pa[,model4$tg.sd.sorted.newlyn$ix[2]], Pa3=model4$Pa[,model4$tg.sd.sorted.newlyn$ix[3]],
  Pa4=model4$Pa[,model4$tg.sd.sorted.newlyn$ix[4]], Pa5=model4$Pa[,model4$tg.sd.sorted.newlyn$ix[5]],
  Pa6=model4$Pa[,model4$tg.sd.sorted.newlyn$ix[6]], Pa7=model4$Pa[,model4$tg.sd.sorted.newlyn$ix[7]],
  Pa8=model4$Pa[,model4$tg.sd.sorted.newlyn$ix[8]], Pa9=model4$Pa[,model4$tg.sd.sorted.newlyn$ix[9]])

model4$tg.lmRob.newlyn <- lmRob(msl ~ t + Pa1 + Pa2 + Pa3 + Pa4 + Pa5 + Pa6 + Pa7 + Pa8 + Pa9, x=T, y=T,
                               data=model4$data.newlyn, control=lmRobControl)

# Variance reduction
#> (model4$tg.lmRob.newlyn$r.sq-model4$tg.lmRob.newlynTotalP$r.sq)
#[1] 0.3638209
#> model4$tg.lmRob.newlyn$r.sq
#[1] 0.7181428
#> c(coef(model4$tg.lmRob.brest)[2]/10.045, coef(model4$tg.lmRob.newlyn)[2]/10.045)
#       t        t 
#2.079806 1.544290
#> c(sqrt(diag(model4$tg.lmRob.brest$cov))[2]/10.045,  sqrt(diag(model4$tg.lmRob.newlyn$cov))[2]/10.045)
#        t         t 
#0.7195811 0.6130437

#############
## Model 5 ##
#############
# Follow Thompson (1980) and see if we can improve
# Uses pressures from HadSLP2r
# Calculate over 1916-2008 period

model5<-new.env()

# Just use annual pressures over the periods of the model data (Pa) ending in 2000 (end of model pressure) 
model5$newlynPa <-
  met$hadSLP2rAnnualPressureArray[intersect(which(met$hadSLP2rAnnualTime>=1916), which(met$hadSLP2rAnnualTime<=2008)),1]
model5$brestPa <-
  met$hadSLP2rAnnualPressureArray[intersect(which(met$hadSLP2rAnnualTime>=1916), which(met$hadSLP2rAnnualTime<=2008)),2]

model5$time <- c(1916:2008)
tg_annual$brest1608 <- tg_annual$brest[110:202,]
# Add 22mm to Brest post-war to see what effect this has
tg_annual$brest1608_22 <- tg_annual$brest1608
tg_annual$brest1608_22[38:93,2] <- tg_annual$brest1608_22[38:93,2]+22

tg_annual$newlyn1608 <- tg_annual$newlyn[1:93,]

# Calculate the total pressure
model5$brestTotalP <- 1025*9.8*tg_annual$brest1608[,2]/1000 + model5$brestPa
model5$brestTotalP_22<- 1025*9.8*tg_annual$brest1608_22[,2]/1000 + model5$brestPa
model5$newlynTotalP <- 1025*9.8*tg_annual$newlyn1608[,2]/1000 + model5$newlynPa

model5$tg.lmRob.brestTotalP <- lmRob(model5$brestTotalP ~ model5$time, control=lmRobControl)
model5$tg.lmRob.brestTotalP_22<- lmRob(model5$brestTotalP_22 ~ model5$time, control=lmRobControl)
model5$tg.lmRob.newlynTotalP <- lmRob(model5$newlynTotalP ~ model5$time, x=T, control=lmRobControl)

model5$Pa <- array(NA,dim=c(93,(met$hadNUniq*3)))
# Zero lag
for(i in 1:met$hadNUniq){
  model5$Pa[,i] <- as.vector(met$hadSLP2rAnnualPressureArray[67:159,met$hadSLP2rUniq[i]])
}
# Lag 1 year
for(i in (met$hadNUniq+1):(met$hadNUniq*2)){
  model5$Pa[,i] <- as.vector(met$hadSLP2rAnnualPressureArray[66:158,met$hadSLP2rUniq[(i-met$hadNUniq)]])
}
# Lag 2 years....
for(i in ((met$hadNUniq*2)+1):(met$hadNUniq*3)){
  model5$Pa[,i] <- as.vector(met$hadSLP2rAnnualPressureArray[65:157,met$hadSLP2rUniq[(i-(met$hadNUniq*2))]])
}

#########
# Brest #
#########

# What we are going to do now is calculate the sd of the residuals as each of the
# components is added in turn 
model5$tg.sd.brest <- vector(mode="numeric", length=(met$hadNUniq*3))

for(i in 1:(met$hadNUniq*3)){
  model5$data.brest <- data.frame(msl=model5$brestTotalP, t=seq(from=1916,to=2008), Pa=model5$Pa[,i])
  model5$tg.lmRob.brest <- lmRob(msl ~ t + Pa, data=model5$data.brest, control=lmRobControl)
  
  model5$tg.sd.brest[i] <- model5$tg.lmRob.brest$scale
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model5$tg.sd.sorted.brest <- sort(model5$tg.sd.brest, index=TRUE)

# Now construct models of the top 9 components
model5$data.brest <- data.frame(msl=model5$brestTotalP,
  t=seq(from=1916,to=2008),
  Pa1=model5$Pa[,model5$tg.sd.sorted.brest$ix[1]], 
  Pa2=model5$Pa[,model5$tg.sd.sorted.brest$ix[2]], Pa3=model5$Pa[,model5$tg.sd.sorted.brest$ix[3]],
  Pa4=model5$Pa[,model5$tg.sd.sorted.brest$ix[4]], Pa5=model5$Pa[,model5$tg.sd.sorted.brest$ix[5]],
  Pa6=model5$Pa[,model5$tg.sd.sorted.brest$ix[6]], Pa7=model5$Pa[,model5$tg.sd.sorted.brest$ix[7]],
  Pa8=model5$Pa[,model5$tg.sd.sorted.brest$ix[8]], Pa9=model5$Pa[,model5$tg.sd.sorted.brest$ix[9]])

model5$tg.lmRob.brest <- lmRob(msl ~ t + Pa1 + Pa2 + Pa3 + Pa4 + Pa5 + Pa6 + Pa7 + Pa8 + Pa9, x=T, y=T,
                               data=model5$data.brest, control=lmRobControl)
model5$tg.lmRob.brest.fitted.lm <- lmRob(model5$tg.lmRob.brest$fitted  ~ model5$tg.lmRob.brest$x[,2], control=lmRobControl)

# Variance reduction
#> (model5$tg.lmRob.brest$r.sq-model5$tg.lmRob.brestTotalP$r.sq)
#[1] 0.1072223
#> model5$tg.lmRob.brest$r.sq
#[1] 0.6555045

#################
# Brest plus 22 #
#################

# What we are going to do now is calculate the sd of the residuals as each of the
# components is added in turn 
model5$tg.sd.brest_22 <- vector(mode="numeric", length=(met$hadNUniq*3))

for(i in 1:(met$hadNUniq*3)){
  model5$data.brest <- data.frame(msl=model5$brestTotalP_22, t=seq(from=1916,to=2008), Pa=model5$Pa[,i])
  model5$tg.lmRob.brest_22 <- lmRob(msl ~ t + Pa, data=model5$data.brest, control=lmRobControl)
  
  model5$tg.sd.brest_22[i] <- model5$tg.lmRob.brest_22$scale
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model5$tg.sd.sorted.brest_22<- sort(model5$tg.sd.brest_22, index=TRUE)

# Now construct models of the top 9 components
model5$data.brest <- data.frame(msl=model5$brestTotalP_22,
  t=seq(from=1916,to=2008),
  Pa1=model5$Pa[,model5$tg.sd.sorted.brest_22$ix[1]], 
  Pa2=model5$Pa[,model5$tg.sd.sorted.brest_22$ix[2]], Pa3=model5$Pa[,model5$tg.sd.sorted.brest_22$ix[3]],
  Pa4=model5$Pa[,model5$tg.sd.sorted.brest_22$ix[4]], Pa5=model5$Pa[,model5$tg.sd.sorted.brest_22$ix[5]],
  Pa6=model5$Pa[,model5$tg.sd.sorted.brest_22$ix[6]], Pa7=model5$Pa[,model5$tg.sd.sorted.brest_22$ix[7]],
  Pa8=model5$Pa[,model5$tg.sd.sorted.brest_22$ix[8]], Pa9=model5$Pa[,model5$tg.sd.sorted.brest_22$ix[9]])

model5$tg.lmRob.brest_22 <- lmRob(msl ~ t + Pa1 + Pa2 + Pa3 + Pa4 + Pa5 + Pa6 + Pa7 + Pa8 + Pa9, x=T, y=T,
                               data=model5$data.brest, control=lmRobControl)

# Variance reduction
#> (model5$tg.lmRob.brest_22$r.sq-model5$tg.lmRob.brestTotalP_22$r.sq)
#[1] 0.1072223
#> model5$tg.lmRob.brest_22$r.sq
#[1] 0.6555045

##########
# Newlyn #
##########

# What we are going to do now is calculate the sd of the residuals as each of the
# components is added in turn 
model5$tg.sd.newlyn <- vector(mode="numeric", length=(met$hadNUniq*3))

for(i in 1:(met$hadNUniq*3)){
  model5$data.newlyn <- data.frame(msl=model5$newlynTotalP, t=seq(from=1916, to=2008), Pa=model5$Pa[,i])
  model5$tg.lmRob.newlyn <- lmRob(msl ~ t + Pa, data=model5$data.newlyn, control=lmRobControl)
  
  model5$tg.sd.newlyn[i] <- model5$tg.lmRob.newlyn$scale
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model5$tg.sd.sorted.newlyn <- sort(model5$tg.sd.newlyn, index=TRUE)

# Now construct models of the top 9 components
model5$data.newlyn <- data.frame(msl=model5$newlynTotalP,
  t=seq(from=1916,to=2008),
  Pa1=model5$Pa[,model5$tg.sd.sorted.newlyn$ix[1]], 
  Pa2=model5$Pa[,model5$tg.sd.sorted.newlyn$ix[2]], Pa3=model5$Pa[,model5$tg.sd.sorted.newlyn$ix[3]],
  Pa4=model5$Pa[,model5$tg.sd.sorted.newlyn$ix[4]], Pa5=model5$Pa[,model5$tg.sd.sorted.newlyn$ix[5]],
  Pa6=model5$Pa[,model5$tg.sd.sorted.newlyn$ix[6]], Pa7=model5$Pa[,model5$tg.sd.sorted.newlyn$ix[7]],
  Pa8=model5$Pa[,model5$tg.sd.sorted.newlyn$ix[8]], Pa9=model5$Pa[,model5$tg.sd.sorted.newlyn$ix[9]])

model5$tg.lmRob.newlyn <- lmRob(msl ~ t + Pa1 + Pa2 + Pa3 + Pa4 + Pa5 + Pa6 + Pa7 + Pa8 + Pa9, x=T, y=T,
                               data=model5$data.newlyn, control=lmRobControl)
model5$tg.lmRob.newlyn.fitted.lm <- lmRob(model5$tg.lmRob.newlyn$fitted  ~ model5$tg.lmRob.newlyn$x[,2], control=lmRobControl)
# Variance reduction
#> (model5$tg.lmRob.newlyn$r.sq-model5$tg.lmRob.newlynTotalP$r.sq)
#[1] 0.1392065
#> model5$tg.lmRob.newlyn$r.sq
#[1] 0.7895451
#c(coef(model5$tg.lmRob.brest)[2]/10.045, coef(model5$tg.lmRob.newlyn)[2]/10.045)
#       t        t 
#1.471297 1.816946
#c(sqrt(diag(model5$tg.lmRob.brest$cov))[2]/10.045,  sqrt(diag(model5$tg.lmRob.newlyn$cov))[2]/10.045)
#        t         t 
#0.1200115 0.1068158

#############
## Model 6 ##
#############
# Follow Thompson (1980) and see if we can improve
# Uses pressures from HadSLP2r
# Calculate over 1850-1943 period, just for Brest ;)

model6<-new.env()

# Just use annual pressures over the periods of the model data (Pa) ending in 2000 (end of model pressure) 
model6$brestPa <-
  met$hadSLP2rAnnualPressureArray[intersect(which(met$hadSLP2rAnnualTime>=1850), which(met$hadSLP2rAnnualTime<=1943)),2]
model6$time <- c(1850:1943)
tg_annual$brest5043 <- tg_annual$brest[44:137,]

# Calculate the total pressure
model6$brestTotalP <- 1025*9.8*tg_annual$brest5043[,2]/1000 + model6$brestPa

model6$tg.lmRob.brestTotalP <- lmRob(model6$brestTotalP ~ model6$time, control=lmRobControl)

model6$Pa <- array(NA,dim=c(94,(met$hadNUniq*3)))
# Zero lag
for(i in 1:met$hadNUniq){
  model6$Pa[,i] <- as.vector(met$hadSLP2rAnnualPressureArray[1:94,met$hadSLP2rUniq[i]])
}
# Lag 1 year
for(i in (met$hadNUniq+1):(met$hadNUniq*2)){
  model6$Pa[2:94,i] <- as.vector(met$hadSLP2rAnnualPressureArray[1:93,met$hadSLP2rUniq[(i-met$hadNUniq)]])
  model6$Pa[1,i] <- mean(as.vector(met$hadSLP2rAnnualPressureArray[1:94,met$hadSLP2rUniq[(i-(met$hadNUniq*2))]]), na.rm=T)
}
# Lag 2 years....
for(i in ((met$hadNUniq*2)+1):(met$hadNUniq*3)){
  model6$Pa[3:94,i] <- as.vector(met$hadSLP2rAnnualPressureArray[1:92,met$hadSLP2rUniq[(i-(met$hadNUniq*2))]])
  model6$Pa[1:2,i] <- mean(as.vector(met$hadSLP2rAnnualPressureArray[1:94,met$hadSLP2rUniq[(i-(met$hadNUniq*2))]]), na.rm=T)
}

#########
# Brest #
#########

# What we are going to do now is calculate the sd of the residuals as each of the
# components is added in turn 
model6$tg.sd.brest <- vector(mode="numeric", length=(met$hadNUniq*3))

for(i in 1:(met$hadNUniq*3)){
  model6$data.brest <- data.frame(msl=model6$brestTotalP, t=seq(from=1850,to=1943), Pa=model6$Pa[,i])
  model6$tg.lmRob.brest <- lmRob(msl ~ t + Pa, data=model6$data.brest, control=lmRobControl)
  
  model6$tg.sd.brest[i] <- model6$tg.lmRob.brest$scale
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model6$tg.sd.sorted.brest <- sort(model6$tg.sd.brest, index=TRUE)

# Now construct models of the top 9 components
model6$data.brest <- data.frame(msl=model6$brestTotalP,
  t=seq(from=1850,to=1943),
  Pa1=model6$Pa[,model6$tg.sd.sorted.brest$ix[1]], 
  Pa2=model6$Pa[,model6$tg.sd.sorted.brest$ix[2]], Pa3=model6$Pa[,model6$tg.sd.sorted.brest$ix[3]],
  Pa4=model6$Pa[,model6$tg.sd.sorted.brest$ix[4]], Pa5=model6$Pa[,model6$tg.sd.sorted.brest$ix[5]],
  Pa6=model6$Pa[,model6$tg.sd.sorted.brest$ix[6]], Pa7=model6$Pa[,model6$tg.sd.sorted.brest$ix[7]],
  Pa8=model6$Pa[,model6$tg.sd.sorted.brest$ix[8]], Pa9=model6$Pa[,model6$tg.sd.sorted.brest$ix[9]])

model6$tg.lmRob.brest <- lmRob(msl ~ t + Pa1 + Pa2 + Pa3 + Pa4 + Pa5 + Pa6 + Pa7 + Pa8 + Pa9, x=T, y=T,
                               data=model6$data.brest, control=lmRobControl)

# Variance reduction
#> (model6$tg.lmRob.brest$r.sq-model6$tg.lmRob.brestTotalP$r.sq)
#[1] 0.1658905
#> model6$tg.lmRob.brest$r.sq
#[1] 0.5441145
#>sqrt(diag(model6$tg.lmRob.brest$cov))[2]/10.045
#        t 
#0.2692527 
#> coef(model6$tg.lmRob.brest)[2]/10.045
#       t 
#1.126984

#############
## Model 7 ##
#############
# Follow Thompson (1980) and see if we can improve
# Uses pressures from HadSLP2r
# Calculate over 1888-1943 period, just for Brest - 56 years same as 1953-2008

model7<-new.env()

# Just use annual pressures over the periods of the model data (Pa) ending in 2000 (end of model pressure) 
model7$brestPa <-
  met$hadSLP2rAnnualPressureArray[intersect(which(met$hadSLP2rAnnualTime>=1888), which(met$hadSLP2rAnnualTime<=1943)),2]
model7$time <- c(1888:1943)
tg_annual$brest8843 <- tg_annual$brest[82:137,]

# Calculate the total pressure
model7$brestTotalP <- 1025*9.8*tg_annual$brest8843[,2]/1000 + model7$brestPa

model7$tg.lmRob.brestTotalP <- lmRob(model7$brestTotalP ~ model7$time, control=lmRobControl)

model7$Pa <- array(NA,dim=c(56,(met$hadNUniq*3)))
# Zero lag
for(i in 1:met$hadNUniq){
  model7$Pa[,i] <- as.vector(met$hadSLP2rAnnualPressureArray[39:94,met$hadSLP2rUniq[i]])
}
# Lag 1 year
for(i in (met$hadNUniq+1):(met$hadNUniq*2)){
  model7$Pa[,i] <- as.vector(met$hadSLP2rAnnualPressureArray[38:93,met$hadSLP2rUniq[(i-met$hadNUniq)]])
}
# Lag 2 years....
for(i in ((met$hadNUniq*2)+1):(met$hadNUniq*3)){
  model7$Pa[,i] <- as.vector(met$hadSLP2rAnnualPressureArray[37:92,met$hadSLP2rUniq[(i-(met$hadNUniq*2))]])
}

#########
# Brest #
#########

# What we are going to do now is calculate the sd of the residuals as each of the
# components is added in turn 
model7$tg.sd.brest <- vector(mode="numeric", length=(met$hadNUniq*3))

for(i in 1:(met$hadNUniq*3)){
  model7$data.brest <- data.frame(msl=model7$brestTotalP, t=seq(from=1888,to=1943), Pa=model7$Pa[,i])
  model7$tg.lmRob.brest <- lmRob(msl ~ t + Pa, data=model7$data.brest, control=lmRobControl)
  
  model7$tg.sd.brest[i] <- model7$tg.lmRob.brest$scale
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model7$tg.sd.sorted.brest <- sort(model7$tg.sd.brest, index=TRUE)

# Now construct models of the top 9 components
model7$data.brest <- data.frame(msl=model7$brestTotalP,
  t=seq(from=1888,to=1943),
  Pa1=model7$Pa[,model7$tg.sd.sorted.brest$ix[1]], 
  Pa2=model7$Pa[,model7$tg.sd.sorted.brest$ix[2]], Pa3=model7$Pa[,model7$tg.sd.sorted.brest$ix[3]],
  Pa4=model7$Pa[,model7$tg.sd.sorted.brest$ix[4]], Pa5=model7$Pa[,model7$tg.sd.sorted.brest$ix[5]],
  Pa6=model7$Pa[,model7$tg.sd.sorted.brest$ix[6]], Pa7=model7$Pa[,model7$tg.sd.sorted.brest$ix[7]],
  Pa8=model7$Pa[,model7$tg.sd.sorted.brest$ix[8]], Pa9=model7$Pa[,model7$tg.sd.sorted.brest$ix[9]])

model7$tg.lmRob.brest <- lmRob(msl ~ t + Pa1 + Pa2 + Pa3 + Pa4 + Pa5 + Pa6 + Pa7 + Pa8 + Pa9, x=T, y=T,
                               data=model7$data.brest, control=lmRobControl)

# Variance reduction
#> (model7$tg.lmRob.brest$r.sq-model7$tg.lmRob.brestTotalP$r.sq)
#[1] 0.1658905
#> model7$tg.lmRob.brest$r.sq
#[1] 0.5441145
#>sqrt(diag(model7$tg.lmRob.brest$cov))[2]/10.045
#        t 
#0.4507884 
#> coef(model7$tg.lmRob.brest)[2]/10.045
#       t 
#2.336161 

#############
## Model 8 ##
#############
# Follow Thompson (1980) and see if we can improve
# Uses pressures from HadSLP2r
# Calculate over 1850-2008 period, just for Brest ;)

model8<-new.env()

# Just use annual pressures over the periods of the model data (Pa) ending in 2000 (end of model pressure) 
model8$brestPa <-
  met$hadSLP2rAnnualPressureArray[intersect(which(met$hadSLP2rAnnualTime>=1850), which(met$hadSLP2rAnnualTime<=2008)),2]
model8$time <- c(1850:2008)
tg_annual$brest5008 <- tg_annual$brest[44:202,]

# Calculate the total pressure
model8$brestTotalP <- 1025*9.8*tg_annual$brest5008[,2]/1000 + model8$brestPa

model8$tg.lmRob.brestTotalP <- lmRob(model8$brestTotalP ~ model8$time, control=lmRobControl)

model8$Pa <- array(NA,dim=c(159,(met$hadNUniq*3)))
# Zero lag
for(i in 1:met$hadNUniq){
  model8$Pa[,i] <- as.vector(met$hadSLP2rAnnualPressureArray[1:159,met$hadSLP2rUniq[i]])
}
# Lag 1 year
for(i in (met$hadNUniq+1):(met$hadNUniq*2)){
  model8$Pa[2:159,i] <- as.vector(met$hadSLP2rAnnualPressureArray[1:158,met$hadSLP2rUniq[(i-met$hadNUniq)]])
  model8$Pa[1,i] <- mean(as.vector(met$hadSLP2rAnnualPressureArray[1:159,met$hadSLP2rUniq[(i-(met$hadNUniq*2))]]), na.rm=T)
}
# Lag 2 years....
for(i in ((met$hadNUniq*2)+1):(met$hadNUniq*3)){
  model8$Pa[3:159,i] <- as.vector(met$hadSLP2rAnnualPressureArray[1:157,met$hadSLP2rUniq[(i-(met$hadNUniq*2))]])
  model8$Pa[1:2,i] <- mean(as.vector(met$hadSLP2rAnnualPressureArray[1:159,met$hadSLP2rUniq[(i-(met$hadNUniq*2))]]), na.rm=T)
}

#########
# Brest #
#########

# What we are going to do now is calculate the sd of the residuals as each of the
# components is added in turn 
model8$tg.sd.brest <- vector(mode="numeric", length=(met$hadNUniq*3))

for(i in 1:(met$hadNUniq*3)){
  model8$data.brest <- data.frame(msl=model8$brestTotalP, t=seq(from=1850,to=2008), Pa=model8$Pa[,i])
  model8$tg.lmRob.brest <- lmRob(msl ~ t + Pa, data=model8$data.brest, control=lmRobControl)
  
  model8$tg.sd.brest[i] <- model8$tg.lmRob.brest$scale
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model8$tg.sd.sorted.brest <- sort(model8$tg.sd.brest, index=TRUE)

# Now construct models of the top 9 components
model8$data.brest <- data.frame(msl=model8$brestTotalP,
  t=seq(from=1850,to=2008),
  Pa1=model8$Pa[,model8$tg.sd.sorted.brest$ix[1]], 
  Pa2=model8$Pa[,model8$tg.sd.sorted.brest$ix[2]], Pa3=model8$Pa[,model8$tg.sd.sorted.brest$ix[3]],
  Pa4=model8$Pa[,model8$tg.sd.sorted.brest$ix[4]], Pa5=model8$Pa[,model8$tg.sd.sorted.brest$ix[5]],
  Pa6=model8$Pa[,model8$tg.sd.sorted.brest$ix[6]], Pa7=model8$Pa[,model8$tg.sd.sorted.brest$ix[7]],
  Pa8=model8$Pa[,model8$tg.sd.sorted.brest$ix[8]], Pa9=model8$Pa[,model8$tg.sd.sorted.brest$ix[9]])

model8$tg.lmRob.brest <- lmRob(msl ~ t + Pa1 + Pa2 + Pa3 + Pa4 + Pa5 + Pa6 + Pa7 + Pa8 + Pa9, x=T, y=T,
                               data=model8$data.brest, control=lmRobControl)

# Variance reduction
#> (model8$tg.lmRob.brest$r.sq-model8$tg.lmRob.brestTotalP$r.sq)
#[1] 0.1658905
#> model8$tg.lmRob.brest$r.sq
#[1] 0.5441145
#>sqrt(diag(model8$tg.lmRob.brest$cov))[2]/10.045
#        t 
#0.2692527 
#> coef(model8$tg.lmRob.brest)[2]/10.045
#       t 
#1.126984
