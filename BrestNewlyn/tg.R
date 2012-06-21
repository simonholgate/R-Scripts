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
brest <- read.table("~/Dropbox/brestNewlynData/analysis/paper/brestAnnual.rlrdata", col.names=c("Year", "Height", "Flag", "QC"), colClasses=c("numeric","numeric", "character", "character"), sep=";", na.string="-99999")
newlyn <- read.table("~/Dropbox/brestNewlynData/analysis/paper/newlynAnnual.rlrdata", col.names=c("Year", "Height", "Flag", "QC"), colClasses=c("numeric","numeric", "character", "character"), sep=";", na.string="-99999")
keywest <- read.table("~/Dropbox/brestNewlynData/analysis/paper/keyWestAnnual.rlrdata", col.names=c("Year", "Height", "Flag", "QC"), colClasses=c("numeric","numeric", "character", "character"), sep=";", na.string="-99999")
sanfran <- read.table("~/Dropbox/brestNewlynData/analysis/paper/sanFranciscoAnnual.rlrdata", col.names=c("Year", "Height", "Flag", "QC"), colClasses=c("numeric","numeric", "character", "character"), sep=";", na.string="-99999")

# Delete the 1916 value from both time-series so that all the series
# are zero in 1916
newlyn$Height <- newlyn$Height - newlyn$Height[1]
brest$Height <- brest$Height - brest$Height[which(brest$Year==1916)]
keywest$Height <- keywest$Height - keywest$Height[which(keywest$Year==1916)]
sanfran$Height <- sanfran$Height - sanfran$Height[which(sanfran$Year==1916)]

# Use robust methods
library(robust)
lmRobControl <- lmRob.control(mxr=100,mxf=100,trace=F)

## Calculate rate over entire Newlyn period (1916-2008)
newlyn.lmRob <- lmRob(Height ~ Year, data=newlyn, x=T)
#> summary(newlyn.lmRob)
#
#Call: lmRob(formula = Height ~ Year, data = newlyn)
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

newlyn.lm <- lm(Height ~ Year, data=newlyn, x=T)

newlynGap <- newlyn
newlynGap[29:37,2] <- NA
newlynGap.lmRob <- lmRob(Height ~ Year, data=newlynGap, x=T)
#> summary(newlynGap.lmRob)
#
#Call: lmRob(formula = Height ~ Year, data = newlynGap, x = T)
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

newlyn1920 <- newlyn[4:93,]
newlyn1920.lmRob <- lmRob(Height ~ Year, data=newlyn1920, x=T)
#> summary(newlyn1920.lmRob)
#
#Call: lmRob(formula = Height ~ Year, data = newlyn1920, x = T)
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
brest1916 <- brest[110:202,]
brest1916.lmRob <- lmRob(Height ~ Year, data=brest1916, x=T)
#> summary(brest1916.lmRob)
#
#Call: lmRob(formula = Height ~ Year, data = brest1916, x = T)
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

brest1916.lm <- lm(Height ~ Year, data=brest1916, x=T)

brest1920 <- brest[114:202,]
brest1920.lmRob <- lmRob(Height ~ Year, data=brest1920, x=T)
#> summary(brest1920.lmRob)
#
#Call: lmRob(formula = Height ~ Year, data = brest1920, x = T)
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
newlyn1st <- newlyn[1:28,]
newlyn1st.lmRob <- lmRob(Height ~ Year, data=newlyn1st, x=T)
#> summary(newlyn1st.lmRob)
#
#Call: lmRob(formula = Height ~ Year, data = newlyn1st)
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

newlyn2nd <- newlyn[38:93,]
newlyn2nd.lmRob <- lmRob(Height ~ Year, data=newlyn2nd, x=T)
#> summary(newlyn2nd.lmRob)
#
#Call: lmRob(formula = Height ~ Year, data = newlyn2nd)
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

# Introduce a third section of 1992-2008
newlyn3rd <- newlyn[77:93,]
newlyn3rd.lmRob <- lmRob(Height ~ Year, data=newlyn3rd, x=T)
#> summary(newlyn3rd.lmRob)

## Call: lmRob(formula = Height ~ Year, data = newlyn3rd, x = T)

## Residuals:
##        Min         1Q     Median         3Q        Max 
## -34.525632 -12.250767   1.958521  14.016175  64.700598 

## Coefficients:
##             Value         Std. Error    t value       Pr(>|t|)     
## (Intercept) -9.358338e+03  3.921977e+03 -2.386128e+00  3.169953e-02
## Year         4.742077e+00  1.961542e+00  2.417525e+00  2.984918e-02

## Residual standard error: 24.8797 on 14 degrees of freedom
## 1 observation deleted due to missingness 
## Multiple R-Squared: 0.41861 

## Test for Bias:
##             statistic   p-value
## M-estimate  0.6616453 0.7183326
## LS-estimate 0.5308272 0.7668887

# Do the same for Brest
brest1st <- brest[110:137,]
brest1st.lmRob <- lmRob(Height ~ Year, data=brest1st, x=T)
#>  summary(brest1st.lmRob) 
#
#Call: lmRob(formula = Height ~ Year, data = brest1st)
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

brest2nd <- brest[147:202,]
brest2nd.lmRob <- lmRob(Height ~ Year, data=brest2nd, x=T)
#> summary(brest2nd.lmRob)
#
#Call: lmRob(formula = Height ~ Year, data = brest2nd)
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

# Introduce a third section of 1992-2008
brest3rd <- brest[186:202,]
brest3rd.lmRob <- lmRob(Height ~ Year, data=brest3rd, x=T)
#> summary(newlyn3rd.lmRob)

## Call: lmRob(formula = Height ~ Year, data = newlyn3rd, x = T)

## Residuals:
##        Min         1Q     Median         3Q        Max 
## -34.525632 -12.250767   1.958521  14.016175  64.700598 

## Coefficients:
##             Value         Std. Error    t value       Pr(>|t|)     
## (Intercept) -9.358338e+03  3.921977e+03 -2.386128e+00  3.169953e-02
## Year         4.742077e+00  1.961542e+00  2.417525e+00  2.984918e-02

## Residual standard error: 24.8797 on 14 degrees of freedom
## 1 observation deleted due to missingness 
## Multiple R-Squared: 0.41861 

## Test for Bias:
##             statistic   p-value
## M-estimate  0.6616453 0.7183326
## LS-estimate 0.5308272 0.7668887

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


# Delfzijl
tg_monthly <- new.env()
tg_monthly$delfzijl <- read.table("~/Dropbox/brestNewlynData/analysis/paper/delfzijlMonthly.metdata", col.names=c("Year", "Height", "Flag", "QC"), colClasses=c("numeric","numeric", "character", "character"), sep=";", na.string="-99999")

## Create annual means
tmp<-new.env()
tmp$delfzijl <-  tg_monthly$delfzijl$Height
dim(tmp$delfzijl) <- c(12,143)
tmp$delfzijlAnnual <- colMeans(tmp$delfzijl)
delfzijl <- data.frame(c(1865:2007), tmp$delfzijlAnnual)
colnames(delfzijl) <- c("Year", "Height")
rm(tmp)
rm(tg_monthly)

delfzijl1916 <- delfzijl[52:143,]
delfzijl1916.lmRob <- lmRob(Height ~ Year, data=delfzijl1916, x=T)
#> summary(delfzijl1916.lmRob)
#
#Call: lmRob(formula = Height ~ Year, data = delfzijl1916, 
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

save(file="tg.RData", list=ls())
