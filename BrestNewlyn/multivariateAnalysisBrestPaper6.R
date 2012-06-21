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

# Introduce a third section of 1992-2008
tg_annual$newlyn3rd <- tg_annual$newlyn[77:93,]
tg_annual$newlyn3rd.lmRob <- lmRob(Height ~ Year, data=tg_annual$newlyn3rd, x=T)
#> summary(tg_annual$newlyn3rd.lmRob)

## Call: lmRob(formula = Height ~ Year, data = tg_annual$newlyn3rd, x = T)

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

# Introduce a third section of 1992-2008
tg_annual$brest3rd <- tg_annual$brest[186:202,]
tg_annual$brest3rd.lmRob <- lmRob(Height ~ Year, data=tg_annual$brest3rd, x=T)
#> summary(tg_annual$newlyn3rd.lmRob)

## Call: lmRob(formula = Height ~ Year, data = tg_annual$newlyn3rd, x = T)

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



tmp <- new.env()
# HadSLP2r data is 1850-2009 = 160 years = 1920 months. 37 cols - Newlyn, Brest and 35 grid points
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

# CRU data is 1873-2000 = 128 years = 1536 months. 
met$cru.brest<-read.table("~/Dropbox/brestNewlynData/pressure/cru.brest.1873.2000.txt", col.names=c("Year","SLP"))


lmRobControl <- lmRob.control(mxr=100,mxf=100,trace=F)

#######################
## Model 3 1953-2008 ##
#######################
# Follow Thompson (1980) and see if we can improve
# Uses HadSLP2r
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

model3$tg.lmRob.brestTotalP <- lmRob(model3$brestTotalP ~ model3$time, control=lmRobControl, x=T)
model3$tg.lmRob.newlynTotalP <- lmRob(model3$newlynTotalP ~ model3$time, control=lmRobControl, x=T)

model3$Pa <- array(NA,dim=c(56,met$hadNUniq))
# Zero lag
for(i in 1:met$hadNUniq){
  model3$Pa[,i] <- as.vector(met$hadSLP2rAnnualPressureArray[104:159,met$hadSLP2rUniq[i]])
}

#########
# Brest #
#########

# What we are going to do now is calculate the sd of the residuals as each of the
# components is added in turn 
model3$tg.sd.brest <- vector(mode="numeric", length=met$hadNUniq)

for(i in 1:met$hadNUniq){
  model3$data.brest <- data.frame(msl=model3$brestTotalP, t=seq(from=1953, to=2008), Pa=model3$Pa[,i])
  model3$tg.lmRob.brest <- lmRob(msl ~ t + Pa, data=model3$data.brest, control=lmRobControl, na.action=na.exclude)
  model3$tg.sd.brest[i] <- model3$tg.lmRob.brest$scale
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model3$tg.sd.sorted.brest <- sort(model3$tg.sd.brest, index=TRUE)

# Now construct models of the top 9 components
model3$data.brest <- data.frame(msl=model3$brestTotalP,
  t=seq(from=1953,to=2008), model3$Pa[,model3$tg.sd.sorted.brest$ix[1:9]])

model3$tg.lmRob.brest <- lmRob(msl ~ ., x=T, y=T, data=model3$data.brest, control=lmRobControl, na.action=na.exclude)

# Test model prediction skill
# Need to build new set of top 9 components for partial set
# What we are going to do now is calculate the sd of the residuals as each of the
# components is added in turn 
model3$tg.sd.brest.partial <- vector(mode="numeric", length=met$hadNUniq)

for(i in 1:met$hadNUniq){
  model3$data.brest.partial <- data.frame(msl=model3$brestTotalP[1:37], t=seq(from=1953,to=1989), Pa=model3$Pa[1:37,i])
  model3$tg.lmRob.brest.partial <- lmRob(msl ~ ., data=model3$data.brest.partial, control=lmRobControl, na.action=na.exclude)
  
  model3$tg.sd.brest.partial[i] <- model3$tg.lmRob.brest.partial$scale
}

# Now we can sort these sds. 
model3$tg.sd.sorted.brest.partial<- sort(model3$tg.sd.brest.partial, index=TRUE)

model3$data.brest.partial <- data.frame(msl=model3$brestTotalP[1:37],
  t=seq(from=1953,to=1989), model3$Pa[1:37,model3$tg.sd.sorted.brest.partial$ix[1:9]])

model3$tg.lmRob.brest.partial <- lmRob(msl ~ ., x=T, y=T,
                               data=model3$data.brest.partial, control=lmRobControl, na.action=na.exclude)

model3$data.brest.new<- data.frame(t=seq(from=1990,to=2008),
  model3$Pa[38:56,model3$tg.sd.sorted.brest.partial$ix[1:9]])

model3$tg.lmRob.brest.pred<- predict.lmRob(model3$tg.lmRob.brest.partial, se.fit=T, newdata=model3$data.brest.new, interval='confidence')

# Calculate correlation array. What is the relationship between each of the pressure series and the sea level time-series?
# Columns are Brest, Newlyn then Pa
model3$corr.array <- array(NA,dim=c(104,104))
model3$corr.data <- cbind(model3$brestTotalP,model3$newlynTotalP,model3$Pa)
for(i in 1:104){
  for(j in 1:104){
    model3$corr.array[i,j] <- cor.test(model3$corr.data[,i], model3$corr.data[,j], method="p", alternative="t", na.action=na.fail)$estimate
  }
}

### How does a model based on the most correlated 9 components, compare?
model3$tg.cor.sorted.brest <- sort(abs(model3$corr.array[1,3:37]), index=T, decreasing=T)$ix[1:9]
model3$data.brest.partial <- data.frame(msl=model3$brestTotalP[1:37], t=seq(from=1953,to=1989), Pa=model3$Pa[1:37,model3$tg.cor.sorted.brest])
model3$data.brest.new <- data.frame(t=seq(from=1990,to=2008), Pa=model3$Pa[38:56,3:37])
model3$tg.cor.lmRob.brest.partial <- lmRob(msl ~ ., data=model3$data.brest.partial)
model3$tg.cor.lmRob.brest.pred <- predict.lmRob(model3$tg.cor.lmRob.brest.partial, se.fit=T, newdata=model3$data.brest.new, interval='confidence')

# The predictions from the most correlated pressures are worse than from those with the R^2 technique

##########
# Newlyn #
##########

# What we are going to do now is calculate the sd of the residuals as each of the
# components is added in turn 
model3$tg.sd.newlyn <- vector(mode="numeric", length=met$hadNUniq)

for(i in 1:met$hadNUniq){
  model3$data.newlyn <- data.frame(msl=model3$newlynTotalP, t=seq(from=1953, to=2008), Pa=model3$Pa[,i])
  model3$tg.lmRob.newlyn <- lmRob(msl ~ t + Pa, data=model3$data.newlyn, control=lmRobControl, na.action=na.exclude)
  
  model3$tg.sd.newlyn[i] <- model3$tg.lmRob.newlyn$scale
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model3$tg.sd.sorted.newlyn <- sort(model3$tg.sd.newlyn, index=TRUE)

# Now construct models of the top 9 components
model3$data.newlyn <- data.frame(msl=model3$newlynTotalP,
  t=seq(from=1953,to=2008),
  model3$Pa[,model3$tg.sd.sorted.newlyn$ix[1:9]])

model3$tg.lmRob.newlyn <- lmRob(msl ~ ., x=T, y=T, data=model3$data.newlyn, control=lmRobControl, na.action=na.exclude)

# Test model prediction skill
# Need to build new set of top 9 components for partial set
# What we are going to do now is calculate the sd of the residuals as each of the
# components is added in turn 
model3$tg.sd.newlyn.partial <- vector(mode="numeric", length=met$hadNUniq)

for(i in 1:met$hadNUniq){
  model3$data.newlyn.partial <- data.frame(msl=model3$newlynTotalP[1:37], t=seq(from=1953,to=1989), Pa=model3$Pa[1:37,i])
  model3$tg.lmRob.newlyn.partial <- lmRob(msl ~ ., data=model3$data.newlyn.partial, control=lmRobControl, na.action=na.exclude)
  
  model3$tg.sd.newlyn.partial[i] <- model3$tg.lmRob.newlyn.partial$scale
}

# Now we can sort these sds. 
model3$tg.sd.sorted.newlyn.partial<- sort(model3$tg.sd.newlyn.partial, index=TRUE)

model3$data.newlyn.partial <- data.frame(msl=model3$newlynTotalP[1:37],
  t=seq(from=1953,to=1989), model3$Pa[1:37,model3$tg.sd.sorted.newlyn.partial$ix[1:9]])

model3$tg.lmRob.newlyn.partial <- lmRob(msl ~ ., x=T, y=T,
                               data=model3$data.newlyn.partial, control=lmRobControl, na.action=na.exclude)

model3$data.newlyn.new<- data.frame(t=seq(from=1990,to=2008),
  model3$Pa[38:56,model3$tg.sd.sorted.newlyn.partial$ix[1:9]])

model3$tg.lmRob.newlyn.pred<- predict.lmRob(model3$tg.lmRob.newlyn.partial, se.fit=T, newdata=model3$data.newlyn.new, interval='confidence')

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

model4$tg.lmRob.brestTotalP <- lmRob(model4$brestTotalP ~ model4$time, control=lmRobControl, na.action=na.exclude)
model4$tg.lmRob.newlynTotalP <- lmRob(model4$newlynTotalP ~ model4$time, x=T, control=lmRobControl, na.action=na.exclude)

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

# Calculate correlation array. What is the relationship between each of the pressure series and the sea level time-series?
# Columns are Brest, Newlyn then Pa
model4$corr.array <- array(NA,dim=c(104,104))
model4$corr.data <- cbind(model4$brestTotalP,model4$newlynTotalP,model4$Pa)
for(i in 1:104){
  for(j in 1:104){
    model4$corr.array[i,j] <- cor.test(model4$corr.data[,i], model4$corr.data[,j], method="p", alternative="t", na.action=na.fail)$estimate
  }
}
                

#########
# Brest #
#########

# What we are going to do now is calculate the sd of the residuals as each of the
# components is added in turn 
model4$tg.sd.brest <- vector(mode="numeric", length=(met$hadNUniq*3))

for(i in 1:(met$hadNUniq*3)){
  model4$data.brest <- data.frame(msl=model4$brestTotalP, t=seq(from=1916,to=1943), Pa=model4$Pa[,i])
  model4$tg.lmRob.brest <- lmRob(msl ~ t + Pa, data=model4$data.brest, control=lmRobControl, na.action=na.exclude)
  
  model4$tg.sd.brest[i] <- model4$tg.lmRob.brest$scale
  model4$tg.rsq.brest[i] <- model4$tg.lmRob.brest$r.sq
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model4$tg.sd.sorted.brest <- sort(model4$tg.sd.brest, index=TRUE)
model4$tg.rsq.sorted.brest <- sort(model4$tg.rsq.brest, index=TRUE, decreasing=TRUE)

# Now construct models of the top 9 components
model4$data.brest <- data.frame(msl=model4$brestTotalP,
  t=seq(from=1916,to=1943), model4$Pa[,model4$tg.sd.sorted.brest$ix[1:9]]])

model4$tg.lmRob.brest <- lmRob(msl ~ ., x=T, y=T, data=model4$data.brest, control=lmRobControl, na.action=na.exclude)

# Test model prediction skill
model4$data.brest.partial <- data.frame(msl=model4$brestTotalP[1:18],
  t=seq(from=1916,to=1933), model4$Pa[1:18,model4$tg.sd.sorted.brest$ix[1:9]])

model4$tg.lmRob.brest.partial <- lmRob(msl ~ ., x=T, y=T,
                               data=model4$data.brest.partial, control=lmRobControl, na.action=na.exclude)

model4$data.brest.new <- data.frame(t=seq(from=1934,to=1943),
  model4$Pa[19:28,model4$tg.sd.sorted.brest$ix[1:9]])

model4$tg.lmRob.brest.partial <- lmRob(msl ~ ., data=model4$data.brest.partial)
model4$tg.lmRob.brest.pred <- predict.lmRob(model4$tg.lmRob.brest.partial, se.fit=T, newdata=model4$data.brest.new, interval='confidence')


model4$data.brest <- data.frame(msl=model4$brestTotalP,
  t=seq(from=1916,to=1943),
  model4$Pa[,model4$tg.rsq.sorted.brest$ix[1:9]])

model4$tg.lmRob.rsq.brest <- lmRob(msl ~ ., x=T, y=T, data=model4$data.brest, control=lmRobControl, na.action=na.exclude)

# Test model prediction skill
model4$data.brest.partial <- data.frame(msl=model4$brestTotalP[1:18],
  t=seq(from=1916,to=1933),
  model4$Pa[1:18,model4$tg.rsq.sorted.brest$ix[1:9]])

model4$tg.lmRob.brest.partial <- lmRob(msl ~ ., x=T, y=T, data=model4$data.brest.partial, control=lmRobControl, na.action=na.exclude)

model4$data.brest.new <- data.frame(t=seq(from=1934,to=1943),
  model4$Pa[19:28,model4$tg.rsq.sorted.brest$ix[1:9]])

model4$tg.lmRob.brest.partial <- lmRob(msl ~ ., data=model4$data.brest.partial)
model4$tg.lmRob.brest.pred <- predict.lmRob(model4$tg.lmRob.brest.partial, se.fit=T, newdata=model4$data.brest.new, interval='confidence')
# Again sd is better for predictor than rsq

# Variance reduction
#> (model4$tg.lmRob.brest$r.sq-model4$tg.lmRob.brestTotalP$r.sq)
#[1] 0.3982849
#> model4$tg.lmRob.brest$r.sq
#[1] 0.6360687

#> model4$tg.lmRob.rsq.brest$r.sq
#[1] 0.5790682

##########
# Newlyn #
##########

# What we are going to do now is calculate the sd of the residuals as each of the
# components is added in turn 
model4$tg.sd.newlyn <- vector(mode="numeric", length=(met$hadNUniq*3))

for(i in 1:(met$hadNUniq*3)){
  model4$data.newlyn <- data.frame(msl=model4$newlynTotalP, t=seq(from=1916, to=1943), Pa=model4$Pa[,i])
  model4$tg.lmRob.newlyn <- lmRob(msl ~ t + Pa, data=model4$data.newlyn, control=lmRobControl, na.action=na.exclude)
  
  model4$tg.sd.newlyn[i] <- model4$tg.lmRob.newlyn$scale
  model4$tg.rsq.newlyn[i] <- model4$tg.lmRob.newlyn$r.sq
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model4$tg.sd.sorted.newlyn <- sort(model4$tg.sd.newlyn, index=TRUE)
model4$tg.rsq.sorted.newlyn <- sort(model4$tg.rsq.newlyn, index=TRUE, decreasing=TRUE)

# Now construct models of the top 9 components
model4$data.newlyn <- data.frame(msl=model4$newlynTotalP,
  t=seq(from=1916,to=1943), model4$Pa[,model4$tg.sd.sorted.newlyn$ix[1:9]])

model4$tg.lmRob.newlyn <- lmRob(msl ~ ., x=T, y=T, data=model4$data.newlyn, control=lmRobControl, na.action=na.exclude)

# Test model prediction skill
model4$data.newlyn.partial <- data.frame(msl=model4$newlynTotalP[1:18],
  t=seq(from=1916,to=1933),
  model4$Pa[1:18,model4$tg.sd.sorted.newlyn$ix[1:9]])

model4$tg.lmRob.newlyn.partial <- lmRob(msl ~ ., x=T, y=T,
                               data=model4$data.newlyn.partial, control=lmRobControl, na.action=na.exclude)

model4$data.newlyn.new <- data.frame(t=seq(from=1934,to=1943),
  model4$Pa[19:28,model4$tg.sd.sorted.newlyn$ix[1:9]])

model4$tg.lmRob.newlyn.partial <- lmRob(msl ~ ., data=model4$data.newlyn.partial)
model4$tg.lmRob.newlyn.pred <- predict.lmRob(model4$tg.lmRob.newlyn.partial, se.fit=T, newdata=model4$data.newlyn.new, interval='confidence')

# Now construct models of the top 9 components
model4$data.newlyn <- data.frame(msl=model4$newlynTotalP,
  t=seq(from=1916,to=1943),
  model4$Pa[,model4$tg.rsq.sorted.newlyn$ix[1:9]])

model4$tg.lmRob.rsq.newlyn <- lmRob(msl ~ ., x=T, y=T, data=model4$data.newlyn, control=lmRobControl, na.action=na.exclude)
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

#> model4$tg.lmRob.rsq.newlyn$r.sq
#[1] 0.742686
#> c(coef(model4$tg.lmRob.rsq.brest)[2]/10.045, coef(model4$tg.lmRob.rsq.newlyn)[2]/10.045)
#       t        t 
#2.295370 2.122536
#> c(sqrt(diag(model4$tg.lmRob.rsq.brest$cov))[2]/10.045,  sqrt(diag(model4$tg.lmRob.rsq.newlyn$cov))[2]/10.045)
#        t         t 
#0.3551371 0.3473025

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
# Subtract 22mm to Brest pre-war to see what effect this has
tg_annual$brest1608_22_2 <- tg_annual$brest1608
tg_annual$brest1608_22_2[1:37,2] <- tg_annual$brest1608_22[1:37,2]-22

tg_annual$newlyn1608 <- tg_annual$newlyn[1:93,]

# Calculate the total pressure
model5$brestTotalP <- 1025*9.8*tg_annual$brest1608[,2]/1000 + model5$brestPa
model5$brestTotalP_22<- 1025*9.8*tg_annual$brest1608_22[,2]/1000 + model5$brestPa
model5$brestTotalP_22_2<- 1025*9.8*tg_annual$brest1608_22_2[,2]/1000 + model5$brestPa
model5$newlynTotalP <- 1025*9.8*tg_annual$newlyn1608[,2]/1000 + model5$newlynPa
model5$delfzijlTotalP <- 1025*9.8*tg_annual$delfzijl1916[,2]/1000 + model5$newlynPa[1:92]

model5$tg.lmRob.brestTotalP <- lmRob(model5$brestTotalP ~ model5$time, x=T, y=T, control=lmRobControl, na.action=na.exclude)
model5$tg.lmRob.brestTotalP_22<- lmRob(model5$brestTotalP_22 ~ model5$time, control=lmRobControl, na.action=na.exclude)
model5$tg.lmRob.newlynTotalP <- lmRob(model5$newlynTotalP ~ model5$time, x=T, y=T, control=lmRobControl, na.action=na.exclude)

model5$Pa <- array(NA,dim=c(93,met$hadNUniq*))
# Zero lag
for(i in 1:met$hadNUniq){
  model5$Pa[,i] <- as.vector(met$hadSLP2rAnnualPressureArray[67:159,met$hadSLP2rUniq[i]])
}


#########
# Brest #
#########

# What we are going to do now is calculate the sd of the residuals as each of the
# components is added in turn 
model5$tg.sd.brest <- vector(mode="numeric", length=(met$hadNUniq))

for(i in 1:(met$hadNUniq)){
  model5$data.brest <- data.frame(msl=model5$brestTotalP, t=seq(from=1916,to=2008), Pa=model5$Pa[,i])
  model5$tg.lmRob.brest <- lmRob(msl ~ t + Pa, data=model5$data.brest, control=lmRobControl, na.action=na.exclude)
  
  model5$tg.sd.brest[i] <- model5$tg.lmRob.brest$scale
  model5$tg.rsq.brest[i] <- model5$tg.lmRob.brest$r.sq
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model5$tg.sd.sorted.brest <- sort(model5$tg.sd.brest, index=TRUE)
model5$tg.rsq.sorted.brest <- sort(model5$tg.rsq.brest, index=TRUE, decreasing=TRUE)

# Now construct models of the top 9 components
model5$data.brest <- data.frame(msl=model5$brestTotalP,
  t=seq(from=1916,to=2008),
  model5$Pa[,model5$tg.sd.sorted.brest$ix[1:9]])

model5$tg.lmRob.brest <- lmRob(msl ~ ., x=T, y=T, data=model5$data.brest, control=lmRobControl, na.action=na.exclude)
model5$tg.lmRob.brest.fitted.lm <- lmRob(model5$tg.lmRob.brest$fitted  ~ model5$tg.lmRob.brest$x[,2], control=lmRobControl, na.action=na.exclude)

model5$data.brest <- data.frame(msl=model5$brestTotalP,
  t=seq(from=1916,to=2008),
  model5$Pa[,model5$tg.rsq.sorted.brest$ix[1:9]])

model5$tg.lmRob.rsq.brest <- lmRob(msl ~ ., x=T, y=T, data=model5$data.brest, control=lmRobControl, na.action=na.exclude)

# Test model prediction skill
# Need to build new set of top 9 components for partial set
# What we are going to do now is calculate the sd of the residuals as each of the
# components is added in turn 
model5$tg.sd.brest.partial <- vector(mode="numeric", length=(met$hadNUniq*3))

for(i in 1:(met$hadNUniq*3)){
  model5$data.brest.partial <- data.frame(msl=model5$brestTotalP[38:93], t=seq(from=1953,to=2008), Pa=model5$Pa[38:93,i])
  model5$tg.lmRob.brest.partial <- lmRob(msl ~ ., data=model5$data.brest.partial, control=lmRobControl, na.action=na.exclude)
  
  model5$tg.sd.brest.partial[i] <- model5$tg.lmRob.brest.partial$scale
}

# Now we can sort these sds. 
model5$tg.sd.sorted.brest.partial<- sort(model5$tg.sd.brest.partial, index=TRUE)

model5$data.brest.partial <- data.frame(msl=model5$brestTotalP[38:93],
  t=seq(from=1953,to=2008),
  model5$Pa[38:93,model5$tg.sd.sorted.brest.partial$ix[1:9]])

model5$tg.lmRob.brest.partial <- lmRob(msl ~ ., x=T, y=T,
                               data=model5$data.brest.partial, control=lmRobControl, na.action=na.exclude)

model5$data.brest.new<- data.frame(t=seq(from=1916,to=1952),
  model5$Pa[1:37,model5$tg.sd.sorted.brest.partial$ix[1:9]])

model5$tg.lmRob.brest.pred<- predict.lmRob(model5$tg.lmRob.brest.partial, se.fit=T, newdata=model5$data.brest.new, interval='confidence')

# Variance reduction
#> (model5$tg.lmRob.brest$r.sq-model5$tg.lmRob.brestTotalP$r.sq)
#[1] 0.1072223
#> model5$tg.lmRob.brest$r.sq
#[1] 0.6555045
#> model5$tg.lmRob.rsq.brest$r.sq
#[1] 0.6398969

#################
# Brest plus 22 #
#################

# What we are going to do now is calculate the sd of the residuals as each of the
# components is added in turn 
model5$tg.sd.brest_22 <- vector(mode="numeric", length=(met$hadNUniq))

for(i in 1:(met$hadNUniq)){
  model5$data.brest_22<- data.frame(msl=model5$brestTotalP_22, t=seq(from=1916,to=2008), Pa=model5$Pa[,i])
  model5$tg.lmRob.brest_22 <- lmRob(msl ~ t + Pa, data=model5$data.brest_22, control=lmRobControl, na.action=na.exclude)
  
  model5$tg.sd.brest_22[i] <- model5$tg.lmRob.brest_22$scale
  model5$tg.rsq.brest_22[i] <- model5$tg.lmRob.brest_22$r.sq
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model5$tg.sd.sorted.brest_22<- sort(model5$tg.sd.brest_22, index=TRUE)
model5$tg.rsq.sorted.brest_22<- sort(model5$tg.rsq.brest_22, index=TRUE, decreasing=TRUE)

# Now construct models of the top 9 components
model5$data.brest_22<- data.frame(msl=model5$brestTotalP_22,
  t=seq(from=1916,to=2008),
  model5$Pa[,model5$tg.sd.sorted.brest_22$ix[1:9]])

model5$tg.lmRob.brest_22 <- lmRob(msl ~ ., x=T, y=T, data=model5$data.brest_22, control=lmRobControl, na.action=na.exclude)

model5$data.brest <- data.frame(msl=model5$brestTotalP_22,
  t=seq(from=1916,to=2008),
  model5$Pa[,model5$tg.rsq.sorted.brest_22$ix[1:9]])

model5$tg.lmRob.rsq.brest_22 <- lmRob(msl ~ ., x=T, y=T, data=model5$data.brest_22, control=lmRobControl, na.action=na.exclude)

# Test model prediction skill
# Need to build new set of top 9 components for partial set
# What we are going to do now is calculate the sd of the residuals as each of the
# components is added in turn 
model5$tg.sd.brest.partial_22 <- vector(mode="numeric", length=(met$hadNUniq*3))

for(i in 1:(met$hadNUniq*3)){
  model5$data.brest.partial_22 <- data.frame(msl=model5$brestTotalP_22[38:93], t=seq(from=1953,to=2008), Pa=model5$Pa[38:93,i])
  model5$tg.lmRob.brest.partial_22 <- lmRob(msl ~ ., data=model5$data.brest.partial_22, control=lmRobControl, na.action=na.exclude)
  
  model5$tg.sd.brest.partial_22[i] <- model5$tg.lmRob.brest.partial_22$scale
}

# Now we can sort these sds. 
model5$tg.sd.sorted.brest.partial_22<- sort(model5$tg.sd.brest.partial_22, index=TRUE)

model5$data.brest.partial_22 <- data.frame(msl=model5$brestTotalP_22[38:93],
  t=seq(from=1953,to=2008),
  model5$Pa[38:93,model5$tg.sd.sorted.brest.partial_22$ix[1:9]])

model5$tg.lmRob.brest.partial_22 <- lmRob(msl ~ ., x=T, y=T,
                               data=model5$data.brest.partial_22, control=lmRobControl, na.action=na.exclude)

model5$data.brest.new_22<- data.frame(t=seq(from=1916,to=1952),
  model5$Pa[1:37,model5$tg.sd.sorted.brest.partial_22$ix[1:9]])

model5$tg.lmRob.brest.pred_22<- predict.lmRob(model5$tg.lmRob.brest.partial_22, se.fit=T, newdata=model5$data.brest.new_22, interval='confidence')

# Residual sum of squares
sum((model5$brestTotalP[1:28]-model5$tg.lmRob.brest.pred_22$fit[1:28])^2)
#[1] 1746609 No 22: 3107529 Newlyn: 974070
sum((model5$brestTotalP[38:93]-model5$tg.lmRob.brest.partial_22$fitted)^2)
#[1] 1657279 No 22: 4085193 Newlyn: 1210126


# Variance reduction
#> (model5$tg.lmRob.brest_22$r.sq-model5$tg.lmRob.brestTotalP_22$r.sq)
#[1] 0.113057
#> model5$tg.lmRob.brest_22$r.sq
#[1] 0.7397892
#> model5$tg.lmRob.rsq.brest_22$r.sq
#[1] 0.7429975

#################
# Brest minus 22 #
#################

# What we are going to do now is calculate the sd of the residuals as each of the
# components is added in turn 
model5$tg.sd.brest_22_2 <- vector(mode="numeric", length=(met$hadNUniq))

for(i in 1:(met$hadNUniq)){
  model5$data.brest_22_2<- data.frame(msl=model5$brestTotalP_22_2, t=seq(from=1916,to=2008), Pa=model5$Pa[,i])
  model5$tg.lmRob.brest_22_2 <- lmRob(msl ~ t + Pa, data=model5$data.brest_22_2, control=lmRobControl, na.action=na.exclude)
  
  model5$tg.sd.brest_22_2[i] <- model5$tg.lmRob.brest_22_2$scale
  model5$tg.rsq.brest_22_2[i] <- model5$tg.lmRob.brest_22_2$r.sq
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model5$tg.sd.sorted.brest_22_2<- sort(model5$tg.sd.brest_22_2, index=TRUE)
model5$tg.rsq.sorted.brest_22_2<- sort(model5$tg.rsq.brest_22_2, index=TRUE, decreasing=TRUE)

# Now construct models of the top 9 components
model5$data.brest_22_2<- data.frame(msl=model5$brestTotalP_22_2,
  t=seq(from=1916,to=2008),
  model5$Pa[,model5$tg.sd.sorted.brest_22_2$ix[1:9]])

model5$tg.lmRob.brest_22_2 <- lmRob(msl ~ ., x=T, y=T, data=model5$data.brest_22_2, control=lmRobControl, na.action=na.exclude)

model5$data.brest <- data.frame(msl=model5$brestTotalP_22_2,
  t=seq(from=1916,to=2008),
  model5$Pa[,model5$tg.rsq.sorted.brest_22_2$ix[1:9]])

model5$tg.lmRob.rsq.brest_22_2 <- lmRob(msl ~ ., x=T, y=T, data=model5$data.brest_22_2, control=lmRobControl, na.action=na.exclude)

# Test model prediction skill
# Need to build new set of top 9 components for partial set
# What we are going to do now is calculate the sd of the residuals as each of the
# components is added in turn 
model5$tg.sd.brest.partial_22_2 <- vector(mode="numeric", length=(met$hadNUniq*3))

for(i in 1:(met$hadNUniq*3)){
  model5$data.brest.partial_22_2 <- data.frame(msl=model5$brestTotalP_22_2[38:93], t=seq(from=1953,to=2008), Pa=model5$Pa[38:93,i])
  model5$tg.lmRob.brest.partial_22_2 <- lmRob(msl ~ ., data=model5$data.brest.partial_22_2, control=lmRobControl, na.action=na.exclude)
  
  model5$tg.sd.brest.partial_22_2[i] <- model5$tg.lmRob.brest.partial_22_2$scale
}

# Now we can sort these sds. 
model5$tg.sd.sorted.brest.partial_22_2<- sort(model5$tg.sd.brest.partial_22_2, index=TRUE)

model5$data.brest.partial_22_2 <- data.frame(msl=model5$brestTotalP_22_2[38:93],
  t=seq(from=1953,to=2008),
  model5$Pa[38:93,model5$tg.sd.sorted.brest.partial_22_2$ix[1:9]])

model5$tg.lmRob.brest.partial_22_2 <- lmRob(msl ~ ., x=T, y=T,
                               data=model5$data.brest.partial_22_2, control=lmRobControl, na.action=na.exclude)

model5$data.brest.new_22_2<- data.frame(t=seq(from=1916,to=1952),
  model5$Pa[1:37,model5$tg.sd.sorted.brest.partial_22_2$ix[1:9]])

model5$tg.lmRob.brest.pred_22_2<- predict.lmRob(model5$tg.lmRob.brest.partial_22_2, se.fit=T, newdata=model5$data.brest.new_22_2, interval='confidence')

# Residual sum of squares
sum((model5$data.brest$msl[1:28]-model5$tg.lmRob.brest.pred_22_2$fit[1:28])^2)
#[1] 1746609 No 22: 3107529 Newlyn: 974070
sum((model5$data.brest$msl[38:93]-model5$tg.lmRob.brest.partial_22_2$fitted)^2)
#[1] 1657279 No 22: 4085193 Newlyn: 1210126


# Variance reduction
#> (model5$tg.lmRob.brest_22_2$r.sq-model5$tg.lmRob.brestTotalP_22_2$r.sq)
#[1] 0.113057
#> model5$tg.lmRob.brest_22_2$r.sq
#[1] 0.7397892
#> model5$tg.lmRob.rsq.brest_22_2$r.sq
#[1] 0.7429975

##########
# Newlyn #
##########

# What we are going to do now is calculate the sd of the residuals as each of the
# components is added in turn 
model5$tg.sd.newlyn <- vector(mode="numeric", length=(met$hadNUniq*3))

for(i in 1:(met$hadNUniq*3)){
  model5$data.newlyn <- data.frame(msl=model5$newlynTotalP, t=seq(from=1916, to=2008), Pa=model5$Pa[,i])
  model5$tg.lmRob.newlyn <- lmRob(msl ~ t + Pa, data=model5$data.newlyn, control=lmRobControl, na.action=na.exclude)
  
  model5$tg.sd.newlyn[i] <- model5$tg.lmRob.newlyn$scale
  model5$tg.rsq.newlyn[i] <- model5$tg.lmRob.newlyn$r.sq
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model5$tg.sd.sorted.newlyn <- sort(model5$tg.sd.newlyn, index=TRUE)
model5$tg.rsq.sorted.newlyn <- sort(model5$tg.rsq.newlyn, index=TRUE, decreasing=TRUE)

# Now construct models of the top 9 components
model5$data.newlyn <- data.frame(msl=model5$newlynTotalP, t=seq(from=1916,to=2008), model5$Pa[,model5$tg.sd.sorted.newlyn$ix[1:9]])

model5$tg.lmRob.newlyn <- lmRob(msl ~ ., x=T, y=T, data=model5$data.newlyn, control=lmRobControl, na.action=na.exclude)
model5$tg.lmRob.newlyn.fitted.lm <- lmRob(model5$tg.lmRob.newlyn$fitted  ~ model5$tg.lmRob.newlyn$x[,2], control=lmRobControl, na.action=na.exclude)

# Test model prediction skill
# Need to build new set of top 9 components for partial set
# What we are going to do now is calculate the sd of the residuals as each of the
# components is added in turn 
model5$tg.sd.newlyn.partial <- vector(mode="numeric", length=(met$hadNUniq*3))

for(i in 1:(met$hadNUniq*3)){
  model5$data.newlyn.partial <- data.frame(msl=model5$newlynTotalP[38:91], t=seq(from=1953,to=2006), Pa=model5$Pa[38:91,i])
  model5$tg.lmRob.newlyn.partial <- lmRob(msl ~ ., data=model5$data.newlyn.partial, control=lmRobControl, na.action=na.exclude)
  
  model5$tg.sd.newlyn.partial[i] <- model5$tg.lmRob.newlyn.partial$scale
}

# Now we can sort these sds. 
model5$tg.sd.sorted.newlyn.partial <- sort(model5$tg.sd.newlyn.partial, index=TRUE)

model5$data.newlyn.partial <- data.frame(msl=model5$newlynTotalP[38:91], t=seq(from=1953,to=2006),
  model5$Pa[38:91,model5$tg.sd.sorted.newlyn.partial$ix[1:9]])

model5$tg.lmRob.newlyn.partial <- lmRob(msl ~ ., x=T, y=T, data=model5$data.newlyn.partial, control=lmRobControl, na.action=na.exclude)

model5$data.newlyn.new <- data.frame(t=seq(from=1916,to=1952),
  model5$Pa[1:37,model5$tg.sd.sorted.newlyn.partial$ix[1:9]])

model5$tg.lmRob.newlyn.partial <- lmRob(msl ~ ., data=model5$data.newlyn.partial)
model5$tg.lmRob.newlyn.pred <- predict.lmRob(model5$tg.lmRob.newlyn.partial, se.fit=T, newdata=model5$data.newlyn.new, interval='confidence')

# Residual sum of squares
sum((model5$data.newlyn$msl[1:37]-model5$tg.lmRob.newlyn.pred$fit[1:37])^2)
#[1] 974070.3
sum((model5$data.newlyn$msl[38:91]-model5$tg.lmRob.newlyn.partial$fitted)^2)
#[1] 1210126


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

############
# Delfzijl #
############

# What we are going to do now is calculate the sd of the residuals as each of the
# components is added in turn 
model5$tg.sd.delfzijl <- vector(mode="numeric", length=(met$hadNUniq*3))

for(i in 1:(met$hadNUniq*3)){
  model5$data.delfzijl <- data.frame(msl=model5$delfzijlTotalP, t=seq(from=1916, to=2007), Pa=model5$Pa[1:92,i])
  model5$tg.lmRob.delfzijl <- lmRob(msl ~ t + Pa, data=model5$data.delfzijl, control=lmRobControl, na.action=na.exclude)
  
  model5$tg.sd.delfzijl[i] <- model5$tg.lmRob.delfzijl$scale
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model5$tg.sd.sorted.delfzijl <- sort(model5$tg.sd.delfzijl, index=TRUE)

# Now construct models of the top 9 components
model5$data.delfzijl <- data.frame(msl=model5$delfzijlTotalP,
  t=seq(from=1916,to=2007),
  model5$Pa[1:92,model5$tg.sd.sorted.delfzijl$ix[1:9]])

model5$tg.lmRob.delfzijl <- lmRob(msl ~ ., x=T, y=T, data=model5$data.delfzijl, control=lmRobControl, na.action=na.exclude)

# Test model prediction skill
model5$data.delfzijl.partial <- data.frame(msl=model5$delfzijlTotalP[38:92],
  t=seq(from=1953,to=2007),
  model5$Pa[38:92,model5$tg.sd.sorted.delfzijl$ix[1:9]])

model5$tg.lmRob.delfzijl.partial <- lmRob(msl ~ ., x=T, y=T, data=model5$data.delfzijl.partial, control=lmRobControl, na.action=na.exclude)

model5$data.delfzijl.new <- data.frame(t=seq(from=1916,to=1952),
  model5$Pa[1:37,model5$tg.sd.sorted.delfzijl$ix[1:9]])

model5$tg.lmRob.delfzijl.partial <- lmRob(msl ~ ., data=model5$data.delfzijl.partial)
model5$tg.lmRob.delfzijl.pred <- predict.lmRob(model5$tg.lmRob.delfzijl.partial, se.fit=T, newdata=model5$data.delfzijl.new, interval='confidence')

# Residual sum of squares
sum((model5$data.delfzijl$msl[1:28]-model5$tg.lmRob.delfzijl.pred$fit[1:28])^2)
#[1] 1918570 Brest: 3107529 Newlyn: 974070
sum((model5$data.delfzijl$msl[38:92]-model5$tg.lmRob.delfzijl.partial$fitted)^2)
#[1] 3514912 Brest: 4085193 Newlyn: 1210126



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

model6$tg.lmRob.brestTotalP <- lmRob(model6$brestTotalP ~ model6$time, control=lmRobControl, na.action=na.exclude)

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
  model6$tg.lmRob.brest <- lmRob(msl ~ t + Pa, data=model6$data.brest, control=lmRobControl, na.action=na.exclude)
  
  model6$tg.sd.brest[i] <- model6$tg.lmRob.brest$scale
  model6$tg.rsq.brest[i] <- model6$tg.lmRob.brest$r.sq
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model6$tg.sd.sorted.brest <- sort(model6$tg.sd.brest, index=TRUE)
model6$tg.rsq.sorted.brest <- sort(model6$tg.rsq.brest, index=TRUE, decreasing=TRUE)

# Now construct models of the top 9 components
model6$data.brest <- data.frame(msl=model6$brestTotalP,
  t=seq(from=1850,to=1943),
  Pa1=model6$Pa[,model6$tg.sd.sorted.brest$ix[1]], 
  Pa2=model6$Pa[,model6$tg.sd.sorted.brest$ix[2]], Pa3=model6$Pa[,model6$tg.sd.sorted.brest$ix[3]],
  Pa4=model6$Pa[,model6$tg.sd.sorted.brest$ix[4]], Pa5=model6$Pa[,model6$tg.sd.sorted.brest$ix[5]],
  Pa6=model6$Pa[,model6$tg.sd.sorted.brest$ix[6]], Pa7=model6$Pa[,model6$tg.sd.sorted.brest$ix[7]],
  Pa8=model6$Pa[,model6$tg.sd.sorted.brest$ix[8]], Pa9=model6$Pa[,model6$tg.sd.sorted.brest$ix[9]])

model6$tg.lmRob.brest <- lmRob(msl ~ t + Pa1 + Pa2 + Pa3 + Pa4 + Pa5 + Pa6 + Pa7 + Pa8 + Pa9, x=T, y=T,
                               data=model6$data.brest, control=lmRobControl, na.action=na.exclude)

model6$data.brest <- data.frame(msl=model6$brestTotalP,
  t=seq(from=1850,to=1943),
  Pa1=model6$Pa[,model6$tg.rsq.sorted.brest$ix[1]], 
  Pa2=model6$Pa[,model6$tg.rsq.sorted.brest$ix[2]], Pa3=model6$Pa[,model6$tg.rsq.sorted.brest$ix[3]],
  Pa4=model6$Pa[,model6$tg.rsq.sorted.brest$ix[4]], Pa5=model6$Pa[,model6$tg.rsq.sorted.brest$ix[5]],
  Pa6=model6$Pa[,model6$tg.rsq.sorted.brest$ix[6]], Pa7=model6$Pa[,model6$tg.rsq.sorted.brest$ix[7]],
  Pa8=model6$Pa[,model6$tg.rsq.sorted.brest$ix[8]], Pa9=model6$Pa[,model6$tg.rsq.sorted.brest$ix[9]])

model6$tg.lmRob.rsq.brest <- lmRob(msl ~ t + Pa1 + Pa2 + Pa3 + Pa4 + Pa5 + Pa6 + Pa7 + Pa8 + Pa9, x=T, y=T,
                               data=model6$data.brest, control=lmRobControl, na.action=na.exclude)

# Test model prediction skill
model6$data.brest.partial <- data.frame(msl=model6$brestTotalP[12:62],
  t=seq(from=1861,to=1911),
  model6$Pa[12:62,model6$tg.sd.sorted.brest$ix[1:9]])

model6$tg.lmRob.brest.partial <- lmRob(msl ~ ., x=T, y=T,
                               data=model6$data.brest.partial, control=lmRobControl, na.action=na.exclude)

model6$data.brest.new <- data.frame(t=seq(from=1912,to=1943),
  model6$Pa[63:94,model6$tg.sd.sorted.brest$ix[1:9]])

#model6$tg.lmRob.brest.partial <- lmRob(msl ~ ., data=model6$data.brest.partial)
model6$tg.lmRob.brest.pred <- predict.lmRob(model6$tg.lmRob.brest.partial, se.fit=T, newdata=model6$data.brest.new, interval='confidence')

# Residual sum of squares
sum((model6$data.brest$msl[1:37]-model6$tg.lmRob.brest.pred$fit[1:37])^2)
#[1] 974070.3
sum((model6$data.brest$msl[38:91]-model6$tg.lmRob.brest.partial$fitted)^2)
#[1] 1210126


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

#> model6$tg.lmRob.rsq.brest$r.sq
#[1] 0.4983623
#> coef(model6$tg.lmRob.rsq.brest)[2]/10.045
#       t 
#1.051769

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

model7$tg.lmRob.brestTotalP <- lmRob(model7$brestTotalP ~ model7$time, control=lmRobControl, na.action=na.exclude)

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
  model7$tg.lmRob.brest <- lmRob(msl ~ t + Pa, data=model7$data.brest, control=lmRobControl, na.action=na.exclude)
  
  model7$tg.sd.brest[i] <- model7$tg.lmRob.brest$scale
  model7$tg.rsq.brest[i] <- model7$tg.lmRob.brest$r.sq
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model7$tg.sd.sorted.brest <- sort(model7$tg.sd.brest, index=TRUE)
model7$tg.rsq.sorted.brest <- sort(model7$tg.rsq.brest, index=TRUE, decreasing=TRUE)

# Now construct models of the top 9 components
model7$data.brest <- data.frame(msl=model7$brestTotalP,
  t=seq(from=1888,to=1943),
  Pa1=model7$Pa[,model7$tg.sd.sorted.brest$ix[1]], 
  Pa2=model7$Pa[,model7$tg.sd.sorted.brest$ix[2]], Pa3=model7$Pa[,model7$tg.sd.sorted.brest$ix[3]],
  Pa4=model7$Pa[,model7$tg.sd.sorted.brest$ix[4]], Pa5=model7$Pa[,model7$tg.sd.sorted.brest$ix[5]],
  Pa6=model7$Pa[,model7$tg.sd.sorted.brest$ix[6]], Pa7=model7$Pa[,model7$tg.sd.sorted.brest$ix[7]],
  Pa8=model7$Pa[,model7$tg.sd.sorted.brest$ix[8]], Pa9=model7$Pa[,model7$tg.sd.sorted.brest$ix[9]])

model7$tg.lmRob.brest <- lmRob(msl ~ t + Pa1 + Pa2 + Pa3 + Pa4 + Pa5 + Pa6 + Pa7 + Pa8 + Pa9, x=T, y=T,
                               data=model7$data.brest, control=lmRobControl, na.action=na.exclude)

model7$data.brest <- data.frame(msl=model7$brestTotalP,
  t=seq(from=1888,to=1943),
  Pa1=model7$Pa[,model7$tg.rsq.sorted.brest$ix[1]], 
  Pa2=model7$Pa[,model7$tg.rsq.sorted.brest$ix[2]], Pa3=model7$Pa[,model7$tg.rsq.sorted.brest$ix[3]],
  Pa4=model7$Pa[,model7$tg.rsq.sorted.brest$ix[4]], Pa5=model7$Pa[,model7$tg.rsq.sorted.brest$ix[5]],
  Pa6=model7$Pa[,model7$tg.rsq.sorted.brest$ix[6]], Pa7=model7$Pa[,model7$tg.rsq.sorted.brest$ix[7]],
  Pa8=model7$Pa[,model7$tg.rsq.sorted.brest$ix[8]], Pa9=model7$Pa[,model7$tg.rsq.sorted.brest$ix[9]])

model7$tg.lmRob.rsq.brest <- lmRob(msl ~ t + Pa1 + Pa2 + Pa3 + Pa4 + Pa5 + Pa6 + Pa7 + Pa8 + Pa9, x=T, y=T,
                               data=model7$data.brest, control=lmRobControl, na.action=na.exclude)

# Variance reduction
#> (model7$tg.lmRob.brest$r.sq-model7$tg.lmRob.brestTotalP$r.sq)
#[1] 0.1658905
#> model7$tg.lmRob.brest$r.sq
#[1] 0.6909427
#>sqrt(diag(model7$tg.lmRob.brest$cov))[2]/10.045
#        t 
#0.4507884 
#> coef(model7$tg.lmRob.brest)[2]/10.045
#       t 
#2.336161

#> model7$tg.lmRob.rsq.brest$r.sq
#[1] 0.6821647
#> coef(model7$tg.lmRob.rsq.brest)[2]/10.045
#       t 
#2.468026 
#> sqrt(diag(model7$tg.lmRob.rsq.brest$cov))[2]/10.045
#        t 
#0.2279729 

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

model8$tg.lmRob.brestTotalP <- lmRob(model8$brestTotalP ~ model8$time, control=lmRobControl, na.action=na.exclude)

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
  model8$tg.lmRob.brest <- lmRob(msl ~ t + Pa, data=model8$data.brest, control=lmRobControl, na.action=na.exclude)
  
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
                               data=model8$data.brest, control=lmRobControl, na.action=na.exclude)

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

#############
## Model 9 ##
#############
# Scaled IB comparison with Brest as in Douglas 2008
# Take 1880-2008
model9<-new.env()

# Just use annual pressures over the periods of the model data (Pa) ending in 2000 (end of model pressure) 
model9$brestPa <-
  met$hadSLP2rAnnualPressureArray[intersect(which(met$hadSLP2rAnnualTime>=1880), which(met$hadSLP2rAnnualTime<=2008)),2]
model9$time <- c(1880:2008)
tg_annual$brest8008 <- tg_annual$brest[74:202,]

# Remove linear trend from both pressure and sea level
model9$tg.lmRob.brest <- lmRob(tg_annual$brest8008[,2] ~ model9$time, control=lmRobControl, x=T, y=T)
model9$pa.lmRob.brest <- lmRob(model9$brestPa ~ model9$time, control=lmRobControl, x=T, y=T)

# normalise both time series
model9$tg.norm.brest <- model9$tg.lmRob.brest$resid/max(abs(model9$tg.lmRob.brest$resid), na.rm=T)
model9$pa.norm.brest <- model9$pa.lmRob.brest$resid/max(abs(model9$pa.lmRob.brest$resid), na.rm=T)

# invert pressure?
model9$pa.norm.brest <- model9$pa.norm.brest*(-1)

# Find scale of pressure to sea level with lowest residual sd (best fit)
model9$pa.coperiod.brest.resid <- c(model9$pa.lmRob.brest$resid[1:64],model9$pa.lmRob.brest$resid[74:129])*-10
model9$pa.scaled.brest <- lmRob(model9$tg.lmRob.brest$resid ~ model9$pa.coperiod.brest.resid)
## summary(model9$pa.scaled.brest)

## Call: lmRob(formula = model9$tg.lmRob.brest$resid ~ model9$pa.coperiod.brest.resid)

## Residuals:
##        Min         1Q     Median         3Q        Max 
## -64.706984 -19.730825   1.093687  18.400327  85.948380 

## Coefficients:
##                                Value        Std. Error   t value     
## (Intercept)                    9.233435e-01 2.443849e+00 3.778235e-01
## model9$pa.coperiod.brest.resid 1.842819e+00 2.300845e-01 8.009312e+00
##                                Pr(>|t|)    
## (Intercept)                    7.062405e-01
## model9$pa.coperiod.brest.resid 9.143797e-13

## Residual standard error: 28.0952 on 118 degrees of freedom
## Multiple R-Squared: 0.328645 

## Test for Bias:
##              statistic   p-value
## M-estimate   0.1215612 0.9410297
## LS-estimate -2.3815424 1.0000000

#############
## Model 10 ##
#############
# Follow Thompson (1980) and see if we can improve
# Uses pressures from HadSLP2r
# Calculate over altimetry period 1992-2008

model10<-new.env()

# Just use annual pressures over the periods of the model data (Pa) ending in 2000 (end of model pressure) 
model10$newlynPa <-
  met$hadSLP2rAnnualPressureArray[intersect(which(met$hadSLP2rAnnualTime>=1992), which(met$hadSLP2rAnnualTime<=2008)),1]
model10$brestPa <-
  met$hadSLP2rAnnualPressureArray[intersect(which(met$hadSLP2rAnnualTime>=1992), which(met$hadSLP2rAnnualTime<=2008)),2]

model10$time <- c(1992:2008)
tg_annual$brest9208 <- tg_annual$brest[186:202,]
tg_annual$newlyn9208 <- tg_annual$newlyn[77:93,]

# Calculate the total pressure
model10$brestTotalP <- 1025*9.8*tg_annual$brest9208[,2]/1000 + model10$brestPa
model10$newlynTotalP <- 1025*9.8*tg_annual$newlyn9208[,2]/1000 + model10$newlynPa

model10$tg.lmRob.brestTotalP <- lmRob(model10$brestTotalP ~ model10$time, x=T, y=T, control=lmRobControl, na.action=na.exclude)
model10$tg.lmRob.newlynTotalP <- lmRob(model10$newlynTotalP ~ model10$time, x=T, y=T, control=lmRobControl, na.action=na.exclude)

model10$Pa <- array(NA,dim=c(17,(met$hadNUniq*3)))
# Zero lag
for(i in 1:met$hadNUniq){
  model10$Pa[,i] <- as.vector(met$hadSLP2rAnnualPressureArray[143:159,met$hadSLP2rUniq[i]])
}
# Lag 1 year
for(i in (met$hadNUniq+1):(met$hadNUniq*2)){
  model10$Pa[,i] <- as.vector(met$hadSLP2rAnnualPressureArray[142:158,met$hadSLP2rUniq[(i-met$hadNUniq)]])
}
# Lag 2 years....
for(i in ((met$hadNUniq*2)+1):(met$hadNUniq*3)){
  model10$Pa[,i] <- as.vector(met$hadSLP2rAnnualPressureArray[141:157,met$hadSLP2rUniq[(i-(met$hadNUniq*2))]])
}

#########
# Brest #
#########

# What we are going to do now is calculate the sd of the residuals as each of the
# components is added in turn 
model10$tg.sd.brest <- vector(mode="numeric", length=(met$hadNUniq*3))

for(i in 1:(met$hadNUniq*3)){
  model10$data.brest <- data.frame(msl=model10$brestTotalP, t=seq(from=1992,to=2008), Pa=model10$Pa[,i])
  model10$tg.lmRob.brest <- lmRob(msl ~ t + Pa, data=model10$data.brest, control=lmRobControl, na.action=na.exclude)
  
  model10$tg.sd.brest[i] <- model10$tg.lmRob.brest$scale
  model10$tg.rsq.brest[i] <- model10$tg.lmRob.brest$r.sq
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model10$tg.sd.sorted.brest <- sort(model10$tg.sd.brest, index=TRUE)
model10$tg.rsq.sorted.brest <- sort(model10$tg.rsq.brest, index=TRUE, decreasing=TRUE)

# Now construct models of the top 9 components
model10$data.brest <- data.frame(msl=model10$brestTotalP,
  t=seq(from=1992,to=2008),
  Pa1=model10$Pa[,model10$tg.sd.sorted.brest$ix[1]], 
  Pa2=model10$Pa[,model10$tg.sd.sorted.brest$ix[2]], Pa3=model10$Pa[,model10$tg.sd.sorted.brest$ix[3]],
  Pa4=model10$Pa[,model10$tg.sd.sorted.brest$ix[4]], Pa5=model10$Pa[,model10$tg.sd.sorted.brest$ix[5]],
  Pa6=model10$Pa[,model10$tg.sd.sorted.brest$ix[6]], Pa7=model10$Pa[,model10$tg.sd.sorted.brest$ix[7]],
  Pa8=model10$Pa[,model10$tg.sd.sorted.brest$ix[8]], Pa9=model10$Pa[,model10$tg.sd.sorted.brest$ix[9]])

model10$tg.lmRob.brest <- lmRob(msl ~ t + Pa1 + Pa2 + Pa3 + Pa4 + Pa5 + Pa6 + Pa7 + Pa8 + Pa9, x=T, y=T,
                               data=model10$data.brest, control=lmRobControl, na.action=na.exclude)
model10$tg.lmRob.brest.fitted.lm <- lmRob(model10$tg.lmRob.brest$fitted  ~ model10$tg.lmRob.brest$x[,2], control=lmRobControl, na.action=na.exclude)

model10$data.brest <- data.frame(msl=model10$brestTotalP,
  t=seq(from=1992,to=2008),
  Pa1=model10$Pa[,model10$tg.rsq.sorted.brest$ix[1]], 
  Pa2=model10$Pa[,model10$tg.rsq.sorted.brest$ix[2]], Pa3=model10$Pa[,model10$tg.rsq.sorted.brest$ix[3]],
  Pa4=model10$Pa[,model10$tg.rsq.sorted.brest$ix[4]], Pa5=model10$Pa[,model10$tg.rsq.sorted.brest$ix[5]],
  Pa6=model10$Pa[,model10$tg.rsq.sorted.brest$ix[6]], Pa7=model10$Pa[,model10$tg.rsq.sorted.brest$ix[7]],
  Pa8=model10$Pa[,model10$tg.rsq.sorted.brest$ix[8]], Pa9=model10$Pa[,model10$tg.rsq.sorted.brest$ix[9]])

model10$tg.lmRob.rsq.brest <- lmRob(msl ~ t + Pa1 + Pa2 + Pa3 + Pa4 + Pa5 + Pa6 + Pa7 + Pa8 + Pa9, x=T, y=T,
                               data=model10$data.brest, control=lmRobControl, na.action=na.exclude)

# Variance reduction
#> (model10$tg.lmRob.brest$r.sq-model10$tg.lmRob.brestTotalP$r.sq)
#[1] 0.7041578
#> model10$tg.lmRob.brest$r.sq
#[1] 0.8119005
#> model10$tg.lmRob.rsq.brest$r.sq
#[1] 0.8156477


##########
# Newlyn #
##########

# What we are going to do now is calculate the sd of the residuals as each of the
# components is added in turn 
model10$tg.sd.newlyn <- vector(mode="numeric", length=(met$hadNUniq*3))

for(i in 1:(met$hadNUniq*3)){
  model10$data.newlyn <- data.frame(msl=model10$newlynTotalP, t=seq(from=1992, to=2008), Pa=model10$Pa[,i])
  model10$tg.lmRob.newlyn <- lmRob(msl ~ t + Pa, data=model10$data.newlyn, control=lmRobControl, na.action=na.exclude)
  
  model10$tg.sd.newlyn[i] <- model10$tg.lmRob.newlyn$scale
  model10$tg.rsq.newlyn[i] <- model10$tg.lmRob.newlyn$r.sq
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model10$tg.sd.sorted.newlyn <- sort(model10$tg.sd.newlyn, index=TRUE)
model10$tg.rsq.sorted.newlyn <- sort(model10$tg.rsq.newlyn, index=TRUE, decreasing=TRUE)

# Now construct models of the top 9 components
model10$data.newlyn <- data.frame(msl=model10$newlynTotalP,
  t=seq(from=1992,to=2008),
  Pa1=model10$Pa[,model10$tg.sd.sorted.newlyn$ix[1]], 
  Pa2=model10$Pa[,model10$tg.sd.sorted.newlyn$ix[2]], Pa3=model10$Pa[,model10$tg.sd.sorted.newlyn$ix[3]],
  Pa4=model10$Pa[,model10$tg.sd.sorted.newlyn$ix[4]], Pa5=model10$Pa[,model10$tg.sd.sorted.newlyn$ix[5]],
  Pa6=model10$Pa[,model10$tg.sd.sorted.newlyn$ix[6]], Pa7=model10$Pa[,model10$tg.sd.sorted.newlyn$ix[7]],
  Pa8=model10$Pa[,model10$tg.sd.sorted.newlyn$ix[8]], Pa9=model10$Pa[,model10$tg.sd.sorted.newlyn$ix[9]])

model10$tg.lmRob.newlyn <- lmRob(msl ~ t + Pa1 + Pa2 + Pa3 + Pa4 + Pa5 + Pa6 + Pa7 + Pa8 + Pa9, x=T, y=T,
                               data=model10$data.newlyn, control=lmRobControl, na.action=na.exclude)
model10$tg.lmRob.newlyn.fitted.lm <- lmRob(model10$tg.lmRob.newlyn$fitted  ~ model10$tg.lmRob.newlyn$x[,2], control=lmRobControl, na.action=na.exclude)

model10$data.newlyn <- data.frame(msl=model10$newlynTotalP,
  t=seq(from=1992,to=2008),
  Pa1=model10$Pa[,model10$tg.rsq.sorted.newlyn$ix[1]], 
  Pa2=model10$Pa[,model10$tg.rsq.sorted.newlyn$ix[2]], Pa3=model10$Pa[,model10$tg.rsq.sorted.newlyn$ix[3]],
  Pa4=model10$Pa[,model10$tg.rsq.sorted.newlyn$ix[4]], Pa5=model10$Pa[,model10$tg.rsq.sorted.newlyn$ix[5]],
  Pa6=model10$Pa[,model10$tg.rsq.sorted.newlyn$ix[6]], Pa7=model10$Pa[,model10$tg.rsq.sorted.newlyn$ix[7]],
  Pa8=model10$Pa[,model10$tg.rsq.sorted.newlyn$ix[8]], Pa9=model10$Pa[,model10$tg.rsq.sorted.newlyn$ix[9]])

model10$tg.lmRob.rsq.newlyn <- lmRob(msl ~ t + Pa1 + Pa2 + Pa3 + Pa4 + Pa5 + Pa6 + Pa7 + Pa8 + Pa9, x=T, y=T,
                               data=model10$data.newlyn, control=lmRobControl, na.action=na.exclude)
# Variance reduction
#> (model10$tg.lmRob.newlyn$r.sq-model10$tg.lmRob.newlynTotalP$r.sq)
#[1] 0.3892037
#> model10$tg.lmRob.newlyn$r.sq
#[1] 0.8096695
#c(coef(model10$tg.lmRob.brest)[2]/10.045, coef(model10$tg.lmRob.newlyn)[2]/10.045)
#       t        t 
#2.560348 4.363938
#c(sqrt(diag(model10$tg.lmRob.brest$cov))[2]/10.045,  sqrt(diag(model10$tg.lmRob.newlyn$cov))[2]/10.045)
#        t         t 
#0.6579228 0.2529815

#> model10$tg.lmRob.rsq.newlyn$r.sq
#[1] 0.7694162
#c(coef(model10$tg.lmRob.rsq.brest)[2]/10.045, coef(model10$tg.lmRob.rsq.newlyn)[2]/10.045)
#       t        t 
#1.482763 3.274131
#c(sqrt(diag(model10$tg.lmRob.rsq.brest$cov))[2]/10.045, sqrt(diag(model10$tg.lmRob.rsq.newlyn$cov))[2]/10.045)
#        t         t
#0.4943032 0.4423348

##############
## Model 11 ##
##############
# Follow Thompson (1980) and see if we can improve
# Uses pressures from HadSLP2r
# Calculate over Rossiter period of 1894-1961, just for Brest ;)

model11<-new.env()

# Just use annual pressures over the periods of the model data (Pa) ending in 2000 (end of model pressure) 
model11$brestPa <-
  met$hadSLP2rAnnualPressureArray[intersect(which(met$hadSLP2rAnnualTime>=1894), which(met$hadSLP2rAnnualTime<=1961)),2]
model11$time <- c(1894:1961)
tg_annual$brest9461 <- tg_annual$brest[88:155,]

# Calculate the total pressure
model11$brestTotalP <- 1025*9.8*tg_annual$brest9461[,2]/1000 + model11$brestPa

model11$tg.lmRob.brestTotalP <- lmRob(model11$brestTotalP ~ model11$time, control=lmRobControl, na.action=na.exclude)
model11$tg.lm.brestTotalP <- lm(model11$brestTotalP ~ model11$time)

model11$Pa <- array(NA,dim=c(68,(met$hadNUniq*3)))
# Zero lag
for(i in 1:met$hadNUniq){
  model11$Pa[,i] <- as.vector(met$hadSLP2rAnnualPressureArray[45:112,met$hadSLP2rUniq[i]])
}
# Lag 1 year
for(i in (met$hadNUniq+1):(met$hadNUniq*2)){
  model11$Pa[,i] <- as.vector(met$hadSLP2rAnnualPressureArray[44:111,met$hadSLP2rUniq[(i-met$hadNUniq)]])
}
# Lag 2 years....
for(i in ((met$hadNUniq*2)+1):(met$hadNUniq*3)){
  model11$Pa[,i] <- as.vector(met$hadSLP2rAnnualPressureArray[43:110,met$hadSLP2rUniq[(i-(met$hadNUniq*2))]])
}

#########
# Brest #
#########

# What we are going to do now is calculate the sd of the residuals as each of the
# components is added in turn 
model11$tg.sd.brest <- vector(mode="numeric", length=(met$hadNUniq*3))

for(i in 1:(met$hadNUniq*3)){
  model11$data.brest <- data.frame(msl=model11$brestTotalP, t=seq(from=1894,to=1961), Pa=model11$Pa[,i])
  model11$tg.lmRob.brest <- lmRob(msl ~ t + Pa, data=model11$data.brest, control=lmRobControl, na.action=na.exclude)
  
  model11$tg.sd.brest[i] <- model11$tg.lmRob.brest$scale
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model11$tg.sd.sorted.brest <- sort(model11$tg.sd.brest, index=TRUE)

# Now construct models of the top 9 components
model11$data.brest <- data.frame(msl=model11$brestTotalP,
  t=seq(from=1894,to=1961),
  Pa1=model11$Pa[,model11$tg.sd.sorted.brest$ix[1]], 
  Pa2=model11$Pa[,model11$tg.sd.sorted.brest$ix[2]], Pa3=model11$Pa[,model11$tg.sd.sorted.brest$ix[3]],
  Pa4=model11$Pa[,model11$tg.sd.sorted.brest$ix[4]], Pa5=model11$Pa[,model11$tg.sd.sorted.brest$ix[5]],
  Pa6=model11$Pa[,model11$tg.sd.sorted.brest$ix[6]], Pa7=model11$Pa[,model11$tg.sd.sorted.brest$ix[7]],
  Pa8=model11$Pa[,model11$tg.sd.sorted.brest$ix[8]], Pa9=model11$Pa[,model11$tg.sd.sorted.brest$ix[9]])

model11$tg.lmRob.brest <- lmRob(msl ~ t + Pa1 + Pa2 + Pa3 + Pa4 + Pa5 + Pa6 + Pa7 + Pa8 + Pa9, x=T, y=T,
                               data=model11$data.brest, control=lmRobControl, na.action=na.exclude)

# Variance reduction
#> (model11$tg.lmRob.brest$r.sq-model11$tg.lmRob.brestTotalP$r.sq)
#[1] 0.2202793
#> model11$tg.lmRob.brest$r.sq
#[1] 0.6940627
#>sqrt(diag(model11$tg.lmRob.brest$cov))[2]/10.045
#        t 
#0.3093761 
#> coef(model11$tg.lmRob.brest)[2]/10.045
#       t 
#1.911946

#############
## Model 12 ##
#############
# Scaled IB comparison with Brest as in Douglas 2008 using CRU data
# Take 1880-2008
model12<-new.env()

# Just use annual pressures over the periods of the model data (Pa) ending in 2000 (end of model pressure) 
model12$brestPa <-
  met$cru.brest$SLP[intersect(which(met$cru.brest$Year>=1880), which(met$cru.brest$Year<=2000))]
model12$time <- c(1880:2000)
tg_annual$brest8000 <- tg_annual$brest[74:194,]

# Remove linear trend from both pressure and sea level
model12$tg.lmRob.brest <- lmRob(tg_annual$brest8000[,2] ~ model12$time, control=lmRobControl, x=T, y=T)
model12$pa.lmRob.brest <- lmRob(model12$brestPa ~ model12$time, control=lmRobControl, x=T, y=T)

# normalise both time series
model12$tg.norm.brest <- model12$tg.lmRob.brest$resid/max(abs(model12$tg.lmRob.brest$resid), na.rm=T)
model12$pa.norm.brest <- model12$pa.lmRob.brest$resid/max(abs(model12$pa.lmRob.brest$resid), na.rm=T)
# Actually, Bruce uses the ratio of SDs between the two time series as his nomralisation, not the absolute range. With robust statistics,
# this is the "scale" parameter
model12$tg.norm.scale.brest <- model12$tg.lmRob.brest$resid/model12$tg.lmRob.brest$scale
# In addition to this scaling we need an additional factor of 10 to convert pressure in Pascals to mm of sea level
model12$pa.norm.scale.brest <- model12$pa.lmRob.brest$resid*10/model12$pa.lmRob.brest$scale

# invert pressure?
model12$pa.norm.brest <- model12$pa.norm.brest*(-1)
model12$pa.norm.scale.brest <- model12$pa.norm.scale.brest*(-1)

# Find scale of pressure to sea level with lowest residual sd (best fit)
model12$pa.coperiod.brest.resid <- c(model12$pa.lmRob.brest$resid[1:64],model12$pa.lmRob.brest$resid[74:121])*-10
model12$pa.scaled.brest <- lmRob(model12$tg.lmRob.brest$resid ~ model12$pa.coperiod.brest.resid)
## summary(model12$pa.scaled.brest)

## Call: lmRob(formula = model12$tg.lmRob.brest$resid ~ model12$pa.coperiod.brest.resid)
##
##Residuals:
##       Min         1Q     Median         3Q        Max 
##-70.324366 -20.254001   3.726153  19.478394  83.294073 
##
##Coefficients:
##                                Value        Std. Error   t value     
##(Intercept)                     9.797562e-02 2.714828e+00 3.608907e-02
##model12$pa.coperiod.brest.resid 1.433893e+00 2.053923e-01 6.981237e+00
##                                Pr(>|t|)    
##(Intercept)                     9.712768e-01
##model12$pa.coperiod.brest.resid 2.330798e-10
##
##Residual standard error: 30.5129 on 110 degrees of freedom
##Multiple R-Squared: 0.28764 
##
##Test for Bias:
##            statistic   p-value
##M-estimate   1.652864 0.4376079
##LS-estimate -6.329955 1.0000000

##############
## Model 13 ##
##############
# Follow Thompson (1979) with N Atl EOFs and see if we can improve
# Uses pressures from HadSLP2r and EOFs calculated in slpEOF.R
# Calculate over 1916-2008 period

model13<-new.env()

# Calculate correlation coefficients with EOFs.

# Just use annual pressures over the periods of the model data (Pa) ending in 2000 (end of model pressure) 
model13$newlynPa <-
  met$hadSLP2rAnnualPressureArray[intersect(which(met$hadSLP2rAnnualTime>=1916), which(met$hadSLP2rAnnualTime<=2008)),1]
model13$brestPa <-
  met$hadSLP2rAnnualPressureArray[intersect(which(met$hadSLP2rAnnualTime>=1916), which(met$hadSLP2rAnnualTime<=2008)),2]

model13$time <- c(1916:2008)

# Calculate the total pressure
model13$brestTotalP <- 1025*9.8*tg_annual$brest1608[,2]/1000 + model13$brestPa
model13$brestTotalP_22<- 1025*9.8*tg_annual$brest1608_22[,2]/1000 + model13$brestPa
model13$newlynTotalP <- 1025*9.8*tg_annual$newlyn1608[,2]/1000 + model13$newlynPa

model13$tg.lmRob.brestTotalP <- lmRob(model13$brestTotalP ~ model13$time, x=T, y=T, control=lmRobControl, na.action=na.exclude)
model13$tg.lmRob.brestTotalP_22<- lmRob(model13$brestTotalP_22 ~ model13$time, control=lmRobControl, na.action=na.exclude)
model13$tg.lmRob.newlynTotalP <- lmRob(model13$newlynTotalP ~ model13$time, x=T, y=T, control=lmRobControl, na.action=na.exclude)

# First 14 EOFs explain 96.8% of the annual variance

model13$EOF <- array(NA,dim=c(93,14))
# Zero lag
for(i in 1:14){
  model13$EOF[,i] <- (slp$dmSlpAnnual %*% slp$svdDmSlpAnnual$v[,i])[67:159]
}
# Lag 1 year
#for(i in (met$hadNUniq+1):(met$hadNUniq*2)){
#  model13$Pa[,i] <- as.vector(met$hadSLP2rAnnualPressureArray[66:158,met$hadSLP2rUniq[(i-met$hadNUniq)]])
#}
# Lag 2 years....
#for(i in ((met$hadNUniq*2)+1):(met$hadNUniq*3)){
#  model13$Pa[,i] <- as.vector(met$hadSLP2rAnnualPressureArray[65:157,met$hadSLP2rUniq[(i-(met$hadNUniq*2))]])
#}
# Calculate correlation array. What is the relationship between each of the pressure series and the sea level time-series?
# Columns are Brest, Newlyn then EOFs
model13$corr.array1608 <- array(NA,dim=c(16,16))
model13$corr.data1608 <- cbind(model13$brestTotalP,model13$newlynTotalP,model13$EOF)
for(i in 1:16){
  for(j in 1:16){
    model13$corr.array1608[i,j] <- cor.test(model13$corr.data1608[,i], model13$corr.data1608[,j], method="p", alternative="t", na.action=na.fail)$estimate
  }
}
# First half
model13$corr.array1643 <- array(NA,dim=c(16,16))
# Model 4 has the total pressures at Brest and Newlyn for 1916-1943
model13$corr.data1643 <- cbind(model4$brestTotalP,model4$newlynTotalP,model13$EOF[1:28,])
for(i in 1:16){
  for(j in 1:16){
    model13$corr.array1643[i,j] <- cor.test(model13$corr.data1643[,i], model13$corr.data1643[,j], method="p", alternative="t", na.action=na.fail)$estimate
  }
}
# Second half
model13$corr.array5308 <- array(NA,dim=c(16,16))
# Model 3 has the total pressures at Brest and Newlyn for 1953-2008
model13$corr.data5308 <- cbind(model3$brestTotalP,model3$newlynTotalP,model13$EOF[38:93,])
for(i in 1:16){
  for(j in 1:16){
    model13$corr.array5308[i,j] <- cor.test(model13$corr.data5308[,i], model13$corr.data5308[,j], method="p", alternative="t", na.action=na.fail)$estimate
  }
}

#########
# Brest #
#########

# Construct models of the 14 components
model13$data.brest <- data.frame(msl=model13$brestTotalP,
  t=seq(from=1916,to=2008),
  EOF1=model13$EOF[,1], EOF2=model13$EOF[,2], EOF3=model13$EOF[,3], EOF4=model13$EOF[,4], EOF5=model13$EOF[,5],
  EOF6=model13$EOF[,6], EOF7=model13$EOF[,7], EOF8=model13$EOF[,8], EOF9=model13$EOF[,9], EOF10=model13$EOF[,10],
  EOF11=model13$EOF[,11], EOF12=model13$EOF[,12], EOF13=model13$EOF[,13], EOF14=model13$EOF[,14])

model13$tg.lmRob.brest <- lmRob(msl ~ t + EOF1 + EOF2 + EOF3 + EOF4 + EOF5 + EOF6 + EOF7 + EOF8 + EOF9 + EOF10
                                + EOF11 + EOF12 + EOF13 + EOF14,
                                x=T, y=T, data=model13$data.brest, control=lmRobControl, na.action=na.exclude)
model13$tg.lmRob.brest.fitted.lm <- lmRob(model13$tg.lmRob.brest$fitted  ~ model13$tg.lmRob.brest$x[,2], control=lmRobControl, na.action=na.exclude)

# Test model prediction skill
model13$data.brest.partial <- data.frame(msl=model13$brestTotalP[38:93],
  t=seq(from=1953,to=2008),
  model13$EOF[38:93,1:14])

model13$tg.lmRob.brest.partial <- lmRob(msl ~ ., x=T, y=T, data=model13$data.brest.partial, control=lmRobControl, na.action=na.exclude)

model13$data.brest.new <- data.frame(t=seq(from=1916,to=1952),
  model13$EOF[1:37,1:14])

model13$tg.lmRob.brest.partial <- lmRob(msl ~ ., data=model13$data.brest.partial)
model13$tg.lmRob.brest.pred <- predict.lmRob(model13$tg.lmRob.brest.partial, se.fit=T, newdata=model13$data.brest.new, interval='confidence')

# Residual sum of squares
sum((model13$data.brest$msl[1:28]-model13$tg.lmRob.brest.pred$fit[1:28])^2)
#[1] 4399047 Newlyn: EOF: 2142474 Pa: 974070
sum((model13$data.brest$msl[38:93]-model13$tg.lmRob.brest.partial$fitted)^2)
#[1] 1153239 Newlyn: EOF: 1299343 Pa: 1210126

# Variance reduction
#> (model13$tg.lmRob.brest$r.sq-model13$tg.lmRob.brestTotalP$r.sq)
#[1] 0.1222487
#> model13$tg.lmRob.brest$r.sq
#[1] 0.670531


#################
# Brest plus 22 #
#################

# Construct models of the 14 components
model13$data.brest <- data.frame(msl=model13$brestTotalP_22,
  t=seq(from=1916,to=2008),
  EOF1=model13$EOF[,1], EOF2=model13$EOF[,2], EOF3=model13$EOF[,3], EOF4=model13$EOF[,4], EOF5=model13$EOF[,5],
  EOF6=model13$EOF[,6], EOF7=model13$EOF[,7], EOF8=model13$EOF[,8], EOF9=model13$EOF[,9], EOF10=model13$EOF[,10],
  EOF11=model13$EOF[,11], EOF12=model13$EOF[,12], EOF13=model13$EOF[,13], EOF14=model13$EOF[,14])

model13$tg.lmRob.brest_22<- lmRob(msl ~ t + EOF1 + EOF2 + EOF3 + EOF4 + EOF5 + EOF6 + EOF7 + EOF8 + EOF9 + EOF10
                                + EOF11 + EOF12 + EOF13 + EOF14,
                                x=T, y=T, data=model13$data.brest, control=lmRobControl, na.action=na.exclude)

# Variance reduction
#> (model13$tg.lmRob.brest_22$r.sq-model13$tg.lmRob.brestTotalP_22$r.sq)
#[1] 0.0847231
#> model13$tg.lmRob.brest_22$r.sq
#[1] 0.7114546

##########
# Newlyn #
##########

# Construct models of the 14 components
model13$data.newlyn <- data.frame(msl=model13$newlynTotalP,
  t=seq(from=1916,to=2008),
  model13$EOF[,1:14])

model13$tg.lmRob.newlyn<- lmRob(msl ~ ., x=T, y=T, data=model13$data.newlyn, control=lmRobControl, na.action=na.exclude)

model13$tg.lmRob.newlyn.fitted.lm <- lmRob(model13$tg.lmRob.newlyn$fitted  ~ model13$tg.lmRob.newlyn$x[,2], control=lmRobControl, na.action=na.exclude)

## Test model prediction skill
#model13$data.newlyn.partial <- data.frame(msl=model13$newlynTotalP[38:91],
#  t=seq(from=1953,to=2006),
#  model13$EOF[38:91,1:14])
#
#model13$tg.lmRob.newlyn.partial <- lmRob(msl ~ ., x=T, y=T, data=model13$data.newlyn.partial, control=lmRobControl,
#  na.action=na.exclude)
#
#model13$data.newlyn.new <- data.frame(t=seq(from=1916,to=1952),
#  model13$EOF[1:37,1:14])

#model13$tg.lmRob.newlyn.partial <- lmRob(msl ~ ., data=model13$data.newlyn.partial)
#model13$tg.lmRob.newlyn.pred <- predict.lmRob(model13$tg.lmRob.newlyn.partial, se.fit=T,
#  newdata=model13$data.newlyn.new, interval='confidence')

# Test model prediction skill
# Need to build new set of top 9 components for partial set
# What we are going to do now is calculate the sd of the residuals as each of the
# components is added in turn 
model13$tg.sd.newlyn.partial <- vector(mode="numeric", length=(met$hadNUniq*3))

for(i in 1:(met$hadNUniq*3)){
  model13$data.newlyn.partial <- data.frame(msl=model13$newlynTotalP[38:91], t=seq(from=1953,to=2006), Pa=model13$EOF[38:91,i])
  model13$tg.lmRob.newlyn.partial <- lmRob(msl ~ ., data=model13$data.newlyn.partial, control=lmRobControl, na.action=na.exclude)
  
  model13$tg.sd.newlyn.partial[i] <- model13$tg.lmRob.newlyn.partial$scale
}

# Now we can sort these sds. 
model13$tg.sd.sorted.newlyn.partial <- sort(model13$tg.sd.newlyn.partial, index=TRUE)

model13$data.newlyn.partial <- data.frame(msl=model13$newlynTotalP[38:91], t=seq(from=1953,to=2006),
  model13$Pa[38:91,model13$tg.sd.sorted.newlyn.partial$ix[1:9]])

model13$tg.lmRob.newlyn.partial <- lmRob(msl ~ ., x=T, y=T, data=model13$data.newlyn.partial, control=lmRobControl, na.action=na.exclude)

model13$data.newlyn.new <- data.frame(t=seq(from=1916,to=1952),
  model13$Pa[1:37,model13$tg.sd.sorted.newlyn.partial$ix[1:9]])

model13$tg.lmRob.newlyn.partial <- lmRob(msl ~ ., data=model13$data.newlyn.partial)
model13$tg.lmRob.newlyn.pred <- predict.lmRob(model13$tg.lmRob.newlyn.partial, se.fit=T, newdata=model13$data.newlyn.new, interval='confidence')

# Residual sum of squares
sum((model13$data.newlyn$msl[1:37]-model13$tg.lmRob.newlyn.pred$fit[1:37])^2)
#[1] 2142474 Pa: 974070
sum((model13$data.newlyn$msl[38:91]-model13$tg.lmRob.newlyn.partial$fitted)^2)
#[1] 1299343 Pa: 1210126

# Variance reduction
#> (model13$tg.lmRob.newlyn$r.sq-model13$tg.lmRob.newlynTotalP$r.sq)
#[1] 0.08409328
#> model13$tg.lmRob.newlyn$r.sq
#[1] 0.7344319
#c(coef(model13$tg.lmRob.brest)[2]/10.045, coef(model13$tg.lmRob.newlyn)[2]/10.045)
#       t        t 
#1.423981 1.590407
#c(sqrt(diag(model13$tg.lmRob.brest$cov))[2]/10.045,  sqrt(diag(model13$tg.lmRob.newlyn$cov))[2]/10.045)
#        t         t 
#0.2496875 0.1107558
