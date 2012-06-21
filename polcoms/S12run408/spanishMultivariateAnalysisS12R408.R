###############################################################################
#
# spanishMultivariateAnalysis.R
# 
# Implementation of Thompson's methods (Thompson, 1980; 1986) to the tide gauge
# time series at Spanish Atlantic stations
#
# Author: simonh
###############################################################################
model<-new.env()
load("~/diskx/polcoms/S12run408/spain/spanishModelMonthlyMean.RData", envir=model)

interp<-new.env()
load("~/diskx/polcoms/spain/spanishMetMonthlyMeanInterp.RData", envir=interp)

# Model to be fitted
# p*g*n + pa = a*tx + b*ty + c11*cos(w1*t) + c12*sin(w1*t) 
# 				+ c21*cos(w2*t) + c22*sin(w2*t) + e
# where w1 = 2*pi/12 month and w2 = 2*pi/6 month, pa = local air pressure, 
# (p*g*n + pa) is the total pressure, (tx, ty) is influence of local wind stress
# modelled by (a*tx + b*ty) where (a,b) are regression coefficients to be
# determined from the data

model$numYears <- 45
model$time <- seq.Date(from=as.Date("1960/1/15"), to=as.Date("2004/12/15"), by="1 month")
model$t <- seq(from=0,to=((model$numYears*12)-1))

# Generate annual means from model data
model$annualYrs <- seq.Date(from=as.Date("1960/7/1"), to=as.Date("2004/7/1"), by="1 year")
model$annualMeans <- array(NA,dim=c(model$numYears,12))
for (i in 1:12){
  model$tmp <- model$spanishMonthlyMean[,i]
  dim(model$tmp) <- c(12,model$numYears)
  model$annualMeans[,i] <- colMeans(model$tmp)*1000
}
rm(tmp, envir=model)

x11()
par(family="HersheySans")
plot(model$annualYrs, model$annualMeans[,11]-mean(model$annualMeans[,11], na.rm=T)-100, type='l', col='blue', ann=F, ylim=c(-300,200), lwd=2)
lines(model$annualYrs, model$annualMeans[,8]-mean(model$annualMeans[,8], na.rm=T)-50, col='red', lwd=2)
lines(model$annualYrs, model$annualMeans[,7]-mean(model$annualMeans[,7], na.rm=T), col='magenta', lwd=2)
lines(model$annualYrs, model$annualMeans[,3]-mean(model$annualMeans[,3], na.rm=T)+50, col='orange', lwd=2) 
legend(x=1980, y=-150,
  legend=c('Santander I','Coruna I','Coruna II', 'Vigo'),
  col=c('blue', 'red', 'magenta', 'orange'), lwd=2)
title(ylab='Height [mm]', xlab='Year')
#dev2bitmap(file="modelAnnualDataAtlSpain.png", res=150)

# Calculate running 10 year means for the Atlantic stations
model$tenYrTrends <- array(NA,dim=c((model$numYears-9),12))
model$midpointYrs <- seq.Date(from=as.Date("1965/1/1"), length=(model$numYears-9), by="1 year")
for (i in 1:(model$numYears-9)) {
  for ( j in 1:12) {
     if (length(which(is.finite(model$annualMeans[i:(i+9),j])))>=7){
       model$junk <- lm(model$annualMeans[i:(i+9),j] ~ c(1:10), na.action = na.exclude)
       model$tenYrTrends[i,j] <- model$junk$coef[2]
     } else {
       model$tenYrTrends[i,j] <- NA
     }
  }
}

model$meanTenYrTrend <- rowMeans(model$tenYrTrends, na.rm=T)
model$meanTenYrTrendLongRecs <- rowMeans(model$tenYrTrends[,c(11,8,7,3)], na.rm=T)


x11()
#par(family="HersheySans")

plot(model$midpointYrs, model$tenYrTrends[,11]-15, type='l', col='blue',
  ann=F, ylim=c(-40,30), lwd=2)
title(ylab='Rate [mm/yr]', xlab='Year')
lines(model$midpointYrs, model$tenYrTrends[,8]-5, col='red', lwd=2)
lines(model$midpointYrs, model$tenYrTrends[,7], col='magenta', lwd=2)
lines(model$midpointYrs, model$tenYrTrends[,3]+10, col='orange', lwd=2)
lines(model$midpointYrs, model$meanTenYrTrendLongRecs, lwd=2, lty=4)

legend(x=as.Date("1965/7/1"), y=-22,
  legend=c('Santander I','Coruna I','Coruna II', 'Vigo', 'Mean'),
  col=c('blue', 'red', 'magenta', 'orange', 'black'), lwd=2,
  lty=c(1,1,1,1,4))
#dev2bitmap(file="modelDecadalratesAtlSpain.png", res=150)


model$totalP <- array(NA, dim=c(length(model$spanishMonthlyMean[,1]),12))
for (i in 1:12){
  model$totalP[,i] <- 
    1025*9.8*model$spanishMonthlyMean[,i] + as.vector(interp$spanishMetMonthlyMeanInterp[,1,,i])
}

w1 <- 2*pi/12
w2 <- 2*pi/6
  
rainbow_colours <- rainbow(12) 

x11()
par(family="HersheySans")

model1<-new.env()
model1$resid <- array(NA,dim=c(480,12))
model1$ds.totalP <- array(NA,dim=c(480,12))

model$sd.totalP <- vector(mode='numeric', length=12)
model1$sd.ds.totalP <- vector(mode='numeric', length=12)
model1$sd.resid <- vector(mode='numeric', length=12)

for (i in 1:12){
  model1$data <- data.frame(totP=model$totalP[,i], 
    tauX=as.vector(interp$spanishMetMonthlyMeanInterp[,4,,i]), 
    tauY=as.vector(interp$spanishMetMonthlyMeanInterp[,5,,i]), 
    t=model$t, w1=w1, w2=w2)

  model1$lm <- lm(totP ~ t + tauX + tauY + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)), 
    data=model1$data)
  model1$resid[,i] <- model1$lm$resid

# Deseasonalise the data
  model1$ds.totalP[,i] <- model$totalP[,i]-
    (model1$lm$coef[4]*cos(w1 * model$t)+model1$lm$coef[5]*sin(w1 * model$t))-
    (model1$lm$coef[6]*cos(w2 * model$t)+model1$lm$coef[7]*sin(w2 * model$t))+
    model1$lm$coef[1] 

  model$sd.totalP[i] <- sd(model$totalP[,i]/(1025*9.8), na.rm=T)
  model1$sd.ds.totalP[i] <- sd(model1$ds.totalP[,12]/(1025*9.8), na.rm=T)
  model1$sd.resid[i] <- sd(model1$resid[,i]/(1025*9.8), na.rm=T)
  
  if (i == 1){
    plot(model$time, model1$resid[,i], type='l', col=rainbow_colours[i])
  } else {
    lines(model$time, model1$resid[,i], col=rainbow_colours[i])
  }
}


#> range((model1$lm$coef[4]*cos(w1 * model1$data$t)+model1$lm$coef[5]*sin(w1 * model1$data$t)))[2]/(1025*9.8)
#[1] 0.0453652
# i.e. the mean seasonal cycle at Vigo II has an amplitude of 4.5cm cf 4cm at Newlyn in Thompson 1986
# sd(model$totalP[,12]/(1025*9.8), na.rm=T)
#[1] 0.0635112
# sd(model1$ds.totalP[,12]/(1025*9.8), na.rm=T)
# [1] 0.0481164 (i.e. 6.2cm cf 7cm in Thompson, 1986)
# sd(model1$resid[,12]/(1025*9.8), na.rm=T)
#[1] 0.04056846 (i.e. 5.4cm cf 3.1cm in Thompson, 1986)

# Compare with TG data
tg<-new.env()
load("~/diskx/polcoms/spain/spanishData.RData", envir=tg)
load("~/diskx/polcoms/spain/atlanticSpain.RData", envir=tg)

tg$resid <- array(NA,dim=c(480,12))
tg$ds.totalP <- array(NA,dim=c(480,12))

model$sd.totalP <- vector(mode='numeric', length=12)
tg$sd.ds.totalP <- vector(mode='numeric', length=12)
tg$sd.resid <- vector(mode='numeric', length=12)

tg$time <- seq.Date(from=as.Date("1944/1/15"), to=as.Date("2007/12/15"), by="1 month")
# Compare just over the period of the model
tg$mp <- array(NA,dim=c(480,12))
# Correlations
tg$cor <- vector("numeric", 12)
for (i in 1:12){
 tg$mp[,i] <- tg$spanishMonthly[intersect(which(tg$time>=as.Date("1960/1/15")), which(tg$time<=as.Date("1999/12/15"))),i]
 tg$junk <- cor.test(tg$mp[,i],(model$spanishMonthlyMean[,i]*1000), alternative="greater", na.action="na.fail")
 tg$cor[i] <- tg$junk$estimate
}

x11()
par(family="HersheySans")
plot(model$time, tg$mp[,11], type='l', col='orange')
lines(model$time,tg$mp[,12], col='blue')

x11()
par(family="HersheySans")
plot(tg$mp[,11],(model$spanishMonthlyMean[,11]*1000))

x11()
par(family="HersheySans")
plot(tg$mp[,7],(model$spanishMonthlyMean[,7]*1000))

# Calculate sds
tg$sd <- vector(mode='numeric', length=12)
model$sd <- vector(mode='numeric', length=12)

for (i in 1:12){
  tg$sd[i] <- sd(tg$mp[,i], na.rm=T)
  model$sd[i] <- sd(model$spanishMonthlyMean[,i], na.rm=T)
}

tg$totalP <- array(NA, dim=dim(tg$mp))

x11()
par(family="HersheySans")

for (i in 1:12){
  tg$totalP[,i] <- 1025*9.8*tg$mp[,i]/1000 + as.vector(interp$spanishMetMonthlyMeanInterp[,1,,i])

  tg$data <- data.frame(totP=tg$totalP[,i], tauX=as.vector(interp$spanishMetMonthlyMeanInterp[,4,,i]), 
    tauY=as.vector(interp$spanishMetMonthlyMeanInterp[,5,,i]), t=model$t, w1=w1, w2=w2)
  
  tg$lm <- lm(totP ~ t + tauX + tauY + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)), 
    data=tg$data, na.action=na.exclude, x=T)
  tg$resid[(as.vector(tg$lm$x[,2])+1),i] <- tg$lm$resid

# Deseasonalise the data
  tg$ds.totalP[,i] <- model$totalP[,i]-
    (tg$lm$coef[4]*cos(w1 * model$t)+tg$lm$coef[5]*sin(w1 * model$t))-
    (tg$lm$coef[6]*cos(w2 * model$t)+tg$lm$coef[7]*sin(w2 * model$t))+
    tg$lm$coef[1] 

  tg$sd.ds.totalP[i] <- sd(tg$ds.totalP[,1]/(1025*9.8), na.rm=T)
  tg$sd.resid[i] <- sd(tg$resid[,1], na.rm=T)
  
  if (i == 1){
    plot(model$time, tg$resid[,i], type='l', col=rainbow_colours[i])
  } else {
    lines(model$time, tg$resid[,i], col=rainbow_colours[i])
  }

}

# Annual signal amplitude
# range(tg$lm.brest$coef[4]*cos(w1 * tg$data.newlyn$t)+tg$lm.brest$coef[5]*sin(w1 * tg$data.newlyn$t))
# 392.3397/(1025*9.8)
#[1] 0.03905821
# i.e. the mean seasonal cycle has an amplitude of 3.9cm cf 4cm at Newlyn in Thompson 1986

# sd(tg$totalPNewlyn/(1025*9.8), na.rm=T)
#[1] 0.0817145
# sd(tg$ds.totalPNewlyn/(1025*9.8), na.rm=T)
#[1] 0.06992447 (i.e. 7cm cf 7cm in Thompson, 1986)
# sd(tg$resid.newlyn/(1025*9.8), na.rm=T)
#[1] 0.05948289 (i.e. 5.9cm cf 3.1cm in Thompson, 1986)

# Thompson's 4th model is rather harder to implement
# 
model4 <- new.env() 

model4$m.resid <- array(NA,dim=c(480,12))
model4$tg.resid <- array(NA,dim=c(480,12))
model4$ds.totalP <- array(NA,dim=c(480,12))

model$sd.totalP <- vector(mode='numeric', length=12)
model4$sd.ds.totalP <- vector(mode='numeric', length=12)
model4$sd.resid <- vector(mode='numeric', length=12)

# Again, just use pressures over the periods of the model data (Pm)

# In Thompson's model 4 he also introduces lags of 0, 1 & 2 months which in this
# case gives a total of 51 components. However, to facilitate comparison with
# model 3 he only uses the first 9 components which he chooses by stepwise
# regression. We achieve this by adding in each component in turn and calculating
# the variance reduction, which we can then rank and analyze.
model4$Pm <- array(NA,dim=c(480,51))
# Zero lag
for(i in 1:17){
  model4$Pm[,i] <- as.vector(interp$extraMetMonthlyMeanInterp[,1,,i])
}
# Lag 1
for(i in 18:34){
  model4$Pm[2:480,i] <- as.vector(interp$extraMetMonthlyMeanInterp[,1,,(i-17)])[1:479]
}
# Lag 2
for(i in 35:51){
  model4$Pm[3:480,i] <- as.vector(interp$extraMetMonthlyMeanInterp[,1,,(i-34)])[1:478]
}

# What we are going to do now is calculate the sd of the residuals as each of the
# 51 components is added in turn 
model4$m.sd <- array(NA, dim=c(51,12))
model4$tg.sd <- array(NA, dim=c(51,12))

for(i in 1:51){
  for (j in 1:12){
    model4$data <- data.frame(mmsl=model$spanishMonthlyMean[,j], msl=tg$mp[,j],
      t=model$t, w1=w1, w2=w2, Pm=model4$Pm[,i])

    model4$m.lm <- lm(mmsl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) +
      I(sin(w2*t)) + Pm, data=model4$data, x=T)
    model4$tg.lm <- lm(msl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + 
      I(sin(w2*t)) + Pm, data=model4$data, x=T)
  
    model4$m.resid[(as.vector(model4$m.lm$x[,2])+1),j] <- model4$m.lm$resid
    model4$tg.resid[(as.vector(model4$tg.lm$x[,2])+1),j] <- model4$tg.lm$resid
  
    model4$m.sd[i,j] <- sd(model4$m.resid[,j], na.rm=T)
    model4$tg.sd[i,j] <- sd(model4$tg.resid[,j], na.rm=T)
  }
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
x11()
par(family="HersheySans")

for (i in 1:12){
  model4$m.sd.sorted <- sort(model4$m.sd[,i], index=TRUE)
  model4$tg.sd.sorted <- sort(model4$tg.sd[,i], index=TRUE)


# Now construct models of the top 9 components
  model4$data <- data.frame(mmsl=model$spanishMonthlyMean[,i], msl=tg$mp[,i],
    t=model$t, w1=w1, w2=w2, Pm1=model4$Pm[,model4$m.sd.sorted$ix[1]], 
    Pm2=model4$Pm[,model4$m.sd.sorted$ix[2]], Pm3=model4$Pm[,model4$m.sd.sorted$ix[3]],
    Pm4=model4$Pm[,model4$m.sd.sorted$ix[4]], Pm5=model4$Pm[,model4$m.sd.sorted$ix[5]],
    Pm6=model4$Pm[,model4$m.sd.sorted$ix[6]], Pm7=model4$Pm[,model4$m.sd.sorted$ix[7]],
    Pm8=model4$Pm[,model4$m.sd.sorted$ix[8]], Pm9=model4$Pm[,model4$m.sd.sorted$ix[9]])

  model4$m.lm <- lm(mmsl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)) + 
    Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9,  data=model4$data)

  model4$data <- data.frame(mmsl=model$spanishMonthlyMean[,i], msl=tg$mp[,i],
    t=model$t, w1=w1, w2=w2, Pm1=model4$Pm[,model4$tg.sd.sorted$ix[1]], 
    Pm2=model4$Pm[,model4$tg.sd.sorted$ix[2]], Pm3=model4$Pm[,model4$tg.sd.sorted$ix[3]],
    Pm4=model4$Pm[,model4$tg.sd.sorted$ix[4]], Pm5=model4$Pm[,model4$tg.sd.sorted$ix[5]],
    Pm6=model4$Pm[,model4$tg.sd.sorted$ix[6]], Pm7=model4$Pm[,model4$tg.sd.sorted$ix[7]],
    Pm8=model4$Pm[,model4$tg.sd.sorted$ix[8]], Pm9=model4$Pm[,model4$tg.sd.sorted$ix[9]])

  model4$tg.lm <- lm(msl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)) + 
    Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9,  data=model4$data, x=T)

  model4$m.resid[,i] <- model4$m.lm$resid
  model4$tg.resid[(as.vector(model4$tg.lm$x[,2])+1),i] <- model4$tg.lm$resid
  
# Deseasonalise the data
  model4$ds.totalP[,i] <- model$totalP[,i]-
    (model4$m.lm$coef[4]*cos(w1 * model$t)+model4$m.lm$coef[5]*sin(w1 * model$t))-
    (model4$m.lm$coef[6]*cos(w2 * model$t)+model4$m.lm$coef[7]*sin(w2 * model$t))+
    model4$m.lm$coef[1] 

  model4$sd.ds.totalP[i] <- sd(model4$ds.totalP[,i]/(1025*9.8), na.rm=T)
  model4$sd.resid[i] <- sd(model4$m.resid[,i], na.rm=T)
  
  if (i == 1){
    plot(model$time, model4$m.resid[,i], type='l', col=rainbow_colours[i])
  } else {
    lines(model$time, model4$m.resid[,i], col=rainbow_colours[i])
  }
}

## After choosing i=7
#save(file="corunaI.RData",data,envir=model4)

#> sd(model4$data.newlyn$mmsl, na.rm=T)
#[1] 0.0676091
#> sd(model4$m.resid.newlyn, na.rm=T)
#[1] 0.02422649

#> sd(model4$data.newlyn$msl, na.rm=T)
#[1] 82.10498
#> sd(model4$tg.resid.newlyn, na.rm=T)
#[1] 25.16201
# No signficant improvement

## Now the new stuff. Add deep sea boundary condition.
#model5 <- new.env()
#
#model5$data.newlyn <- data.frame(mmsl=newlynMonthlyMean, msl=tg$mp.newlyn,
#  t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2, Pm1=model4$Pm[,model4$m.sd.sorted$ix[1]], 
#  Pm2=model4$Pm[,model4$m.sd.sorted$ix[2]], Pm3=model4$Pm[,model4$m.sd.sorted$ix[3]],
#  Pm4=model4$Pm[,model4$m.sd.sorted$ix[4]], Pm5=model4$Pm[,model4$m.sd.sorted$ix[5]],
#  Pm6=model4$Pm[,model4$m.sd.sorted$ix[6]], Pm7=model4$Pm[,model4$m.sd.sorted$ix[7]],
#  Pm8=model4$Pm[,model4$m.sd.sorted$ix[8]], Pm9=model4$Pm[,model4$m.sd.sorted$ix[9]],
#  dssl=deepSeaMonthlyMean)
#
#model5$m.lm.newlyn <- lm(mmsl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)) + 
#  Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 + dssl,  data=model5$data.newlyn)
#
#model5$m.resid.newlyn <- model5$m.lm.newlyn$resid
#model5$tg.resid.newlyn <- model5$tg.lm.newlyn$resid
#
##> sd(model5$data.newlyn$mmsl, na.rm=T)
##[1] 0.0676091
##> sd(model5$m.resid.newlyn, na.rm=T)
##[1] 0.02422649
#
##> sd(model5$data.newlyn$msl, na.rm=T)
##[1] 82.10498
##> sd(model5$tg.resid.newlyn, na.rm=T)
##[1] 25.16201

#############
## Model 6 ##
#############
# Focus on Coruna I as the longest record

# Use winds from the model instead of pressures
# Treat the 18 (corunaI + 17 others) wind points in the same way as we did the pressure
# points with 3 lags for the 17 distant stations (=51 combinations) but
# just use zero lag for corunaI
# Es is the E wind stress, Ns is the N wind stress
# 
model6 <- new.env() 
# Again, just use winds over the periods of the model data (Pm)

# As with the pressures we use only the first 9 components which we choose by stepwise
# regression. Again, we achieve this by adding in each component in turn and calculating
# the variance reduction, which we can then rank and analyze.
model6$Ns <- array(NA,dim=c(480,51))
model6$Es <- array(NA,dim=c(480,51))

model6$corunaIEs <- as.vector(interp$spanishMetMonthlyMeanInterp[,2,,7])
model6$corunaINs <- as.vector(interp$spanishMetMonthlyMeanInterp[,3,,7])
model6$corunaIPm <- as.vector(interp$spanishMetMonthlyMeanInterp[,1,,7])

# Zero lag
for(i in 1:17){
  model6$Es[,i] <- as.vector(interp$extraMetMonthlyMeanInterp[,2,,i])
  model6$Ns[,i] <- as.vector(interp$extraMetMonthlyMeanInterp[,3,,i])
}
# Lag 1
for(i in 18:34){
  model6$Es[2:480,i] <- as.vector(interp$extraMetMonthlyMeanInterp[,2,,(i-17)])[1:479]
  model6$Ns[2:480,i] <- as.vector(interp$extraMetMonthlyMeanInterp[,3,,(i-17)])[1:479]
}
# Lag 2
for(i in 35:51){
  model6$Es[3:480,i] <- as.vector(interp$extraMetMonthlyMeanInterp[,2,,(i-34)])[1:478]
  model6$Ns[3:480,i] <- as.vector(interp$extraMetMonthlyMeanInterp[,3,,(i-34)])[1:478]
}

# What we are going to do now is calculate the sd of the residuals as each of the
# 51 components is added in turn 
model6$m.sd <- vector(mode="numeric", length=51)
model6$tg.sd <- vector(mode="numeric", length=51)

for(i in 1:51){
  model6$data.corunaI <- data.frame(mmsl=model$spanishMonthlyMean[,7], msl=tg$mp[,7],
    t=model$t, w1=w1, w2=w2, Pm=model6$corunaIPm,
    NEs=model6$corunaIEs, NNs=model6$corunaINs, Es=model6$Es[,i], Ns=model6$Ns[,i])

  model6$m.lm.corunaI <- lm(mmsl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) +
    I(sin(w2*t)) + Es + Ns + NEs + NNs + Pm, data=model6$data.corunaI)
  model6$tg.lm.corunaI <- lm(msl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + 
    I(sin(w2*t)) + Es + Ns + NEs + NNs + Pm, data=model6$data.corunaI)
  
  model6$m.resid.corunaI <- model6$m.lm.corunaI$resid
  model6$tg.resid.corunaI <- model6$tg.lm.corunaI$resid
  
  model6$m.sd[i] <- sd(model6$m.resid.corunaI, na.rm=T)
  model6$tg.sd[i] <- sd(model6$tg.resid.corunaI, na.rm=T)
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model6$m.sd.sorted <- sort(model6$m.sd, index=TRUE)
model6$tg.sd.sorted <- sort(model6$tg.sd, index=TRUE)

# Now construct models of the top 9 components
#model6$t <- seq(from=0,to=((40*12)-1))
model6$data.corunaI <- data.frame(mmsl=model$spanishMonthlyMean[,7], 
  msl=tg$mp[,7], t=model$t, w1=w1, w2=w2, 
  Pm1=model4$Pm[,model4$m.sd.sorted$ix[1]], 
  Pm2=model4$Pm[,model4$m.sd.sorted$ix[2]], Pm3=model4$Pm[,model4$m.sd.sorted$ix[3]],
  Pm4=model4$Pm[,model4$m.sd.sorted$ix[4]], Pm5=model4$Pm[,model4$m.sd.sorted$ix[5]],
  Pm6=model4$Pm[,model4$m.sd.sorted$ix[6]], Pm7=model4$Pm[,model4$m.sd.sorted$ix[7]],
  Pm8=model4$Pm[,model4$m.sd.sorted$ix[8]], Pm9=model4$Pm[,model4$m.sd.sorted$ix[9]],
  Es1=model6$Es[,model6$m.sd.sorted$ix[1]], 
  Es2=model6$Es[,model6$m.sd.sorted$ix[2]], Es3=model6$Es[,model6$m.sd.sorted$ix[3]],
  Es4=model6$Es[,model6$m.sd.sorted$ix[4]], Es5=model6$Es[,model6$m.sd.sorted$ix[5]],
  Es6=model6$Es[,model6$m.sd.sorted$ix[6]], Es7=model6$Es[,model6$m.sd.sorted$ix[7]],
  Es8=model6$Es[,model6$m.sd.sorted$ix[8]], Es9=model6$Es[,model6$m.sd.sorted$ix[9]],
  Ns1=model6$Ns[,model6$m.sd.sorted$ix[1]], 
  Ns2=model6$Ns[,model6$m.sd.sorted$ix[2]], Ns3=model6$Ns[,model6$m.sd.sorted$ix[3]],
  Ns4=model6$Ns[,model6$m.sd.sorted$ix[4]], Ns5=model6$Ns[,model6$m.sd.sorted$ix[5]],
  Ns6=model6$Ns[,model6$m.sd.sorted$ix[6]], Ns7=model6$Ns[,model6$m.sd.sorted$ix[7]],
  Ns8=model6$Ns[,model6$m.sd.sorted$ix[8]], Ns9=model6$Ns[,model6$m.sd.sorted$ix[9]],
  NEs=model6$corunaIEs, NNs=model6$corunaINs, Pm=model6$corunaIPm,
  w1cos=cos(w1*model$t), w1sin=sin(w1*model$t), 
  w2cos=cos(w2*model$t), w2sin=sin(w2*model$t))

model6$m.lm.corunaI <- lm(mmsl ~ t + w1cos + w1sin + w2cos + w2sin + 
  Es1 + Es2 + Es3 + Es4 + Es5 + Es6 + Es7 + Es8 + Es9 +
  Ns1 + Ns2 + Ns3 + Ns4 + Ns5 + Ns6 + Ns7 + Ns8 + Ns9 +
  NEs + NNs + Pm,  data=model6$data.corunaI)

model6$data.corunaI <- data.frame(mmsl=model$spanishMonthlyMean[,7], msl=tg$mp[,7],
  t=model$t, w1=w1, w2=w2, 
  Pm1=model4$Pm[,model4$tg.sd.sorted$ix[1]], 
  Pm2=model4$Pm[,model4$tg.sd.sorted$ix[2]], Pm3=model4$Pm[,model4$tg.sd.sorted$ix[3]],
  Pm4=model4$Pm[,model4$tg.sd.sorted$ix[4]], Pm5=model4$Pm[,model4$tg.sd.sorted$ix[5]],
  Pm6=model4$Pm[,model4$tg.sd.sorted$ix[6]], Pm7=model4$Pm[,model4$tg.sd.sorted$ix[7]],
  Pm8=model4$Pm[,model4$tg.sd.sorted$ix[8]], Pm9=model4$Pm[,model4$tg.sd.sorted$ix[9]],
  Es1=model6$Es[,model6$tg.sd.sorted$ix[1]], 
  Es2=model6$Es[,model6$tg.sd.sorted$ix[2]], Es3=model6$Es[,model6$tg.sd.sorted$ix[3]],
  Es4=model6$Es[,model6$tg.sd.sorted$ix[4]], Es5=model6$Es[,model6$tg.sd.sorted$ix[5]],
  Es6=model6$Es[,model6$tg.sd.sorted$ix[6]], Es7=model6$Es[,model6$tg.sd.sorted$ix[7]],
  Es8=model6$Es[,model6$tg.sd.sorted$ix[8]], Es9=model6$Es[,model6$tg.sd.sorted$ix[9]],
  Ns1=model6$Ns[,model6$tg.sd.sorted$ix[1]], 
  Ns2=model6$Ns[,model6$tg.sd.sorted$ix[2]], Ns3=model6$Ns[,model6$tg.sd.sorted$ix[3]],
  Ns4=model6$Ns[,model6$tg.sd.sorted$ix[4]], Ns5=model6$Ns[,model6$tg.sd.sorted$ix[5]],
  Ns6=model6$Ns[,model6$tg.sd.sorted$ix[6]], Ns7=model6$Ns[,model6$tg.sd.sorted$ix[7]],
  Ns8=model6$Ns[,model6$tg.sd.sorted$ix[8]], Ns9=model6$Ns[,model6$tg.sd.sorted$ix[9]],
  NEs=model6$corunaIEs, NNs=model6$corunaINs, Pm=model6$corunaIPm,
  w1cos=cos(w1*model$t), w1sin=sin(w1*model$t), 
  w2cos=cos(w2*model$t), w2sin=sin(w2*model$t))

model6$tg.lm.corunaI <- lm(msl ~ t + w1cos + w1sin + w2cos + w2sin + 
  Es1 + Es2 + Es3 + Es4 + Es5 + Es6 + Es7 + Es8 + Es9 +
  Ns1 + Ns2 + Ns3 + Ns4 + Ns5 + Ns6 + Ns7 + Ns8 + Ns9 +
  NEs + NNs + Pm,  data=model6$data.corunaI, x=T, y=T)

model6$m.resid.corunaI <- model6$m.lm.corunaI$resid
model6$tg.resid.corunaI <- model6$tg.lm.corunaI$resid

#> sd(model6$data.corunaI$mmsl, na.rm=T)
#[1] 0.06452842
#> sd(model6$m.resid.corunaI, na.rm=T)
#[1] 0.02891333

#> sd(model6$data.corunaI$msl, na.rm=T)
#[1] 84.9157
#> sd(model6$tg.resid.corunaI, na.rm=T)
#[1] 40.14192

#> model6$tg.lm.corunaI$coef[2]*12
#       t 
#2.368088 

#> lm(tg$spanishAnnual[17:56,7] ~ c(1:40), na.action=na.exclude)
#
#Call:
#lm(formula = tg$spanishAnnual[17:56, 7] ~ c(1:40), na.action = na.exclude)
#
#Coefficients:
#(Intercept)      c(1:40)  
#   6998.529        2.017  

# Calculate running 10 year means for the Atlantic stations
# Put back the missing vals
model6$tg.resid.full <- vector('numeric',480)
model6$tg.resid.full[as.vector(model6$tg.lm.corunaI$x[,2])] <- model6$tg.resid.corunaI
model6$tg.fit.full <- vector('numeric',480)
model6$tg.fit.full[as.vector(model6$tg.lm.corunaI$x[,2])] <- model6$tg.lm.corunaI$fit
model6$tg.fit.full[which(model6$tg.fit.full==0)] <- NA
# Make annual means
dim(model6$tg.resid.full)<-c(12,40)
model6$tg.resid.annual <- colMeans(model6$tg.resid.full, na.rm=T)

model6$tenYrTrends <- vector('numeric',length=(40-9))
model6$midpointYrs <- seq.Date(from=as.Date("1965/1/15"), length=(40-9), by="1 year")
for (i in 1:(40-9)) {
   if (length(which(is.finite(model6$tg.resid.annual[i:(i+9)])))>=7){
     model6$junk <- lm(model6$tg.resid.annual[i:(i+9)] ~ c(1:10), na.action = na.exclude)
     model6$tenYrTrends[i] <- model6$junk$coef[2]
   } else {
     model6$tenYrTrends[i] <- NA
   }
}

x11()
par(family="HersheySans")

plot(model$time, model6$data.corunaI$msl - 
  mean(model6$data.corunaI$msl, na.rm=T), type='l', col='blue',
  ann=F, lwd=2)
lines(model$time, 
  as.vector(model6$tg.fit.full) - mean(as.vector(model6$tg.fit.full), 
  na.rm=T), col='red', lwd=2)
lines(model$time, as.vector(model6$tg.resid.full), col='magenta', lwd=2)

title(ylab='Height [mm]', xlab='Year')
legend(x=as.Date("1972/1/1"), y=400,
  legend=c('Original, sd=84.9 mm','Fit','Residual, sd=40.1 mm'),
  col=c('blue', 'red','magenta'), lwd=2)
#dev2bitmap(file="corunaIMonthlyMVFit.png", res=150)

x11()
par(family="HersheySans")

plot(model6$midpointYrs, model6$tenYrTrends, type='l', col='blue',
  ann=F, ylim=c(-10,15), lwd=2)
lines(tg$midpointYrs, tg$tenYrTrends[,7], col='red', lwd=2)

title(ylab='Rate [mm/yr]', xlab='Year')
legend(x=1972, y=15,
  legend=c('Residual','Original'),
  col=c('blue', 'red'), lwd=2)
#dev2bitmap(file="corunaIDecadalRatesAtlSpain.png", res=150)

 
save(file="corunaI.RData",list=ls(model6), envir=model6)

#############
## Model 7 ##
#############
# Focus on Santander I as another long record

# Use winds from the model instead of pressures
# Treat the 18 (santanderI + 17 others) wind points in the same way as we did the pressure
# points with 3 lags for the 17 distant stations (=51 combinations) but
# just use zero lag for santanderI
# Es is the E wind stress, Ns is the N wind stress
# 
model7 <- new.env() 
# Again, just use winds over the periods of the model data (Pm)

# As with the pressures we use only the first 9 components which we choose by stepwise
# regression. Again, we achieve this by adding in each component in turn and calculating
# the variance reduction, which we can then rank and analyze.
model7$Ns <- array(NA,dim=c(480,51))
model7$Es <- array(NA,dim=c(480,51))

model7$santanderIEs <- as.vector(interp$spanishMetMonthlyMeanInterp[,2,,3])
model7$santanderINs <- as.vector(interp$spanishMetMonthlyMeanInterp[,3,,3])
model7$santanderIPm <- as.vector(interp$spanishMetMonthlyMeanInterp[,1,,3])

# Zero lag
for(i in 1:17){
  model7$Es[,i] <- as.vector(interp$extraMetMonthlyMeanInterp[,2,,i])
  model7$Ns[,i] <- as.vector(interp$extraMetMonthlyMeanInterp[,3,,i])
}
# Lag 1
for(i in 18:34){
  model7$Es[2:480,i] <- as.vector(interp$extraMetMonthlyMeanInterp[,2,,(i-17)])[1:479]
  model7$Ns[2:480,i] <- as.vector(interp$extraMetMonthlyMeanInterp[,3,,(i-17)])[1:479]
}
# Lag 2
for(i in 35:51){
  model7$Es[3:480,i] <- as.vector(interp$extraMetMonthlyMeanInterp[,2,,(i-34)])[1:478]
  model7$Ns[3:480,i] <- as.vector(interp$extraMetMonthlyMeanInterp[,3,,(i-34)])[1:478]
}

# What we are going to do now is calculate the sd of the residuals as each of the
# 51 components is added in turn 
model7$m.sd <- vector(mode="numeric", length=51)
model7$tg.sd <- vector(mode="numeric", length=51)

for(i in 1:51){
  model7$data.santanderI <- data.frame(mmsl=model$spanishMonthlyMean[,3], msl=tg$mp[,3],
    t=model$t, w1=w1, w2=w2, Pm=model7$santanderIPm,
    NEs=model7$santanderIEs, NNs=model7$santanderINs, Es=model7$Es[,i], Ns=model7$Ns[,i])

  model7$m.lm.santanderI <- lm(mmsl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) +
    I(sin(w2*t)) + Es + Ns + NEs + NNs + Pm, data=model7$data.santanderI)
  model7$tg.lm.santanderI <- lm(msl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + 
    I(sin(w2*t)) + Es + Ns + NEs + NNs + Pm, data=model7$data.santanderI)
  
  model7$m.resid.santanderI <- model7$m.lm.santanderI$resid
  model7$tg.resid.santanderI <- model7$tg.lm.santanderI$resid
  
  model7$m.sd[i] <- sd(model7$m.resid.santanderI, na.rm=T)
  model7$tg.sd[i] <- sd(model7$tg.resid.santanderI, na.rm=T)
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model7$m.sd.sorted <- sort(model7$m.sd, index=TRUE)
model7$tg.sd.sorted <- sort(model7$tg.sd, index=TRUE)

# Now construct models of the top 9 components
#model7$t <- seq(from=0,to=((40*12)-1))
model7$data.santanderI <- data.frame(mmsl=model$spanishMonthlyMean[,3], 
  msl=tg$mp[,3], t=model$t, w1=w1, w2=w2, 
  Pm1=model4$Pm[,model4$m.sd.sorted$ix[1]], 
  Pm2=model4$Pm[,model4$m.sd.sorted$ix[2]], Pm3=model4$Pm[,model4$m.sd.sorted$ix[3]],
  Pm4=model4$Pm[,model4$m.sd.sorted$ix[4]], Pm5=model4$Pm[,model4$m.sd.sorted$ix[5]],
  Pm6=model4$Pm[,model4$m.sd.sorted$ix[6]], Pm7=model4$Pm[,model4$m.sd.sorted$ix[7]],
  Pm8=model4$Pm[,model4$m.sd.sorted$ix[8]], Pm9=model4$Pm[,model4$m.sd.sorted$ix[9]],
  Es1=model7$Es[,model7$m.sd.sorted$ix[1]], 
  Es2=model7$Es[,model7$m.sd.sorted$ix[2]], Es3=model7$Es[,model7$m.sd.sorted$ix[3]],
  Es4=model7$Es[,model7$m.sd.sorted$ix[4]], Es5=model7$Es[,model7$m.sd.sorted$ix[5]],
  Es6=model7$Es[,model7$m.sd.sorted$ix[6]], Es7=model7$Es[,model7$m.sd.sorted$ix[7]],
  Es8=model7$Es[,model7$m.sd.sorted$ix[8]], Es9=model7$Es[,model7$m.sd.sorted$ix[9]],
  Ns1=model7$Ns[,model7$m.sd.sorted$ix[1]], 
  Ns2=model7$Ns[,model7$m.sd.sorted$ix[2]], Ns3=model7$Ns[,model7$m.sd.sorted$ix[3]],
  Ns4=model7$Ns[,model7$m.sd.sorted$ix[4]], Ns5=model7$Ns[,model7$m.sd.sorted$ix[5]],
  Ns6=model7$Ns[,model7$m.sd.sorted$ix[6]], Ns7=model7$Ns[,model7$m.sd.sorted$ix[7]],
  Ns8=model7$Ns[,model7$m.sd.sorted$ix[8]], Ns9=model7$Ns[,model7$m.sd.sorted$ix[9]],
  NEs=model7$santanderIEs, NNs=model7$santanderINs, Pm=model7$santanderIPm,
  w1cos=cos(w1*model$t), w1sin=sin(w1*model$t), 
  w2cos=cos(w2*model$t), w2sin=sin(w2*model$t))

model7$m.lm.santanderI <- lm(mmsl ~ t + w1cos + w1sin + w2cos + w2sin + 
  Es1 + Es2 + Es3 + Es4 + Es5 + Es6 + Es7 + Es8 + Es9 +
  Ns1 + Ns2 + Ns3 + Ns4 + Ns5 + Ns6 + Ns7 + Ns8 + Ns9 +
  NEs + NNs + Pm,  data=model7$data.santanderI)

model7$data.santanderI <- data.frame(mmsl=model$spanishMonthlyMean[,3], msl=tg$mp[,3],
  t=model$t, w1=w1, w2=w2, 
  Pm1=model4$Pm[,model4$tg.sd.sorted$ix[1]], 
  Pm2=model4$Pm[,model4$tg.sd.sorted$ix[2]], Pm3=model4$Pm[,model4$tg.sd.sorted$ix[3]],
  Pm4=model4$Pm[,model4$tg.sd.sorted$ix[4]], Pm5=model4$Pm[,model4$tg.sd.sorted$ix[5]],
  Pm6=model4$Pm[,model4$tg.sd.sorted$ix[6]], Pm7=model4$Pm[,model4$tg.sd.sorted$ix[7]],
  Pm8=model4$Pm[,model4$tg.sd.sorted$ix[8]], Pm9=model4$Pm[,model4$tg.sd.sorted$ix[9]],
  Es1=model7$Es[,model7$tg.sd.sorted$ix[1]], 
  Es2=model7$Es[,model7$tg.sd.sorted$ix[2]], Es3=model7$Es[,model7$tg.sd.sorted$ix[3]],
  Es4=model7$Es[,model7$tg.sd.sorted$ix[4]], Es5=model7$Es[,model7$tg.sd.sorted$ix[5]],
  Es6=model7$Es[,model7$tg.sd.sorted$ix[6]], Es7=model7$Es[,model7$tg.sd.sorted$ix[7]],
  Es8=model7$Es[,model7$tg.sd.sorted$ix[8]], Es9=model7$Es[,model7$tg.sd.sorted$ix[9]],
  Ns1=model7$Ns[,model7$tg.sd.sorted$ix[1]], 
  Ns2=model7$Ns[,model7$tg.sd.sorted$ix[2]], Ns3=model7$Ns[,model7$tg.sd.sorted$ix[3]],
  Ns4=model7$Ns[,model7$tg.sd.sorted$ix[4]], Ns5=model7$Ns[,model7$tg.sd.sorted$ix[5]],
  Ns6=model7$Ns[,model7$tg.sd.sorted$ix[6]], Ns7=model7$Ns[,model7$tg.sd.sorted$ix[7]],
  Ns8=model7$Ns[,model7$tg.sd.sorted$ix[8]], Ns9=model7$Ns[,model7$tg.sd.sorted$ix[9]],
  NEs=model7$santanderIEs, NNs=model7$santanderINs, Pm=model7$santanderIPm,
  w1cos=cos(w1*model$t), w1sin=sin(w1*model$t), 
  w2cos=cos(w2*model$t), w2sin=sin(w2*model$t))

model7$tg.lm.santanderI <- lm(msl ~ t + w1cos + w1sin + w2cos + w2sin + 
  Es1 + Es2 + Es3 + Es4 + Es5 + Es6 + Es7 + Es8 + Es9 +
  Ns1 + Ns2 + Ns3 + Ns4 + Ns5 + Ns6 + Ns7 + Ns8 + Ns9 +
  NEs + NNs + Pm,  data=model7$data.santanderI, x=T, y=T)

model7$m.resid.santanderI <- model7$m.lm.santanderI$resid
model7$tg.resid.santanderI <- model7$tg.lm.santanderI$resid
model7$tg.fit.full <- vector('numeric',480)
model7$tg.fit.full[as.vector(model7$tg.lm.santanderI$x[,2])] <- model7$tg.lm.santanderI$fit
model7$tg.fit.full[which(model7$tg.fit.full==0)] <- NA

#> sd(model7$data.santanderI$mmsl, na.rm=T)
#[1] 0.05944343
#> sd(model7$m.resid.santanderI, na.rm=T)
#[1] 0.02726072

#> sd(model7$data.santanderI$msl, na.rm=T)
#[1] 78.91515
#> sd(model7$tg.resid.santanderI, na.rm=T)
#[1] 37.74835

#> model7$tg.lm.santanderI$coef[2]*12
#       t 
#1.62323 

#> lm(tg$spanishAnnual[17:56,7] ~ c(1:40), na.action=na.exclude)
#
#Call:
#lm(formula = tg$spanishAnnual[17:56, 7] ~ c(1:40), na.action = na.exclude)
#
#Coefficients:
#(Intercept)      c(1:40)  
#   6998.529        2.017  

# Calculate running 10 year means for the Atlantic stations
# Put back the missing vals
model7$tg.resid.full <- vector('numeric',480)
model7$tg.resid.full[as.vector(model7$tg.lm.santanderI$x[,2])] <- model7$tg.resid.santanderI
# Make annual means
dim(model7$tg.resid.full)<-c(12,40)
model7$tg.resid.annual <- colMeans(model7$tg.resid.full, na.rm=T)

model7$tenYrTrends <- vector('numeric',length=(40-9))
model7$midpointYrs <- seq.Date(from=as.Date("1965/1/15"), length=(40-9), by="1 year")
for (i in 1:(40-9)) {
   if (length(which(is.finite(model7$tg.resid.annual[i:(i+9)])))>=7){
     model7$junk <- lm(model7$tg.resid.annual[i:(i+9)] ~ c(1:10), na.action = na.exclude)
     model7$tenYrTrends[i] <- model7$junk$coef[2]
   } else {
     model7$tenYrTrends[i] <- NA
   }
}

x11()
par(family="HersheySans")

plot(model$time, model7$data.santanderI$msl - 
  mean(model7$data.santanderI$msl, na.rm=T), type='l', col='blue',
  ann=F, lwd=2, ylim=c(-200,320))
lines(model$time, 
  as.vector(model7$tg.fit.full) - mean(as.vector(model7$tg.fit.full), 
  na.rm=T), col='red', lwd=2)
lines(model$time, as.vector(model7$tg.resid.full), col='magenta', lwd=2)

title(ylab='Height [mm]', xlab='Year')
legend(x=as.Date("1960/1/1"), y=320,
  legend=c('Original, sd=78.9 mm','Fit','Residual, sd=37.7 mm'),
  col=c('blue', 'red','magenta'), lwd=2)
#dev2bitmap(file="santanderIMonthlyMVFit.png", res=150)

x11()
par(family="HersheySans")

plot(model7$midpointYrs, model7$tenYrTrends, type='l', col='blue',
  ann=F, ylim=c(-40,30), lwd=2)
lines(tg$midpointYrs, tg$tenYrTrends[,3], col='red', lwd=2)
title(ylab='Rate [mm/yr]', xlab='Year')
legend(x=1975, y=-10,
  legend=c('Residual','Original'),
  col=c('blue', 'red'), lwd=2)

#dev2bitmap(file="santanderIDecadalRatesAtlSpain.png", res=150)

save(file="santanderI.RData", list=ls(model7), envir=model7)

# So wind stresses (or just winds) actually do worse than just static pressures. Why?
# Dynamic height in the model is not the same as sea level so I should diagnose and use a
# linear combination of diagnostic heights first.
# Jason also suggests using the currents. Although this isn't a very useful tool
# for reconstructing sea level (unless we have a model) it is useful as a diagnostic tool
# to tell us what is actually controlling SL at the coast.
