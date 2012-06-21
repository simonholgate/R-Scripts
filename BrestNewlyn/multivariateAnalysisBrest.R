# multivariateAnalysisBrest.R
# 
# Implementation of Thompson's methods (Thompson, 1980; 1986) to the tide gauge
# time series at Brest. 
# 
# Make sure that mulivariateAnalysisNewlyn.R has been sourced first.
#
# Author: simonh
###############################################################################

brestTime <- seq.Date(from=as.Date("1914/1/15"), to=as.Date("2006/12/15"), by="1 month")

# Model to be fitted
# p*g*n + pa = a*tx + b*ty + c11*cos(w1*t) + c12*sin(w1*t) 
# 				+ c21*cos(w2*t) + c22*sin(w2*t) + e
# where w1 = 2*pi/12 month and w2 = 2*pi/6 month, pa = local air pressure, 
# (p*g*n + pa) is the total pressure, (tx, ty) is influence of local wind stress
# modelled by (a*tx + b*ty) where (a,b) are regression coefficients to be
# determined from the data

totalPBrest <- 1025*9.8*brestMonthlyMean + as.vector(brestMetMonthlyMean[,1,])

data.brest <- data.frame(totP=totalPBrest, tauX=as.vector(brestMetMonthlyMean[,4,]), 
  tauY=as.vector(brestMetMonthlyMean[,5,]), t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2)
w1 <- 2*pi/12
w2 <- 2*pi/6

lm.brest <- lm(totP ~ tauX + tauY + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)), 
  data=data.brest)
resid.brest <- lm.brest$resid

x11()
par(family="HersheySans")
plot(time, data.brest$totP, type='l', col='orange')
lines(time, (2.737833e+02*cos(w1 * data.brest$t)+-4.358475e+01*sin(w1 * data.brest$t))-2.161031e+03, col='magenta', lwd=2)
lines(time, (-1.045618e+02*cos(w2 * data.brest$t)+4.918339e+01*sin(w2 * data.brest$t))-2.161031e+03, col='cyan', lwd=2)
lines(time, -4.992530e+06*data.brest$tauX - 2.161031e+03, col='red', lwd=2)
lines(time, 2.355104e+07*data.brest$tauY - 2.161031e+03, col='green', lwd=2)
#> range((2.737833e+02*cos(w1 * data.brest$t)+-4.358475e+01*sin(w1 * data.brest$t)))
#> 273.7833/(1025*9.8)
#[1] 0.02725568
# i.e. the mean seasonal cycle has an amplitude of 2.7cm cf 4cm at brest in Thompson 1986
# Deseasonalise the brest data
ds.totalPbrest <- totalPBrest-(lm.brest$coef[4]*cos(w1 * data.brest$t)+lm.brest$coef[5]*sin(w1 * data.brest$t))-(lm.brest$coef[6]*cos(w2 * data.brest$t)+lm.brest$coef[7]*sin(w2 * data.brest$t))+lm.brest$coef[1] 
# sd(totalPBrest/(1025*9.8), na.rm=T)
#[1] 0.0705519
# sd(ds.totalPbrest/(1025*9.8), na.rm=T)
# [1] 0.06383467 (i.e. 6.4cm)
# sd(resid.brest/(1025*9.8), na.rm=T)
#[1] 0.05502322 (i.e. 5.5cm)

# Calculate demeaned model data and convert to mm
dm.mm.brestMonthlyMean <- (brestMonthlyMean-mean(brestMonthlyMean))*1000

# Compare just over the period of the model
tg$mp.brest <- tg$brestMonthly[intersect(which(tg$time>=as.Date("1960/1/15")), which(tg$time<=as.Date("1999/12/15")))]
tg$mp.time <- tg$time[intersect(which(tg$time>=as.Date("1960/1/15")), which(tg$time<=as.Date("1999/12/15")))]

x11()
par(family="HersheySans")
plot(time, tg$mp.newlyn, type='l', col='orange')
lines(time,tg$mp.brest, col='blue')

x11()
par(family="HersheySans")
plot(tg$mp.brest, brestMonthlyMean)
cor.brest <- cor.test(tg$mp.brest, brestMonthlyMean, alternative="greater", na.action="na.fail")

# Demean tg data
tg$dm.mp.brest <- tg$mp.brest-mean(tg$mp.brest, na.rm=T)

# Detrend tg data
tg$data.brest <- data.frame(dm=tg$dm.mp.brest, t=seq(from=0,to=((40*12)-1)))
tg$lm.brest <- lm(dm ~ t, data=tg$data.brest)
tg$dt.mp.brest <- tg$dm.mp.brest - tg$lm.brest$coef[2]*tg$data.brest$t + tg$lm.brest$coef[1]

x11()
par(family="HersheySans")
plot(tg$dt.mp.brest, dm.mm.brestMonthlyMean)
cor.dt.brest <- cor.test(tg$dt.mp.brest, dm.mm.brestMonthlyMean, alternative="greater", na.action="na.fail")

tg$totalPNewlyn <- 1025*9.8*tg$mp.newlyn/1000 + as.vector(newlynMetMonthlyMean[,1,])
tg$totalPBrest <- 1025*9.8*tg$mp.brest/1000 + as.vector(brestMetMonthlyMean[,1,])

tg$data.newlyn <- data.frame(totP=tg$totalPNewlyn, tauX=as.vector(newlynMetMonthlyMean[,4,]), 
  tauY=as.vector(newlynMetMonthlyMean[,5,]), t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2)
tg$data.brest <- data.frame(totP=tg$totalPBrest, tauX=as.vector(brestMetMonthlyMean[,4,]), 
  tauY=as.vector(brestMetMonthlyMean[,5,]), t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2)
  
tg$lm.newlyn <- lm(totP ~ tauX + tauY + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)), 
  data=tg$data.newlyn)
tg$resid.newlyn <- tg$lm.newlyn$resid

tg$lm.brest <- lm(totP ~ tauX + tauY + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)), 
  data=tg$data.brest)
tg$resid.brest <- tg$lm.brest$resid

x11()
par(family="HersheySans")
plot(time, tg$data.newlyn$totP, type='l', col='orange')
lines(time, tg$data.brest$totP, col='blue')
lines(time, (tg$lm.brest$coef[4]*cos(w1 * tg$data.newlyn$t)+tg$lm.brest$coef[5]*sin(w1 * tg$data.newlyn$t))+tg$lm.brest$coef[1], col='magenta', lwd=2)
lines(time, (tg$lm.brest$coef[6]*cos(w2 * tg$data.newlyn$t)+tg$lm.brest$coef[7]*sin(w2 * tg$data.newlyn$t))+tg$lm.brest$coef[1], col='cyan', lwd=2)
lines(time, tg$lm.brest$coef[2]*tg$data.newlyn$tauX +tg$lm.brest$coef[1], col='red', lwd=2)
# Annual signal amplitude
# range(tg$lm.brest$coef[4]*cos(w1 * tg$data.newlyn$t)+tg$lm.brest$coef[5]*sin(w1 * tg$data.newlyn$t))
# 392.3397/(1025*9.8)
#[1] 0.03905821
# i.e. the mean seasonal cycle has an amplitude of 3.9cm cf 4cm at Newlyn in Thompson 1986

# Deseasonalise the Newlyn data
tg$ds.totalPNewlyn <- tg$totalPNewlyn-(tg$lm.newlyn$coef[4]*cos(w1 * tg$data.newlyn$t)+tg$lm.newlyn$coef[5]*sin(w1 * tg$data.newlyn$t))-(tg$lm.newlyn$coef[6]*cos(w2 * tg$data.newlyn$t)+tg$lm.newlyn$coef[7]*sin(w2 * tg$data.newlyn$t))+tg$lm.newlyn$coef[1] 
# sd(tg$totalPNewlyn/(1025*9.8), na.rm=T)
#[1] 0.0817145
# sd(tg$ds.totalPNewlyn/(1025*9.8), na.rm=T)
#[1] 0.06992447 (i.e. 7cm cf 7cm in Thompson, 1986)
# sd(tg$resid.newlyn/(1025*9.8), na.rm=T)
#[1] 0.05948289 (i.e. 5.9cm cf 3.1cm in Thompson, 1986)

# Check against interpolated data
interp<-new.env()
load("~/diskx/polcoms/brestNewlyn/met_data/brestNewlynMetMonthlyMeanInterp.RData", envir=interp)

interp$totalPNewlyn <- 1025*9.8*newlynMonthlyMean + as.vector(interp$newlynMetMonthlyMeanInterp[,1,])
interp$totalPBrest <- 1025*9.8*brestMonthlyMean + as.vector(interp$brestMetMonthlyMeanInterp[,1,])

interp$data.newlyn <- data.frame(totP=interp$totalPNewlyn, tauX=as.vector(interp$newlynMetMonthlyMeanInterp[,4,]), 
  tauY=as.vector(interp$newlynMetMonthlyMeanInterp[,5,]), t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2)
interp$data.brest <- data.frame(totP=interp$totalPBrest, tauX=as.vector(interp$brestMetMonthlyMeanInterp[,4,]), 
  tauY=as.vector(interp$brestMetMonthlyMeanInterp[,5,]), t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2)
  
interp$lm.newlyn <- lm(totP ~ t + tauX + tauY + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)), 
  data=interp$data.newlyn)
interp$resid.newlyn <- interp$lm.newlyn$resid

interp$lm.brest <- lm(totP ~ t + tauX + tauY + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)), 
  data=interp$data.brest)
interp$resid.brest <- interp$lm.brest$resid

x11()
par(family="HersheySans")
plot(time, interp$data.newlyn$totP, type='l', col='orange')
lines(time, interp$data.brest$totP, col='blue')
lines(time, (interp$lm.brest$coef[4]*cos(w1 * interp$data.brest$t)+interp$lm.brest$coef[5]*sin(w1 * interp$data.brest$t))+interp$lm.brest$coef[1], col='magenta', lwd=2)
lines(time, (interp$lm.brest$coef[6]*cos(w2 * interp$data.brest$t)+interp$lm.brest$coef[7]*sin(w2 * interp$data.brest$t))+interp$lm.brest$coef[1], col='cyan', lwd=2)
lines(time, interp$lm.brest$coef[3]*interp$data.brest$tauY +interp$lm.brest$coef[1], col='green', lwd=2)
lines(time, interp$lm.brest$coef[2]*interp$data.brest$tauX +interp$lm.brest$coef[1], col='red', lwd=2)
#> sd(interp$resid.newlyn/(1025*9.8), na.rm=T)
#[1] 0.05487602

#############
## Model 1 ##
#############
# Follow Thompson (1980) and see if we can improve
plymouth<-new.env()
load("~/diskx/polcoms/brestNewlyn/pressure/plymouthPressures.RData", envir=plymouth)

# Just use pressures over the periods of the model data (Pm) and TG data (Ptg)
plymouth$Pm <- plymouth$pressure[intersect(which(plymouth$time>=as.Date("1960/1/15")), which(plymouth$time<=as.Date("1999/12/15")))]
plymouth$Ptg <- plymouth$pressure[intersect(which(plymouth$time>=as.Date("1914/1/15")), which(plymouth$time<=as.Date("2006/12/15")))]

plymouth$data.newlyn <- data.frame(mmsl=newlynMonthlyMean, msl=tg$mp.newlyn,
  tauX=as.vector(interp$newlynMetMonthlyMeanInterp[,4,]), 
  tauY=as.vector(interp$newlynMetMonthlyMeanInterp[,5,]), t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2, 
  Pm=plymouth$Pm)

plymouth$m.lm.newlyn <- lm(mmsl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)) + Pm, 
  data=plymouth$data.newlyn)
plymouth$tg.lm.newlyn <- lm(msl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)) + Pm, 
  data=plymouth$data.newlyn)
  
plymouth$m.resid.newlyn <- plymouth$m.lm.newlyn$resid
plymouth$tg.resid.newlyn <- plymouth$tg.lm.newlyn$resid
#> sd(plymouth$data.newlyn$mmsl, na.rm=T)
#[1] 0.0676091
#> sd(plymouth$m.resid.newlyn, na.rm=T)
#[1] 0.03332927

#> sd(plymouth$data.newlyn$msl, na.rm=T)
#[1] 82.10498
#> sd(plymouth$tg.resid.newlyn, na.rm=T)
#[1] 41.12733

x11()
par(family="HersheySans")
plot(time, plymouth$data.newlyn$mmsl, type='l', col='orange')
lines(time, (plymouth$m.lm.newlyn$coef[3]*cos(w1 * plymouth$data.newlyn$t)+plymouth$m.lm.newlyn$coef[4]*sin(w1 * plymouth$data.newlyn$t))+plymouth$m.lm.newlyn$coef[1], col='magenta', lwd=2)
lines(time, (plymouth$m.lm.newlyn$coef[5]*cos(w2 * plymouth$data.newlyn$t)+plymouth$m.lm.newlyn$coef[6]*sin(w2 * plymouth$data.newlyn$t))+plymouth$m.lm.newlyn$coef[1], col='cyan', lwd=2)
lines(time, plymouth$m.lm.newlyn$coef[7]*plymouth$data.newlyn$Pm +plymouth$m.lm.newlyn$coef[1], col='green', lwd=2)
lines(time, plymouth$m.lm.newlyn$coef[2]*plymouth$data.newlyn$t +plymouth$m.lm.newlyn$coef[1], col='red', lwd=2)

plymouth$m.var.newlyn <- (plymouth$m.lm.newlyn$coef[3]*cos(w1 * plymouth$data.newlyn$t)+plymouth$m.lm.newlyn$coef[4]*sin(w1 * plymouth$data.newlyn$t)) +
                         (plymouth$m.lm.newlyn$coef[5]*cos(w2 * plymouth$data.newlyn$t)+plymouth$m.lm.newlyn$coef[6]*sin(w2 * plymouth$data.newlyn$t)) +
                         plymouth$m.lm.newlyn$coef[7]*plymouth$data.newlyn$Pm +
                         plymouth$m.lm.newlyn$coef[2]*plymouth$data.newlyn$t +
                         plymouth$m.lm.newlyn$coef[1]
plymouth$tg.var.newlyn <- (plymouth$tg.lm.newlyn$coef[3]*cos(w1 * plymouth$data.newlyn$t)+plymouth$tg.lm.newlyn$coef[4]*sin(w1 * plymouth$data.newlyn$t)) +
                         (plymouth$tg.lm.newlyn$coef[5]*cos(w2 * plymouth$data.newlyn$t)+plymouth$tg.lm.newlyn$coef[6]*sin(w2 * plymouth$data.newlyn$t)) +
                         plymouth$tg.lm.newlyn$coef[7]*plymouth$data.newlyn$Pm +
                         plymouth$tg.lm.newlyn$coef[2]*plymouth$data.newlyn$t +
                         plymouth$tg.lm.newlyn$coef[1]

# This is a big improvement. 

# Now look at all the Newlyn data. Also look at the Rossiter period (to 1962) and 
# the Cartwright period (to 1980)
# Plymouth pressures only run to 12/2004.
plymouth$PtgRossiter <- 
  plymouth$pressure[intersect(which(plymouth$time>=as.Date("1914/1/15")), 
  which(plymouth$time<=as.Date("1962/12/15")))]
plymouth$PtgCartwright <- 
  plymouth$pressure[intersect(which(plymouth$time>=as.Date("1914/1/15")), 
  which(plymouth$time<=as.Date("1980/12/15")))]
  
plymouth$tott=seq(from=0,to=((91*12)-1))
plymouth$Rot=seq(from=0,to=((49*12)-1))
plymouth$Cat=seq(from=0,to=((67*12)-1))

plymouth$data.newlyn.tot <- data.frame(msl=tg$newlynMonthly[1:1092], t=plymouth$tott, 
  Pm=plymouth$Ptg,  w1cos=cos(w1*plymouth$tott), w1sin=sin(w1*plymouth$tott), 
  w2cos=cos(w2*plymouth$tott), w2sin=sin(w2*plymouth$tott))
  
plymouth$data.newlyn.Ro <- data.frame(mslRo=tg$newlynMonthly[1:588], 
  PtgRo=plymouth$PtgRossiter, tRo=plymouth$Rot, 
  w1cosRo=cos(w1*plymouth$Rot), w1sinRo=sin(w1*plymouth$Rot), 
  w2cosRo=cos(w2*plymouth$Rot), w2sinRo=sin(w2*plymouth$Rot))
  
plymouth$data.newlyn.Ca <- data.frame(mslCa=tg$newlynMonthly[1:804],
  PtgCa=plymouth$PtgCartwright, tCa=plymouth$Cat,
  w1cosCa=cos(w1*plymouth$Cat), w1sinCa=sin(w1*plymouth$Cat), 
  w2cosCa=cos(w2*plymouth$Cat), w2sinCa=sin(w2*plymouth$Cat))
  
plymouth$tg.lm.newlyn.tot <- lm(msl ~ t + w1cos + w1sin + w2cos + w2sin + Pm,  
  data=plymouth$data.newlyn.tot)
plymouth$tg.resid.newlyn.tot <- plymouth$tg.lm.newlyn.tot$resid
#> sd(plymouth$data.newlyn.tot$msl, na.rm=T)
#[1] 91.69582
#> sd(plymouth$tg.resid.newlyn.tot, na.rm=T)
#[1] 40.46192
#> plymouth$tg.lm.newlyn.tot$coef[2]*12
#       t 
#1.822942 

plymouth$tg.lm.newlyn.Ro <- lm(mslRo ~ tRo + w1cosRo + w1sinRo + w2cosRo + w2sinRo + PtgRo,  
  data=plymouth$data.newlyn.Ro)
plymouth$tg.resid.newlyn.Ro <- plymouth$tg.lm.newlyn.Ro$resid
#> sd(plymouth$data.newlyn.Ro$msl, na.rm=T)
#[1] 83.08406
#> sd(plymouth$tg.resid.newlyn.Ro, na.rm=T)
#[1] 39.05382
#> plymouth$tg.lm.newlyn.Ro$coef[2]*12
#       t 
#2.273566
# cf 2.2 mm/yr in Rossiter (1972). This is effectively the same.

plymouth$tg.lm.newlyn.Ca <- lm(mslCa ~ tCa + w1cosCa + w1sinCa + w2cosCa + w2sinCa + PtgCa,  
  data=plymouth$data.newlyn.Ca)
plymouth$tg.resid.newlyn.Ca <- plymouth$tg.lm.newlyn.Ca$resid
#> sd(plymouth$data.newlyn.Ca$msl, na.rm=T)
#[1] 85.26331
#> sd(plymouth$tg.resid.newlyn.Ca, na.rm=T)
#[1] 40.21784
#> plymouth$tg.lm.newlyn.Ca$coef[2]*12
#       t 
#1.824541
# cf 1.34 mm/yr in Cartwright (1983). This is quite different!

#############
## Model 2 ##
#############  
# See if we can achieve the same thing with this model 2 of Thompson (1980) and
# static pressures from the model        
model2<-new.env()
# Just use pressures over the periods of the model data (Pm)
model2$Pm <- as.vector(interp$newlynMetMonthlyMeanInterp[,1,])

model2$data.newlyn <- data.frame(mmsl=newlynMonthlyMean, msl=tg$mp.newlyn,
  t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2, Pm=model2$Pm)

model2$m.lm.newlyn <- lm(mmsl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)) + Pm, 
  data=model2$data.newlyn)
model2$tg.lm.newlyn <- lm(msl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)) + Pm, 
  data=model2$data.newlyn)
  
model2$m.resid.newlyn <- model2$m.lm.newlyn$resid
model2$tg.resid.newlyn <- model2$tg.lm.newlyn$resid
#> sd(model2$data.newlyn$mmsl, na.rm=T)
#[1] 0.0676091
#> sd(model2$m.resid.newlyn, na.rm=T)
#[1] 0.0292884

#> sd(model2$data.newlyn$msl, na.rm=T)
#[1] 82.10498
#> sd(model2$tg.resid.newlyn, na.rm=T)
#[1] 35.91815

#############
## Model 3 ##
#############
# The next extension in Thompson is to use a linear combination of 9
# air pressures from around the UK. We first need to extract those air
# pressures 
model3 <- new.env() 
# Just use pressures over the periods of the model data (Pm)

model3$Pm1 <- as.vector(interp$extraMetMonthlyMeanInterp[,1,,1])
model3$Pm2 <- as.vector(interp$extraMetMonthlyMeanInterp[,1,,2])
model3$Pm3 <- as.vector(interp$extraMetMonthlyMeanInterp[,1,,3])
model3$Pm4 <- as.vector(interp$extraMetMonthlyMeanInterp[,1,,4])
model3$Pm5 <- as.vector(interp$extraMetMonthlyMeanInterp[,1,,5])
model3$Pm6 <- as.vector(interp$extraMetMonthlyMeanInterp[,1,,6])
model3$Pm7 <- as.vector(interp$extraMetMonthlyMeanInterp[,1,,7])
model3$Pm8 <- as.vector(interp$extraMetMonthlyMeanInterp[,1,,8])
model3$Pm9 <- as.vector(interp$extraMetMonthlyMeanInterp[,1,,9])

model3$data.newlyn <- data.frame(mmsl=newlynMonthlyMean, msl=tg$mp.newlyn,
  t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2, Pm1=model3$Pm1, Pm2=model3$Pm2,
  Pm3=model3$Pm3, Pm4=model3$Pm4, Pm5=model3$Pm5, Pm6=model3$Pm6, Pm7=model3$Pm7,
  Pm8=model3$Pm8, Pm9=model3$Pm9)

model3$m.lm.newlyn <- lm(mmsl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)) + 
  Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9,  data=model3$data.newlyn)
model3$tg.lm.newlyn <- lm(msl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)) + 
  Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9, data=model3$data.newlyn)
  
model3$m.resid.newlyn <- model3$m.lm.newlyn$resid
model3$tg.resid.newlyn <- model3$tg.lm.newlyn$resid
#> sd(model3$data.newlyn$mmsl, na.rm=T)
#[1] 0.0676091
#> sd(model3$m.resid.newlyn, na.rm=T)
#[1] 0.02449552

#> sd(model3$data.newlyn$msl, na.rm=T)
#[1] 82.10498
#> sd(model3$tg.resid.newlyn, na.rm=T)
#[1] 27.06582

# This is also an improvement!

#############
## Model 4 ##
#############
# Thompson's 4th model is rather harder to implement
# 
model4 <- new.env() 
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
model4$m.sd <- vector(mode="numeric", length=51)
model4$tg.sd <- vector(mode="numeric", length=51)

for(i in 1:51){
  model4$data.newlyn <- data.frame(mmsl=newlynMonthlyMean, msl=tg$mp.newlyn,
    t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2, Pm=model4$Pm[,i])

  model4$m.lm.newlyn <- lm(mmsl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) +
    I(sin(w2*t)) + Pm, data=model4$data.newlyn)
  model4$tg.lm.newlyn <- lm(msl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + 
    I(sin(w2*t)) + Pm, data=model4$data.newlyn)
  
  model4$m.resid.newlyn <- model4$m.lm.newlyn$resid
  model4$tg.resid.newlyn <- model4$tg.lm.newlyn$resid
  
  model4$m.sd[i] <- sd(model4$m.resid.newlyn, na.rm=T)
  model4$tg.sd[i] <- sd(model4$tg.resid.newlyn, na.rm=T)
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model4$m.sd.sorted <- sort(model4$m.sd, index=TRUE)
model4$tg.sd.sorted <- sort(model4$tg.sd, index=TRUE)

# Now construct models of the top 9 components
model4$data.newlyn <- data.frame(mmsl=newlynMonthlyMean, msl=tg$mp.newlyn,
  t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2, Pm=model2$Pm,
  Pm1=model4$Pm[,model4$m.sd.sorted$ix[1]], 
  Pm2=model4$Pm[,model4$m.sd.sorted$ix[2]], Pm3=model4$Pm[,model4$m.sd.sorted$ix[3]],
  Pm4=model4$Pm[,model4$m.sd.sorted$ix[4]], Pm5=model4$Pm[,model4$m.sd.sorted$ix[5]],
  Pm6=model4$Pm[,model4$m.sd.sorted$ix[6]], Pm7=model4$Pm[,model4$m.sd.sorted$ix[7]],
  Pm8=model4$Pm[,model4$m.sd.sorted$ix[8]], Pm9=model4$Pm[,model4$m.sd.sorted$ix[9]])

model4$m.lm.newlyn <- lm(mmsl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)) + 
  Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 + Pm,  data=model4$data.newlyn)

model4$data.newlyn <- data.frame(mmsl=newlynMonthlyMean, msl=tg$mp.newlyn,
  t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2, Pm=model2$Pm,
  Pm1=model4$Pm[,model4$tg.sd.sorted$ix[1]], 
  Pm2=model4$Pm[,model4$tg.sd.sorted$ix[2]], Pm3=model4$Pm[,model4$tg.sd.sorted$ix[3]],
  Pm4=model4$Pm[,model4$tg.sd.sorted$ix[4]], Pm5=model4$Pm[,model4$tg.sd.sorted$ix[5]],
  Pm6=model4$Pm[,model4$tg.sd.sorted$ix[6]], Pm7=model4$Pm[,model4$tg.sd.sorted$ix[7]],
  Pm8=model4$Pm[,model4$tg.sd.sorted$ix[8]], Pm9=model4$Pm[,model4$tg.sd.sorted$ix[9]])

model4$tg.lm.newlyn <- lm(msl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)) + 
  Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 + Pm,  data=model4$data.newlyn)

#model4$m.resid.newlyn <- model4$m.lm.newlyn$resid
#model4$tg.resid.newlyn <- model4$tg.lm.newlyn$resid

model4$m.resid.newlyn <- model4$data.newlyn$mmsl - model4$m.lm.newlyn$coef[2]*model4$data.newlyn$t + model4$m.lm.newlyn$coef[1]
model4$tg.resid.newlyn <- model4$data.newlyn$msl - model4$tg.lm.newlyn$coef[2]*model4$data.newlyn$t + model4$tg.lm.newlyn$coef[1]

#> sd(model4$data.newlyn$mmsl, na.rm=T)
#[1] 0.0676091
#> sd(model4$m.resid.newlyn, na.rm=T)
#[1] 0.02422649

#> sd(model4$data.newlyn$msl, na.rm=T)
#[1] 82.10498
#> sd(model4$tg.resid.newlyn, na.rm=T)
#[1] 25.16201
# No signficant improvement

#############
## Model 5 ##
#############
# Now the new stuff. Add deep sea boundary conditions.
# Treat the 6 deep sea points in the same way as we did the pressure
# points with 3 lags (=18 combinations)
model5 <- new.env()
model5$Ds <- array(NA,dim=c(480,18))
# Zero lag
for(i in 1:6){
  model5$Ds[,i] <- deepSeaMonthlyMean[,i]
}
# Lag 1
for(i in 7:12){
  model5$Ds[2:480,i] <- deepSeaMonthlyMean[1:479,(i-6)]
}
# Lag 2
for(i in 13:18){
  model5$Ds[3:480,i] <- deepSeaMonthlyMean[1:478,(i-12)]
}

# What we are going to do now is calculate the sd of the residuals as each of the
# 18 components is added in turn 
model5$m.sd <- vector(mode="numeric", length=18)
model5$tg.sd <- vector(mode="numeric", length=18)

for(i in 1:18){
  model5$data.newlyn <- data.frame(mmsl=newlynMonthlyMean, msl=tg$mp.newlyn,
    t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2, Pm1=model4$Pm[,model4$m.sd.sorted$ix[1]], 
    Pm2=model4$Pm[,model4$m.sd.sorted$ix[2]], Pm3=model4$Pm[,model4$m.sd.sorted$ix[3]],
    Pm4=model4$Pm[,model4$m.sd.sorted$ix[4]], Pm5=model4$Pm[,model4$m.sd.sorted$ix[5]],
    Pm6=model4$Pm[,model4$m.sd.sorted$ix[6]], Pm7=model4$Pm[,model4$m.sd.sorted$ix[7]],
    Pm8=model4$Pm[,model4$m.sd.sorted$ix[8]], Pm9=model4$Pm[,model4$m.sd.sorted$ix[9]],
    dssl=model5$Ds[,i])

  model5$m.lm.newlyn <- lm(mmsl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)) + 
    Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 + dssl,  data=model5$data.newlyn)
    
  model5$tg.lm.newlyn <- lm(msl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)) + 
    Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 + dssl,  data=model5$data.newlyn)
  
  model5$m.resid.newlyn <- model5$m.lm.newlyn$resid
  model5$tg.resid.newlyn <- model5$tg.lm.newlyn$resid
  
  model5$m.sd[i] <- sd(model5$m.resid.newlyn, na.rm=T)
  model5$tg.sd[i] <- sd(model5$tg.resid.newlyn, na.rm=T)
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model5$m.sd.sorted <- sort(model5$m.sd, index=TRUE)
model5$tg.sd.sorted <- sort(model5$tg.sd, index=TRUE)

model5$t <- seq(from=0,to=((40*12)-1))

# Now construct models of the top 9 components
model5$data.newlyn <- data.frame(mmsl=newlynMonthlyMean, msl=tg$mp.newlyn,
  t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2, Pm1=model4$Pm[,model4$m.sd.sorted$ix[1]], 
  Pm2=model4$Pm[,model4$m.sd.sorted$ix[2]], Pm3=model4$Pm[,model4$m.sd.sorted$ix[3]],
  Pm4=model4$Pm[,model4$m.sd.sorted$ix[4]], Pm5=model4$Pm[,model4$m.sd.sorted$ix[5]],
  Pm6=model4$Pm[,model4$m.sd.sorted$ix[6]], Pm7=model4$Pm[,model4$m.sd.sorted$ix[7]],
  Pm8=model4$Pm[,model4$m.sd.sorted$ix[8]], Pm9=model4$Pm[,model4$m.sd.sorted$ix[9]],
  Ds1=model5$Ds[,model5$m.sd.sorted$ix[1]], 
  Ds2=model5$Ds[,model5$m.sd.sorted$ix[2]], Ds3=model5$Ds[,model5$m.sd.sorted$ix[3]],
  Ds4=model5$Ds[,model5$m.sd.sorted$ix[4]], Ds5=model5$Ds[,model5$m.sd.sorted$ix[5]],
  Ds6=model5$Ds[,model5$m.sd.sorted$ix[6]], Ds7=model5$Ds[,model5$m.sd.sorted$ix[7]],
  Ds8=model5$Ds[,model5$m.sd.sorted$ix[8]], Ds9=model5$Ds[,model5$m.sd.sorted$ix[9]])

model5$m.lm.newlyn <- lm(mmsl ~ t + w1cos + w1sin + w2cos + w2sin + 
  Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 +
  Ds1 + Ds2 + Ds3 + Ds4 + Ds5 + Ds6 + Ds7 + Ds8 + Ds9,  data=model5$data.newlyn)

model5$data.newlyn <- data.frame(mmsl=newlynMonthlyMean, msl=tg$mp.newlyn,
  t=model5$t, w1cos=cos(w1*model5$t), w1sin=sin(w1*model5$t), w2cos=cos(w2*model5$t), w2sin=sin(w2*model5$t), 
  Pm1=model4$Pm[,model4$tg.sd.sorted$ix[1]], 
  Pm2=model4$Pm[,model4$tg.sd.sorted$ix[2]], Pm3=model4$Pm[,model4$tg.sd.sorted$ix[3]],
  Pm4=model4$Pm[,model4$tg.sd.sorted$ix[4]], Pm5=model4$Pm[,model4$tg.sd.sorted$ix[5]],
  Pm6=model4$Pm[,model4$tg.sd.sorted$ix[6]], Pm7=model4$Pm[,model4$tg.sd.sorted$ix[7]],
  Pm8=model4$Pm[,model4$tg.sd.sorted$ix[8]], Pm9=model4$Pm[,model4$tg.sd.sorted$ix[9]],
  Ds1=model5$Ds[,model5$tg.sd.sorted$ix[1]], 
  Ds2=model5$Ds[,model5$tg.sd.sorted$ix[2]], Ds3=model5$Ds[,model5$tg.sd.sorted$ix[3]],
  Ds4=model5$Ds[,model5$tg.sd.sorted$ix[4]], Ds5=model5$Ds[,model5$tg.sd.sorted$ix[5]],
  Ds6=model5$Ds[,model5$tg.sd.sorted$ix[6]], Ds7=model5$Ds[,model5$tg.sd.sorted$ix[7]],
  Ds8=model5$Ds[,model5$tg.sd.sorted$ix[8]], Ds9=model5$Ds[,model5$tg.sd.sorted$ix[9]])

model5$tg.lm.newlyn <- lm(msl ~ t + w1cos + w1sin + w2cos + w2sin + 
  Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 +
  Ds1 + Ds2 + Ds3 + Ds4 + Ds5 + Ds6 + Ds7 + Ds8 + Ds9,  data=model5$data.newlyn)
  
#model5$m.resid.newlyn <- model5$m.lm.newlyn$resid
#model5$tg.resid.newlyn <- model5$tg.lm.newlyn$resid

model5$m.fitted.newlyn <- model5$m.lm.newlyn$coef[1]
for (i in 2:24){
	model5$m.fitted.newlyn <- model5$m.fitted.newlyn + model5$m.lm.newlyn$coef[i]*model5$data.newlyn[[i+1]]
}
model5$m.resid.newlyn <- model5$data.newlyn$mmsl - model5$m.fitted.newlyn

x11()
plot(tg$mp.time,dm.mm.newlynMonthlyMean, type='l', col='blue')
lines(tg$mp.time, (model5$m.fitted.newlyn-mean(model5$m.fitted.newlyn, na.rm=T))*1000, col='red')
lines(tg$mp.time, model5$m.resid.newlyn*1000, col='magenta', lwd=2)


model5$tg.fitted.newlyn <- model5$tg.lm.newlyn$coef[1]
for (i in 2:24){
	model5$tg.fitted.newlyn <- model5$tg.fitted.newlyn + model5$tg.lm.newlyn$coef[i]*model5$data.newlyn[[i+1]]
}
model5$tg.resid.newlyn <- model5$data.newlyn$msl - model5$tg.fitted.newlyn

x11()
plot(tg$mp.time,tg$dm.mp.newlyn, type='l', col='blue')
lines(tg$mp.time, model5$tg.fitted.newlyn-mean(model5$tg.fitted.newlyn, na.rm=T), col='red')
lines(tg$mp.time, model5$tg.resid.newlyn, col='magenta', lwd=2)

#> sd(model5$data.newlyn$mmsl, na.rm=T)
#[1] 0.0676091
#> sd(model5$m.resid.newlyn, na.rm=T)
#[1] 0.01615011

#> sd(model5$data.newlyn$msl, na.rm=T)
#[1] 82.10498
#> sd(model5$tg.resid.newlyn, na.rm=T)
#[1] 23.78966

# Generate annual means from regression data
model5$numYrs <- 40

model5$annualYrs <- seq.Date(from=as.Date("1960/7/1"), to=as.Date("1999/7/1"), by="1 year")

model5$tmp <- tg$dm.mp.newlyn
dim(model5$tmp) <- c(12,model5$numYrs)
model5$tg.annualMeans.newlyn <- colMeans(model5$tmp)

model5$tmp <- model5$tg.fitted.newlyn-mean(model5$tg.fitted.newlyn, na.rm=T)
dim(model5$tmp) <- c(12,model5$numYrs)
model5$tg.fitted.annualMeans.newlyn <- colMeans(model5$tmp)

model5$tmp <- model5$tg.resid.newlyn
dim(model5$tmp) <- c(12,model5$numYrs)
model5$tg.resid.annualMeans.newlyn <- colMeans(model5$tmp)

rm(tmp, envir=model5)

x11()
plot(model5$annualYrs, model5$tg.annualMeans.newlyn, type='l', col='blue')
lines(model5$annualYrs, model5$tg.fitted.annualMeans.newlyn, col='red')
lines(model5$annualYrs, model5$tg.resid.annualMeans.newlyn, col='magenta', lwd=2)
#var(model5$tg.annualMeans.newlyn, na.rm=T)
#[1] 807.7804
#var(model5$tg.resid.annualMeans.newlyn, na.rm=T)
#[1] 176.7980
#var(model5$tg.resid.annualMeans.newlyn, na.rm=T)/var(model5$tg.annualMeans.newlyn, na.rm=T)*100
#[1] 21.88689

# Produce decadal means to examine reduction in decadal variance
model5$tg.tenYrTrends.newlyn <- vector(mode='numeric',length=(model5$numYrs-9))
model5$tg.resid.tenYrTrends.newlyn <- model5$tg.tenYrTrends.newlyn 
model5$tg.fitted.tenYrTrends.newlyn <- model5$tg.tenYrTrends.newlyn

model5$midpointYrs <- seq.Date(from=as.Date("1965/1/1"), length=(model5$numYrs-9), by="1 year")
for (i in 1:(model5$numYrs-9)) {
	if (length(which(is.finite(model5$tg.annualMeans.newlyn[i:(i+9)])))>=7){
		model5$junk <- lm(model5$tg.annualMeans.newlyn[i:(i+9)] ~ c(1:10), na.action = na.exclude)
		model5$tg.tenYrTrends.newlyn[i] <- model5$junk$coef[2]
	} else {
		model5$tg.tenYrTrends.newlyn[i] <- NA
	}
	if (length(which(is.finite(model5$tg.fitted.annualMeans.newlyn[i:(i+9)])))>=7){
		model5$junk <- lm(model5$tg.fitted.annualMeans.newlyn[i:(i+9)] ~ c(1:10), na.action = na.exclude)
		model5$tg.fitted.tenYrTrends.newlyn[i] <- model5$junk$coef[2]
	} else {
		model5$tg.fitted.tenYrTrends.newlyn[i] <- NA
	}
	if (length(which(is.finite(model5$tg.resid.annualMeans.newlyn[i:(i+9)])))>=7){
		model5$junk <- lm(model5$tg.resid.annualMeans.newlyn[i:(i+9)] ~ c(1:10), na.action = na.exclude)
		model5$tg.resid.tenYrTrends.newlyn[i] <- model5$junk$coef[2]
	} else {
		model5$tg.resid.tenYrTrends.newlyn[i] <- NA
	}
	
}
rm(junk, envir=model5)

x11()
plot(model5$midpointYrs, model5$tg.tenYrTrends.newlyn, type='l', col='blue')
lines(model5$midpointYrs, model5$tg.fitted.tenYrTrends.newlyn, col='red')
lines(model5$midpointYrs, model5$tg.resid.tenYrTrends.newlyn, col='magenta', lwd=2)
#var(model5$tg.tenYrTrends.newlyn, na.rm=T)
#[1] 15.47344
#var(model5$tg.resid.tenYrTrends.newlyn, na.rm=T)
#[1] 5.570029
#var(model5$tg.resid.tenYrTrends.newlyn, na.rm=T)/var(model5$tg.tenYrTrends.newlyn, na.rm=T)*100
#[1] 35.99736

#############
## Model 6 ##
#############
# Use winds from the model instead of pressures
# Treat the 18 (Newlyn + 17 others) wind points in the same way as we did the pressure
# points with 3 lags for the 17 distant stations (=51 combinations) but
# just use zero lag for Newlyn
# Es is the E wind stress, Ns is the N wind stress
# 
model6 <- new.env() 
# Again, just use winds over the periods of the model data (Pm)

# As with the pressures we use only the first 9 components which we choose by stepwise
# regression. Again, we achieve this by adding in each component in turn and calculating
# the variance reduction, which we can then rank and analyze.
model6$Ns <- array(NA,dim=c(480,51))
model6$Es <- array(NA,dim=c(480,51))

model6$NewlynEs <- as.vector(interp$newlynMetMonthlyMeanInterp[,2,])
model6$NewlynNs <- as.vector(interp$newlynMetMonthlyMeanInterp[,3,])

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
  model6$data.newlyn <- data.frame(mmsl=newlynMonthlyMean, msl=tg$mp.newlyn,
    t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2, Pm=model2$Pm,
    NEs=model6$NewlynEs, NNs=model6$NewlynNs, Es=model6$Es[,i], Ns=model6$Ns[,i])

  model6$m.lm.newlyn <- lm(mmsl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) +
    I(sin(w2*t)) + Es + Ns + NEs + NNs + Pm, data=model6$data.newlyn)
  model6$tg.lm.newlyn <- lm(msl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + 
    I(sin(w2*t)) + Es + Ns + NEs + NNs + Pm, data=model6$data.newlyn)
  
  model6$m.resid.newlyn <- model6$m.lm.newlyn$resid
  model6$tg.resid.newlyn <- model6$tg.lm.newlyn$resid
  
  model6$m.sd[i] <- sd(model6$m.resid.newlyn, na.rm=T)
  model6$tg.sd[i] <- sd(model6$tg.resid.newlyn, na.rm=T)
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model6$m.sd.sorted <- sort(model6$m.sd, index=TRUE)
model6$tg.sd.sorted <- sort(model6$tg.sd, index=TRUE)

# Now construct models of the top 9 components
model6$t <- seq(from=0,to=((40*12)-1))
model6$data.newlyn <- data.frame(mmsl=newlynMonthlyMean, msl=tg$mp.newlyn,
  t=model6$t, w1=w1, w2=w2, 
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
  NEs=model6$NewlynEs, NNs=model6$NewlynNs, Pm=model2$Pm,
  w1cos=cos(w1*model6$t), w1sin=sin(w1*model6$t), 
  w2cos=cos(w2*model6$t), w2sin=sin(w2*model6$t))

model6$m.lm.newlyn <- lm(mmsl ~ t + w1cos + w1sin + w2cos + w2sin + 
  Es1 + Es2 + Es3 + Es4 + Es5 + Es6 + Es7 + Es8 + Es9 +
  Ns1 + Ns2 + Ns3 + Ns4 + Ns5 + Ns6 + Ns7 + Ns8 + Ns9 +
  NEs + NNs + Pm,  data=model6$data.newlyn)

model6$data.newlyn <- data.frame(mmsl=newlynMonthlyMean, msl=tg$mp.newlyn,
  t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2, 
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
  NEs=model6$NewlynEs, NNs=model6$NewlynNs, Pm=model2$Pm,
  w1cos=cos(w1*model6$t), w1sin=sin(w1*model6$t), 
  w2cos=cos(w2*model6$t), w2sin=sin(w2*model6$t))

model6$tg.lm.newlyn <- lm(msl ~ t + w1cos + w1sin + w2cos + w2sin + 
  Es1 + Es2 + Es3 + Es4 + Es5 + Es6 + Es7 + Es8 + Es9 +
  Ns1 + Ns2 + Ns3 + Ns4 + Ns5 + Ns6 + Ns7 + Ns8 + Ns9 +
  NEs + NNs + Pm,  data=model6$data.newlyn)

model6$m.resid.newlyn <- model6$m.lm.newlyn$resid
model6$tg.resid.newlyn <- model6$tg.lm.newlyn$resid

#> sd(model6$data.newlyn$mmsl, na.rm=T)
#[1] 0.0676091
#> sd(model6$m.resid.newlyn, na.rm=T)
#[1] 0.02276734

#> sd(model6$data.newlyn$msl, na.rm=T)
#[1] 82.10498
#> sd(model6$tg.resid.newlyn, na.rm=T)
#[1] 22.99741

# So wind stresses (or just winds) actually do worse than just static pressures. Why?
# Dynamic height in the model is not the same as sea level so I should diagnose and use a
# linear combination of diagnostic heights first.
# Jason also suggests using the currents. Although this isn't a very useful tool
# for reconstructing sea level (unless we have a model) it is useful as a diagnostic tool
# to tell us what is actually controlling SL at the coast.
