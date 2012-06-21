# multivariateAnalysisBrest3.R
# 
# Implementation of Thompson's methods (Thompson, 1980; 1986) to the tide gauge
# time series at Brest
# These are the same methods implemented in multivariateAnalysisNewlyn.R
#
# Version 3 is changed so that we use total pressure throughout, as in
# Thompson '86, rather than having it as a forcing term on the RHS (as
# in Thompson '80). This avoids the problem of aliasing.
# We also include HadSLP2r data for the first time.
#
# Author: simonh
###############################################################################

#load("~/diskx/polcoms/brestNewlyn/newlynBrestTG.RData")
load("~/diskx/polcoms/brestNewlyn/brestNewlynModelMonthlyMean.RData")
load("~/diskx/polcoms/brestNewlyn/met_data/brestNewlynMetMonthlyMeanInterp.RData")
load("~/diskx/polcoms/brestNewlyn/met_data/brestNewlynMetMonthlyMean.RData")
#rm(junk, drv,con,i, selectString,res)

# Model to be fitted
# p*g*n + pa = a*tx + b*ty + c11*cos(w1*t) + c12*sin(w1*t) 
# 				+ c21*cos(w2*t) + c22*sin(w2*t) + e
# where w1 = 2*pi/12 month and w2 = 2*pi/6 month, pa = local air pressure, 
# (p*g*n + pa) is the total pressure, (tx, ty) is influence of local wind stress
# modelled by (a*tx + b*ty) where (a,b) are regression coefficients to be
# determined from the data

# Determined from model data
w1 <- 2*pi/12
w2 <- 2*pi/6
time <- seq.Date(from=as.Date("1960/1/15"), to=as.Date("1999/12/15"), by="1 month")
totalPBrest <- 1025*9.8*brestMonthlyMean + as.vector(brestMetMonthlyMean[,1,])
data.brest <- data.frame(totP=totalPBrest, tauX=as.vector(brestMetMonthlyMean[,4,]), 
  tauY=as.vector(brestMetMonthlyMean[,5,]), t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2)

lm.brest <- lm(totP ~ tauX + tauY + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)), 
  data=data.brest)
resid.brest <- lm.brest$resid

bc <- coef(lm.brest)

x11()
par(family="HersheySans")
plot(time, data.brest$totP, type='l', col='orange')
lines(time, (bc[6]*cos(w2 * data.brest$t)+bc[7]*sin(w2 * data.brest$t))+bc[1], col='cyan', lwd=2)
lines(time, bc[2]*data.brest$tauX + bc[1], col='red', lwd=2)
lines(time, bc[3]*data.brest$tauY + bc[1], col='green', lwd=2)
#> range((bc[6]*cos(w1 * data.brest$t)+bc[7]*sin(w1 * data.brest$t)))
#> 171.7752/(1025*9.8)
#[1] 0.01710057
# i.e. the mean seasonal cycle has an amplitude of 1.7cm cf 2.7 cm at Newlyn in analysis and 4cm at Newlyn in Thompson 1986
# Deseasonalise the Brest data
ds.totalPBrest <- totalPBrest-(lm.brest$coef[4]*cos(w1 * data.brest$t)+lm.brest$coef[5]*sin(w1 * data.brest$t))-(lm.brest$coef[6]*cos(w2 * data.brest$t)+lm.brest$coef[7]*sin(w2 * data.brest$t))+lm.brest$coef[1] 
# sd(totalPBrest/(1025*9.8), na.rm=T)
#[1] 0.0705519 cf 0.06716031 a Newlyn
# sd(ds.totalPBrest/(1025*9.8), na.rm=T)
# [1] 0.06383467 cf 0.06170118 at Newlyn (i.e. 6.2cm cf 7cm in Thompson, 1986)
# sd(resid.brest/(1025*9.8), na.rm=T)
#[1] 0.05502322 cf 0.05399438 at Newlyn (i.e. 5.4cm cf 3.1cm in Thompson, 1986)

# Calculate demeaned model data and convert to mm
dm.mm.brestMonthlyMean <- (brestMonthlyMean-mean(brestMonthlyMean))*1000
dm.mm.newlynMonthlyMean <- (newlynMonthlyMean-mean(newlynMonthlyMean))*1000

x11()
plot(time, dm.mm.brestMonthlyMean, type='l', col='blue', ylim=c(-200,300))
lines(time, dm.mm.newlynMonthlyMean, col='red')
lines(time, dm.mm.brestMonthlyMean-dm.mm.newlynMonthlyMean, col='magenta', lwd=2)


# Add model data from S12run408 with realistic deep ocean boundary condition 
# which varies over the 45 year of the run (1960-2004) unlike the previous model
# S12run405
S12R408<-new.env()
load("~/diskx/polcoms/S12run408/brestNewlyn/brestNewlynMonthlyS12R408.RData", envir=S12R408)
S12R408$time <- S12R408$monthsArray
x11()
plot(S12R408$time, S12R408$brestMonthlyMean, type='l', col='orange')
# Compare with TG data
tg<-new.env()
load("~/diskx/polcoms/brestNewlyn/newlynBrestTG.RData", envir=tg)
tg$time <- seq.Date(from=as.Date("1914/1/15"), to=as.Date("2006/12/15"), by="1 month")
# Compare just over the period of the model
tg$mp.brest <- tg$brestMonthly[intersect(which(tg$time>=as.Date("1960/1/15")), which(tg$time<=as.Date("1999/12/15")))]
tg$mp.time <- tg$time[intersect(which(tg$time>=as.Date("1960/1/15")), which(tg$time<=as.Date("1999/12/15")))]
tg$mp.brest.orig <- tg$mp.brest
tg$mp.brest[which(is.na(tg$mp.brest))] <- mean(tg$mp.brest, na.rm=T)

tg$mp.newlyn <- tg$newlynMonthly[intersect(which(tg$time>=as.Date("1960/1/15")), which(tg$time<=as.Date("1999/12/15")))]
tg$mp.newlyn.orig <- tg$mp.newlyn
tg$mp.newlyn[which(is.na(tg$mp.newlyn))] <- mean(tg$mp.newlyn, na.rm=T)

x11()
par(family="HersheySans")
plot(time, tg$mp.brest, type='l', col='orange')

x11()
par(family="HersheySans")
plot(tg$mp.brest, brestMonthlyMean)
cor.brest <- cor.test(tg$mp.brest, brestMonthlyMean, alternative="greater", na.action="na.fail")

# Demean tg data
tg$dm.mp.brest <- tg$mp.brest-mean(tg$mp.brest, na.rm=T)
tg$dm.mp.newlyn <- tg$mp.newlyn-mean(tg$mp.newlyn, na.rm=T)

resid <- tg$dm.mp.brest-dm.mm.brestMonthlyMean

x11()
plot(tg$mp.time, tg$dm.mp.brest, type='l', col='blue', ylim=c(-200,300))
lines(tg$mp.time, tg$dm.mp.newlyn, col='red')
lines(tg$mp.time, tg$dm.mp.brest-tg$dm.mp.newlyn, col='magenta', lwd=2)

x11()
spectrum(tg$dm.mp.newlyn, kernel=kernel("daniell", c(11,7,3)), col='red')
x11()
spectrum(tg$dm.mp.brest, kernel=kernel("daniell", c(11,7,3)), col='blue')
x11()
spectrum((tg$dm.mp.brest-tg$dm.mp.newlyn), kernel=kernel("daniell", c(11,7,3)), col='blue')


#> sd(resid, na.rm=T)/sd(tg$dm.mp.brest, na.rm=T)*100
#[1] 62.03784 Newlyn: 59.66754
#> sd(resid, na.rm=T)
#[1] 53.96208 Newlyn: 48.99002
#> sd(tg$dm.mp.brest, na.rm=T)
#[1] 86.98253 Newlyn: 82.10498
#> sd(dm.mm.brestMonthlyMean, na.rm=T)
#[1] 70.98115 Newlyn: 67.6091

# Detrend tg data
tg$data.brest <- data.frame(dm=tg$dm.mp.brest, t=seq(from=0,to=((40*12)-1)))
tg$lm.brest.time <- lm(dm ~ t, data=tg$data.brest)
#> tg$lm.brest.time$coef[2]*12
#t 
#1.363057 mm/yr over 1960-2000

tg$dt.mp.brest <- tg$dm.mp.brest - tg$lm.brest$coef[2]*tg$data.brest$t + tg$lm.brest$coef[1]

x11()
par(family="HersheySans")
plot(tg$dt.mp.brest, dm.mm.brestMonthlyMean)
cor.dt.brest <- cor.test(tg$dt.mp.brest, dm.mm.brestMonthlyMean, alternative="greater", na.action="na.fail")

tg$totalPBrest <- 1025*9.8*tg$mp.brest/1000 + as.vector(brestMetMonthlyMean[,1,])

tg$data.brest <- data.frame(totP=tg$totalPBrest, tauX=as.vector(brestMetMonthlyMean[,4,]), 
  tauY=as.vector(brestMetMonthlyMean[,5,]), t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2)

tg$lm.brest <- lm(totP ~ tauX + tauY + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)), 
  data=tg$data.brest)
tg$resid.brest <- tg$lm.brest$resid

x11()
par(family="HersheySans")
plot(time, tg$data.brest$totP, type='l', col='orange')
lines(time, (tg$lm.brest$coef[4]*cos(w1 * tg$data.brest$t)+tg$lm.brest$coef[5]*sin(w1 * tg$data.brest$t))+tg$lm.brest$coef[1], col='magenta', lwd=2)
lines(time, (tg$lm.brest$coef[6]*cos(w2 * tg$data.brest$t)+tg$lm.brest$coef[7]*sin(w2 * tg$data.brest$t))+tg$lm.brest$coef[1], col='cyan', lwd=2)
lines(time, tg$lm.brest$coef[2]*tg$data.brest$tauX +tg$lm.brest$coef[1], col='red', lwd=2)
# Annual signal amplitude
# range(tg$lm.brest$coef[4]*cos(w1 * tg$data.brest$t)+tg$lm.brest$coef[5]*sin(w1 * tg$data.brest$t))
# 392.3397/(1025*9.8)
#[1] 0.03905821
# i.e. the mean seasonal cycle has an amplitude of 3.9cm cf 5cm at Newlyn and 4cm at Newlyn in Thompson 1986

# Semi-annual signal
#> range(tg$lm.brest$coef[6]*cos(w1 * tg$data.brest$t)+tg$lm.brest$coef[7]*sin(w1 * tg$data.brest$t))
#[1] -237.3465  237.3465
#> 237.3465/(1025*9.8)
#[1] 0.02362832 Newlyn: 0.01894187

# Deseasonalise the Brest data
tg$ds.totalPBrest <- tg$totalPBrest-(tg$lm.brest$coef[4]*cos(w1 * tg$data.brest$t)+tg$lm.brest$coef[5]*sin(w1 * tg$data.brest$t))-(tg$lm.brest$coef[6]*cos(w2 * tg$data.brest$t)+tg$lm.brest$coef[7]*sin(w2 * tg$data.brest$t))+tg$lm.brest$coef[1] 
# sd(tg$totalPBrest/(1025*9.8), na.rm=T)
#[1] 0.08678313 cf 0.0817145 at Newlyn
# sd(tg$ds.totalPBrest/(1025*9.8), na.rm=T)
#[1] 0.07774023 cf 0.06992447 at Newlyn (i.e. 7cm cf 7cm in Thompson, 1986)
# sd(tg$resid.brest/(1025*9.8), na.rm=T)
#[1] 0.06585905 cf 0.05948289 at Newlyn (i.e. 5.9cm cf 3.1cm in Thompson, 1986)

# Check against interpolated data
interp<-new.env()
load("~/diskx/polcoms/brestNewlyn/met_data/brestNewlynMetMonthlyMeanInterp.RData", envir=interp)

interp$totalPBrest <- 1025*9.8*brestMonthlyMean + as.vector(interp$brestMetMonthlyMeanInterp[,1,])

interp$data.brest <- data.frame(totP=interp$totalPBrest, tauX=as.vector(interp$brestMetMonthlyMeanInterp[,4,]), 
  tauY=as.vector(interp$brestMetMonthlyMeanInterp[,5,]), t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2)
  
interp$lm.brest <- lm(totP ~ t + tauX + tauY + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)), 
  data=interp$data.brest)
interp$resid.brest <- interp$lm.brest$resid

x11()
par(family="HersheySans")
plot(time, interp$data.brest$totP, type='l', col='orange')
lines(time, (interp$lm.brest$coef[5]*cos(w1 * interp$data.brest$t)+interp$lm.brest$coef[6]*sin(w1 * interp$data.brest$t))+interp$lm.brest$coef[1], col='magenta', lwd=2)
lines(time, (interp$lm.brest$coef[7]*cos(w2 * interp$data.brest$t)+interp$lm.brest$coef[8]*sin(w2 * interp$data.brest$t))+interp$lm.brest$coef[1], col='cyan', lwd=2)
lines(time, interp$lm.brest$coef[4]*interp$data.brest$tauY +interp$lm.brest$coef[1], col='green', lwd=2)
lines(time, interp$lm.brest$coef[3]*interp$data.brest$tauX +interp$lm.brest$coef[1], col='red', lwd=2)
#> sd(interp$resid.brest/(1025*9.8), na.rm=T)
#[1] 0.05462631 cf 0.05487602 at Newlyn

#############
## Model 1 ##
#############
# Follow Thompson (1980) and see if we can improve
plymouth<-new.env()
load("~/diskx/polcoms/brestNewlyn/pressure/plymouthPressures.RData", envir=plymouth)

# Just use pressures over the periods of the model data (Pm) and TG data (Ptg)
plymouth$Pm <- plymouth$pressure[intersect(which(plymouth$time>=as.Date("1960/1/15")), which(plymouth$time<=as.Date("1999/12/15")))]
plymouth$Ptg <- plymouth$pressure[intersect(which(plymouth$time>=as.Date("1914/1/15")), which(plymouth$time<=as.Date("2006/12/15")))]

plymouth$data.brest <- data.frame(mmsl=brestMonthlyMean, msl=tg$mp.brest,
  tauX=as.vector(interp$brestMetMonthlyMeanInterp[,4,]), 
  tauY=as.vector(interp$brestMetMonthlyMeanInterp[,5,]), t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2, 
  Pm=plymouth$Pm)

plymouth$m.lm.brest <- lm(mmsl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)) + Pm, 
  data=plymouth$data.brest)
plymouth$tg.lm.brest <- lm(msl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)) + Pm, 
  data=plymouth$data.brest)
  
plymouth$m.resid.brest <- plymouth$m.lm.brest$resid
plymouth$tg.resid.brest <- plymouth$tg.lm.brest$resid
#> sd(plymouth$data.brest$mmsl, na.rm=T)
#[1] 0.07098115 cf 0.0676091 at Newlyn
#> sd(plymouth$m.resid.brest, na.rm=T)
#[1] 0.03311669 cf 0.03332927 at Newlyn

#> sd(plymouth$data.brest$msl, na.rm=T)
#[1] 87.1647 cf 82.10498 at Newlyn
#> sd(plymouth$tg.resid.brest, na.rm=T)
#[1] 44.49036 cf 41.12733 at Newlyn

x11()
par(family="HersheySans")
plot(time, plymouth$data.brest$mmsl, type='l', col='orange')
lines(time, (plymouth$m.lm.brest$coef[3]*cos(w1 * plymouth$data.brest$t)+plymouth$m.lm.brest$coef[4]*sin(w1 * plymouth$data.brest$t))+plymouth$m.lm.brest$coef[1], col='magenta', lwd=2)
lines(time, (plymouth$m.lm.brest$coef[5]*cos(w2 * plymouth$data.brest$t)+plymouth$m.lm.brest$coef[6]*sin(w2 * plymouth$data.brest$t))+plymouth$m.lm.brest$coef[1], col='cyan', lwd=2)
lines(time, plymouth$m.lm.brest$coef[7]*plymouth$data.brest$Pm +plymouth$m.lm.brest$coef[1], col='green', lwd=2)

plymouth$m.var.brest <- (plymouth$m.lm.brest$coef[3]*cos(w1 * plymouth$data.brest$t)+plymouth$m.lm.brest$coef[4]*sin(w1 * plymouth$data.brest$t)) +
                         (plymouth$m.lm.brest$coef[5]*cos(w2 * plymouth$data.brest$t)+plymouth$m.lm.brest$coef[6]*sin(w2 * plymouth$data.brest$t)) +
                         plymouth$m.lm.brest$coef[7]*plymouth$data.brest$Pm +
                         plymouth$m.lm.brest$coef[2]*plymouth$data.brest$t +
                         plymouth$m.lm.brest$coef[1]
plymouth$tg.var.brest <- (plymouth$tg.lm.brest$coef[3]*cos(w1 * plymouth$data.brest$t)+plymouth$tg.lm.brest$coef[4]*sin(w1 * plymouth$data.brest$t)) +
                         (plymouth$tg.lm.brest$coef[5]*cos(w2 * plymouth$data.brest$t)+plymouth$tg.lm.brest$coef[6]*sin(w2 * plymouth$data.brest$t)) +
                         plymouth$tg.lm.brest$coef[7]*plymouth$data.brest$Pm +
                         plymouth$tg.lm.brest$coef[2]*plymouth$data.brest$t +
                         plymouth$tg.lm.brest$coef[1]

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

plymouth$data.brest.tot <- data.frame(msl=tg$brestMonthly[1:1092], t=plymouth$tott, 
  Pm=plymouth$Ptg,  w1cos=cos(w1*plymouth$tott), w1sin=sin(w1*plymouth$tott), 
  w2cos=cos(w2*plymouth$tott), w2sin=sin(w2*plymouth$tott))
  
plymouth$data.brest.Ro <- data.frame(mslRo=tg$brestMonthly[1:588], 
  PtgRo=plymouth$PtgRossiter, tRo=plymouth$Rot, 
  w1cosRo=cos(w1*plymouth$Rot), w1sinRo=sin(w1*plymouth$Rot), 
  w2cosRo=cos(w2*plymouth$Rot), w2sinRo=sin(w2*plymouth$Rot))
  
plymouth$data.brest.Ca <- data.frame(mslCa=tg$brestMonthly[1:804],
  PtgCa=plymouth$PtgCartwright, tCa=plymouth$Cat,
  w1cosCa=cos(w1*plymouth$Cat), w1sinCa=sin(w1*plymouth$Cat), 
  w2cosCa=cos(w2*plymouth$Cat), w2sinCa=sin(w2*plymouth$Cat))
  
plymouth$tg.lm.brest.tot <- lm(msl ~ t + w1cos + w1sin + w2cos + w2sin + Pm,  
  data=plymouth$data.brest.tot)
plymouth$tg.resid.brest.tot <- plymouth$tg.lm.brest.tot$resid
#> sd(plymouth$data.brest.tot$msl, na.rm=T)
#[1] 92.76088 Newlyn: 91.69582
#> sd(plymouth$tg.resid.brest.tot, na.rm=T)
#[1] 48.26179 Newlyn: 40.46192
#> plymouth$tg.lm.brest.tot$coef[2]*12
#       t 
#1.264667 Newlyn: 1.822942 

plymouth$tg.lm.brest.Ro <- lm(mslRo ~ tRo + w1cosRo + w1sinRo + w2cosRo + w2sinRo + PtgRo,  
  data=plymouth$data.brest.Ro)
plymouth$tg.resid.brest.Ro <- plymouth$tg.lm.brest.Ro$resid
#> sd(plymouth$data.brest.Ro$msl, na.rm=T)
#[1] 86.38129 Newlyn: 83.08406
#> sd(plymouth$tg.resid.brest.Ro, na.rm=T)
#[1] 50.58848 Newlyn: 39.05382
#> plymouth$tg.lm.brest.Ro$coef[2]*12
#       t 
#1.368652 Newlyn: 2.273566
# cf 2.2 mm/yr in Rossiter (1972). This is effectively the same.

plymouth$tg.lm.brest.Ca <- lm(mslCa ~ tCa + w1cosCa + w1sinCa + w2cosCa + w2sinCa + PtgCa,  
  data=plymouth$data.brest.Ca)
plymouth$tg.resid.brest.Ca <- plymouth$tg.lm.brest.Ca$resid
#> sd(plymouth$data.brest.Ca$msl, na.rm=T)
#[1] 85.26331
#> sd(plymouth$tg.resid.brest.Ca, na.rm=T)
#[1] 40.21784
#> plymouth$tg.lm.brest.Ca$coef[2]*12
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
model2$Pm <- as.vector(interp$brestMetMonthlyMeanInterp[,1,])

model2$data.brest <- data.frame(mmsl=brestMonthlyMean, msl=tg$mp.brest,
  t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2, Pm=model2$Pm)

model2$m.lm.brest <- lm(mmsl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)) + Pm, 
  data=model2$data.brest)
model2$tg.lm.brest <- lm(msl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)) + Pm, 
  data=model2$data.brest)
  
model2$m.resid.brest <- model2$m.lm.brest$resid
model2$tg.resid.brest <- model2$tg.lm.brest$resid
#> sd(model2$data.brest$mmsl, na.rm=T)
#[1] 0.07098115 Newlyn: 0.0676091
#> sd(model2$m.resid.brest, na.rm=T)
#[1] 0.02752077 Newlyn: 0.0292884

#> sd(model2$data.brest$msl, na.rm=T)
#[1] 87.1647 Newlyn: 82.10498
#> sd(model2$tg.resid.brest, na.rm=T)
#[1] 39.22745 Newlyn: 35.91815

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

model3$data.brest <- data.frame(mmsl=brestMonthlyMean, msl=tg$mp.brest,
  t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2, Pm1=model3$Pm1, Pm2=model3$Pm2,
  Pm3=model3$Pm3, Pm4=model3$Pm4, Pm5=model3$Pm5, Pm6=model3$Pm6, Pm7=model3$Pm7,
  Pm8=model3$Pm8, Pm9=model3$Pm9)

model3$m.lm.brest <- lm(mmsl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)) + 
  Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9,  data=model3$data.brest)
model3$tg.lm.brest <- lm(msl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)) + 
  Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9, data=model3$data.brest)
  
model3$m.resid.brest <- model3$m.lm.brest$resid
model3$tg.resid.brest <- model3$tg.lm.brest$resid
#> sd(model3$data.brest$mmsl, na.rm=T)
#[1] 0.07098115 Newlyn: 0.0676091
#> sd(model3$m.resid.brest, na.rm=T)
#[1] 0.0241897 Newlyn: 0.02449552

#> sd(model3$data.brest$msl, na.rm=T)
#[1] 87.1647 Newlyn: 82.10498
#> sd(model3$tg.resid.brest, na.rm=T)
#[1] 35.15331 Newlyn: 27.06582

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
  model4$Pm[1,i] <- mean(as.vector(interp$extraMetMonthlyMeanInterp[,1,,(i-17)])[1:479], na.rm=T)
}
# Lag 2
for(i in 35:51){
  model4$Pm[3:480,i] <- as.vector(interp$extraMetMonthlyMeanInterp[,1,,(i-34)])[1:478]
  model4$Pm[1:2,i] <- mean(as.vector(interp$extraMetMonthlyMeanInterp[,1,,(i-34)])[1:478], na.rm=T)
}

# What we are going to do now is calculate the sd of the residuals as each of the
# 51 components is added in turn 
model4$m.sd <- vector(mode="numeric", length=51)
model4$tg.sd <- vector(mode="numeric", length=51)

for(i in 1:51){
  model4$data.brest <- data.frame(mmsl=brestMonthlyMean, msl=tg$mp.brest,
    t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2, Pm=model4$Pm[,i])

  model4$m.lm.brest <- lm(mmsl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) +
    I(sin(w2*t)) + Pm, data=model4$data.brest)
  model4$tg.lm.brest <- lm(msl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + 
    I(sin(w2*t)) + Pm, data=model4$data.brest)
  
  model4$m.resid.brest <- model4$m.lm.brest$resid
  model4$tg.resid.brest <- model4$tg.lm.brest$resid
  
  model4$m.sd[i] <- sd(model4$m.resid.brest, na.rm=T)
  model4$tg.sd[i] <- sd(model4$tg.resid.brest, na.rm=T)
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model4$m.sd.sorted <- sort(model4$m.sd, index=TRUE)
model4$tg.sd.sorted <- sort(model4$tg.sd, index=TRUE)

# Now construct models of the top 9 components
model4$data.brest <- data.frame(mmsl=brestMonthlyMean, msl=tg$mp.brest,
  t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2, Pm=model2$Pm,
  Pm1=model4$Pm[,model4$m.sd.sorted$ix[1]], 
  Pm2=model4$Pm[,model4$m.sd.sorted$ix[2]], Pm3=model4$Pm[,model4$m.sd.sorted$ix[3]],
  Pm4=model4$Pm[,model4$m.sd.sorted$ix[4]], Pm5=model4$Pm[,model4$m.sd.sorted$ix[5]],
  Pm6=model4$Pm[,model4$m.sd.sorted$ix[6]], Pm7=model4$Pm[,model4$m.sd.sorted$ix[7]],
  Pm8=model4$Pm[,model4$m.sd.sorted$ix[8]], Pm9=model4$Pm[,model4$m.sd.sorted$ix[9]])

model4$m.lm.brest <- lm(mmsl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)) + 
  Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 + Pm,  data=model4$data.brest)

model4$data.brest <- data.frame(mmsl=brestMonthlyMean, msl=tg$mp.brest,
  t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2, Pm=model2$Pm,
  Pm1=model4$Pm[,model4$tg.sd.sorted$ix[1]], 
  Pm2=model4$Pm[,model4$tg.sd.sorted$ix[2]], Pm3=model4$Pm[,model4$tg.sd.sorted$ix[3]],
  Pm4=model4$Pm[,model4$tg.sd.sorted$ix[4]], Pm5=model4$Pm[,model4$tg.sd.sorted$ix[5]],
  Pm6=model4$Pm[,model4$tg.sd.sorted$ix[6]], Pm7=model4$Pm[,model4$tg.sd.sorted$ix[7]],
  Pm8=model4$Pm[,model4$tg.sd.sorted$ix[8]], Pm9=model4$Pm[,model4$tg.sd.sorted$ix[9]])

model4$tg.lm.brest <- lm(msl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)) + 
  Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 + Pm,  data=model4$data.brest)

model4$m.resid.brest <- model4$m.lm.brest$resid
model4$tg.resid.brest <- model4$tg.lm.brest$resid

#model4$m.resid.brest <- model4$data.brest$mmsl - model4$m.lm.brest$coef[2]*model4$data.brest$t + model4$m.lm.brest$coef[1]
#model4$tg.resid.brest <- model4$data.brest$msl - model4$tg.lm.brest$coef[2]*model4$data.brest$t + model4$tg.lm.brest$coef[1]

#> sd(model4$data.brest$mmsl, na.rm=T)
#[1] 0.07098115 Newlyn: 0.0676091
#> sd(model4$m.resid.brest, na.rm=T)
#[1] 0.02322259 Newlyn: 0.02422649

# Variance reduction
# Brest
#> (0.07098115-0.02322259)/0.07098115
#[1] 0.6728344
# Newlyn
#> (0.0676091-0.02422649)/0.0676091
#[1] 0.6416682

#> sd(model4$data.brest$msl, na.rm=T)
#[1] 87.1647 Newlyn: 82.10498
#> sd(model4$tg.resid.brest, na.rm=T)
#[1] 33.23035 Newlyn: 25.16201

# Variance reduction
# Brest
#> (87.1647-33.23035)/87.1647
#[1] 0.6187637
# Newlyn
#> (82.10498-25.16201)/82.10498
#[1] 0.6935386

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
  model5$Ds[1,i] <- mean(deepSeaMonthlyMean[1:479,(i-6)])
}
# Lag 2
for(i in 13:18){
  model5$Ds[3:480,i] <- deepSeaMonthlyMean[1:478,(i-12)]
  model5$Ds[1:2,i] <- mean(deepSeaMonthlyMean[1:478,(i-12)])
}

# What we are going to do now is calculate the sd of the residuals as each of the
# 18 components is added in turn 
model5$m.sd <- vector(mode="numeric", length=18)
model5$tg.sd <- vector(mode="numeric", length=18)

for(i in 1:18){
  model5$data.brest <- data.frame(mmsl=brestMonthlyMean, msl=tg$mp.brest,
    t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2, Pm1=model4$Pm[,model4$m.sd.sorted$ix[1]], 
    Pm2=model4$Pm[,model4$m.sd.sorted$ix[2]], Pm3=model4$Pm[,model4$m.sd.sorted$ix[3]],
    Pm4=model4$Pm[,model4$m.sd.sorted$ix[4]], Pm5=model4$Pm[,model4$m.sd.sorted$ix[5]],
    Pm6=model4$Pm[,model4$m.sd.sorted$ix[6]], Pm7=model4$Pm[,model4$m.sd.sorted$ix[7]],
    Pm8=model4$Pm[,model4$m.sd.sorted$ix[8]], Pm9=model4$Pm[,model4$m.sd.sorted$ix[9]],
    dssl=model5$Ds[,i])

  model5$m.lm.brest <- lm(mmsl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)) + 
    Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 + dssl,  data=model5$data.brest)
    
  model5$tg.lm.brest <- lm(msl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)) + 
    Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 + dssl,  data=model5$data.brest)
  
  model5$m.resid.brest <- model5$m.lm.brest$resid
  model5$tg.resid.brest <- model5$tg.lm.brest$resid
  
  model5$m.sd[i] <- sd(model5$m.resid.brest, na.rm=T)
  model5$tg.sd[i] <- sd(model5$tg.resid.brest, na.rm=T)
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model5$m.sd.sorted <- sort(model5$m.sd, index=TRUE)
model5$tg.sd.sorted <- sort(model5$tg.sd, index=TRUE)

model5$t <- seq(from=0,to=((40*12)-1))

# Now construct models of the top 9 components
model5$data.brest <- data.frame(mmsl=brestMonthlyMean, msl=tg$mp.brest,
  t=seq(from=0,to=((40*12)-1)), w1cos=cos(w1*model5$t), w1sin=sin(w1*model5$t), w2cos=cos(w2*model5$t), w2sin=sin(w2*model5$t),
  Pm1=model4$Pm[,model4$m.sd.sorted$ix[1]], 
  Pm2=model4$Pm[,model4$m.sd.sorted$ix[2]], Pm3=model4$Pm[,model4$m.sd.sorted$ix[3]],
  Pm4=model4$Pm[,model4$m.sd.sorted$ix[4]], Pm5=model4$Pm[,model4$m.sd.sorted$ix[5]],
  Pm6=model4$Pm[,model4$m.sd.sorted$ix[6]], Pm7=model4$Pm[,model4$m.sd.sorted$ix[7]],
  Pm8=model4$Pm[,model4$m.sd.sorted$ix[8]], Pm9=model4$Pm[,model4$m.sd.sorted$ix[9]],
  Ds1=model5$Ds[,model5$m.sd.sorted$ix[1]], 
  Ds2=model5$Ds[,model5$m.sd.sorted$ix[2]], Ds3=model5$Ds[,model5$m.sd.sorted$ix[3]],
  Ds4=model5$Ds[,model5$m.sd.sorted$ix[4]], Ds5=model5$Ds[,model5$m.sd.sorted$ix[5]],
  Ds6=model5$Ds[,model5$m.sd.sorted$ix[6]], Ds7=model5$Ds[,model5$m.sd.sorted$ix[7]],
  Ds8=model5$Ds[,model5$m.sd.sorted$ix[8]], Ds9=model5$Ds[,model5$m.sd.sorted$ix[9]])

model5$m.lm.brest <- lm(mmsl ~ t + w1cos + w1sin + w2cos + w2sin + 
  Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 +
  Ds1 + Ds2 + Ds3 + Ds4 + Ds5 + Ds6 + Ds7 + Ds8 + Ds9,  data=model5$data.brest)

model5$data.brest <- data.frame(mmsl=brestMonthlyMean, msl=tg$mp.brest,
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

model5$tg.lm.brest <- lm(msl ~ t + w1cos + w1sin + w2cos + w2sin + 
  Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 +
  Ds1 + Ds2 + Ds3 + Ds4 + Ds5 + Ds6 + Ds7 + Ds8 + Ds9,  data=model5$data.brest)
  
model5$m.resid.brest <- model5$m.lm.brest$resid
model5$tg.resid.brest <- model5$tg.lm.brest$resid

model5$m.fitted.brest <- model5$m.lm.brest$fitted

#model5$m.fitted.brest <- model5$m.lm.brest$coef[1]
#for (i in 2:24){
#	model5$m.fitted.brest <- model5$m.fitted.brest + model5$m.lm.brest$coef[i]*model5$data.brest[[i+1]]
#}

x11()
plot(tg$mp.time,dm.mm.brestMonthlyMean, type='l', col='blue')
lines(tg$mp.time, (model5$m.fitted.brest-mean(model5$m.fitted.brest, na.rm=T))*1000, col='red')
lines(tg$mp.time, model5$m.resid.brest*1000, col='magenta', lwd=2)

model5$tg.fitted.brest <- model5$tg.lm.brest$fitted
#model5$tg.fitted.brest <- model5$tg.lm.brest$coef[1]
#for (i in 2:24){
#	model5$tg.fitted.brest <- model5$tg.fitted.brest + model5$tg.lm.brest$coef[i]*model5$data.brest[[i+1]]
#}

x11()
plot(tg$mp.time,tg$dm.mp.brest, type='l', col='blue')
lines(tg$mp.time, model5$tg.fitted.brest-mean(model5$tg.fitted.brest, na.rm=T), col='red')
lines(tg$mp.time, model5$tg.resid.brest, col='magenta', lwd=2)

#> sd(model5$data.brest$mmsl, na.rm=T)
#[1] 0.07098115 Newlyn: 0.0676091
#> sd(model5$m.resid.brest, na.rm=T)
#[1] 0.01703225 Newlyn: 0.01615011

#> sd(model5$data.brest$msl, na.rm=T)
#[1] 82.10498
#> sd(model5$tg.resid.brest, na.rm=T)
#[1] 30.08910 Newlyn: 23.78966

# Generate annual means from regression data
model5$numYrs <- 40

model5$annualYrs <- seq.Date(from=as.Date("1960/7/1"), to=as.Date("1999/7/1"), by="1 year")

tg$dm.mp.brest[which(is.na(tg$dm.mp.brest))] <- mean(tg$dm.mp.brest, na.rm=T)
model5$tmp <- tg$dm.mp.brest
dim(model5$tmp) <- c(12,model5$numYrs)
model5$tg.annualMeans.brest <- colMeans(model5$tmp)

model5$tmp <- model5$tg.fitted.brest-mean(model5$tg.fitted.brest, na.rm=T)
dim(model5$tmp) <- c(12,model5$numYrs)
model5$tg.fitted.annualMeans.brest <- colMeans(model5$tmp)

model5$tmp <- model5$tg.resid.brest
dim(model5$tmp) <- c(12,model5$numYrs)
model5$tg.resid.annualMeans.brest <- colMeans(model5$tmp)

# Model
model5$tmp <- dm.mm.brestMonthlyMean
dim(model5$tmp) <- c(12,model5$numYrs)
model5$m.annualMeans.brest <- colMeans(model5$tmp)

model5$tmp <- model5$m.fitted.brest-mean(model5$m.fitted.brest, na.rm=T)
dim(model5$tmp) <- c(12,model5$numYrs)
model5$m.fitted.annualMeans.brest <- colMeans(model5$tmp)*1000

model5$tmp <- model5$m.resid.brest
dim(model5$tmp) <- c(12,model5$numYrs)
model5$m.resid.annualMeans.brest <- colMeans(model5$tmp)*1000

rm(tmp, envir=model5)

x11()
plot(model5$annualYrs, model5$tg.annualMeans.brest, type='l', col='blue', lwd=2,
xlab="Year", ylab="Brest Tide Gauge SL [mm]")
lines(model5$annualYrs, model5$tg.fitted.annualMeans.brest, col='red', lwd=2)
lines(model5$annualYrs, model5$tg.resid.annualMeans.brest, col='magenta', lwd=2)
#var(model5$tg.annualMeans.brest, na.rm=T)
#[1] 1115.831 Newlyn: 807.7804
#var(model5$tg.resid.annualMeans.brest, na.rm=T)
#[1] 282.7684 Newlyn: 176.7980
#var(model5$tg.resid.annualMeans.brest, na.rm=T)/var(model5$tg.annualMeans.brest, na.rm=T)*100
#[1] 25.34151 Newlyn: 21.88689

x11()
plot(model5$annualYrs, model5$m.annualMeans.brest, type='l', col='blue', lwd=2,
		xlab="Year", ylab="Brest Model SL [mm]")
lines(model5$annualYrs, model5$m.fitted.annualMeans.brest, col='red', lwd=2)
lines(model5$annualYrs, model5$m.resid.annualMeans.brest, col='magenta', lwd=2)
#var(model5$m.annualMeans.brest, na.rm=T)
#[1] 441.1876 Newlyn: 470.1199
#var(model5$m.resid.annualMeans.brest, na.rm=T)
#[1] 88.4361 Newlyn: 90.13054
#var(model5$m.resid.annualMeans.brest, na.rm=T)/var(model5$m.annualMeans.brest, na.rm=T)*100
#[1] 20.04501 Newlyn: 19.17182

# Produce decadal means to examine reduction in decadal variance
model5$tg.tenYrTrends.brest <- vector(mode='numeric',length=(model5$numYrs-9))
model5$tg.resid.tenYrTrends.brest <- model5$tg.tenYrTrends.brest 
model5$tg.fitted.tenYrTrends.brest <- model5$tg.tenYrTrends.brest

model5$midpointYrs <- seq.Date(from=as.Date("1965/1/1"), length=(model5$numYrs-9), by="1 year")
for (i in 1:(model5$numYrs-9)) {
	if (length(which(is.finite(model5$tg.annualMeans.brest[i:(i+9)])))>=7){
		model5$junk <- lm(model5$tg.annualMeans.brest[i:(i+9)] ~ c(1:10), na.action = na.exclude)
		model5$tg.tenYrTrends.brest[i] <- model5$junk$coef[2]
	} else {
		model5$tg.tenYrTrends.brest[i] <- NA
	}
	if (length(which(is.finite(model5$tg.fitted.annualMeans.brest[i:(i+9)])))>=7){
		model5$junk <- lm(model5$tg.fitted.annualMeans.brest[i:(i+9)] ~ c(1:10), na.action = na.exclude)
		model5$tg.fitted.tenYrTrends.brest[i] <- model5$junk$coef[2]
	} else {
		model5$tg.fitted.tenYrTrends.brest[i] <- NA
	}
	if (length(which(is.finite(model5$tg.resid.annualMeans.brest[i:(i+9)])))>=7){
		model5$junk <- lm(model5$tg.resid.annualMeans.brest[i:(i+9)] ~ c(1:10), na.action = na.exclude)
		model5$tg.resid.tenYrTrends.brest[i] <- model5$junk$coef[2]
	} else {
		model5$tg.resid.tenYrTrends.brest[i] <- NA
	}
	
}

# Model
model5$m.tenYrTrends.brest <- vector(mode='numeric',length=(model5$numYrs-9))
model5$m.resid.tenYrTrends.brest <- model5$m.tenYrTrends.brest 
model5$m.fitted.tenYrTrends.brest <- model5$m.tenYrTrends.brest

for (i in 1:(model5$numYrs-9)) {
	if (length(which(is.finite(model5$m.annualMeans.brest[i:(i+9)])))>=7){
		model5$junk <- lm(model5$m.annualMeans.brest[i:(i+9)] ~ c(1:10), na.action = na.exclude)
		model5$m.tenYrTrends.brest[i] <- model5$junk$coef[2]
	} else {
		model5$m.tenYrTrends.brest[i] <- NA
	}
	if (length(which(is.finite(model5$m.fitted.annualMeans.brest[i:(i+9)])))>=7){
		model5$junk <- lm(model5$m.fitted.annualMeans.brest[i:(i+9)] ~ c(1:10), na.action = na.exclude)
		model5$m.fitted.tenYrTrends.brest[i] <- model5$junk$coef[2]
	} else {
		model5$m.fitted.tenYrTrends.brest[i] <- NA
	}
	if (length(which(is.finite(model5$m.resid.annualMeans.brest[i:(i+9)])))>=7){
		model5$junk <- lm(model5$m.resid.annualMeans.brest[i:(i+9)] ~ c(1:10), na.action = na.exclude)
		model5$m.resid.tenYrTrends.brest[i] <- model5$junk$coef[2]
	} else {
		model5$m.resid.tenYrTrends.brest[i] <- NA
	}
	
}

rm(junk, envir=model5)

x11()
plot(model5$midpointYrs, model5$tg.tenYrTrends.brest, type='l', col='blue', lwd=2,
		xlab="Year", ylab="Brest Tide Gauge SL [mm]")
lines(model5$midpointYrs, model5$tg.fitted.tenYrTrends.brest, col='red', lwd=2)
lines(model5$midpointYrs, model5$tg.resid.tenYrTrends.brest, col='magenta', lwd=2)
#var(model5$tg.tenYrTrends.brest, na.rm=T)
#[1] 19.10664 Newlyn: 15.47344
#var(model5$tg.resid.tenYrTrends.brest, na.rm=T)
#[1] 8.113846 Newlyn: 5.570029
#var(model5$tg.resid.tenYrTrends.brest, na.rm=T)/var(model5$tg.tenYrTrends.brest, na.rm=T)*100
#[1] 42.4661 Newlyn: 35.99736

x11()
plot(model5$midpointYrs, model5$m.tenYrTrends.brest, type='l', col='blue', lwd=2,
		xlab="Year", ylab="Brest Model SL [mm]")
lines(model5$midpointYrs, model5$m.fitted.tenYrTrends.brest, col='red', lwd=2)
lines(model5$midpointYrs, model5$m.resid.tenYrTrends.brest, col='magenta', lwd=2)
#var(model5$m.tenYrTrends.brest, na.rm=T)
#[1] 6.492875 Newlyn: 7.017346
#var(model5$m.resid.tenYrTrends.brest, na.rm=T)
#[1] 1.196288 Newlyn: 1.048865
#var(model5$m.resid.tenYrTrends.brest, na.rm=T)/var(model5$m.tenYrTrends.brest, na.rm=T)*100
#[1] 18.42463 Newlyn: 14.94675

##############
## Model 4b ##
##############
# As with model 4 but using POLCOMS run S12R408
# S12R408 was run for longer, so we need to truncate it at 480 months
# 
model4b <- new.env() 
# Again, just use pressures over the periods of the model data (Pm)

# In Thompson's model 4 he also introduces lags of 0, 1 & 2 months which in this
# case gives a total of 51 components. However, to facilitate comparison with
# model 3 he only uses the first 9 components which he chooses by stepwise
# regression. We achieve this by adding in each component in turn and calculating
# the variance reduction, which we can then rank and analyze.
model4b$Pm <- array(NA,dim=c(480,51))
# Zero lag
for(i in 1:17){
	model4b$Pm[,i] <- as.vector(interp$extraMetMonthlyMeanInterp[,1,,i])
}
# Lag 1
for(i in 18:34){
	model4b$Pm[2:480,i] <- as.vector(interp$extraMetMonthlyMeanInterp[,1,,(i-17)])[1:479]
	model4b$Pm[1,i] <- mean(as.vector(interp$extraMetMonthlyMeanInterp[,1,,(i-17)])[1:479], na.rm=T)
}
# Lag 2
for(i in 35:51){
	model4b$Pm[3:480,i] <- as.vector(interp$extraMetMonthlyMeanInterp[,1,,(i-34)])[1:478]
	model4b$Pm[1:2,i] <- mean(as.vector(interp$extraMetMonthlyMeanInterp[,1,,(i-34)])[1:478], na.rm=T)
}

# What we are going to do now is calculate the sd of the residuals as each of the
# 51 components is added in turn 
model4b$m.sd <- vector(mode="numeric", length=51)
model4b$tg.sd <- vector(mode="numeric", length=51)

for(i in 1:51){
	model4b$data.brest <- data.frame(mmsl=S12R408$brestMonthlyMean[1:480,1], msl=tg$mp.brest,
			t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2, Pm=model4b$Pm[,i])
	
	model4b$m.lm.brest <- lm(mmsl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) +
					I(sin(w2*t)) + Pm, data=model4b$data.brest)
	model4b$tg.lm.brest <- lm(msl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + 
					I(sin(w2*t)) + Pm, data=model4b$data.brest)
	
	model4b$m.resid.brest <- model4b$m.lm.brest$resid
	model4b$tg.resid.brest <- model4b$tg.lm.brest$resid
	
	model4b$m.sd[i] <- sd(model4b$m.resid.brest, na.rm=T)
	model4b$tg.sd[i] <- sd(model4b$tg.resid.brest, na.rm=T)
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model4b$m.sd.sorted <- sort(model4b$m.sd, index=TRUE)
model4b$tg.sd.sorted <- sort(model4b$tg.sd, index=TRUE)

# Now construct models of the top 9 components
model4b$data.brest <- data.frame(mmsl=S12R408$brestMonthlyMean[1:480,1], msl=tg$mp.brest,
		t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2, Pm=model2$Pm,
		Pm1=model4b$Pm[,model4b$m.sd.sorted$ix[1]], 
		Pm2=model4b$Pm[,model4b$m.sd.sorted$ix[2]], Pm3=model4b$Pm[,model4b$m.sd.sorted$ix[3]],
		Pm4=model4b$Pm[,model4b$m.sd.sorted$ix[4]], Pm5=model4b$Pm[,model4b$m.sd.sorted$ix[5]],
		Pm6=model4b$Pm[,model4b$m.sd.sorted$ix[6]], Pm7=model4b$Pm[,model4b$m.sd.sorted$ix[7]],
		Pm8=model4b$Pm[,model4b$m.sd.sorted$ix[8]], Pm9=model4b$Pm[,model4b$m.sd.sorted$ix[9]])

model4b$m.lm.brest <- lm(mmsl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)) + 
				Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 + Pm,  data=model4b$data.brest)

model4b$data.brest <- data.frame(mmsl=S12R408$brestMonthlyMean[1:480,1], msl=tg$mp.brest,
		t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2, Pm=model2$Pm,
		Pm1=model4b$Pm[,model4b$tg.sd.sorted$ix[1]], 
		Pm2=model4b$Pm[,model4b$tg.sd.sorted$ix[2]], Pm3=model4b$Pm[,model4b$tg.sd.sorted$ix[3]],
		Pm4=model4b$Pm[,model4b$tg.sd.sorted$ix[4]], Pm5=model4b$Pm[,model4b$tg.sd.sorted$ix[5]],
		Pm6=model4b$Pm[,model4b$tg.sd.sorted$ix[6]], Pm7=model4b$Pm[,model4b$tg.sd.sorted$ix[7]],
		Pm8=model4b$Pm[,model4b$tg.sd.sorted$ix[8]], Pm9=model4b$Pm[,model4b$tg.sd.sorted$ix[9]])

model4b$tg.lm.brest <- lm(msl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)) + 
				Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 + Pm,  data=model4b$data.brest)

model4b$m.resid.brest <- model4b$m.lm.brest$resid
model4b$tg.resid.brest <- model4b$tg.lm.brest$resid

#model4b$m.resid.brest <- model4b$data.brest$mmsl - model4b$m.lm.brest$coef[2]*model4b$data.brest$t + model4b$m.lm.brest$coef[1]
#model4b$tg.resid.brest <- model4b$data.brest$msl - model4b$tg.lm.brest$coef[2]*model4b$data.brest$t + model4b$tg.lm.brest$coef[1]

#> sd(model4b$data.brest$mmsl, na.rm=T)
#[1] 0.07098115 Newlyn: [1] 0.0676091 (model4)
#[1] 0.07975017 Newlyn: [1] 0.07649279 (model4b)

#> sd(model4b$m.resid.brest, na.rm=T)
#[1] 0.02322259 Newlyn: [1] 0.02422649 (model4)
#[1] 0.02352657 Newlyn: [1] 0.06616694 (model4b)

#> sd(model4b$data.brest$msl, na.rm=T)
#[1] 86.98253 Newlyn: 82.10498
#> sd(model4b$tg.resid.brest, na.rm=T)
#[1] 33.23898 Newlyn: 25.16201
# No signficant improvement

#############
## Model 5b ##
#############
# Same as model 5 but use POLCOMS run S12R408. 
model5b <- new.env()
model5b$Ds <- array(NA,dim=c(480,18))
# Zero lag
for(i in 1:6){
	model5b$Ds[,i] <- S12R408$deepSeaMonthlyMean[1:480,i]
}
# Lag 1
for(i in 7:12){
	model5b$Ds[2:480,i] <- S12R408$deepSeaMonthlyMean[1:479,(i-6)]
	model5b$Ds[1,i] <- mean(S12R408$deepSeaMonthlyMean[1:479,(i-6)], na.rm=T)
}
# Lag 2
for(i in 13:18){
	model5b$Ds[3:480,i] <- S12R408$deepSeaMonthlyMean[1:478,(i-12)]
	model5b$Ds[1:2,i] <- mean(S12R408$deepSeaMonthlyMean[1:478,(i-12)], na.rm=T)
}

# What we are going to do now is calculate the sd of the residuals as each of the
# 18 components is added in turn 
model5b$m.sd <- vector(mode="numeric", length=18)
model5b$tg.sd <- vector(mode="numeric", length=18)

for(i in 1:18){
	model5b$data.brest <- data.frame(mmsl=S12R408$brestMonthlyMean[1:480,1], msl=tg$mp.brest,
			t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2, Pm1=model4b$Pm[,model4b$m.sd.sorted$ix[1]], 
			Pm2=model4b$Pm[,model4b$m.sd.sorted$ix[2]], Pm3=model4b$Pm[,model4b$m.sd.sorted$ix[3]],
			Pm4=model4b$Pm[,model4b$m.sd.sorted$ix[4]], Pm5=model4b$Pm[,model4b$m.sd.sorted$ix[5]],
			Pm6=model4b$Pm[,model4b$m.sd.sorted$ix[6]], Pm7=model4b$Pm[,model4b$m.sd.sorted$ix[7]],
			Pm8=model4b$Pm[,model4b$m.sd.sorted$ix[8]], Pm9=model4b$Pm[,model4b$m.sd.sorted$ix[9]],
			dssl=model5b$Ds[,i])
	
	model5b$m.lm.brest <- lm(mmsl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)) + 
					Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 + dssl,  data=model5b$data.brest)
	
	model5b$tg.lm.brest <- lm(msl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)) + 
					Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 + dssl,  data=model5b$data.brest)
	
	model5b$m.resid.brest <- model5b$m.lm.brest$resid
	model5b$tg.resid.brest <- model5b$tg.lm.brest$resid
	
	model5b$m.sd[i] <- sd(model5b$m.resid.brest, na.rm=T)
	model5b$tg.sd[i] <- sd(model5b$tg.resid.brest, na.rm=T)
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model5b$m.sd.sorted <- sort(model5b$m.sd, index=TRUE)
model5b$tg.sd.sorted <- sort(model5b$tg.sd, index=TRUE)

model5b$t <- seq(from=0,to=((40*12)-1))

# Now construct models of the top 9 components
model5b$data.brest <- data.frame(mmsl=S12R408$brestMonthlyMean[1:480,1], msl=tg$mp.brest,
		t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2, Pm1=model4b$Pm[,model4b$m.sd.sorted$ix[1]], 
		Pm2=model4b$Pm[,model4b$m.sd.sorted$ix[2]], Pm3=model4b$Pm[,model4b$m.sd.sorted$ix[3]],
		Pm4=model4b$Pm[,model4b$m.sd.sorted$ix[4]], Pm5=model4b$Pm[,model4b$m.sd.sorted$ix[5]],
		Pm6=model4b$Pm[,model4b$m.sd.sorted$ix[6]], Pm7=model4b$Pm[,model4b$m.sd.sorted$ix[7]],
		Pm8=model4b$Pm[,model4b$m.sd.sorted$ix[8]], Pm9=model4b$Pm[,model4b$m.sd.sorted$ix[9]],
		Ds1=model5b$Ds[,model5b$m.sd.sorted$ix[1]], 
		Ds2=model5b$Ds[,model5b$m.sd.sorted$ix[2]], Ds3=model5b$Ds[,model5b$m.sd.sorted$ix[3]],
		Ds4=model5b$Ds[,model5b$m.sd.sorted$ix[4]], Ds5=model5b$Ds[,model5b$m.sd.sorted$ix[5]],
		Ds6=model5b$Ds[,model5b$m.sd.sorted$ix[6]], Ds7=model5b$Ds[,model5b$m.sd.sorted$ix[7]],
		Ds8=model5b$Ds[,model5b$m.sd.sorted$ix[8]], Ds9=model5b$Ds[,model5b$m.sd.sorted$ix[9]],
		w1cos=cos(w1*model5b$t), w1sin=sin(w1*model5b$t), 
		w2cos=cos(w2*model5b$t), w2sin=sin(w2*model5b$t))

model5b$m.lm.brest <- lm(mmsl ~ t + w1cos + w1sin + w2cos + w2sin + 
				Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 +
				Ds1 + Ds2 + Ds3 + Ds4 + Ds5 + Ds6 + Ds7 + Ds8 + Ds9,  data=model5b$data.brest)

model5b$data.brest <- data.frame(mmsl=S12R408$brestMonthlyMean[1:480,1], msl=tg$mp.brest,
		t=model5b$t, w1cos=cos(w1*model5b$t), w1sin=sin(w1*model5b$t), w2cos=cos(w2*model5b$t), w2sin=sin(w2*model5b$t), 
		Pm1=model4b$Pm[,model4b$tg.sd.sorted$ix[1]], 
		Pm2=model4b$Pm[,model4b$tg.sd.sorted$ix[2]], Pm3=model4b$Pm[,model4b$tg.sd.sorted$ix[3]],
		Pm4=model4b$Pm[,model4b$tg.sd.sorted$ix[4]], Pm5=model4b$Pm[,model4b$tg.sd.sorted$ix[5]],
		Pm6=model4b$Pm[,model4b$tg.sd.sorted$ix[6]], Pm7=model4b$Pm[,model4b$tg.sd.sorted$ix[7]],
		Pm8=model4b$Pm[,model4b$tg.sd.sorted$ix[8]], Pm9=model4b$Pm[,model4b$tg.sd.sorted$ix[9]],
		Ds1=model5b$Ds[,model5b$tg.sd.sorted$ix[1]], 
		Ds2=model5b$Ds[,model5b$tg.sd.sorted$ix[2]], Ds3=model5b$Ds[,model5b$tg.sd.sorted$ix[3]],
		Ds4=model5b$Ds[,model5b$tg.sd.sorted$ix[4]], Ds5=model5b$Ds[,model5b$tg.sd.sorted$ix[5]],
		Ds6=model5b$Ds[,model5b$tg.sd.sorted$ix[6]], Ds7=model5b$Ds[,model5b$tg.sd.sorted$ix[7]],
		Ds8=model5b$Ds[,model5b$tg.sd.sorted$ix[8]], Ds9=model5b$Ds[,model5b$tg.sd.sorted$ix[9]],
		w1cos=cos(w1*model5b$t), w1sin=sin(w1*model5b$t), 
		w2cos=cos(w2*model5b$t), w2sin=sin(w2*model5b$t))

model5b$tg.lm.brest <- lm(msl ~ t + w1cos + w1sin + w2cos + w2sin + 
				Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 +
				Ds1 + Ds2 + Ds3 + Ds4 + Ds5 + Ds6 + Ds7 + Ds8 + Ds9,  data=model5b$data.brest)

model5b$m.resid.brest <- model5b$m.lm.brest$resid
model5b$tg.resid.brest <- model5b$tg.lm.brest$resid

model5b$m.fitted.brest <- model5b$m.lm.brest$fitted
#model5b$m.fitted.brest <- model5b$m.lm.brest$coef[1]
#for (i in 2:24){
#	model5b$m.fitted.brest <- model5b$m.fitted.brest + model5b$m.lm.brest$coef[i]*model5b$data.brest[[i+1]]
#}
#model5b$m.resid.brest <- model5b$data.brest$mmsl - model5b$m.fitted.brest

x11()
plot(tg$mp.time,dm.mm.brestMonthlyMean, type='l', col='blue')
lines(tg$mp.time, (model5b$m.fitted.brest-mean(model5b$m.fitted.brest, na.rm=T))*1000, col='red')
lines(tg$mp.time, model5b$m.resid.brest*1000, col='magenta', lwd=2)

model5b$tg.fitted.brest <- model5b$tg.lm.brest$fitted
#model5b$tg.fitted.brest <- model5b$tg.lm.brest$coef[1]
#for (i in 2:24){
#	model5b$tg.fitted.brest <- model5b$tg.fitted.brest + model5b$tg.lm.brest$coef[i]*model5b$data.brest[[i+1]]
#}
#model5b$tg.resid.brest <- model5b$data.brest$msl - model5b$tg.fitted.brest

x11()
plot(tg$mp.time,tg$dm.mp.brest, type='l', col='blue')
lines(tg$mp.time, model5b$tg.fitted.brest-mean(model5b$tg.fitted.brest, na.rm=T), col='red')
lines(tg$mp.time, model5b$tg.resid.brest, col='magenta', lwd=2)

#> sd(model5b$data.brest$mmsl, na.rm=T)
#[1] 0.07098115 Newlyn: 0.0676091 (model5)
#[1] 0.07975017 Newlyn: 0.07649279 (model5b)
#> sd(model5b$m.resid.brest, na.rm=T)
#[1] 0.01703225 Newlyn: 0.01615011 (model5)
#[1] 0.01717341 Newlyn: 0.2276618 (model5b)

#> sd(model5b$data.brest$msl, na.rm=T)
#[1] 86.98253 Newlyn: 82.10498
#> sd(model5b$tg.resid.brest, na.rm=T)
#[1] 30.08910 Newlyn: 23.78966 (model5)
#[1] 30.37137 Newlyn: 23.35431 (model5b)

# Generate annual means from regression data
model5b$numYrs <- 40

model5b$annualYrs <- seq.Date(from=as.Date("1960/7/1"), to=as.Date("1999/7/1"), by="1 year")

model5b$tmp <- tg$dm.mp.brest
dim(model5b$tmp) <- c(12,model5b$numYrs)
model5b$tg.annualMeans.brest <- colMeans(model5b$tmp)

model5b$tmp <- model5b$tg.fitted.brest-mean(model5b$tg.fitted.brest, na.rm=T)
dim(model5b$tmp) <- c(12,model5b$numYrs)
model5b$tg.fitted.annualMeans.brest <- colMeans(model5b$tmp)

model5b$tmp <- model5b$tg.resid.brest
dim(model5b$tmp) <- c(12,model5b$numYrs)
model5b$tg.resid.annualMeans.brest <- colMeans(model5b$tmp)

rm(tmp, envir=model5b)

x11()
plot(model5b$annualYrs, model5b$tg.annualMeans.brest, type='l', col='blue')
lines(model5b$annualYrs, model5b$tg.fitted.annualMeans.brest, col='red')
lines(model5b$annualYrs, model5b$tg.resid.annualMeans.brest, col='magenta', lwd=2)
#var(model5b$tg.annualMeans.brest, na.rm=T)
#[1] 1115.831 Newlyn: 807.7804
#var(model5b$tg.resid.annualMeans.brest, na.rm=T)
#[1] 282.7684 Newlyn: 176.7980 (model5)
#[1] 281.4484 Newlyn: 159.8537 (model5b)
#var(model5b$tg.resid.annualMeans.brest, na.rm=T)/var(model5b$tg.annualMeans.brest, na.rm=T)*100
#[1] 25.34151 Newlyn: 21.88689 (model5)
#[1] 25.22321 Newlyn: 19.78925 (model5b)

# Produce decadal means to examine reduction in decadal variance
model5b$tg.tenYrTrends.brest <- vector(mode='numeric',length=(model5b$numYrs-9))
model5b$tg.resid.tenYrTrends.brest <- model5b$tg.tenYrTrends.brest 
model5b$tg.fitted.tenYrTrends.brest <- model5b$tg.tenYrTrends.brest

model5b$midpointYrs <- seq.Date(from=as.Date("1965/1/1"), length=(model5b$numYrs-9), by="1 year")
for (i in 1:(model5b$numYrs-9)) {
	if (length(which(is.finite(model5b$tg.annualMeans.brest[i:(i+9)])))>=7){
		model5b$junk <- lm(model5b$tg.annualMeans.brest[i:(i+9)] ~ c(1:10), na.action = na.exclude)
		model5b$tg.tenYrTrends.brest[i] <- model5b$junk$coef[2]
	} else {
		model5b$tg.tenYrTrends.brest[i] <- NA
	}
	if (length(which(is.finite(model5b$tg.fitted.annualMeans.brest[i:(i+9)])))>=7){
		model5b$junk <- lm(model5b$tg.fitted.annualMeans.brest[i:(i+9)] ~ c(1:10), na.action = na.exclude)
		model5b$tg.fitted.tenYrTrends.brest[i] <- model5b$junk$coef[2]
	} else {
		model5b$tg.fitted.tenYrTrends.brest[i] <- NA
	}
	if (length(which(is.finite(model5b$tg.resid.annualMeans.brest[i:(i+9)])))>=7){
		model5b$junk <- lm(model5b$tg.resid.annualMeans.brest[i:(i+9)] ~ c(1:10), na.action = na.exclude)
		model5b$tg.resid.tenYrTrends.brest[i] <- model5b$junk$coef[2]
	} else {
		model5b$tg.resid.tenYrTrends.brest[i] <- NA
	}
	
}
rm(junk, envir=model5b)

x11()
plot(model5b$midpointYrs, model5b$tg.tenYrTrends.brest, type='l', col='blue')
lines(model5b$midpointYrs, model5b$tg.fitted.tenYrTrends.brest, col='red')
lines(model5b$midpointYrs, model5b$tg.resid.tenYrTrends.brest, col='magenta', lwd=2)
#var(model5b$tg.tenYrTrends.brest, na.rm=T)
#[1] 19.10664 Newlyn: 15.47344
#var(model5b$tg.resid.tenYrTrends.brest, na.rm=T)
#[1] 8.113846 Newlyn: 5.570029 (model5)
#[1] 6.452357 Newlyn: 3.963082 (model5b)
#var(model5b$tg.resid.tenYrTrends.brest, na.rm=T)/var(model5b$tg.tenYrTrends.brest, na.rm=T)*100
#[1] 42.4661 Newlyn: 35.99736 (model5)
#[1] 33.7702 Newlyn: 25.61216 (model5b)

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

model6$NewlynEs <- as.vector(interp$brestMetMonthlyMeanInterp[,2,])
model6$NewlynNs <- as.vector(interp$brestMetMonthlyMeanInterp[,3,])

# Zero lag
for(i in 1:17){
  model6$Es[,i] <- as.vector(interp$extraMetMonthlyMeanInterp[,2,,i])
  model6$Ns[,i] <- as.vector(interp$extraMetMonthlyMeanInterp[,3,,i])
}
# Lag 1
for(i in 18:34){
  model6$Es[2:480,i] <- as.vector(interp$extraMetMonthlyMeanInterp[,2,,(i-17)])[1:479]
  model6$Ns[2:480,i] <- as.vector(interp$extraMetMonthlyMeanInterp[,3,,(i-17)])[1:479]
  model6$Es[1,i] <- mean(as.vector(interp$extraMetMonthlyMeanInterp[,2,,(i-17)])[1:479], na.rm=T)
  model6$Ns[1,i] <- mean(as.vector(interp$extraMetMonthlyMeanInterp[,3,,(i-17)])[1:479], na.rm=T)
  
}
# Lag 2
for(i in 35:51){
  model6$Es[3:480,i] <- as.vector(interp$extraMetMonthlyMeanInterp[,2,,(i-34)])[1:478]
  model6$Ns[3:480,i] <- as.vector(interp$extraMetMonthlyMeanInterp[,3,,(i-34)])[1:478]
  model6$Es[1:2,i] <- mean(as.vector(interp$extraMetMonthlyMeanInterp[,2,,(i-34)])[1:478], na.rm=T)
  model6$Ns[1:2,i] <- mean(as.vector(interp$extraMetMonthlyMeanInterp[,3,,(i-34)])[1:478], na.rm=T)
}

# What we are going to do now is calculate the sd of the residuals as each of the
# 51 components is added in turn 
model6$m.sd <- vector(mode="numeric", length=51)
model6$tg.sd <- vector(mode="numeric", length=51)

for(i in 1:51){
  model6$data.brest <- data.frame(mmsl=brestMonthlyMean, msl=tg$mp.brest,
    t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2, Pm=model2$Pm,
    NEs=model6$NewlynEs, NNs=model6$NewlynNs, Es=model6$Es[,i], Ns=model6$Ns[,i])

  model6$m.lm.brest <- lm(mmsl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) +
    I(sin(w2*t)) + Es + Ns + NEs + NNs + Pm, data=model6$data.brest)
  model6$tg.lm.brest <- lm(msl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + 
    I(sin(w2*t)) + Es + Ns + NEs + NNs + Pm, data=model6$data.brest)
  
  model6$m.resid.brest <- model6$m.lm.brest$resid
  model6$tg.resid.brest <- model6$tg.lm.brest$resid
  
  model6$m.sd[i] <- sd(model6$m.resid.brest, na.rm=T)
  model6$tg.sd[i] <- sd(model6$tg.resid.brest, na.rm=T)
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model6$m.sd.sorted <- sort(model6$m.sd, index=TRUE)
model6$tg.sd.sorted <- sort(model6$tg.sd, index=TRUE)

# Now construct models of the top 9 components
model6$t <- seq(from=0,to=((40*12)-1))
model6$data.brest <- data.frame(mmsl=brestMonthlyMean, msl=tg$mp.brest,
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

model6$m.lm.brest <- lm(mmsl ~ t + w1cos + w1sin + w2cos + w2sin + 
  Es1 + Es2 + Es3 + Es4 + Es5 + Es6 + Es7 + Es8 + Es9 +
  Ns1 + Ns2 + Ns3 + Ns4 + Ns5 + Ns6 + Ns7 + Ns8 + Ns9 +
  NEs + NNs + Pm,  data=model6$data.brest)

model6$data.brest <- data.frame(mmsl=brestMonthlyMean, msl=tg$mp.brest,
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

model6$tg.lm.brest <- lm(msl ~ t + w1cos + w1sin + w2cos + w2sin + 
  Es1 + Es2 + Es3 + Es4 + Es5 + Es6 + Es7 + Es8 + Es9 +
  Ns1 + Ns2 + Ns3 + Ns4 + Ns5 + Ns6 + Ns7 + Ns8 + Ns9 +
  NEs + NNs + Pm,  data=model6$data.brest)

model6$m.resid.brest <- model6$m.lm.brest$resid
model6$tg.resid.brest <- model6$tg.lm.brest$resid

#> sd(model6$data.brest$mmsl, na.rm=T)
#[1] 0.07098115 Newlyn: 0.0676091
#> sd(model6$m.resid.brest, na.rm=T)
#[1] 0.02187177 Newlyn: 0.02276734

#> sd(model6$data.brest$msl, na.rm=T)
#[1] 86.98253 Newlyn: 82.10498
#> sd(model6$tg.resid.brest, na.rm=T)
#[1] 30.78034 Newlyn: 22.99741

# So wind stresses (or just winds) actually do worse than just static pressures in the
# model only domain (though not in the "real" domain). Why?
# Dynamic height in the model is not the same as sea level so I should diagnose and use a
# linear combination of diagnostic heights first.
# Jason also suggests using the currents. Although this isn't a very useful tool
# for reconstructing sea level (unless we have a model) it is useful as a diagnostic tool
# to tell us what is actually controlling SL at the coast.

#> model6$tg.lm.brest$coef[2]*12
#t 
#1.841724 mm/yr 1960-2000 after correction

##########
## CS3X ##
##########
# Look at CS3X data and see if it is any different to POLCOMS
# Columns of brest, brest and deep sea (ds) are z,u,v
# The 3rd dimension of ds covers the 6 points:
# [1,50] [1,100] [1,150] [50,50] [50,100] [50,150]
#
# CS3X data runs from 15/1/58 to 15/12/01 (44yrs, 528 months)

cs3x <- new.env() 
load("~/diskx/cs3x/brestNewlynCS3X.RData", envir=cs3x)

# CS3X is over 44 years so we need a different TG data set
# Compare just over the period of the model
tg$cs3x.brest <- tg$brestMonthly[intersect(which(tg$time>=as.Date("1958/1/15")), which(tg$time<=as.Date("2001/12/15")))]
tg$cs3x.brest <- tg$brestMonthly[intersect(which(tg$time>=as.Date("1958/1/15")), which(tg$time<=as.Date("2001/12/15")))]
tg$cs3x.time <- tg$time[intersect(which(tg$time>=as.Date("1958/1/15")), which(tg$time<=as.Date("2001/12/15")))]


cs3x$Ds <- array(NA,dim=c(528,18))
# Elevation: zero lag
for(i in 1:6){
	cs3x$Ds[,i] <- cs3x$ds[,1,i]
}
# Lag 1
for(i in 7:12){
	cs3x$Ds[2:528,i] <- cs3x$ds[1:527,1,(i-6)]
}
# Lag 2
for(i in 13:18){
	cs3x$Ds[3:528,i] <- cs3x$ds[1:526,1,(i-12)]
}

# What we are going to do now is calculate the sd of the residuals as each of the
# 18 components is added in turn 
cs3x$m.sd <- vector(mode="numeric", length=18)
cs3x$tg.sd <- vector(mode="numeric", length=18)

for(i in 1:18){
	cs3x$data.brest <- data.frame(mmsl=brestMonthlyMean, msl=tg$mp.brest,
			t=seq(from=0,to=((44*12)-1)), w1=w1, w2=w2, 
			w1cos=cos(w1*cs3x$t), w1sin=sin(w1*cs3x$t), w2cos=cos(w2*cs3x$t), w2sin=sin(w2*cs3x$t),
			Pm1=model4$Pm[,model4$m.sd.sorted$ix[1]], 
			Pm2=model4$Pm[,model4$m.sd.sorted$ix[2]], Pm3=model4$Pm[,model4$m.sd.sorted$ix[3]],
			Pm4=model4$Pm[,model4$m.sd.sorted$ix[4]], Pm5=model4$Pm[,model4$m.sd.sorted$ix[5]],
			Pm6=model4$Pm[,model4$m.sd.sorted$ix[6]], Pm7=model4$Pm[,model4$m.sd.sorted$ix[7]],
			Pm8=model4$Pm[,model4$m.sd.sorted$ix[8]], Pm9=model4$Pm[,model4$m.sd.sorted$ix[9]],
			dssl=cs3x$Ds[,i])
	
	cs3x$m.lm.brest <- lm(mmsl ~ t + w1cos + w1sin + w2cos + w2sin + 
					Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 + dssl,  data=cs3x$data.brest)
	
	cs3x$tg.lm.brest <- lm(msl ~ t + w1cos + w1sin + w2cos + w2sin + 
					Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 + dssl,  data=cs3x$data.brest)
	
	cs3x$m.resid.brest <- cs3x$m.lm.brest$resid
	cs3x$tg.resid.brest <- cs3x$tg.lm.brest$resid
	
	cs3x$m.sd[i] <- sd(cs3x$m.resid.brest, na.rm=T)
	cs3x$tg.sd[i] <- sd(cs3x$tg.resid.brest, na.rm=T)
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
cs3x$m.sd.sorted <- sort(cs3x$m.sd, index=TRUE)
cs3x$tg.sd.sorted <- sort(cs3x$tg.sd, index=TRUE)

cs3x$t <- seq(from=0,to=((40*12)-1))

# Now construct models of the top 9 components
cs3x$data.brest <- data.frame(mmsl=brestMonthlyMean, msl=tg$mp.brest,
		t=seq(from=0,to=((40*12)-1)), w1cos=cos(w1*cs3x$t), w1sin=sin(w1*cs3x$t), w2cos=cos(w2*cs3x$t), w2sin=sin(w2*cs3x$t),
		w1=w1, w2=w2, 
		Pm1=model4$Pm[,model4$m.sd.sorted$ix[1]], 
		Pm2=model4$Pm[,model4$m.sd.sorted$ix[2]], Pm3=model4$Pm[,model4$m.sd.sorted$ix[3]],
		Pm4=model4$Pm[,model4$m.sd.sorted$ix[4]], Pm5=model4$Pm[,model4$m.sd.sorted$ix[5]],
		Pm6=model4$Pm[,model4$m.sd.sorted$ix[6]], Pm7=model4$Pm[,model4$m.sd.sorted$ix[7]],
		Pm8=model4$Pm[,model4$m.sd.sorted$ix[8]], Pm9=model4$Pm[,model4$m.sd.sorted$ix[9]],
		Ds1=cs3x$Ds[,cs3x$m.sd.sorted$ix[1]], 
		Ds2=cs3x$Ds[,cs3x$m.sd.sorted$ix[2]], Ds3=cs3x$Ds[,cs3x$m.sd.sorted$ix[3]],
		Ds4=cs3x$Ds[,cs3x$m.sd.sorted$ix[4]], Ds5=cs3x$Ds[,cs3x$m.sd.sorted$ix[5]],
		Ds6=cs3x$Ds[,cs3x$m.sd.sorted$ix[6]], Ds7=cs3x$Ds[,cs3x$m.sd.sorted$ix[7]],
		Ds8=cs3x$Ds[,cs3x$m.sd.sorted$ix[8]], Ds9=cs3x$Ds[,cs3x$m.sd.sorted$ix[9]])

cs3x$m.lm.brest <- lm(mmsl ~ t + w1cos + w1sin + w2cos + w2sin + 
				Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 +
				Ds1 + Ds2 + Ds3 + Ds4 + Ds5 + Ds6 + Ds7 + Ds8 + Ds9,  data=cs3x$data.brest)

cs3x$data.brest <- data.frame(mmsl=brestMonthlyMean, msl=tg$mp.brest,
		t=cs3x$t, w1cos=cos(w1*cs3x$t), w1sin=sin(w1*cs3x$t), w2cos=cos(w2*cs3x$t), w2sin=sin(w2*cs3x$t), 
		Pm1=model4$Pm[,model4$tg.sd.sorted$ix[1]], 
		Pm2=model4$Pm[,model4$tg.sd.sorted$ix[2]], Pm3=model4$Pm[,model4$tg.sd.sorted$ix[3]],
		Pm4=model4$Pm[,model4$tg.sd.sorted$ix[4]], Pm5=model4$Pm[,model4$tg.sd.sorted$ix[5]],
		Pm6=model4$Pm[,model4$tg.sd.sorted$ix[6]], Pm7=model4$Pm[,model4$tg.sd.sorted$ix[7]],
		Pm8=model4$Pm[,model4$tg.sd.sorted$ix[8]], Pm9=model4$Pm[,model4$tg.sd.sorted$ix[9]],
		Ds1=cs3x$Ds[,cs3x$tg.sd.sorted$ix[1]], 
		Ds2=cs3x$Ds[,cs3x$tg.sd.sorted$ix[2]], Ds3=cs3x$Ds[,cs3x$tg.sd.sorted$ix[3]],
		Ds4=cs3x$Ds[,cs3x$tg.sd.sorted$ix[4]], Ds5=cs3x$Ds[,cs3x$tg.sd.sorted$ix[5]],
		Ds6=cs3x$Ds[,cs3x$tg.sd.sorted$ix[6]], Ds7=cs3x$Ds[,cs3x$tg.sd.sorted$ix[7]],
		Ds8=cs3x$Ds[,cs3x$tg.sd.sorted$ix[8]], Ds9=cs3x$Ds[,cs3x$tg.sd.sorted$ix[9]])

cs3x$tg.lm.brest <- lm(msl ~ t + w1cos + w1sin + w2cos + w2sin + 
				Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 +
				Ds1 + Ds2 + Ds3 + Ds4 + Ds5 + Ds6 + Ds7 + Ds8 + Ds9,  data=cs3x$data.brest)

#cs3x$m.resid.brest <- cs3x$m.lm.brest$resid
#cs3x$tg.resid.brest <- cs3x$tg.lm.brest$resid

cs3x$m.fitted.brest <- cs3x$m.lm.brest$coef[1]
for (i in 2:24){
	cs3x$m.fitted.brest <- cs3x$m.fitted.brest + cs3x$m.lm.brest$coef[i]*cs3x$data.brest[[i+1]]
}
cs3x$m.resid.brest <- cs3x$data.brest$mmsl - cs3x$m.fitted.brest

x11()
plot(tg$mp.time,dm.mm.brestMonthlyMean, type='l', col='blue')
lines(tg$mp.time, (cs3x$m.fitted.brest-mean(cs3x$m.fitted.brest, na.rm=T))*1000, col='red')
lines(tg$mp.time, cs3x$m.resid.brest*1000, col='magenta', lwd=2)


cs3x$tg.fitted.brest <- cs3x$tg.lm.brest$coef[1]
for (i in 2:24){
	cs3x$tg.fitted.brest <- cs3x$tg.fitted.brest + cs3x$tg.lm.brest$coef[i]*cs3x$data.brest[[i+1]]
}
cs3x$tg.resid.brest <- cs3x$data.brest$msl - cs3x$tg.fitted.brest

x11()
plot(tg$mp.time,tg$dm.mp.brest, type='l', col='blue')
lines(tg$mp.time, cs3x$tg.fitted.brest-mean(cs3x$tg.fitted.brest, na.rm=T), col='red')
lines(tg$mp.time, cs3x$tg.resid.brest, col='magenta', lwd=2)

#> sd(cs3x$data.brest$mmsl, na.rm=T)
#[1] 0.0676091
#> sd(cs3x$m.resid.brest, na.rm=T)
#[1] 0.01615011

#> sd(cs3x$data.brest$msl, na.rm=T)
#[1] 82.10498
#> sd(cs3x$tg.resid.brest, na.rm=T)
#[1] 23.78966

# Generate annual means from regression data
cs3x$numYrs <- 44

cs3x$annualYrs <- seq.Date(from=as.Date("1960/7/1"), to=as.Date("1999/7/1"), by="1 year")

cs3x$tmp <- tg$dm.mp.brest
dim(cs3x$tmp) <- c(12,cs3x$numYrs)
cs3x$tg.annualMeans.brest <- colMeans(cs3x$tmp)

cs3x$tmp <- cs3x$tg.fitted.brest-mean(cs3x$tg.fitted.brest, na.rm=T)
dim(cs3x$tmp) <- c(12,cs3x$numYrs)
cs3x$tg.fitted.annualMeans.brest <- colMeans(cs3x$tmp)

cs3x$tmp <- cs3x$tg.resid.brest
dim(cs3x$tmp) <- c(12,cs3x$numYrs)
cs3x$tg.resid.annualMeans.brest <- colMeans(cs3x$tmp)

rm(tmp, envir=cs3x)

x11()
plot(cs3x$annualYrs, cs3x$tg.annualMeans.brest, type='l', col='blue')
lines(cs3x$annualYrs, cs3x$tg.fitted.annualMeans.brest, col='red')
lines(cs3x$annualYrs, cs3x$tg.resid.annualMeans.brest, col='magenta', lwd=2)
#var(cs3x$tg.annualMeans.brest, na.rm=T)
#[1] 807.7804
#var(cs3x$tg.resid.annualMeans.brest, na.rm=T)
#[1] 176.7980
#var(cs3x$tg.resid.annualMeans.brest, na.rm=T)/var(cs3x$tg.annualMeans.brest, na.rm=T)*100
#[1] 21.88689

# Produce decadal means to examine reduction in decadal variance
cs3x$tg.tenYrTrends.brest <- vector(mode='numeric',length=(cs3x$numYrs-9))
cs3x$tg.resid.tenYrTrends.brest <- cs3x$tg.tenYrTrends.brest 
cs3x$tg.fitted.tenYrTrends.brest <- cs3x$tg.tenYrTrends.brest

cs3x$midpointYrs <- seq.Date(from=as.Date("1965/1/1"), length=(cs3x$numYrs-9), by="1 year")
for (i in 1:(cs3x$numYrs-9)) {
	if (length(which(is.finite(cs3x$tg.annualMeans.brest[i:(i+9)])))>=7){
		cs3x$junk <- lm(cs3x$tg.annualMeans.brest[i:(i+9)] ~ c(1:10), na.action = na.exclude)
		cs3x$tg.tenYrTrends.brest[i] <- cs3x$junk$coef[2]
	} else {
		cs3x$tg.tenYrTrends.brest[i] <- NA
	}
	if (length(which(is.finite(cs3x$tg.fitted.annualMeans.brest[i:(i+9)])))>=7){
		cs3x$junk <- lm(cs3x$tg.fitted.annualMeans.brest[i:(i+9)] ~ c(1:10), na.action = na.exclude)
		cs3x$tg.fitted.tenYrTrends.brest[i] <- cs3x$junk$coef[2]
	} else {
		cs3x$tg.fitted.tenYrTrends.brest[i] <- NA
	}
	if (length(which(is.finite(cs3x$tg.resid.annualMeans.brest[i:(i+9)])))>=7){
		cs3x$junk <- lm(cs3x$tg.resid.annualMeans.brest[i:(i+9)] ~ c(1:10), na.action = na.exclude)
		cs3x$tg.resid.tenYrTrends.brest[i] <- cs3x$junk$coef[2]
	} else {
		cs3x$tg.resid.tenYrTrends.brest[i] <- NA
	}
	
}
rm(junk, envir=cs3x)

x11()
plot(cs3x$midpointYrs, cs3x$tg.tenYrTrends.brest, type='l', col='blue')
lines(cs3x$midpointYrs, cs3x$tg.fitted.tenYrTrends.brest, col='red')
lines(cs3x$midpointYrs, cs3x$tg.resid.tenYrTrends.brest, col='magenta', lwd=2)
#var(cs3x$tg.tenYrTrends.brest, na.rm=T)
#[1] 15.47344
#var(cs3x$tg.resid.tenYrTrends.brest, na.rm=T)
#[1] 5.570029
#var(cs3x$tg.resid.tenYrTrends.brest, na.rm=T)/var(cs3x$tg.tenYrTrends.brest, na.rm=T)*100
#[1] 35.99736

# Use HadSLP2r pressure data
hadley<-new.env()
load("~/diskx/polcoms/brestNewlyn/pressure/brestNewlynHadSLP2.RData", envir=hadley)
# We need to work in terms of total pressure rather than in sea level as in Thompson '86
hadley$totp.tg.mp.brest <- 1025*9.8*tg$mp.brest/1000
hadley$totp.tg.brest <- 1025*9.8*tg$brestMonthly/1000
hadley$totp.m.brest <- 1025*9.8*S12R408$brestMonthlyMean[1:480,1]
#
colnames(hadley$slpNewlynStns) <- c("Newlyn", "Brest", "DS1", "DS2", "DS3", "DS4", "DS5", "DS6", "StMawgan","ThorneyIsland","Shoeburyness","Gorleston","Kilnsea","Eskdalemuir","Kinloss","Ronaldsway","Mumbles","EM1","EM2","EM3","EM4","EM5","EM6","EM7","EM8")
# Some of these columns are duplicates due to the large (5deg) grid size of HadSLP2r. The unique columns are:
hadley$uniq <- which(!duplicated(hadley$slpNewlynStns[1,]))
# [1]  1  3  4  5  6  7  8 10 13 14 15 20 21 22 23 24
hadley$slpStns<-hadley$slpNewlynStns[,which(!duplicated(hadley$slpNewlynStns[1,]))]

# HadSLP2r runs from Jan 1850 to Dec 2007
hadley$time <- seq.Date(from=as.Date("1850/1/15"), to=as.Date("2007/12/15"), by="1 month")

# Compare over the period of the model
hadley$mp.slp <- hadley$slpStns[intersect(which(hadley$time>=as.Date("1960/1/15")), which(hadley$time<=as.Date("1999/12/15"))),]
hadley$mp.slpNewlyn <- hadley$slpNewlynStns[intersect(which(hadley$time>=as.Date("1960/1/15")), which(hadley$time<=as.Date("1999/12/15"))),]
# Compare over the period of Newlyn
hadley$n.slp <- hadley$slpStns[intersect(which(hadley$time>=as.Date("1914/1/15")), which(hadley$time<=as.Date("2006/12/15"))),]

hadley$totp.tg.mp.brest <- hadley$totp.tg.mp.brest + hadley$mp.slp[,1]
hadley$totp.tg.brest <- hadley$totp.tg.brest + hadley$n.slp[,1]
hadley$totp.m.brest <- hadley$totp.m.brest + hadley$mp.slp[,1]

##############
## Model 4c ##
##############
# As with model 4b (using S12R408) but using HadSLP2r data.
# Also, put local pressure term on the LHS
# to avoid aliasing and use total pressure.
# S12R408 was run for longer, so we need to truncate it at 480 months
# 
model4c <- new.env() 
# Again, just use pressures over the periods of the model data (Pm)

# In Thompson's model 4 he introduces lags of 0, 1 & 2 months. In this
# case, due to duplication of pressure points (leaving local + 6 deep sea points and 9
# extra met points) gives a total of 27 components (Extra Met only).
# However, to facilitate comparison with
# model 3 he only uses the first 9 components which he chooses by stepwise
# regression. We achieve this by adding in each component in turn and calculating
# the variance reduction, which we can then rank and analyze.
model4c$Pm <- array(NA,dim=c(480,27))
# Zero lag
for(i in 1:9){
	model4c$Pm[,i] <- as.vector(hadley$mp.slp[,(i+7)])
}
# Lag 1
for(i in 10:18){
	model4c$Pm[2:480,i] <- as.vector(hadley$mp.slp[,(i-2)])[1:479]
	model4c$Pm[1,i] <- mean(as.vector(hadley$mp.slp[,(i-2)])[1:479], na.rm=T)
}
# Lag 2
for(i in 19:27){
	model4c$Pm[3:480,i] <- as.vector(hadley$mp.slp[,(i-11)])[1:478]
	model4c$Pm[1:2,i] <- mean(as.vector(hadley$mp.slp[,(i-11)])[1:478], na.rm=T)
}

# What we are going to do now is calculate the sd of the residuals as each of the
# 27 components is added in turn 
model4c$m.sd <- vector(mode="numeric", length=27)
model4c$tg.sd <- vector(mode="numeric", length=27)
model4c$t <- seq(from=0,to=((40*12)-1))

for(i in 1:27){
	model4c$data.brest <- data.frame(mmsl=hadley$totp.m.brest, msl=hadley$totp.tg.mp.brest,
			t=model4c$t, w1cos=cos(w1*model4c$t), w1sin=sin(w1*model4c$t), 
		        w2cos=cos(w2*model4c$t), w2sin=sin(w2*model4c$t), Pm=model4c$Pm[,i])
	
	model4c$m.lm.brest <- lm(mmsl ~ t + w1cos + w1sin + w2cos + w2sin + Pm, data=model4c$data.brest)
	model4c$tg.lm.brest <- lm(msl ~ t + w1cos + w1sin + w2cos + w2sin + Pm, data=model4c$data.brest)
	
	model4c$m.resid.brest <- model4c$m.lm.brest$resid
	model4c$tg.resid.brest <- model4c$tg.lm.brest$resid
	
	model4c$m.sd[i] <- sd(model4c$m.resid.brest, na.rm=T)
	model4c$tg.sd[i] <- sd(model4c$tg.resid.brest, na.rm=T)
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model4c$m.sd.sorted <- sort(model4c$m.sd, index=TRUE)
model4c$tg.sd.sorted <- sort(model4c$tg.sd, index=TRUE)

# Now construct models of the top 9 components
model4c$data.brest <- data.frame(mmsl=hadley$totp.m.brest, msl=hadley$totp.tg.mp.brest,
			t=model4c$t, w1cos=cos(w1*model4c$t), w1sin=sin(w1*model4c$t), 
		        w2cos=cos(w2*model4c$t), w2sin=sin(w2*model4c$t), Pm1=model4c$Pm[,model4c$m.sd.sorted$ix[1]], 
                        Pm2=model4c$Pm[,model4c$m.sd.sorted$ix[2]], Pm3=model4c$Pm[,model4c$m.sd.sorted$ix[3]],
		        Pm4=model4c$Pm[,model4c$m.sd.sorted$ix[4]], Pm5=model4c$Pm[,model4c$m.sd.sorted$ix[5]],
		        Pm6=model4c$Pm[,model4c$m.sd.sorted$ix[6]], Pm7=model4c$Pm[,model4c$m.sd.sorted$ix[7]],
		        Pm8=model4c$Pm[,model4c$m.sd.sorted$ix[8]], Pm9=model4c$Pm[,model4c$m.sd.sorted$ix[9]])

model4c$m.lm.brest <- lm(mmsl ~ t + w1cos + w1sin + w2cos + w2sin + 
				Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9,  data=model4c$data.brest)

model4c$data.brest <- data.frame(mmsl=hadley$totp.m.brest, msl=hadley$totp.tg.mp.brest,
			t=model4c$t, w1cos=cos(w1*model4c$t), w1sin=sin(w1*model4c$t), 
		        w2cos=cos(w2*model4c$t), w2sin=sin(w2*model4c$t), Pm1=model4c$Pm[,model4c$tg.sd.sorted$ix[1]], 
		        Pm2=model4c$Pm[,model4c$tg.sd.sorted$ix[2]], Pm3=model4c$Pm[,model4c$tg.sd.sorted$ix[3]],
		        Pm4=model4c$Pm[,model4c$tg.sd.sorted$ix[4]], Pm5=model4c$Pm[,model4c$tg.sd.sorted$ix[5]],
		        Pm6=model4c$Pm[,model4c$tg.sd.sorted$ix[6]], Pm7=model4c$Pm[,model4c$tg.sd.sorted$ix[7]],
		        Pm8=model4c$Pm[,model4c$tg.sd.sorted$ix[8]], Pm9=model4c$Pm[,model4c$tg.sd.sorted$ix[9]])

model4c$tg.lm.brest <- lm(msl ~ t + w1cos + w1sin + w2cos + w2sin + 
				Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9,  data=model4c$data.brest)

model4c$m.resid.brest <- model4c$m.lm.brest$resid
model4c$tg.resid.brest <- model4c$tg.lm.brest$resid

#model4c$m.resid.brest <- model4c$data.brest$mmsl - model4c$m.lm.brest$coef[2]*model4c$data.brest$t + model4c$m.lm.brest$coef[1]
#model4c$tg.resid.brest <- model4c$data.brest$msl - model4c$tg.lm.brest$coef[2]*model4c$data.brest$t + model4c$tg.lm.brest$coef[1]

#> sd(model4c$data.brest$mmsl/(1025*9.8), na.rm=T)
#[1] 0.07098115 Newlyn: [1] 0.0676091 (model4)
#[1] 0.07972603 Newlyn: [1] 0.07649279 (model4c)

#> sd(model4c$m.resid.brest/(1025*9.8), na.rm=T)
#[1] 0.02322259 Newlyn: [1] 0.02422649 (model4)
#[1] 0.0655804 Newlyn: [1] 0.06616694 (model4c)

#> sd(model4c$data.brest$msl/(1025*9.8), na.rm=T)
#[1] 0.08690541 Newlyn: 82.10498
#> sd(model4c$tg.resid.brest/(1025*9.8), na.rm=T)
#[1] 0.07343209 Newlyn: 25.16201
# No signficant improvement

#############
## Model 5c ##
#############
# Same as model 5 but use POLCOMS run S12R408. 
model5c <- new.env()
model5c$Ds <- array(NA,dim=c(480,18))
# Zero lag
for(i in 1:6){
	model5c$Ds[,i] <- S12R408$deepSeaMonthlyMean[1:480,i]
}
# Lag 1
for(i in 7:12){
	model5c$Ds[2:480,i] <- S12R408$deepSeaMonthlyMean[1:479,(i-6)]
	model5c$Ds[1,i] <- mean(S12R408$deepSeaMonthlyMean[1:479,(i-6)], na.rm=T)
}
# Lag 2
for(i in 13:18){
	model5c$Ds[3:480,i] <- S12R408$deepSeaMonthlyMean[1:478,(i-12)]
	model5c$Ds[1:2,i] <- mean(S12R408$deepSeaMonthlyMean[1:478,(i-12)], na.rm=T)
}

# What we are going to do now is calculate the sd of the residuals as each of the
# 18 components is added in turn 
model5c$m.sd <- vector(mode="numeric", length=18)
model5c$tg.sd <- vector(mode="numeric", length=18)

for(i in 1:18){
	model5c$data.brest <- data.frame(mmsl=S12R408$brestMonthlyMean[1:480,1], msl=tg$mp.brest,
			t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2, Pm1=model4c$Pm[,model4c$m.sd.sorted$ix[1]], 
			Pm2=model4c$Pm[,model4c$m.sd.sorted$ix[2]], Pm3=model4c$Pm[,model4c$m.sd.sorted$ix[3]],
			Pm4=model4c$Pm[,model4c$m.sd.sorted$ix[4]], Pm5=model4c$Pm[,model4c$m.sd.sorted$ix[5]],
			Pm6=model4c$Pm[,model4c$m.sd.sorted$ix[6]], Pm7=model4c$Pm[,model4c$m.sd.sorted$ix[7]],
			Pm8=model4c$Pm[,model4c$m.sd.sorted$ix[8]], Pm9=model4c$Pm[,model4c$m.sd.sorted$ix[9]],
			dssl=model5c$Ds[,i])
	
	model5c$m.lm.brest <- lm(mmsl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)) + 
					Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 + dssl,  data=model5c$data.brest)
	
	model5c$tg.lm.brest <- lm(msl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)) + 
					Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 + dssl,  data=model5c$data.brest)
	
	model5c$m.resid.brest <- model5c$m.lm.brest$resid
	model5c$tg.resid.brest <- model5c$tg.lm.brest$resid
	
	model5c$m.sd[i] <- sd(model5c$m.resid.brest, na.rm=T)
	model5c$tg.sd[i] <- sd(model5c$tg.resid.brest, na.rm=T)
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model5c$m.sd.sorted <- sort(model5c$m.sd, index=TRUE)
model5c$tg.sd.sorted <- sort(model5c$tg.sd, index=TRUE)

model5c$t <- seq(from=0,to=((40*12)-1))

# Now construct models of the top 9 components
model5c$data.brest <- data.frame(mmsl=S12R408$brestMonthlyMean[1:480,1], msl=tg$mp.brest,
		t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2, Pm1=model4c$Pm[,model4c$m.sd.sorted$ix[1]], 
		Pm2=model4c$Pm[,model4c$m.sd.sorted$ix[2]], Pm3=model4c$Pm[,model4c$m.sd.sorted$ix[3]],
		Pm4=model4c$Pm[,model4c$m.sd.sorted$ix[4]], Pm5=model4c$Pm[,model4c$m.sd.sorted$ix[5]],
		Pm6=model4c$Pm[,model4c$m.sd.sorted$ix[6]], Pm7=model4c$Pm[,model4c$m.sd.sorted$ix[7]],
		Pm8=model4c$Pm[,model4c$m.sd.sorted$ix[8]], Pm9=model4c$Pm[,model4c$m.sd.sorted$ix[9]],
		Ds1=model5c$Ds[,model5c$m.sd.sorted$ix[1]], 
		Ds2=model5c$Ds[,model5c$m.sd.sorted$ix[2]], Ds3=model5c$Ds[,model5c$m.sd.sorted$ix[3]],
		Ds4=model5c$Ds[,model5c$m.sd.sorted$ix[4]], Ds5=model5c$Ds[,model5c$m.sd.sorted$ix[5]],
		Ds6=model5c$Ds[,model5c$m.sd.sorted$ix[6]], Ds7=model5c$Ds[,model5c$m.sd.sorted$ix[7]],
		Ds8=model5c$Ds[,model5c$m.sd.sorted$ix[8]], Ds9=model5c$Ds[,model5c$m.sd.sorted$ix[9]],
		w1cos=cos(w1*model5c$t), w1sin=sin(w1*model5c$t), 
		w2cos=cos(w2*model5c$t), w2sin=sin(w2*model5c$t))

model5c$m.lm.brest <- lm(mmsl ~ t + w1cos + w1sin + w2cos + w2sin + 
				Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 +
				Ds1 + Ds2 + Ds3 + Ds4 + Ds5 + Ds6 + Ds7 + Ds8 + Ds9,  data=model5c$data.brest)

model5c$data.brest <- data.frame(mmsl=S12R408$brestMonthlyMean[1:480,1], msl=tg$mp.brest,
		t=model5c$t, w1cos=cos(w1*model5c$t), w1sin=sin(w1*model5c$t), w2cos=cos(w2*model5c$t), w2sin=sin(w2*model5c$t), 
		Pm1=model4c$Pm[,model4c$tg.sd.sorted$ix[1]], 
		Pm2=model4c$Pm[,model4c$tg.sd.sorted$ix[2]], Pm3=model4c$Pm[,model4c$tg.sd.sorted$ix[3]],
		Pm4=model4c$Pm[,model4c$tg.sd.sorted$ix[4]], Pm5=model4c$Pm[,model4c$tg.sd.sorted$ix[5]],
		Pm6=model4c$Pm[,model4c$tg.sd.sorted$ix[6]], Pm7=model4c$Pm[,model4c$tg.sd.sorted$ix[7]],
		Pm8=model4c$Pm[,model4c$tg.sd.sorted$ix[8]], Pm9=model4c$Pm[,model4c$tg.sd.sorted$ix[9]],
		Ds1=model5c$Ds[,model5c$tg.sd.sorted$ix[1]], 
		Ds2=model5c$Ds[,model5c$tg.sd.sorted$ix[2]], Ds3=model5c$Ds[,model5c$tg.sd.sorted$ix[3]],
		Ds4=model5c$Ds[,model5c$tg.sd.sorted$ix[4]], Ds5=model5c$Ds[,model5c$tg.sd.sorted$ix[5]],
		Ds6=model5c$Ds[,model5c$tg.sd.sorted$ix[6]], Ds7=model5c$Ds[,model5c$tg.sd.sorted$ix[7]],
		Ds8=model5c$Ds[,model5c$tg.sd.sorted$ix[8]], Ds9=model5c$Ds[,model5c$tg.sd.sorted$ix[9]],
		w1cos=cos(w1*model5c$t), w1sin=sin(w1*model5c$t), 
		w2cos=cos(w2*model5c$t), w2sin=sin(w2*model5c$t))

model5c$tg.lm.brest <- lm(msl ~ t + w1cos + w1sin + w2cos + w2sin + 
				Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 +
				Ds1 + Ds2 + Ds3 + Ds4 + Ds5 + Ds6 + Ds7 + Ds8 + Ds9,  data=model5c$data.brest)

model5c$m.resid.brest <- model5c$m.lm.brest$resid
model5c$tg.resid.brest <- model5c$tg.lm.brest$resid

model5c$m.fitted.brest <- model5c$m.lm.brest$fitted
#model5c$m.fitted.brest <- model5c$m.lm.brest$coef[1]
#for (i in 2:24){
#	model5c$m.fitted.brest <- model5c$m.fitted.brest + model5c$m.lm.brest$coef[i]*model5c$data.brest[[i+1]]
#}
#model5c$m.resid.brest <- model5c$data.brest$mmsl - model5c$m.fitted.brest

x11()
plot(tg$mp.time,dm.mm.brestMonthlyMean, type='l', col='blue')
lines(tg$mp.time, (model5c$m.fitted.brest-mean(model5c$m.fitted.brest, na.rm=T))*1000, col='red')
lines(tg$mp.time, model5c$m.resid.brest*1000, col='magenta', lwd=2)

model5c$tg.fitted.brest <- model5c$tg.lm.brest$fitted
#model5c$tg.fitted.brest <- model5c$tg.lm.brest$coef[1]
#for (i in 2:24){
#	model5c$tg.fitted.brest <- model5c$tg.fitted.brest + model5c$tg.lm.brest$coef[i]*model5c$data.brest[[i+1]]
#}
#model5c$tg.resid.brest <- model5c$data.brest$msl - model5c$tg.fitted.brest

x11()
plot(tg$mp.time,tg$dm.mp.brest, type='l', col='blue')
lines(tg$mp.time, model5c$tg.fitted.brest-mean(model5c$tg.fitted.brest, na.rm=T), col='red')
lines(tg$mp.time, model5c$tg.resid.brest, col='magenta', lwd=2)

#> sd(model5c$data.brest$mmsl, na.rm=T)
#[1] 0.07098115 Newlyn: 0.0676091 (model5)
#[1] 0.07975017 Newlyn: 0.07649279 (model5c)
#> sd(model5c$m.resid.brest, na.rm=T)
#[1] 0.01703225 Newlyn: 0.01615011 (model5)
#[1] 0.01717341 Newlyn: 0.2276618 (model5c)

#> sd(model5c$data.brest$msl, na.rm=T)
#[1] 86.98253 Newlyn: 82.10498
#> sd(model5c$tg.resid.brest, na.rm=T)
#[1] 30.08910 Newlyn: 23.78966 (model5)
#[1] 30.37137 Newlyn: 23.35431 (model5c)

# Generate annual means from regression data
model5c$numYrs <- 40

model5c$annualYrs <- seq.Date(from=as.Date("1960/7/1"), to=as.Date("1999/7/1"), by="1 year")

model5c$tmp <- tg$dm.mp.brest
dim(model5c$tmp) <- c(12,model5c$numYrs)
model5c$tg.annualMeans.brest <- colMeans(model5c$tmp)

model5c$tmp <- model5c$tg.fitted.brest-mean(model5c$tg.fitted.brest, na.rm=T)
dim(model5c$tmp) <- c(12,model5c$numYrs)
model5c$tg.fitted.annualMeans.brest <- colMeans(model5c$tmp)

model5c$tmp <- model5c$tg.resid.brest
dim(model5c$tmp) <- c(12,model5c$numYrs)
model5c$tg.resid.annualMeans.brest <- colMeans(model5c$tmp)

rm(tmp, envir=model5c)

x11()
plot(model5c$annualYrs, model5c$tg.annualMeans.brest, type='l', col='blue')
lines(model5c$annualYrs, model5c$tg.fitted.annualMeans.brest, col='red')
lines(model5c$annualYrs, model5c$tg.resid.annualMeans.brest, col='magenta', lwd=2)
#var(model5c$tg.annualMeans.brest, na.rm=T)
#[1] 1115.831 Newlyn: 807.7804
#var(model5c$tg.resid.annualMeans.brest, na.rm=T)
#[1] 282.7684 Newlyn: 176.7980 (model5)
#[1] 281.4484 Newlyn: 159.8537 (model5c)
#var(model5c$tg.resid.annualMeans.brest, na.rm=T)/var(model5c$tg.annualMeans.brest, na.rm=T)*100
#[1] 25.34151 Newlyn: 21.88689 (model5)
#[1] 25.22321 Newlyn: 19.78925 (model5c)

# Produce decadal means to examine reduction in decadal variance
model5c$tg.tenYrTrends.brest <- vector(mode='numeric',length=(model5c$numYrs-9))
model5c$tg.resid.tenYrTrends.brest <- model5c$tg.tenYrTrends.brest 
model5c$tg.fitted.tenYrTrends.brest <- model5c$tg.tenYrTrends.brest

model5c$midpointYrs <- seq.Date(from=as.Date("1965/1/1"), length=(model5c$numYrs-9), by="1 year")
for (i in 1:(model5c$numYrs-9)) {
	if (length(which(is.finite(model5c$tg.annualMeans.brest[i:(i+9)])))>=7){
		model5c$junk <- lm(model5c$tg.annualMeans.brest[i:(i+9)] ~ c(1:10), na.action = na.exclude)
		model5c$tg.tenYrTrends.brest[i] <- model5c$junk$coef[2]
	} else {
		model5c$tg.tenYrTrends.brest[i] <- NA
	}
	if (length(which(is.finite(model5c$tg.fitted.annualMeans.brest[i:(i+9)])))>=7){
		model5c$junk <- lm(model5c$tg.fitted.annualMeans.brest[i:(i+9)] ~ c(1:10), na.action = na.exclude)
		model5c$tg.fitted.tenYrTrends.brest[i] <- model5c$junk$coef[2]
	} else {
		model5c$tg.fitted.tenYrTrends.brest[i] <- NA
	}
	if (length(which(is.finite(model5c$tg.resid.annualMeans.brest[i:(i+9)])))>=7){
		model5c$junk <- lm(model5c$tg.resid.annualMeans.brest[i:(i+9)] ~ c(1:10), na.action = na.exclude)
		model5c$tg.resid.tenYrTrends.brest[i] <- model5c$junk$coef[2]
	} else {
		model5c$tg.resid.tenYrTrends.brest[i] <- NA
	}
	
}
rm(junk, envir=model5c)

x11()
plot(model5c$midpointYrs, model5c$tg.tenYrTrends.brest, type='l', col='blue')
lines(model5c$midpointYrs, model5c$tg.fitted.tenYrTrends.brest, col='red')
lines(model5c$midpointYrs, model5c$tg.resid.tenYrTrends.brest, col='magenta', lwd=2)
#var(model5c$tg.tenYrTrends.brest, na.rm=T)
#[1] 19.10664 Newlyn: 15.47344
#var(model5c$tg.resid.tenYrTrends.brest, na.rm=T)
#[1] 8.113846 Newlyn: 5.570029 (model5)
#[1] 6.452357 Newlyn: 3.963082 (model5c)
#var(model5c$tg.resid.tenYrTrends.brest, na.rm=T)/var(model5c$tg.tenYrTrends.brest, na.rm=T)*100
#[1] 42.4661 Newlyn: 35.99736 (model5)
#[1] 33.7702 Newlyn: 25.61216 (model5c)

