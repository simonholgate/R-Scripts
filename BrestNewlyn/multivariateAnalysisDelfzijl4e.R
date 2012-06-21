# multivariateAnalysisDelfzijl4e.R
# 
# Implementation of Thompson's methods (Thompson, 1980; 1986) to the tide gauge
# time series at Delfzijl
# These are the same methods implemented in multivariateAnalysisDelfzijl.R
#
# Version 4e reverts to using pressure as a forcing term on the RHS for
# comparison with earlier model results.
# We also include all the HadSLP2r data and don't remove duplicates due to the
# lower resolution (unlike 4c and 4d)
#
# Author: simonh
###############################################################################
##############
## Model 4e ##
##############
# As with model 4 but using POLCOMS run S12R408
# S12R408 was run for longer, so we need to truncate it at 480 months
#
# Assume that environment model4e already exists as it was created for Brest 
#model4e <- new.env() 
# Again, just use pressures over the periods of the model data (Pm)

# In Thompson's model 4 he also introduces lags of 0, 1 & 2 months which in this
# case gives a total of 51 components. However, to facilitate comparison with
# model 3 he only uses the first 9 components which he chooses by stepwise
# regression. We achieve this by adding in each component in turn and calculating
# the variance reduction, which we can then rank and analyze.
load("~/diskx/polcoms/S12run408/brestNewlyn/delfzijlMonthlyS12R408.RData", envir=S12R408)
tg$dm.mp.delfzijl <- tg$mp.delfzijl-mean(tg$mp.delfzijl, na.rm=T)

x11()
plot(tg$mp.time, tg$dm.mp.brest, type='l', col='blue', ylim=c(-200,300))
lines(tg$mp.time, tg$dm.mp.newlyn, col='red')
lines(tg$mp.time, tg$dm.mp.delfzijl, col='magenta')
x11()
plot(tg$dm.mp.delfzijl, tg$dm.mp.brest, col='blue')
points(tg$dm.mp.delfzijl, tg$dm.mp.newlyn, col='red')
# There is no correlation between Delfzijl and Brest or Newlyn at monthly timescales.
# Try 5 year running mean as per Bruce. 5 years is 60 months = 61 so odd numbered
trian61ptFilter<-c(1:31,seq(from=30,to=1,by=-1))
trian61ptFilter<-trian61ptFilter/sum(trian61ptFilter)
tg$filt.dm.mp.delfzijl<-filter(tg$dm.mp.delfzijl,trian61ptFilter,method="c",sides=2)
tg$filt.dm.mp.brest<-filter(tg$dm.mp.brest,trian61ptFilter,method="c",sides=2)
tg$filt.dm.mp.newlyn<-filter(tg$dm.mp.newlyn,trian61ptFilter,method="c",sides=2)
x11()
plot(tg$mp.time, tg$filt.dm.mp.brest, type='l', col='blue', ylim=c(-200,300))
lines(tg$mp.time, tg$filt.dm.mp.newlyn, col='red')
lines(tg$mp.time, tg$filt.dm.mp.delfzijl, col='magenta')

#model4e$Pm <- array(NA,dim=c(480,51))
# Zero lag
#for(i in 1:17){
#	model4e$Pm[,i] <- as.vector(hadley$mp.slpDelfzijl[,i])
#}
# Lag 1
#for(i in 18:34){
#	model4e$Pm[2:480,i] <- as.vector(hadley$mp.slpDelfzijl[,(i-17)])[1:479]
#	model4e$Pm[1,i] <- mean(as.vector(hadley$mp.slpDelfzijl[,(i-17)])[1:479], na.rm=T)
#}
# Lag 2
#for(i in 35:51){
#	model4e$Pm[3:480,i] <- as.vector(hadley$mp.slpDelfzijl[,(i-34)])[1:478]
#	model4e$Pm[1:2,i] <- mean(as.vector(hadley$mp.slpDelfzijl[,(i-34)])[1:478], na.rm=T)
#}

# What we are going to do now is calculate the sd of the residuals as each of the
# 51 components is added in turn 
model4e$m.delfzijl.sd <- vector(mode="numeric", length=51)
model4e$tg.delfzijl.sd <- vector(mode="numeric", length=51)

for(i in 1:51){
  model4e$data.delfzijl <- data.frame(mmsl=S12R408$delfzijlMonthlyMean[1:480,1], msl=tg$mp.delfzijl,
    t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2, Pm=model4e$Pm[,i])
	
  model4e$m.lm.delfzijl <- lm(mmsl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) +
    I(sin(w2*t)) + Pm, data=model4e$data.delfzijl)
  model4e$tg.lm.delfzijl <- lm(msl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + 
    I(sin(w2*t)) + Pm, data=model4e$data.delfzijl)
	
  model4e$m.resid.delfzijl <- model4e$m.lm.delfzijl$resid
  model4e$tg.resid.delfzijl <- model4e$tg.lm.delfzijl$resid
	
  model4e$m.delfzijl.sd[i] <- sd(model4e$m.resid.delfzijl, na.rm=T)
  model4e$tg.delfzijl.sd[i] <- sd(model4e$tg.resid.delfzijl, na.rm=T)
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model4e$m.delfzijl.sd.sorted <- sort(model4e$m.delfzijl.sd, index=TRUE)
model4e$tg.delfzijl.sd.sorted <- sort(model4e$tg.delfzijl.sd, index=TRUE)

# Now construct models of the top 9 components
model4e$data.delfzijl <- data.frame(mmsl=S12R408$delfzijlMonthlyMean[1:480,1], msl=tg$mp.delfzijl,
  t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2, Pm=model2$Pm,
  Pm1=model4e$Pm[,model4e$m.delfzijl.sd.sorted$ix[1]], 
  Pm2=model4e$Pm[,model4e$m.delfzijl.sd.sorted$ix[2]], Pm3=model4e$Pm[,model4e$m.delfzijl.sd.sorted$ix[3]],
  Pm4=model4e$Pm[,model4e$m.delfzijl.sd.sorted$ix[4]], Pm5=model4e$Pm[,model4e$m.delfzijl.sd.sorted$ix[5]],
  Pm6=model4e$Pm[,model4e$m.delfzijl.sd.sorted$ix[6]], Pm7=model4e$Pm[,model4e$m.delfzijl.sd.sorted$ix[7]],
  Pm8=model4e$Pm[,model4e$m.delfzijl.sd.sorted$ix[8]], Pm9=model4e$Pm[,model4e$m.delfzijl.sd.sorted$ix[9]])

model4e$m.lm.delfzijl <- lm(mmsl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)) + 
  Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 + Pm,  data=model4e$data.delfzijl)

model4e$data.delfzijl <- data.frame(mmsl=S12R408$delfzijlMonthlyMean[1:480,1], msl=tg$mp.delfzijl,
  t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2, Pm=model2$Pm,
  Pm1=model4e$Pm[,model4e$tg.delfzijl.sd.sorted$ix[1]], 
  Pm2=model4e$Pm[,model4e$tg.delfzijl.sd.sorted$ix[2]], Pm3=model4e$Pm[,model4e$tg.delfzijl.sd.sorted$ix[3]],
  Pm4=model4e$Pm[,model4e$tg.delfzijl.sd.sorted$ix[4]], Pm5=model4e$Pm[,model4e$tg.delfzijl.sd.sorted$ix[5]],
  Pm6=model4e$Pm[,model4e$tg.delfzijl.sd.sorted$ix[6]], Pm7=model4e$Pm[,model4e$tg.delfzijl.sd.sorted$ix[7]],
  Pm8=model4e$Pm[,model4e$tg.delfzijl.sd.sorted$ix[8]], Pm9=model4e$Pm[,model4e$tg.delfzijl.sd.sorted$ix[9]])

model4e$tg.lm.delfzijl <- lm(msl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)) + 
  Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 + Pm,  data=model4e$data.delfzijl)

model4e$m.resid.delfzijl <- model4e$m.lm.delfzijl$resid
model4e$tg.resid.delfzijl <- model4e$tg.lm.delfzijl$resid

#model4e$m.resid.delfzijl <- model4e$data.delfzijl$mmsl - model4e$m.lm.delfzijl$coef[2]*model4e$data.delfzijl$t + model4e$m.lm.delfzijl$coef[1]
#model4e$tg.resid.delfzijl <- model4e$data.delfzijl$msl - model4e$tg.lm.delfzijl$coef[2]*model4e$data.delfzijl$t + model4e$tg.lm.delfzijl$coef[1]

#> sd(model4e$data.delfzijl$mmsl, na.rm=T)
#[1] 0.07098115 Delfzijl: [1] 0.0676091 (model4)
#[1] 0.07975017 Delfzijl: [1] 0.07649279 (model4e)

#> sd(model4e$m.resid.delfzijl, na.rm=T)
#[1] 0.02322259 Delfzijl: [1] 0.02422649 (model4)
#[1] 0.02532422 Delfzijl: [1] 0.02690714 (model4e)

#> sd(model4e$data.delfzijl$msl, na.rm=T)
#[1] 86.98253 Delfzijl: 142.6104
#> sd(model4e$tg.resid.delfzijl, na.rm=T)
#[1] 36.48843 Delfzijl: 72.62975
# This is slightly worse than using the CS3X forcing data or the POLCOMS
# forcing data (suggesting that in fact this system is sensitive to
# resolution). However, it is encourgaing that this lower resolution set is
# good enough that we might be able to extend this approach to the full
# co-period of both Delfzijl and Brest. 

