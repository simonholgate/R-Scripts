# multivariateAnalysisBrest4e.R
# 
# Implementation of Thompson's methods (Thompson, 1980; 1986) to the tide gauge
# time series at Brest
# These are the same methods implemented in multivariateAnalysisNewlyn.R
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
model4e <- new.env() 
# Again, just use pressures over the periods of the model data (Pm)

# In Thompson's model 4 he also introduces lags of 0, 1 & 2 months which in this
# case gives a total of 51 components. However, to facilitate comparison with
# model 3 he only uses the first 9 components which he chooses by stepwise
# regression. We achieve this by adding in each component in turn and calculating
# the variance reduction, which we can then rank and analyze.
model4e$Pm <- array(NA,dim=c(480,51))
# Zero lag
for(i in 1:17){
	model4e$Pm[,i] <- as.vector(hadley$mp.slpNewlyn[,i])
}
# Lag 1
for(i in 18:34){
	model4e$Pm[2:480,i] <- as.vector(hadley$mp.slpNewlyn[,(i-17)])[1:479]
	model4e$Pm[1,i] <- mean(as.vector(hadley$mp.slpNewlyn[,(i-17)])[1:479], na.rm=T)
}
# Lag 2
for(i in 35:51){
	model4e$Pm[3:480,i] <- as.vector(hadley$mp.slpNewlyn[,(i-34)])[1:478]
	model4e$Pm[1:2,i] <- mean(as.vector(hadley$mp.slpNewlyn[,(i-34)])[1:478], na.rm=T)
}

# What we are going to do now is calculate the sd of the residuals as each of the
# 51 components is added in turn 
model4e$m.sd <- vector(mode="numeric", length=51)
model4e$tg.sd <- vector(mode="numeric", length=51)

for(i in 1:51){
  model4e$data.brest <- data.frame(mmsl=S12R408$brestMonthlyMean[1:480,1], msl=tg$mp.brest,
    t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2, Pm=model4e$Pm[,i])
	
  model4e$m.lm.brest <- lm(mmsl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) +
    I(sin(w2*t)) + Pm, data=model4e$data.brest)
  model4e$tg.lm.brest <- lm(msl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + 
    I(sin(w2*t)) + Pm, data=model4e$data.brest)
	
  model4e$m.resid.brest <- model4e$m.lm.brest$resid
  model4e$tg.resid.brest <- model4e$tg.lm.brest$resid
	
  model4e$m.sd[i] <- sd(model4e$m.resid.brest, na.rm=T)
  model4e$tg.sd[i] <- sd(model4e$tg.resid.brest, na.rm=T)
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model4e$m.sd.sorted <- sort(model4e$m.sd, index=TRUE)
model4e$tg.sd.sorted <- sort(model4e$tg.sd, index=TRUE)

# Now construct models of the top 9 components
model4e$data.brest <- data.frame(mmsl=S12R408$brestMonthlyMean[1:480,1], msl=tg$mp.brest,
  t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2, Pm=model2$Pm,
  Pm1=model4e$Pm[,model4e$m.sd.sorted$ix[1]], 
  Pm2=model4e$Pm[,model4e$m.sd.sorted$ix[2]], Pm3=model4e$Pm[,model4e$m.sd.sorted$ix[3]],
  Pm4=model4e$Pm[,model4e$m.sd.sorted$ix[4]], Pm5=model4e$Pm[,model4e$m.sd.sorted$ix[5]],
  Pm6=model4e$Pm[,model4e$m.sd.sorted$ix[6]], Pm7=model4e$Pm[,model4e$m.sd.sorted$ix[7]],
  Pm8=model4e$Pm[,model4e$m.sd.sorted$ix[8]], Pm9=model4e$Pm[,model4e$m.sd.sorted$ix[9]])

model4e$m.lm.brest <- lm(mmsl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)) + 
  Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 + Pm,  data=model4e$data.brest)

model4e$data.brest <- data.frame(mmsl=S12R408$brestMonthlyMean[1:480,1], msl=tg$mp.brest,
  t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2, Pm=model2$Pm,
  Pm1=model4e$Pm[,model4e$tg.sd.sorted$ix[1]], 
  Pm2=model4e$Pm[,model4e$tg.sd.sorted$ix[2]], Pm3=model4e$Pm[,model4e$tg.sd.sorted$ix[3]],
  Pm4=model4e$Pm[,model4e$tg.sd.sorted$ix[4]], Pm5=model4e$Pm[,model4e$tg.sd.sorted$ix[5]],
  Pm6=model4e$Pm[,model4e$tg.sd.sorted$ix[6]], Pm7=model4e$Pm[,model4e$tg.sd.sorted$ix[7]],
  Pm8=model4e$Pm[,model4e$tg.sd.sorted$ix[8]], Pm9=model4e$Pm[,model4e$tg.sd.sorted$ix[9]])

model4e$tg.lm.brest <- lm(msl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)) + 
  Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 + Pm,  data=model4e$data.brest)

model4e$m.resid.brest <- model4e$m.lm.brest$resid
model4e$tg.resid.brest <- model4e$tg.lm.brest$resid

#model4e$m.resid.brest <- model4e$data.brest$mmsl - model4e$m.lm.brest$coef[2]*model4e$data.brest$t + model4e$m.lm.brest$coef[1]
#model4e$tg.resid.brest <- model4e$data.brest$msl - model4e$tg.lm.brest$coef[2]*model4e$data.brest$t + model4e$tg.lm.brest$coef[1]

#> sd(model4e$data.brest$mmsl, na.rm=T)
#[1] 0.07098115 Newlyn: [1] 0.0676091 (model4)
#[1] 0.07975017 Newlyn: [1] 0.07649279 (model4e)

#> sd(model4e$m.resid.brest, na.rm=T)
#[1] 0.02322259 Newlyn: [1] 0.02422649 (model4)
#[1] 0.02352657 Newlyn: [1] 0.06616694 (model4e)

#> sd(model4e$data.brest$msl, na.rm=T)
#[1] 86.98253 Newlyn: 82.10498
#> sd(model4e$tg.resid.brest, na.rm=T)
#[1] 36.48843 Newlyn: 25.16201
# This is slightly worse than using the CS3X forcing data or the POLCOMS
# forcing data (suggesting that in fact this system is sensitive to
# resolution). However, it is encourgaing that this lower resolution set is
# good enough that we might be able to extend this approach to the full
# co-period of both Newlyn and Brest. 

