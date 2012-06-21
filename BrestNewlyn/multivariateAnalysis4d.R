# multivariateAnalysisBrest4d.R
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
##############
## Model 4d ##
##############
# As with model 4 but using POLCOMS run S12R408
# S12R408 was run for longer, so we need to truncate it at 480 months
# 
model4d <- new.env() 
# Again, just use pressures over the periods of the model data (Pm)

# In Thompson's model 4 he also introduces lags of 0, 1 & 2 months which in this
# case gives a total of 51 components. However, to facilitate comparison with
# model 3 he only uses the first 9 components which he chooses by stepwise
# regression. We achieve this by adding in each component in turn and calculating
# the variance reduction, which we can then rank and analyze.
model4d$Pm <- array(NA,dim=c(480,27))
# Zero lag
for(i in 1:9){
  model4d$Pm[,i] <- as.vector(hadley$mp.slp[,(i+7)])
}
# Lag 1
for(i in 10:18){
  model4d$Pm[2:480,i] <- as.vector(hadley$mp.slp[,(i-2)])[1:479]
  model4d$Pm[1,i] <- mean(as.vector(hadley$mp.slp[,(i-2)])[1:479], na.rm=T)
}
# Lag 2
for(i in 19:27){
  model4d$Pm[3:480,i] <- as.vector(hadley$mp.slp[,(i-11)])[1:478]
  model4d$Pm[1:2,i] <- mean(as.vector(hadley$mp.slp[,(i-11)])[1:478], na.rm=T)
}

# What we are going to do now is calculate the sd of the residuals as each of the
# 27 components is added in turn 
model4d$m.sd <- vector(mode="numeric", length=27)
model4d$tg.sd <- vector(mode="numeric", length=27)

for(i in 1:27){
  model4d$data.brest <- data.frame(mmsl=S12R408$brestMonthlyMean[1:480,1], msl=tg$mp.brest,
    t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2, Pm=model4d$Pm[,i])
	
  model4d$m.lm.brest <- lm(mmsl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) +
    I(sin(w2*t)) + Pm, data=model4d$data.brest)
  model4d$tg.lm.brest <- lm(msl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + 
    I(sin(w2*t)) + Pm, data=model4d$data.brest)
	
  model4d$m.resid.brest <- model4d$m.lm.brest$resid
  model4d$tg.resid.brest <- model4d$tg.lm.brest$resid
	
  model4d$m.sd[i] <- sd(model4d$m.resid.brest, na.rm=T)
  model4d$tg.sd[i] <- sd(model4d$tg.resid.brest, na.rm=T)
}

# Now we can sort these sds. Thompson took the first 9 and we will do the same for
# comparison. 
model4d$m.sd.sorted <- sort(model4d$m.sd, index=TRUE)
model4d$tg.sd.sorted <- sort(model4d$tg.sd, index=TRUE)

# Now construct models of the top 9 components
model4d$data.brest <- data.frame(mmsl=S12R408$brestMonthlyMean[1:480,1], msl=tg$mp.brest,
  t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2, Pm=model2$Pm,
  Pm1=model4d$Pm[,model4d$m.sd.sorted$ix[1]], 
  Pm2=model4d$Pm[,model4d$m.sd.sorted$ix[2]], Pm3=model4d$Pm[,model4d$m.sd.sorted$ix[3]],
  Pm4=model4d$Pm[,model4d$m.sd.sorted$ix[4]], Pm5=model4d$Pm[,model4d$m.sd.sorted$ix[5]],
  Pm6=model4d$Pm[,model4d$m.sd.sorted$ix[6]], Pm7=model4d$Pm[,model4d$m.sd.sorted$ix[7]],
  Pm8=model4d$Pm[,model4d$m.sd.sorted$ix[8]], Pm9=model4d$Pm[,model4d$m.sd.sorted$ix[9]])

model4d$m.lm.brest <- lm(mmsl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)) + 
  Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 + Pm,  data=model4d$data.brest)

model4d$data.brest <- data.frame(mmsl=S12R408$brestMonthlyMean[1:480,1], msl=tg$mp.brest,
  t=seq(from=0,to=((40*12)-1)), w1=w1, w2=w2, Pm=model2$Pm,
  Pm1=model4d$Pm[,model4d$tg.sd.sorted$ix[1]], 
  Pm2=model4d$Pm[,model4d$tg.sd.sorted$ix[2]], Pm3=model4d$Pm[,model4d$tg.sd.sorted$ix[3]],
  Pm4=model4d$Pm[,model4d$tg.sd.sorted$ix[4]], Pm5=model4d$Pm[,model4d$tg.sd.sorted$ix[5]],
  Pm6=model4d$Pm[,model4d$tg.sd.sorted$ix[6]], Pm7=model4d$Pm[,model4d$tg.sd.sorted$ix[7]],
  Pm8=model4d$Pm[,model4d$tg.sd.sorted$ix[8]], Pm9=model4d$Pm[,model4d$tg.sd.sorted$ix[9]])

model4d$tg.lm.brest <- lm(msl ~ t + I(cos(w1*t)) + I(sin(w1*t)) + I(cos(w2*t)) + I(sin(w2*t)) + 
  Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 + Pm,  data=model4d$data.brest)

model4d$m.resid.brest <- model4d$m.lm.brest$resid
model4d$tg.resid.brest <- model4d$tg.lm.brest$resid

#model4d$m.resid.brest <- model4d$data.brest$mmsl - model4d$m.lm.brest$coef[2]*model4d$data.brest$t + model4d$m.lm.brest$coef[1]
#model4d$tg.resid.brest <- model4d$data.brest$msl - model4d$tg.lm.brest$coef[2]*model4d$data.brest$t + model4d$tg.lm.brest$coef[1]

#> sd(model4d$data.brest$mmsl, na.rm=T)
#[1] 0.07098115 Newlyn: [1] 0.0676091 (model4)
#[1] 0.07975017 Newlyn: [1] 0.07649279 (model4d)

#> sd(model4d$m.resid.brest, na.rm=T)
#[1] 0.02322259 Newlyn: [1] 0.02422649 (model4)
#[1] 0.02352657 Newlyn: [1] 0.06616694 (model4d)

#> sd(model4d$data.brest$msl, na.rm=T)
#[1] 86.98253 Newlyn: 82.10498
#> sd(model4d$tg.resid.brest, na.rm=T)
#[1] 33.23898 Newlyn: 25.16201
# No signficant improvement

