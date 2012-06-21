###############################################################################
# neuralNetwork.R
# 
# First attamepts at using neural networks to see if they can do a better job 
# than a linear model
#
# Author: simonh
###############################################################################

library(nnet)

time <- seq.Date(from=as.Date("1960/1/15"), to=as.Date("1999/12/15"), by="1 month")

corunaI<-new.env()
santanderI<-new.env()

load("~/diskx/polcoms/spain/analysis/corunaI.RData", envir=corunaI)
load("~/diskx/polcoms/spain/analysis/santanderI.RData", envir=santanderI)

corunaI$nmin <- apply(corunaI$data.corunaI, 2, min)
corunaI$nmax <- apply(corunaI$data.corunaI, 2, max)
corunaI$nrnge <- corunaI$nmax-corunaI$nmin

scaled.data.corunaI <- scale(corunaI$data.corunaI, center=corunaI$nmin, scale=corunaI$nrnge)

#

santanderI$nmin <- apply(santanderI$data.santanderI, 2, min)
santanderI$nmax <- apply(santanderI$data.santanderI, 2, max)
santanderI$nrnge <- santanderI$nmax-santanderI$nmin

scaled.data.santanderI <- scale(santanderI$data.santanderI, center=santanderI$nmin, scale=santanderI$nrnge)

#############
## Model 1 ##
#############

nn.corunaI<- nnet(mmsl ~ Pm + NEs + NNs,  
  data=scaled.data.corunaI, size=5, subset=c(1:382), decay = 5e-4, maxit = 200, rang = 0.1)

pred.corunaI <- predict(nn.corunaI, scaled.data.corunaI[383:480,], type='raw')

x11()
par(family="HersheySans")
plot(time[1:382],scaled.data.corunaI[1:382,1]-mean(scaled.data.corunaI[,1]), type='l')
lines(time[1:382],nn.corunaI$fitted.values-mean(nn.corunaI$fitted.values),col='red')
lines(time[1:382],nn.corunaI$residuals-mean(nn.corunaI$residuals),col='blue')

x11()
par(family="HersheySans")
plot(time[383:480],scaled.data.corunaI[383:480,1]-mean(scaled.data.corunaI[383:480,1]), type='l')
lines(time[383:480],pred.corunaI-mean(pred.corunaI),col='red')

#############
## Model 2 ##
#############

nn2.corunaI <- nnet(mmsl ~ t + w1cos + w1sin + w2cos + w2sin + Pm + NEs + NNs,  
  data=scaled.data.corunaI, size=5, subset=c(1:382), decay = 5e-4, maxit = 200, rang = 0.1)

pred2.corunaI <- predict(nn2.corunaI, scaled.data.corunaI[383:480,], type='raw')

x11()
par(family="HersheySans")
plot(time[1:382],scaled.data.corunaI[1:382,1]-mean(scaled.data.corunaI[,1]), type='l')
lines(time[1:382],nn2.corunaI$fitted.values-mean(nn2.corunaI$fitted.values),col='red')
lines(time[1:382],nn2.corunaI$residuals-mean(nn2.corunaI$residuals),col='blue')

x11()
par(family="HersheySans")
plot(time[383:480],scaled.data.corunaI[383:480,1]-mean(scaled.data.corunaI[383:480,1]), type='l')
lines(time[383:480],pred2.corunaI-mean(pred2.corunaI),col='red')

#############
## Model 3 ##
#############

nn3.corunaI <- nnet(mmsl ~ t + w1cos + w1sin + w2cos + w2sin +
  Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 +
  NEs + NNs + Pm,  data=scaled.data.corunaI, size=5, subset=c(1:382),
  decay = 5e-4, maxit = 1000, rang = 0.1)

pred3.corunaI <- predict(nn3.corunaI, scaled.data.corunaI[383:480,], type='raw')
resid3.corunaI <- scaled.data.corunaI[383:480,1]-mean(scaled.data.corunaI[383:480,1])-pred3.corunaI-mean(pred3.corunaI)

x11()
par(family="HersheySans")
plot(time[1:382],scaled.data.corunaI[1:382,1]-mean(scaled.data.corunaI[,1]), type='l', ann=F)
lines(time[1:382],nn3.corunaI$fitted.values-mean(nn3.corunaI$fitted.values),col='red')
lines(time[1:382],nn3.corunaI$residuals-mean(nn3.corunaI$residuals),col='blue')

x11()
par(family="HersheySans")
plot(time[383:480],scaled.data.corunaI[383:480,1]-mean(scaled.data.corunaI[383:480,1]), type='l', ann=F)
lines(time[383:480],pred3.corunaI-mean(pred3.corunaI),col='red')
lines(time[383:480],resid3.corunaI-mean(resid3.corunaI), col='magenta')

x11()
par(family="HersheySans")
plot((scaled.data.corunaI[383:480,1]-mean(scaled.data.corunaI[,1])),resid3.corunaI-mean(resid3.corunaI), pch=19)
corr.nn3.corunaI <- cor.test((scaled.data.corunaI[383:480,1]-mean(scaled.data.corunaI[,1])),resid3.corunaI-mean(resid3.corunaI))

nn3.santanderI <- nnet(mmsl ~ t + w1cos + w1sin + w2cos + w2sin +
  Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 +
  NEs + NNs + Pm,  data=scaled.data.santanderI, size=5, subset=c(1:382),
  decay = 5e-4, maxit = 1000, rang = 0.1)

pred3.santanderI <- predict(nn3.santanderI, scaled.data.santanderI[383:480,], type='raw')
resid3.santanderI <- scaled.data.santanderI[383:480,1]-mean(scaled.data.santanderI[383:480,1])-pred3.santanderI-mean(pred3.santanderI)

x11()
par(family="HersheySans")
plot(time[1:382],scaled.data.santanderI[1:382,1]-mean(scaled.data.santanderI[,1]), type='l', ann=F)
lines(time[1:382],nn3.santanderI$fitted.values-mean(nn3.santanderI$fitted.values),col='red')
lines(time[1:382],nn3.santanderI$residuals-mean(nn3.santanderI$residuals),col='blue')
#dev2bitmap(file="santanderIDecadalRatesAtlSpainNeuralFit.png", res=150)

x11()
par(family="HersheySans")
plot(time[383:480],scaled.data.santanderI[383:480,1]-mean(scaled.data.santanderI[383:480,1]), 
  type='l', ann=F, lwd=2, col='blue')
lines(time[383:480],pred3.santanderI-mean(pred3.santanderI),
  col='red', lwd=2)
lines(time[383:480],resid3.santanderI-mean(resid3.santanderI), 
  col='magenta', lwd=2)
title(ylab='Rate [mm/yr]', xlab='Year')
legend(x=1975, y=-10,
  legend=c('Original','Fit','Residual'),
  col=c('blue', 'red', 'magenta'), lwd=2)

#dev2bitmap(file="santanderIDecadalRatesAtlSpainNeuralPred.png", res=150)

x11()
par(family="HersheySans")
plot((scaled.data.santanderI[383:480,1]-mean(scaled.data.santanderI[,1])),resid3.santanderI-mean(resid3.santanderI), pch=19)
corr.nn3.santanderI <- cor.test((scaled.data.santanderI[383:480,1]-mean(scaled.data.santanderI[,1])),resid3.santanderI-mean(resid3.santanderI))
#dev2bitmap(file="santanderIDecadalRatesAtlSpainNeuralBias.png", res=150)



#############
## Model 4 ##
#############

nn4 <- nnet(mmsl ~ t + w1cos + w1sin + w2cos + w2sin + Pm1 + Pm2 +
  Es1 + Es2 + Ns1 + Ns2 + 
  NEs + NNs + Pm,  data=scaled.data, size=5, subset=c(1:382),
  decay = 5e-4, maxit = 500, rang = 0.1) 

pred4 <- predict(nn4, scaled.data[383:480,], type='raw')

x11()
par(family="HersheySans")
plot(time[1:382],scaled.data[1:382,1]-mean(scaled.data[,1]), type='l')
lines(time[1:382],nn4$fitted.values-mean(nn4$fitted.values),col='red')
lines(time[1:382],nn4$residuals-mean(nn4$residuals),col='blue')

x11()
par(family="HersheySans")
plot(time[383:480],scaled.data[383:480,1]-mean(scaled.data[383:480,1]), type='l')
lines(time[383:480],pred4-mean(pred4),col='red')

#############
## Model 4a ##
#############
ids <- sample(seq(1,480),400)
nn4a <- nnet(mmsl[ids] ~ t[ids] + w1cos[ids] + w1sin[ids] + w2cos[ids] + 
  w2sin[ids] + Pm1[ids] + Pm2[ids] +
  Es1[ids] + Es2[ids] + Ns1[ids] + Ns2[ids] + 
  NEs[ids] + NNs[ids] + Pm[ids],  data=scaled.data, size=5,
  decay = 5e-4, maxit = 500, rang = 0.1) 

pred4a <- predict(nn4a, scaled.data[-ids,], type='raw')

x11()
par(family="HersheySans")
plot(time[ids],scaled.data[ids,1]-mean(scaled.data[,1]), type='l')
lines(time[ids],nn4a$fitted.values-mean(nn4a$fitted.values),col='red')
lines(time[ids],nn4a$residuals-mean(nn4a$residuals),col='blue')

x11()
par(family="HersheySans")
plot(time[-ids],scaled.data[-ids,1]-mean(scaled.data[-ids,1]), type='l')
lines(time[-ids],pred4a-mean(pred4a),col='red')

##############
## Model 4b ##
##############
samp <- c(sample(1:50,25), sample(51:100,25), sample(101:150,25), sample(151:200,25),
  sample(201:250,25), sample(251:300,25), sample(301:350,25), sample(351:400,25))

nn4b.corunaI <- nnet(mmsl[samp] ~ t[samp] + w1cos[samp] + w1sin[samp] + w2cos[samp] + 
  w2sin[samp] + Pm1[samp] + Pm2[samp] +
  Es1[samp] + Es2[samp] + Ns1[samp] + Ns2[samp] + 
  NEs[samp] + NNs[samp] + Pm[samp],  data=scaled.data.corunaI, size=5,
  decay = 5e-4, maxit = 500, rang = 0.1) 

pred4b.corunaI <- predict(nn4b.corunaI, scaled.data.corunaI[401:480,], type='raw')

x11()
par(family="HersheySans")
plot(time[samp],scaled.data.corunaI[samp,1]-mean(scaled.data.corunaI[,1]), pch=19)
points(time[samp],nn4b.corunaI$fitted.values-mean(nn4b.corunaI$fitted.values),col='red', pch=19)
points(time[samp],nn4b.corunaI$residuals-mean(nn4b.corunaI$residuals),col='blue', pch=19)

x11()
par(family="HersheySans")
plot(time[401:480],scaled.data.corunaI[401:480,1]-mean(scaled.data.corunaI[401:480,1]), pch=19)
points(time[401:480],pred4b.corunaI-mean(pred4b.corunaI, na.rm=T),col='red', pch=19)

x11()
par(family="HersheySans")
plot((scaled.data.corunaI[samp,1]-mean(scaled.data.corunaI[,1]))-(nn4b.corunaI$residuals-mean(nn4b.corunaI$residuals)),nn4b.corunaI$residuals-mean(nn4b.corunaI$residuals), pch=19)
corr.corunaI <- cor.test((scaled.data.corunaI[samp,1]-mean(scaled.data.corunaI[,1]))-(nn4b.corunaI$residuals-mean(nn4b.corunaI$residuals)),nn4b.corunaI$residuals-mean(nn4b.corunaI$residuals))

nn4b.santander1 <- nnet(mmsl[samp] ~ t[samp] + w1cos[samp] + w1sin[samp] + w2cos[samp] + 
  w2sin[samp] + Pm1[samp] + Pm2[samp] +
  Es1[samp] + Es2[samp] + Ns1[samp] + Ns2[samp] + 
  NEs[samp] + NNs[samp] + Pm[samp],  data=scaled.data.santander1, size=5,
  decay = 5e-4, maxit = 500, rang = 0.1) 

pred4b.santander1 <- predict(nn4b.santander1, scaled.data.santanderI[401:480,], type='raw')

x11()
par(family="HersheySans")
plot(time[samp],scaled.data.santanderI[samp,1]-mean(scaled.data.santanderI[,1]), pch=19)
points(time[samp],nn4b.santanderI$fitted.values-mean(nn4b.santanderI$fitted.values),col='red', pch=19)
points(time[samp],nn4b.santanderI$residuals-mean(nn4b.santanderI$residuals),col='blue', pch=19)

x11()
par(family="HersheySans")
plot(time[401:480],scaled.data.santanderI[401:480,1]-mean(scaled.data.santanderI[401:480,1]), pch=19)
points(time[401:480],pred4b.santanderI-mean(pred4b.santanderI, na.rm=T),col='red', pch=19)

x11()
par(family="HersheySans")
plot((scaled.data.santanderI[samp,1]-mean(scaled.data.santanderI[,1]))-(nn4b.santanderI$residuals-mean(nn4b.santanderI$residuals)),nn4b.santanderI$residuals-mean(nn4b.santanderI$residuals), pch=19)
corr.santanderI <- cor.test((scaled.data.santanderI[samp,1]-mean(scaled.data.santanderI[,1]))-(nn4b.santanderI$residuals-mean(nn4b.santanderI$residuals)),nn4b.santanderI$residuals-mean(nn4b.santanderI$residuals))

#############
## Model 5 ##
#############
#samp <- c(sample(1:50,25), sample(51:100,25), sample(101:150,25), sample(151:200,25),
#  sample(201:250,25), sample(251:300,25), sample(301:350,25), sample(351:400,25))
ids <- sample(seq(1,480),400)

nn5 <- nnet(mmsl[ids] ~ t[ids] + w1cos[ids] + w1sin[ids] + w2cos[ids] + w2sin[ids] +
   Pm1[ids] + Pm2[ids] + Pm3[ids] + Pm4[ids] + Pm5[ids] + Pm6[ids] + Pm7[ids] + Pm8[ids] + Pm9[ids] +
   NEs[ids] + NNs[ids] + Pm[ids],  data=scaled.data, size=10,
  decay = 5e-3, maxit = 500, rang = 0.5, abstol=1e-3)
#    Pm1[ids] + Pm2[ids] + Pm3[ids] + Pm4[ids] + Pm5[ids] + Pm6[ids] + Pm7[ids] + Pm8[ids] + Pm9[ids] +
#  Es1[ids] + Es2[ids] + Es3[ids] + Es4[ids] + Es5[ids] + Es6[ids] + Es7[ids] + Es8[ids] + Es9[ids] +
#  Ns1[ids] + Ns2[ids] + Ns3[ids] + Ns4[ids] + Ns5[ids] + Ns6[ids] + Ns7[ids] + Ns8[ids] + Ns9[ids] +

#  Pm1[ids] + Es1[ids] + Ns1[ids] +
pred5 <- predict(nn5, scaled.data[-ids,], type='raw')

x11()
par(family="HersheySans")
plot(time[ids],scaled.data[ids,1]-mean(scaled.data[,1]), pch=19)
points(time[ids],nn5$fitted.values-mean(nn5$fitted.values),col='red', pch=19)
points(time[ids],nn5$residuals-mean(nn5$residuals),col='blue', pch=19)

x11()
par(family="HersheySans")
plot(time[-ids],scaled.data[-ids,1]-mean(scaled.data[-ids,1]), pch=19)
points(time[-ids],pred5-mean(pred5, na.rm=T),col='red',pch=19)

x11()
par(family="HersheySans")
plot((scaled.data[ids,1]-mean(scaled.data[,1]))-(nn5$residuals-mean(nn5$residuals)),nn5$residuals-mean(nn5$residuals), pch=19)
corr <- cor.test((scaled.data[ids,1]-mean(scaled.data[,1]))-(nn5$residuals-mean(nn5$residuals)),nn5$residuals-mean(nn5$residuals))

#############
## Model 6 ##
#############
ids <- sample(seq(1,480),400)

nn6.corunaI <- nnet(mmsl[ids] ~ t[ids] + w1cos[ids] + w1sin[ids] + w2cos[ids] + w2sin[ids] +
   NEs[ids] + NNs[ids] + Pm[ids],  data=scaled.data.corunaI, size=10,
  decay = 5e-3, maxit = 500, rang = 0.5)
#    Pm1[ids] + Pm2[ids] + Pm3[ids] + Pm4[ids] + Pm5[ids] + Pm6[ids] + Pm7[ids] + Pm8[ids] + Pm9[ids] +
#  Es1[ids] + Es2[ids] + Es3[ids] + Es4[ids] + Es5[ids] + Es6[ids] + Es7[ids] + Es8[ids] + Es9[ids] +
#  Ns1[ids] + Ns2[ids] + Ns3[ids] + Ns4[ids] + Ns5[ids] + Ns6[ids] + Ns7[ids] + Ns8[ids] + Ns9[ids] +

pred6 <- predict(nn6, scaled.data[-ids,], type='raw')

x11()
par(family="HersheySans")
plot(time[ids],scaled.data[ids,1]-mean(scaled.data[,1]), pch=19)
points(time[ids],nn6$fitted.values-mean(nn6$fitted.values),col='red', pch=19)
points(time[ids],nn6$residuals-mean(nn6$residuals),col='blue', pch=19)

x11()
par(family="HersheySans")
plot(time[-ids],scaled.data[-ids,1]-mean(scaled.data[-ids,1]), pch=19)
points(time[-ids],pred6-mean(pred6, na.rm=T),col='red',pch=19)

x11()
par(family="HersheySans")
plot((scaled.data[ids,1]-mean(scaled.data[,1]))-(nn6$residuals-mean(nn6$residuals)),nn6$residuals-mean(nn6$residuals), pch=19)
corr <- cor.test((scaled.data[ids,1]-mean(scaled.data[,1]))-(nn6$residuals-mean(nn6$residuals)),nn6$residuals-mean(nn6$residuals))

#############
## Model 7 ##
#############

# At this point we need to understand what is contributing to the last 27% of variance. 
# As with the linear regression we need to 
