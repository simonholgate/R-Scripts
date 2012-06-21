# neuralNetwork.R
# 
# First attempts at using neural networks to see if they can do a better job 
# than a linear model
#
# Author: simonh
###############################################################################

library(nnet)

load("~/diskx/polcoms/brestNewlyn/analysis/model6.RData")

nmin <- apply(data.newlyn, 2, min)
nmax <- apply(data.newlyn, 2, max)
nrnge <- nmax-nmin

scaled.data <- scale(data.newlyn, center=nmin, scale=nrnge)
#scaled.data <- scale(data.newlyn)

#############
## Model 1 ##
#############

nn <- nnet(mmsl ~ Pm + NEs + NNs,  
  data=scaled.data, size=5, subset=c(1:382), decay = 5e-4, maxit = 200, rang = 0.1)

pred <- predict(nn, scaled.data[383:480,], type='raw')

x11()
par(family="HersheySans")
plot(time[1:382],scaled.data[1:382,1]-mean(scaled.data[,1]), type='l')
lines(time[1:382],nn$fitted.values-mean(nn$fitted.values),col='red')
lines(time[1:382],nn$residuals-mean(nn$residuals),col='blue')

x11()
par(family="HersheySans")
plot(time[383:480],scaled.data[383:480,1]-mean(scaled.data[383:480,1]), type='l')
lines(time[383:480],pred-mean(pred),col='red')

#############
## Model 2 ##
#############

nn2 <- nnet(mmsl ~ t + w1cos + w1sin + w2cos + w2sin + Pm + NEs + NNs,  
  data=scaled.data, size=2, subset=c(1:382), decay = 5e-4, maxit = 200, rang = 0.1)

pred2 <- predict(nn2, scaled.data[383:480,], type='raw')

x11()
par(family="HersheySans")
plot(time[1:382],scaled.data[1:382,1]-mean(scaled.data[,1]), type='l')
lines(time[1:382],nn2$fitted.values-mean(nn2$fitted.values),col='red')
lines(time[1:382],nn2$residuals-mean(nn2$residuals),col='blue')

x11()
par(family="HersheySans")
plot(time[383:480],scaled.data[383:480,1]-mean(scaled.data[383:480,1]), type='l')
lines(time[383:480],pred2-mean(pred2),col='red')

#############
## Model 3 ##
#############

#nn3 <- nnet(mmsl ~ t + w1cos + w1sin + w2cos + w2sin +
#  Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 +
#  Es1 + Es2 + Es3 + Es4 + Es5 + Es6 + Es7 + Es8 + Es9 +
#  Ns1 + Ns2 + Ns3 + Ns4 + Ns5 + Ns6 + Ns7 + Ns8 + Ns9 +
#  NEs + NNs + Pm,  data=scaled.data, size=5, subset=c(1:382),
#  decay = 5e-4, maxit = 500, rang = 0.1)
nn3 <- nnet(mmsl ~ t + w1cos + w1sin + w2cos + w2sin +
  Pm1 + Pm2 + Pm3 + Pm4 + Pm5 + Pm6 + Pm7 + Pm8 + Pm9 +
  NEs + NNs + Pm,  data=scaled.data, size=5, subset=c(1:382),
  decay = 5e-4, maxit = 500, rang = 0.1)
#nn3 <- nnet(mmsl ~ t + w1cos + w1sin + w2cos + w2sin +
#  Pm1 + Pm,  
#  data=scaled.data, size=5, subset=c(1:382),
#  decay = 5e-4, maxit = 200, rang = 0.1)  

pred3 <- predict(nn3, scaled.data[383:480,], type='raw')

x11()
par(family="HersheySans")
plot(time[1:382],scaled.data[1:382,1]-mean(scaled.data[,1]), type='l')
lines(time[1:382],nn3$fitted.values-mean(nn3$fitted.values),col='red')
lines(time[1:382],nn3$residuals-mean(nn3$residuals),col='blue')

x11()
par(family="HersheySans")
plot(time[383:480],scaled.data[383:480,1]-mean(scaled.data[383:480,1]), type='l')
lines(time[383:480],pred3-mean(pred3),col='red')

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

nn4b <- nnet(mmsl[samp] ~ t[samp] + w1cos[samp] + w1sin[samp] + w2cos[samp] + 
  w2sin[samp] + Pm1[samp] + Pm2[samp] +
  Es1[samp] + Es2[samp] + Ns1[samp] + Ns2[samp] + 
  NEs[samp] + NNs[samp] + Pm[samp],  data=scaled.data, size=5,
  decay = 5e-4, maxit = 500, rang = 0.1) 

pred4b <- predict(nn4b, scaled.data[401:480,], type='raw')

x11()
par(family="HersheySans")
plot(time[samp],scaled.data[samp,1]-mean(scaled.data[,1]), pch=19)
points(time[samp],nn4b$fitted.values-mean(nn4b$fitted.values),col='red', pch=19)
points(time[samp],nn4b$residuals-mean(nn4b$residuals),col='blue', pch=19)

x11()
par(family="HersheySans")
plot(time[401:480],scaled.data[401:480,1]-mean(scaled.data[401:480,1]), pch=19)
points(time[401:480],pred4b-mean(pred4b, na.rm=T),col='red', pch=19)

x11()
par(family="HersheySans")
plot((scaled.data[samp,1]-mean(scaled.data[,1]))-(nn4b$residuals-mean(nn4b$residuals)),nn4b$residuals-mean(nn4b$residuals), pch=19)
corr <- cor.test((scaled.data[samp,1]-mean(scaled.data[,1]))-(nn4b$residuals-mean(nn4b$residuals)),nn4b$residuals-mean(nn4b$residuals))

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

nn6 <- nnet(mmsl[ids] ~ t[ids] + w1cos[ids] + w1sin[ids] + w2cos[ids] + w2sin[ids] +
   NEs[ids] + NNs[ids] + Pm[ids],  data=scaled.data, size=10,
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