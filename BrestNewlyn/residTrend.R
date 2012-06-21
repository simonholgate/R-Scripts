plot(1916:2008, newlyn.tot.p.full, type='l', col='blue')
lines(1916:2008, brest22.tot.p.full, type='l', col='red')
plot(1916:2008, newlyn.tot.p.full, type='l', col='blue', ylim=c(1e5,1.04e5))
lines(1916:2008, brest22.tot.p.full, type='l', col='red')
lines(1916:2008, tg.lmRob.newlyn.full$fitted, type='l', col='cyan')
lines(tg.lmRob.newlyn.full$x[,2], tg.lmRob.newlyn.full$fitted, type='l',
      col='cyan')
(var(newlyn.tot.p.full,na.rm=T)-var(tg.lmRob.newlyn.full$fitted,na.rm=T))/
  var(newlyn.tot.p.full,na.rm=T)*100
(var(newlyn.tot.p.full,na.rm=T)-var(tg.lmRob.newlyn.full$resid,na.rm=T))/
  var(newlyn.tot.p.full,na.rm=T)*100
dev.set(3)
x11()
plot(tg.lmRob.newlyn.full$x[,2], tg.lmRob.newlyn.full$resid, type='l',
     col='blue')
summary(tg.lmRob.newlyn.full)
newlyn.clean<-(tg.lmRob.newlyn.full$x[,2]*1.597182e+01 +
               tg.lmRob.newlyn.full$resid)
plot(tg.lmRob.newlyn.full$x[,2], newlyn.clean, type='l',
     col='blue')
newlyn.clean<-(tg.lmRob.newlyn.full$x[,2]*1.597182e+01 +
               tg.lmRob.newlyn.full$resid + 7.210212e+04)
plot(tg.lmRob.newlyn.full$x[,2], newlyn.clean, type='l', col='blue')
plot(1916:2008, newlyn.tot.p.full, type='l', col='blue', ylim=c(1e5,1.04e5))
lines(tg.lmRob.newlyn.full$x[,2], newlyn.clean, type='l', col='cyan')
lines(tg.lmRob.newlyn.full$x[,2],
      tg.lmRob.newlyn.full$x[,2]*1.597182e+01 + 7.210212e+04,
      type='l', col='magenta')
1916*1.597182e+01 
1916*1.597182e+01 +  7.210212e+04
lmRob(newlyn.tot.p.full~c(1916:2008)) 
x11()
plot(tg.lmRob.newlyn.full$x[,2], tg.lmRob.newlyn.full$resid,
     type='l', col='blue')
plot(tg.lmRob.newlyn.full$x[,2], tg.lmRob.newlyn.full$resid/(1025*9.81),
     type='l', col='blue')
x11()
hist(tg.lmRob.newlyn.full$resid/(1025*9.81))
hist(tg.lmRob.newlyn.full$resid/(1025*9.81), n=12)
x11()
plot(1916:2008, newlyn.tot.p.full-mean(newlyn.tot.p.full, na.rm=T),
     type='l', col='blue')
lines(tg.lmRob.newlyn.full$x[,2], newlyn.clean-mean(newlyn.clean, na.rm=T),
      type='l', col='cyan')

plot(1916:2008, newlyn.tot.p.full-mean(newlyn.tot.p.full, na.rm=T),
     type='l', col='blue')
lines(tg.lmRob.newlyn.full$x[,2], newlyn.clean-mean(newlyn.clean, na.rm=T),
      type='l', col='cyan')
lines(tg.lmRob.newlyn.full$x[,2], tg.lmRob.newlyn.full$fitted, col='magenta')
lines(tg.lmRob.newlyn.full$x[,2], tg.lmRob.newlyn.full$fitted-
      mean(tg.lmRob.newlyn.full$fitted, na.rm=T), col='magenta')
(var(newlyn.tot.p.full,na.rm=T)-var(tg.lmRob.newlyn.full$resid,na.rm=T))/
  var(newlyn.tot.p.full,na.rm=T)*100
(var(newlyn.tot.p.full,na.rm=T)-var(newlyn.clean,na.rm=T))/
  var(newlyn.tot.p.full,na.rm=T)*100
(var(newlyn.tot.p.full,na.rm=T)-var(tg.lmRob.newlyn.full$fitted,na.rm=T))/
  var(newlyn.tot.p.full,na.rm=T)*100
newlyn.detrended <- detrend(newlyn.tot.p.full)
var(newlyn.detrended, na.rm =T)
(var(newlyn.detrended, na.rm =T)-var(tg.lmRob.newlyn.full$resid,na.rm=T))/
  var(newlyn.detrended,na.rm=T)*100
#[1] 36.96341
newlyn.fitted.detrended <- detrend(tg.lmRob.newlyn.full$fitted)
(var(newlyn.detrended, na.rm =T)-var(newlyn.fitted.detrended,na.rm=T))/
  var(newlyn.tot.p.full,na.rm=T)*100
#[1] 34.66034
(var(newlyn.tot.p.full, na.rm =T)-var(tg.lmRob.newlyn.full$resid,na.rm=T))/
  var(newlyn.tot.p.full,na.rm=T)*100
#[1] 85.60281
