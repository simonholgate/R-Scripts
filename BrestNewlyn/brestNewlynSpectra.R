ts.model.newlyn <- ts(dm.mm.newlynMonthlyMean, start=c(1960,1), frequency=12)
ts.model.brest <- ts(dm.mm.brestMonthlyMean, start=c(1960,1), frequency=12)
# TG timeseries has gaps making spectra difficult. Fill single gaps with
# spline fit
tg$spl.newlyn <- spline(c(1:1116),tg$dm.newlyn,xout=c(871,1059))
tg$spl.dm.newlyn <- tg$dm.newlyn
tg$spl.dm.newlyn[c(871,1059)] <- tg$spl.newlyn$y

tg$spl.brest <- spline(c(1:1116),tg$dm.brest,xout=c(793,821))
tg$spl.dm.brest <- tg$dm.brest
tg$spl.dm.brest[c(793,821)] <- tg$spl.brest$y

ts.tg.newlyn <- ts(tg$spl.dm.newlyn[466:1116],start=c(1953,1), frequency=12)
ts.tg.brest <- ts(tg$spl.dm.brest[466:1116],start=c(1953,1), frequency=12)
xspec.tg.brest.newlyn <- spec.pgram(ts.union(ts.tg.brest,ts.tg.newlyn), spans=c(3,3), main="Cross-spectrum of TG Brest & Newlyn")
xspec.model.brest.newlyn <- spec.pgram(ts.union(ts.model.brest,ts.model.newlyn), spans=c(3,3), main="Cross-spectrum of Model Brest & Newlyn")
x11()
plot.spec.coherency(xspec.tg.brest.newlyn)
