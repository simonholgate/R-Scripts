# Script to plot the annual cycles a Newlyn and the North European Shelf along
# with the sea level re-created from hydrographic data
postscript(file="annualCycles.ps")
plot((annualCycleNES-mean(annualCycleNES))/1000,type="b",ann=FALSE,col="bl
ue")
history()
x11();plot(monthlyMeanSeaLevels$month,monthlyMeanSeaLevels$mean-mean(monthlyMean
SeaLevels$mean),type='b')
x11();plot(monthlyMeanSeaLevels$month,monthlyMeanSeaLevels$mean-mean(monthlyMean
SeaLevels$mean),type='b')
x11();plot(c(1:12),newlynMeanAnnualCycle,type="b")

