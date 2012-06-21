# Plots for Brest-Newlyn paper
postscript(file="brestMinusNewlynModel.eps", colormodel="cmyk",
           width = 4.0, height = 3.0,
           horizontal = FALSE, onefile = FALSE, paper = "special",
           family = "ComputerModern", encoding = "TeXtext.enc")
plot(time, dm.mm.brestMonthlyMean, type='l', col='blue', ylim=c(-200,300), ann=F, lwd=2)
lines(time, dm.mm.newlynMonthlyMean, col='red', lwd=2, lty='dotdash')
lines(time, dm.mm.brestMonthlyMean-dm.mm.newlynMonthlyMean, col='magenta', lwd=2)
title(xlab="Year", ylab="Sea level deviation [mm]")
dev.off()

#######

postscript(file="brestMinusNewlynTG.eps", colormodel="cmyk",
           width = 4.0, height = 3.0,
           horizontal = FALSE, onefile = FALSE, paper = "special",
           family = "ComputerModern", encoding = "TeXtext.enc")
plot(tg$mp.time, tg$dm.mp.brest, type='l', col='blue', ylim=c(-200,300), ann=F, lwd=2)
lines(tg$mp.time, tg$dm.mp.newlyn, col='red', lwd=2, lty='dotdash')
lines(tg$mp.time, tg$dm.mp.brest-tg$dm.mp.newlyn, col='magenta', lwd=2)
title(xlab="Year", ylab="Sea level deviation [mm]")
dev.off()
#> sd(tg$dm.mp.brest-tg$dm.mp.newlyn)
#[1] 30.67184
#> sd(tg$dm.mp.brest)
#[1] 86.98253
#> sd(tg$dm.mp.newlyn)
#[1] 82.01923

#######

postscript(file="brestMinusNewlynTGFullPeriod.eps", colormodel="cmyk",
           width = 4.0, height = 3.0,
           horizontal = FALSE, onefile = FALSE, paper = "special",
           family = "ComputerModern", encoding = "TeXtext.enc")
png(file="brestMinusNewlynTGFullPeriod.png")
plot(tg$time, tg$dm.brest, type='l', col='blue', ylim=c(-330,330), ann=F, lwd=2)
lines(tg$time, tg$dm.newlyn, col='red', lwd=2, lty='dotdash')
lines(tg$time, tg$dm.resid, col='magenta', lwd=2)
lines(tg$time, tg$dm.mean, col='black', lwd=2)
legend(as.Date("1980/1/1"),-200, legend=c("Brest","Newlyn","Mean","Difference"), col=c("blue","red","black","magenta"),
       lty=c(1,4,1,1), lwd=2)

title(xlab="Year", ylab="Sea level deviation [mm]")
dev.off()
#> sd(tg$dm.brest-tg$dm.newlyn, na.rm=T)
#[1] 38.04088
#> sd(tg$dm.brest, na.rm=T)
#[1] 92.71704
#> sd(tg$dm.newlyn, na.rm=T)
#[1] 92.26965

#######

postscript(file="brestNewlynTGAnnualMeans.eps", colormodel="cmyk",
           width = 4.0, height = 3.0,
           horizontal = FALSE, onefile = FALSE, paper = "special",
           family = "ComputerModern", encoding = "TeXtext.enc")
png(file="brestNewlynTGAnnualMeans.png")
plot(tg$annualYrs, tg$annualMeans.brest, type='l', col='blue', ylim=c(-120,120), ann=F, lwd=2)
lines(tg$annualYrs, tg$annualMeans.newlyn, col='red', lwd=2, lty='dotdash')
lines(tg$annualYrs, tg$annualMeans.resid, col='magenta', lwd=2)
lines(tg$annualYrs, tg$annualMeans.mean, col='black', lwd=2)
legend(as.Date("1980/1/1"),-70, legend=c("Brest","Newlyn","Mean","Difference"), col=c("blue","red","black","magenta"), lty=c(1,4,1,1), lwd=2)

title(xlab="Year", ylab="Sea level deviation [mm]")
dev.off()
#> sd(tg$annualMeans.brest, na.rm=T)
#[1] 48.09981
#> sd(tg$annualMeans.newlyn, na.rm=T)
#[1] 51.10194

#######

postscript(file="brestNewlynTGDecadalTrends.eps", colormodel="cmyk",
           width = 4.0, height = 3.0,
           horizontal = FALSE, onefile = FALSE, paper = "special",
           family = "ComputerModern", encoding = "TeXtext.enc")
png(file="brestNewlynTGDecadalTrends.png")
plot(tg$midpointYrs, tg$tenYrTrends.brest, type='l', col='blue', ylim=c(-15,15), ann=F, lwd=2)
lines(tg$midpointYrs, tg$tenYrTrends.newlyn, col='red', lwd=2, lty='dotdash')
lines(tg$midpointYrs, tg$tenYrTrends.resid, col='magenta', lwd=2)
lines(tg$midpointYrs, tg$tenYrTrends.mean, col='black', lwd=2)

title(xlab="Year", ylab="Sea level deviation [mm]")
legend(as.Date("1925/1/1"),-8, legend=c("Brest","Newlyn","Mean","Difference"), col=c("blue","red","black","magenta"), lty=c(1,4,1,1), lwd=2)
dev.off()
#> sd(tg$annualMeans.brest, na.rm=T)
#[1] 48.09981
#> sd(tg$annualMeans.newlyn, na.rm=T)
#[1] 51.10194

#######

postscript(file="brestNewlynTGDecadalMeans.eps", colormodel="cmyk",
           width = 4.0, height = 3.0,
           horizontal = FALSE, onefile = FALSE, paper = "special",
           family = "ComputerModern", encoding = "TeXtext.enc")
png(file="brestNewlynTGDecadalMeans.png")
plot(tg$midpointYrs, tg$tenYrMeans.brest, type='l', col='blue', ylim=c(-100,100), ann=F, lwd=2)
lines(tg$midpointYrs, tg$tenYrMeans.newlyn, col='red', lwd=2, lty='dotdash')
lines(tg$midpointYrs, tg$tenYrMeans.resid, col='magenta', lwd=2)
lines(tg$midpointYrs, tg$tenYrMeans.mean, col='black', lwd=2)

title(xlab="Year", ylab="Sea level deviation [mm]")
legend(as.Date("1960/1/1"),-50, legend=c("Brest","Newlyn","Mean","Difference"), col=c("blue","red","black","magenta"), lty=c(1,4,1,1), lwd=2)
dev.off()
#> sd(tg$annualMeans.brest, na.rm=T)
#[1] 48.09981
#> sd(tg$annualMeans.newlyn, na.rm=T)
#[1] 51.10194
