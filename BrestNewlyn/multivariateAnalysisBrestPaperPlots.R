# Plots for paper
#quartz()
png("seaLevelChangeBrestNewlyn.png")
plot(tg_annual$newlyn$Year, tg_annual$newlyn$Height, type='l', col='blue', ylim=c(-120,210), ann=F, lwd=2)
lines(tg_annual$brest$Year, tg_annual$brest$Height, col='red')
#lines(tg_annual$delfzijl$Year, tg_annual$delfzijl$Height, col='maroon')
lines(tg_annual$newlyn.lmRob$x[,2], tg_annual$newlyn.lmRob$fitted, col='blue', lty="dashed", lwd=2)
#lines(tg_annual$newlynGap.lmRob$x[,2], tg_annual$newlynGap.lmRob$fitted, col='darkgreen', lty='dotted')
#lines(tg_annual$delfzijl1916.lmRob$x[,2], tg_annual$delfzijl1916.lmRob$fitted, col='darkblue')
#lines(tg_annual$newlyn.lm$x[,2], tg_annual$newlyn.lm$fitted, col='green', lty='dashed')
lines(tg_annual$newlyn1st.lmRob$x[,2], tg_annual$newlyn1st.lmRob$fitted, col='blue', lty="dashed", lwd=2)
lines(tg_annual$newlyn2nd.lmRob$x[,2], tg_annual$newlyn2nd.lmRob$fitted, col='blue', lty="dashed", lwd=2)
lines(tg_annual$brest1916.lmRob$x[,2], tg_annual$brest1916.lmRob$fitted, col='red', lty="dashed", lwd=2)
#lines(tg_annual$brest1916.lm$x[,2], tg_annual$brest1916.lm$fitted, col='red', lty='dashed')
lines(tg_annual$brest1st.lmRob$x[,2], tg_annual$brest1st.lmRob$fitted, col='red', lty="dashed", lwd=2)
lines(tg_annual$brest2nd.lmRob$x[,2], tg_annual$brest2nd.lmRob$fitted, col='red', lty="dashed", lwd=2)
#title(main="Sea level change at Brest and Newlyn", ylab="Sea Level [mm]", xlab="Year")
title(ylab="Sea Level [mm]", xlab="Year")
dev.off()

# Scatterplot of residuals
#x11()
#par(family="HersheySans")
png("scatterplotBrestNewlynMonthly19532008.png")
plot(tg_monthly$newlyn5308.lmRob$resid, tg_monthly$brest5308.lmRob$resid, type='p', col='blue', ann=F, ylim=c(-200,300), xlim=c(-200,300))
title(xlab="Residual Sea Level at Newlyn [mm]", ylab="Residual Sea Level at Brest [mm]")
lines(c(-200,300),c(-200,300), col="black", lwd=2)
dev.off()

png("scatterplotBrestNewlynAnnual19532008.png")
plot(model3$tg.lmRob.newlynTotalP$resid, model3$tg.lmRob.brestTotalP$resid[as.integer(rownames(model3$tg.lmRob.newlynTotalP$x))], type='p', col='blue', ann=F, ylim=c(-600,850), xlim=c(-600,850))
title(xlab="Residual Sea Level at Newlyn [mm]", ylab="Residual Sea Level at Brest [mm]")
lines(c(-600,850),c(-600,850), col="black", lwd=2)
dev.off()
# Plot of residuals after subtracting Brest from Newlyn
#x11()
#par(family="HersheySans")

png("scatterplotSubtractingBrestFromNewlyn.png")
plot(tg_monthly$mp.time, tg_monthly$mp.resid.diff, type='l', col='red')
#title("Plot of residuals after subtracting Brest from Newlyn")
title(xlab="Year", ylab="Residual [mm]")
dev.off()

# Model 1
plot(model1$time, model1$data.brest$msl, col='blue', ann=F)
points(model1$tg.lmRob.brest$x[,2], model1$tg.lmRob.brest$fitted, col='red')
lines(model1$time, model1$data.brest$msl, col='blue')
lines(model1$time, model1$tg.lmRob.brest$fitted, col='red')
#title(main="Brest", xlab="Year", ylab="Total Pressure [hPa]")
title(xlab="Year", ylab="Total Pressure [hPa]")

plot(model1$time, model1$data.newlyn$msl, col='blue', ann=F)
points(model1$tg.lmRob.newlyn$x[,2], model1$tg.lmRob.newlyn$fitted, col='red')
lines(model1$time, model1$data.newlyn$msl, col='blue')
lines(model1$tg.lmRob.newlyn$x[,2], model1$tg.lmRob.newlyn$fitted, col='red')
#title(main="Newlyn", xlab="Year", ylab="Total Pressure [hPa]")
title(xlab="Year", ylab="Total Pressure [hPa]")

# Model 2
png("model2-1960-1999.png")
split.screen(c(2,1))
screen(1)
plot(model2$time, model2$data.brest$msl, col='blue', ann=F)
points(model2$tg.lmRob.brest$x[,2], model2$tg.lmRob.brest$fitted, col='red')
lines(model2$time, model2$data.brest$msl, col='blue')
lines(model2$time, model2$tg.lmRob.brest$fitted, col='red')
lines(model1$tg.lmRob.brest$x[,2], model1$tg.lmRob.brest$fitted, col='magenta')
#title(main="Model 2: Brest (HadSLP2r) 1960-1999", xlab="Year", ylab="Total Pressure [hPa]")
title(xlab="Year", ylab="Total Pressure [hPa]")

screen(2)
plot(model2$time, model2$data.newlyn$msl, col='blue', ann=F)
points(model2$tg.lmRob.newlyn$x[,2], model2$tg.lmRob.newlyn$fitted, col='red')
lines(model2$time, model2$data.newlyn$msl, col='blue')
lines(model2$tg.lmRob.newlyn$x[,2], model2$tg.lmRob.newlyn$fitted, col='red')
lines(model1$tg.lmRob.newlyn$x[,2], model1$tg.lmRob.newlyn$fitted, col='magenta')
#title(main="Model 2: Newlyn (HadSLP2r) 1960-1999", xlab="Year", ylab="Total Pressure [hPa]")
title(xlab="Year", ylab="Total Pressure [hPa]")
close.screen(all = TRUE)
dev.off()

# Model 3
png("model3-1953-2008.png")
split.screen(c(2,1))
screen(1)
plot(model3$time, model3$data.brest$msl, col='blue', ann=F)
points(model3$tg.lmRob.brest$x[,2], model3$tg.lmRob.brest$fitted, col='red')
lines(model3$time, model3$data.brest$msl, col='blue')
lines(model3$time, model3$tg.lmRob.brest$fitted, col='red')
lines(model1$tg.lmRob.brest$x[,2], model1$tg.lmRob.brest$fitted, col='magenta')
#title(main="Model 3: Brest (HadSLP2r) 1953-2008", xlab="Year", ylab="Total Pressure [hPa]")
title(xlab="Year", ylab="Total Pressure [hPa]")

screen(2)
plot(model3$time, model3$data.newlyn$msl, col='blue', ann=F)
points(model3$tg.lmRob.newlyn$x[,2], model3$tg.lmRob.newlyn$fitted, col='red')
lines(model3$time, model3$data.newlyn$msl, col='blue')
lines(model3$tg.lmRob.newlyn$x[,2], model3$tg.lmRob.newlyn$fitted, col='red')
lines(model1$tg.lmRob.newlyn$x[,2], model1$tg.lmRob.newlyn$fitted, col='magenta')
#title(main="Model 3: Newlyn (HadSLP2r) 1953-1999", xlab="Year", ylab="Total Pressure [hPa]")
title(xlab="Year", ylab="Total Pressure [hPa]")
close.screen(all = TRUE)
dev.off()

# Model 4
png("model4-1916-1943.png")
split.screen(c(2,1))
screen(1)
plot(model4$time, model4$data.brest$msl, col='blue', ann=F)
points(model4$tg.lmRob.brest$x[,2], model4$tg.lmRob.brest$fitted, col='red')
lines(model4$time, model4$data.brest$msl, col='blue')
lines(model4$time, model4$tg.lmRob.brest$fitted, col='red')
#title(main="Model 4: Brest (HadSLP2r) 1916-1943", xlab="Year", ylab="Total Pressure [hPa]")
title(xlab="Year", ylab="Total Pressure [hPa]")

screen(2)
plot(model4$time, model4$data.newlyn$msl, col='blue', ann=F)
points(model4$tg.lmRob.newlyn$x[,2], model4$tg.lmRob.newlyn$fitted, col='red')
lines(model4$time, model4$data.newlyn$msl, col='blue')
lines(model4$tg.lmRob.newlyn$x[,2], model4$tg.lmRob.newlyn$fitted, col='red')
#title(main="Model 4: Newlyn (HadSLP2r) 1916-1943", xlab="Year", ylab="Total Pressure [hPa]")
title(xlab="Year", ylab="Total Pressure [hPa]")
close.screen(all = TRUE)
dev.off()

# Model 5
png("model5-1916-2008.png")
split.screen(c(2,1))
screen(1)
plot(model5$time, model5$data.brest$msl, lwd=2, ann=F)
points(model5$tg.lmRob.brest$x[,2], model5$tg.lmRob.brest$fitted, lty='dotdash')
lines(model5$time, model5$data.brest$msl, lwd=2)
lines(model5$tg.lmRob.brest$x[,2], model5$tg.lmRob.brest$fitted, lty='dotdash')
lines(model5$tg.lmRob.brest_22$x[,2], model5$tg.lmRob.brest_22$fitted, lty='dotted', col='red')
#title(main="Model 5: Brest (HadSLP2r) 1916-2008", xlab="Year", ylab="Total Pressure [hPa]")
title(xlab="Year", ylab="Total Pressure [hPa]")
screen(2)
plot(model5$time, model5$data.newlyn$msl, lwd=2, ann=F)
points(model5$tg.lmRob.newlyn$x[,2], model5$tg.lmRob.newlyn$fitted, lty='dotdash')
lines(model5$time, model5$data.newlyn$msl, lwd=2)
lines(model5$tg.lmRob.newlyn$x[,2], model5$tg.lmRob.newlyn$fitted, lty='dotdash')
#lines(model5$tg.lmRob.newlyn$x[,2], model5$tg.lmRob.newlyn$x[,2]*coef(model5$tg.lmRob.newlyn)[2]+coef(model5$tg.lmRob.newlyn.fitted.lm)[1] , col='red', lty='dashed')
#title(main="Model 5: Newlyn (HadSLP2r) 1916-1999", xlab="Year", ylab="Total Pressure [hPa]")
title(xlab="Year", ylab="Total Pressure [hPa]")
close.screen(all = TRUE)
dev.off()

# Map
png("mapOfPressureStations.png")
library(maps)
library(mapdata)
map("worldHires", xlim=c(-25,15), ylim=c(35,65), interior=F, fill=F, col="grey50", resolution=0)
map.axes()
points(expand.grid(c(-20, -15,-10, -5, 0 , 5, 10), c(40,45,50,55,60)), pch=19, col="grey20")
points( y=50.1, x=-5.55, pch='*', col="grey20", cex=3)
points( y=48.38, x=-4.5, pch='*', col="grey20", cex=3)
dev.off()

# Rates
library(gplots)
png("ratesAndErrorBars.png")
plotCI(c(2.08,2.34,1.67,1.82,1.44, 1.43),col='blue', uiw=c(0.72,0.45,0.60,0.26,0.14,0.14), pch=c(21,22,22,23,24,25), ann=F,
       ylim=c(1,3.3), xaxt='n')
plotCI(c(1.54,-1,1.76,1.8,1.71, 1.71),col='red', uiw=c(0.61,0,0.37,0.29,0.12,0.12), pch=c(21,22,22,23,24,25), add=T, ann=F)
text(1,3.1,"1916 - 1943", srt=90, cex=1.2, font=2)
text(2,3.1,"1888 - 1943", srt=90, cex=1.2, font=2)
text(3,2.7,"1953 - 1999", srt=90, cex=1.2, font=2)
text(4,2.5,"1953 - 2008", srt=90, cex=1.2, font=2)
text(5,2.2,"1916 - 1999", srt=90, cex=1.2, font=2)
text(6,2.2,"1916 - 2008", srt=90, cex=1.2, font=2)
#title(main="Rates of Change with Error Bars", ylab="Robust Standard Error [mm/yr]")
title(ylab="Robust Standard Error [mm/yr]")
dev.off()

png("comparativePlotsBrest1850And1888.png")
split.screen(c(2,1))
screen(1)
# Model 6
plot(model6$time, model6$data.brest$msl, lwd=2, ann=F)
points(model6$tg.lmRob.brest$x[,2], model6$tg.lmRob.brest$fitted, lty='dotdash')
lines(model6$time, model6$data.brest$msl, lwd=2)
lines(model6$tg.lmRob.brest$x[,2], model6$tg.lmRob.brest$fitted, lty='dotdash')
#title(main="Model 6: Brest (HadSLP2r) 1850-1943", xlab="Year", ylab="Total Pressure [hPa]")
title(xlab="Year", ylab="Total Pressure [hPa]")
screen(2)
# Model 7
plot(model7$time, model7$data.brest$msl, lwd=2, ann=F)
points(model7$tg.lmRob.brest$x[,2], model7$tg.lmRob.brest$fitted, lty='dotdash')
lines(model7$time, model7$data.brest$msl, lwd=2)
lines(model7$tg.lmRob.brest$x[,2], model7$tg.lmRob.brest$fitted, lty='dotdash')
#title(main="Model 7: Brest (HadSLP2r) 1888-1943", xlab="Year", ylab="Total Pressure [hPa]")
title(xlab="Year", ylab="Total Pressure [hPa]")
close.screen(all = TRUE)
dev.off()

# Acceleration?
png("brestAcceleration18502008.png")
plot(model8$tg.lmRob.brest$x[,2], model8$tg.lmRob.brest$resid, lwd=2, ann=F, type='b', col='red')
lines(model5$tg.lmRob.newlyn$x[,2], model5$tg.lmRob.newlyn$resid, ann=F, type='b', lty='dotdash', col='blue')
#title(main="Model 8: Brest Residuals (HadSLP2r) 1850-2008", xlab="Year", ylab="Total Pressure [hPa]")
title(xlab="Year", ylab="Total Pressure [hPa]")
dev.off()

# Douglas IB scaling?
png("brestIBscaling18802008.png")
plot(model9$pa.lmRob.brest$x[,2], model9$pa.scaled.brest, ann=F, type='l', lty='dashed', col='magenta')
lines(model9$tg.lmRob.brest$x[,2], model9$tg.norm.brest, lwd=2, ann=F, type='l', col='red')
lines(model9$pa.lmRob.brest$x[,2], model9$pa.norm.brest, ann=F, type='l', col='blue')
#title(main="Model 9: Brest Normalised and Scaled Residuals Compared with Pressure 1880-2008", xlab="Year", ylab="Normalised sea level")
title(xlab="Year", ylab="Normalised sea level")
dev.off()
