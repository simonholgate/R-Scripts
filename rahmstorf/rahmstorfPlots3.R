# Script to produce plots for comment on Rahmstorf (ScienceExpress, 2006) 

# Global temperature
postscript("globalTemp.ps",horizontal = FALSE,
  onefile = FALSE, paper = "special", width=7, height=7)
#postscript("figure1.ps",horizontal = FALSE,
#  onefile = FALSE, paper = "special", width=5, height=12)

#def.par <- par(no.readonly = TRUE) # save default, for resetting...
#nf <- layout(matrix(c(1,2,3),3,1), widths=lcm(c(10,10,10)), 
#  heights=lcm(c(10,10,10)), respect=TRUE)
#layout.show(nf)

plot(data.all$years,data.all$tanom, col='blue',type='l', lwd=2, ann=FALSE)
lines(data.all$years,filtAnnMeanTemp,lwd=2, col='red')
#title("GISS global surface temperature from Hansen et al", xlab="Year",
#  ylab="Temperature [C]")
title(xlab="Year", ylab=expression(paste("Temperature anomaly [",degree,"C]")))
legend("topleft",c("Annual mean global temperature",
  "Temperature filtered with 15 yr running mean"), 
  col=c("blue","red"), lty=1, lwd=2, bty='n')

dev.off()

# Global sea level
postscript("globalSeaLevel.ps",horizontal = FALSE,
  onefile = FALSE, paper = "special", width=7, height=7)

plot(data.all$years,data.all$sl, type='l', lwd=2, col='blue', ann=FALSE)
lines(data.all$years,filtCWAnnual, lwd=2, col='red')
#title("Annual mean sea level from Church and White", 
title(xlab="Year", ylab="Global mean sea level [mm]")
legend("topleft",c("Annual mean sea level",
  "Sea level filtered with 15 yr running mean"), 
  col=c("blue","red"), lty=1, lwd=2, bty='n')

dev.off()

# Look at patterns in the scatter
# C+W graph suggests that 1930 is about the break point in the 2 sticks of the
# curve. Lets look at the correlation above and below 1930
#postscript("globalTempSeaLevelRateByEpoch.ps",horizontal = FALSE, 
postscript("figure1.ps",horizontal = FALSE,
  onefile = FALSE, paper = "special", width=7, height=7)
#> which(cwYears==1930.5)
#[1] 61
plot(filtAnnMeanTemp[1:41],diffFiltCWAnnual[1:41], ann=FALSE, pch=25, 
  ylim=c(0,3.5), xlim=c(-0.5,0.5))
points(filtAnnMeanTemp[42:53],diffFiltCWAnnual[42:53], pch=19)
points(filtAnnMeanTemp[54:98],diffFiltCWAnnual[54:98], pch=23)
points(filtAnnMeanTemp[99:121],diffFiltCWAnnual[99:121], pch=22, bg='black')
#points(filtAnnMeanTemp5Yr,diffCWAnn5Yr, pch=19, bg='red')
# Add lines to show Rahmstorf fit and my fit to first half of data
lines(c(-0.6,0.6), 3.4*(c(-0.6,0.6)+0.5), lwd=2)
#lines(c(-0.6,0), meanSlope1st*c(-0.6,0)+mean(myintercepts1st, na.rm=TRUE), 
#  lwd=2, lty='dotdash')
title(
#main="Relationship between rate of sea level change and global temperature", 
  xlab=expression(paste("Temperature anomaly [",degree,"C]")), 
  ylab="Rate of sea level change [mm/yr]")
legend("bottomright", c("1880-1920", "1921-1932", "1933-1977", "1978-2000"),
  pch=c(25,19,23,22), pt.bg=c('white','black','white','black'), bty='n')

par <- def.par
dev.off()
#
postscript("figure2.ps",horizontal = FALSE, 
  onefile = FALSE, paper = "special", width=7, height=7)

junkx <- c(data.all$years[8:115], rev(data.all$years[8:115]))
junkFy <- c(cumsum(predF)-(sum(predF[1:27]) -
 mean(data.all$sl[8:61]))+ConfInt1st,
  rev(cumsum(predF)-(sum(predF[1:27]) -
 mean(data.all$sl[8:61]))-ConfInt1st))
junkSy <- c(cumsum(predS)-(sum(predS[1:81]) -
 mean(data.all$sl[62:115]))+ConfInt2nd,
  rev(cumsum(predS)-(sum(predS[1:81]) -
 mean(data.all$sl[62:115]))-ConfInt2nd))

plot(data.all$years[7:114],(cumsum(predR[7:114])-(sum(predR[7:61]) -
  mean(c(data.all$sl[59],data.all$sl[62]))) - 80)/10, lwd=2,
  type='l',  ann=FALSE, xlim=c(1880,2000), ylim=c(-19, 5))
#polygon(junkx,junkFy,col="grey80", border=NA)
#polygon(junkx,junkSy,col="grey60", border=NA)
#lines(data.all$years[7:114],cumsum(predR[7:114])-(sum(predR[7:61]) -
#  0), lwd=2)
points(data.all$years,(data.all$sl - 80)/10, pch=20)
lines(data.all$years[8:115],(cumsum(predF)-(sum(predF[1:27]) -
  mean(data.all$sl[8:61])) - 80)/10, lwd=2, 
  lty='dotdash', col='blue')
lines(data.all$years[8:115],(cumsum(predS)-(sum(predS[1:81]) -
  mean(data.all$sl[62:115])) - 80)/10, lwd=2, 
  lty='dashed', col='red')

#title(main="Predictive skill of the linear fit", 
title(xlab="Year", 
  ylab="Sea Level [cm]")
legend("topleft", 
c("Prediction based on 1880-1940", 
  "Prediction based on all data 1880-2000", 
  "Prediction based on 1941-2000"),
  lty=c("dotdash","solid","dashed"), col=c('blue','black','red'),
  lwd=2, bty='n')
legend(1878.5,2.8, "Sea level reconstruction (from [4])", pch=20, col="black",
  bty="n", adj=-0.03)

dev.off()

