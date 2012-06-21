# Script to read in Church and White updated data (to 2007) for comaprison with my curve
cw <- read.table("church_white_new_gmsl.lis", col.names=c("Year","GMSL","Error"),
                 colClasses=c("numeric","numeric","numeric"))

gmsl <- cumsum(s177GlobalMean[1:58])
cw1952 <- which(cw$Year==1952.5)
cw1992 <- which(cw$Year==1992.5)
gmsl1992 <- which(s177DecadeMidPoints==1992.5)

x11()
plot(s177DecadeMidPoints[1:58], gmsl, type="l",ann=FALSE, lwd=2, col='red', ylim=c(0,120), xlim=c(1951,2009))
lines(s177DecadeMidPoints[1:58], gmsl+s177GlobalRS[1:58,2], lwd=2, col='red', lty='dotted')
lines(s177DecadeMidPoints[1:58], gmsl-s177GlobalRS[1:58,2], lwd=2, col='red', lty='dotted')
lines(cw$Year, cw$GMSL-cw$GMSL[cw1952], col='blue', lwd=2)
lines(cw$Year, cw$GMSL-cw$GMSL[cw1952]+cw$Error, lty="dotted", col='blue', lwd=2)
lines(cw$Year, cw$GMSL-cw$GMSL[cw1952]-cw$Error, lty="dotted", col='blue', lwd=2)

grid(col="black", lwd=1)
title(main=
"Integral of the decadal rates of SLR from the 177 records & CW",
xlab="Year", ylab="Sea level [mm]")
legend(x=1955, y=117, 
		legend=c('C&W','177 Stns'),
		col=c('blue','red'), lwd=2)


library(robust)
df.gmsl <- data.frame(year=s177DecadeMidPoints[1:58], gmsl=gmsl)
lm.gmsl <- lmRob(gmsl ~ year, data=df.gmsl)

df.cw <- data.frame(year=cw$Year[cw1952:138], gmsl=cw$GMSL[cw1952:138])
lm.cw <- lmRob(gmsl ~ year, data=df.cw)

# Altimetry period
df.gmsl1992 <- data.frame(year=s177DecadeMidPoints[gmsl1992:58], gmsl=gmsl[gmsl1992:58])
lm.gmsl1992 <- lmRob(gmsl ~ year, data=df.gmsl1992)

df.cw1992 <- data.frame(year=cw$Year[cw1992:138], gmsl=cw$GMSL[cw1992:138])
lm.cw1992 <- lmRob(gmsl ~ year, data=df.cw1992)
