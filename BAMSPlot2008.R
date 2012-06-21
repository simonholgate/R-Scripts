postscript("BAMSPlot2008.ps",horizontal = FALSE, onefile = FALSE,
paper = "special", width=7, height=7)

plot(X[1:48],yearMean[1:48],col="blue",lwd=2,
  type="n", lty=3, ann=FALSE, ylim=c(-4,8))
lines(X[1:48],yearMean[1:48],col="blue",lwd=2,
  type="l")
grid(col="black")

title(xlab="Decade mid-point [year]", ylab="Rate of sea level change [mm/yr]")
dev.off()


