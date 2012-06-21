# Script to read in the monthly mean arrays created by makeZetDailyMeans.r and
# extract the data at each station location

#######################
# _*S2 dimensions:    #
#                     #
# *_Xmin:   -19.83333 #
# Xres:   1/6         #
# Xn:     198         #
#                     #
# Ymin:   40.11111    #
# Yres:   1/9         #
# Yn:     224         #
#######################

# Newlyn 50 06 N  05 33 W => 50.1 -5.55
#> (-5.55-xmin)/xres
#[1] 85.69998
#> (50.1-ymin)/yres
#[1] 89.90001
# Newlyn is at x=86, y=90 
# intersect(which(isea==86),which(jsea==90)) => ip=17701

xmin <- -19.83333
xres <- 1/6
ymin <- 40.11111
yres <- 1/9

####################
# Read header info #
####################
inNpseaIseaJsea <- file("npseaiseajsea.dat", "rb")
  npsea<-readBin(inNpseaIseaJsea, what="integer", n = 1, size = NA)

  isea<-array(NA,dim=npsea)
  jsea<-array(NA,dim=npsea)

  isea<-readBin(inNpseaIseaJsea, what="integer", n = npsea, size = NA)
  jsea<-readBin(inNpseaIseaJsea, what="integer", n = npsea, size = NA)
close(inNpseaIseaJsea)

yearsArray<-array(data=
  c(1960,1961,1962,1963,1964,1965,1966,1967,1968,1969,
  1970,1971,1972,1973,1974,1975,1976,1977,1978,1979,
  1980,1981,1982,1983,1984,1985,1986,1987,1988,1989,
  1990,1991,1992,1993,1994,1995,1996,1997,1998,1999), dim=c(40,1))

numYears <- length(yearsArray)

########################
# Read in station list #
########################
stns <- read.table(file="~/diskx/modelComparison/177StationsLatLons.txt",
  col.names=c("lat","lon","code"))

# Calculate nearest grid point in model, based on station lat/lon
numStns <- length(stns$lat)
stnPoints <- vector(mode="integer", length=numStns)
iseai <- vector(mode="integer", length=numStns)
jseai <- vector(mode="integer", length=numStns)

for (i in 1:numStns){
# Convert latitudes to -180:+180
  if (stns$lon[i] > 180) {
    ilon <- stns$lon[i] - 360
  } else {
    ilon <- stns$lon[i]
  }
  iseai[i] <- round((ilon - xmin)/xres)
  jseai[i] <- round((stns$lat[i] - ymin)/yres)
  ipi <- intersect(which(isea==iseai[i]),which(jsea==jseai[i]))
  stnPoints[i] <- max(ipi, na.rm=TRUE)
}

numModelStns <- length(which(is.finite(stnPoints)))
modelStnCodes <- stns$code[which(is.finite(stnPoints))]
modelStnsLat <- stns$lat[which(is.finite(stnPoints))]
modelStnsLon <- stns$lon[which(is.finite(stnPoints))]
modelStnPoints <- stnPoints[which(is.finite(stnPoints))]

stationMonthlyMeans <- array(NA, dim=c(numModelStns,numYears*12))

for (jj in 1:numYears) {
  inConn <- file(paste("zett.mm.",yearsArray[jj],".dat",sep=""), "rb") 
  zetMonthlyMeans <-readBin(inConn, what="numeric", n = npsea*12, 
    size = 4)
  close(inConn)

  message(paste("Processing year:",as.character(yearsArray[jj])))

  dim(zetMonthlyMeans) <- c(npsea,12)

  zetMMs <- array(NA,dim=c(numModelStns,12))
  for (mm in 1:12){
    zetMMs[,mm] <- zetMonthlyMeans[modelStnPoints,mm]
  }

  stationMonthlyMeans[,((jj*12)-11):(jj*12)] <- zetMMs

}

save(file="stationMonthlyMeans.RData", stationMonthlyMeans)
