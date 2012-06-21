# Script to read in the monthly mean arrays created by makeZetDailyMeans.r and
# make global monthly means from these
yearsArray<-array(data=
  c(1960,1961,1962,1963,1964,1965,1966,1967,1968,1969,
  1970,1971,1972,1973,1974,1975,1976,1977,1978,1979,
  1980,1981,1982,1983,1984,1985,1986,1987,1988,1989,
  1990,1991,1992,1993,1994,1995,1996,1997,1998,1999), dim=c(40,1))

numYears <- length(yearsArray)
globalZetMonthlyMean <- vector(mode="numeric", length=(numYears*12))

for (jj in 1:numYears) {
  inConn <- file(paste("zett.mm.",yearsArray[jj],".dat",sep=""), "rb") 
  zetMonthlyMeans <-readBin(inConn, what="numeric", n = 28420*12, 
    size = NA, endian = "big")
  close(inConn)

  message(paste("Processing year:",as.character(yearsArray[jj])))

  dim(monthlyMeans) <- c(28420,12)

  zetMMs <- vector(mode="numeric", length=12)
  for (mm in 1:12){
    zetMMs[mm] <- mean(zetMonthlyMean[,mm])
  }

  globalZetMonthlyMean[((jj*12)-11):(jj*12)] <- zetMMs

}
