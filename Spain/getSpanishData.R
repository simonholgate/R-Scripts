library(ROracle)

# Load in list of RLR sites to get data for
stations <- read.table("spanish_RLR_station_codes.txt", 
  col.names=c("ccode","scode","sname"), 
  colClasses=c("character", "character", "character"), sep=",")
  
drv <- dbDriver("Oracle")
con <- dbConnect(drv,"psmsl2/psmsl2@BIA")

numStations <- length(stations$ccode)
# Earliest data is from 1944
numYrs <- length(1944:2007)

annualYrs <- array(data=1944:2007,dim=c(numYrs,1))
monthlyYrs <- seq.Date(from=as.Date("1944/1/15"),to=as.Date("2007/12/15"),
  by="1 month")
spanishMonthlyArr <- array(dim=c(numStations,numYrs,12))
spanishMonthly <- array(dim=c((numYrs*12),numStations))
spanishAnnual <- array(dim=c(numYrs,numStations))

months <- array(data=
c('JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'),
dim=c(12,1))

for (i in 1:numStations){

# Monthly data
  for (j in 1:12){
	selectString <- paste("SELECT yr,",months[j],",rlrfac FROM data",
	" WHERE ccode='",stations$ccode[i],"' AND scode='",stations$scode[i],"'",
        " ORDER BY yr",sep="")
	res <- dbGetQuery(con,selectString)
	indices <- which(is.element(res[,1],annualYrs))
	dataIndices <- which(is.element(annualYrs,res[indices,1]))
	spanishMonthlyArr[i,dataIndices,j] <- res[indices,2]+res[indices,3]
#	dbClearResult(res)
  }
  junk <- t(spanishMonthlyArr[i,,])
  dim(junk) <- c((12*numYrs),1)
  spanishMonthly[,i] <- junk

# Annual data 
  selectString <- paste("SELECT yr, ann ,rlrfac FROM data",
	" WHERE ccode='",stations$ccode[i],"' AND scode='",stations$scode[i],"'",
        " ORDER BY yr",sep="")
  res <- dbGetQuery(con,selectString)
  indices <- which(is.element(res[,1],annualYrs))
  dataIndices <- which(is.element(annualYrs,res[indices,1]))
  spanishAnnual[dataIndices,i] <- res[indices,2]+res[indices,3]
#  dbClearResult(res)
  
}

dbDisconnect(con)
dbUnloadDriver(drv)

save(file="spanishData.RData", 
  monthlyYrs, spanishMonthly, stations, numStations, spanishAnnual, annualYrs)

x11()
par(family="HersheySans")
plot(monthlyYrs,spanishMonthly[,1], type='l', col='red')
lines(monthlyYrs,spanishMonthly[,27], type='l', col='blue')


x11()
par(family="HersheySans")
plot(annualYrs,spanishAnnual[,1], type='l', col='red')
lines(annualYrs,spanishAnnual[,27], type='l', col='blue')