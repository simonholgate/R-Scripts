# Script to get monthly mean data for Leriwck, Wick, 
# Invergordon, Stornoway and Torshavn from the database.

library(ROracle)
user <- ""
pass <- ""
db <- ""
drv <- dbDriver("Oracle")
con <- dbConnect(drv,paste(user, pass "@", db, sep=""))
nYrs <- 51
nStns <- 6
monthlyData <- array(dim=c(nYrs,12,nStns))

stationCodes <- c('170','001','170','005','170','007','170','251','170','011', '015', '011')
dim(stationCodes)<-c(2,nStns)
annualYrs <- array(data=1957:(1957+nYrs-1),dim=c(nYrs,1))
months <- array(data=
c('JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'),
dim=c(12,1))
for (j in 1:nStns){
  for (i in 1:12){
	selectString <- paste("SELECT yr,",months[i],",rlrfac FROM data",
	" WHERE ccode='",stationCodes[1,j],"' AND scode='",
    stationCodes[2,j],"'",
        " ORDER BY yr",sep="")
	res <- dbGetQuery(con,selectString)
	indices <- which(is.element(res[,1],annualYrs))
	dataIndices <- which(is.element(annualYrs,res[indices,1]))
	monthlyData[dataIndices,i,j] <- res[indices,2]+res[indices,3]
  }
}

lerwickMonthly <- t(monthlyData[,,1])
dim(lerwickMonthly) <- c((12*nYrs),1)

wickMonthly <- t(monthlyData[,,2])
dim(wickMonthly) <- c((12*nYrs),1)

invergordonMonthly <- t(monthlyData[,,3])
dim(invergordonMonthly) <- c((12*nYrs),1)

stornowayMonthly <- t(monthlyData[,,4])
dim(stornowayMonthly) <- c((12*nYrs),1)

aberdeenMonthly <- t(monthlyData[,,5])
dim(aberdeenMonthly) <- c((12*nYrs),1)

torshavnMonthly <- t(monthlyData[,,6])
dim(torshavnMonthly) <- c((12*nYrs),1)

# Get annual data
annualData <- array(dim=c(nYrs,nStns))
for (i in 1:nStns) {
  selectString <- paste("SELECT yr,ann,rlrfac FROM data",
  " WHERE ccode='",stationCodes[1,i],"' AND scode='",
  stationCodes[2,i],"' ORDER BY yr",sep="")
  res <- dbGetQuery(con,selectString)
  indices <- which(is.element(res[,1],annualYrs))
  dataIndices <- which(is.element(annualYrs,res[indices,1]))
  annualData[dataIndices,i] <- res[indices,2]+res[indices,3]
}

lerwickAnnual <- annualData[,1]
wickAnnual <- annualData[,2]
invergordonAnnual <- annualData[,3]
stornowayAnnual <- annualData[,4]
aberdeenAnnual <- annualData[,5]
torshavnAnnual <- annualData[,6]

dbDisconnect(con)
dbUnloadDriver(drv)
