library(ROracle)
drv <- dbDriver("Oracle")
con <- dbConnect(drv,"psmsl2/psmsl2@BIA")
newlynMonthly <- array(dim=c(91,12))
annualYrs <- array(data=1914:2004,dim=c(91,1))
months <- array(data=
c('JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'),
dim=c(12,1))
for (i in 1:12){
	selectString <- paste("SELECT yr,",months[i],",rlrfac FROM data",
	" WHERE ccode='170' AND scode='161'",
        " ORDER BY yr",sep="")
	res <- dbGetQuery(con,selectString)
	indices <- which(is.element(res[,1],annualYrs))
	dataIndices <- which(is.element(annualYrs,res[indices,1]))
	newlynMonthly[dataIndices,i] <- res[indices,2]+res[indices,3]
}
dbDisconnect(con)
dbUnloadDriver(drv)
