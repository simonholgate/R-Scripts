library(ROracle)
user <- ""
pass <- ""
db <- ""
drv <- dbDriver("Oracle")
con <- dbConnect(drv,paste(user, pass "@", db, sep=""))
lerwickMonthly <- array(dim=c(48,12))
annualYrs <- array(data=1957:2004,dim=c(48,1))
months <- array(data=
c('JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'),
dim=c(12,1))
for (i in 1:12){
	selectString <- paste("SELECT yr,",months[i],",rlrfac FROM data",
	" WHERE ccode='170' AND scode='001'",
        " ORDER BY yr",sep="")
	res <- dbGetQuery(con,selectString)
	indices <- which(is.element(res[,1],annualYrs))
	dataIndices <- which(is.element(annualYrs,res[indices,1]))
	lerwickMonthly[dataIndices,i] <- res[indices,2]+res[indices,3]
}
lerwickMonthly <- t(lerwickMonthly)
dim(lerwickMonthly) <- c((12*48),1)
# Get annual data
lerwickAnnual <- array(dim=c(48,1))
selectString <- paste("SELECT yr,ann,rlrfac FROM data",
" WHERE ccode='170' AND scode='001' ORDER BY yr",sep="")
res <- dbGetQuery(con,selectString)
indices <- which(is.element(res[,1],annualYrs))
dataIndices <- which(is.element(annualYrs,res[indices,1]))
lerwickAnnual[dataIndices] <- res[indices,2]+res[indices,3]

dbDisconnect(con)
dbUnloadDriver(drv)
