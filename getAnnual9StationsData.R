library(ROracle)
user <- ""
pass <- ""
db <- ""
drv <- dbDriver("Oracle")
con <- dbConnect(drv,paste(user, pass "@", db, sep=""))
nineStationsAnnual <- array(dim=c(87,9))
annualYrs <- array(data=1914:2000,dim=c(87,1))
for (i in 1:9){
	selectString <- paste("SELECT yr,ann,rlrfac FROM data",
	" WHERE ccode='",nineStations[i,1],"' AND scode='",
	nineStations[i,2],"' ORDER BY yr",sep="")
	res <- dbGetQuery(con,selectString)
	indices <- which(is.element(res[,1],annualYrs))
	dataIndices <- which(is.element(annualYrs,res[indices,1]))
	nineStationsAnnual[dataIndices,i] <- res[indices,2]+res[indices,3]
}
dbDisconnect(con)
dbUnloadDriver(drv)
