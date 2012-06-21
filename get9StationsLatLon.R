library(ROracle)
user <- ""
pass <- ""
db <- ""
drv <- dbDriver("Oracle")
con <- dbConnect(drv,paste(user, pass "@", db, sep=""))
nineStationsLatLon <- array(NA,dim=c(9,2))
for (i in 1:9){
	selectString <- paste("SELECT slatdec,slondec FROM simon_sttn_copy",
	" WHERE ccode='",nineStations[i,1],"' AND scode='",
	nineStations[i,2],"'",sep="")
	res <- dbGetQuery(con,selectString)
	nineStationsLatLon[i,1] <- res[[1]][1]
	nineStationsLatLon[i,2] <- res[[2]][1]
}
dbDisconnect(con)
dbUnloadDriver(drv)
