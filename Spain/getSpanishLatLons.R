library(ROracle)

# Load in list of RLR sites to get data for
stations <- read.table("spanish_RLR_station_codes.txt", 
  col.names=c("ccode","scode","sname"), 
  colClasses=c("character", "character", "character"), sep=",")
  
drv <- dbDriver("Oracle")
con <- dbConnect(drv,"psmsl2/psmsl2@BIA")

numStations <- length(stations$ccode)
latLons <- array(NA,dim=c(numStations,2))

for (i in 1:numStations){

  selectString <- paste("SELECT slatdec, slondec FROM simon_sttn_copy",
	" WHERE ccode='",stations$ccode[i],"' AND scode='",stations$scode[i],"'",sep="")
  res <- dbGetQuery(con,selectString)
  latLons[i,1] <- res[1,1]
  latLons[i,2] <- res[1,2]
}

dbDisconnect(con)
dbUnloadDriver(drv)

save(file="spanishLatLons.RData", stations, latLons)