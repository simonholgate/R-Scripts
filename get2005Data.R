library(ROracle)
user <- ""
pass <- ""
db <- ""
drv <- dbDriver("Oracle")
conn <- dbConnect(drv,paste(user, pass "@", db, sep=""))

selectString <- paste("SELECT ccode,scode,datent FROM data",
" WHERE datent >= to_date('01012005','ddmmyyyy')",
" AND datent < to_date('01012006','ddmmyyyy')",
" ORDER BY ccode,scode",sep="")
res <- dbGetQuery(conn,selectString)

ccode <-res[,1]
numStations <- length(ccode)

junk<-array(NA,dim=c(numStations,3))
junk[,1] <-res[,1]
junk[,2] <-res[,2]
junk[,3] <-res[,3]

j <- which(junk[,1]=="A")
junk[j,1]<-"999"
cscode<-as.integer(junk[,1])+as.integer(junk[,2])/1000
cscodeUnique<-unique(cscode)

numStations <- length(cscodeUnique)
cscodeArray <- array(NA,dim=c(numStations,2))

#lonLatArray <- array(NA,dim=c(numStations,2))
fid<-file("tideGaugeLonLat2005.gmt","w")

for (j in 1:numStations){
  cscodeArray[j,1] <- substr(as.character(cscodeUnique[j]),1,3)
  if (cscodeArray[j,1]=="999") {
    cscodeArray[j,1] <- "A"
  }
  cscodeArray[j,2] <- substr(as.character(cscodeUnique[j]),5,7)

  selectString <- paste("SELECT slondec,slatdec FROM simon_sttn_copy",
  " WHERE ccode='",cscodeArray[j,1],"' AND scode='",cscodeArray[j,2],"'",sep="")
  res <- dbGetQuery(conn,selectString)
#  lonLatArray[j,1] <- res[,1]
#  lonLatArray[j,2] <- res[,2]
  writeLines(paste(res[,1],res[,2],sep="\t"),con=fid)
#  dbClearResult(res)
}

dbDisconnect(conn)
dbUnloadDriver(drv)
close(fid)

