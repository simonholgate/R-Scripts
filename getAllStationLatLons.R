library(ROracle)
user <- ""
pass <- ""
db <- ""
drv <- dbDriver("Oracle")
conn <- dbConnect(drv,paste(user, pass "@", db, sep=""))

fid<-file("tideGaugeLonLatAll.gmt","w")
selectString <- paste("SELECT slondec,slatdec FROM simon_sttn_copy",sep="")
res <- dbGetQuery(conn,selectString)
writeLines(paste(res[,1],res[,2],sep="\t"),con=fid)

dbDisconnect(conn)
dbUnloadDriver(drv)
close(fid)

