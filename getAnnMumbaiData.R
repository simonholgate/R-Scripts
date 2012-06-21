# Setup new oracle (9.2.0) first!
library(ROracle)
user <- ""
pass <- ""
db <- ""
drv <- dbDriver("Oracle")
con <- dbConnect(drv,paste(user, pass "@", db, sep=""))
mumbaiAnnual <- vector(mode='numeric', length=90)
annualYrs <- 1878:1994
selectString <- paste("SELECT yr,ann,rlrfac FROM data",
  " WHERE ccode='500' AND scode='041' ORDER BY yr",sep="")
res <- dbGetQuery(con,selectString)
dbDisconnect(con)
dbUnloadDriver(drv)
for (i in 1:117){
  index <- which(annualYrs==res[[1]][i])
  mumbaiAnnual[index] <- res[[2]][index]+res[[3]][index]
}
