# Setup new oracle (9.2.0) first!
library(ROracle)
drv <- dbDriver("Oracle")
con <- dbConnect(drv,"psmsl2/psmsl2@BIA")
newlynAnnual <- vector(mode='numeric', length=90)
annualYrs <- 1915:2004
selectString <- paste("SELECT yr,ann,rlrfac FROM data",
  " WHERE ccode='170' AND scode='161' ORDER BY yr",sep="")
res <- dbGetQuery(con,selectString)
dbDisconnect(con)
dbUnloadDriver(drv)
for (i in 1:90){
  index <- which(annualYrs==res[[1]][i])
  newlynAnnual[index] <- res[[2]][index]+res[[3]][index]
}
