library(ROracle)
user <- ""
pass <- ""
db <- ""
drv <- dbDriver("Oracle")
con <- dbConnect(drv,paste(user, pass "@", db, sep=""))
rs <- dbSendQuery(con,statement = paste(
 "SELECT sname,ccode,scode", 
 "FROM sttn WHERE ccode='170' and scode='214'"))
data <- fetch(rs, n = -1)   # extract all rows
dim(data)
dbDisconnect(con)
