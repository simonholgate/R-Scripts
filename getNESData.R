# This gets the North European Shelf data from the database
library(ROracle)
user <- ""
pass <- ""
db <- ""
drv <- dbDriver("Oracle")
con <- dbConnect(drv,paste(user, pass "@", db, sep=""))
selectString <- paste("SELECT ccode,scode,sname FROM sttn",
" WHERE ccode>='140' AND ccode<='180' and iyrlr IS NOT NULL",
" ORDER BY ccode,scode",sep="")
res <- dbGetQuery(con,selectString)
ccode <-res[,1]
scode <-res[,2]
sname <-res[,3]

numStations <- length(ccode)
NESMonthly <- array(dim=c(91*numStations,12))
annualYrs <- array(data=1914:2004,dim=c(91,1))
months <- array(data=
c('JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'),
dim=c(12,1))
for (j in 1:numStations){
  for (i in 1:12){
	selectString <- paste("SELECT yr,",months[i],",rlrfac FROM data",
	" WHERE ccode='",ccode[j],"' AND scode='",scode[j],"'",
        " ORDER BY yr",sep="")
	res <- dbGetQuery(con,selectString)
	indices <- which(is.element(res[,1],annualYrs))
	dataIndices <- which(is.element(annualYrs,res[indices,1]))
        dataIndices <- dataIndices + (j-1)*91
	NESMonthly[dataIndices,i] <- res[indices,2]+res[indices,3]
  }
}
dbDisconnect(con)
dbUnloadDriver(drv)

annualCycleNES<-mean(as.data.frame(NESMonthly),na.rm=TRUE)

plot((annualCycleNES-mean(annualCycleNES))/1000,type="b",ann=FALSE,col="blue")
lines((annualCycleNewlyn-mean(annualCycleNewlyn))/1000,type="b",ann=FALSE,col="red")
title(main="Mean Annual Sea Level Cycle 1914-2004 from Tide Gauges",xlab="Month number",ylab="Sea level height [m]")
grid(col="black")
legend(2,0.07,c("North European Shelf","Newlyn"),col=c("blue","red"),lty=c(1,1),bg="gray90")
