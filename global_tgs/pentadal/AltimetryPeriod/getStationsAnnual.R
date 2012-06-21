#library(RMySQL)
library(ROracle)
#drv <- dbDriver("MySQL")
drv <- dbDriver("Oracle")
con <- dbConnect(drv, "psmsl2/psmsl2")
stnsAnnual <- array(NA,dim=c(lenYrs,numStns))
stnsAnnualYrs <- array(data=startYr:endYr,dim=c(lenYrs,1))
for (i in 1:numStns){
	selectString <- paste("SELECT yr,ann,rlrfac,missdays FROM data",
	" WHERE ccode='",stns$ccode[i],"' AND scode='",
	stns$scode[i],"' ORDER BY yr",sep="")
	res <- dbGetQuery(con,selectString)
	indices <- which(is.element(res[,1],stnsAnnualYrs))
# For comparison with jabber output, ignore years that have missing months
        junk<-which(substr(res[,4],25,26)=="XX")
        indices <- setdiff(indices,junk)
	dataIndices <- which(is.element(stnsAnnualYrs,res[indices,1]))
        if (stns[i,1] == "150"){
	  stnsAnnual[dataIndices,i] <- res[indices,2]
        } else {
	  stnsAnnual[dataIndices,i] <- res[indices,2]+res[indices,3]
        }
}
dbDisconnect(con)
dbUnloadDriver(drv)
