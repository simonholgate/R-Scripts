#library(RMySQL)
library(ROracle)
#drv <- dbDriver("MySQL")
drv <- dbDriver("Oracle")
con <- dbConnect(drv, "psmsl2/psmsl2")
stns177Annual <- array(NA,dim=c(lenYrs,177))
stns177AnnualYrs <- array(data=startYr:endYr,dim=c(lenYrs,1))
for (i in 1:177){
	selectString <- paste("SELECT yr,ann,rlrfac,missdays FROM data",
	" WHERE ccode='",stns177[i,1],"' AND scode='",
	stns177[i,2],"' ORDER BY yr",sep="")
	res <- dbGetQuery(con,selectString)
	indices <- which(is.element(res[,1],stns177AnnualYrs))
# For comparison with jabber output, ignore years that have missing months
        junk<-which(substr(res[,4],25,26)=="XX")
        indices <- setdiff(indices,junk)
	dataIndices <- which(is.element(stns177AnnualYrs,res[indices,1]))
        if (stns177[i,1] == "150"){
	  stns177Annual[dataIndices,i] <- res[indices,2]
        } else {
	  stns177Annual[dataIndices,i] <- res[indices,2]+res[indices,3]
        }
}
dbDisconnect(con)
dbUnloadDriver(drv)
