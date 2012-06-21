#library(RMySQL)
library(ROracle)
#drv <- dbDriver("MySQL")
drv <- dbDriver("Oracle")
con <- dbConnect(drv, "psmsl2/psmsl2")
stnsMonthly <- array(NA,dim=c(12,lenYrs,numStns))
stnsAnnualYrs <- array(data=startYr:endYr,dim=c(lenYrs,1))
for (i in 1:numStns){
	selectString <- paste("SELECT yr,jan,feb,mar,apr,may,jun,jul,",
	"aug,sep,oct,nov,dec,rlrfac,missdays FROM data",
	" WHERE ccode='",stns$ccode[i],"' AND scode='",
	stns$scode[i],"' ORDER BY yr",sep="")
	res <- dbGetQuery(con,selectString)
	indices <- which(is.element(res[,1],stnsAnnualYrs))
	lenIndices <- length(indices)
	dataIndices <- which(is.element(stnsAnnualYrs,res[indices,1]))
        if (stns[i,1] == "150"){
          for(j in 1:lenIndices){
            for(k in 1:12){
	          stnsMonthly[k,dataIndices[j],i] <- res[indices[j],(k+1)]
	        }
	      }
        } else {
          for(j in 1:lenIndices){
            for(k in 1:12){
	          stnsMonthly[k,dataIndices[j],i] <- res[indices[j],(k+1)]+res[indices[j],14]
	        }
	      }
        }
}
stnsMonthlyArray <- stnsMonthly
dim(stnsMonthly) <- c(lenMons,numStns)
dbDisconnect(con)
dbUnloadDriver(drv)
