library(ROracle)
user <- ""
pass <- ""
db <- ""
drv <- dbDriver("Oracle")
con <- dbConnect(drv,paste(user, pass "@", db, sep=""))
newlynMonthly <- array(dim=c(93,12))
brestMonthly <- array(dim=c(93,12))
annualYrs <- array(data=1914:2006,dim=c(93,1))
months <- array(data=
c('JAN','FEB','MAR','APR','MAY','JUN','JUL','AUG','SEP','OCT','NOV','DEC'),
dim=c(12,1))
for (i in 1:12){
# Newlyn
	selectString <- paste("SELECT yr,",months[i],",rlrfac FROM data",
	" WHERE ccode='170' AND scode='161'",
        " ORDER BY yr",sep="")
	res <- dbGetQuery(con,selectString)
	indices <- which(is.element(res[,1],annualYrs))
	dataIndices <- which(is.element(annualYrs,res[indices,1]))
	newlynMonthly[dataIndices,i] <- res[indices,2]+res[indices,3]
# Brest
		selectString <- paste("SELECT yr,",months[i],",rlrfac FROM data",
	" WHERE ccode='190' AND scode='091'",
        " ORDER BY yr",sep="")
	res <- dbGetQuery(con,selectString)
	indices <- which(is.element(res[,1],annualYrs))
	dataIndices <- which(is.element(annualYrs,res[indices,1]))
	brestMonthly[dataIndices,i] <- res[indices,2]+res[indices,3]
	
}

newlynMonthly <- t(newlynMonthly)
dim(newlynMonthly) <- c((12*93),1)

brestMonthly <- t(brestMonthly)
dim(brestMonthly) <- c((12*93),1)

#junk <- lm(brestMonthly ~ c(1:1116))
#junk$coef[2]*12
#c(1:1116)
# 1.235732

#junk <- lm(newlynMonthly ~ c(1:1116))
#junk$coef[2]*12
#c(1:1116)
# 1.724267


dbDisconnect(con)
dbUnloadDriver(drv)
