## Load pressure data
tmp <- new.env()
## HadSLP2r data is 1850-2009 = 160 years = 1920 months. 37 cols - Newlyn, Brest and 35 grid points
load("~/Dropbox/brestNewlynData/pressure/brestNewlynHadSLP2Paper2.RData", envir=tmp)
hadnstns <- 37
hadSLP2rAnnualPressureArray <- array(NA, dim=c(160,hadnstns))
for(i in 1:hadnstns){
  tmp$hadSLP2rMonthlyPressureArray <- tmp$slpNewlynStns[,i]
  dim(tmp$hadSLP2rMonthlyPressureArray) <- c(12,160)
  hadSLP2rAnnualPressureArray[,i] <- colMeans(tmp$hadSLP2rMonthlyPressureArray, na.rm=T)
}
hadSLP2rAnnualTime <- c(1850:2009)
rm(tmp)

## CRU data is 1873-2000 = 128 years = 1536 months. 
cru.brest<-read.table("~/Dropbox/brestNewlynData/pressure/cru.brest.1873.2000.txt", col.names=c("Year","SLP"))

save(file="met.RData", list=ls())
