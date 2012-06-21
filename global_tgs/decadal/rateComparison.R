# Script to compare the rates of the first and second halves of the records
# and plot the difference as a function of latitude
library(RMySQL)
drv <- dbDriver("MySQL")
con <- dbConnect(drv)
slat <- vector(mode="character", length=177)
for (i in 1:177){
  selectString <- paste("select slat from sttn where ccode=",
    as.character(stns177[i,1]), " and scode=",
    as.character(stns177[i,2]), sep='')
  slat[i] <- dbGetQuery(con,selectString)
}
dbDisconnect(con)
dbUnloadDriver(drv)
#
for (i in 1:177){
  if (substr(slat[[i]],8,8)=="S"){
    if (stns177Lat[i] == abs(stns177Lat[i])){
      stns177Lat[i] <- -stns177Lat[i]
    }
  }
}
#
diffS177PgrCorrRates <- array(NA,dim=c(3,177))
for ( i in 1:177 ){
  diffS177PgrCorrRates[1,i] <- mean(s177PgrCorrectedRates[1:22,i], na.rm=TRUE)
  diffS177PgrCorrRates[2,i] <- mean(s177PgrCorrectedRates[23:44,i], na.rm=TRUE)
  diffS177PgrCorrRates[3,i] <- 
    diffS177PgrCorrRates[2,i] - diffS177PgrCorrRates[1,i]
}
# Zonal means
zones <- c(-90,-30,-20,-10,0,10,20,30,40,50,60,90)
zones2 <- c(-35,-25,-15,-5,5,15,25,35,45,55,65)
zonalMeans <- vector(mode="numeric",length=11)
for ( i in 1:11) {
  j <- intersect(which(stns177Lat< zones[i]), which(stns177Lat< zones[i+1]))
  zonalMeans[i] <- mean(diffS177PgrCorrRates[3,j], na.rm=TRUE)
}
#plot(stns177Lat,diffS177PgrCorrRates[3,], col='blue')
plot(zones2,zonalMeans, col='blue')
