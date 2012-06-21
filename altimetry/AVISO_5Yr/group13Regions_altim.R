# Script to group the stations into 13 regions and calculate their means and
# standard deviations etc.

s177CSCode<-as.integer(stntable$ccode)+as.integer(stntable$scode)/1000
lenRates <- length(trends[,1])-1
lenRatesP1 <- lenRates + 1
stns177Lon <- lon
stns177Lat <- lat

jCounter<-0
jArray<-0

#
# 1. Scandianvia
#
j <- intersect(which(s177CSCode >= 25), which(s177CSCode < 131))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

stns177LonNeg <- stns177Lon
stns177LonNeg[j]<- -1*(360-stns177Lon[j])

s177ScandinaviaS <- trends[lenRatesP1,j]
s177Scandinavia <- trends[1:lenRates,j]
s177ScandinaviaLat <- stns177Lat[j]
s177ScandinaviaLon <- stns177LonNeg[j]
s177ScandinaviaNames <- stntable$sname[j]

s177ScandinaviaMean <- array(NA,dim=c(lenRatesP1,1))
s177ScandinaviaSD <- array(NA,dim=c(lenRatesP1,1))
s177ScandinaviaSE <- array(NA,dim=c(lenRatesP1,1))
for ( i in 1:lenRates){
  s177ScandinaviaMean[i] <- mean(s177Scandinavia[i,], na.rm=TRUE)
  s177ScandinaviaSD[i] <- sd(s177Scandinavia[i,], na.rm=TRUE)
  s177ScandinaviaSE[i] <- s177ScandinaviaSD[i]/
    sqrt(length(which(is.finite(s177Scandinavia[i,]))))
}
s177ScandinaviaMean[lenRatesP1] <- mean(s177ScandinaviaS, na.rm=TRUE)
s177ScandinaviaSD[lenRatesP1] <- sd(s177ScandinaviaS, na.rm=TRUE)
s177ScandinaviaSE[lenRatesP1] <- s177ScandinaviaSD[lenRatesP1]/
  sqrt(length(which(is.finite(s177ScandinaviaS))))

#
# 2. N Europe
#
j <- intersect(which(s177CSCode >= 140), which(s177CSCode < 200))
j <- union(j, which(s177CSCode == 210))
j <- union(j, which(s177CSCode == 360))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177NEuropeS <- trends[lenRatesP1,j]
s177NEurope <- trends[1:lenRates,j]
s177NEuropeLat <- stns177Lat[j]
s177NEuropeLon <- stns177LonNeg[j]
s177NEuropeNames <- stntable$sname[j]

s177NEuropeMean <- array(NA,dim=c(lenRatesP1,1))
s177NEuropeSD <- array(NA,dim=c(lenRatesP1,1))
s177NEuropeSE <- array(NA,dim=c(lenRatesP1,1))
for ( i in 1:lenRates){
  s177NEuropeMean[i] <- mean(s177NEurope[i,], na.rm=TRUE)
  s177NEuropeSD[i] <- sd(s177NEurope[i,], na.rm=TRUE)
  s177NEuropeSE[i] <- s177NEuropeSD[i]/
    sqrt(length(which(is.finite(s177NEurope[i,]))))
}
s177NEuropeMean[lenRatesP1] <- mean(s177NEuropeS, na.rm=TRUE)
s177NEuropeSD[lenRatesP1] <- sd(s177NEuropeS, na.rm=TRUE)
s177NEuropeSE[lenRatesP1] <- s177NEuropeSD[lenRatesP1]/
  sqrt(length(which(is.finite(s177NEuropeS))))

#
# 3. E Atlantic
#
j <- intersect(which(s177CSCode >= 340), which(s177CSCode <= 447))
j <- union(j, intersect(which(s177CSCode >= 200), which(s177CSCode < 201)))
j <- union(j, intersect(which(s177CSCode >= 220), which(s177CSCode <= 220.011)))
j <- union(j, intersect(which(s177CSCode >= 430), which(s177CSCode <= 430.081)))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177EAtlanticS <- trends[lenRatesP1,j]
s177EAtlantic <- trends[1:lenRates,j]
s177EAtlanticLat <- stns177Lat[j]
s177EAtlanticLon <- stns177LonNeg[j]
s177EAtlanticNames <- stntable$sname[j]

s177EAtlanticMean <- array(NA,dim=c(lenRatesP1,1))
s177EAtlanticSD <- array(NA,dim=c(lenRatesP1,1))
s177EAtlanticSE <- array(NA,dim=c(lenRatesP1,1))
for ( i in 1:lenRates){
  s177EAtlanticMean[i] <- mean(s177EAtlantic[i,], na.rm=TRUE)
  s177EAtlanticSD[i] <- sd(s177EAtlantic[i,], na.rm=TRUE)
  s177EAtlanticSE[i] <- s177EAtlanticSD[i]/
    sqrt(length(which(is.finite(s177EAtlantic[i,]))))
}
s177EAtlanticMean[lenRatesP1] <- mean(s177EAtlanticS, na.rm=TRUE)
s177EAtlanticSD[lenRatesP1] <- sd(s177EAtlanticS, na.rm=TRUE)
s177EAtlanticSE[lenRatesP1] <- s177EAtlanticSD[lenRatesP1]/
  sqrt(length(which(is.finite(s177EAtlanticS))))

#
# 4. Mediterranean
#
j <- intersect(which(s177CSCode >= 225), which(s177CSCode <= 290))
j <- union(j, intersect(which(s177CSCode > 220.011), which(s177CSCode < 221)))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177MediterraneanS <- trends[lenRatesP1,j]
s177Mediterranean <- trends[1:lenRates,j]
s177MediterraneanLat <- stns177Lat[j]
s177MediterraneanLon <- stns177LonNeg[j]
s177MediterraneanNames <- stntable$sname[j]

s177MediterraneanMean <- array(NA,dim=c(lenRatesP1,1))
s177MediterraneanSD <- array(NA,dim=c(lenRatesP1,1))
s177MediterraneanSE <- array(NA,dim=c(lenRatesP1,1))
for ( i in 1:lenRates){
  s177MediterraneanMean[i] <- mean(s177Mediterranean[i,], na.rm=TRUE)
  s177MediterraneanSD[i] <- sd(s177Mediterranean[i,], na.rm=TRUE)
  s177MediterraneanSE[i] <- s177MediterraneanSD[i]/
    sqrt(length(which(is.finite(s177Mediterranean[i,]))))
}
s177MediterraneanMean[lenRatesP1] <- mean(s177MediterraneanS, na.rm=TRUE)
s177MediterraneanSD[lenRatesP1] <- sd(s177MediterraneanS, na.rm=TRUE)
s177MediterraneanSE[lenRatesP1] <- s177MediterraneanSD[lenRatesP1]/
  sqrt(length(which(is.finite(s177MediterraneanS))))

#
# 5. SE Asia
#
j <- intersect(which(s177CSCode >= 530), which(s177CSCode <= 648))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177SEAsiaS <- trends[lenRatesP1,j]
s177SEAsia <- trends[1:lenRates,j]
s177SEAsiaLat <- stns177Lat[j]
s177SEAsiaLon <- stns177LonNeg[j]
s177SEAsiaNames <- stntable$sname[j]

s177SEAsiaMean <- array(NA,dim=c(lenRatesP1,1))
s177SEAsiaSD <- array(NA,dim=c(lenRatesP1,1))
s177SEAsiaSE <- array(NA,dim=c(lenRatesP1,1))
for ( i in 1:lenRates){
  s177SEAsiaMean[i] <- mean(s177SEAsia[i,], na.rm=TRUE)
  s177SEAsiaSD[i] <- sd(s177SEAsia[i,], na.rm=TRUE)
  s177SEAsiaSE[i] <- s177SEAsiaSD[i]/
    sqrt(length(which(is.finite(s177SEAsia[i,]))))
}
s177SEAsiaMean[lenRatesP1] <- mean(s177SEAsiaS, na.rm=TRUE)
s177SEAsiaSD[lenRatesP1] <- sd(s177SEAsiaS, na.rm=TRUE)
s177SEAsiaSE[lenRatesP1] <- s177SEAsiaSD[lenRatesP1]/sqrt(length(which(is.finite(s177SEAsiaS))))

# 
# 6. Australasia
#
j <- intersect(which(s177CSCode >= 680), which(s177CSCode <= 700))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177AustralasiaS <- trends[lenRatesP1,j]
s177Australasia <- trends[1:lenRates,j]
s177AustralasiaLat <- stns177Lat[j]
s177AustralasiaLon <- stns177LonNeg[j]
s177AustralasiaNames <- stntable$sname[j]

s177AustralasiaMean <- array(NA,dim=c(lenRatesP1,1))
s177AustralasiaSD <- array(NA,dim=c(lenRatesP1,1))
s177AustralasiaSE <- array(NA,dim=c(lenRatesP1,1))
for ( i in 1:lenRates){
  s177AustralasiaMean[i] <- mean(s177Australasia[i,], na.rm=TRUE)
  s177AustralasiaSD[i] <- sd(s177Australasia[i,], na.rm=TRUE)
  s177AustralasiaSE[i] <- s177AustralasiaSD[i]/ 
    sqrt(length(which(is.finite(s177Australasia[i,]))))
}
s177AustralasiaMean[lenRatesP1] <- mean(s177AustralasiaS, na.rm=TRUE)
s177AustralasiaSD[lenRatesP1] <- sd(s177AustralasiaS, na.rm=TRUE)
s177AustralasiaSE[lenRatesP1] <- s177AustralasiaSD[lenRatesP1]/
  sqrt(length(which(is.finite(s177AustralasiaS))))

#
# 7. Pacific Islands
#
j <- intersect(which(s177CSCode >= 701), which(s177CSCode < 823))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177PacificIslandsS <- trends[lenRatesP1,j]
s177PacificIslands <- trends[1:lenRates,j]
s177PacificIslandsLat <- stns177Lat[j]
s177PacificIslandsLon <- stns177LonNeg[j]
s177PacificIslandsNames <- stntable$sname[j]

s177PacificIslandsMean <- array(NA,dim=c(lenRatesP1,1))
s177PacificIslandsSD <- array(NA,dim=c(lenRatesP1,1))
s177PacificIslandsSE <- array(NA,dim=c(lenRatesP1,1))
for ( i in 1:lenRates){
  s177PacificIslandsMean[i] <- mean(s177PacificIslands[i,], na.rm=TRUE)
  s177PacificIslandsSD[i] <- sd(s177PacificIslands[i,], na.rm=TRUE)
  s177PacificIslandsSE[i] <- s177PacificIslandsSD[i]/
    sqrt(length(which(is.finite(s177PacificIslands[i,]))))
}
s177PacificIslandsMean[lenRatesP1] <- mean(s177PacificIslandsS, na.rm=TRUE)
s177PacificIslandsSD[lenRatesP1] <- sd(s177PacificIslandsS, na.rm=TRUE)
s177PacificIslandsSE[lenRatesP1] <- s177PacificIslandsSD[lenRatesP1]/
  sqrt(length(which(is.finite(s177PacificIslandsS))))

#
# 8. E Pacific
#
j <- intersect(which(s177CSCode >= 823), which(s177CSCode <= 836))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177EPacificS <- trends[lenRatesP1,j]
s177EPacific <- trends[1:lenRates,j]
s177EPacificLat <- stns177Lat[j]
s177EPacificLon <- stns177LonNeg[j]
s177EPacificNames <- stntable$sname[j]

s177EPacificMean <- array(NA,dim=c(lenRatesP1,1))
s177EPacificSD <- array(NA,dim=c(lenRatesP1,1))
s177EPacificSE <- array(NA,dim=c(lenRatesP1,1))
for ( i in 1:lenRates){
  s177EPacificMean[i] <- mean(s177EPacific[i,], na.rm=TRUE)
  s177EPacificSD[i] <- sd(s177EPacific[i,], na.rm=TRUE)
  s177EPacificSE[i] <- s177EPacificSD[i]/
    sqrt(length(which(is.finite(s177EPacific[i,]))))
}
s177EPacificMean[lenRatesP1] <- mean(s177EPacificS, na.rm=TRUE)
s177EPacificSD[lenRatesP1] <- sd(s177EPacificS, na.rm=TRUE)
s177EPacificSE[lenRatesP1] <- s177EPacificSD[lenRatesP1]/
  sqrt(length(which(is.finite(s177EPacificS))))

#
# 9. SE Pacific
#
j <- intersect(which(s177CSCode >= 840), which(s177CSCode < 851))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177SEPacificS <- trends[lenRatesP1,j]
s177SEPacific <- trends[1:lenRates,j]
s177SEPacificLat <- stns177Lat[j]
s177SEPacificLon <- stns177LonNeg[j]
s177SEPacificNames <- stntable$sname[j]

s177SEPacificMean <- array(NA,dim=c(lenRatesP1,1))
s177SEPacificSD <- array(NA,dim=c(lenRatesP1,1))
s177SEPacificSE <- array(NA,dim=c(lenRatesP1,1))
for ( i in 1:lenRates){
  s177SEPacificMean[i] <- mean(s177SEPacific[i,], na.rm=TRUE)
  s177SEPacificSD[i] <- sd(s177SEPacific[i,], na.rm=TRUE)
  s177SEPacificSE[i] <- s177SEPacificSD[i]/
    sqrt(length(which(is.finite(s177SEPacific[i,]))))
}
s177SEPacificMean[lenRatesP1] <- mean(s177SEPacificS, na.rm=TRUE)
s177SEPacificSD[lenRatesP1] <- sd(s177SEPacificS, na.rm=TRUE)
s177SEPacificSE[lenRatesP1] <- s177SEPacificSD[lenRatesP1]/
  sqrt(length(which(is.finite(s177SEPacificS))))

#
# 10. SW Atlantic
#
j <- intersect(which(s177CSCode >= 860), which(s177CSCode < 875))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177SWAtlanticS <- trends[lenRatesP1,j]
s177SWAtlantic <- trends[1:lenRates,j]
s177SWAtlanticLat <- stns177Lat[j]
s177SWAtlanticLon <- stns177LonNeg[j]
s177SWAtlanticNames <- stntable$sname[j]

s177SWAtlanticMean <- array(NA,dim=c(lenRatesP1,1))
s177SWAtlanticSD <- array(NA,dim=c(lenRatesP1,1))
s177SWAtlanticSE <- array(NA,dim=c(lenRatesP1,1))
for ( i in 1:lenRates){
  s177SWAtlanticMean[i] <- mean(s177SWAtlantic[i,], na.rm=TRUE)
  s177SWAtlanticSD[i] <- sd(s177SWAtlantic[i,], na.rm=TRUE)
  s177SWAtlanticSE[i] <- s177SWAtlanticSD[i]/
    sqrt(length(which(is.finite(s177SWAtlantic[i,]))))
}
s177SWAtlanticMean[lenRatesP1] <- mean(s177SWAtlanticS, na.rm=TRUE)
s177SWAtlanticSD[lenRatesP1] <- sd(s177SWAtlanticS, na.rm=TRUE)
s177SWAtlanticSE[lenRatesP1] <- s177SWAtlanticSD[lenRatesP1]/
  sqrt(length(which(is.finite(s177SWAtlanticS))))

#
# 11. US Gulf
#
j <- intersect(which(s177CSCode >= 938), which(s177CSCode < 943))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177USGulfS <- trends[lenRatesP1,j]
s177USGulf <- trends[1:lenRates,j]
s177USGulfLat <- stns177Lat[j]
s177USGulfLon <- stns177LonNeg[j]
s177USGulfNames <- stntable$sname[j]

s177USGulfMean <- array(NA,dim=c(lenRatesP1,1))
s177USGulfSD <- array(NA,dim=c(lenRatesP1,1))
s177USGulfSE <- array(NA,dim=c(lenRatesP1,1))
for ( i in 1:lenRates){
  s177USGulfMean[i] <- mean(s177USGulf[i,], na.rm=TRUE)
  s177USGulfSD[i] <- sd(s177USGulf[i,], na.rm=TRUE)
  s177USGulfSE[i] <- s177USGulfSD[i]/
    sqrt(length(which(is.finite(s177USGulf[i,]))))
}
s177USGulfMean[lenRatesP1] <- mean(s177USGulfS, na.rm=TRUE)
s177USGulfSD[lenRatesP1] <- sd(s177USGulfS, na.rm=TRUE)
s177USGulfSE[lenRatesP1] <- s177USGulfSD[lenRatesP1]/sqrt(length(which(is.finite(s177USGulfS))))

#
# 12. W Atlantic
#
j <- intersect(which(s177CSCode >= 950), which(s177CSCode < 961))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177WAtlanticS <- trends[lenRatesP1,j]
s177WAtlantic <- trends[1:lenRates,j]
s177WAtlanticLat <- stns177Lat[j]
s177WAtlanticLon <- stns177LonNeg[j]
s177WAtlanticNames <- stntable$sname[j]

s177WAtlanticMean <- array(NA,dim=c(lenRatesP1,1))
s177WAtlanticSD <- array(NA,dim=c(lenRatesP1,1))
s177WAtlanticSE <- array(NA,dim=c(lenRatesP1,1))
for ( i in 1:lenRates){
  s177WAtlanticMean[i] <- mean(s177WAtlantic[i,], na.rm=TRUE)
  s177WAtlanticSD[i] <- sd(s177WAtlantic[i,], na.rm=TRUE)
  s177WAtlanticSE[i] <- s177WAtlanticSD[i]/
    sqrt(length(which(is.finite(s177WAtlantic[i,]))))
}
s177WAtlanticMean[lenRatesP1] <- mean(s177WAtlanticS, na.rm=TRUE)
s177WAtlanticSD[lenRatesP1] <- sd(s177WAtlanticS, na.rm=TRUE)
s177WAtlanticSE[lenRatesP1] <- s177WAtlanticSD[lenRatesP1]/
  sqrt(length(which(is.finite(s177WAtlanticS))))

#
# 13. NW Atlantic
#
j <- which(s177CSCode >= 970)
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177NWAtlanticS <- trends[lenRatesP1,j]
s177NWAtlantic <- trends[1:lenRates,j]
s177NWAtlanticLat <- stns177Lat[j]
s177NWAtlanticLon <- stns177LonNeg[j]
s177NWAtlanticNames <- stntable$sname[j]

s177NWAtlanticMean <- array(NA,dim=c(lenRatesP1,1))
s177NWAtlanticSD <- array(NA,dim=c(lenRatesP1,1))
s177NWAtlanticSE <- array(NA,dim=c(lenRatesP1,1))
for ( i in 1:lenRates){
  s177NWAtlanticMean[i] <- mean(s177NWAtlantic[i,], na.rm=TRUE)
  s177NWAtlanticSD[i] <- sd(s177NWAtlantic[i,], na.rm=TRUE)
  s177NWAtlanticSE[i] <- s177NWAtlanticSD[i]/
    sqrt(length(which(is.finite(s177NWAtlantic[i,]))))
}
s177NWAtlanticMean[lenRatesP1] <- mean(s177NWAtlanticS, na.rm=TRUE)
s177NWAtlanticSD[lenRatesP1] <- sd(s177NWAtlanticS, na.rm=TRUE)
s177NWAtlanticSE[lenRatesP1] <- s177NWAtlanticSD[lenRatesP1]/
  sqrt(length(which(is.finite(s177NWAtlanticS))))

#
# Global
#
s177GlobalArray <- array(c(s177ScandinaviaMean, s177NEuropeMean, 
  s177EAtlanticMean, s177MediterraneanMean, s177SEAsiaMean, 
  s177AustralasiaMean, s177PacificIslandsMean, s177EPacificMean, 
  s177SEPacificMean, s177SWAtlanticMean, s177USGulfMean, s177WAtlanticMean,
  s177NWAtlanticMean),dim=c(lenRatesP1,13))
s177GlobalMean <- array(NA,dim=c(lenRatesP1,1))
s177GlobalMeanNo6 <- array(NA,dim=c(lenRatesP1,1))
s177GlobalSD <- array(NA,dim=c(lenRatesP1,1))
s177GlobalSE <- array(NA,dim=c(lenRatesP1,1))

for ( i in 1:lenRatesP1){
  s177GlobalMean[i] <- mean(s177GlobalArray[i,], na.rm=TRUE)
  s177GlobalMeanNo6[i] <- mean(s177GlobalArray[i,c(1:5,7:13)], na.rm=TRUE)
  s177GlobalSD[i] <- sd(s177GlobalArray[i,], na.rm=TRUE)
  s177GlobalSE[i] <- s177GlobalSD[i]/
    sqrt(length(which(is.finite(s177GlobalArray[i,]))))
}
# Write out station names
region_names <- ls(pat="Names")
for ( i in 1:13 ){
write.table(file=paste(region_names[i],".txt",sep=""), get(region_names[i]))
}

x11()
plot(trendMidPoints, s177GlobalMean, col='black', lwd=2, type='l', 
  lty='dashed', ylim=range(s177GlobalArray, na.rm=T))
lines(trendMidPoints, s177GlobalArray[,1], col='red', lwd=2)
lines(trendMidPoints, s177GlobalArray[,2], col='blue', lwd=2)
lines(trendMidPoints, s177GlobalArray[,3], col='magenta', lwd=2)
lines(trendMidPoints, s177GlobalArray[,4], col='cyan', lwd=2)
lines(trendMidPoints, s177GlobalArray[,5], col='green', lwd=2)
lines(trendMidPoints, s177GlobalArray[,6], col='orange', lwd=2)
lines(trendMidPoints, s177GlobalArray[,7], col='darkblue', lwd=2, lty='dotdash')
lines(trendMidPoints, s177GlobalArray[,8], col='grey', lwd=2, lty='dotdash')
lines(trendMidPoints, s177GlobalArray[,9], col='yellow', lwd=2, lty='dotdash') 
lines(trendMidPoints, s177GlobalArray[,10], col='violet', lwd=2, lty='dotdash')
lines(trendMidPoints, s177GlobalArray[,11], col='brown', lwd=2, lty='dotdash') 
lines(trendMidPoints, s177GlobalArray[,12], col='salmon', lwd=2, lty='dotdash')
lines(trendMidPoints, s177GlobalArray[,13], col='mintcream', lwd=2, lty='dotdash')
lines(trendMidPoints, s177GlobalArray[,13], col='turquoise', lwd=2, lty='dotdash')

save(file="aviso_stns_5yrs.RData", trendMidPoints, trends, s177GlobalMean, s177GlobalMeanNo6, s177GlobalArray)

#
# Check completeness of data
#
completeRegions <- vector(mode="integer", length=lenRatesP1)
completeStations <- vector(mode="integer", length=lenRatesP1)
for (i in 1:lenRatesP1) {
  completeRegions[i] <- length(which(is.finite(s177GlobalArray[i,])))
  completeStations[i] <- length(which(is.finite(trends[i,])))
}
x11()
plot(trendMidPoints, completeRegions, type="s")
x11()
plot(trendMidPoints, completeStations, type="s")
#
