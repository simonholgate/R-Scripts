# Script to group the stations into 13 regions and calculate their means and
# standard deviations etc.

cSCode<-as.integer(stns$ccode)+as.integer(stns$scode)/1000
jCounter<-0
jArray<-0

# From Matlab m-file
#    j = intersect(find(shortCCode>=25), find(shortCCode<=130));
j <- intersect(which(cSCode >= 25), which(cSCode < 131))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

stnsLonNeg <- stnsLon
stnsLonNeg[j]<- -1*(360-stnsLon[j])

scandinaviaS <- pgrCorrectedRates[lenRatesP1,j]
scandinavia <- pgrCorrectedRates[1:lenRates,j]
scandinaviaLat <- stnsLat[j]
scandinaviaLon <- stnsLonNeg[j]
scandinaviaNames <- stns$sname[j]

scandinaviaMean <- array(NA,dim=c(lenRatesP1,1))
scandinaviaSD <- array(NA,dim=c(lenRatesP1,1))
scandinaviaSE <- array(NA,dim=c(lenRatesP1,1))
for ( i in 1:lenRates){
  scandinaviaMean[i] <- mean(scandinavia[i,], na.rm=TRUE)
  scandinaviaSD[i] <- sd(scandinavia[i,], na.rm=TRUE)
  scandinaviaSE[i] <- scandinaviaSD[i]/
    sqrt(length(which(is.finite(scandinavia[i,]))))
}
scandinaviaMean[lenRatesP1] <- mean(scandinaviaS, na.rm=TRUE)
scandinaviaSD[lenRatesP1] <- sd(scandinaviaS, na.rm=TRUE)
scandinaviaSE[lenRatesP1] <- scandinaviaSD[lenRatesP1]/
  sqrt(length(which(is.finite(scandinaviaS))))

#    j = intersect(find(shortCCode>=140), find(shortCCode<200));
#    j = [j;find(shortCCode==210);find(shortCCode==360)];
j <- intersect(which(cSCode >= 140), which(cSCode < 200))
j <- union(j, which(cSCode == 210))
j <- union(j, which(cSCode == 360))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

nEuropeS <- pgrCorrectedRates[lenRatesP1,j]
nEurope <- pgrCorrectedRates[1:lenRates,j]
nEuropeLat <- stnsLat[j]
nEuropeLon <- stnsLonNeg[j]
nEuropeNames <- stns$sname[j]

nEuropeMean <- array(NA,dim=c(lenRatesP1,1))
nEuropeSD <- array(NA,dim=c(lenRatesP1,1))
nEuropeSE <- array(NA,dim=c(lenRatesP1,1))
for ( i in 1:lenRates){
  nEuropeMean[i] <- mean(nEurope[i,], na.rm=TRUE)
  nEuropeSD[i] <- sd(nEurope[i,], na.rm=TRUE)
  nEuropeSE[i] <- nEuropeSD[i]/
    sqrt(length(which(is.finite(nEurope[i,]))))
}
nEuropeMean[lenRatesP1] <- mean(nEuropeS, na.rm=TRUE)
nEuropeSD[lenRatesP1] <- sd(nEuropeS, na.rm=TRUE)
nEuropeSE[lenRatesP1] <- nEuropeSD[lenRatesP1]/
  sqrt(length(which(is.finite(nEuropeS))))

#    j = intersect(find(shortCCode>=340), find(shortCCode<=447));
#    j = [j; find(shortCCode==200)];
#    j = [j;intersect(find(shortCSCode>=220), find(shortCSCode<=220.011))];
#    j = [j;intersect(find(shortCSCode>=430), find(shortCSCode<=430.081))];
j <- intersect(which(cSCode >= 340), which(cSCode <= 447))
j <- union(j, intersect(which(cSCode >= 200), which(cSCode < 201)))
j <- union(j, intersect(which(cSCode >= 220), which(cSCode <= 220.011)))
j <- union(j, intersect(which(cSCode >= 430), which(cSCode <= 430.081)))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

eAtlanticS <- pgrCorrectedRates[lenRatesP1,j]
eAtlantic <- pgrCorrectedRates[1:lenRates,j]
eAtlanticLat <- stnsLat[j]
eAtlanticLon <- stnsLonNeg[j]
eAtlanticNames <- stns$sname[j]

eAtlanticMean <- array(NA,dim=c(lenRatesP1,1))
eAtlanticSD <- array(NA,dim=c(lenRatesP1,1))
eAtlanticSE <- array(NA,dim=c(lenRatesP1,1))
for ( i in 1:lenRates){
  eAtlanticMean[i] <- mean(eAtlantic[i,], na.rm=TRUE)
  eAtlanticSD[i] <- sd(eAtlantic[i,], na.rm=TRUE)
  eAtlanticSE[i] <- eAtlanticSD[i]/
    sqrt(length(which(is.finite(eAtlantic[i,]))))
}
eAtlanticMean[lenRatesP1] <- mean(eAtlanticS, na.rm=TRUE)
eAtlanticSD[lenRatesP1] <- sd(eAtlanticS, na.rm=TRUE)
eAtlanticSE[lenRatesP1] <- eAtlanticSD[lenRatesP1]/
  sqrt(length(which(is.finite(eAtlanticS))))

#    j = intersect(find(shortCCode>=225), find(shortCCode<=290));
#    j = [j;intersect(find(shortCSCode>220.011), find(shortCSCode<221))];
j <- intersect(which(cSCode >= 225), which(cSCode <= 290))
j <- union(j, intersect(which(cSCode > 220.011), which(cSCode < 221)))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

mediterraneanS <- pgrCorrectedRates[lenRatesP1,j]
mediterranean <- pgrCorrectedRates[1:lenRates,j]
mediterraneanLat <- stnsLat[j]
mediterraneanLon <- stnsLonNeg[j]
mediterraneanNames <- stns$sname[j]

mediterraneanMean <- array(NA,dim=c(lenRatesP1,1))
mediterraneanSD <- array(NA,dim=c(lenRatesP1,1))
mediterraneanSE <- array(NA,dim=c(lenRatesP1,1))
for ( i in 1:lenRates){
  mediterraneanMean[i] <- mean(mediterranean[i,], na.rm=TRUE)
  mediterraneanSD[i] <- sd(mediterranean[i,], na.rm=TRUE)
  mediterraneanSE[i] <- mediterraneanSD[i]/
    sqrt(length(which(is.finite(mediterranean[i,]))))
}
mediterraneanMean[lenRatesP1] <- mean(mediterraneanS, na.rm=TRUE)
mediterraneanSD[lenRatesP1] <- sd(mediterraneanS, na.rm=TRUE)
mediterraneanSE[lenRatesP1] <- mediterraneanSD[lenRatesP1]/
  sqrt(length(which(is.finite(mediterraneanS))))

#    j = intersect(find(shortCCode>=530), find(shortCCode<=648));
j <- intersect(which(cSCode >= 530), which(cSCode <= 648))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

sEAsiaS <- pgrCorrectedRates[lenRatesP1,j]
sEAsia <- pgrCorrectedRates[1:lenRates,j]
sEAsiaLat <- stnsLat[j]
sEAsiaLon <- stnsLonNeg[j]
sEAsiaNames <- stns$sname[j]

sEAsiaMean <- array(NA,dim=c(lenRatesP1,1))
sEAsiaSD <- array(NA,dim=c(lenRatesP1,1))
sEAsiaSE <- array(NA,dim=c(lenRatesP1,1))
for ( i in 1:lenRates){
  sEAsiaMean[i] <- mean(sEAsia[i,], na.rm=TRUE)
  sEAsiaSD[i] <- sd(sEAsia[i,], na.rm=TRUE)
  sEAsiaSE[i] <- sEAsiaSD[i]/
    sqrt(length(which(is.finite(sEAsia[i,]))))
}
sEAsiaMean[lenRatesP1] <- mean(sEAsiaS, na.rm=TRUE)
sEAsiaSD[lenRatesP1] <- sd(sEAsiaS, na.rm=TRUE)
sEAsiaSE[lenRatesP1] <- sEAsiaSD[lenRatesP1]/sqrt(length(which(is.finite(sEAsiaS))))

#    j = intersect(find(shortCCode>=680), find(shortCCode<=700));
j <- intersect(which(cSCode >= 680), which(cSCode <= 700))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

australasiaS <- pgrCorrectedRates[lenRatesP1,j]
australasia <- pgrCorrectedRates[1:lenRates,j]
australasiaLat <- stnsLat[j]
australasiaLon <- stnsLonNeg[j]
australasiaNames <- stns$sname[j]

australasiaMean <- array(NA,dim=c(lenRatesP1,1))
australasiaSD <- array(NA,dim=c(lenRatesP1,1))
australasiaSE <- array(NA,dim=c(lenRatesP1,1))
for ( i in 1:lenRates){
  australasiaMean[i] <- mean(australasia[i,], na.rm=TRUE)
  australasiaSD[i] <- sd(australasia[i,], na.rm=TRUE)
  australasiaSE[i] <- australasiaSD[i]/ 
    sqrt(length(which(is.finite(australasia[i,]))))
}
australasiaMean[lenRatesP1] <- mean(australasiaS, na.rm=TRUE)
australasiaSD[lenRatesP1] <- sd(australasiaS, na.rm=TRUE)
australasiaSE[lenRatesP1] <- australasiaSD[lenRatesP1]/
  sqrt(length(which(is.finite(australasiaS))))

#    j = intersect(find(shortCCode>=701), find(shortCCode<823));
j <- intersect(which(cSCode >= 701), which(cSCode < 823))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

pacificIslandsS <- pgrCorrectedRates[lenRatesP1,j]
pacificIslands <- pgrCorrectedRates[1:lenRates,j]
pacificIslandsLat <- stnsLat[j]
pacificIslandsLon <- stnsLonNeg[j]
pacificIslandsNames <- stns$sname[j]

pacificIslandsMean <- array(NA,dim=c(lenRatesP1,1))
pacificIslandsSD <- array(NA,dim=c(lenRatesP1,1))
pacificIslandsSE <- array(NA,dim=c(lenRatesP1,1))
for ( i in 1:lenRates){
  pacificIslandsMean[i] <- mean(pacificIslands[i,], na.rm=TRUE)
  pacificIslandsSD[i] <- sd(pacificIslands[i,], na.rm=TRUE)
  pacificIslandsSE[i] <- pacificIslandsSD[i]/
    sqrt(length(which(is.finite(pacificIslands[i,]))))
}
pacificIslandsMean[lenRatesP1] <- mean(pacificIslandsS, na.rm=TRUE)
pacificIslandsSD[lenRatesP1] <- sd(pacificIslandsS, na.rm=TRUE)
pacificIslandsSE[lenRatesP1] <- pacificIslandsSD[lenRatesP1]/
  sqrt(length(which(is.finite(pacificIslandsS))))

#    j = intersect(find(shortCCode>=823), find(shortCCode<=836));
j <- intersect(which(cSCode >= 823), which(cSCode <= 836))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

ePacificS <- pgrCorrectedRates[lenRatesP1,j]
ePacific <- pgrCorrectedRates[1:lenRates,j]
ePacificLat <- stnsLat[j]
ePacificLon <- stnsLonNeg[j]
ePacificNames <- stns$sname[j]

ePacificMean <- array(NA,dim=c(lenRatesP1,1))
ePacificSD <- array(NA,dim=c(lenRatesP1,1))
ePacificSE <- array(NA,dim=c(lenRatesP1,1))
for ( i in 1:lenRates){
  ePacificMean[i] <- mean(ePacific[i,], na.rm=TRUE)
  ePacificSD[i] <- sd(ePacific[i,], na.rm=TRUE)
  ePacificSE[i] <- ePacificSD[i]/
    sqrt(length(which(is.finite(ePacific[i,]))))
}
ePacificMean[lenRatesP1] <- mean(ePacificS, na.rm=TRUE)
ePacificSD[lenRatesP1] <- sd(ePacificS, na.rm=TRUE)
ePacificSE[lenRatesP1] <- ePacificSD[lenRatesP1]/
  sqrt(length(which(is.finite(ePacificS))))

#    j = intersect(find(shortCCode>=840), find(shortCCode<=850));
j <- intersect(which(cSCode >= 840), which(cSCode < 851))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

nWAtlanticS <- pgrCorrectedRates[lenRatesP1,j]
nWAtlantic <- pgrCorrectedRates[1:lenRates,j]
nWAtlanticLat <- stnsLat[j]
nWAtlanticLon <- stnsLonNeg[j]
nWAtlanticNames <- stns$sname[j]

nWAtlanticMean <- array(NA,dim=c(lenRatesP1,1))
nWAtlanticSD <- array(NA,dim=c(lenRatesP1,1))
nWAtlanticSE <- array(NA,dim=c(lenRatesP1,1))

# SE Pacific only has one station in this period
if(length(j)>1){
  for ( i in 1:lenRates){
    nWAtlanticMean[i] <- mean(nWAtlantic[i,], na.rm=TRUE)
    nWAtlanticSD[i] <- sd(nWAtlantic[i,], na.rm=TRUE)
    nWAtlanticSE[i] <- nWAtlanticSD[i]/
      sqrt(length(which(is.finite(nWAtlantic[i,]))))
  }
  nWAtlanticMean[lenRatesP1] <- mean(nWAtlanticS, na.rm=TRUE)
  nWAtlanticSD[lenRatesP1] <- sd(nWAtlanticS, na.rm=TRUE)
  nWAtlanticSE[lenRatesP1] <- nWAtlanticSD[lenRatesP1]/
    sqrt(length(which(is.finite(nWAtlanticS))))
} else {
  nWAtlanticMean <- c(nWAtlantic, nWAtlanticS)
  nWAtlanticSD <- nWAtlanticSD*NA
  nWAtlanticSE <- nWAtlanticSE*NA
}
#    j = intersect(find(shortCCode>=860), find(shortCCode<=874));
j <- intersect(which(cSCode >= 860), which(cSCode < 875))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

sWAtlanticS <- pgrCorrectedRates[lenRatesP1,j]
sWAtlantic <- pgrCorrectedRates[1:lenRates,j]
sWAtlanticLat <- stnsLat[j]
sWAtlanticLon <- stnsLonNeg[j]
sWAtlanticNames <- stns$sname[j]

sWAtlanticMean <- array(NA,dim=c(lenRatesP1,1))
sWAtlanticSD <- array(NA,dim=c(lenRatesP1,1))
sWAtlanticSE <- array(NA,dim=c(lenRatesP1,1))
for ( i in 1:lenRates){
  sWAtlanticMean[i] <- mean(sWAtlantic[i,], na.rm=TRUE)
  sWAtlanticSD[i] <- sd(sWAtlantic[i,], na.rm=TRUE)
  sWAtlanticSE[i] <- sWAtlanticSD[i]/
    sqrt(length(which(is.finite(sWAtlantic[i,]))))
}
sWAtlanticMean[lenRatesP1] <- mean(sWAtlanticS, na.rm=TRUE)
sWAtlanticSD[lenRatesP1] <- sd(sWAtlanticS, na.rm=TRUE)
sWAtlanticSE[lenRatesP1] <- sWAtlanticSD[lenRatesP1]/
  sqrt(length(which(is.finite(sWAtlanticS))))

#    j = intersect(find(shortCCode>=938), find(shortCCode<=940));
j <- intersect(which(cSCode >= 938), which(cSCode < 943))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

uSGulfS <- pgrCorrectedRates[lenRatesP1,j]
uSGulf <- pgrCorrectedRates[1:lenRates,j]
uSGulfLat <- stnsLat[j]
uSGulfLon <- stnsLonNeg[j]
uSGulfNames <- stns$sname[j]

uSGulfMean <- array(NA,dim=c(lenRatesP1,1))
uSGulfSD <- array(NA,dim=c(lenRatesP1,1))
uSGulfSE <- array(NA,dim=c(lenRatesP1,1))
for ( i in 1:lenRates){
  uSGulfMean[i] <- mean(uSGulf[i,], na.rm=TRUE)
  uSGulfSD[i] <- sd(uSGulf[i,], na.rm=TRUE)
  uSGulfSE[i] <- uSGulfSD[i]/
    sqrt(length(which(is.finite(uSGulf[i,]))))
}
uSGulfMean[lenRatesP1] <- mean(uSGulfS, na.rm=TRUE)
uSGulfSD[lenRatesP1] <- sd(uSGulfS, na.rm=TRUE)
uSGulfSE[lenRatesP1] <- uSGulfSD[lenRatesP1]/sqrt(length(which(is.finite(uSGulfS))))

#    j = intersect(find(shortCCode>=950), find(shortCCode<=960));
j <- intersect(which(cSCode >= 950), which(cSCode < 961))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

wAtlanticS <- pgrCorrectedRates[lenRatesP1,j]
wAtlantic <- pgrCorrectedRates[1:lenRates,j]
wAtlanticLat <- stnsLat[j]
wAtlanticLon <- stnsLonNeg[j]
wAtlanticNames <- stns$sname[j]

wAtlanticMean <- array(NA,dim=c(lenRatesP1,1))
wAtlanticSD <- array(NA,dim=c(lenRatesP1,1))
wAtlanticSE <- array(NA,dim=c(lenRatesP1,1))
for ( i in 1:lenRates){
  wAtlanticMean[i] <- mean(wAtlantic[i,], na.rm=TRUE)
  wAtlanticSD[i] <- sd(wAtlantic[i,], na.rm=TRUE)
  wAtlanticSE[i] <- wAtlanticSD[i]/
    sqrt(length(which(is.finite(wAtlantic[i,]))))
}
wAtlanticMean[lenRatesP1] <- mean(wAtlanticS, na.rm=TRUE)
wAtlanticSD[lenRatesP1] <- sd(wAtlanticS, na.rm=TRUE)
wAtlanticSE[lenRatesP1] <- wAtlanticSD[lenRatesP1]/
  sqrt(length(which(is.finite(wAtlanticS))))

#    j = find(shortCCode>=970);
j <- which(cSCode >= 970)
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

nWAtlanticS <- pgrCorrectedRates[lenRatesP1,j]
nWAtlantic <- pgrCorrectedRates[1:lenRates,j]
nWAtlanticLat <- stnsLat[j]
nWAtlanticLon <- stnsLonNeg[j]
nWAtlanticNames <- stns$sname[j]

nWAtlanticMean <- array(NA,dim=c(lenRatesP1,1))
nWAtlanticSD <- array(NA,dim=c(lenRatesP1,1))
nWAtlanticSE <- array(NA,dim=c(lenRatesP1,1))
# NW Atlantic only has one station in this period
if(length(j)>1){
  for ( i in 1:lenRates){
    nWAtlanticMean[i] <- mean(nWAtlantic[i,], na.rm=TRUE)
    nWAtlanticSD[i] <- sd(nWAtlantic[i,], na.rm=TRUE)
    nWAtlanticSE[i] <- nWAtlanticSD[i]/
      sqrt(length(which(is.finite(nWAtlantic[i,]))))
  }
  nWAtlanticMean[lenRatesP1] <- mean(nWAtlanticS, na.rm=TRUE)
  nWAtlanticSD[lenRatesP1] <- sd(nWAtlanticS, na.rm=TRUE)
  nWAtlanticSE[lenRatesP1] <- nWAtlanticSD[lenRatesP1]/
    sqrt(length(which(is.finite(nWAtlanticS))))
} else {
  nWAtlanticMean <- c(nWAtlantic, nWAtlanticS)
  nWAtlanticSD <- nWAtlanticSD*NA
  nWAtlanticSE <- nWAtlanticSE*NA
}


globalArray <- array(c(scandinaviaMean, nEuropeMean, 
  eAtlanticMean, mediterraneanMean, sEAsiaMean, 
  australasiaMean, pacificIslandsMean, ePacificMean, 
  nWAtlanticMean, sWAtlanticMean, uSGulfMean, wAtlanticMean,
  nWAtlanticMean),dim=c(lenRatesP1,13))
globalMean <- array(NA,dim=c(lenRatesP1,1))
#globalMeanNo9 <- array(NA,dim=c(lenRatesP1,1))
globalSD <- array(NA,dim=c(lenRatesP1,1))
globalSE <- array(NA,dim=c(lenRatesP1,1))

for ( i in 1:lenRatesP1){
  globalMean[i] <- mean(globalArray[i,], na.rm=TRUE)
#  globalMeanNo9[i] <- mean(globalArray[i,c(1:8,10:13)], na.rm=TRUE)
  globalSD[i] <- sd(globalArray[i,], na.rm=TRUE)
  globalSE[i] <- globalSD[i]/
    sqrt(length(which(is.finite(globalArray[i,]))))
}
# Write out station names
region_names <- ls(pat="Names")
for ( i in 1:13 ){
write.table(file=paste(region_names[i],".txt",sep=""), get(region_names[i]))
}

