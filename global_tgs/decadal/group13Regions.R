# Script to group the stations into 13 regions and calculate their means and
# standard deviations etc.

s177CSCode<-as.integer(stns177$ccode)+as.integer(stns177$scode)/1000
jCounter<-0
jArray<-0

# From Matlab m-file
#    j = intersect(find(shortCCode>=25), find(shortCCode<=130));
j <- intersect(which(s177CSCode >= 25), which(s177CSCode < 131))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

stns177LonNeg <- stns177Lon
stns177LonNeg[j]<- -1*(360-stns177Lon[j])

s177ScandinaviaS <- s177PgrCorrectedRates[42,j]
s177Scandinavia <- s177PgrCorrectedRates[1:41,j]
s177ScandinaviaLat <- stns177Lat[j]
s177ScandinaviaLon <- stns177LonNeg[j]

s177ScandinaviaMean <- array(NA,dim=c(42,1))
s177ScandinaviaSD <- array(NA,dim=c(42,1))
s177ScandinaviaSE <- array(NA,dim=c(42,1))
for ( i in 1:41){
  s177ScandinaviaMean[i] <- mean(s177Scandinavia[i,], na.rm=TRUE)
  s177ScandinaviaSD[i] <- sd(s177Scandinavia[i,], na.rm=TRUE)
  s177ScandinaviaSE[i] <- s177ScandinaviaSD[i]/
    sqrt(length(which(is.finite(s177Scandinavia[i,]))))
}
s177ScandinaviaMean[42] <- mean(s177ScandinaviaS, na.rm=TRUE)
s177ScandinaviaSD[42] <- sd(s177ScandinaviaS, na.rm=TRUE)
s177ScandinaviaSE[42] <- s177ScandinaviaSD[42]/
  sqrt(length(which(is.finite(s177ScandinaviaS))))

#    j = intersect(find(shortCCode>=140), find(shortCCode<200));
#    j = [j;find(shortCCode==210);find(shortCCode==360)];
j <- intersect(which(s177CSCode >= 140), which(s177CSCode < 200))
j <- union(j, which(s177CSCode == 210))
j <- union(j, which(s177CSCode == 360))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177NEuropeS <- s177PgrCorrectedRates[42,j]
s177NEurope <- s177PgrCorrectedRates[1:41,j]
s177NEuropeLat <- stns177Lat[j]
s177NEuropeLon <- stns177LonNeg[j]

s177NEuropeMean <- array(NA,dim=c(42,1))
s177NEuropeSD <- array(NA,dim=c(42,1))
s177NEuropeSE <- array(NA,dim=c(42,1))
for ( i in 1:41){
  s177NEuropeMean[i] <- mean(s177NEurope[i,], na.rm=TRUE)
  s177NEuropeSD[i] <- sd(s177NEurope[i,], na.rm=TRUE)
  s177NEuropeSE[i] <- s177NEuropeSD[i]/
    sqrt(length(which(is.finite(s177NEurope[i,]))))
}
s177NEuropeMean[42] <- mean(s177NEuropeS, na.rm=TRUE)
s177NEuropeSD[42] <- sd(s177NEuropeS, na.rm=TRUE)
s177NEuropeSE[42] <- s177NEuropeSD[42]/
  sqrt(length(which(is.finite(s177NEuropeS))))

#    j = intersect(find(shortCCode>=340), find(shortCCode<=427));
#    j = [j; find(shortCCode==200)];
#    j = [j;intersect(find(shortCSCode>=220), find(shortCSCode<=220.011))];
#    j = [j;intersect(find(shortCSCode>=430), find(shortCSCode<=430.081))];
j <- intersect(which(s177CSCode >= 340), which(s177CSCode <= 427))
j <- union(j, intersect(which(s177CSCode >= 200), which(s177CSCode < 201)))
j <- union(j, intersect(which(s177CSCode >= 220), which(s177CSCode <= 220.011)))
j <- union(j, intersect(which(s177CSCode >= 430), which(s177CSCode <= 430.081)))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177EAtlanticS <- s177PgrCorrectedRates[42,j]
s177EAtlantic <- s177PgrCorrectedRates[1:41,j]
s177EAtlanticLat <- stns177Lat[j]
s177EAtlanticLon <- stns177LonNeg[j]

s177EAtlanticMean <- array(NA,dim=c(42,1))
s177EAtlanticSD <- array(NA,dim=c(42,1))
s177EAtlanticSE <- array(NA,dim=c(42,1))
for ( i in 1:41){
  s177EAtlanticMean[i] <- mean(s177EAtlantic[i,], na.rm=TRUE)
  s177EAtlanticSD[i] <- sd(s177EAtlantic[i,], na.rm=TRUE)
  s177EAtlanticSE[i] <- s177EAtlanticSD[i]/
    sqrt(length(which(is.finite(s177EAtlantic[i,]))))
}
s177EAtlanticMean[42] <- mean(s177EAtlanticS, na.rm=TRUE)
s177EAtlanticSD[42] <- sd(s177EAtlanticS, na.rm=TRUE)
s177EAtlanticSE[42] <- s177EAtlanticSD[42]/
  sqrt(length(which(is.finite(s177EAtlanticS))))

#    j = intersect(find(shortCCode>=225), find(shortCCode<=290));
#    j = [j;intersect(find(shortCSCode>220.011), find(shortCSCode<221))];
j <- intersect(which(s177CSCode >= 225), which(s177CSCode <= 290))
j <- union(j, intersect(which(s177CSCode > 220.011), which(s177CSCode < 221)))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177MediterraneanS <- s177PgrCorrectedRates[42,j]
s177Mediterranean <- s177PgrCorrectedRates[1:41,j]
s177MediterraneanLat <- stns177Lat[j]
s177MediterraneanLon <- stns177LonNeg[j]

s177MediterraneanMean <- array(NA,dim=c(42,1))
s177MediterraneanSD <- array(NA,dim=c(42,1))
s177MediterraneanSE <- array(NA,dim=c(42,1))
for ( i in 1:41){
  s177MediterraneanMean[i] <- mean(s177Mediterranean[i,], na.rm=TRUE)
  s177MediterraneanSD[i] <- sd(s177Mediterranean[i,], na.rm=TRUE)
  s177MediterraneanSE[i] <- s177MediterraneanSD[i]/
    sqrt(length(which(is.finite(s177Mediterranean[i,]))))
}
s177MediterraneanMean[42] <- mean(s177MediterraneanS, na.rm=TRUE)
s177MediterraneanSD[42] <- sd(s177MediterraneanS, na.rm=TRUE)
s177MediterraneanSE[42] <- s177MediterraneanSD[42]/
  sqrt(length(which(is.finite(s177MediterraneanS))))

#    j = intersect(find(shortCCode>=530), find(shortCCode<=648));
j <- intersect(which(s177CSCode >= 530), which(s177CSCode <= 648))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177SEAsiaS <- s177PgrCorrectedRates[42,j]
s177SEAsia <- s177PgrCorrectedRates[1:41,j]
s177SEAsiaLat <- stns177Lat[j]
s177SEAsiaLon <- stns177LonNeg[j]

s177SEAsiaMean <- array(NA,dim=c(42,1))
s177SEAsiaSD <- array(NA,dim=c(42,1))
s177SEAsiaSE <- array(NA,dim=c(42,1))
for ( i in 1:41){
  s177SEAsiaMean[i] <- mean(s177SEAsia[i,], na.rm=TRUE)
  s177SEAsiaSD[i] <- sd(s177SEAsia[i,], na.rm=TRUE)
  s177SEAsiaSE[i] <- s177SEAsiaSD[i]/
    sqrt(length(which(is.finite(s177SEAsia[i,]))))
}
s177SEAsiaMean[42] <- mean(s177SEAsiaS, na.rm=TRUE)
s177SEAsiaSD[42] <- sd(s177SEAsiaS, na.rm=TRUE)
s177SEAsiaSE[42] <- s177SEAsiaSD[42]/sqrt(length(which(is.finite(s177SEAsiaS))))

#    j = intersect(find(shortCCode>=680), find(shortCCode<=700));
j <- intersect(which(s177CSCode >= 680), which(s177CSCode <= 700))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177AustralasiaS <- s177PgrCorrectedRates[42,j]
s177Australasia <- s177PgrCorrectedRates[1:41,j]
s177AustralasiaLat <- stns177Lat[j]
s177AustralasiaLon <- stns177LonNeg[j]

s177AustralasiaMean <- array(NA,dim=c(42,1))
s177AustralasiaSD <- array(NA,dim=c(42,1))
s177AustralasiaSE <- array(NA,dim=c(42,1))
for ( i in 1:41){
  s177AustralasiaMean[i] <- mean(s177Australasia[i,], na.rm=TRUE)
  s177AustralasiaSD[i] <- sd(s177Australasia[i,], na.rm=TRUE)
  s177AustralasiaSE[i] <- s177AustralasiaSD[i]/ 
    sqrt(length(which(is.finite(s177Australasia[i,]))))
}
s177AustralasiaMean[42] <- mean(s177AustralasiaS, na.rm=TRUE)
s177AustralasiaSD[42] <- sd(s177AustralasiaS, na.rm=TRUE)
s177AustralasiaSE[42] <- s177AustralasiaSD[42]/
  sqrt(length(which(is.finite(s177AustralasiaS))))

#    j = intersect(find(shortCCode>=701), find(shortCCode<823));
j <- intersect(which(s177CSCode >= 701), which(s177CSCode < 823))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177PacificIslandsS <- s177PgrCorrectedRates[42,j]
s177PacificIslands <- s177PgrCorrectedRates[1:41,j]
s177PacificIslandsLat <- stns177Lat[j]
s177PacificIslandsLon <- stns177LonNeg[j]

s177PacificIslandsMean <- array(NA,dim=c(42,1))
s177PacificIslandsSD <- array(NA,dim=c(42,1))
s177PacificIslandsSE <- array(NA,dim=c(42,1))
for ( i in 1:41){
  s177PacificIslandsMean[i] <- mean(s177PacificIslands[i,], na.rm=TRUE)
  s177PacificIslandsSD[i] <- sd(s177PacificIslands[i,], na.rm=TRUE)
  s177PacificIslandsSE[i] <- s177PacificIslandsSD[i]/
    sqrt(length(which(is.finite(s177PacificIslands[i,]))))
}
s177PacificIslandsMean[42] <- mean(s177PacificIslandsS, na.rm=TRUE)
s177PacificIslandsSD[42] <- sd(s177PacificIslandsS, na.rm=TRUE)
s177PacificIslandsSE[42] <- s177PacificIslandsSD[42]/
  sqrt(length(which(is.finite(s177PacificIslandsS))))

#    j = intersect(find(shortCCode>=823), find(shortCCode<=836));
j <- intersect(which(s177CSCode >= 823), which(s177CSCode <= 836))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177EPacificS <- s177PgrCorrectedRates[42,j]
s177EPacific <- s177PgrCorrectedRates[1:41,j]
s177EPacificLat <- stns177Lat[j]
s177EPacificLon <- stns177LonNeg[j]

s177EPacificMean <- array(NA,dim=c(42,1))
s177EPacificSD <- array(NA,dim=c(42,1))
s177EPacificSE <- array(NA,dim=c(42,1))
for ( i in 1:41){
  s177EPacificMean[i] <- mean(s177EPacific[i,], na.rm=TRUE)
  s177EPacificSD[i] <- sd(s177EPacific[i,], na.rm=TRUE)
  s177EPacificSE[i] <- s177EPacificSD[i]/
    sqrt(length(which(is.finite(s177EPacific[i,]))))
}
s177EPacificMean[42] <- mean(s177EPacificS, na.rm=TRUE)
s177EPacificSD[42] <- sd(s177EPacificS, na.rm=TRUE)
s177EPacificSE[42] <- s177EPacificSD[42]/
  sqrt(length(which(is.finite(s177EPacificS))))

#    j = intersect(find(shortCCode>=840), find(shortCCode<=850));
j <- intersect(which(s177CSCode >= 840), which(s177CSCode < 851))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177SEPacificS <- s177PgrCorrectedRates[42,j]
s177SEPacific <- s177PgrCorrectedRates[1:41,j]
s177SEPacificLat <- stns177Lat[j]
s177SEPacificLon <- stns177LonNeg[j]

s177SEPacificMean <- array(NA,dim=c(42,1))
s177SEPacificSD <- array(NA,dim=c(42,1))
s177SEPacificSE <- array(NA,dim=c(42,1))
for ( i in 1:41){
  s177SEPacificMean[i] <- mean(s177SEPacific[i,], na.rm=TRUE)
  s177SEPacificSD[i] <- sd(s177SEPacific[i,], na.rm=TRUE)
  s177SEPacificSE[i] <- s177SEPacificSD[i]/
    sqrt(length(which(is.finite(s177SEPacific[i,]))))
}
s177SEPacificMean[42] <- mean(s177SEPacificS, na.rm=TRUE)
s177SEPacificSD[42] <- sd(s177SEPacificS, na.rm=TRUE)
s177SEPacificSE[42] <- s177SEPacificSD[42]/
  sqrt(length(which(is.finite(s177SEPacificS))))

#    j = intersect(find(shortCCode>=860), find(shortCCode<=874));
j <- intersect(which(s177CSCode >= 860), which(s177CSCode < 875))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177SWAtlanticS <- s177PgrCorrectedRates[42,j]
s177SWAtlantic <- s177PgrCorrectedRates[1:41,j]
s177SWAtlanticLat <- stns177Lat[j]
s177SWAtlanticLon <- stns177LonNeg[j]

s177SWAtlanticMean <- array(NA,dim=c(42,1))
s177SWAtlanticSD <- array(NA,dim=c(42,1))
s177SWAtlanticSE <- array(NA,dim=c(42,1))
for ( i in 1:41){
  s177SWAtlanticMean[i] <- mean(s177SWAtlantic[i,], na.rm=TRUE)
  s177SWAtlanticSD[i] <- sd(s177SWAtlantic[i,], na.rm=TRUE)
  s177SWAtlanticSE[i] <- s177SWAtlanticSD[i]/
    sqrt(length(which(is.finite(s177SWAtlantic[i,]))))
}
s177SWAtlanticMean[42] <- mean(s177SWAtlanticS, na.rm=TRUE)
s177SWAtlanticSD[42] <- sd(s177SWAtlanticS, na.rm=TRUE)
s177SWAtlanticSE[42] <- s177SWAtlanticSD[42]/
  sqrt(length(which(is.finite(s177SWAtlanticS))))

#    j = intersect(find(shortCCode>=938), find(shortCCode<=940));
j <- intersect(which(s177CSCode >= 938), which(s177CSCode < 941))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177USGulfS <- s177PgrCorrectedRates[42,j]
s177USGulf <- s177PgrCorrectedRates[1:41,j]
s177USGulfLat <- stns177Lat[j]
s177USGulfLon <- stns177LonNeg[j]

s177USGulfMean <- array(NA,dim=c(42,1))
s177USGulfSD <- array(NA,dim=c(42,1))
s177USGulfSE <- array(NA,dim=c(42,1))
for ( i in 1:41){
  s177USGulfMean[i] <- mean(s177USGulf[i,], na.rm=TRUE)
  s177USGulfSD[i] <- sd(s177USGulf[i,], na.rm=TRUE)
  s177USGulfSE[i] <- s177USGulfSD[i]/
    sqrt(length(which(is.finite(s177USGulf[i,]))))
}
s177USGulfMean[42] <- mean(s177USGulfS, na.rm=TRUE)
s177USGulfSD[42] <- sd(s177USGulfS, na.rm=TRUE)
s177USGulfSE[42] <- s177USGulfSD[42]/sqrt(length(which(is.finite(s177USGulfS))))

#    j = intersect(find(shortCCode>=950), find(shortCCode<=960));
j <- intersect(which(s177CSCode >= 950), which(s177CSCode < 961))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177WAtlanticS <- s177PgrCorrectedRates[42,j]
s177WAtlantic <- s177PgrCorrectedRates[1:41,j]
s177WAtlanticLat <- stns177Lat[j]
s177WAtlanticLon <- stns177LonNeg[j]

s177WAtlanticMean <- array(NA,dim=c(42,1))
s177WAtlanticSD <- array(NA,dim=c(42,1))
s177WAtlanticSE <- array(NA,dim=c(42,1))
for ( i in 1:41){
  s177WAtlanticMean[i] <- mean(s177WAtlantic[i,], na.rm=TRUE)
  s177WAtlanticSD[i] <- sd(s177WAtlantic[i,], na.rm=TRUE)
  s177WAtlanticSE[i] <- s177WAtlanticSD[i]/
    sqrt(length(which(is.finite(s177WAtlantic[i,]))))
}
s177WAtlanticMean[42] <- mean(s177WAtlanticS, na.rm=TRUE)
s177WAtlanticSD[42] <- sd(s177WAtlanticS, na.rm=TRUE)
s177WAtlanticSE[42] <- s177WAtlanticSD[42]/
  sqrt(length(which(is.finite(s177WAtlanticS))))

#    j = find(shortCCode>=970);
j <- which(s177CSCode >= 970)
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177NWAtlanticS <- s177PgrCorrectedRates[42,j]
s177NWAtlantic <- s177PgrCorrectedRates[1:41,j]
s177NWAtlanticLat <- stns177Lat[j]
s177NWAtlanticLon <- stns177LonNeg[j]

s177NWAtlanticMean <- array(NA,dim=c(42,1))
s177NWAtlanticSD <- array(NA,dim=c(42,1))
s177NWAtlanticSE <- array(NA,dim=c(42,1))
for ( i in 1:41){
  s177NWAtlanticMean[i] <- mean(s177NWAtlantic[i,], na.rm=TRUE)
  s177NWAtlanticSD[i] <- sd(s177NWAtlantic[i,], na.rm=TRUE)
  s177NWAtlanticSE[i] <- s177NWAtlanticSD[i]/
    sqrt(length(which(is.finite(s177NWAtlantic[i,]))))
}
s177NWAtlanticMean[42] <- mean(s177NWAtlanticS, na.rm=TRUE)
s177NWAtlanticSD[42] <- sd(s177NWAtlanticS, na.rm=TRUE)
s177NWAtlanticSE[42] <- s177NWAtlanticSD[42]/
  sqrt(length(which(is.finite(s177NWAtlanticS))))

s177GlobalArray <- array(c(s177ScandinaviaMean, s177NEuropeMean, 
  s177EAtlanticMean, s177MediterraneanMean, s177SEAsiaMean, 
  s177AustralasiaMean, s177PacificIslandsMean, s177EPacificMean, 
  s177SEPacificMean, s177SWAtlanticMean, s177USGulfMean, s177WAtlanticMean,
  s177NWAtlanticMean),dim=c(42,13))
s177GlobalMean <- array(NA,dim=c(42,1))
s177GlobalSD <- array(NA,dim=c(42,1))
s177GlobalSE <- array(NA,dim=c(42,1))

for ( i in 1:42){
  s177GlobalMean[i] <- mean(s177GlobalArray[i,], na.rm=TRUE)
  s177GlobalSD[i] <- sd(s177GlobalArray[i,], na.rm=TRUE)
  s177GlobalSE[i] <- s177GlobalSD[i]/
    sqrt(length(which(is.finite(s177GlobalArray[i,]))))
}

