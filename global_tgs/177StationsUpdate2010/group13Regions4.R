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

s177ScandinaviaS <- s177PgrCorrectedRates[58,j]
s177Scandinavia <- s177PgrCorrectedRates[1:57,j]
s177ScandinaviaAlt <- s177PgrCorrectedRates[59,j]

s177ScandinaviaLat <- stns177Lat[j]
s177ScandinaviaLon <- stns177LonNeg[j]
s177ScandinaviaNames <- stns177$sname[j]

s177ScandinaviaMean <- array(NA,dim=c(59,1))
s177ScandinaviaSD <- array(NA,dim=c(59,1))
s177ScandinaviaSE <- array(NA,dim=c(59,1))
for ( i in 1:57){
  s177ScandinaviaMean[i] <- mean(s177Scandinavia[i,], na.rm=TRUE)
  s177ScandinaviaSD[i] <- sd(s177Scandinavia[i,], na.rm=TRUE)
  s177ScandinaviaSE[i] <- s177ScandinaviaSD[i]/
    sqrt(length(which(is.finite(s177Scandinavia[i,]))))
}
s177ScandinaviaMean[58] <- mean(s177ScandinaviaS, na.rm=TRUE)
s177ScandinaviaSD[58] <- sd(s177ScandinaviaS, na.rm=TRUE)
s177ScandinaviaSE[58] <- s177ScandinaviaSD[58]/
  sqrt(length(which(is.finite(s177ScandinaviaS))))

# Altimetry period
s177ScandinaviaMean[59] <- mean(s177ScandinaviaAlt, na.rm=TRUE)
s177ScandinaviaSD[59] <- sd(s177ScandinaviaAlt, na.rm=TRUE)
s177ScandinaviaSE[59] <- s177ScandinaviaSD[59]/
  sqrt(length(which(is.finite(s177ScandinaviaAlt))))

#    j = intersect(find(shortCCode>=140), find(shortCCode<200));
#    j = [j;find(shortCCode==210);find(shortCCode==360)];
j <- intersect(which(s177CSCode >= 140), which(s177CSCode < 200))
j <- union(j, which(s177CSCode == 210))
j <- union(j, which(s177CSCode == 360))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177NEuropeS <- s177PgrCorrectedRates[58,j]
s177NEurope <- s177PgrCorrectedRates[1:57,j]
s177NEuropeAlt <- s177PgrCorrectedRates[59,j]

s177NEuropeLat <- stns177Lat[j]
s177NEuropeLon <- stns177LonNeg[j]
s177NEuropeNames <- stns177$sname[j]

s177NEuropeMean <- array(NA,dim=c(58,1))
s177NEuropeSD <- array(NA,dim=c(58,1))
s177NEuropeSE <- array(NA,dim=c(58,1))
for ( i in 1:57){
  s177NEuropeMean[i] <- mean(s177NEurope[i,], na.rm=TRUE)
  s177NEuropeSD[i] <- sd(s177NEurope[i,], na.rm=TRUE)
  s177NEuropeSE[i] <- s177NEuropeSD[i]/
    sqrt(length(which(is.finite(s177NEurope[i,]))))
}
s177NEuropeMean[58] <- mean(s177NEuropeS, na.rm=TRUE)
s177NEuropeSD[58] <- sd(s177NEuropeS, na.rm=TRUE)
s177NEuropeSE[58] <- s177NEuropeSD[58]/
  sqrt(length(which(is.finite(s177NEuropeS))))


# Altimetry period
s177NEuropeMean[59] <- mean(s177NEuropeAlt, na.rm=TRUE)
s177NEuropeSD[59] <- sd(s177NEuropeAlt, na.rm=TRUE)
s177NEuropeSE[59] <- s177NEuropeSD[59]/
  sqrt(length(which(is.finite(s177NEuropeAlt))))

#    j = intersect(find(shortCCode>=340), find(shortCCode<=447));
#    j = [j; find(shortCCode==200)];
#    j = [j;intersect(find(shortCSCode>=220), find(shortCSCode<=220.011))];
#    j = [j;intersect(find(shortCSCode>=430), find(shortCSCode<=430.081))];
j <- intersect(which(s177CSCode >= 340), which(s177CSCode <= 447))
j <- union(j, intersect(which(s177CSCode >= 200), which(s177CSCode < 201)))
j <- union(j, intersect(which(s177CSCode >= 220), which(s177CSCode <= 220.011)))
j <- union(j, intersect(which(s177CSCode >= 430), which(s177CSCode <= 430.081)))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177EAtlanticS <- s177PgrCorrectedRates[58,j]
s177EAtlantic <- s177PgrCorrectedRates[1:57,j]
s177EAtlanticAlt <- s177PgrCorrectedRates[59,j]

s177EAtlanticLat <- stns177Lat[j]
s177EAtlanticLon <- stns177LonNeg[j]
s177EAtlanticNames <- stns177$sname[j]

s177EAtlanticMean <- array(NA,dim=c(59,1))
s177EAtlanticSD <- array(NA,dim=c(59,1))
s177EAtlanticSE <- array(NA,dim=c(59,1))
for ( i in 1:57){
  s177EAtlanticMean[i] <- mean(s177EAtlantic[i,], na.rm=TRUE)
  s177EAtlanticSD[i] <- sd(s177EAtlantic[i,], na.rm=TRUE)
  s177EAtlanticSE[i] <- s177EAtlanticSD[i]/
    sqrt(length(which(is.finite(s177EAtlantic[i,]))))
}
s177EAtlanticMean[58] <- mean(s177EAtlanticS, na.rm=TRUE)
s177EAtlanticSD[58] <- sd(s177EAtlanticS, na.rm=TRUE)
s177EAtlanticSE[58] <- s177EAtlanticSD[58]/
  sqrt(length(which(is.finite(s177EAtlanticS))))


# Altimetry period
s177EAtlanticMean[59] <- mean(s177EAtlanticAlt, na.rm=TRUE)
s177EAtlanticSD[59] <- sd(s177EAtlanticAlt, na.rm=TRUE)
s177EAtlanticSE[59] <- s177EAtlanticSD[59]/
  sqrt(length(which(is.finite(s177EAtlanticAlt))))

#    j = intersect(find(shortCCode>=225), find(shortCCode<=290));
#    j = [j;intersect(find(shortCSCode>220.011), find(shortCSCode<221))];
j <- intersect(which(s177CSCode >= 225), which(s177CSCode <= 290))
j <- union(j, intersect(which(s177CSCode > 220.011), which(s177CSCode < 221)))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177MediterraneanS <- s177PgrCorrectedRates[58,j]
s177Mediterranean <- s177PgrCorrectedRates[1:57,j]
s177MediterraneanAlt <- s177PgrCorrectedRates[59,j]

s177MediterraneanLat <- stns177Lat[j]
s177MediterraneanLon <- stns177LonNeg[j]
s177MediterraneanNames <- stns177$sname[j]

s177MediterraneanMean <- array(NA,dim=c(59,1))
s177MediterraneanSD <- array(NA,dim=c(59,1))
s177MediterraneanSE <- array(NA,dim=c(59,1))
for ( i in 1:57){
  s177MediterraneanMean[i] <- mean(s177Mediterranean[i,], na.rm=TRUE)
  s177MediterraneanSD[i] <- sd(s177Mediterranean[i,], na.rm=TRUE)
  s177MediterraneanSE[i] <- s177MediterraneanSD[i]/
    sqrt(length(which(is.finite(s177Mediterranean[i,]))))
}
s177MediterraneanMean[58] <- mean(s177MediterraneanS, na.rm=TRUE)
s177MediterraneanSD[58] <- sd(s177MediterraneanS, na.rm=TRUE)
s177MediterraneanSE[58] <- s177MediterraneanSD[58]/
  sqrt(length(which(is.finite(s177MediterraneanS))))

# Altimetry period
s177MediterraneanMean[59] <- mean(s177MediterraneanAlt, na.rm=TRUE)
s177MediterraneanSD[59] <- sd(s177MediterraneanAlt, na.rm=TRUE)
s177MediterraneanSE[59] <- s177MediterraneanSD[59]/
  sqrt(length(which(is.finite(s177MediterraneanAlt))))

#    j = intersect(find(shortCCode>=530), find(shortCCode<=648));
j <- intersect(which(s177CSCode >= 530), which(s177CSCode <= 648))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177SEAsiaS <- s177PgrCorrectedRates[58,j]
s177SEAsia <- s177PgrCorrectedRates[1:57,j]
s177SEAsiaAlt <- s177PgrCorrectedRates[59,j]

s177SEAsiaLat <- stns177Lat[j]
s177SEAsiaLon <- stns177LonNeg[j]
s177SEAsiaNames <- stns177$sname[j]

s177SEAsiaMean <- array(NA,dim=c(59,1))
s177SEAsiaSD <- array(NA,dim=c(59,1))
s177SEAsiaSE <- array(NA,dim=c(59,1))
for ( i in 1:57){
  s177SEAsiaMean[i] <- mean(s177SEAsia[i,], na.rm=TRUE)
  s177SEAsiaSD[i] <- sd(s177SEAsia[i,], na.rm=TRUE)
  s177SEAsiaSE[i] <- s177SEAsiaSD[i]/
    sqrt(length(which(is.finite(s177SEAsia[i,]))))
}
s177SEAsiaMean[58] <- mean(s177SEAsiaS, na.rm=TRUE)
s177SEAsiaSD[58] <- sd(s177SEAsiaS, na.rm=TRUE)
s177SEAsiaSE[58] <- s177SEAsiaSD[58]/sqrt(length(which(is.finite(s177SEAsiaS))))

# Altimetry period
s177SEAsiaMean[59] <- mean(s177SEAsiaAlt, na.rm=TRUE)
s177SEAsiaSD[59] <- sd(s177SEAsiaAlt, na.rm=TRUE)
s177SEAsiaSE[59] <- s177SEAsiaSD[59]/
  sqrt(length(which(is.finite(s177SEAsiaAlt))))

#    j = intersect(find(shortCCode>=680), find(shortCCode<=700));
j <- intersect(which(s177CSCode >= 680), which(s177CSCode <= 700))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177AustralasiaS <- s177PgrCorrectedRates[58,j]
s177Australasia <- s177PgrCorrectedRates[1:57,j]
s177AustralasiaAlt <- s177PgrCorrectedRates[59,j]

s177AustralasiaLat <- stns177Lat[j]
s177AustralasiaLon <- stns177LonNeg[j]
s177AustralasiaNames <- stns177$sname[j]

s177AustralasiaMean <- array(NA,dim=c(59,1))
s177AustralasiaSD <- array(NA,dim=c(59,1))
s177AustralasiaSE <- array(NA,dim=c(59,1))
for ( i in 1:57){
  s177AustralasiaMean[i] <- mean(s177Australasia[i,], na.rm=TRUE)
  s177AustralasiaSD[i] <- sd(s177Australasia[i,], na.rm=TRUE)
  s177AustralasiaSE[i] <- s177AustralasiaSD[i]/ 
    sqrt(length(which(is.finite(s177Australasia[i,]))))
}
s177AustralasiaMean[58] <- mean(s177AustralasiaS, na.rm=TRUE)
s177AustralasiaSD[58] <- sd(s177AustralasiaS, na.rm=TRUE)
s177AustralasiaSE[58] <- s177AustralasiaSD[58]/
  sqrt(length(which(is.finite(s177AustralasiaS))))

# Altimetry period
s177AustralasiaMean[59] <- mean(s177AustralasiaAlt, na.rm=TRUE)
s177AustralasiaSD[59] <- sd(s177AustralasiaAlt, na.rm=TRUE)
s177AustralasiaSE[59] <- s177AustralasiaSD[59]/
  sqrt(length(which(is.finite(s177AustralasiaAlt))))

#    j = intersect(find(shortCCode>=701), find(shortCCode<823));
j <- intersect(which(s177CSCode >= 701), which(s177CSCode < 823))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177PacificIslandsS <- s177PgrCorrectedRates[58,j]
s177PacificIslands <- s177PgrCorrectedRates[1:57,j]
s177PacificIslandsAlt <- s177PgrCorrectedRates[59,j]

s177PacificIslandsLat <- stns177Lat[j]
s177PacificIslandsLon <- stns177LonNeg[j]
s177PacificIslandsNames <- stns177$sname[j]

s177PacificIslandsMean <- array(NA,dim=c(59,1))
s177PacificIslandsSD <- array(NA,dim=c(59,1))
s177PacificIslandsSE <- array(NA,dim=c(59,1))
for ( i in 1:57){
  s177PacificIslandsMean[i] <- mean(s177PacificIslands[i,], na.rm=TRUE)
  s177PacificIslandsSD[i] <- sd(s177PacificIslands[i,], na.rm=TRUE)
  s177PacificIslandsSE[i] <- s177PacificIslandsSD[i]/
    sqrt(length(which(is.finite(s177PacificIslands[i,]))))
}
s177PacificIslandsMean[58] <- mean(s177PacificIslandsS, na.rm=TRUE)
s177PacificIslandsSD[58] <- sd(s177PacificIslandsS, na.rm=TRUE)
s177PacificIslandsSE[58] <- s177PacificIslandsSD[58]/
  sqrt(length(which(is.finite(s177PacificIslandsS))))

# Altimetry period
s177PacificIslandsMean[59] <- mean(s177PacificIslandsAlt, na.rm=TRUE)
s177PacificIslandsSD[59] <- sd(s177PacificIslandsAlt, na.rm=TRUE)
s177PacificIslandsSE[59] <- s177PacificIslandsSD[59]/
  sqrt(length(which(is.finite(s177PacificIslandsAlt))))

#    j = intersect(find(shortCCode>=823), find(shortCCode<=836));
j <- intersect(which(s177CSCode >= 823), which(s177CSCode <= 836))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177EPacificS <- s177PgrCorrectedRates[58,j]
s177EPacific <- s177PgrCorrectedRates[1:57,j]
s177EPacificAlt <- s177PgrCorrectedRates[59,j]

s177EPacificLat <- stns177Lat[j]
s177EPacificLon <- stns177LonNeg[j]
s177EPacificNames <- stns177$sname[j]

s177EPacificMean <- array(NA,dim=c(59,1))
s177EPacificSD <- array(NA,dim=c(59,1))
s177EPacificSE <- array(NA,dim=c(59,1))
for ( i in 1:57){
  s177EPacificMean[i] <- mean(s177EPacific[i,], na.rm=TRUE)
  s177EPacificSD[i] <- sd(s177EPacific[i,], na.rm=TRUE)
  s177EPacificSE[i] <- s177EPacificSD[i]/
    sqrt(length(which(is.finite(s177EPacific[i,]))))
}
s177EPacificMean[58] <- mean(s177EPacificS, na.rm=TRUE)
s177EPacificSD[58] <- sd(s177EPacificS, na.rm=TRUE)
s177EPacificSE[58] <- s177EPacificSD[58]/
  sqrt(length(which(is.finite(s177EPacificS))))

# Altimetry period
s177EPacificMean[59] <- mean(s177EPacificAlt, na.rm=TRUE)
s177EPacificSD[59] <- sd(s177EPacificAlt, na.rm=TRUE)
s177EPacificSE[59] <- s177EPacificSD[59]/
  sqrt(length(which(is.finite(s177EPacificAlt))))

#    j = intersect(find(shortCCode>=840), find(shortCCode<=850));
j <- intersect(which(s177CSCode >= 840), which(s177CSCode < 851))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177SEPacificS <- s177PgrCorrectedRates[58,j]
s177SEPacific <- s177PgrCorrectedRates[1:57,j]
s177SEPacificAlt <- s177PgrCorrectedRates[59,j]

s177SEPacificLat <- stns177Lat[j]
s177SEPacificLon <- stns177LonNeg[j]
s177SEPacificNames <- stns177$sname[j]

s177SEPacificMean <- array(NA,dim=c(59,1))
s177SEPacificSD <- array(NA,dim=c(59,1))
s177SEPacificSE <- array(NA,dim=c(59,1))
for ( i in 1:57){
  s177SEPacificMean[i] <- mean(s177SEPacific[i,], na.rm=TRUE)
  s177SEPacificSD[i] <- sd(s177SEPacific[i,], na.rm=TRUE)
  s177SEPacificSE[i] <- s177SEPacificSD[i]/
    sqrt(length(which(is.finite(s177SEPacific[i,]))))
}
s177SEPacificMean[58] <- mean(s177SEPacificS, na.rm=TRUE)
s177SEPacificSD[58] <- sd(s177SEPacificS, na.rm=TRUE)
s177SEPacificSE[58] <- s177SEPacificSD[58]/
  sqrt(length(which(is.finite(s177SEPacificS))))

# Altimetry period
s177SEPacificMean[59] <- mean(s177SEPacificAlt, na.rm=TRUE)
s177SEPacificSD[59] <- sd(s177SEPacificAlt, na.rm=TRUE)
s177SEPacificSE[59] <- s177SEPacificSD[59]/
  sqrt(length(which(is.finite(s177SEPacificAlt))))

#    j = intersect(find(shortCCode>=860), find(shortCCode<=874));
j <- intersect(which(s177CSCode >= 860), which(s177CSCode < 875))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177SWAtlanticS <- s177PgrCorrectedRates[58,j]
s177SWAtlantic <- s177PgrCorrectedRates[1:57,j]
s177SWAtlanticAlt <- s177PgrCorrectedRates[59,j]

s177SWAtlanticLat <- stns177Lat[j]
s177SWAtlanticLon <- stns177LonNeg[j]
s177SWAtlanticNames <- stns177$sname[j]

s177SWAtlanticMean <- array(NA,dim=c(59,1))
s177SWAtlanticSD <- array(NA,dim=c(59,1))
s177SWAtlanticSE <- array(NA,dim=c(59,1))
for ( i in 1:57){
  s177SWAtlanticMean[i] <- mean(s177SWAtlantic[i,], na.rm=TRUE)
  s177SWAtlanticSD[i] <- sd(s177SWAtlantic[i,], na.rm=TRUE)
  s177SWAtlanticSE[i] <- s177SWAtlanticSD[i]/
    sqrt(length(which(is.finite(s177SWAtlantic[i,]))))
}
s177SWAtlanticMean[58] <- mean(s177SWAtlanticS, na.rm=TRUE)
s177SWAtlanticSD[58] <- sd(s177SWAtlanticS, na.rm=TRUE)
s177SWAtlanticSE[58] <- s177SWAtlanticSD[58]/
  sqrt(length(which(is.finite(s177SWAtlanticS))))

# Altimetry period
s177SWAtlanticMean[59] <- mean(s177SWAtlanticAlt, na.rm=TRUE)
s177SWAtlanticSD[59] <- sd(s177SWAtlanticAlt, na.rm=TRUE)
s177SWAtlanticSE[59] <- s177SWAtlanticSD[59]/
  sqrt(length(which(is.finite(s177SWAtlanticAlt))))


#    j = intersect(find(shortCCode>=938), find(shortCCode<=940));
j <- intersect(which(s177CSCode >= 938), which(s177CSCode < 943))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177USGulfS <- s177PgrCorrectedRates[58,j]
s177USGulf <- s177PgrCorrectedRates[1:57,j]
s177USGulfAlt <- s177PgrCorrectedRates[59,j]

s177USGulfLat <- stns177Lat[j]
s177USGulfLon <- stns177LonNeg[j]
s177USGulfNames <- stns177$sname[j]

s177USGulfMean <- array(NA,dim=c(59,1))
s177USGulfSD <- array(NA,dim=c(59,1))
s177USGulfSE <- array(NA,dim=c(59,1))
for ( i in 1:57){
  s177USGulfMean[i] <- mean(s177USGulf[i,], na.rm=TRUE)
  s177USGulfSD[i] <- sd(s177USGulf[i,], na.rm=TRUE)
  s177USGulfSE[i] <- s177USGulfSD[i]/
    sqrt(length(which(is.finite(s177USGulf[i,]))))
}
s177USGulfMean[58] <- mean(s177USGulfS, na.rm=TRUE)
s177USGulfSD[58] <- sd(s177USGulfS, na.rm=TRUE)
s177USGulfSE[58] <- s177USGulfSD[58]/sqrt(length(which(is.finite(s177USGulfS))))

# Altimetry period
s177USGulfMean[59] <- mean(s177USGulfAlt, na.rm=TRUE)
s177USGulfSD[59] <- sd(s177USGulfAlt, na.rm=TRUE)
s177USGulfSE[59] <- s177USGulfSD[59]/
  sqrt(length(which(is.finite(s177USGulfAlt))))


#    j = intersect(find(shortCCode>=950), find(shortCCode<=960));
j <- intersect(which(s177CSCode >= 950), which(s177CSCode < 961))
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177WAtlanticS <- s177PgrCorrectedRates[58,j]
s177WAtlantic <- s177PgrCorrectedRates[1:57,j]
s177WAtlanticAlt <- s177PgrCorrectedRates[59,j]

s177WAtlanticLat <- stns177Lat[j]
s177WAtlanticLon <- stns177LonNeg[j]
s177WAtlanticNames <- stns177$sname[j]

s177WAtlanticMean <- array(NA,dim=c(59,1))
s177WAtlanticSD <- array(NA,dim=c(59,1))
s177WAtlanticSE <- array(NA,dim=c(59,1))
for ( i in 1:57){
  s177WAtlanticMean[i] <- mean(s177WAtlantic[i,], na.rm=TRUE)
  s177WAtlanticSD[i] <- sd(s177WAtlantic[i,], na.rm=TRUE)
  s177WAtlanticSE[i] <- s177WAtlanticSD[i]/
    sqrt(length(which(is.finite(s177WAtlantic[i,]))))
}
s177WAtlanticMean[58] <- mean(s177WAtlanticS, na.rm=TRUE)
s177WAtlanticSD[58] <- sd(s177WAtlanticS, na.rm=TRUE)
s177WAtlanticSE[58] <- s177WAtlanticSD[58]/
  sqrt(length(which(is.finite(s177WAtlanticS))))

# Altimetry period
s177WAtlanticMean[59] <- mean(s177WAtlanticAlt, na.rm=TRUE)
s177WAtlanticSD[59] <- sd(s177WAtlanticAlt, na.rm=TRUE)
s177WAtlanticSE[59] <- s177WAtlanticSD[59]/
  sqrt(length(which(is.finite(s177WAtlanticAlt))))


#    j = find(shortCCode>=970);
j <- which(s177CSCode >= 970)
jCounter <- jCounter+length(j)
jArray<-c(jArray,j)

s177NWAtlanticS <- s177PgrCorrectedRates[58,j]
s177NWAtlantic <- s177PgrCorrectedRates[1:57,j]
s177NWAtlanticAlt <- s177PgrCorrectedRates[59,j]

s177NWAtlanticLat <- stns177Lat[j]
s177NWAtlanticLon <- stns177LonNeg[j]
s177NWAtlanticNames <- stns177$sname[j]

s177NWAtlanticMean <- array(NA,dim=c(59,1))
s177NWAtlanticSD <- array(NA,dim=c(59,1))
s177NWAtlanticSE <- array(NA,dim=c(59,1))
for ( i in 1:57){
  s177NWAtlanticMean[i] <- mean(s177NWAtlantic[i,], na.rm=TRUE)
  s177NWAtlanticSD[i] <- sd(s177NWAtlantic[i,], na.rm=TRUE)
  s177NWAtlanticSE[i] <- s177NWAtlanticSD[i]/
    sqrt(length(which(is.finite(s177NWAtlantic[i,]))))
}
s177NWAtlanticMean[58] <- mean(s177NWAtlanticS, na.rm=TRUE)
s177NWAtlanticSD[58] <- sd(s177NWAtlanticS, na.rm=TRUE)
s177NWAtlanticSE[58] <- s177NWAtlanticSD[58]/
  sqrt(length(which(is.finite(s177NWAtlanticS))))

# Altimetry period
s177NWAtlanticMean[59] <- mean(s177NWAtlanticAlt, na.rm=TRUE)
s177NWAtlanticSD[59] <- sd(s177NWAtlanticAlt, na.rm=TRUE)
s177NWAtlanticSE[59] <- s177NWAtlanticSD[59]/
  sqrt(length(which(is.finite(s177NWAtlanticAlt))))


s177GlobalArray <- array(c(s177ScandinaviaMean, s177NEuropeMean, 
  s177EAtlanticMean, s177MediterraneanMean, s177SEAsiaMean, 
  s177AustralasiaMean, s177PacificIslandsMean, s177EPacificMean, 
  s177SEPacificMean, s177SWAtlanticMean, s177USGulfMean, s177WAtlanticMean,
  s177NWAtlanticMean),dim=c(59,13))
s177GlobalMean <- array(NA,dim=c(59,1))
s177GlobalSD <- array(NA,dim=c(59,1))
s177GlobalSE <- array(NA,dim=c(59,1))
# Robust scale and median
s177GlobalRS <- array(NA,dim=c(59,2))
s177GlobalRS_MAD <- array(NA,dim=c(59,2))

for ( i in 1:59){
	finiteStns<-which(is.finite(s177GlobalArray[i,]))
  s177GlobalMean[i] <- mean(s177GlobalArray[i,], na.rm=TRUE)
  s177GlobalSD[i] <- sd(s177GlobalArray[i,], na.rm=TRUE)
  s177GlobalSE[i] <- s177GlobalSD[i]/
    sqrt(length(finiteStns))
	RS <- s_Sn(s177GlobalArray[i,finiteStns], mu.too = T)
	Med <- median(s177GlobalArray[i,finiteStns])
	MAD <- mad(s177GlobalArray[i,finiteStns])
	s177GlobalRS[i,]<-RS
	s177GlobalRS_MAD[i,]<-c(Med,MAD)
}

# Write out station names

#region_names <- ls(pat="Names")
#for ( i in 1:13 ){
#write.table(file=paste(region_names[i],".txt",sep=""), get(region_names[i]))
#}

