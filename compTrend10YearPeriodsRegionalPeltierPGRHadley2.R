#% compTrend10YearPeriodsRegionalPeltierPGR.m
#% MATLAB script to calculate the pgr corrected sea level over 10 year periods 
#% where there is at least 50 of the last 60 years of data for a selection of stations
#%
#
startYr <- 1948
endYr <- 1992
#
#% Required years of data in shorter section 
nReqYrShort <- 7
#
#% Required years data in longer section 
nReqYrLong <- 38   
longFile = 'hadleyAPC2/temp27.1948-2002.hadley2.apc.sed'
pgrFile = 'pgr.rates.stns.out.sed'
#
longDataFile <- read.table(file=longFile)
longYears <- longDataFile[,1]
#% delete stations with less than nReqYrLong in the 1948-2002 period
longEnough <- which(longYears >= nReqYrLong)
longDataFile <- longDataFile[longEnough,]
#
longYears <- longDataFile[,1]
longCCode <- longDataFile[,2]
longSCode <- longDataFile[,3]
longCSCode <- longCCode+longSCode/1000
#
longLon <- longDataFile[,4]
longLat <- longDataFile[,5]
longTrend <- longDataFile[,6]
#
pgrDataFile <- read.table(pgrFile)
pgrCSCode <- pgrDataFile[,1]
#
#% Regions: nEurope, scandinavia, iceland, mediterranean, eAtlantic, seAsia,
#% australasia, pacificIslands, ePacific, swAtlantic, usGulf, wAtlantic,
#% nwAtlantic
#nEuropeS <- numeric(0) 
#scandinaviaS <- numeric(0) 
#mediterraneanS <- numeric(0) 
#eAtlanticS <- numeric(0)
#seAsiaS <- numeric(0) 
#australasiaS <- numeric(0) 
#pacificIslandsS <- numeric(0) 
#ePacificS <- numeric(0)
#sePacificS <- numeric(0) 
#usGulfS <- numeric(0) 
#wAtlanticS <- numeric(0) 
#nwAtlanticS <- numeric(0)
#swAtlanticS <- numeric(0) 
#icelandS <- numeric(0) 
#antarcticS <- numeric(0) 
#indOcnS <- numeric(0) 
#blackSeaS <- numeric(0)
#
dataList <- list()
#
mnE <- numeric(0) 
mSc <- numeric(0) 
mMe <- numeric(0) 
mseA <- numeric(0) 
mAus <- numeric(0) 
mpI <- numeric(0) 
meP <- numeric(0)
musG <- numeric(0) 
mwA <- numeric(0) 
mnwA <- numeric(0) 
mBlSea <- numeric(0) 
meA <- numeric(0) 
mIc <- numeric(0) 
mswA <- numeric(0)
mInd <- numeric(0) 
mAnt <- numeric(0) 
mseP <- numeric(0)
#
nECt <- numeric(0) 
ScCt <- numeric(0) 
MeCt <- numeric(0) 
seACt <- numeric(0) 
AusCt <- numeric(0) 
pICt <- numeric(0) 
ePCt <- numeric(0)
usGCt <- numeric(0) 
wACt <- numeric(0) 
nwACt <- numeric(0) 
BlSeaCt <- numeric(0) 
eACt <- numeric(0) 
icCt <- numeric(0) 
swACt <- numeric(0) 
antCt <- numeric(0) 
indCt <- numeric(0) 
sePCt <- numeric(0)
#
yearMean<- numeric(0)
yearCount<- numeric(0)
yearSD<- numeric(0)
stationSet<- numeric(0)
#
antarcticLatLon<- numeric(0) 
icelandLatLon<- numeric(0) 
scandinaviaLatLon<- numeric(0) 
nEuropeLatLon<- numeric(0)
eAtlanticLatLon<- numeric(0) 
mediterraneanLatLon<- numeric(0) 
blackSeaLatLon<- numeric(0) 
indOcnLatLon<- numeric(0)
seAsiaLatLon<- numeric(0) 
australasiaLatLon<- numeric(0) 
pacificIslandsLatLon<- numeric(0) 
ePacificLatLon<- numeric(0)
sePacificLatLon<- numeric(0) 
swAtlanticLatLon<- numeric(0) 
usGulfLatLon<- numeric(0) 
wAtlanticLatLon<- numeric(0)
nwAtlanticLatLon<- numeric(0) 
# 
diffCount <- 0
diffTrend<- numeric(0)
diffTrendArray<- numeric(0)
#
for (i in startYr:endYr) {
  shortTrendFile = paste('hadleyAPC2/temp27.',as.character(i),'-',
    as.character(i+9),'.hadley2.apc.sed', sep="")
#
#nEurope <- numeric(0) 
#scandinavia <- numeric(0) 
#mediterranean <- numeric(0) 
#eAtlantic <- numeric(0)
#seAsia <- numeric(0) 
#australasia <- numeric(0) 
#pacificIslands <- numeric(0) 
#ePacific <- numeric(0)
#sePacific <- numeric(0) 
#usGulf <- numeric(0) 
#wAtlantic <- numeric(0) 
#nwAtlantic <- numeric(0)
#swAtlantic <- numeric(0) 
#iceland <- numeric(0) 
#antarctic <- numeric(0) 
#indOcn <- numeric(0) 
#blackSea <- numeric(0)
# 
  rm(nEurope, scandinavia, mediterranean, eAtlantic, 
  seAsia, australasia, pacificIslands, ePacific,
  sePacific, usGulf, wAtlantic, nwAtlantic,
  swAtlantic, iceland, antarctic, indOcn, blackSea)
#% Open file of stations to get data from
  shortDataFile <- read.table(shortTrendFile)
  shortYears <- shortDataFile[,1]
#
#% Ensure that there are at least nReqYrShort years in the record
  longEnough <- which(shortYears>=nReqYrShort)
  shortDataFile <- shortDataFile[longEnough,]
#  
  shortCSCode <- shortDataFile[,2]+shortDataFile[,3]/1000
  shortCCode <- shortDataFile[,2]
#
#% Skip various stations on basis of visual inspection
#% Barentsburg, Murmansk, Narvik, Cuxhaven 2, Malaga
  skipStations <- c(025.001, 030.018, 040.081, 140.012, 220.031,
                  270.054, seq(from=310.000,to=310.036,length=37), 410.001,
                  seq(from=940.000,to=940.040,length=41),970.078)
  skipCountries <- c(010, seq(from=296,to=309,length=14), 
                     seq(from=431,to=500,length=70), 612,
                   620, 625, 630, 660, 820, 821, 822, 999)
#  
  deleteStations <- numeric(0)
  for (j in 1:length(skipStations)){
    foundStation <- which(shortCSCode == skipStations[j])
    deleteStations <- c(deleteStations,foundStation)
  }
  for (j in 1:length(skipCountries)){
    foundStation <- which(shortCCode == skipCountries[j])
    deleteStations <- c(deleteStations,foundStation)
  }
#
  nonDeleteStations <- setdiff(1:length(shortCCode), deleteStations)
  shortDataFile <- shortDataFile[nonDeleteStations,]
#
  shortCCode <- shortDataFile[,2]
  shortSCode <- shortDataFile[,3]
  shortCSCode <- shortCCode+shortSCode/1000;
#
#% Only keep shortData stations which have entry in longData
  keepStations <- numeric(0)
  for (j in 1:length(shortCSCode)){
    foundStation <- which(longCSCode == shortCSCode[j])
    if (length(foundStation)>0){
        keepStations <- c(keepStations,j)
    }
  }
#
  shortDataFile <- shortDataFile[keepStations,]
#
  shortCCode <- shortDataFile[,2]
  shortSCode <- shortDataFile[,3]
  shortCSCode <- shortCCode+shortSCode/1000
  stationSet <- union(stationSet,shortCSCode)
#
  dataArray <- array(NA,dim=c(length(shortCSCode), 7))
  shortLon <- shortDataFile[,4]
  shortLat <- shortDataFile[,5]
  shortTrend <- as.numeric(as.character(shortDataFile[,6]))
  dataArray[,1] <- shortCSCode
  dataArray[,2] <- shortLon
  dataArray[,3] <- shortLat
  dataArray[,4] <- shortDataFile[,1]
  dataArray[,5] <- shortTrend
#
#% Now check there is a corresponding record in the pgr file and delete the
#% others from the pgr file
  keepStations<-numeric(0)
  for (j in 1:length(shortCSCode)) {
    foundStation <- which(pgrCSCode == shortCSCode[j])
    if (length(foundStation)>0) {
        keepStations <- c(keepStations,foundStation)
    } else {
      if (shortCSCode[j]==690.012){
        foundStation <- which(pgrCSCode==690.011)
        keepStations <- c(keepStations,foundStation)
      } else {
        stop(paste('shortCSCode ',as.character(shortCSCode[j]),
          ' not found in pgr file', sep=""))
      }
    }
  }
#        
  pgr <- pgrDataFile[keepStations,2]
  dataArray[,6] <- pgr
  diffTrend <- shortTrend - pgr
#
  if (i==endYr) {
    diffTrendArray <- c(diffTrendArray,c(shortLon,shortLat,diffTrend))
  }
#
  diffCount <- diffCount + length(pgr)
#
#% Regional trend analysis
  if (i==endYr) {
#%    j = find(shortCCode == 999);
#%    antarcticS = diffTrend(j);
#%    antarcticLatLon = [antarcticLatLon;[shortLat(j),shortLon(j)]];
#
#%    j = find(shortCCode == 10);
#%    icelandS = diffTrend(j);
#%    icelandLatLon = [icelandLatLon;[shortLat(j),shortLon(j)]];
#    
    j <- intersect(which(shortCCode>=25), which(shortCCode<=130))
    scandinaviaS <- diffTrend[j]
      scandinaviaLatLon <- c(scandinaviaLatLon,c(shortLat[j],shortLon[j]))
# Scandianvia is region 1
    dataArray[j,7] <- 1
#
    j <- intersect(which(shortCCode>=140), which(shortCCode<200))
    j <- c(j,which(shortCCode==210),which(shortCCode==360))
    nEuropeS <- diffTrend[j]
      nEuropeLatLon <- c(nEuropeLatLon,c(shortLat[j],shortLon[j]))
# nEurope is region 2
    dataArray[j,7] <- 2
#
    j <- intersect(which(shortCCode>=340), which(shortCCode<=427))
    j <- c(j, which(shortCCode==200))
    j <- c(j,intersect(which(shortCSCode>=220), which(shortCSCode<=220.011)))
    j <- c(j,intersect(which(shortCSCode>=430), which(shortCSCode<=430.081)))
    eAtlanticS <- diffTrend[j] 
      eAtlanticLatLon <- c(eAtlanticLatLon,c(shortLat[j],shortLon[j]))
# eAtlantic is region 3
    dataArray[j,7] <- 3
#
    j <- intersect(which(shortCCode>=225), which(shortCCode<=290))
    j <- c(j,intersect(which(shortCSCode>220.011), which(shortCSCode<221)))
    mediterraneanS <- diffTrend[j] 
      mediterraneanLatLon <- c(mediterraneanLatLon, c(shortLat[j],shortLon[j]))
# mediterranean is region 4
    dataArray[j,7] <- 4
#
#%    j = intersect(find(shortCCode>295), find(shortCCode<310));
#%    j = [j;intersect(find(shortCSCode>=310), find(shortCSCode<=310.036))];
#%    blackSeaS = diffTrend(j); 
#%    blackSeaLatLon = [blackSeaLatLon;[shortLat(j),shortLon(j)]];
#
#%    j = intersect(find(shortCCode>430), find(shortCCode<=500));
#%    indOcnS = diffTrend(j); 
#%    indOcnLatLon = [indOcnLatLon;[shortLat(j),shortLon(j)]];
#    
    j <- intersect(which(shortCCode>=530), which(shortCCode<=648))
    seAsiaS <- diffTrend[j]
      seAsiaLatLon <- c(seAsiaLatLon,c(shortLat[j],shortLon[j]))
# seAsia is region 5
    dataArray[j,7] <- 5

# 
    j <- intersect(which(shortCCode>=680), which(shortCCode<=700))
    australasiaS <- diffTrend[j]
      australasiaLatLon <- c(australasiaLatLon,c(shortLat[j],shortLon[j]))
# australasia is region 6
    dataArray[j,7] <- 6
#    
    j <- intersect(which(shortCCode>=701), which(shortCCode<823))
    pacificIslandsS <- diffTrend[j]
      pacificIslandsLatLon <- c(pacificIslandsLatLon,c(shortLat[j],shortLon[j]))
# pacificIslands is region 7
    dataArray[j,7] <- 7
#    
    j <- intersect(which(shortCCode>=823), which(shortCCode<=836))
    ePacificS <- diffTrend[j]
      ePacificLatLon <- c(ePacificLatLon,c(shortLat[j],shortLon[j]))
# ePacific is region 8
    dataArray[j,7] <- 8
#    
    j <- intersect(which(shortCCode>=840), which(shortCCode<=850))
    sePacificS <- diffTrend[j] 
      sePacificLatLon <- c(sePacificLatLon,c(shortLat[j],shortLon[j]))
# sePacific is region 9
    dataArray[j,7] <- 9
#
    j <- intersect(which(shortCCode>=860), which(shortCCode<=874))
    swAtlanticS <- diffTrend[j] 
      swAtlanticLatLon <- c(swAtlanticLatLon,c(shortLat[j],shortLon[j]))
# swAtlantic is region 10
    dataArray[j,7] <- 10
#
    j <- intersect(which(shortCCode>=938), which(shortCCode<=940))
    usGulfS <- diffTrend[j]
      usGulfLatLon <- c(usGulfLatLon,c(shortLat[j],shortLon[j]))
# usGulf is region 11
    dataArray[j,7] <- 11
#
    j <- intersect(which(shortCCode>=950), which(shortCCode<=960))
    wAtlanticS <- diffTrend[j] 
      wAtlanticLatLon <- c(wAtlanticLatLon,c(shortLat[j],shortLon[j]))
# wAtlantic is region 12
    dataArray[j,7] <- 12
#
    j <- which(shortCCode>=970)
    nwAtlanticS <- diffTrend[j]
      nwAtlanticLatLon <- c(nwAtlanticLatLon,c(shortLat[j],shortLon[j]))
# nwAtlantic is region 13
    dataArray[j,7] <- 13
#
   } else {
# 
    j <- intersect(which(shortCCode>=25), which(shortCCode<=130))
    scandinavia <- diffTrend[j]
      scandinaviaLatLon <- c(scandinaviaLatLon,c(shortLat[j],shortLon[j]))
    dataArray[j,7] <- 1
#
    j <- intersect(which(shortCCode>=140), which(shortCCode<200))
    j <- c(j,which(shortCCode==210),which(shortCCode==360))
    nEurope <- diffTrend[j]
      nEuropeLatLon <- c(nEuropeLatLon,c(shortLat[j],shortLon[j]))
    dataArray[j,7] <- 2
#
    j <- intersect(which(shortCCode>=340), which(shortCCode<=427))
    j <- c(j, which(shortCCode==200))
    j <- c(j,intersect(which(shortCSCode>=220), which(shortCSCode<=220.011)))
    j <- c(j,intersect(which(shortCSCode>=430), which(shortCSCode<=430.081)))
    eAtlantic <- diffTrend[j] 
      eAtlanticLatLon <- c(eAtlanticLatLon,c(shortLat[j],shortLon[j]))
    dataArray[j,7] <- 3
#
    j <- intersect(which(shortCCode>=225), which(shortCCode<=290))
    j <- c(j,intersect(which(shortCSCode>220.011), which(shortCSCode<221)))
    mediterranean <- diffTrend[j] 
      mediterraneanLatLon <- c(mediterraneanLatLon, c(shortLat[j],shortLon[j]))
    dataArray[j,7] <- 4
#
    j <- intersect(which(shortCCode>=530), which(shortCCode<=648))
    seAsia <- diffTrend[j]
      seAsiaLatLon <- c(seAsiaLatLon,c(shortLat[j],shortLon[j]))
    dataArray[j,7] <- 5
# 
    j <- intersect(which(shortCCode>=680), which(shortCCode<=700))
    australasia <- diffTrend[j]
      australasiaLatLon <- c(australasiaLatLon,c(shortLat[j],shortLon[j]))
    dataArray[j,7] <- 6
#    
    j <- intersect(which(shortCCode>=701), which(shortCCode<823))
    pacificIslands <- diffTrend[j]
      pacificIslandsLatLon <- c(pacificIslandsLatLon,c(shortLat[j],shortLon[j]))
    dataArray[j,7] <- 7
#    
    j <- intersect(which(shortCCode>=823), which(shortCCode<=836))
    ePacific <- diffTrend[j]
      ePacificLatLon <- c(ePacificLatLon,c(shortLat[j],shortLon[j]))
    dataArray[j,7] <- 8
#    
    j <- intersect(which(shortCCode>=840), which(shortCCode<=850))
    sePacific <- diffTrend[j] 
      sePacificLatLon <- c(sePacificLatLon,c(shortLat[j],shortLon[j]))
    dataArray[j,7] <- 9
#
    j <- intersect(which(shortCCode>=860), which(shortCCode<=874))
    swAtlantic <- diffTrend[j] 
      swAtlanticLatLon <- c(swAtlanticLatLon,c(shortLat[j],shortLon[j]))
    dataArray[j,7] <- 10
#
    j <- intersect(which(shortCCode>=938), which(shortCCode<=940))
    usGulf <- diffTrend[j]
      usGulfLatLon <- c(usGulfLatLon,c(shortLat[j],shortLon[j]))
    dataArray[j,7] <- 11
#
    j <- intersect(which(shortCCode>=950), which(shortCCode<=960))
    wAtlantic <- diffTrend[j] 
      wAtlanticLatLon <- c(wAtlanticLatLon,c(shortLat[j],shortLon[j]))
    dataArray[j,7] <- 12
#
    j <- which(shortCCode>=970)
    nwAtlantic <- diffTrend[j] 
      nwAtlanticLatLon <- c(nwAtlanticLatLon,c(shortLat[j],shortLon[j]))
    dataArray[j,7] <- 13
  } 
#
  junkArray <- numeric(0)
#
  if (i!=endYr) {
#
    junk <- mean(nEurope)
    if (is.finite(junk)){
      mnE <- c(mnE, junk)
      junkArray <- c(junkArray, junk)
      nECt <- c(nECt, (i-startYr + 1))
    }
    junk <- mean(scandinavia)
    if (is.finite(junk)){
      mSc <- c(mSc, junk)
      ScCt <- c(ScCt, (i-startYr + 1)) 
      junkArray <- c(junkArray, junk)
    }
    junk <- mean(mediterranean)
    if (is.finite(junk)){
      mMe <- c(mMe, junk)
      MeCt <- c(MeCt, (i-startYr + 1))
      junkArray <- c(junkArray, junk)
    }
    junk <- mean(eAtlantic)
    if (is.finite(junk)){
      meA <- c(meA, junk)
      eACt <- c(eACt, (i-startYr + 1))
      junkArray <- c(junkArray, junk)
    }
    junk <- mean(seAsia)
    if (is.finite(junk)){
      mseA <- c(mseA, junk)
      seACt <- c(seACt, (i-startYr + 1))
      junkArray <- c(junkArray, junk)
    }
    junk <- mean(australasia)
    if (is.finite(junk)){
      mAus <- c(mAus, junk)
      AusCt <- c(AusCt, (i-startYr + 1))
      junkArray <- c(junkArray, junk)
    }
    junk <- mean(pacificIslands)
    if (is.finite(junk)){
      mpI <- c(mpI, junk)
      pICt <- c(pICt, (i-startYr + 1))
      junkArray <- c(junkArray, junk)
    }
    junk <- mean(ePacific)
    if (is.finite(junk)){
      meP <- c(meP, junk)
      ePCt <- c(ePCt, (i-startYr + 1))
      junkArray <- c(junkArray, junk)
    }
    junk <- mean(usGulf)
    if (is.finite(junk)){ 
      musG <- c(musG, junk)
      usGCt <- c(usGCt, (i-startYr + 1))
      junkArray <- c(junkArray, junk)
    }
    junk <- mean(wAtlantic)
    if (is.finite(junk)){
      mwA <- c(mwA, junk)
      wACt <- c(wACt, (i-startYr + 1))
      junkArray <- c(junkArray, junk)
    }
    junk <- mean(nwAtlantic)
    if (is.finite(junk)){
      mnwA <- c(mnwA, junk)
      nwACt <- c(nwACt, (i-startYr + 1))
      junkArray <- c(junkArray, junk)
    }
    junk <- mean(swAtlantic)
    if (is.finite(junk)){
      mswA <- c(mswA, junk)
      swACt <- c(swACt, (i-startYr + 1))
      junkArray <- c(junkArray, junk)
    }
#%    junk = mean(iceland);
#%    if isfinite(junk), mIc = [mIc, junk];
#%      icCt = [icCt, i-startYr + 1]; 
#%      junkArray = [junkArray, junk];
#%    end
#%    junk = mean(indOcn);
#%    if isfinite(junk), mInd = [mInd, junk];
#%      indCt = [indCt, i-startYr + 1]; 
#%      junkArray = [junkArray, junk];
#%    end
#%    junk = mean(antarctic);
#%    if isfinite(junk), mAnt = [mAnt, junk];
#%      antCt = [antCt, i-startYr + 1]; 
#%      junkArray = [junkArray, junk];
#%    end
#%    junk = mean(blackSea);
#%    if isfinite(junk), mBlSea = [mBlSea, junk];
#%      BlSeaCt = [BlSeaCt, i-startYr + 1]; 
#%      junkArray = [junkArray, junk];
#%    end
    junk <- mean(sePacific)
    if (is.finite(junk)){
      mseP <- c(mseP, junk)
      sePCt <- c(sePCt, (i-startYr + 1))
      junkArray <- c(junkArray, junk)
    }

    yearMean <- c(yearMean,mean(junkArray))
    yearCount <- c(yearCount,length(junkArray))
    yearSD <- c(yearSD,sd(junkArray))
  }

# 
  dataList[(i-startYr+1)] <- list(data=dataArray)

} # for loop
#
#% Write file of all sites included in the study
#% Produce output file of all sites with > longYears data
#FID5 <- file('psmslSitesListR.out','wt');
#FID6 <- file('GMT/globalMapInputR.out','wt');
#[cSC,sName] <- read.table(file='ncepAPC/temp26.1948-2002.ncep.mod',sep=':');
#stationNames <- strvcat(sName{:});
#longLat = num2str(longLat,'%8.2f');
#longLon = num2str(longLon,'%8.2f');
#longYears = num2str(longYears);
#longCCode = num2str(longCCode,'%03.0f');
#longSCode = num2str(longSCode,'%03.0f');
#
#nameCSCode = zeros(length(cSC),1);
#for i = 1:length(cSC)
#  nameCSCode(i) = str2num(cSC{i});
#end
#for i = 1:length(stationSet)
#  j = find(longCSCode==stationSet(i));
#  fwrite(FID6,['    ',longLon(j,:),'    ',longLat(j,:),'    ',...
#        '1.0',' ',longYears(j,:),' ',longCCode(j,:),' ',...
#        longSCode(j,:),010], 'uchar');
#% Find the equivalent line in the names file
#  l = find(nameCSCode==stationSet(i));
#  if l
#    fwrite(FID5,[cSC{l},'  ',stationNames(l,:),' ', longLon(j,:),'    ',...
#           longLat(j,:),'     ',longYears(j,:),010], 'uchar');
#  else
#    disp(['Name not found: ',num2str(stationSet(i))])
#  end
#end
#
#fclose(FID5);
#fclose(FID6);
#
#save (file="compTrend10YearPeriodsRegionalPeltierPGR.RData")
#save -V6 dataV6.mat
#load compTrend10YearPeriodsRegionalPeltierPGR
        
message('Ignored: Venice 270/054 and Takoradi 410/001')
message('Ignored: Taiwan, N and S Korea and Pacific Russia, Phillipines -')
message('612 620 625 630 660')
message('Ignored: Aleutians, Alaska and Canadian Pacific - 820 821 822')
message('Ignored: Gulf Coast (incl. Galveston) before Pensacola (940/041)')

message(paste('Total pgr corrected records: ', as.character(diffCount),
  sep=""));

#% Calculate mean difference of selections
message(paste('Mean trend of pgr corrected records is ',   
  as.character(mean(diffTrend)), sep=""))
message('Regional means population: ')
mmnE <- mean(mnE) 
message(paste('nEurope: ',as.character(mmnE), sep=""))
mmSc <- mean(mSc)
message(paste('scandinavia: ',as.character(mmSc), sep=""))
mmMe <- mean(mMe)
message(paste('mediterranean: ', as.character(mmMe),sep=""))
mmeA <- mean(meA)
message(paste('eAtlantic: ',as.character(mmeA),sep=""))
mmseA <- mean(mseA)
message(paste('seAsia: ',as.character(mmseA), sep=""))
mmAus <- mean(mAus)
message(paste('australasia: ', as.character(mmAus), sep=""))
mmpI <- mean(mpI)
message(paste('pacificIslands: ',as.character(mmpI), sep=""))
mmeP <- mean(meP)
message(paste('ePacific: ',as.character(mmeP), sep=""))
mmseP <- mean(mseP)
message(paste('sePacific: ',as.character(mmseP), sep=""))
mmusG <- mean(musG)
message(paste('usGulf: ',as.character(mmusG), sep=""))
mmwA <- mean(mwA)
message(paste('wAtlantic: ',as.character(mmwA), sep=""))
mmnwA <- mean(mnwA)
message(paste('nwAtlantic: ', as.character(mmnwA), sep=""))
mmswA <- mean(mswA)
message(paste('swAtlantic: ', as.character(mmswA), sep=""))
#mmIc = mean(mIc); message(['iceland: ', num2str(mmIc)])
#mmInd = mean(mInd); message(['indOcn: ', num2str(mmInd)])
#mmAnt = mean(mAnt); message(['antarctic: ', num2str(mmAnt)])
#mmBlSea = mean(mBlSea); message(['blackSea: ', num2str(mmBlSea)])
#
message('Regional means sample: ')
mnES <- mean(nEuropeS)
message(paste('nEuropeS: ',as.character(mnES)))
mScS <- mean(scandinaviaS)
message(paste('scandinaviaS: ',as.character(mScS)))
mMeS <- mean(mediterraneanS)
message(paste('mediterraneanS: ', as.character(mMeS)))
meAS <- mean(eAtlanticS) 
message(paste('eAtlanticS: ',as.character(meAS)))
mseAS <- mean(seAsiaS)
message(paste('seAsiaS: ',as.character(mseAS)))
mAusS <- mean(australasiaS)
message(paste('australasiaS: ', as.character(mAusS)))
mpIS <- mean(pacificIslandsS)
message(paste('pacificIslandsS: ',as.character(mpIS)))
mePS <- mean(ePacificS)
message(paste('ePacificS: ',as.character(mePS)))
msePS <- mean(sePacificS)
message(paste('sePacificS: ',as.character(msePS)))
musGS <- mean(usGulfS)
message(paste('usGulfS: ',as.character(musGS)))
mwAS <- mean(wAtlanticS)
message(paste('wAtlanticS: ',as.character(mwAS)))
mnwAS <- mean(nwAtlanticS)
message(paste('nwAtlanticS: ', as.character(mnwAS)))
mswAS <- mean(swAtlanticS)
message(paste('swAtlanticS: ', as.character(mswAS)))
#mIcS = mean(icelandS); message(['icelandS: ', num2str(mIcS)])
#mIndS = mean(indOcnS); message(['indOcnS: ', num2str(mIndS)])
#mAntS = mean(antarcticS); message(['antarcticS: ', num2str(mAntS)])
#mBlSeaS = mean(blackSeaS); message(['blackSeaS: ', num2str(mBlSeaS)])
#
regArray <- c(mmnE,mmSc,mmMe,mmeA,mmseA,mmAus,mmpI,mmeP,mmseP,mmusG,
            mmwA,mmnwA,mmswA)
#[mmnE;mmSc;mmMe;mmeA;mmseA;mmAus;mmpI;mmeP;mmseP;mmusG;mmwA;mmnwA;mmswA];
mRegArray <- mean(regArray)
sdRegArray <- sd(regArray)
message(paste('mean of regions: ', as.character(mRegArray)))
message(paste('sd  of regions: ', as.character(sdRegArray)))
#
regArrayS <- c(mnES,mScS,mMeS,meAS,mseAS,mAusS,mpIS,mePS,msePS,musGS,mwAS,
             mnwAS,mswAS)
#[mnES;mScS;mMeS;meAS;mseAS;mAusS;mpIS;mePS;msePS;musGS;mwAS;mswAS];
mRegArrayS <- mean(regArrayS)
sdRegArrayS <- sd(regArrayS)
message(paste('mean of regions (sample): ', as.character(mRegArrayS)))
message(paste('sd  of regions (sample): ', as.character(sdRegArrayS)))
#
popArray <- c(mnE,mSc,mMe,meA,mseA,mAus,mpI,meP,mseP,musG,mwA,mnwA,mswA)
#[mnE,mSc,mMe,meA,mseA,mAus,mpI,meP,mseP,musG,mwA,mswA];
sampleArray <- c(mnES,mScS,mMeS,meAS,mseAS,mAusS,mpIS,mePS,msePS,musGS,
               mwAS,mnwAS,mswAS)
#[mnES;mScS;mMeS;meAS;mseAS;mAusS;mpIS;mePS;msePS;musGS;mwAS;mswAS];
yearMean <- c(yearMean,mean(sampleArray))
yearCount <- c(yearCount,length(sampleArray))
yearSD <- c(yearSD,sd(sampleArray))
yearSE <- yearSD/sqrt(yearCount)
# 
# Non parametric test
sigLevel <- 0.95
message(paste('Significance at the ',as.character(sigLevel*100),'% level...'))
ranksum <- wilcox.test(x=popArray,y=sampleArray,conf.level=sigLevel)
if (ranksum$p.value>0.95) {
  message('Populations are significantly different')
} else {
  message('Populations are NOT significantly different')
}
#
# Compare regions with the mean: which region is the most highly correlated?
X <- 1953:1997
data <- cbind(t(X),t(yearMean),t(c(mSc,mScS)),t(c(mnE,mnES)),t(c(mMe,mMeS)),
        t(c(meA,meAS)), t(c(musG,musGS)),t(c(mwA,mwAS)),t(c(mswA,mswAS)),
        t(c(mseA,mseAS)),t(c(mAus,mAusS)), t(c(mpI,mpIS)),t(c(meP,mePS)), 
        t(c(mseP,msePS)), t(c(mnwA,mnwAS)))
#        [mpI,mpIS]',[meP,mePS]', [mseP,msePS]'];
R<-cor(as.matrix(data))
#j<-which(P>0.05)
#Rnan<-R
#Rnan[j]<-na
#R0<-R
#R0[j]<-0
#max(R0[,2])
#message(as.character(which(R0[,2]==max(R0[,2]))))
#
#figure
#imagesc(Rnan)
#colorbar;
#%
#figure
#H=plot(X,yearMean,linspace(1950,2000,length(X)),linspace(1.7,1.7,length(X)),'-.');
#set(H,'linewidth',2);
#hold on 
#H=plot(X,yearMean,'r.');
#set(H,'markersize',24)
#grid on
#%title('Global mean Peltier-GIA corrected sea level trends during overlapping 10 year periods.')
#xlabel('Year')
#%xlabel('Year','fontsize',48)
#ylabel('GIA corrected sea level trend [mm yr^{-1}]')
#%
#figure
#set(gcf,'PaperOrientation','portrait','papertype','a4','PaperPosition', [0.634517 0.634517 19.715 28.4084])
#H1=subplot(3,1,1);
#pos1=[0.13 0.701222 0.775 0.283778];
#set(H1,'Position', pos1)
#H=plot(X,data(:,3)+15,'k',X,data(:,4)+10,'k--',X,data(:,5)+5,'k:',X,data(:,6),'k-.');
#set(H,'linewidth',2);
#set(gca,'YLim',[-10 55], 'Xticklabel',[],'XLim',[1950 2000]);
#legh=legend('Scandinavia +15','N. Europe +10','Mediterranean +5','E Atlantic',2);
#set(legh,'fontsize',16,'color','none','box','off','xaxislocation','top');
#set(legh,'position',[0.13 pos1(2)+pos1(4)-0.178533 0.27625 0.178533])
#grid on
#H2=subplot(3,1,2);
#pos2=[0.13 0.405611 0.775 0.283778];
#set(H2,'Position', pos2)
#H=plot(X,data(:,7)+15,'k',X,data(:,8)+10,'k--',X,data(:,9)+5,'k:',...
#X,data(:,15),'k-.');
#set(H,'linewidth',2);
#set(gca,'YLim',[-10 55], 'Xticklabel',[],'XLim',[1950 2000]);
#ylabel('GIA corrected sea level trend [mm yr^{-1}]')
#legh=legend('Caribbean +15','W. Atlantic +10','S.W. Atlantic +5','N.W. Atlantic',2);
#set(legh,'fontsize',16,'color','none','box','off','xaxislocation','top');
#set(legh,'position',[0.13 pos2(2)+pos2(4)-0.178533 0.27625 0.178533])
#grid on
#H3=subplot(3,1,3);
#pos3=[0.13 0.11 0.775 0.283778];
#set(H3,'Position', pos3)
#H=plot(X,data(:,10)+20,'k',X,data(:,11)+15,'k--',X,data(:,12)+10,'k:',...
#X,data(:,13)+5,'k-.',X,data(:,14),'k');
#set(H,'linewidth',2);
#set(gca,'YLim',[-10 55],'XLim',[1950 2000]);
#xlabel('Year')
#legh=legend('S.E. Asia +20','Australasia +15','Central Pacfic +10','N.E. Pacific +5','S.E. Pacific',2);
#set(legh,'fontsize',16,'color','none','box','off','xaxislocation','top');
#set(legh,'position',[0.13 pos3(2)+pos3(4)-0.178533 0.27625 0.178533])
#grid on
#%
#%cd GMT
#%save -ascii -tabs antarcticLatLon.out antarcticLatLon;
#%save -ascii -tabs icelandLatLon.out icelandLatLon; 
#%save -ascii -tabs scandinaviaLatLon.out scandinaviaLatLon; 
#%save -ascii -tabs nEuropeLatLon.out nEuropeLatLon;
#%save -ascii -tabs eAtlanticLatLon.out eAtlanticLatLon; 
#%save -ascii -tabs mediterraneanLatLon.out mediterraneanLatLon; 
#%save -ascii -tabs blackSeaLatLon.out blackSeaLatLon; 
#%save -ascii -tabs indOcnLatLon.out indOcnLatLon;
#%save -ascii -tabs seAsiaLatLon.out seAsiaLatLon; 
#%save -ascii -tabs australasiaLatLon.out australasiaLatLon; 
#%save -ascii -tabs pacificIslandsLatLon.out pacificIslandsLatLon; 
#%save -ascii -tabs ePacificLatLon.out ePacificLatLon;
#%save -ascii -tabs sePacificLatLon.out sePacificLatLon; 
#%save -ascii -tabs swAtlanticLatLon.out swAtlanticLatLon; 
#%save -ascii -tabs usGulfLatLon.out usGulfLatLon; 
#%save -ascii -tabs wAtlanticLatLon.out wAtlanticLatLon;
#%save -ascii -tabs nwAtlanticLatLon.out nwAtlanticLatLon;
#%save regions antarcticLatLon icelandLatLon scandinaviaLatLon nEuropeLatLon eAtlanticLatLon mediterraneanLatLon blackSeaLatLon indOcnLatLon seAsiaLatLon australasiaLatLon ePacificLatLon sePacificLatLon swAtlanticLatLon usGulfLatLon wAtlanticLatLon nwAtlanticLatLon pacificIslandsLatLon
#%cd .. 
#%
message(paste('Mean of lower and upper halves: ',as.character(mean(yearMean[1:23])),' ',as.character(mean(yearMean[23:4]))))
message(paste('Mean total: ',as.character(mean(yearMean))))
#figure
#%subplot(2,1,1)
#H=errorbar(X,cumsum(yearMean),yearSE);
#set(H,'linewidth',2);
#hold on 
#H=plot(X,cumsum(yearMean),'r.');
#set(H,'markersize',16)
#grid on
#%title('Integral of 10 year average sea level rise rates 1948-2002')
#xlabel('Year')
#ylabel('GIA corrected sea level rise [mm]')
#set(gca,'YLim',[0 80])
#%set(gca,'Xticklabel',[]);
#%subplot(2,1,2)
#%intData=cumsum(data);
#%H=plot(X,intData(:,3)+60,X,intData(:,4)+55,X,intData(:,5)+50,X,intData(:,6)+45,...
#%X,intData(:,7)+40,X,intData(:,8)+35,X,intData(:,9)+30,X,intData(:,10)+25,...
#%X,intData(:,11)+20,X,intData(:,12)+15,X,intData(:,13)+10,X,intData(:,14)+5,...
#%X,intData(:,15));
#%xlabel('Year')
#%ylabel('GIA corrected sea level rise [mm]')
#%grid on;
#%set(gca,'YLim',[-20 190]);
