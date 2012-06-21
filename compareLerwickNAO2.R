# Script showing exploration of the rate of sea level at Lerwick and it's
# relation to the NAO and to model data
#source('~/workspace/RScripts/getLerwickData2.R')
load('lerwickData.RData')
annualYrs <- seq.Date(from=as.Date("1957/7/1"),
 to=as.Date("2007/7/1"), by="1 year")
monthlyYrs <- seq.Date(from=as.Date("1957/1/15"), 
 length.out=612, by="1 month")
#monthlyYrs<-seq(from=(1957+0.5/12),to=(2004+11.5/12),by=1/12)

x11()
par(family='HersheySans')
plot(monthlyYrs,lerwickMonthly-mean(lerwickMonthly, na.rm=TRUE),type='l',col='red', 
  lwd=2, ann=F, ylim=c(-300,430))
title(ylab="Sea Level [mm]", xlab="Year")
lines(monthlyYrs,wickMonthly-mean(wickMonthly, na.rm=TRUE),col='blue', lwd=2)> lines(monthlyYrs,aberdeenMonthly-mean(aberdeenMonthly, na.rm=TRUE),col='magenta', lwd=2)
lines(monthlyYrs,invergordonMonthly-mean(invergordonMonthly, na.rm=TRUE),col='green', lwd=2)
lines(monthlyYrs,stornowayMonthly-mean(stornowayMonthly, na.rm=TRUE),col='orange', lwd=2)
lines(monthlyYrs,torshavnMonthly-mean(torshavnMonthly, na.rm=TRUE),col='cyan', lwd=2)
lines(monthlyYrs,lerwickMonthly-mean(lerwickMonthly, na.rm=TRUE),col='red', lwd=2)
legend(x=as.Date("1958/1/1"), y=430, 
  legend=c('Lerwick','Aberdeen I','Invergordon','Stornoway','Wick','Torshavn'),
  col=c('red','magenta','green', 'orange', 'blue', 'cyan'), lwd=2, pch=21)
  
x11()
par(family='HersheySans')
plot(annualYrs,lerwickAnnual-mean(lerwickAnnual, na.rm=TRUE),
	type='l',col='red', lwd=2, ylim=c(-80,100), ann=F)
title(ylab="Sea Level [mm]", xlab="Year", main="Tide Gauges Around N Scotland")	
lines(annualYrs,wickAnnual-mean(wickAnnual, na.rm=TRUE),col='blue', lwd=2)
lines(annualYrs,aberdeenAnnual-mean(aberdeenAnnual, na.rm=TRUE),col='magenta', lwd=2)
lines(annualYrs,invergordonAnnual-mean(invergordonAnnual, na.rm=TRUE),col='green', lwd=2)
lines(annualYrs,stornowayAnnual-mean(stornowayAnnual, na.rm=TRUE),col='orange', lwd=2)
lines(annualYrs,torshavnAnnual-mean(torshavnAnnual, na.rm=TRUE),col='cyan', lwd=2)
lines(annualYrs,lerwickAnnual-mean(lerwickAnnual, na.rm=TRUE),col='red', lwd=2)
legend(x=as.Date("1965/1/1"), y=95, 
  legend=c('Lerwick','Aberdeen I','Invergordon','Stornoway','Wick','Torshavn'),
  col=c('red','magenta','green', 'orange', 'blue', 'cyan'), lwd=2, pch=21)
#
load("lerwickStnsSlpRates.RData")
# Calculate linear trends from the data
seqMonYrs <- c(1:length(monthlyYrs))
tgMonthlyRateArray <- vector(length=6, mode='numeric')
junk <- lm(lerwickMonthly ~ seqMonYrs)
tgMonthlyRateArray[1] <- junk$coef[2]
junk <- lm(wickMonthly ~ seqMonYrs)
tgMonthlyRateArray[2] <- junk$coef[2]
junk <- lm(invergordonMonthly ~ seqMonYrs)
tgMonthlyRateArray[3] <- junk$coef[2]
junk <- lm(stornowayMonthly ~ seqMonYrs)
tgMonthlyRateArray[4] <- junk$coef[2]
junk <- lm(aberdeenMonthly ~ seqMonYrs)
tgMonthlyRateArray[5] <- junk$coef[2]
junk <- lm(torshavnMonthly ~ seqMonYrs)
tgMonthlyRateArray[6] <- junk$coef[2]
# Correct for SLP
tgMonthlyRateArray <- tgMonthlyRateArray*12 + stnsSlpRates
#> tgMonthlyRateArray (uncorrected)
#[1] -0.320573  1.376047 -3.030858  2.044843  1.440282  1.925399
#> tgMonthlyRateArray (corrected)
#[1] -0.3216091  1.3753486 -3.0227314  2.0464070  1.4400261  1.9239595

postAltim <- which(monthlyYrs>=as.Date("1992/10/14"))
seqAltimYrs <- c(1:length(monthlyYrs[postAltim]))
tgAltMonthlyRateArray <- vector(length=6, mode='numeric')
junk <- lm(lerwickMonthly[postAltim] ~ seqAltimYrs)
tgAltMonthlyRateArray[1] <- junk$coef[2]
junk <- lm(wickMonthly[postAltim] ~ seqAltimYrs)
tgAltMonthlyRateArray[2] <- junk$coef[2]
junk <- lm(invergordonMonthly[postAltim] ~ seqAltimYrs)
tgAltMonthlyRateArray[3] <- junk$coef[2]
junk <- lm(stornowayMonthly[postAltim] ~ seqAltimYrs)
tgAltMonthlyRateArray[4] <- junk$coef[2]
junk <- lm(aberdeenMonthly[postAltim] ~ seqAltimYrs)
tgAltMonthlyRateArray[5] <- junk$coef[2]
junk <- lm(torshavnMonthly[postAltim] ~ seqAltimYrs)
tgAltMonthlyRateArray[6] <- junk$coef[2]
#> tgAltMonthlyRateArray*12
#[1] 3.697524 5.422544       NA 6.328178 5.911570 5.584658

# Correct for SLP
tgAltMonthlyRateArray <- tgAltMonthlyRateArray*12 + stnsSlpRates
#> tgAltMonthlyRateArray
#[1] 3.696488 5.421846       NA 6.329742 5.911314 5.583219

post1965 <- which(monthlyYrs>=as.Date("1965/1/1"))
seq1965Yrs <- c(1:length(monthlyYrs[post1965]))
tg1965MonthlyRateArray <- vector(length=6, mode='numeric')
junk <- lm(lerwickMonthly[post1965] ~ seq1965Yrs)
tg1965MonthlyRateArray[1] <- junk$coef[2]
junk <- lm(wickMonthly[post1965] ~ seq1965Yrs)
tg1965MonthlyRateArray[2] <- junk$coef[2]
junk <- lm(invergordonMonthly[post1965] ~ seq1965Yrs)
tg1965MonthlyRateArray[3] <- junk$coef[2]
junk <- lm(stornowayMonthly[post1965] ~ seq1965Yrs)
tg1965MonthlyRateArray[4] <- junk$coef[2]
junk <- lm(aberdeenMonthly[post1965] ~ seq1965Yrs)
tg1965MonthlyRateArray[5] <- junk$coef[2]
junk <- lm(torshavnMonthly[post1965] ~ seq1965Yrs)
tg1965MonthlyRateArray[6] <- junk$coef[2]
#> tg1965MonthlyRateArray*12
#[1] 0.3109622 1.3760466 2.3043650 2.0448427 1.9608385 2.0040454
# Correct for SLP
tg1965MonthlyRateArray <- tg1965MonthlyRateArray*12 + stnsSlpRates
#> tg1965MonthlyRateArray
#[1] 0.3099260 1.3753486 2.3124913 2.0464070 1.9605831 2.0026064

seqAnnYrs <- c(1:length(annualYrs))
tgAnnualRateArray <- vector(length=6, mode='numeric')
junk <- lm(lerwickAnnual ~ seqAnnYrs)
tgAnnualRateArray[1] <- junk$coef[2]
junk <- lm(wickAnnual ~ seqAnnYrs)
tgAnnualRateArray[2] <- junk$coef[2]
junk <- lm(invergordonAnnual ~ seqAnnYrs)
tgAnnualRateArray[3] <- junk$coef[2]
junk <- lm(stornowayAnnual ~ seqAnnYrs)
tgAnnualRateArray[4] <- junk$coef[2]
junk <- lm(aberdeenAnnual ~ seqAnnYrs)
tgAnnualRateArray[5] <- junk$coef[2]
junk <- lm(torshavnAnnual ~ seqAnnYrs)
tgAnnualRateArray[6] <- junk$coef[2]
#> tgAnnualRateArray
#[1] -0.5451691  1.3839711 -3.9230769  1.8337018  1.4446302  1.7885586
# Correct for SLP
tgAnnualRateArray <- tgAnnualRateArray + stnsSlpRates
#> tgAnnualRateArray
#[1] -0.5462052  1.3832732 -3.9149506  1.8352661  1.4443747  1.7871196

# GIA correction
pgrFile = '~/data/cabanes/accel/pgr.rates.stns.out.sed'
pgrDataFile <- read.table(pgrFile)
pgrCSCode <- pgrDataFile[,1]
# Extract relevant PGR rates (ICE-4G VM2)
stationCSCodes=vector(length=6, mode="numeric")
stnPGRrates=vector(length=6, mode="numeric")
for (i in 1:nStns){
  stationCSCodes[i] <- as.numeric(paste(stationCodes[1,i],'.',stationCodes[2,i], sep=""))
  junk <- which(as.numeric(pgrCSCode)==stationCSCodes[i])
  stnPGRrates[i] <- pgrDataFile[junk,2]
}
#> stnPGRrates - uses PLW's rate.f to extract numbers
#[1] 0.793 0.318 0.095 0.365 0.044 1.399
# From Peltier's book
#    0.47  -0.32 -0.65 -0.34 -0.60 1.34

#> tgAnnualRateArray-stnPGRrates
#[1] -1.3381691  1.0659711 -4.0180769  1.4687018  1.4006302  0.3895586
#> tgMonthlyRateArray*12-stnPGRrates
#[1] -1.1135730  1.0580466 -3.1258577  1.6798427  1.3962816  0.5263985
#> tg1965MonthlyRateArray*12-stnPGRrates
#[1] -0.4820378  1.0580466  2.2093650  1.6798427  1.9168385  0.6050454
#> tgAltMonthlyRateArray*12-stnPGRrates
#[1] 2.904524 5.104544       NA 5.963178 5.867570 4.185658

# Extract relevant PGR rates (ICE-5G VM2)
ice5gFile = '~/diskx/polcoms/lerwick/drsl.PSMSL.ICE5G_VM2_L90.mod'
ice5gDataFile <- read.table(ice5gFile)
ice5gCSCode <- as.numeric(ice5gDataFile[,3])/1000
stn5Grates=vector(length=6, mode="numeric")
for (i in 1:nStns){
  junk <- which(ice5gCSCode==stationCSCodes[i])
  stn5Grates[i] <- ice5gDataFile[junk,6]
}
#> tgAnnualRateArray-stn5Grates
#[1] -0.4251691  2.2639711 -2.8330769  2.5537018  2.3146302  0.4785586
#> tgMonthlyRateArray*12-stn5Grates
#[1] -0.2005730  2.2560466 -1.9408577  2.7648427  2.3102816  0.6153985

# Using numbers from Teferle et al 2006

# Lambeck
#> round(tgMonthlyRateArray[1],1)--1.8
#Lerwick: 1.5
#> round(tgMonthlyRateArray[5],1)-0
#Aberdeen: 1.4

# Peltier (book numbers)
#> round(tgMonthlyRateArray[1],1)--0.5
#Lerwick: 0.2
#> round(tgMonthlyRateArray[5],1)-0.6
#Aberdeen: 0.8

# CGPS/AG - CGPS doesn't seem right for Lerwick and there is no AG value for Aberdeen
# Assume AG is correct for Lerwick, and CGPS is correct for Aberdeen. Obtain proxy
# AG value for Aberdeen by assuming Newlyn CGPS-AG (=-0.5) is correct offset. Applying
# same offset to Aberdeen CGPS => Aberdeen AG value of +0.5 mm/yr
#> round(tgMonthlyRateArray[1],1)--1.1
#Lerwick: 0.8
#> round(tgMonthlyRateArray[5],1)-0.5
#Aberdeen: 0.9

# Using numbers from Sarah Bradley
#> round(tgMonthlyRateArray[1],1)--1.07
#Lerwick: 0.77
#> round(tgMonthlyRateArray[2],1)-0.13
#Wick: 1.27
#> round(tgMonthlyRateArray[3],1)-0.59
#Invergordon: -3.59
#> round(tgMonthlyRateArray[4],1)--0.08
#Stornoway: 2.08
#> round(tgMonthlyRateArray[5],1)-0.32
#Aberdeen: 1.08
#> round(tgMonthlyRateArray[6],1)--1.37
#Torshavn 3.27


nao<-read.table(file='~/data/NAO/nao.dat',na.string=-99.99)
#naoMonthlyYrs<-seq(from=(1821+0.5/12),to=(2000+11.5/12),by=1/12)
naoMonthlyYrs <- seq.Date(from=as.Date("1821/1/15"), 
 to=as.Date("2000/12/15", by="1 month")
naoIndex<-nao[,2:13]
naoIndex<-t(naoIndex)
dim(naoIndex)<-c(180*12,1)

j<-which(naoMonthlyYrs>as.Date("1957/1/1"))

x11()
par(family='HersheySans')
plot(naoMonthlyYrs[j],naoIndex[j],type='l',col='blue')
k<-which(nao[,1]>=1957)
lines(as.Date(paste(nao[k,1],"/1/1",sep='')),nao[k,14],col='red',lwd=2)
grid(col='black',lwd=1)

naoDJFM<-read.table(file='~/data/NAO/naodjfmindex.1864-2004.asc')
naoDJFMyrs <- seq.Date(from=as.Date("1864/7/1"), 
	to=as.Date("2004/7/1", by="1 year")
l<-which(naoDJFMyrs>=as.Date("1957/1/1"))
demeanLerwickMonthly<-lerwickMonthly-mean(lerwickMonthly,na.rm=TRUE)

x11()
par(family='HersheySans')
plot(naoDJFMyrs[l],naoDJFM[l,2]/abs(max(naoDJFM[l,2])),type='l',col='blue',lwd=2,
ylim=c(-1,1))
lines(monthlyYrs,demeanLerwickMonthly/abs(max(demeanLerwickMonthly,na.rm=TRUE)),
col='red')

demeanLerwickAnnual<-lerwickAnnual-mean(lerwickAnnual,na.rm=TRUE)
normDMLerwickAnnual<-demeanLerwickAnnual/max(abs(demeanLerwickAnnual),
na.rm=TRUE)

x11()
par(family='HersheySans')
plot(naoDJFM[l,1],naoDJFM[l,2]/abs(max(naoDJFM[l,2])),type='l',col='blue',lwd=2,
ylim=c(-1,1))
lines(annualYrs,demeanLerwickAnnual/max(abs(demeanLerwickAnnual),na.rm=TRUE),col
='red',lwd=2)

library(stats)
corNaoDJFMLerwickAnnual<-cor.test(demeanLerwickAnnual,naoDJFM[l,2])
corNaoDJFMLerwickAnnual

acf(demeanLerwickAnnual,na.action=na.pass)
acf(naoDJFM[l, 2],na.action=na.pass)

trian5ptFilter<-c(1,2,3,2,1)
trian5ptFilter<-trian5ptFilter/sum(trian5ptFilter)
filtDemeanLerwickAnnual<-filter(demeanLerwickAnnual,trian5ptFilter,method="c",  
  sides=2)
filtLerwickNaoDJFM<-filter(lerwickNaoDJFM,trian5ptFilter,method="c",sides=2)

x11()
par(family='HersheySans')
plot(filtDemeanLerwickAnnual/max(abs(filtDemeanLerwickAnnual),na.rm=TRUE),
  type='b',ylim=c(-1,1))
lines(filtLerwickNaoDJFM/max(abs(filtLerwickNaoDJFM),na.rm=TRUE),col='red')

corFiltNaoDJFMLerwickAnnual<-cor.test(filtDemeanLerwickAnnual,filtLerwickNaoDJFM ,alternative='less')

# POLCOMS data read with readLerwickModelTimeSeries.R
# Rename POLCOMS data
pcMonthlyMean<-lerwickMonthlyMean

# demean the Lerwick model data
dmPCMonthlyMean <- dmLerwickMonthlyMean
#rm(lerwickMonthlyMean,dmLerwickMonthlyMean)

s<-intersect(which(naoMonthlyYrs>1960),which(naoMonthlyYrs<1990))
u<-intersect(which(monthlyYrs>1960),which(monthlyYrs<1990))
acf(demeanLerwickMonthly[u],na.action=na.pass,lag.max=50)
acf(lerwickMonthly[u],na.action=na.pass,lag.max=75)
cor.test(lerwickMonthly[u],pcMonthlyMean)

normDMPCMonthlyMean<-pcMonthlyMean-mean(pcMonthlyMean)
normDMLerwickMonthly<-demeanLerwickMonthly/max(abs(demeanLerwickMonthly),
na.rm=TRUE)
normDMPCMonthlyMean<-dmPCMonthlyMean/max(abs(dmPCMonthlyMean))

x11()
par(family='HersheySans')
plot(naoMonthlyYrs[s],naoIndex[s]/max(abs(naoIndex[s]),na.rm=TRUE),type='l',col= 'blue',ylim=c(-1,1))
lines(naoMonthlyYrs[s],normDMPCMonthlyMean,col='red')
lines(monthlyYrs[u],normDMLerwickMonthly[u],col='magenta')

# Produce annual means of POLCOMS data
pcMonthlyMeanArray<-pcMonthlyMean
dim(pcMonthlyMeanArray)<-c(12,30)
pcAnnualMean<-array(NA,dim=c(30,1))
for (i in 1:30){
 pcAnnualMean[i]<-mean(pcMonthlyMeanArray[,i])
}
rm(pcMonthlyMeanArray)
# Produce DJFM means of POLCOMS data
pcMonthlyMeanArray<-pcMonthlyMean[12:(length(pcMonthlyMean)-1)]
dim(pcMonthlyMeanArray)<-c(12,29)
pcDJFMMean<-array(NA,dim=c(29,1))
for (i in 1:29){
 pcDJFMMean[i]<-mean(pcMonthlyMeanArray[1:4,i])
}
rm(pcMonthlyMeanArray)

dmPCAnnualMean<-pcAnnualMean-mean(pcAnnualMean)
dmPCDJFMMean<-pcDJFMMean-mean(pcDJFMMean)
filtDMPCAnnualMean<-filter(dmPCAnnualMean,trian5ptFilter,method="c",sides=2)
filtDMPCDJFMMean<-filter(dmPCDJFMMean,trian5ptFilter,method="c",sides=2)

x11()
par(family='HersheySans')
plot(annualYrs,filtDemeanLerwickAnnual/max(abs(filtDemeanLerwickAnnual),na.rm=TRUE),type='b',ylim=c(-1,1))
lines(annualYrs,filtLerwickNaoDJFM/max(abs(filtLerwickNaoDJFM),na.rm=TRUE),col='red')
lines(annualYrs[5:33],filtDMPCDJFMMean/max(abs(filtDMPCDJFMMean),na.rm=TRUE),col='magenta')
lines(annualYrs[4:33],filtDMPCAnnualMean/max(abs(filtDMPCAnnualMean),na.rm=TRUE),col='blue')

fiveYrRunMeanDMLerwickAnnual<-array(NA,c(43,1))
for (i in 1:43){
  fiveYrRunMeanDMLerwickAnnual[i]<-mean(demeanLerwickAnnual[i:(i+4)],na.rm=TRUE)
}

fiveYrRunMeanLerwickNaoDJFM<-array(NA,c(43,1))
for (i in 1:43){
  fiveYrRunMeanLerwickNaoDJFM[i]<-mean(lerwickNaoDJFM[i:(i+4)],na.rm=TRUE)
}

fiveYrRunMeanDMPCDJFMMean<-array(NA,c(25,1))
for (i in 1:25){
  fiveYrRunMeanDMPCDJFMMean[i]<-mean(dmPCDJFMMean[i:(i+4)],na.rm=TRUE)
}

fiveYrRunMeanDMPCAnnualMean<-array(NA,c(26,1))
for (i in 1:26){
  fiveYrRunMeanDMPCAnnualMean[i]<-mean(dmPCAnnualMean[i:(i+4)],na.rm=TRUE)
}

x11()
par(family='HersheySans')
plot(annualYrs[3:45],fiveYrRunMeanDMLerwickAnnual/max(abs(fiveYrRunMeanDMLerwickAnnual)),type='l',ylim=c(-1,1), ann=F)
lines(annualYrs[3:45],fiveYrRunMeanLerwickNaoDJFM/max(abs(fiveYrRunMeanLerwickNaoDJFM)),col='red')

lines(annualYrs[7:31],fiveYrRunMeanDMPCDJFMMean/max(abs(fiveYrRunMeanDMPCDJFMMean)),col='blue')
lines(annualYrs[6:31],fiveYrRunMeanDMPCAnnualMean/max(abs(fiveYrRunMeanDMPCAnnualMean)),col='magenta')

title(main="Five Year Running Means of Lerwick Data", xlab="Year", ylab="Normalised Values")
grid(col='gray4')
legend(1964,0.98,c("Lerwick TG","NAO DJFM", "POLCOMS Annual", "POLCOMS DJFM"),col=c("black","red","blue","magenta"),lty=c(1,1,1,1),bg="gray90")

corFiveYrDMLerwickNAODJFM <- cor.test(fiveYrRunMeanDMLerwickAnnual,fiveYrRunMeanLerwickNaoDJFM,alternative='less')
corFiveYrDMPCDJFMNAODJFM <- cor.test(fiveYrRunMeanDMPCDJFMMean,fiveYrRunMeanLerwickNaoDJFM[4:28],alternative='greater')
corFiveYrDMPCAnnualNAODJFM <- cor.test(fiveYrRunMeanDMPCAnnualMean,fiveYrRunMeanLerwickNaoDJFM[3:28],alternative='greater')
corFiveYrDMLerwickPCAnnual <- cor.test(fiveYrRunMeanDMLerwickAnnual[3:28],fiveYrRunMeanDMPCAnnualMean,alternative='less')
corFiveYrDMLerwickPCDJFM <- cor.test(fiveYrRunMeanDMLerwickAnnual[4:28],fiveYrRunMeanDMPCDJFMMean,alternative='less')

acf(fiveYrRunMeanDMLerwickAnnual)
acf(fiveYrRunMeanLerwickNaoDJFM)
acf(fiveYrRunMeanDMPCDJFMMean)
acf(fiveYrRunMeanDMPCAnnualMean)

annualMeans<-array(NA,dim=c(48,1))
for (i in 1957:2004){
  j <- intersect(which(monthlyYrs<(i+1)), which(monthlyYrs>=i))
  annualMeans[i-1956] <- mean(lerwickMonthly[j],na.rm=TRUE)
}

# Create NAO annual values
naoAnnual <- array(NA,dim=c(48,1))
for ( i in 1957:2004) {
  l <- intersect(which(naoMonthlyYrs>=i), which(naoMonthlyYrs<(i+1)))
  naoAnnual[i-1956] <- mean(naoIndex[l],na.rm=TRUE)
}
fiveYrRunMeanNAOAnnual<-array(NA,c(43,1))
for (i in 1:43){
  fiveYrRunMeanNAOAnnual[i]<-mean(naoAnnual[i:(i+4)],na.rm=TRUE)
}
cor.test(lerwickAnnual,naoAnnual)
cor.test(fiveYrRunMeanDMLerwickAnnual,fiveYrRunMeanNAOAnnual)

invergordon_ann <-
read.table(file='/home/simonh/diskx/polcoms/lerwick/invergordon_ann.txt',
  fill=TRUE,col.names=c("Yrs","SL","XX"))
wick_ann <-
read.table(file='/home/simonh/diskx/polcoms/lerwick/wick_ann.txt',
  fill=TRUE,col.names=c("Yrs","SL","XX"))
stornoway_ann <-
read.table(file='/home/simonh/diskx/polcoms/lerwick/stornoway_ann.txt',
  fill=TRUE,col.names=c("Yrs","SL","XX"))

dmInvergordonAnn<-invergordon_ann$SL-mean(invergordon_ann$SL,na.rm=TRUE)
dmWickAnn<-wick_ann$SL-mean(wick_ann$SL,na.rm=TRUE)
dmStornowayAnn<-stornoway_ann$SL-mean(stornoway_ann$SL,na.rm=TRUE)

normDMInvergordonAnn<-dmInvergordonAnn/max(abs(dmInvergordonAnn),na.rm=TRUE)
normDMWickAnn<-dmWickAnn/max(abs(dmWickAnn),na.rm=TRUE)
normDMStornowayAnn<-dmStornowayAnn/max(abs(dmStornowayAnn),na.rm=TRUE)

x11()
par(family='HersheySans')
plot(invergordon_ann$Yrs,normDMInvergordonAnn,type='l',col='blue',
xlim=c(1957,2004),ylim=c(-1,1))
lines(wick_ann$Yrs,normDMWickAnn,col='red')
lines(stornoway_ann$Yrs,normDMStornowayAnn,col='magenta')
lines(annualYrs,normDMLerwickAnnual)

normDMWickAnn40<-array(NA,dim=c(40,1))
for (i in 1:40){
  l<-which(wick_ann$Yrs==annualYrs[k[i]])
  if(length(l)>0){
    normDMWickAnn40[i]<-normDMWickAnn[l]
  }
}

cor.test(normDMLerwickAnnual[k],normDMWickAnn40)
#
#        Pearson's product-moment correlation
#
#data:  normDMLerwickAnnual[k] and normDMWickAnn40 
#t = 5.0017, df = 25, p-value = 3.711e-05
#alternative hypothesis: true correlation is not equal to 0 
#95 percent confidence interval:
# 0.4474731 0.8569341 
#sample estimates:
#      cor 
#0.7072261

cor.test(lerwickNaoDJFM[k],normDMWickAnn40)
#
#        Pearson's product-moment correlation
#
#data:  lerwickNaoDJFM[k] and normDMWickAnn40 
#t = 1.8637, df = 32, p-value = 0.07156
#alternative hypothesis: true correlation is not equal to 0 
#95 percent confidence interval:
# -0.02823931  0.58877761 
#sample estimates:
#      cor 
#0.3129146

cor.test(naoAnnual[k],normDMWickAnn40)
#
#        Pearson's product-moment correlation
#
#data:  naoAnnual[k] and normDMWickAnn40 
#t = 2.665, df = 29, p-value = 0.01245
#alternative hypothesis: true correlation is not equal to 0 
#95 percent confidence interval:
# 0.1058244 0.6895088 
# sample estimates:
#     cor 
#0.443533

stornoway_ann$Yrs[1]
#[1] 1977
k<-which(annualYrs>=stornoway_ann$Yrs[1])
setdiff(annualYrs[k],stornoway_ann$Yrs)
#[1] 1982 1983 1984 1985 1995 2004
length(annualYrs[k])
#[1] 28

normDMStornowayAnn28<-array(NA,dim=c(28,1))
for (i in 1:28){
  l<-which(stornoway_ann$Yrs==annualYrs[k[i]])
  if(length(l)>0){
    normDMStornowayAnn28[i]<-normDMStornowayAnn[l]
  }
}
cor.test(normDMLerwickAnnual[k],normDMStornowayAnn28)
#
#        Pearson's product-moment correlation
#
#data:  normDMLerwickAnnual[k] and normDMStornowayAnn28
#t = 2.3239, df = 13, p-value = 0.03698
#alternative hypothesis: true correlation is not equal to 0 
#95 percent confidence interval:
# 0.04082616 0.82505108 
#sample estimates:
#      cor 
#0.5417588

cor.test(lerwickNaoDJFM[k],normDMStornowayAnn28)
#
#        Pearson's product-moment correlation
#
#data:  lerwickNaoDJFM[k] and normDMStornowayAnn28 
#t = 0.8296, df = 20, p-value = 0.4166
#alternative hypothesis: true correlation is not equal to 0 
#95 percent confidence interval:
# -0.2591504  0.5608660 
#sample estimates:
#      cor 
#0.1823855

cor.test(naoAnnual[k],normDMStornowayAnn28)
#
#        Pearson's product-moment correlation
#
#data:  naoAnnual[k] and normDMStornowayAnn28
#t = 0.2992, df = 17, p-value = 0.7684
#alternative hypothesis: true correlation is not equal to 0
#95 percent confidence interval:
# -0.3948027  0.5098337
#sample estimates:
#       cor
#0.07238724
 
invergordon_ann$Yrs[1]
#[1] 1960
k<-which(annualYrs>=invergordon_ann$Yrs[1])
setdiff(annualYrs[k],invergordon_ann$Yrs)
#[1] 1982 1983 1984 1985 1995 2004
#[1] 1972 1973 1974 1975 1976 1977 1978 1979 1980 1981 1982 1983 1984 1985
#1986
#[16] 1987 1988 1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 1999 2000
#2001
#[31] 2002 2003 2004

length(annualYrs[k])
#[1] 45

normDMInvergordonAnn45<-array(NA,dim=c(45,1))
for (i in 1:45){
  l<-which(invergordon_ann$Yrs==annualYrs[k[i]])
  if(length(l)>0){
    normDMInvergordonAnn45[i]<-normDMInvergordonAnn[l]
  }
}
cor.test(normDMLerwickAnnual[k],normDMInvergordonAnn45)
#        Pearson's product-moment correlation
#
#data:  normDMLerwickAnnual[k] and normDMInvergordonAnn45 
#t = 3.4934, df = 10, p-value = 0.00579
#alternative hypothesis: true correlation is not equal to 0 
#95 percent confidence interval:
# 0.2914806 0.9226896 
#sample estimates:
#      cor 
#0.7413657

cor.test(lerwickNaoDJFM[k],normDMInvergordonAnn45)
#
#        Pearson's product-moment correlation
#
#data:  lerwickNaoDJFM[k] and normDMInvergordonAnn45 
#t = 1.1847, df = 10, p-value = 0.2635
#alternative hypothesis: true correlation is not equal to 0 
#95 percent confidence interval:
# -0.2793146  0.7697460 
#sample estimates:
#      cor 
#0.3508238

cor.test(naoAnnual[k],normDMInvergordonAnn45)
#
#        Pearson's product-moment correlation
#
#data:  naoAnnual[k] and normDMInvergordonAnn45 
#t = 1.0207, df = 10, p-value = 0.3314
#alternative hypothesis: true correlation is not equal to 0 
#95 percent confidence interval:
# -0.3238102  0.7490314 
#sample estimates:
#      cor 
#0.3071754

fiveYrRunMeanDMWickAnn<-array(NA,c(36,1))
for (i in 1:36){
  fiveYrRunMeanDMWickAnn[i]<-mean(normDMWickAnn40[i:(i+4)],na.rm=TRUE)
}

fiveYrRunMeanDMStornowayAnn<-array(NA,c(24,1))
for (i in 1:24){
  fiveYrRunMeanDMStornowayAnn[i]<-mean(normDMStornowayAnn28[i:(i+4)],na.rm=TRUE)
}

fiveYrRunMeanDMInvergordonAnn<-array(NA,c(41,1))
for (i in 1:41){
  fiveYrRunMeanDMInvergordonAnn[i]<-mean(normDMInvergordonAnn45[i:(i+4)],na.rm=TRUE)
}

