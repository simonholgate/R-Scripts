# Calculate EOFs over the N Atlantic to give some feel for the data
# Just focus on N Atlantic 100W-15E and 5S-80N

library(fields)

slp <- new.env()

slp$nlon<-72
slp$nlat<-37
slp$nyr<-160
slp$nmon<-12
slp$xlon<-seq(from=-180,by=5,length=slp$nlon)
slp$ylat<-seq(from=90,by=-5,length=slp$nlat)

slp$nAtlLon <- match(seq(from=-100,to=15,by=5), slp$xlon)
slp$nAtlLat <- match(seq(from=-5,to=80, by=5), slp$ylat)
slp$nNAtlLon <- length(slp$nAtlLon)
slp$nNAtlLat <- length(slp$nAtlLat)
                     
slp$slpMonthlyArray<-array(NA,dim=c(slp$nlon,slp$nlat,slp$nmon,slp$nyr))

slp$con<-file("~/diskx/HadSLP2r/hadSLP2_kij_1850-2009.bin","rb")

for (i in 1:slp$nlat){
  for (j in 1:slp$nlon){
    for (k in 1:slp$nyr){
      slp$slp<-readBin(slp$con,what="numeric", size=4, n=slp$nmon, endian='little')
      slp$slpMonthlyArray[j,i,,k]<-slp$slp
    }
  }
}

close(slp$con)

dim(slp$slpMonthlyArray)<-c(slp$nlon,slp$nlat,slp$nmon*slp$nyr)
slp$slpNAtlantic <- slp$slpMonthlyArray[slp$nAtlLon,slp$nAtlLat,]

slp$monthlyYears <-  seq(as.Date("1850/1/15"), by="month", length.out=1920)

# Map
library(maps)
library(mapdata)

x11()
map("worldHires", xlim=c(-100,15), ylim=c(-5,80), interior=F, fill=F, col="grey50", resolution=0)
map.axes()
contour(slp$xlon[slp$nAtlLon],slp$ylat[slp$nAtlLat],slp$slpNAtlantic[,,1], add=T)

# Now produce EOFs of the monthly data
slp$dmSlpMonthly <- slp$slpNAtlantic
slp$nstns <- (slp$nNAtlLon*slp$nNAtlLat)
dim(slp$dmSlpMonthly) <- c(slp$nstns,1920)

# Rotate so that there are nstns columns
slp$dmSlpMonthly <- t(slp$dmSlpMonthly)

# Remove column means
slp$cmeans <- colMeans(slp$dmSlpMonthly, na.rm=T)
for(i in 1:slp$nstns){
  slp$dmSlpMonthly[,i] <- slp$dmSlpMonthly[,i]-slp$cmeans[i]
}

# Form the co-variance matrix
slp$covDmSlpMonthly <- t(slp$dmSlpMonthly) %*% slp$dmSlpMonthly

# Calculate eigenvalues and eigenvectors
slp$eigCovDmSlpMonthly <- eigen(slp$covDmSlpMonthly)

# % variance explained
#slp$eigCovDmSlpMonthly$values/sum(slp$eigCovDmSlpMonthly$values)*100

# EOF1
slp$EOF1 <- slp$eigCovDmSlpMonthly$vectors[,1]
dim(slp$EOF1) <- c(slp$nNAtlLon, slp$nNAtlLat)

# EOF2
slp$EOF2 <- slp$eigCovDmSlpMonthly$vectors[,2]
dim(slp$EOF2) <- c(slp$nNAtlLon, slp$nNAtlLat)

# Find expansion coefficients
slp$C <- diag(x=slp$eigCovDmSlpMonthly$values,nrow=slp$nstns, ncol=slp$nstns)
#for (i in 1:nstns){
 slp$PC <- slp$dmSlpMonthly %*% slp$C[,1]

#quartz()
x11()
map("worldHires", xlim=c(-100,15), ylim=c(-5,80), interior=F, fill=F, col="grey50", resolution=0)
map.axes()
#EOF1
contour(slp$xlon[slp$nAtlLon],slp$ylat[slp$nAtlLat], slp$EOF1, lwd=2, add=T, col='blue')
#EOF2
contour(slp$xlon[slp$nAtlLon],slp$ylat[slp$nAtlLat], slp$EOF2, lty='dotdash', add=T, col='green')

# Compare with built-in R method
# Spatial loading pattern
slp$pcaSlpMonthly <- princomp(slp$dmSlpMonthly, scale=TRUE)

#PCA1
slp$PCA1 <- slp$pcaSlpMonthly$loadings[,1]
dim(slp$PCA1) <- c(slp$nNAtlLon, slp$nNAtlLat)

#PCA2
slp$PCA2 <- slp$pcaSlpMonthly$loadings[,2]
dim(slp$PCA2) <- c(slp$nNAtlLon, slp$nNAtlLat)

x11()
map("worldHires", xlim=c(-100,15), ylim=c(-5,80), interior=F, fill=F, col="grey50", resolution=0)
map.axes()
#PCA1
contour(slp$xlon[slp$nAtlLon],slp$ylat[slp$nAtlLat], slp$PCA1, lwd=2, add=T, col='red')
#PCA2
contour(slp$xlon[slp$nAtlLon],slp$ylat[slp$nAtlLat], slp$PCA2, lty='dotdash', add=T, col='blue')


x11()
map("worldHires", xlim=c(-25,15), ylim=c(35,65), interior=F, fill=F, col="grey50", resolution=0)
map.axes()
#PCA1
contour(slp$xlon[slp$nAtlLon],slp$ylat[slp$nAtlLat], slp$PCA1, lwd=2, add=T, col='red')
#PCA2
contour(slp$xlon[slp$nAtlLon],slp$ylat[slp$nAtlLat], slp$PCA2, lty='dotdash', add=T, col='blue')

# Time series pattern
#quartz()
x11()
plot(slp$monthlyYears,slp$pcaSlpMonthly$scores[,1], type='l', col='red')
lines(slp$monthlyYears,slp$pcaSlpMonthly$scores[,2], col='blue')
lines(slp$monthlyYears,slp$pcaSlpMonthly$scores[,3], col='magenta')
lines(slp$monthlyYears,(slp$PC/50000), col='cyan', lty='dotted')

# Annual data
tmp <- new.env()
slp$hadSLP2rAnnualPressureArray <- array(NA, dim=c(160,432))
k <- 1
for(i in 1:slp$nNAtlLat){
  for(j in 1:slp$nNAtlLon){
    tmp$hadSLP2rMonthlyPressureArray <- slp$slpNAtlantic[j,i,]
    dim(tmp$hadSLP2rMonthlyPressureArray) <- c(12,160)
    slp$hadSLP2rAnnualPressureArray[,k] <- colMeans(tmp$hadSLP2rMonthlyPressureArray, na.rm=T)
    k <- k+1
  }
}

slp$annualYears <- seq(as.Date("1850/6/15"), as.Date("2009/6/15"), "years")
rm(tmp)

# Spatial loading pattern
#slp$pcaSlpAnnual <- princomp(slp$hadSLP2rAnnualPressureArray, scale=TRUE)

#summary(slp$pcaSlpAnnual)

# EOF1
#slp$PCAann1 <- slp$pcaSlpAnnual$loadings[,1]
#dim(slp$PCAann1) <- c(24,18)

# EOF2
#slp$PCAann2 <- slp$pcaSlpAnnual$loadings[,2]
#dim(slp$PCAann2) <- c(24,18)

# EOF3
#slp$PCAann3 <- slp$pcaSlpAnnual$loadings[,3]
#dim(slp$PCAann3) <- c(24,18)

# Remove column means
slp$dmSlpAnnual <- slp$hadSLP2rAnnualPressureArray
slp$cmeans <- colMeans(slp$dmSlpAnnual, na.rm=T)
for(i in 1:slp$nstns){
  slp$dmSlpAnnual[,i] <- slp$dmSlpAnnual[,i]-slp$cmeans[i]
}

# SVD method
slp$svdDmSlpAnnual <- svd(slp$dmSlpAnnual)

# Amount of variance explained is in d^2
slp$svdVarExp <- slp$svdDmSlpAnnual$d^2/sum(slp$svdDmSlpAnnual$d^2)*100

#PCA1
slp$PCAann1 <- slp$svdDmSlpAnnual$v[,1]
dim(slp$PCAann1) <- c(slp$nNAtlLon, slp$nNAtlLat)

#PCA2
slp$PCAann2 <- slp$svdDmSlpAnnual$v[,2]
dim(slp$PCAann2) <- c(slp$nNAtlLon, slp$nNAtlLat)

#PCA3
slp$PCAann3 <- slp$svdDmSlpAnnual$v[,3]
dim(slp$PCAann3) <- c(slp$nNAtlLon, slp$nNAtlLat)

# Find expansion coefficients
slp$ECann1 <- slp$dmSlpAnnual %*% slp$svdDmSlpAnnual$v[,1]
slp$ECann2 <- slp$dmSlpAnnual %*% slp$svdDmSlpAnnual$v[,2]
slp$ECann3 <- slp$dmSlpAnnual %*% slp$svdDmSlpAnnual$v[,3]

# Map
#quartz()
x11()
map("worldHires", xlim=c(-100,15), ylim=c(-5,80), interior=F, fill=F, col="grey50", resolution=0)
map.axes()
#PCA1
contour(slp$xlon[slp$nAtlLon],slp$ylat[slp$nAtlLat], slp$PCAann1, lwd=2, add=T, col='red')
#PCA2
contour(slp$xlon[slp$nAtlLon],slp$ylat[slp$nAtlLat], slp$PCAann2, lty='dotdash', add=T, col='blue')
#PCA3
contour(slp$xlon[slp$nAtlLon],slp$ylat[slp$nAtlLat], slp$PCAann3, lty='dotted', add=T, col='green')

# Time series pattern
#quartz()
x11()
plot(slp$annualYears,slp$ECann1, type='l', col='red')
lines(slp$annualYears,slp$ECann2, col='blue')
lines(slp$annualYears,slp$ECann3, col='magenta')


# Eigenvector method
# Form the co-variance matrix
slp$covDmSlpAnnual <- t(slp$dmSlpAnnual) %*% slp$dmSlpAnnual

# Calculate eigenvalues and eigenvectors
slp$eigCovDmSlpAnnual <- eigen(slp$covDmSlpAnnual)

# % variance explained
slp$eigCovDmSlpAnnual$values/sum(slp$eigCovDmSlpAnnual$values)*100

# EOF1
slp$EOFann1 <- slp$eigCovDmSlpAnnual$vectors[,1]
dim(slp$EOFann1) <- c(slp$nNAtlLon, slp$nNAtlLat)

# EOF2
slp$EOFann2 <- slp$eigCovDmSlpAnnual$vectors[,2]
dim(slp$EOFann2) <- c(slp$nNAtlLon, slp$nNAtlLat)

# EOF3
slp$EOFann3 <- slp$eigCovDmSlpAnnual$vectors[,3]
dim(slp$EOFann3) <- c(slp$nNAtlLon, slp$nNAtlLat)

# Find expansion coefficients
slp$Cann <- diag(x=slp$eigCovDmSlpAnnual$values,nrow=slp$nstns, ncol=slp$nstns)
#for (i in 1:nstns){
slp$PCann1 <- slp$dmSlpAnnual %*% slp$Cann[,1]
slp$PCann2 <- slp$dmSlpAnnual %*% slp$Cann[,2]
slp$PCann3 <- slp$dmSlpAnnual %*% slp$Cann[,3]

# Map
#quartz()
x11()
map("worldHires", xlim=c(-100,15), ylim=c(-5,80), interior=F, fill=F, col="grey50", resolution=0)
map.axes()
#PCA1
contour(slp$xlon[slp$nAtlLon],slp$ylat[slp$nAtlLat], slp$EOFann1, lwd=2, add=T, col='red')
#PCA2
contour(slp$xlon[slp$nAtlLon],slp$ylat[slp$nAtlLat], slp$EOFann2, lty='dotdash', add=T, col='blue')
#PCA3
contour(slp$xlon[slp$nAtlLon],slp$ylat[slp$nAtlLat], slp$EOFann3, lty='dotted', add=T, col='green')

# Time series pattern
#quartz()
x11()
plot(slp$annualYears,slp$PCann1, type='l', col='red')
lines(slp$annualYears,slp$PCann2, col='blue')
lines(slp$annualYears,slp$PCann3, col='magenta')


# Save the data
#save(slp$slpNAtlantic, file="nAtlHadSLP2.RData")

