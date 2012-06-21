###############################################################################
# Read the monthly mean temperatures, salinities, two components of
# velocity and two components of stress calculated from POLCOMS model output. 
# 
# The arrays were written from FORTRAN with the readmonthmean_TSUV_series
# program and are big endian real*4 numbers. 
#
# Simon Holgate, July 2009
#
###############################################################################

#######################
# _*S2 dimensions:    #
#                     #
# *_Xmin:   -19.83333 #
# Xres:   1/6         #
# Xn:     198         #
#                     #
# Ymin:   40.11111    #
# Yres:   1/9         #
# Yn:     224         #
#######################

# Newlyn 50 06 N  05 33 W => 50.1 -5.55
# Newlyn is at x=86, y=90
#> (-5.55-xmin)/xres
#[1] 85.69998
#> (50.1-ymin)/yres
#[1] 89.90001
# Brest 48 23 N  04 30 W
# Brest is at x=92, y=74
#> (-4.5-min)/xres
#[1] 91.99998
#> (48.38-ymin)/yres
#[1] 74.42001

# Also choose 6 boundary points to reflect the deep ocean boundary condition
# Choose points at:
# x=1, y=50, x=1, y=100, x=1, y=150
# x=50, y=50, x=50, y=100, x=50, y=150,

load("~/diskx/polcoms/iseajseanpsea.Rdata")
load("~/diskx/polcoms/iuseajuseanpusea.RData")

library(fields)
# Use jet colors for images
source('~/bin/RScripts/jet.colors.R')
lat<-seq(from=40.11111,by=1/9,length=224)
lon<-seq(from=-19.83333,by=1/6,length=198)

l<-198
m<-224
lm<-l*m
n<-40
nmonths<-n*12

nseries <- 8

yearsArray <- c(1960,1961,1962,1963,1964,1965,1966,1967,1968,1969,
    1970,1971,1972,1973,1974,1975,1976,1977,1978,1979,
    1980,1981,1982,1983,1984,1985,1986,1987,1988,1989,
    1990,1991,1992,1993,1994,1995,1996,1997,1998,1999)

time <- seq.Date(from=as.Date("1960/1/15"), to=as.Date("1999/12/15"), by="1 month")

monthlyMeanV<-array(NA,dim=c(nseries,nmonths,3))
monthlyMeanU<-array(NA,dim=c(nseries,nmonths,3))
monthlyMeanTau<-array(NA,dim=c(nseries,nmonths,4))
monthlyMeanT<-array(NA,dim=c(nseries,nmonths,3))
monthlyMeanS<-array(NA,dim=c(nseries,nmonths,3))

monthCount<-1

#
# Barotropic flow
#
monthlyMeanUbVb<-array(NA,dim=c(l,m,2))
inMonthlyMeanUbVb <- 
  file("~/diskx/polcoms/S12run405UVTSFiles/UbVbmonthMeanArray.S12run405.dat", "rb")
junk1<-readBin(inMonthlyMeanUbVb, what="integer", n = 1, size = NA, endian = "big")
ubvb<-readBin(inMonthlyMeanUbVb,n=npusea, what='numeric', size=4, endian='big')
for (ip in 1:npusea) {
  monthlyMeanUbVb[iusea[ip],jusea[ip],1]<-ubvb[ip]
}
junk1<-readBin(inMonthlyMeanUbVb, what="integer", n = 1, size = NA, endian = "big")
ubvb<-readBin(inMonthlyMeanUbVb,n=npusea, what='numeric', size=4, endian='big')
for (ip in 1:npusea) {
  monthlyMeanUbVb[iusea[ip],jusea[ip],2]<-ubvb[ip]
}
close(inMonthlyMeanUbVb)
x11()
image.plot(lon,lat,monthlyMeanUbVb[,,1])
x11()
image.plot(lon,lat,monthlyMeanUbVb[,,2])

# Magnitude of the barotropic flow
magUbVb <- sqrt(monthlyMeanUbVb[,,1]^2 + monthlyMeanUbVb[,,2]^2)
x11()
image.plot(lon,lat,magUbVb,zlim=c(0,0.2))

# Phase of the barotropic flow
phaseUbVb <- array(NA,dim=c(l,m))
for (ip in 1:npusea) {
  vb <- monthlyMeanUbVb[iusea[ip],jusea[ip],2]
  ub <- monthlyMeanUbVb[iusea[ip],jusea[ip],1]
  if ((ub>0)&&(vb>0)){
    phaseUbVb[iusea[ip],jusea[ip]] <- atan(abs(vb)/abs(ub))*180/pi
  }
  if ((ub<0)&&(vb>0)){
	  phaseUbVb[iusea[ip],jusea[ip]] <- 180-atan(abs(vb)/abs(ub))*180/pi
  }
  if ((ub>0)&&(vb<0)){
	  phaseUbVb[iusea[ip],jusea[ip]] <- 360-atan(abs(vb)/abs(ub))*180/pi
  }
  if ((ub<0)&&(vb<0)){
	  phaseUbVb[iusea[ip],jusea[ip]] <- 180+atan(abs(vb)/abs(ub))*180/pi
  }
}
x11()
image.plot(lon,lat,phaseUbVb)

#
# Surface flow
#
monthlyMeanUmVm<-array(NA,dim=c(l,m,2))

inMonthlyMeanUmVm <- 
		file("~/diskx/polcoms/S12run405UVTSFiles/UmVmmonthMeanArray.S12run405.dat", "rb")
junk1<-readBin(inMonthlyMeanUmVm, what="integer", n = 1, size = NA, endian = "big")
umvm<-readBin(inMonthlyMeanUmVm,n=npusea, what='numeric', size=4, endian='big')
for (ip in 1:npusea) {
	monthlyMeanUmVm[iusea[ip],jusea[ip],1]<-ubvb[ip]
}
junk1<-readBin(inMonthlyMeanUmVm, what="integer", n = 1, size = NA, endian = "big")
umvm<-readBin(inMonthlyMeanUmVm,n=npusea, what='numeric', size=4, endian='big')
for (ip in 1:npusea) {
	monthlyMeanUmVm[iusea[ip],jusea[ip],2]<-umvm[ip]
}
close(inMonthlyMeanUmVm)

x11()
image.plot(lon[75:100],lat[50:100],monthlyMeanUmVm[75:100,50:100,1])
x11()
image.plot(lon[75:100],lat[50:100],monthlyMeanUmVm[75:100,50:100,2])

# Magnitude of the surface flow
magUmVm <- sqrt(monthlyMeanUmVm[,,1]^2 + monthlyMeanUmVm[,,2]^2)
x11()
image.plot(lon[75:100],lat[50:100],magUmVm[75:100,50:100],zlim=c(0,0.2))

# Phase of the surface flow
phaseUmVm <- array(NA,dim=c(l,m))
for (ip in 1:npusea) {
	vm <- monthlyMeanUmVm[iusea[ip],jusea[ip],2]
	um <- monthlyMeanUmVm[iusea[ip],jusea[ip],1]
	if ((um>0)&&(vm>0)){
		phaseUmVm[iusea[ip],jusea[ip]] <- atan(abs(vm)/abs(um))*180/pi
	}
	if ((um<0)&&(vm>0)){
		phaseUmVm[iusea[ip],jusea[ip]] <- 180-atan(abs(vm)/abs(um))*180/pi
	}
	if ((um>0)&&(vm<0)){
		phaseUmVm[iusea[ip],jusea[ip]] <- 360-atan(abs(vm)/abs(um))*180/pi
	}
	if ((um<0)&&(vm<0)){
		phaseUmVm[iusea[ip],jusea[ip]] <- 180+atan(abs(vm)/abs(um))*180/pi
	}
}
x11()
image.plot(lon[75:100],lat[50:100],phaseUmVm[75:100,50:100])

#
# Temperature
#

inMonthlyMeanTb <-
  file("~/diskx/polcoms/S12run405UVTSFiles/TbmonthMeanArray.S12run405.dat", "rb")
junk1<-readBin(inMonthlyMeanTb, what="integer", n = 1, size = NA, endian = "big")
Tb<-readBin(inMonthlyMeanTb,n=npsea, what='numeric', size=8, endian='big')
close(inMonthlyMeanTb)
for (ip in 1:npsea) {
  monthlyMeanTb[isea[ip],jsea[ip]]<-Tb[ip]
}
x11()
image.plot(lon,lat,monthlyMeanTb)

mask<-array(NA,dim=c(l,m))
mask[86,90] <- 0
mask[92,74] <- 0
mask[1,50] <- 0
mask[1,100] <- 0
mask[1,150] <- 0
mask[50,150] <- 0
mask[50,100] <- 0
mask[50,50] <- 0
image.plot(lon,lat,mask, add=T, col=grey(c(0,1)))
#
# Time series
#

# i,j indices of the nseries points to be extracted e.g.
# 86 90
# 92 74
#  2 50
#  2 100
#  2 150
# 50  50
# 50 100
# 50 150

#c Tmonthmean.S12run405.dat
#write(66) (tmpb(ipt(is),jpt(is)), is=1,nseries)
#c          write(66) ((tmpm(isea(i),jsea(i),k),i=1,npsea),k=1,n-2)
#write(66) (tmpm(ipt(is),jpt(is),n-2), is=1,nseries)
#write(66) (tmpm(ipt(is),jpt(is),1), is=1,nseries)
#
# monthlyMeanT<-array(NA,dim=c(nseries,12*n,3))
inFile <-
  file("~/diskx/polcoms/S12run405UVTSFiles/fullSeries/Tmonthmean.S12run405.dat", 
    "rb")
for (i in 1:nmonths){
  junk1 <- readBin(inFile, what="integer", n = 1, size = NA, endian = "big")
  data <- readBin(inFile,n=nseries, what='numeric', size=8, endian='big')
  monthlyMeanT[,i,1] <- data
  junk1 <- readBin(inFile, what="integer", n = 1, size = NA, endian = "big")

  junk1 <- readBin(inFile, what="integer", n = 1, size = NA, endian = "big")
  data <- readBin(inFile,n=nseries, what='numeric', size=4, endian='big')
  monthlyMeanT[,i,2] <- data
  junk1 <- readBin(inFile, what="integer", n = 1, size = NA, endian = "big")
  
  junk1 <- readBin(inFile, what="integer", n = 1, size = NA, endian = "big")
  data <- readBin(inFile,n=nseries, what='numeric', size=4, endian='big')
  monthlyMeanT[,i,3] <- data
  junk1 <- readBin(inFile, what="integer", n = 1, size = NA, endian = "big")
}
close(inFile)

x11()
plot(time, monthlyMeanT[1,,2], type='l', col='blue')
lines(time, monthlyMeanT[1,,3], col='magenta')
lines(time, monthlyMeanT[1,,1], col='red')

#c Smonthmean.S12run405.dat
#write(67) (salb(ipt(is),jpt(is)), is=1,nseries)
#c          write(67) ((salm(isea(i),jsea(i),k),i=1,npsea),k=1,n-2)
#write(67) (salm(ipt(is),jpt(is),n-2), is=1,nseries)
#write(67) (salm(ipt(is),jpt(is),1), is=1,nseries)
#
# monthlyMeanS<-array(NA,dim=c(nseries,12*n,3))
inFile <-
		file("~/diskx/polcoms/S12run405UVTSFiles/fullSeries/Smonthmean.S12run405.dat", 
				"rb")
for (i in 1:nmonths){
	junk1 <- readBin(inFile, what="integer", n = 1, size = NA, endian = "big")
	data <- readBin(inFile,n=nseries, what='numeric', size=8, endian='big')
	monthlyMeanS[,i,1] <- data
	junk1 <- readBin(inFile, what="integer", n = 1, size = NA, endian = "big")
	
	junk1 <- readBin(inFile, what="integer", n = 1, size = NA, endian = "big")
	data <- readBin(inFile,n=nseries, what='numeric', size=4, endian='big')
	monthlyMeanS[,i,2] <- data
	junk1 <- readBin(inFile, what="integer", n = 1, size = NA, endian = "big")
	
	junk1 <- readBin(inFile, what="integer", n = 1, size = NA, endian = "big")
	data <- readBin(inFile,n=nseries, what='numeric', size=4, endian='big')
	monthlyMeanS[,i,3] <- data
	junk1 <- readBin(inFile, what="integer", n = 1, size = NA, endian = "big")
}
close(inFile)

x11()
plot(time, monthlyMeanS[1,,2], type='l', col='blue')
lines(time, monthlyMeanS[1,,3], col='magenta')
lines(time, monthlyMeanS[1,,1], col='red')

#c UVmonthmean.S12run405.dat
#write(167) (ub(ipt(is),jpt(is)), is=1,nseries)
#c          write(167) ((um(iusea(i),jusea(i),k),i=1,npusea),k=1,n-2)
#write(167) (um(ipt(is),jpt(is),n-2), is=1,nseries)
#write(167) (um(ipt(is),jpt(is),1), is=1,nseries)
#c           write(6,*) 'w',tau_u(110,185)
#write(167) (vb(ipt(is),jpt(is)), is=1,nseries)
#c          write(167) ((vm(iusea(i),jusea(i),k),i=1,npusea),k=1,n-2)
#write(167) (vm(ipt(is),jpt(is),n-2), is=1,nseries)
#write(167) (vm(ipt(is),jpt(is),1), is=1,nseries)
#monthlyMeanV<-array(NA,dim=c(nseries,nmonths,3))
#monthlyMeanU<-array(NA,dim=c(nseries,nmonths,3))

inFile <-
		file("~/diskx/polcoms/S12run405UVTSFiles/fullSeries/UVmonthmean.S12run405.dat", 
				"rb")
for (i in 1:nmonths){
	junk1 <- readBin(inFile, what="integer", n = 1, size = NA, endian = "big")
	data <- readBin(inFile,n=nseries, what='numeric', size=4, endian='big')
	monthlyMeanU[,i,1] <- data
	junk1 <- readBin(inFile, what="integer", n = 1, size = NA, endian = "big")
	
	junk1 <- readBin(inFile, what="integer", n = 1, size = NA, endian = "big")
	data <- readBin(inFile,n=nseries, what='numeric', size=4, endian='big')
	monthlyMeanU[,i,2] <- data
	junk1 <- readBin(inFile, what="integer", n = 1, size = NA, endian = "big")
	
	junk1 <- readBin(inFile, what="integer", n = 1, size = NA, endian = "big")
	data <- readBin(inFile,n=nseries, what='numeric', size=4, endian='big')
	monthlyMeanU[,i,3] <- data
	junk1 <- readBin(inFile, what="integer", n = 1, size = NA, endian = "big")
	
	junk1 <- readBin(inFile, what="integer", n = 1, size = NA, endian = "big")
	data <- readBin(inFile,n=nseries, what='numeric', size=4, endian='big')
	monthlyMeanV[,i,1] <- data
	junk1 <- readBin(inFile, what="integer", n = 1, size = NA, endian = "big")
	
	junk1 <- readBin(inFile, what="integer", n = 1, size = NA, endian = "big")
	data <- readBin(inFile,n=nseries, what='numeric', size=4, endian='big')
	monthlyMeanV[,i,2] <- data
	junk1 <- readBin(inFile, what="integer", n = 1, size = NA, endian = "big")
	
	junk1 <- readBin(inFile, what="integer", n = 1, size = NA, endian = "big")
	data <- readBin(inFile,n=nseries, what='numeric', size=4, endian='big')
	monthlyMeanV[,i,3] <- data
	junk1 <- readBin(inFile, what="integer", n = 1, size = NA, endian = "big")
}
close(inFile)

x11()
plot(time, monthlyMeanU[1,,2], type='l', col='blue')
lines(time, monthlyMeanU[1,,3], col='magenta')
lines(time, monthlyMeanU[1,,1], col='red')

x11()
plot(time, monthlyMeanV[1,,2], type='l', col='blue')
lines(time, monthlyMeanV[1,,3], col='magenta')
lines(time, monthlyMeanV[1,,1], col='red')

# Magnitude timeseries

monthlyMeanMagUV <- sqrt(monthlyMeanU[,,]^2 + monthlyMeanV[,,]^2)
x11()
plot(time, monthlyMeanMagUV[3,,2], type='l', col='blue')
lines(time, monthlyMeanMagUV[3,,1], col='red')
lines(time, monthlyMeanMagUV[3,,3], col='magenta')

# Phase timeseries

monthlyMeanPhaseUV<-array(NA,dim=c(nseries,nmonths,3))

for (i in 1:nseries) {
  for (j in 1:nmonths){
	for (k in 1:3){
	  vm <- monthlyMeanV[i,j,k]
	  um <- monthlyMeanU[i,j,k]
	  if (is.nan(um)||is.nan(vm)){
		  monthlyMeanPhaseUV[i,j,k] <- NA
	  } else {
	    if ((um>0)&&(vm>0)){
		  monthlyMeanPhaseUV[i,j,k] <- atan(abs(vm)/abs(um))*180/pi
	    }
	    if ((um<0)&&(vm>0)){
		  monthlyMeanPhaseUV[i,j,k] <- 180-atan(abs(vm)/abs(um))*180/pi
	    }
	    if ((um>0)&&(vm<0)){
		  monthlyMeanPhaseUV[i,j,k] <- 360-atan(abs(vm)/abs(um))*180/pi
	    }
	    if ((um<0)&&(vm<0)){
		  monthlyMeanPhaseUV[i,j,k] <- 180+atan(abs(vm)/abs(um))*180/pi
	    }
      }
    }
  }
}
x11()
plot(time, monthlyMeanPhaseUV[1,,2], type='l', col='blue')
lines(time, monthlyMeanPhaseUV[1,,1], col='red')
lines(time, monthlyMeanPhaseUV[1,,3], col='magenta')

#c FBGBmonthmean.S12run405.dat
#write(168) (tau_u(ipt(is),jpt(is)), is=1,nseries)
#write(168) (tau_v(ipt(is),jpt(is)), is=1,nseries)
#write(168) (tau_ub(ipt(is),jpt(is)), is=1,nseries)
#write(168) (tau_vb(ipt(is),jpt(is)), is=1,nseries)
# monthlyMeanTau<-array(NA,dim=c(nseries,nmonths,4))
inFile <-
		file("~/diskx/polcoms/S12run405UVTSFiles/fullSeries/FBGBmonthmean.S12run405.dat", 
				"rb")
for (i in 1:nmonths){
	junk1 <- readBin(inFile, what="integer", n = 1, size = NA, endian = "big")
	data <- readBin(inFile,n=nseries, what='numeric', size=8, endian='big')
	monthlyMeanTau[,i,1] <- data
	junk1 <- readBin(inFile, what="integer", n = 1, size = NA, endian = "big")
	
	junk1 <- readBin(inFile, what="integer", n = 1, size = NA, endian = "big")
	data <- readBin(inFile,n=nseries, what='numeric', size=8, endian='big')
	monthlyMeanTau[,i,2] <- data
	junk1 <- readBin(inFile, what="integer", n = 1, size = NA, endian = "big")
	
	junk1 <- readBin(inFile, what="integer", n = 1, size = NA, endian = "big")
	data <- readBin(inFile,n=nseries, what='numeric', size=8, endian='big')
	monthlyMeanTau[,i,3] <- data
	junk1 <- readBin(inFile, what="integer", n = 1, size = NA, endian = "big")
	
	junk1 <- readBin(inFile, what="integer", n = 1, size = NA, endian = "big")
	data <- readBin(inFile,n=nseries, what='numeric', size=8, endian='big')
	monthlyMeanTau[,i,4] <- data
	junk1 <- readBin(inFile, what="integer", n = 1, size = NA, endian = "big")
}
close(inFile)

x11()
plot(time, monthlyMeanTau[1,,2], type='l', col='blue')
lines(time, monthlyMeanTau[1,,3], col='magenta')
lines(time, monthlyMeanTau[1,,1], col='red')
lines(time, monthlyMeanTau[1,,4], col='cyan')

# Magnitude timeseries

monthlyMeanMagTau <- sqrt(monthlyMeanTau[,,1]^2 + monthlyMeanTau[,,2]^2)
monthlyMeanMagTauB <- sqrt(monthlyMeanTau[,,3]^2 + monthlyMeanTau[,,4]^2)
x11()
plot(time, monthlyMeanMagTau[1,], type='l', col='blue')
lines(time, monthlyMeanMagTauB[1,], col='red')

# Phase timeseries

monthlyMeanPhaseTau<-array(NA,dim=c(nseries,nmonths,2))

for (i in 1:nseries) {
	for (j in 1:nmonths){

		vm <- monthlyMeanTau[i,j,2]
		um <- monthlyMeanTau[i,j,1]
		if (is.nan(um)||is.nan(vm)){
			monthlyMeanPhaseTau[i,j,1] <- NA
		} else {
			if ((um>0)&&(vm>0)){
				monthlyMeanPhaseTau[i,j,1] <- atan(abs(vm)/abs(um))*180/pi
			}
			if ((um<0)&&(vm>0)){
				monthlyMeanPhaseTau[i,j,1] <- 180-atan(abs(vm)/abs(um))*180/pi
			}
			if ((um>0)&&(vm<0)){
				monthlyMeanPhaseTau[i,j,1] <- 360-atan(abs(vm)/abs(um))*180/pi
			}
			if ((um<0)&&(vm<0)){
				monthlyMeanPhaseTau[i,j,1] <- 180+atan(abs(vm)/abs(um))*180/pi
			}
		}
		vm <- monthlyMeanTau[i,j,4]
		um <- monthlyMeanTau[i,j,3]
		if (is.nan(um)||is.nan(vm)){
			monthlyMeanPhaseTau[i,j,2] <- NA
		} else {
			if ((um>0)&&(vm>0)){
				monthlyMeanPhaseTau[i,j,2] <- atan(abs(vm)/abs(um))*180/pi
			}
			if ((um<0)&&(vm>0)){
				monthlyMeanPhaseTau[i,j,2] <- 180-atan(abs(vm)/abs(um))*180/pi
			}
			if ((um>0)&&(vm<0)){
				monthlyMeanPhaseTau[i,j,2] <- 360-atan(abs(vm)/abs(um))*180/pi
			}
			if ((um<0)&&(vm<0)){
				monthlyMeanPhaseTau[i,j,2] <- 180+atan(abs(vm)/abs(um))*180/pi
			}
		}
		
	}
}
x11()
plot(time, monthlyMeanPhaseTau[1,,1], type='l', col='blue')
lines(time, monthlyMeanPhaseTau[1,,2], col='red')


