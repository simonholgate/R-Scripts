# R-Script to read monthly mean sea level data calculated from CS3X files
# and to extract the timeseries at Brest and at Newlyn
# 
# Author: simonh
###############################################################################

#
# File format: Each file is named e.g. 'tid.0000Z01011958.0000Z01011959.mm'
# The '.mm' indicates that the file contains monthly means. 
# The files are unfomratted binary written from FORTRAN. Each contains 12
# months of full arrays (not compact arrays as CS3X outputs) of sea level,
# u velocity and v velocity, in tha order i.e. Jan: sl, uv, vv, Feb: sl, 
# uv, vv etc. The arrays are of size (num_lon, num_lat). 
# The fastest moving index is j (latitude).
#

# Newlyn 50 06 N  05 33 W => 50.1 -5.55
# Newlyn is at x=88, y=91
#> (-5.55-START_LON)/x_res
#[1] 86.7
#> (50.1-START_LAT)/y_res
#[1] 90.0
# Brest 48 23 N  04 30 W
# Brest is at x=92, y=74
#> (-4.5-START_LON)/x_res
#[1] 93
#> (48.38-START_LAT)/y_res
#[1] 75.42

# Also choose 6 boundary points to reflect the deep ocean boundary condition
# Choose points at:
# x=1, y=50, x=1, y=100, x=1, y=150
# x=50, y=50, x=50, y=100, x=50, y=150

# Libraries
library(fields)
source("~/bin/RScripts/matrixMethods.R")

# Parameters

# Boundary
START_LAT <- 40 
START_LON <- -20

# Grid size
num_lon <- 198
num_lat <- 207
itot<-num_lon*num_lat

x_res <- 1/6
y_res <- 1/9

# Cell centres
lat <- seq(from=(40+y_res/2), length=num_lat, by=y_res)
lon <- seq(from=(-20+x_res/2), length=num_lon, by=x_res)

# File reading

# Get a list of all the files we want to operate on - unformatted binary
fileList <- list.files("/home/simonh/diskx/cs3x/", pattern="mm", full.names=T)
# Bathymetry file - text. Format is 2 header lines followed by num_lat *
# num_lon numbers columns giving the depths (in m?). Mask value is -999.99.
# The line length is 10 numbers, so the first latitude row is 19.8 text
# rows long.
bathyFile <- file("/home/simonh/diskx/cs3x/setupcs3x.dat", "r")
bathyArray <- array(NA, dim=c(num_lon, num_lat))
# Read the bathy file in
junk <- readLines(con=bathyFile, n=2) # first 2 lines are junk
for (i in 1:num_lat){
	junk <- readLines(con=bathyFile, n=20)
	junkStr <- ""
	for (j in 1:20){
		junkStr <- paste(junkStr, junk[j])
	}
	bathy <- strsplit(junkStr, split=' +')
	bathyArray[,i] <- as.numeric(unlist(bathy)[2:(num_lon+1)])
}
close(bathyFile)

# Set land to be NA
bathyArray[which(bathyArray==-999.99)] <- NA

num_files <- length(fileList)
num_mnths <- num_files*12
print(paste("Number of months:", num_mnths))

month_count <-1

brest<-array(NA, dim=c(num_mnths,3))
newlyn<-array(NA, dim=c(num_mnths,3))
ds<-array(NA, dim=c(num_mnths,3,6))

for (i in 1: num_files){
#for (i in 1: 1){ # Just read 1 file for now
	dataFile <- file(fileList[i], "rb")
	print(fileList[i])
	
	for(j in 1:12){# monthly means
#	for(j in 1:1){# Just read 1 record for now

# Read 1st integer that Fortran writes
		junk <- readBin(dataFile, what="integer", n=1, endian="little")
		z <- readBin(dataFile, what="numeric", size=4, n=itot, endian="little")
# Read end of record integer from Fortran
		junk <- readBin(dataFile, what="integer", n=1, endian="little")
	
# Read 1st integer that Fortran writes
		junk <- readBin(dataFile, what="integer", n=1, endian="little")
		u <- readBin(dataFile, what="numeric", size=4, n=itot, endian="little")
# Read end of record integer from Fortran
		junk <- readBin(dataFile, what="integer", n=1, endian="little")
	
# Read 1st integer that Fortran writes
		junk <- readBin(dataFile, what="integer", n=1, endian="little")
		v <- readBin(dataFile, what="numeric", size=4, n=itot, endian="little")
# Read end of record integer from Fortran
		junk <- readBin(dataFile, what="integer", n=1, endian="little")
	    
		nas <- which(z < -1e33)
		z[nas] <- NA
		u[nas] <- NA
		v[nas] <- NA
		
		dim(z) <- c(num_lat, num_lon)
		dim(u) <- c(num_lat, num_lon)
		dim(v) <- c(num_lat, num_lon)
		
		newlyn[month_count,1] <- z[88,91]
		newlyn[month_count,2] <- u[88,91]
		newlyn[month_count,3] <- v[88,91]
		
		brest[month_count,1] <- z[92,74]
		brest[month_count,2] <- u[92,74]
		brest[month_count,3] <- v[92,74]
		

		ds[month_count,1,1] <- z[1,50]
		ds[month_count,1,2] <- z[1,100]
		ds[month_count,1,3] <- z[1,150]
		ds[month_count,1,4] <- z[50,50]
		ds[month_count,1,5] <- z[50,100]
		ds[month_count,1,6] <- z[50,150]

		ds[month_count,2,1] <- u[1,50]
		ds[month_count,2,2] <- u[1,100]
		ds[month_count,2,3] <- u[1,150]
		ds[month_count,2,4] <- u[50,50]
		ds[month_count,2,5] <- u[50,100]
		ds[month_count,2,6] <- u[50,150]
		
		ds[month_count,3,1] <- v[1,50]
		ds[month_count,3,2] <- v[1,100]
		ds[month_count,3,3] <- v[1,150]
		ds[month_count,3,4] <- v[50,50]
		ds[month_count,3,5] <- v[50,100]
		ds[month_count,3,6] <- v[50,150]
		
		month_count <- month_count+1
		
	}
	close(dataFile)
}
bathyArray <- mirror.matrix(bathyArray)
z <- rotate90.matrix(mirror.matrix(z))
u <- rotate90.matrix(mirror.matrix(u))
v <- rotate90.matrix(mirror.matrix(v))

x11()
image.plot(lon,lat,z)
#x11()
#image.plot(lon,lat,u)
#x11()
#image.plot(lon,lat,v)
#x11()
#image.plot(lon,lat,bathyArray)

# Magnitude
magnitude <- sqrt(u^2 + v^2)
x11()
image.plot(lon,lat, magnitude)

# Phase
phase<-array(NA,dim=c(num_lon,num_lat))

for (i in 1:num_lon) {
	for (j in 1:num_lat){
		
		vm <- v[i,j]
		um <- u[i,j]
		if (is.na(um)||is.na(vm)){
			phase[i,j] <- NA
		} else {
			if ((um>0)&&(vm>0)){
				phase[i,j] <- atan(abs(vm)/abs(um))*180/pi
			}
			if ((um<0)&&(vm>0)){
				phase[i,j] <- 180-atan(abs(vm)/abs(um))*180/pi
			}
			if ((um>0)&&(vm<0)){
				phase[i,j] <- 360-atan(abs(vm)/abs(um))*180/pi
			}
			if ((um<0)&&(vm<0)){
				phase[i,j] <- 180+atan(abs(vm)/abs(um))*180/pi
			}
		}		
	}
}
x11()
image.plot(lon,lat, phase)

time <- seq(from=as.Date("1958/1/15"),to=as.Date("2001/12/15"),by="1 month")
x11()
plot(time, newlyn[,1]-mean(newlyn[,1]), col='red', type='l')
lines(time, brest[,1]-mean(brest[,1]), col='blue')

#save(file="brestNewlynCS3X.RData", brest, newlyn, ds, time)