# R-Script to read sea level data output files from CS3X
#
# Note that this script needs to access the CS3X files on the main POL network 
# as they are too big to transfer over to livtech59
# 
# Author: simonh
###############################################################################

# Details on the CS3X format from Chris....

#Hi Simon, (cc Kevin)
#
#Perhaps the easiest way to read/reformat the data is to use my routine 
#which converts the model output to netCDF as a basis.
#
#The model output is unformatted direct access, real*4, little-endian.  
#The fiddly thing is that it's in 'compact addressing'.  This is a 
#		hangover from the days of vector computing, where the 2D model arrays 
#		were slotted into 1D arrays with the land removed.  Therefore, you need 
#		the full land mask to transform them back to regular 2D lat-lon.
#		The prog which converts to netCDF contains all the required stuff, and 
#		you have the option of writing netCDF or just hacking that part out and 
#		writing to whatever format you choose.
#		
#		I've put the progs and input files in ~cwi/public/CS3X_ERA40_tools.  
#Have a look at them and I'll explain anything that isn't obvious 
#(possibly quite a lot of the code).
#
#excw_nc.f is the main prog
#expand.fil contains:
#		i. the path to the compact address form of the model output to be 
#converted (e.g. /bank/cwi/CS3X_ERA40/tid.0000Z01011998.0000Z01011999.uda)
#ii. the setup file which contains the total number of compact address 
#points and the start and end points of the 2D array (setupcs3xSGIl_LE.uda)
#iii. the land mask file, since the setup file does not contain info 
#about land between start and end points of each latitude; depths < 5 m 
#are masked (setupcs3x.dat)
#iv. the output file (e.g. /scratch/cwi/tid.0000Z01011998.0000Z01011999.nc)
#
#So, you probably only want to change the first and fourth lines of 
#expand.fil.
#
#excw_nc.f contains the line
#parameter (START_LAT = 40, START_LON = -20)
#which tells you the most southern and western position of the grid 
#boundaries (40 N and 20 W).  The grid is 198 longitudes (1/6 deg res) x 
#207 latitudes (1/9 deg res).
#Therefore the most northern and eastern boundaries are at 63 N and 13 
#E.  Equivalently, the centre of the top-left cell is 62.9444 deg N 
#19.91667 deg W.
#
#Here are some trivial examples of cell centre calculations which I 
#copied from an old email:
#		
#		e.g. longitude of cell (I, J) is -19.91667 + (I-1)/6 and if the number 
#remains negative then it is West of Greenwich
#
#so (100,J) is -19.91667 + 99/6 = -3.417 = 3.417 deg W
#
#similar sort of thing for latitudes. Latitude = 62.9444 - (J-1)/9
#
#Hopefully the model output filenames in /bank/cwi/CS3X_ERA40 will be 
#obvious.  The dataset contains the tidal component of SL (from a run 
#with no met forcing, just tidal forcing at the boundaries) and the total 
#SL (same, but with met forcing also).
#
#This run is forced by a 1 degree lat-lon version of the 6 hourly ERA-40 
#winds and SLP, linearly interpolated to the model grid and timestep (30 
#sec).  In case you're interested, the original ERA-40 run is actually on 
#		an N80 reduced Gaussian Grid.  That's what should be taken as 
#representative of the actual resolution of the re-analysis (see 
#				http://badc.nerc.ac.uk/data/ecmwf-op/grids.html#n80_reduced).
#
#Let me know how you get on.
#
#Chris
#

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
x_res <- 1/6
y_res <- 1/9

# Cell centres
lat <- seq(from=(40+y_res/2), length=num_lat, by=y_res)
lon <- seq(from=(-20+x_res/2), length=num_lon, by=x_res)

# File reading
nrecs <- 8761 # (number of hours in a year)+1. It is the number of records we read

# Get a list of all the files we want to operate on. Use the "total" sea level run
# not the "tide only" run  - unformatted binary
fileList <- list.files("/bank/cwi/CS3X_ERA40", pattern="tot", full.names=T)
# Bathymetry file - text. Format is 2 header lines followed by num_lat *
# num_lon numbers columns giving the depths (in m?). Mask value is -999.99.
# The line length is 10 numbers, so the first latitude row is 19.8 text
# rows long.
bathyFile <- file("/login/cwi/public/CS3X_ERA40_tools/setupcs3x.dat", "r")
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

# Read compact array information from the setup file
# Setup file - unformatted binary
setupFile <- file("/login/cwi/public/CS3X_ERA40_tools/setupcs3xSGIl_LE.uda", "rb")

# Read 1st integer that Fortran writes
junk <- readBin(setupFile, what="integer", n=1, endian="little")
params <- readBin(setupFile, what="integer", n=4, endian="little")
itot <- params[3]
# Read end of record integer from Fortran
junk <- readBin(setupFile, what="integer", n=1, endian="little")

# Read 1st integer that Fortran writes
junk <- readBin(setupFile, what="integer", n=1, endian="little")
ntrns <- readBin(setupFile, what="integer", n=num_lat, endian="little")
# Read end of record integer from Fortran
junk <- readBin(setupFile, what="integer", n=1, endian="little")

# Read 1st integer that Fortran writes
junk <- readBin(setupFile, what="integer", n=1, endian="little")
ntrnt <- readBin(setupFile, what="integer", n=num_lat, endian="little")
# Read end of record integer from Fortran
junk <- readBin(setupFile, what="integer", n=1, endian="little")

# Read 1st integer that Fortran writes
junk <- readBin(setupFile, what="integer", n=1, endian="little")
nsum <- readBin(setupFile, what="integer", n=num_lat, endian="little")
# Read end of record integer from Fortran
junk <- readBin(setupFile, what="integer", n=1, endian="little")

close(setupFile)

num_files <- length(fileList)
z <- array(NA, dim=c(num_lon, num_lat))
u <- z
v <- z

dates <- array(NA, dim=c(nrecs,5))

#for (i in 1: num_files){
for (i in 1: 1){ # Just read 1 file for now
	dataFile <- file(fileList[i], "rb")
	
	for(j in 1:nrecs){# Dumps are hourly
#	for(j in 1:1){# Just read 1 record for now
# Read 1st integer that Fortran writes
		junk <- readBin(dataFile, what="integer", n=1, endian="little")
		dates[j,] <- readBin(dataFile, what="integer", n=5, endian="little")
# Read end of record integer from Fortran
		junk <- readBin(dataFile, what="integer", n=1, endian="little")
	
# Read 1st integer that Fortran writes
		junk <- readBin(dataFile, what="integer", n=1, endian="little")
		zs <- readBin(dataFile, what="numeric", size=4, n=itot, endian="little")
# Read end of record integer from Fortran
		junk <- readBin(dataFile, what="integer", n=1, endian="little")
	
# Read 1st integer that Fortran writes
		junk <- readBin(dataFile, what="integer", n=1, endian="little")
		us <- readBin(dataFile, what="numeric", size=4, n=itot, endian="little")
# Read end of record integer from Fortran
		junk <- readBin(dataFile, what="integer", n=1, endian="little")
	
# Read 1st integer that Fortran writes
		junk <- readBin(dataFile, what="integer", n=1, endian="little")
		vs <- readBin(dataFile, what="numeric", size=4, n=itot, endian="little")
# Read end of record integer from Fortran
		junk <- readBin(dataFile, what="integer", n=1, endian="little")
	
		isum <- 1
		for (k in 1:num_lat){
			ie <- ntrns[k]
			is <- ntrnt[k] + 1
		
			for (l in is:ie){
				isum <- isum+1
				if(is.finite(bathyArray[l,k])){
					if (bathyArray[l,k] >= 5){
						z[l,k] <- zs[isum]
						u[l,k] <- us[isum]
						v[l,k] <- vs[isum]
					}
				}
			}	
		}
	}
	close(dataFile)
}

#x11()
#image.plot(lon,lat,mirror.matrix(z))
#x11()
#image.plot(lon,lat,mirror.matrix(u))
#x11()
#image.plot(lon,lat,mirror.matrix(v))
#x11()
#image.plot(lon,lat,mirror.matrix(bathyArray))