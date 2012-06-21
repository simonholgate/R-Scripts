###############################################################################
# Script to write netcdf files of the 2D daily mean sea 
# level calculated from POLCOMMS model output along the TOPEX/Jason altimetry
# tracks provided by Laurent Testut (see output from readModelAltimTracks.R). 
# 
#
# Simon Holgate, August 2005
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
load("~/diskx/polcoms/brestNewlyn/tracks/modelAltimTracks.RData")

library(ncdf)

JD <- julian(monthsArray, origin = as.Date("1950-01-01"))

#Track 61
x <- altimTrackPoints61[,1]
x[which(is.na(x))] <- 1.0e30
y <- altimTrackPoints61[,2]
y[which(is.na(y))] <- 1.0e30

X <- dim.def.ncdf("Longitude", "degreesE", x)
Y <- dim.def.ncdf("Latitude", "degreesN", y)
T <- dim.def.ncdf("Time", "Julian days since 1950-01-01", JD)

LON <- var.def.ncdf(name="Lon", units="degreesE", dim=X, missval=1.0e30, 
                longname="Longitude", prec="single")
LAT <- var.def.ncdf(name="Lat", units="degreesN", dim=Y, missval=1.0e30, 
                longname="Latitude", prec="single")
SSH <- var.def.ncdf(name="MMSSH", units="millimetres", dim=list(X,T), missval=1.0e30, 
                longname="Monthly Mean Sea Surface Height", prec="single")

# Create a netCDF file with these variables
ncnew <- create.ncdf( "mmssh.track061.nc", list(LON,LAT,SSH))
put.var.ncdf(ncnew, LON, altimTrackPoints61[,1], start=1, count=length(x))
put.var.ncdf(ncnew, LAT, altimTrackPoints61[,2], start=1, count=length(y))

for( i in 1:length(JD)) 
             put.var.ncdf( ncnew, SSH, track61MonthlyMean[i,], start=c(1,i), count=c(-1,1) )
close.ncdf(ncnew)

# Track 70
x <- altimTrackPoints70[,1]
x[which(is.na(x))] <- 1.0e30
y <- altimTrackPoints70[,2]
y[which(is.na(y))] <- 1.0e30

X <- dim.def.ncdf("Longitude", "degreesE", x)
Y <- dim.def.ncdf("Latitude", "degreesN", y)
T <- dim.def.ncdf("Time", "Julian days since 1950-01-01", JD)

LON <- var.def.ncdf(name="Lon", units="degreesE", dim=X, missval=1.0e30, 
                longname="Longitude", prec="single")
LAT <- var.def.ncdf(name="Lat", units="degreesN", dim=Y, missval=1.0e30, 
                longname="Latitude", prec="single")
SSH <- var.def.ncdf(name="MMSSH", units="millimetres", dim=list(X,T), missval=1.0e30, 
                longname="Monthly Mean Sea Surface Height", prec="single")

# Create a netCDF file with these variables
ncnew <- create.ncdf( "mmssh.track070.nc", list(LON,LAT,SSH))
put.var.ncdf(ncnew, LON, altimTrackPoints70[,1], start=1, count=length(x))
put.var.ncdf(ncnew, LAT, altimTrackPoints70[,2], start=1, count=length(y))

for( i in 1:length(JD)) 
             put.var.ncdf( ncnew, SSH, track70MonthlyMean[i,], start=c(1,i), count=c(-1,1) )
close.ncdf(ncnew)

# Track 239
x <- altimTrackPoints239[,1]
x[which(is.na(x))] <- 1.0e30
y <- altimTrackPoints239[,2]
y[which(is.na(y))] <- 1.0e30

X <- dim.def.ncdf("Longitude", "degreesE", x)
Y <- dim.def.ncdf("Latitude", "degreesN", y)
T <- dim.def.ncdf("Time", "Julian days since 1950-01-01", JD)

LON <- var.def.ncdf(name="Lon", units="degreesE", dim=X, missval=1.0e30, 
                longname="Longitude", prec="single")
LAT <- var.def.ncdf(name="Lat", units="degreesN", dim=Y, missval=1.0e30, 
                longname="Latitude", prec="single")
SSH <- var.def.ncdf(name="MMSSH", units="millimetres", dim=list(X,T), missval=1.0e30, 
                longname="Monthly Mean Sea Surface Height", prec="single")

# Create a netCDF file with these variables
ncnew <- create.ncdf( "mmssh.track239.nc", list(LON,LAT,SSH))
put.var.ncdf(ncnew, LON, altimTrackPoints239[,1], start=1, count=length(x))
put.var.ncdf(ncnew, LAT, altimTrackPoints239[,2], start=1, count=length(y))

for( i in 1:length(JD)) 
             put.var.ncdf( ncnew, SSH, track239MonthlyMean[i,], start=c(1,i), count=c(-1,1) )
close.ncdf(ncnew)

# Track 248
x <- altimTrackPoints248[,1]
x[which(is.na(x))] <- 1.0e30
y <- altimTrackPoints248[,2]
y[which(is.na(y))] <- 1.0e30

X <- dim.def.ncdf("Longitude", "degreesE", x)
Y <- dim.def.ncdf("Latitude", "degreesN", y)
T <- dim.def.ncdf("Time", "Julian days since 1950-01-01", JD)

LON <- var.def.ncdf(name="Lon", units="degreesE", dim=X, missval=1.0e30, 
                longname="Longitude", prec="single")
LAT <- var.def.ncdf(name="Lat", units="degreesN", dim=Y, missval=1.0e30, 
                longname="Latitude", prec="single")
SSH <- var.def.ncdf(name="MMSSH", units="millimetres", dim=list(X,T), missval=1.0e30, 
                longname="Monthly Mean Sea Surface Height", prec="single")

# Create a netCDF file with these variables
ncnew <- create.ncdf( "mmssh.track248.nc", list(LON,LAT,SSH))
put.var.ncdf(ncnew, LON, altimTrackPoints248[,1], start=1, count=length(x))
put.var.ncdf(ncnew, LAT, altimTrackPoints248[,2], start=1, count=length(y))

for( i in 1:length(JD)) 
             put.var.ncdf( ncnew, SSH, track248MonthlyMean[i,], start=c(1,i), count=c(-1,1) )

close.ncdf(ncnew)




#save(file="modelAltimTracks.RData", altimTrackPoints, trackMonthlyMean, 
#  altimTrackPoints61, track61MonthlyMean, altimTrackPoints70, track70MonthlyMean,
#  altimTrackPoints239, track239MonthlyMean, altimTrackPoints248, track248MonthlyMean)