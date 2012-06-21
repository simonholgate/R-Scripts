# readBrestPressures.R
# 
# Script to read in and reformat Brest air pressures from the WMO
# since 1951
#
# Author: simonh
###############################################################################

brest <- read.table("/home/simonh/diskx/polcoms/brestNewlyn/pressure/slp_brest1951.dat", 
                na.strings = "-999.9", skip = 5, 
                col.names=c("Year","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
brest <- as.matrix(brest, dimnames=NULL)
# Ignore first column which is the year and transpose
pressure <- t(brest[,2:13])
# Convert to a vector
pressure <- as.vector(pressure)
# Create time dimension
time <- seq.Date(from=as.Date("1951/1/15"), to=as.Date("2002/12/15"), by="1 month")

save(file="brestPressures.RData", time, pressure)
