# readPlymouthPressures.R
# 
# Script to read in and reformat Plymouth air pressures from the Met Office
# via PLW
#
# Author: simonh
###############################################################################

plymouth <- read.table("~simonh/diskx/polcoms/brestNewlyn/pressure/Plymouth_03827_18502004_v2.fts", 
                na.strings = "-999.9", skip = 1, 
                col.names=c("Year","Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"))
plymouth <- as.matrix(plymouth, dimnames=NULL)
# Ignore first column which is the year and transpose
pressure <- t(plymouth[,2:13])
# Convert to a vector
pressure <- as.vector(pressure)
# Create time dimension
time <- seq.Date(from=as.Date("1850/1/15"), to=as.Date("2004/12/15"), by="1 month")

save(file="plymouthPressures.RData", time, pressure)
