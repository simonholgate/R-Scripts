# Read in NCEP netCDF file of monthly values and write out binary file of
# annual values
#
# Note: 
# Refer to http://www.cdc.noaa.gov/cdc/data.ncep.reanalysis.derived.html#surface
#
# These variables are instantaneous values at the reference time.
#
# These variables are averages of instantaneous values at the 4 reference times;
# 0,6,12 and 18z over the averaging period (month, day)
#
# Spatial coverage:
#
#    * 2.5-degree latitude x 2.5-degree longitude global grid with 144x73 points
#    * 90N-90S, 0E-357.5E 
#
# Temporal coverage:
#
#    * 1/1/1948 (1958 for some variables)- close to present (updated near 
#      beginning of month)
#    * Data for the current year from the CDAS program is being made available 
#      along with the historical data. Files are updated through Jan 2002 
#
# Levels:
#
#    * Surface or near the surface (.995 sigma level)
#
#
library(ncdf)
nc <- open.ncdf("/diskx/users/simonh/NCEP/slp.mon.mean.nc")
#> nc$var[[1]]$name
#[1] "slp"
#> nc$var[[1]]$varsize
#[1] 144  73 666
#> nc$var[[1]]$units  
#[1] "millibars"
#> nc$var[[1]]$ndims
#[1] 3
#> nc$var[[1]]$missval
#[1] -9.96921e+36

NCEPslp <- nc$var[[1]]
slp <- get.var.ncdf(nc, NCEPslp)
close.ncdf(nc)

# Write out monthly file
# Use time as the leading dimension
con<-file("ncep.monthlies.le.kij.bin","wb")
 
for (i in 1:144){
  for (j in 1:73){
    temp<-slp[i,j,]
    writeBin(as.vector(temp),con, size=4)
  }
}
close(con)

# Use time as the last dimension
con<-file("ncep.monthlies.le.bin","wb")
for (i in 1:666){
  temp<-slp[,,i]
  writeBin(as.vector(temp),con, size=4)
}
close(con)

## Make annual means (integers as with HadSLP2)
#nyr <- floor(666/12) # 55 years
#nlon <- 144
#nlat <- 73
#slpAnnualArray <- array(NA, dim=c(nlon,nlat,nyr))
#
#for (k in seq(from=1, to=(nyr*12), by=12)){
#  for (i in 1:nlon) {
#    for (j in 1:nlat) {
#      temp <- slp[i,j,(k:k+11)]
## Replace missing values with NAs
#      l <- which(temp<=-9e36)
#      temp[l] <- NA
#      slpAnnualArray[i,j,((k-1)/12+1)] <- mean(temp,na.rm=TRUE)
#    }
#  }
#}
  
