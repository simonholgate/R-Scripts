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
#%%%%%%%%%%%%%%%%%%
# Old Matlab code
#%%%%%%%%%%%%%%%%%%
#
#nc = netcdf('/users/plw/slp.mon.mean.nc', 'nowrite');
#variables = var(nc);
#Lat=variables{1}(:);
#Lon=variables{2}(:);
#Time=variables{3}(:);
#slpData=variables{4}(:);
#% 51 41 S  57 49 W  STANLEY II
#stanley2Lon=intersect(find(Lon>=301),find(Lon<303));
#stanley2Lat=intersect(find(Lat>=-53),find(Lat<-50));
#lengthT=length(Time);
#t=linspace(1948+1/12,1948+lengthT*1/12,lengthT);
#stanleyPressure=squeeze(slpData(:,stanley2Lat,stanley2Lon));
#plot(t,stanleyPressure);
#array=[t',stanleyPressure];
#save -ascii -tabs stanleyPressurePLWNew.txt array
#close(nc);
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

#>      for( i in 1:nc$ndims ) {
#+              print(paste("Here is information about dimension
#number",i,":"))
#+              d <- nc$dim[[i]]
#+              print(paste("    Name  :",d$name))
#+              print(paste("    Units :",d$units))
#+              print(paste("    Length:",d$len))
#+              print(paste("    Unlimited:",d$unlim))
#+              }
#[1] "Here is information about dimension number 1 :"
#[1] "    Name  : lon"
#[1] "    Units : degrees_east"
#[1] "    Length: 144"
#[1] "    Unlimited: FALSE"
#[1] "Here is information about dimension number 2 :"
#[1] "    Name  : lat"
#[1] "    Units : degrees_north"
#[1] "    Length: 73"
#[1] "    Unlimited: FALSE"
#[1] "Here is information about dimension number 3 :"
#[1] "    Name  : time"
#[1] "    Units : hours since 1-1-1 00:00:0.0"
#[1] "    Length: 666"
#[1] "    Unlimited: TRUE"
Lon <- nc$dim[[1]]$vals
Lat <- nc$dim[[2]]$vals
Time <- nc$dim[[3]]$vals
NCEPslp <- nc$var[[1]]
slp <- get.var.ncdf(nc, NCEPslp)
lengthT <- nc$dim[[3]]$len
t <- seq(from=(1948+1/12),to=(1948+lengthT*1/12),length.out=lengthT)

#stanley2Lon <- intersect(which(Lon>=301), which(Lon<303))
#stanley2Lat <- intersect(which(Lat>=-53), which(Lat<(-50)))
#stanleyPressure=slp[stanley2Lon,stanley2Lat,]
#plot(t,stanleyPressure, type='l')

# Get the data for the 177 stations so we can extract the slp at each point
load("/diskx/users/simonh/177Stations/stns177.RData")
slpLons <- stns177$slondec
slpLats <- stns177$slatdec
slpArray <- array(NA, dim=c(lengthT, 177))
for (k in 1:177) {
  rndLon <- round(slpLons[k])
  rndLat <- round(slpLats[k])
  if (rndLon>=359) {
    i <- union(which(Lon>=(rndLon-1)), which(Lon<=(rndLon-360+1)))
  } else {
    i <- intersect(which(Lon>=(rndLon-1)), which(Lon<=(rndLon+1)))
  }
  j <- intersect(which(Lat>=(rndLat-1)), which(Lat<=(rndLat+1)))
  slpArray[,k] <- slp[i,j,]
}

# Calculate the mean series
meanSshNCEP <- vector(mode="numeric", length=lengthT)
for (i in 1:lengthT){
  meanSshNCEP[i] <- mean(slpArray[i,], na.rm=TRUE)
}

# Convert slp to ssh 
# Already in mb so no need to convert from Pa.
# meanSshNCEP <- meanSshNCEP*0.01
# In converting to mm from mb remember that 1mb = 9.948mm of sea surface
# height
meanSshNCEP <- meanSshNCEP*9.948
 
# Compare with the Hadley set
load("/diskx/users/simonh/177Stations/meanSshHadley2.RData")

plot(s177SlpAnnualYrs, meanSshHadley2, type='l', col='red', lwd=2)
lines(t, meanSshNCEP, col='blue', lwd=2)
