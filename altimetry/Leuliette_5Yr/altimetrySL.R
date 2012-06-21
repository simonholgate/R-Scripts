# Calculate 5 year run mean rates of SLR from altimetry data supplied from
# Eric Leuliettes website: http://ibis.grdl.noaa.gov/SAT/slr/
# In publications, presentations, or on web pages based on LSA data the
# following acknowledgment should be included.
# "Altimetry data are provided by NOAA Laboratory for Satellite Altimetry."

# Note that Topex (tx) and Jason (j1) are reversed in this file so I have
# un-reversed the names when creating my own variables to avoid my confusion

library(ncdf)
altimetry <- new.env()

# Data from Topex and Jason
nc1 <-
open.ncdf("~/diskx/177StationsUpdate2009/slr_sla_gbl_link_txj1_90.nc")
# Multi satellite data
nc2 <-
open.ncdf("~/diskx/177StationsUpdate2009/slr_sla_gbl_link_all_66.nc")
#
     print(paste("File",nc1$filename,"contains",nc1$nvars,"variables"))
     for( i in 1:nc1$nvars ) {
             v <- nc1$var[[i]]
             print(paste("Here is information on variable number",i))
             print(paste("   Name: ",v$name))
             print(paste("   Units:",v$units))
             print(paste("   Missing value:",v$missval))
             print(paste("   # dimensions :",v$ndims))
             print(paste("   Variable size:",v$varsize))
             }
#[1] "Here is information on variable number 1"
#[1] "   Name:  sla_tx"
#[1] "   Units: mm"
#[1] "   Missing value: 1e+30"
#[1] "   # dimensions : 1"
#[1] "   Variable size: 440"
#[1] "Here is information on variable number 2"
#[1] "   Name:  sla_j1"
#[1] "   Units: mm"
#[1] "   Missing value: 1e+30"
#[1] "   # dimensions : 1"
#[1] "   Variable size: 246"
     print(paste("File",nc1$filename,"contains",nc1$ndims,"dimensions"))
     for( i in 1:nc1$ndims ) {
             print(paste("Here is information about dimension number",i,":"))
             d <- nc1$dim[[i]]
             print(paste("    Name  :",d$name))
             print(paste("    Units :",d$units))
             print(paste("    Length:",d$len))
             print("    Values:")
             print(d$vals)
             if ( i==1 ){
               altimetry$time_tx <- d$vals
             } else {
               altimetry$time_j1 <- d$vals
             }
             print(paste("    Unlimited:",d$unlim))
             }


     print(paste("File",nc2$filename,"contains",nc2$nvars,"variables"))
     for( i in 1:nc2$nvars ) {
             v <- nc2$var[[i]]
             print(paste("Here is information on variable number",i))
             print(paste("   Name: ",v$name))
             print(paste("   Units:",v$units))
             print(paste("   Missing value:",v$missval))
             print(paste("   # dimensions :",v$ndims))
             print(paste("   Variable size:",v$varsize))
             }
#[1] "Here is information on variable number 1"
#[1] "   Name:  sla_tx"
#[1] "   Units: mm"
#[1] "   Missing value: 1e+30"
#[1] "   # dimensions : 1"
#[1] "   Variable size: 440"
#[1] "Here is information on variable number 2"
#[1] "   Name:  sla_e2"
#[1] "   Units: mm"
#[1] "   Missing value: 1e+30"
#[1] "   # dimensions : 1"
#[1] "   Variable size: 333"
#[1] "Here is information on variable number 3"
#[1] "   Name:  sla_g1"
#[1] "   Units: mm"
#[1] "   Missing value: 1e+30"
#[1] "   # dimensions : 1"
#[1] "   Variable size: 269"
#[1] "Here is information on variable number 4"
#[1] "   Name:  sla_j1"
#[1] "   Units: mm"
#[1] "   Missing value: 1e+30"
#[1] "   # dimensions : 1"
#[1] "   Variable size: 246"
#[1] "Here is information on variable number 5"
#[1] "   Name:  sla_n1"
#[1] "   Units: mm"
#[1] "   Missing value: 1e+30"
#[1] "   # dimensions : 1"
#[1] "   Variable size: 248"
altimetry$tx <- get.var.ncdf(nc1, nc1$var[[1]])
altimetry$j1 <- get.var.ncdf(nc1, nc1$var[[2]])

close.ncdf(nc1)
close.ncdf(nc2)

plot(altimetry$time_j1, altimetry$j1, type='l', col='blue', 
  xlim=c(min(altimetry$time_tx), max(altimetry$time_j1)), 
  ylim=c(min(altimetry$tx),max(altimetry$j1)))
lines(altimetry$time_tx,altimetry$tx, type='l', col='red')

# Calculate 5 year running mean rates
# Average cycle length in this is 10.9 days so there are ~ 167 cycles in 5
# years
altimetry$len_tx <- length(altimetry$tx)
altimetry$rates_tx <- vector(mode="double", length=(altimetry$len_tx-167))
for (i in 1:(altimetry$len_tx-167)){
  fit <- lm(altimetry$tx[i:(i+166)] ~ altimetry$time_tx[i:(i+166)])
  altimetry$rates_tx[i] <- fit$coef[2]
}
altimetry$len_j1 <- length(altimetry$j1)
altimetry$rates_j1 <- vector(mode="double", length=(altimetry$len_j1-167))
for (i in 1:(altimetry$len_j1-167)){
  fit <- lm(altimetry$j1[i:(i+166)] ~ altimetry$time_j1[i:(i+166)])
  altimetry$rates_j1[i] <- fit$coef[2]
}

# Do this again with full year means
# Topex full years: 1993-2004
# Jason full years: 2002-2007
altimetry$tx_ann <- vector(mode="numeric", length=11)
altimetry$tx_midpt <- vector(mode="numeric", length=11)
altimetry$j1_ann <- vector(mode="numeric", length=5)
altimetry$j1_midpt <- vector(mode="numeric", length=5)
for (i in 1:11){
  junk <- intersect(which(altimetry$time_tx>=(i+1992)),
    which(altimetry$time_tx<=(i+1996.999999)))
  altimetry$tx_midpt[i] <- mean(altimetry$time_tx[junk], na.rm=T)
  fit <- lm(altimetry$tx[junk] ~ altimetry$time_tx[junk])
  altimetry$tx_ann[i] <- fit$coef[2]
}
for (i in 1:5){
  junk <- intersect(which(altimetry$time_j1>=(i+2001)),
    which(altimetry$time_j1<=(i+2005.999999)))
  altimetry$j1_midpt[i] <- mean(altimetry$time_j1[junk], na.rm=T)
  fit <- lm(altimetry$j1[junk] ~ altimetry$time_j1[junk])
  altimetry$j1_ann[i] <- fit$coef[2]
}
#
x11()
plot(altimetry$tx_midpt, altimetry$tx_ann, col='red', type='l', 
  ylim=c(0.75, 4.25), xlim=c(1995,2008))
lines(altimetry$j1_midpt, altimetry$j1_ann, col='blue')
# Account for GIA of 0.3 mm/yr
lines(altimetry$tx_midpt, (altimetry$tx_ann-0.3), col='red', lty='dashed')
lines(altimetry$j1_midpt, (altimetry$j1_ann-0.3), col='blue', lty='dashed')

rm(d,fit,i,junk,nc1,nc2,v)

