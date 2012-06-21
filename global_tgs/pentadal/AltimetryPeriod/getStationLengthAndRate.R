# Script to interrogate the database and produce a formatted list of record length
# and annual rate of sea level change, similar to what Jabber produced in Fortran.
# The rates of SLR chnage in Jabber are corrected for inverse barometer effects
# using HadSLP2. Also included in the output are the lats and lons of the station.
#
# Jabber read in a list of country and station codes which it the queried the 
# database about and the same approach will be made here.
#
# Our list will be all RLR stations initially.
#

# Format of output from jabber to reproduce:
#% head temp27.1984-2004.hadley2.apc
# 12 170 184                      355.05  51.92                             7.9
# 19 170 161                      354.45  50.10                             1.8

#% head temp26.1984-2004.hadley2.apc
# 170/184      12  1989 - 2004     7.88 +/-  1.42    22.9   51 55 N  04 57 W   FISHGUARD II                  
# 170/161 241  19  1984 - 2004     1.75 +/-  0.86    21.8   50 06 N  05 33 W   NEWLYN                        

#library(RMySQL)
library(ROracle)

# Read in list of stations to process from a text file
# Format of cmsexec files (R denotes and RLR station is required):
#970162R
#970191R
#970201R
#970211R
#A  001R
#A  003R
#A  005R
#A  024R

stnList <- read.fwf('cmsexec', widths=c(3,3,1), 
  col.names=c('ccode','scode','rlr'), colClasses=c('character'))

startYr <- 1992
endYr <- 2008

# Only replicate temp27 ouput format for now
out27 <- file(paste("temp27.",startYr,"-",endYr,".hadley2.apc", sep=""), open="w")
# We also want a list in the format of includedTideGauges.mod.txt which lists
# just the stations that have at least 70% completeness over this period
outInc <- file(paste("includedTideGaugesAltimPeriod.txt", sep=""), open="w")

pc70 <- round(length(startYr:endYr)*0.7)

numStns <- length(stnList[,1])

#drv <- dbDriver("MySQL")
drv <- dbDriver("Oracle")
con <- dbConnect(drv, "psmsl2/psmsl2")

for (i in 1:numStns){
  if(stnList$rlr[i]=='R'){ # RLR station
	selectString <- paste("SELECT yr,ann,rlrfac,missdays FROM data",
	" WHERE ccode='",stnList$ccode[i],"' AND scode='",
	stnList$scode[i],"' AND RLRFAC IS NOT NULL ", 
	"AND yr>=", startYr," AND yr<=", endYr,
	" ORDER BY yr", sep="")
  } else { # Metric data
  	selectString <- paste("SELECT yr,ann,rlrfac,missdays FROM data",
	" WHERE ccode='",stnList$ccode[i],"' AND scode='",
	stnList$scode[i],"' ",
    "AND yr>=", startYr," AND yr<=", endYr,
	" ORDER BY yr", sep="")
  }
  
  res <- dbGetQuery(con,selectString)

# Check we have any results otherwise skip to next record
  if(dim(res)[1] == 0) next

# For comparison with jabber output, ignore years that have missing months
  indices <- which(substr(res[,4],25,26)!="XX")
# We can include years where the missing days string is NA too
  indices <- c(indices, which(is.na(res[,4])))
  indices <- sort(indices)

# We only want to include years where there is actually an annual value
  indices <- indices[which(is.finite(res[indices,2]))]
# Number of years of data
  numYrs <- length(indices)
# Check that we still have some data left (we need at least 2 years for a trend), 
# else skip to next record...
  if (numYrs<2) next 
  
# Calculate trend
  model <- lm(res[indices,2] ~ res[indices,1])
  trend <- model$coef[2]

# Get stations information from sttn
  selectString <- paste("SELECT slat, slon, sname FROM sttn",
	" WHERE ccode='",stnList$ccode[i],"' AND scode='",
	stnList$scode[i],"'", sep="")
	
  res <- dbGetQuery(con,selectString)
  
# Convert lons and lats to decimal: 
# For lats, N is +ve, S is -ve.
# For lons, values are degrees E of Greenwich

  slat <- strsplit(as.character(res[1]), '[[:blank:]]+')
  slon <- strsplit(as.character(res[2]), '[[:blank:]]+')
  
  slatdec <- as.numeric(slat[[1]][2]) + as.integer(slat[[1]][3])/60
  if(as.character(slat[[1]][4]) == 'S'){
    slatdec <- -slatdec
  }
  
  if (slon[[1]][1] == ""){
    j <- 2
  } else {
    j <- 1
  }
  
  slondec <- as.numeric(slon[[1]][j]) + as.numeric(slon[[1]][j+1])/60
  if(as.character(slon[[1]][j+2]) == 'W'){
    slondec <- 360-slondec
  }
  
# Write it out to file
# 12 170 184                      355.05  51.92                             7.9
  line27 <- sprintf("%3i %03s %03s                      %6.2f  %5.2f                             %3.1f", 
    numYrs, stnList$ccode[i], stnList$scode[i], slondec, slatdec, trend)
  writeLines(line27, con=out27)

# Now do the included stations list  
#970;031;CHARLOTTETOWN;296.88;46.23;42
  if (numYrs>=pc70){
    if (stnList$ccode[i] == 'A  ') stnList$ccode[i] <- '999'
    lineInc <- paste(stnList$ccode[i], stnList$scode[i], res[3], 
      slondec, slatdec, numYrs, sep=";")
    writeLines(lineInc, con=outInc)  
  }
}
dbDisconnect(con)
dbUnloadDriver(drv)

close(out27)
close(outInc)
