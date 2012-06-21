##****************************************************************************************************##
## Calculate annual mean SLP for Key West and San Francisco area using ACRE SLP data
## 
##****************************************************************************************************##

##****************************************************************************************************##
########################################################################################################
## Functions for use below                                                                            ##
########################################################################################################
##****************************************************************************************************##

annual.2d.slp <- function(slpArray){
  ## Convert a 2D array of monthly data into a 2D annual array.
  ## There should be nstns columns and nmon rows
  
  nstns <- dim(slpArray)[2]
  nmon <- dim(slpArray)[1]
  nyr <- nmon/12


  slpAnnualArray <- array(NA, dim=c(nyr,nstns))

  for(i in 1:nstns){
      ## Make a temporray vector of each station
      temp <- slpArray[,i]
      ## Reshape to 3d array with nstns*nyr*12
      dim(temp) <- c(12, nyr)
      ## Place the nyr column means from the temporary array into the column for stn i
      slpAnnualArray[,i] <- colMeans(temp)
  }
  
  slpAnnualArray
}

##*******************************************************************************************************

#########################
## Non-functional part ##
#########################

#library(fields)

nyr<-138

slp.yrs <- c(1871:2008)

##*******************************************************************************************************
## Key West/San Francisco annual pressure

## Monthly data 2 stations. First two are Key West and San Francisco. Variable name is slpKeyWestStns.
load("~/diskx/polcoms/brestNewlyn/analysis/paper/keyWestEofACRE/keyWestACREslp.RData")

# Convert Pa to Mb
slpKeyWestStns <- slpKeyWestStns/100

nstns <- 2

slpAnnArray <- annual.2d.slp(slpKeyWestStns)

save(slpAnnArray, file="keyWestACREannualSLP.RData")
