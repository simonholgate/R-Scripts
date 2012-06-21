###############################################################################
## keyWestFunctionsACRE.R                                                    ##
## Abstracted functions which are common to Key West models                  ##
## Called by control files for each model which define the other variables   ##
## and conduct any other tests                                               ##
###############################################################################

source("~/Dropbox/BrestNewlyn/sorted.scale.R")

##################
##################
## Main program ##
##################
##################

#################################
## Set up additional variables ##
#################################

## Load met data
met <- new.env()
eofatl <- new.env()
eofpac <- new.env()

## Note that Key West and San Francisco have been converted to mb
load("~/Dropbox/brestNewlynData/analysis/paper/keyWestEofACRE/keyWestACREannualSLP.RData", envir=met)

load("~/Dropbox/brestNewlynData/analysis/paper/eofACRE/nAtlACREslpEOF.RData", envir=eofatl) # Whole N Atl
load("~/Dropbox/brestNewlynData/analysis/paper/eofNPacACRE/nPacACREslpEOF.RData", envir=eofpac) # Whole N Pac

met$nstns <- nstns # Number of EOFs to consider

nPa <- met$nstns # All columns are relevant
met$numYrs <- dim(eofatl$dmSlpAnnual)[1]
eofatl$Pa <- array(NA, dim=c(met$numYrs,met$nstns))
eofpac$Pa <- array(NA, dim=c(met$numYrs,met$nstns))

for(i in 1:met$nstns){
  ## Note that pressures have been converted to mb.
  eofatl$Pa[,i] <- (eofatl$dmSlpAnnual %*% eofatl$svdDmSlpAnnual$v[,i])
  eofpac$Pa[,i] <- (eofpac$dmSlpAnnual %*% eofpac$svdDmSlpAnnual$v[,i])
}

met$yrs <- seq(from=1871, to=2008)

pa.start <- which(met$yrs==start) # Start of ACRE pressures for this model
pa.end <- which(met$yrs==end) # End of ACRE pressures for this model

## Use annual pressures over the periods of the model
## In ACRE Key West is column 1 and San Francisco is column 2
## Note that San Fran and Key West have been converted to mb. Multiply by 100 to return to Pa
keywestPa <-
  met$slpAnnArray[intersect(which(met$yrs>=start),
                                            which(met$yrs<=end)),1]*100

sanfranPa <-
  met$slpAnnArray[intersect(which(met$yrs>=start),
                                            which(met$yrs<=end)),2]*100

time <- c(start:end)

p1p.start <- which(time==p1b.start) # Pressure training data start
p1p.end <- which(time==p1b.end) # Pressure training data end
p2p.start <- which(time==p2b.start) # Pressure prediction start
p2p.end <- which(time==p2b.end) # Pressure prediction end


## Load tg data
tg <- new.env()
load("../tg/tg.RData", envir=tg)

tgbs <- which(tg$sanfran$Year==start) # index of start of San Francisco TG data
tgbe <- which(tg$sanfran$Year==end) # index of end of San Francisco TG data

tgns <- which(tg$keywest$Year==start) # index of start of Key West TG data
tgne <- which(tg$keywest$Year==end) # index of end of Key West TG data

tg_sanfran <- tg$sanfran[tgbs:tgbe,]
tg_keywest <- tg$keywest[tgns:tgne,]

## Calculate the total pressure
sanfranTotalP <- 1025*9.8*tg_sanfran[,2]/1000 + sanfranPa
keywestTotalP <- 1025*9.8*tg_keywest[,2]/1000 + keywestPa

tg.lmRob.sanfranTotalP <- lmRob(sanfranTotalP ~ time, control=lmRobControl, x=T)
tg.lmRob.keywestTotalP <- lmRob(keywestTotalP ~ time, control=lmRobControl, x=T)

## Note that Keywest and Sanfran have been converted to mb. Multiply by 100 to return to Pa.
eofatl$Pa <- eofatl$Pa[pa.start:pa.end,]*100
eofpac$Pa <- eofpac$Pa[pa.start:pa.end,]*100

##***********************************************************************************************************

##################################
## Main calculations start here ##
##################################

###################
## San Francisco ##
###################

#tg.sd.sorted.sanfran <- sorted.scale(nPa, sanfranTotalP, start,end, eofpac$Pa, lmRobControl)
#bmax <- 9

data.sanfran <- data.frame(msl=sanfranTotalP, t=seq(from=start,to=end),
                           eofpac$Pa[,1:met$nstns])                           
#                           eofpac$Pa[,tg.sd.sorted.sanfran$ix[1:bmax]])

tg.lmRob.sanfran <- lmRob(msl ~ ., x=T, y=T, data=data.sanfran, control=lmRobControl, na.action=na.exclude)

#################################
## Test model prediction skill ##
#################################
## Need to build new set of top 9 components for partial set

#tg.sd.sorted.sanfran.partial <- sorted.scale(nPa, sanfranTotalP[p1p.start:p1p.end],
#                                           p1b.start,p1b.end, eofpac$Pa[p1p.start:p1p.end,], lmRobControl)

data.sanfran.partial <- data.frame(msl=sanfranTotalP[p1p.start:p1p.end],
  t=seq(from=p1b.start,to=p1b.end), eofpac$Pa[p1p.start:p1p.end,1:met$nstns])

tg.lmRob.sanfran.partial <- lmRob(msl ~ ., x=T, y=T,
                               data=data.sanfran.partial, control=lmRobControl, na.action=na.exclude)

data.sanfran.new<- data.frame(t=seq(from=p2b.start,to=p2b.end),
                            eofpac$Pa[p2p.start:p2p.end,1:met$nstns])

tg.lmRob.sanfran.pred<- predict.lmRob(tg.lmRob.sanfran.partial, se.fit=T,
                                    newdata=data.sanfran.new, interval='confidence')


##############
## Key West ##
##############

#tg.sd.sorted.keywest <- sorted.scale(nPa, keywestTotalP, start,end, eofatl$Pa, lmRobControl)

## Now construct models of the top 9 components
#nmax <- 9

data.keywest <- data.frame(msl=keywestTotalP, t=seq(from=start,to=end), eofatl$Pa[,1:met$nstns])

tg.lmRob.keywest <- lmRob(msl ~ ., x=T, y=T, data=data.keywest, control=lmRobControl, na.action=na.exclude)

#################################
## Test model prediction skill ##
#################################
## Need to build new set of top 9 components for partial set
## lmRob becomes singular and gives an error with more than 4 components in model 4

#tg.sd.sorted.keywest.partial <- sorted.scale(nPa, keywestTotalP[p1p.start:p1p.end],
#                                            p1n.start,p1n.end, eofatl$Pa[p1p.start:p1p.end,], lmRobControl)
#nmax <- 9

data.keywest.partial <- data.frame(msl=keywestTotalP[p1p.start:p1p.end],
  t=seq(from=p1n.start,to=p1n.end), eofatl$Pa[p1p.start:p1p.end,1:met$nstns])

tg.lmRob.keywest.partial <- lmRob(msl ~ ., x=T, y=T,
                               data=data.keywest.partial, control=lmRobControl, na.action=na.exclude)

data.keywest.new<- data.frame(t=seq(from=p2n.start,to=p2n.end),
                             eofatl$Pa[p2p.start:p2p.end,1:met$nstns])

tg.lmRob.keywest.pred<- predict.lmRob(tg.lmRob.keywest.partial, se.fit=T,
                                     newdata=data.keywest.new, interval='confidence')

