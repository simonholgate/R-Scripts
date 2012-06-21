###############################################################################
## brestNewlynFunctions.R                                                    ##
## Abstracted functions which are common to all models                       ##
## Called by control files for each model which define the other variables   ##
## and conduct any other tests                                               ##
###############################################################################

source("~/Dropbox/BrestNewlyn/sorted.scale.R")

##*******************************************************************************************************
## Robustly linearly detrend a vector. Missing values must be included as NA, but are excluded
## in the regression.
detrend <- function(data_vector){
  df.data_vector <- data.frame(p=data_vector,t=seq(from=1,to=length(data_vector)))
  data.lmRob <- lmRob(p ~ t, x=T, y=T, data=df.data_vector, control=lmRobControl,
                              na.action=na.exclude)
## We want the returned vector to be the same length as the original vector so use the $x values
  resid.data_vector <- vector(mode="numeric", length=length(data_vector))
  resid.data_vector[data.lmRob$x[,2]] <- data.lmRob$residuals 

  resid.data_vector
}

##******************************************************************************************************

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

## Note that Newlyn and Brest have been converted to mb
load("~/Dropbox/brestNewlynData/analysis/paper/brestNewlynEofACRE/brestNewlynACREannualSLP.RData", envir=met)

if(model<13){
  met$nstns <- 37
  nPa <- met$nstns-2 # First two columns are Brest and Newlyn
} else {
#  load("~/Dropbox/brestNewlynData/analysis/paper/brestNewlynEofACRE/brestNewlynACREannualSLP.RData", envir=met) # Brest Newlyn region
  load("~/Dropbox/brestNewlynData/analysis/paper/eofACRE/nAtlACREslpEOF.RData", envir=met) # Whole N Atl
  met$nstns <- nstns # Only consider nstns EOFs 
  nPa <- met$nstns # All columns are relevant
  met$numYrs <- dim(met$dmSlpAnnual)[1]
  Pa <- array(NA, dim=c(met$numYrs,met$nstns))
  # Zero lag
  for(i in 1:met$nstns){
    ## Note that Newlyn and Brest have been converted to mb.
    Pa[,i] <- (met$dmSlpAnnual %*% met$svdDmSlpAnnual$v[,i])
  }
}
met$yrs <- seq(from=1871, to=2008)

pa.start <- which(met$yrs==start) # Start of ACRE pressures for this model
pa.end <- which(met$yrs==end) # End of ACRE pressures for this model

## Use annual pressures over the periods of the model
## In ACRE Brest is column 1 and Newlyn is column 2
## Note that Newlyn and Brest have been converted to mb. Multiply by 100 to return to Pa
if(model!=6){
  newlynPa <-
    met$slpAnnArray[intersect(which(met$yrs>=start),
                                            which(met$yrs<=end)),2]*100
}

brestPa <-
  met$slpAnnArray[intersect(which(met$yrs>=start),
                                            which(met$yrs<=end)),1]*100

time <- c(start:end)

p1p.start <- which(time==p1b.start) # Pressure training data start
p1p.end <- which(time==p1b.end) # Pressure training data end
p2p.start <- which(time==p2b.start) # Pressure prediction start
p2p.end <- which(time==p2b.end) # Pressure prediction end


## Load tg data
tg <- new.env()
load("../tg/tg.RData", envir=tg)

tgbs <- which(tg$brest$Year==start) # index of start of Brest TG data
tgbe <- which(tg$brest$Year==end) # index of end of Brest TG data

if (model!=6){
  tgns <- which(tg$newlyn$Year==start) # index of start of Newlyn TG data
  tgne <- which(tg$newlyn$Year==end) # index of end of Newlyn TG data
}

tg_brest <- tg$brest[tgbs:tgbe,]
if((model==5)||(model==6)||(model>=13)){
  tg_brest_22 <- tg$brest[tgbs:tgbe,]
  tg_brest_22[p2p.start:p2p.end,2] <- tg_brest_22[p2p.start:p2p.end,2]-22 # Subtract 22mm pre-war
}

if (model!=6){
  tg_newlyn <- tg$newlyn[tgns:tgne,]
}
## Calculate the total pressure
brestTotalP <- 1025*9.8*tg_brest[,2]/1000 + brestPa
if((model==5)||(model==6)||(model>=13)){
  brestTotalP_22<- 1025*9.8*tg_brest_22[,2]/1000 + brestPa
}

if (model!=6){
  newlynTotalP <- 1025*9.8*tg_newlyn[,2]/1000 + newlynPa
}

tg.lmRob.brestTotalP <- lmRob(brestTotalP ~ time, control=lmRobControl, x=T)
if((model==5)||(model==6)||(model>=13)){
  tg.lmRob.brestTotalP_22<- lmRob(brestTotalP_22 ~ time, control=lmRobControl, x=T)
}
if (model!=6){
  tg.lmRob.newlynTotalP <- lmRob(newlynTotalP ~ time, control=lmRobControl, x=T)
}
## Zero lag
## Note that Newlyn and Brest have been converted to mb. Multiply by 100 to return to Pa.
if(model<13){
  Pa <- as.array(met$slpAnnArray[pa.start:pa.end,3:met$nstns])*100 # 1st two columns are Brest and Newlyn
} else {
  Pa <- Pa[pa.start:pa.end,]*100
}

##***********************************************************************************************************

##################################
## Main calculations start here ##
##################################

###########
## Brest ##
###########
if(model<13) {
  tg.sd.sorted.brest <- sorted.scale(nPa, brestTotalP, start,end, Pa, lmRobControl)
}else{
  tg.sd.sorted.brest <- sort(c(1:nPa), index=T)
}
## Now construct models of the top 9 components.
## lmRob becomes singular and gives an error with more than 7 components in models 3 and 5
if(model<13) {
  bmax <- 9
} else {
  bmax <- nstns
}

if(bmax>1){
  data.brest <- data.frame(msl=brestTotalP, t=seq(from=start,to=end),
                         Pa[,tg.sd.sorted.brest$ix[1:bmax]])
} else {
  data.brest <- data.frame(msl=brestTotalP, t=seq(from=start,to=end), Pa)
}

tg.lmRob.brest <- lmRob(msl ~ ., x=T, y=T, data=data.brest, control=lmRobControl, na.action=na.exclude)

#################################
## Test model prediction skill ##
#################################
## Need to build new set of top 9 components for partial set
if(model<13) {
  tg.sd.sorted.brest.partial <- sorted.scale(nPa, brestTotalP[p1p.start:p1p.end],
                                           p1b.start,p1b.end, Pa[p1p.start:p1p.end,], lmRobControl)
}else{
  tg.sd.sorted.brest.partial <- sort(c(1:nPa), index=T)
}

if(bmax>1){
  data.brest.partial <- data.frame(msl=brestTotalP[p1p.start:p1p.end],
    t=seq(from=p1b.start,to=p1b.end), Pa[p1p.start:p1p.end,tg.sd.sorted.brest.partial$ix[1:bmax]])
} else {
  data.brest.partial <- data.frame(msl=brestTotalP[p1p.start:p1p.end],
    t=seq(from=p1b.start,to=p1b.end), Pa=Pa[p1p.start:p1p.end])
}

tg.lmRob.brest.partial <- lmRob(msl ~ ., x=T, y=T,
                               data=data.brest.partial, control=lmRobControl, na.action=na.exclude)

if (bmax>1){
  data.brest.new<- data.frame(t=seq(from=p2b.start,to=p2b.end),
                            Pa[p2p.start:p2p.end,tg.sd.sorted.brest.partial$ix[1:bmax]])
} else {
  data.brest.new<- data.frame(t=seq(from=p2b.start,to=p2b.end), Pa=Pa[p2p.start:p2p.end])
}

tg.lmRob.brest.pred<- predict.lmRob(tg.lmRob.brest.partial, se.fit=T,
                                    newdata=data.brest.new, interval='confidence')

if((model==5)||(model==6)||(model>=13)){
##############
## Brest 22 ##
##############
if(model<13) {
  tg.sd.sorted.brest_22 <- sorted.scale(nPa, brestTotalP_22, start,end, Pa, lmRobControl)
}else{
  tg.sd.sorted.brest_22 <- sort(c(1:nPa), index=T)
}

  ## Now construct models of the top 9 components.
  ## lmRob becomes singular and gives an error with more than 7 components in model 5
  if(model<13) {
   bmax <- 9
  } else {
   bmax <- nstns
  }

  if(bmax>1){
    data.brest_22 <- data.frame(msl=brestTotalP_22, t=seq(from=start,to=end),
                              Pa[,tg.sd.sorted.brest$ix[1:bmax]])
  } else {
    data.brest_22 <- data.frame(msl=brestTotalP_22, t=seq(from=start,to=end), Pa)
  }

  tg.lmRob.brest_22 <- lmRob(msl ~ ., x=T, y=T, data=data.brest_22, control=lmRobControl,
                             na.action=na.exclude)

#################################
  ## Test model prediction skill ##
#################################
  ## Need to build new set of top 9 components for partial set
if(model<13) {
  tg.sd.sorted.brest.partial_22 <- sorted.scale(nPa, brestTotalP_22[p1p.start:p1p.end],
                                                p1b.start,p1b.end, Pa[p1p.start:p1p.end,], lmRobControl)
}else{
  tg.sd.sorted.brest.partial_22 <- sort(c(1:nPa), index=T)
}

if(bmax>1){
  data.brest.partial_22 <- data.frame(msl=brestTotalP_22[p1p.start:p1p.end],
                                      t=seq(from=p1b.start,to=p1b.end),
                                      Pa[p1p.start:p1p.end,
                                         as.array(tg.sd.sorted.brest.partial_22$ix)[1:bmax]])
} else {
  data.brest.partial_22 <- data.frame(msl=brestTotalP_22[p1p.start:p1p.end],
                                      t=seq(from=p1b.start,to=p1b.end),
                                      Pa=Pa[p1p.start:p1p.end])
}

  tg.lmRob.brest.partial_22 <- lmRob(msl ~ ., x=T, y=T,
                                     data=data.brest.partial_22, control=lmRobControl,
                                     na.action=na.exclude)

  if(bmax>1){
    data.brest.new_22 <- data.frame(t=seq(from=p2b.start,to=p2b.end),
                                  Pa[p2p.start:p2p.end,
                                     as.array(tg.sd.sorted.brest.partial_22$ix)[1:bmax]])
  } else {
    data.brest.new_22 <- data.frame(t=seq(from=p2b.start,to=p2b.end),
                                  Pa=Pa[p2p.start:p2p.end])
  }

  tg.lmRob.brest.pred_22 <- predict.lmRob(tg.lmRob.brest.partial_22, se.fit=T,
                                          newdata=data.brest.new_22, interval='confidence')
}

if(model!=6){
############
## Newlyn ##
############
if(model<13) {
  tg.sd.sorted.newlyn <- sorted.scale(nPa, newlynTotalP, start,end, Pa, lmRobControl)
}else{
  tg.sd.sorted.newlyn <- sort(c(1:nPa), index=T)
}

## Now construct models of the top 9 components
## lmRob becomes singular and gives an error with more than 8 components in model 5
if(model<13) {
  nmax <- 9
} else {
  nmax <- nstns
}

if(nmax>1){
  data.newlyn <- data.frame(msl=newlynTotalP, t=seq(from=start,to=end),
                          Pa[,tg.sd.sorted.newlyn$ix[1:nmax]])
  } else {
    data.newlyn <- data.frame(msl=newlynTotalP, t=seq(from=start,to=end), Pa)
  }

tg.lmRob.newlyn <- lmRob(msl ~ ., x=T, y=T, data=data.newlyn, control=lmRobControl, na.action=na.exclude)

#################################
## Test model prediction skill ##
#################################
## Need to build new set of top 9 components for partial set
## lmRob becomes singular and gives an error with more than 4 components in model 4
if(model<13) {
  tg.sd.sorted.newlyn.partial <- sorted.scale(nPa, newlynTotalP[p1p.start:p1p.end],
                                            p1n.start,p1n.end, Pa[p1p.start:p1p.end,], lmRobControl)
}else{
  tg.sd.sorted.newlyn.partial <- sort(c(1:nPa), index=T)
}

if(model<13) {
  nmax <- 9
} else {
  nmax <- nstns
}

if(nmax>1){
  data.newlyn.partial <- data.frame(msl=newlynTotalP[p1p.start:p1p.end],
    t=seq(from=p1n.start,to=p1n.end), Pa[p1p.start:p1p.end,
                         tg.sd.sorted.newlyn.partial$ix[1:nmax]])
} else {
  data.newlyn.partial <- data.frame(msl=newlynTotalP[p1p.start:p1p.end],
    t=seq(from=p1n.start,to=p1n.end), Pa=Pa[p1p.start:p1p.end])
}

tg.lmRob.newlyn.partial <- lmRob(msl ~ ., x=T, y=T,
                               data=data.newlyn.partial, control=lmRobControl, na.action=na.exclude)

if(nmax>1){
  data.newlyn.new<- data.frame(t=seq(from=p2n.start,to=p2n.end),
                             Pa[p2p.start:p2p.end,
                                tg.sd.sorted.newlyn.partial$ix[1:nmax]])
} else {
  data.newlyn.new<- data.frame(t=seq(from=p2n.start,to=p2n.end),
                             Pa=Pa[p2p.start:p2p.end])
}

tg.lmRob.newlyn.pred<- predict.lmRob(tg.lmRob.newlyn.partial, se.fit=T,
                                     newdata=data.newlyn.new, interval='confidence')
}
