###############################################################################
## brestNewlynFunctions.R                                                    ##
## Abstracted functions which are common to all models                       ##
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
load("../met/met.RData", envir=met)
nPa <- met$hadnstns-2 # First two columns are Brest and Newlyn

pa.start <- which(met$hadSLP2rAnnualTime==start) # Start of Hadley pressures for this model
pa.end <- which(met$hadSLP2rAnnualTime==end) # End of Hadley pressures for this model

## Use annual pressures over the periods of the model
if(model<6){
  newlynPa <-
    met$hadSLP2rAnnualPressureArray[intersect(which(met$hadSLP2rAnnualTime>=start),
                                            which(met$hadSLP2rAnnualTime<=end)),1]
}

brestPa <-
  met$hadSLP2rAnnualPressureArray[intersect(which(met$hadSLP2rAnnualTime>=start),
                                            which(met$hadSLP2rAnnualTime<=end)),2]

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

if (model<6){
  tgns <- which(tg$newlyn$Year==start) # index of start of Newlyn TG data
  tgne <- which(tg$newlyn$Year==end) # index of end of Newlyn TG data
}

tg_brest <- tg$brest[tgbs:tgbe,]
if((model==5)||(model==6)){
  tg_brest_22 <- tg$brest[tgbs:tgbe,]
  tg_brest_22[p2p.start:p2p.end,2] <- tg_brest_22[p2p.start:p2p.end,2]-22 # Subtract 22mm pre-war
}

if (model<6){
  tg_newlyn <- tg$newlyn[tgns:tgne,]
}
## Calculate the total pressure
brestTotalP <- 1025*9.8*tg_brest[,2]/1000 + brestPa
if((model==5)||(model==6)){
  brestTotalP_22<- 1025*9.8*tg_brest_22[,2]/1000 + brestPa
}

if (model<6){
  newlynTotalP <- 1025*9.8*tg_newlyn[,2]/1000 + newlynPa
}

tg.lmRob.brestTotalP <- lmRob(brestTotalP ~ time, control=lmRobControl, x=T)
if((model==5)||(model==6)){
  tg.lmRob.brestTotalP_22<- lmRob(brestTotalP_22 ~ time, control=lmRobControl, x=T)
}
if (model<6){
  tg.lmRob.newlynTotalP <- lmRob(newlynTotalP ~ time, control=lmRobControl, x=T)
}
## Zero lag
Pa <- as.array(met$hadSLP2rAnnualPressureArray[pa.start:pa.end,3:37]) # 1st two columns are Brest and Newlyn

##***********************************************************************************************************

##################################
## Main calculations start here ##
##################################

###########
## Brest ##
###########

tg.sd.sorted.brest <- sorted.scale(nPa, brestTotalP, start,end, Pa, lmRobControl)

## Now construct models of the top 9 components.
## lmRob becomes singular and gives an error with more than 7 components in models 3 and 5
if((model==3)||(model==5)||(model==6)) {
  bmax <- 7
} else {
  bmax <- 9
}

data.brest <- data.frame(msl=brestTotalP, t=seq(from=start,to=end), Pa[,tg.sd.sorted.brest$ix[1:bmax]])

tg.lmRob.brest <- lmRob(msl ~ ., x=T, y=T, data=data.brest, control=lmRobControl, na.action=na.exclude)

#################################
## Test model prediction skill ##
#################################
## Need to build new set of top 9 components for partial set

tg.sd.sorted.brest.partial <- sorted.scale(nPa, brestTotalP[p1p.start:p1p.end],
                                           p1b.start,p1b.end, Pa[p1p.start:p1p.end,], lmRobControl)

data.brest.partial <- data.frame(msl=brestTotalP[p1p.start:p1p.end],
  t=seq(from=p1b.start,to=p1b.end), Pa[p1p.start:p1p.end,tg.sd.sorted.brest.partial$ix[1:bmax]])

tg.lmRob.brest.partial <- lmRob(msl ~ ., x=T, y=T,
                               data=data.brest.partial, control=lmRobControl, na.action=na.exclude)

data.brest.new<- data.frame(t=seq(from=p2b.start,to=p2b.end),
                            Pa[p2p.start:p2p.end,tg.sd.sorted.brest.partial$ix[1:9]])

tg.lmRob.brest.pred<- predict.lmRob(tg.lmRob.brest.partial, se.fit=T,
                                    newdata=data.brest.new, interval='confidence')

if((model==5)||(model==6)){
##############
## Brest 22 ##
##############

  tg.sd.sorted.brest_22 <- sorted.scale(nPa, brestTotalP_22, start,end, Pa, lmRobControl)

  ## Now construct models of the top 9 components.
  ## lmRob becomes singular and gives an error with more than 7 components in model 5
  if((model==5)||(model==6)) {
   bmax <- 7
  } else {
   bmax <- 9
  }

  data.brest_22 <- data.frame(msl=brestTotalP_22, t=seq(from=start,to=end),
                              Pa[,tg.sd.sorted.brest$ix[1:bmax]])

  tg.lmRob.brest_22 <- lmRob(msl ~ ., x=T, y=T, data=data.brest_22, control=lmRobControl,
                             na.action=na.exclude)

#################################
  ## Test model prediction skill ##
#################################
  ## Need to build new set of top 9 components for partial set

  tg.sd.sorted.brest.partial_22 <- sorted.scale(nPa, brestTotalP_22[p1p.start:p1p.end],
                                                p1b.start,p1b.end, Pa[p1p.start:p1p.end,], lmRobControl)

  data.brest.partial_22 <- data.frame(msl=brestTotalP_22[p1p.start:p1p.end],
                                      t=seq(from=p1b.start,to=p1b.end),
                                      Pa[p1p.start:p1p.end,tg.sd.sorted.brest.partial_22$ix[1:bmax]])

  tg.lmRob.brest.partial_22 <- lmRob(msl ~ ., x=T, y=T,
                                     data=data.brest.partial_22, control=lmRobControl,
                                     na.action=na.exclude)

  data.brest.new_22 <- data.frame(t=seq(from=p2b.start,to=p2b.end),
                                  Pa[p2p.start:p2p.end,tg.sd.sorted.brest.partial_22$ix[1:9]])

  tg.lmRob.brest.pred_22 <- predict.lmRob(tg.lmRob.brest.partial_22, se.fit=T,
                                          newdata=data.brest.new_22, interval='confidence')
}

if(model<6){
############
## Newlyn ##
############

tg.sd.sorted.newlyn <- sorted.scale(nPa, newlynTotalP, start,end, Pa, lmRobControl)

## Now construct models of the top 9 components
## lmRob becomes singular and gives an error with more than 8 components in model 5
if(model==5) {
  nmax <- 8
} else {
  nmax <- 9
}
data.newlyn <- data.frame(msl=newlynTotalP, t=seq(from=start,to=end), Pa[,tg.sd.sorted.newlyn$ix[1:nmax]])

tg.lmRob.newlyn <- lmRob(msl ~ ., x=T, y=T, data=data.newlyn, control=lmRobControl, na.action=na.exclude)

#################################
## Test model prediction skill ##
#################################
## Need to build new set of top 9 components for partial set
## lmRob becomes singular and gives an error with more than 4 components in model 4

tg.sd.sorted.newlyn.partial <- sorted.scale(nPa, newlynTotalP[p1p.start:p1p.end],
                                            p1n.start,p1n.end, Pa[p1p.start:p1p.end,], lmRobControl)

if(model==4) {
  nmax <- 4
} else {
  nmax <- 9
}

data.newlyn.partial <- data.frame(msl=newlynTotalP[p1p.start:p1p.end],
  t=seq(from=p1n.start,to=p1n.end), Pa[p1p.start:p1p.end,tg.sd.sorted.newlyn.partial$ix[1:nmax]])

tg.lmRob.newlyn.partial <- lmRob(msl ~ ., x=T, y=T,
                               data=data.newlyn.partial, control=lmRobControl, na.action=na.exclude)

data.newlyn.new<- data.frame(t=seq(from=p2n.start,to=p2n.end),
                             Pa[p2p.start:p2p.end,tg.sd.sorted.newlyn.partial$ix[1:9]])

tg.lmRob.newlyn.pred<- predict.lmRob(tg.lmRob.newlyn.partial, se.fit=T,
                                     newdata=data.newlyn.new, interval='confidence')
}
