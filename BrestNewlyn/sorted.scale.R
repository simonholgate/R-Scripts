###############################
## Locally defined functions ##
###############################

sorted.scale <- function(n, totalP, start, end, pressure, incontrol){
  ## What we are going to do now is calculate the sd of the residuals as each of the
  ## components is added in turn

  require(robust)
  tg.sd <- vector(mode="numeric", length=n)

  for(i in 1:n){
    data.tmp <- data.frame(msl=totalP, t=seq(from=start, to=end), Pa=pressure[,i])
    tg.lmRob <- lmRob(msl ~ t + Pa, data=data.tmp, control=incontrol, na.action=na.exclude)
    tg.sd[i] <- tg.lmRob$scale
  }

  ## Now we can sort these sds. Thompson took the first 9 and we will do the same for
  ## comparison. 
  tg.sd.sorted <- sort(tg.sd, index=TRUE)
  tg.sd.sorted
}
