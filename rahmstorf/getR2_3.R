data.1 <- data.frame(tanom=detrendTemps1st[1:53], slr=diffDtCWAnnual1st)
data.2 <- data.frame(tanom=detrendTemps2nd[1:53], slr=diffDtCWAnnual2nd)

# lets bootstrap (get both R2 and coefs)

niter <- 999
myslopes1st <- rep(NA,niter)
myintercepts1st <- rep(NA,niter)
myslopes2nd <- rep(NA,niter)
myintercepts2nd <- rep(NA,niter)

for(i in 1:niter)
{
	idx <- sample(seq(1:53), replace=T)
	
	cal <- data.frame(tanom=detrendTemps1st[idx],slr=diffDtCWAnnual1st[idx])
	
	fit <- lm(slr ~ tanom, data=cal)
	myslope <- coef(fit)[2]
        myintercept <- coef(fit)[1]
        myx <- -myintercept/myslope
#        if ((myx <= -0.2) && (myx >= -0.6)){
	  myslopes1st[i] <- coef(fit)[2]
	  myintercepts1st[i] <- coef(fit)[1]

#        } else {
#	  myslopes1st[i] <- NA
#	  myintercepts1st[i] <- NA
#        }

}

for(i in 1:niter)
{
	idx <- sample(seq(1:53), replace=T)
	
	cal <- data.frame(tanom=detrendTemps2nd[idx],slr=diffDtCWAnnual2nd[idx])
	
	fit <- lm(slr ~ tanom, data=cal)
	myslope <- coef(fit)[2]
        myintercept <- coef(fit)[1]
        myx <- -myintercept/myslope
#        if ((myx <= -0.2) && (myx >= -0.6)){
	  myslopes2nd[i] <- coef(fit)[2]
	  myintercepts2nd[i] <- coef(fit)[1]

#        } else {
#	  myslopes2nd[i] <- NA
#	  myintercepts2nd[i] <- NA
#        }

}

