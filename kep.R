# Script to read in and plot TIRA file from KEP
kep <-  read.fwf("kepData.new.txt.tira.feb.on", widths=c(6,2,5,4,7,8,8,8,8,8),
col.names = c("seq", "flag", "year", "yd", "time", "pr1", "pr2", "temp",
"tide", "resid"), skip=20, colClasses=c("integer", "integer",
"integer","integer", "numeric", "numeric", "numeric", "numeric", "numeric",
"numeric"), na.strings="9999.99", nrows=2900)
#"integer","integer", "numeric", "numeric", "numeric", "numeric", "character",
#"character"), na.strings="9999.99")

len = max(kep$seq)
time <- vector(mode="character", length=len)
for (i in 1:len){
  time[i] <- paste(kep$year[i],kep$yd[i],as.integer(kep$time[i]), 
         as.integer((kep$time[i]%%1)*60))
}

kep$pr1[which(kep$pr1==9999)] <- NA
kep$pr2[which(kep$pr2==9999)] <- NA
kep$resid[which(kep$resid==9999)] <- NA
kep$tide[which(kep$tide==9999)] <- NA

#> which(strptime(time, format("%Y %j %H %M"))==strptime("2008 32 0 0",
# format("%Y %j %H %M")))
#[1] 17868

pr1 <- kep$pr1-mean(kep$pr1, na.rm=TRUE)
pr2 <- kep$pr2-mean(kep$pr2, na.rm=TRUE)
tide <- kep$tide-mean(kep$tide, na.rm=TRUE)
resid <- kep$resid-mean(kep$resid, na.rm=TRUE)

par(family="HersheySans")
#plot(strptime(time[17868:len], format("%Y %j %H %M")),pr1[17868:len], type='l')
#lines(strptime(time[17868:len], format("%Y %j %H %M")),pr2[17868:len], col='red')
plot(strptime(time, format("%Y %j %H %M")),pr1, type='l')
lines(strptime(time, format("%Y %j %H %M")),pr2, col='red')
x11()
par(family="HersheySans")
plot(strptime(time, format("%Y %j %H %M")),tide, col='blue', type='l')
lines(strptime(time, format("%Y %j %H %M")),resid, col='magenta')

