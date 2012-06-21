# R script to look for 18.6 year nodal cycle in Southampton data
southampton <- read.csv(file="southampton.csv", fill=TRUE)
j <- which(southampton[,1]>1950)
junk <- array(NA,dim=c(length(j),2))
junk[,1] <- southampton[j,1]
junk[,2] <- southampton[j,2]
# Do 5 year running mean on data
len <- length(j)-59
filtered <- vector(mode="numeric", length=len)
filteredYr <- vector(mode="numeric", length=len)
for (i in 1:len){
  filtered[i] <- mean(junk[i:(i+59),2], na.rm=TRUE)
  filteredYr[i] <- mean(junk[i:(i+59),1], na.rm=TRUE)
}
  
