# Load in csv file and plot it
timecsv45 <- read.csv("time45.csv", fill=TRUE, header=FALSE)

lenCol45 <- length(timecsv45[[1]])
radar45 <- vector(mode="numeric", length=lenCol45)

# Columns are time and radar
time <- strptime(timecsv45[[1]], "%d/%m/%Y %H:%M:%S")
radar45[1:lenCol45] <- timecsv45[[2]]

# Find unique values
uni <- !duplicated(time)
time <- time[uni]
radar45 <- radar45[uni]

lenTime <- length(time[[1]])

plot(time, radar45, type='b', col='blue', xaxt="n", ylim=c(0,5000),
ann=FALSE)
r <- as.POSIXct(round(range(time, na.rm=TRUE), "hours"))
axis.POSIXct(1, at=seq(r[1], r[2], by="hour"), format="%d-%H")
title(xlab="Date (day-hour)", ylab="Sea Level [mm]")
# Check for continuity - difference in times should be 1 minute
# Print out the times when the difference is greater
diffTime <- diff(time)

if(length(which(diffTime>1))>0){
  index <- time[which(diffTime>61)]
  index1 <- time[which(diffTime>61)+1]
  message(as.character(paste(time[index], time[index1]," ",
    diffTime[which(diffTime>61)],"\n")))
  points(time[index], radar45[index],col='blue',pch=19)
  points(time[index1], radar45[index1],col='green',pch=3)
  points(time[index1], radar45[index1],col='green',pch=21)
}
