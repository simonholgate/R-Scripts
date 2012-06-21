# Load in csv file and plot it
timecsv <- read.csv("time.csv", fill=TRUE, header=FALSE)
timecsv46 <- read.csv("time46.csv", fill=TRUE, header=FALSE)
timecsv45 <- read.csv("time45.csv", fill=TRUE, header=FALSE)

lenCol <- length(timecsv[[1]])
radar <- vector(mode="numeric", length=lenCol)

lenCol46 <- length(timecsv46[[1]])
radar46 <- vector(mode="numeric", length=lenCol46)
pressure46 <- vector(mode="numeric", length=lenCol46)

lenCol45 <- length(timecsv45[[1]])
radar45 <- vector(mode="numeric", length=lenCol45)

# Columns are time and radar
time <- c(strptime(timecsv[[1]], "%d/%m/%Y %H:%M:%S"), 
  strptime(timecsv45[[1]], "%d/%m/%Y %H:%M:%S"))
radar[1:(lenCol+lenCol45)] <- c(timecsv[[2]], timecsv45[[2]])

time46 <- strptime(timecsv46[,1], "%d/%m/%Y %H:%M:%S")
radar46[1:lenCol46] <- timecsv46[[2]]
pressure46[1:lenCol46] <- timecsv46[[6]]

# Find unique values
uni <- !duplicated(time)
time <- time[uni]
radar <- radar[uni]

uni46 <- !duplicated(time46)
time46 <- time46[uni46]
radar46 <- radar46[uni46]
pressure46 <- pressure46[uni46]

# Find where we go to 1 minute sampling
min1 <- which(time==strptime("2007-09-03 15:07:00", "%Y-%m-%d %H:%M:%S"))
lenTime <- length(time[[1]])
min1 <- c(min1:lenTime)
# convert height in fresh water to height in sea water by multiplying by
# 0.9756
height46 <- pressure46*0.9756098
# We need to invert the radar to put it on the same basis as the pressure
# reading
z <- radar46[324]+height46[324]

plot(time[min1], radar[min1], type='b', col='blue', xaxt="n", ylim=c(0,10000),
ann=FALSE)
points(time46, z - radar46, type='b', col='red', xaxt="n", ylim=c(2,5000))
points(time46, height46, type='b', col='magenta')
#r <- as.POSIXct(round(range(time[min1], na.rm=TRUE), "hours"))
r <- as.POSIXct(round(range(time, na.rm=TRUE), "hours"))
#axis.POSIXct(1, at=seq(r[1], r[2], by="day"), format="%d/%m")
axis.POSIXct(1, at=seq(r[1], r[2], by="hour"), format="%d-%H")
title(xlab="Date (day-hour)", ylab="Sea Level [mm]")
legend(x=time[min1[275]], y=10000, legend=c('Radar','Pressure','Synthetic'),
  col=c('red','magenta','blue'), lwd=2, pch=21)
# Check for continuity - difference in times should be 1 minute
# Print out the times when the difference is greater
diffTime <- diff(time[min1])
index <- min1[which(diffTime>61)]
index1 <- min1[which(diffTime>61)+1]
message(as.character(paste(time[index], time[index1]," ", diffTime[which(diffTime>61)],"\n")))
points(time[index], radar[index],col='blue',pch=19)
points(time[index1], radar[index1],col='green',pch=3)
points(time[index1], radar[index1],col='green',pch=21)
