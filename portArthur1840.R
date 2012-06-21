#
# Script to read in Port Arthur data from the 1840s. The data are high and low
# waters and need to be turned into mean tide levels
#
# Data is in ~simonh/data/portArthur1840/
#

ops <- options()
options("warn" = 1)
library(chron)

daysInMonth <- c(31,28,31,30,31,30,31,31,30,31,30,31)
years <- c(1840,1841,1842)

meanHW <- array(NA, dim=c(12,3))
meanLW <- array(NA, dim=c(12,3))
meanTL <- array(NA, dim=c(12,3))

# Read in data as lists
porta40csv <- read.csv("porta40.csv", fill=TRUE, header=FALSE)
porta41csv <- read.csv("porta41.csv", fill=TRUE, header=FALSE)
porta42csv <- read.csv("porta42.csv", fill=TRUE, header=FALSE)
# Set up data arrays in place of lists
len40 <- length(porta40csv[[1]])
len41 <- length(porta41csv[[1]])
len42 <- length(porta42csv[[1]])

lens <- c(len40,len41,len42)
maxLen <- max(lens)
dataArray <- array(NA, dim=c(maxLen,8,3))
# Fill arrays
# Columns are month, day, hour, min, feet, inches, minutes between readings
for (i in 1:7){
  dataArray[1:len40,i,1] <- porta40csv[[i]]
  dataArray[1:len41,i,2] <- porta41csv[[i]]
  dataArray[1:len42,i,3] <- porta42csv[[i]]
}
# Calculate heights in mm from ft and in and put in column 8
for (i in 1:3) {
  dataArray[,8,i] <- (dataArray[,5,i]*12 + dataArray[,6,i]) * 2.54
}
# Get rid of NAs in dates 
# Look for continuity of data - all days and months should be present and
# there should be no less than 3 and no more than 5 tides
for(k in 1:3) {
  if (leap.year(k+1839)) { 
    daysInMonth[2] <- 29 
  } else {
    daysInMonth[2] <- 28
  } 

  for (i in 1:lens[k]) {
    nextMonth <- dataArray[i,1,k]
    nextDay <- dataArray[i,2,k]
    if(i==1) {
      month <- nextMonth
      day <- nextDay
      lastDay <- nextDay
      tidesInDay <- 1
    }
    if(is.na(nextMonth)) {
      dataArray[i,1,k] <- month
    } else {
      month <- nextMonth
    }
    if(is.na(nextDay)) {
      dataArray[i,2,k] <- day
      tidesInDay <- tidesInDay + 1
    } else {
      if ( i > 1 ) {
        day <- nextDay
# Are days continuous?
        if (day != lastDay + 1) {
          if ((lastDay != daysInMonth[(month - 1)]) && (day != 1)){
            stop("Days are not continuous") 
          }
        }
        lastDay <- day
        if (tidesInDay < 3) { stop("Less than 3 tides in the day") }
        if (tidesInDay == 3) { 
#          warning(paste("3 tides on " , as.character(day), "/", 
#            as.character(month), "/", as.character(k+1839), sep=""))
        }
        if (tidesInDay > 5) { stop("More than 5 tides in the day") }
#        if (tidesInDay == 5) { warning("5 tides in the day") }
        tidesInDay <- 1
      }
    }
  }

# Find HWs and LWs for each month
  for (i in 1:12) {# start in February for 1840
    if ((k==1) && (i==1)) { next }
    month <- which(dataArray[,1,k] == i)
    monthData <- dataArray[month,,k]
# Find all the HWs and LWs in the month
# Define HWs as greater than the mean and LWs as less than the mean
# Keep column 7 as it has information about timing which is needed later
    monthMean <- mean(monthData[,8], na.rm=TRUE)
    HW <- monthData[which(monthData[,8] > monthMean),7:8]
    LW <- monthData[which(monthData[,8] < monthMean),7:8]

# Check HWs and LWs have the same number of values
    lenHW <- length(HW[,1])
    lenLW <- length(LW[,1])
    minLen <- min(lenHW, lenLW)

    if (lenHW > lenLW) {
#      warning(paste("More HWs than LWs in ", as.character(i),
#      "/", as.character(k+1839), sep=""))
    } 
    if (lenHW < lenLW) {
#      warning(paste("More LWs than HWs in ", as.character(i),
#      "/", as.character(k+1839), sep=""))
    }
    HW <- HW[1:minLen,]
    LW <- LW[1:minLen,]
# Delete values where there is less than 200 or more than 600 minutes 
# between readings. Really there should be about 360 minutes
# Check the HW timings first keeping the timings information for looking at
# LWs
    j <- intersect(which(HW[,1]>=200), which(HW[,1]<600))
    if (length(which(HW[,1]>=600))>0) {
      warning(paste("HW timings greater than 600 in ", as.character(i),
      "/", as.character(k+1839), sep=""))
    }
    if (length(which(HW[,1]<=200))>0) {
      warning(paste("HW timings less than 200 in ", as.character(i),
      "/", as.character(k+1839), sep=""))
    }
    HW <- HW[j,]
    LW <- LW[j,]
# Now look at LW timings
    j <- intersect(which(LW[,1]>=200), which(LW[,1]<600))
    if (length(which(LW[,1]>=600))>0) {
      warning(paste("LW timings greater than 600 in ", as.character(i),
      "/", as.character(k+1839), sep=""))
    }
    if (length(which(LW[,1]<=200))>0) {
      warning(paste("LW timings less than 200 in ", as.character(i),
      "/", as.character(k+1839), sep=""))
    }
    HW <- HW[j,]
    LW <- LW[j,]

# At this point timing information can be discarded
    HW <- HW[,2]
    LW <- LW[,2]

    meanHW[i,k] <- mean(HW, na.rm=TRUE)
    meanLW[i,k] <- mean(LW, na.rm=TRUE)

    meanTL[i,k] <- (meanHW[i,k]+meanLW[i,k])/2
  }
}
options(ops)
