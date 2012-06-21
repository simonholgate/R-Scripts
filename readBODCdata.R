# R Script to read BODC formatted data values

# First 11 lines are header. We can ignore these unless we specify that these 
# should be read.

# Arguments are filename to be read (string) and whether or not to print out
# the header lines (TRUE/FALSE)

readBODCdata <- function(filename, readHeader){

  if (readHeader==TRUE) {
    header <- readLines(filename, n=11)
    print(header)
  }

#  bodc <- read.table(filename, 
#    what=c("character","character","character","character","character"),
#  bodc <- scan(filename, 
  bodc <- read.fwf(filename, c(7,20,12,1,12,1),
    colClasses=c("character","character","numeric","character","numeric","character"),
    col.names=c("Cycle","DateTime","SeaLevel","slFlag","Residual","residFlag"), 
    skip=11, na.strings="-99.0000",  nrows=-1, strip.white=TRUE)

#  print(strptime(bodc$DateTime[2],"%Y/%m/%d %H:%M:%S") 
#    - strptime(bodc$DateTime[1],"%Y/%m/%d %H:%M:%S"))

  Date <- as.POSIXlt(bodc$DateTime)
 
  bodc <- cbind(bodc, Date)

  print(bodc$Date[2] - bodc$Date[1])
  bodc
}
# plot(bodc$Date, bodc$SeaLevel, type='l', col='red')
# good <- which(bodc$slFlag!='M')
# plot(bodc$Date[good], bodc$SeaLevel[good], type='l', col='red')
# plot(bodc$Date[1:(24*4)], bodc$SeaLevel[1:(24*4)], type='l', col='red')
# plot(bodc$Date[1:(24*4*14)], bodc$SeaLevel[1:(24*4*14)], type='l', col='red')

