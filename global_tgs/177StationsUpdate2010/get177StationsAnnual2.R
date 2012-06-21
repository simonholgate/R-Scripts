# Use downloaded data here rather than Oracle
# Uses station IDs from downloaded data using new database
# Sept-2010

stns177Annual <- array(NA,dim=c(62,177))
stns177AnnualYrs <- array(data=1948:2009,dim=c(62,1))
for (i in 1:177){
  # Loop through each station id and select annual data file if RLR or monthly data file if Dutch data
  if (s177ids[i,2]==1){
    stns177File <- paste("~/data/rlr_data/140910/rlr_annual/data/",as.character(s177ids[i,1]),".rlrdata", sep="")
    dataFile <- read.table(file=stns177File, sep=";", col.names=c("year", "height", "interp", "flag"),
      colClasses=c("integer","integer","character","integer"), na.strings="-99999")
    numYears <- length(dataFile$year)

    for(j in 1:numYears){
      index <- which(stns177AnnualYrs==dataFile$year[j])
      stns177Annual[index,i] <- dataFile$height[j]
    }
  }else{# Dutch stations
    stns177File <- paste("~/data/rlr_data/140910/met_monthly/data/",as.character(s177ids[i,1]),".metdata", sep="")
    dataFile <- read.table(file=stns177File, sep=";", col.names=c("year", "height", "interp", "flag"),
      colClasses=c("numeric","integer","character","integer"), na.strings="-99999")
    numMonths <- length(dataFile$year)
    
    # We have a problem that the monthly files may not be a whole number of years long.
    # Missing months are at the beginning or end.
    # We'll get around this by putting the data into an array that is a whole number of years long,
    # lining it up and reshaping it
    numYears <- ceiling(numMonths/12)
    
    monthArray <- as.array(dataFile$height)
    yearArray <- array(NA,dim=c(numYears*12,1))
    # Line up monthArray into yearArray
    # Are any months missing?
    if(identical(all.equal((dataFile$year[1] %% 1), (0.5/12), tolerance=0.01), TRUE)){# No months missing
      yearArray[1:numMonths,1] <- monthArray
    } else {# How many months are missing?
      firstMonth <-  round((dataFile$year[1] %% 1)*12+0.5)
      yearArray[firstMonth:(numYears*12),1] <- monthArray
    }
    
    dim(yearArray) <- c(12,numYears)

    # Calculate annual means.
    annualMeans <- colMeans(yearArray, na.rm=T)
    annualYrs <- c(floor(dataFile$year[1]):floor(dataFile$year[numMonths]))
    
    for(j in 1:numYears){
      # Require at least 10 months in annual mean
      if(length(which(is.finite(yearArray[,j]))) < 10) annualMeans[j] <- NA
      
      index <- which(stns177AnnualYrs==annualYrs[j])
      stns177Annual[index,i] <- annualMeans[j]
    }
  }
}
