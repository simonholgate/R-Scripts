# First get new station ids to replace country/station codes
# filelist to convert between systems is in downloaded data
#
rlridlist <- read.table(file="~/data/rlr_data/140910/rlr_annual/filelist.txt", sep=";",
col.names=c("id", "slatdec", "slondec", "sname", "ccode", "scode", "metric"),
colClasses=c("integer","numeric","numeric","character","character","character","character"),
quote=NULL)
rlriddata <- as.numeric(rlridlist$ccode)+(as.numeric(rlridlist$scode)/1000)

# Dutch data (150.xxx) is in metric files
metidlist <- read.table(file="~/data/rlr_data/140910/met_monthly/filelist.txt", sep=";",
col.names=c("id", "slatdec", "slondec", "sname", "ccode", "scode", "metric"),
colClasses=c("integer","numeric","numeric","character","character","character","character"),
quote=NULL)
metiddata <- as.numeric(metidlist$ccode)+(as.numeric(metidlist$scode)/1000)

s177ids <- array(NA, dim=c(177,2))
for (i in 1:177){
  if (floor(dataArray[i,1]) != 150){
    s177ids[i,1] <- rlridlist$id[which(rlriddata == dataArray[i,1])]
    s177ids[i,2] <- 1
  } else {
    s177ids[i,1] <- metidlist$id[which(metiddata == dataArray[i,1])]
    s177ids[i,2] <- 0
  }
}

