###################################
# som on sea level altimetry data #
# previously extracted in         #
# read_altim_nweurope.R           #
###################################

library(som)

# Deseasonalised arrays
matts.n <- normalize(smNWSOMDataArray)

# Calculate SOM
foo <- som(matts.n, xdim=3, ydim=4)

# Get output matrix of values per som map unit
outcode <-cbind(nwSeaPointArray[,1],nwSeaPointArray[,2],t(foo$code))
# Write to file
write.table(outcode,"code.txt",quote=F,row.names=F) 

# Get the time series of BMUs
bmu<-foo$visual$x+(foo$visual$y*3)
outvisual<-cbind(bmu,foo$visual)
write.table(outvisual,"visual.txt",quote=F,row.names=F)

# Get the distribution of som units
unit<-foo$code.sum$x+(foo$code.sum$y*3)
outsum<-cbind(unit,foo$code.sum)
write.table(outsum,"code_sum.txt",quote=F,row.names=F)

