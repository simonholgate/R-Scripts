########################
# som on sealevel data #
# previously extracted #
########################

library(som)

# Deseasonalised arrays
matts.n.dsmma <- normalize(DSMaskedMonthlyArray)
matts.n.dsma <- normalize(DSMonthlyArray)
# Non-deseasonalised masked array
matts.n.mma <- normalize(maskedMonthlyArray)

#
# De-seasonalised masked monthly mean array
#
foo.dsmma <- som(matts.n.dsmma, xdim=3, ydim=5)

# Get output matrix of values per som map unit
outcode.dsmma<-cbind(lon[isea[j]],lat[jsea[j]],t(foo.dsmma$code))
# Write to file
write.table(outcode.dsmma,"code.dsmma.txt",quote=F,row.names=F) 

# Get the time series of BMUs
bmu.dsmma<-foo.dsmma$visual$x+(foo.dsmma$visual$y*3)
outvisual.dsmma<-cbind(bmu.dsmma,foo.dsmma$visual)
write.table(outvisual.dsmma,"visual.dsmma.txt",quote=F,row.names=F)

# Get the distribution of som units
unit.dsmma<-foo.dsmma$code.sum$x+(foo.dsmma$code.sum$y*3)
outsum.dsmma<-cbind(unit.dsmma,foo.dsmma$code.sum)
write.table(outsum.dsmma,"code_sum.dsmma.txt",quote=F,row.names=F)

#
# De-seasonalised monthly mean array
#
foo.dsma <- som(matts.n.dsma, xdim=3, ydim=5)

# Get output matrix of values per som map unit
outcode.dsma<-cbind(lon[isea],lat[jsea],t(foo.dsma$code))
# Write to file
write.table(outcode.dsma,"code.dsma.txt",quote=F,row.names=F) 

# Get the time series of BMUs
bmu.dsma<-foo.dsma$visual$x+(foo.dsma$visual$y*3)
outvisual.dsma<-cbind(bmu.dsma,foo.dsma$visual)
write.table(outvisual.dsma,"visual.dsma.txt",quote=F,row.names=F)

# Get the distribution of som units
unit.dsma<-foo.dsma$code.sum$x+(foo.dsma$code.sum$y*3)
outsum.dsma<-cbind(unit.dsma,foo.dsma$code.sum)
write.table(outsum.dsma,"code_sum.dsma.txt",quote=F,row.names=F)

#
# Masked monthly mean array
#
foo.mma <- som(matts.n.mma, xdim=3, ydim=5)

# Get output matrix of values per som map unit
outcode.mma<-cbind(lon[isea[j]],lat[jsea[j]],t(foo.mma$code))
# Write to file
write.table(outcode.mma,"code.mma.txt",quote=F,row.names=F) 

# Get the time series of BMUs
bmu.mma<-foo.mma$visual$x+(foo.mma$visual$y*3)
outvisual.mma<-cbind(bmu.mma,foo.mma$visual)
write.table(outvisual.mma,"visual.mma.txt",quote=F,row.names=F)

# Get the distribution of som units
unit.mma<-foo.mma$code.sum$x+(foo.mma$code.sum$y*3)
outsum.mma<-cbind(unit.mma,foo.mma$code.sum)
write.table(outsum.mma,"code_sum.mma.txt",quote=F,row.names=F)
