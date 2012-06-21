# Script to read in the monthly mean arrays created by makeZetDailyMeans.r 

#######################
# S2 dimensions:      #
#                     #
# Xmin:   -19.83333   #
# Xres:   1/6         #
# Xn:     198         #
#                     #
# Ymin:   40.11111    #
# Yres:   1/9         #
# Yn:     224         #
#######################

# Newlyn 50 06 N  05 33 W => 50.1 -5.55
#> (-5.55-xmin)/xres
#[1] 85.69998
#> (50.1-ymin)/yres
#[1] 89.90001
# Newlyn is at x=86, y=90 
# intersect(which(isea==86),which(jsea==90)) => ip=17701

xmin <- -19.83333
xres <- 1/6
ymin <- 40.11111
yres <- 1/9

# Set up lat and lon variables for later plotting
lon <- seq(from=-19.833333, by=1/6, length=198)
lat <- seq(from=40.1111, by=1/9, length=224)

####################
# Read header info #
####################

# "npseaiseajsea.dat" just holds header information which describes firstly 
# the number of model points ("sea" points since "land" isn't included), then 2 
# vectors which are the i,j of the "sea" points. 
#
# All values are therefore integers
#
# Files are little endian
#

# Open header file
inNpseaIseaJsea <- file("npseaiseajsea.dat", "rb")
# Read first byte of header information
  npsea<-readBin(inNpseaIseaJsea, what="integer", n = 1, size = NA)
# Construct some arrays to read the data into
  isea<-array(NA,dim=npsea)
  jsea<-array(NA,dim=npsea)
# Read the i,j co-ords of the data into these arrays
  isea<-readBin(inNpseaIseaJsea, what="integer", n = npsea, size = NA)
  jsea<-readBin(inNpseaIseaJsea, what="integer", n = npsea, size = NA)
# close the file
close(inNpseaIseaJsea)


#########################
# Read in monthly means #
#########################

# Sea level values are real*4. Each "year" file contains 12 months of data 
# which are stored sequentially, so to read all the data you need to read
# npsea*12 values (= npsea*12*4 bytes)

# Just read first year - 1960
# Open file
inConn <- file("zett.mm.1960.dat", "rb") 
# Read data
zetMonthlyMeans <-readBin(inConn, what="numeric", n = npsea*12, size = 4)
# Close file
close(inConn)

# Reshape the data into a 2D array which is npsea*12 months in size
dim(zetMonthlyMeans) <- c(npsea,12)

# Take month 1, January, and create a 2D array for visualising
# First create a square array of NAs which therefore includes the land values
jan <- array(NA, dim=c(198,224))

# Now we need to plug the values in at the correct co-ordinates
for (ii in 1:npsea){
  jan[isea[ii],jsea[ii]] <- zetMonthlyMeans[ii,1]
}

# Visualise...
library(fields)
image.plot(lon,lat,jan)
