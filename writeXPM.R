# R function to write out data as a black and white XPM file
# Header looks like this:
#/* XPM */
#static char * TkBW_xpm[] = {
#/* width height num_colors chars_per_pixel */
#"32 32 2 1",
#/* colors */
#"       c #000000",
#".      c #FFFFFF",
# Followed by (in this case) 32 lines of data like:
#/* pixels */
#"          ........              ",
# and ending:
#"                                "};
# Note that comments are C style
#
writeXPM <- function(name,data){

# What are the dimensions of the array?
  yn <- dim(data)[1]
  xn <- dim(data)[2]
# Open file for writing
  conn <- file(paste(name,".xpm", sep=''), open='w')
# Write header
  writeLines("/* XPM */",con=conn)
  writeLines(paste("static char * ",name,"_xpm[] = {", sep=""),con=conn)
  writeLines("/* width height num_colors chars_per_pixel */",con=conn)
# Write dimensions of mask - assumption is that we are using black and white
# and only one character per pixel
  writeLines(paste('"', xn,' ', yn,' 2 1"', sep=""), con=conn)
# Colours
  writeLines("/* colors */", con=conn)
  writeLines('"0      c #000000",', con=conn)
  writeLines('"1      c #FFFFFF",', con=conn)
# Data
  writeLines('#/* pixels */', con=conn)
  for (i in 1:(yn-1)){
    dataLine <- data[i,]
# Sea is finite and is converted white
    index <- which(is.finite(dataLine))
    dataLine[index] <- 0
# Land is NA and is converted black
    index <- which(is.na(dataLine))
    dataLine[index] <- 1
    dataLine <- paste(dataLine, sep="", collapse="")
    writeLines(paste('"', dataLine,'"', sep="", collapse=""), con=conn)
  }
# Last line
#"                                "};
  dataLine <- data[yn,]
# Sea is finite and is converted white
  index <- which(is.finite(dataLine))
  dataLine[index] <- 0
# Land is NaN and is converted black
  index <- which(is.na(dataLine))
  dataLine[index] <- 1
  dataLine <- paste(dataLine, sep="", collapse="")
  writeLines(paste('"', dataLine,'"};', sep="", collapse=""), con=conn)

# Clean up
  close(conn)

}
