# R function to read in data from a black and white XPM file
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
readXPM <- function(name){

# Open file for writing
  conn <- file(paste(name,".xpm", sep=''), open='r')
# Read header
  readLines(con=conn, n=2)
# What are the dimensions of the array?
  text <- readLines(con=conn, n=1)
# Get rid of leading and end chars
  text <- substr(text, 2, nchar(text, "chars")-2)
  dims <- strsplit(text, " ")
  xn <- as.integer(dims[[1]][1])
  yn <- as.integer(dims[[1]][2])
  numColours <- as.integer(dims[[1]][3])

# Set up data array
  data <- array(NA, dim=c(yn,xn))

# Colours
  colours <- vector(length=numColours, mode="integer")
  characters <- vector(length=numColours, mode="character")
  for (i in 1:numColours){
    text <- readLines(con=conn, n=1)
# Get rid of leading and end chars
    text <- substr(text, 2, nchar(text, "chars")-2)
    characters[i] <- substr(text,1,1)
    colours[i] <- i
  }
# Data
  for (i in 1:(yn-1)){
    text <- readLines(con=conn, n=1)
# Get rid of leading and end chars
    text <- substr(text, 2, nchar(text, "chars")-2)
# Substitute colours and characters
    for (j in 1:numColours){
#      text <- gsub(characters[j], as.character(colours[j]), text) 
      text <- gsub(paste("\\",characters[j], sep=''), 
        as.character(colours[j]), text)
    }
    charVector <- strsplit(text, split=NULL)
    for (j in 1:xn){
      data[i,j] <- as.integer(charVector[[1]][j])
    }
  }
# Last line
#"                                "};
  text <- readLines(con=conn, n=1)
# Get rid of leading and extra end chars
  text <- substr(text, 2, nchar(text, "chars")-3)
# Substitute colours and characters
  for (j in 1:numColours){
#    message(paste(as.character(numColours), characters[j], as.character(colours[j])))
    text <- gsub(paste("\\",characters[j], sep=''), 
      as.character(colours[j]), text)
#    message(text)
  }
  charVector <- strsplit(text, split=NULL)
  for (j in 1:xn){
    data[i,j] <- as.integer(charVector[[1]][j])
  }

# Clean up
  close(conn)

  data
}
