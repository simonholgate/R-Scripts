# Write out mask_ew in a FORTRAN readable way
# Integers are natively real*8 in R
outConn <- file("mask.dat", "wb")
# Write first integer that Fortran needs which is the number of bytes in the
# record
numBytes <- 503*533*8
writeBin(numBytes,outConn, endian = "little")
# Replace NAs with missing value of 0
fortranMask <- mask_ew
fortranMask[which(is.na(fortranMask))] <- 0
fortranMask <- as.vector(fortranMask)
# Write mask
writeBin(fortranMask, outConn, endian = "little")
# Write last integer that Fortran needs
writeBin(numBytes,outConn, endian = "little")
close(outConn)
