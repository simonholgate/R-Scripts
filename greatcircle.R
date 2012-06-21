# Script to calculate the great circle distance between two points. We assume
# that the Earth is round with an average radius of 6372.795 km
greatcircle <- function(lat1,lon1,lat2,lon2){
# To calculate the great circle distance we need to know the position of the
# two points on the surface of the sphere. We then need to calculate the
# angular distance between them. 
# Finally, assuming the Earth is spherical we can calculate the arc length.

# First convert degrees to radians
# phi = latitude in radians
# lambda = longitude in radians
  phi1 <- lat1*pi/180
  phi2 <- lat2*pi/180
  lambda1 <- lon1*pi/180
  lambda2 <- lon2*pi/180

# Calculate the anglular between the two points using the formula from
# the Wikipedia, which is more accurate than using the spherical law of
# cosines

# dSigma is the angular difference
  dLambda <- lambda2 - lambda1
  dSigma <- atan2(
    sqrt((cos(phi2)*sin(dLambda))^2 + 
      (cos(phi1)*sin(phi2) - sin(phi1)*cos(phi2)*cos(dLambda))^2)
    /
    (sin(phi1)*sin(phi2)+cos(phi1)*cos(phi2)*cos(dLambda))
    , 1)

# The great circle distance is then r * dSigma, where r is the radius
   r <- 6372795
   gcd <- r*dSigma
   gcd
}
