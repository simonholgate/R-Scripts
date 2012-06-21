specificHeat <- function(s,t,p){
# Function to calculate the specific heat of sea water in j/kg/K
# at salinity, s, temperature, t and pressure, p=0 (in SI units)
# From Gill Appendix 3 based on Millero et al (1973)
  cp0t0 <- 4217.4 - 3.720283*t + 0.1412855 * t^2 - 2.654387e-3 * t^3.0 
           + 2.093236e-5 * t^4.0
  cpst0 <- cp0t0 + s * (-7.6444 + 0.107276 * t - 1.3839e-3 * t^2.0)
           + s^(3.0/2.0) * (0.17709 - 4.0772e-3 * t + 5.3539e-5 * t^2.0)
  return(cpst0)
}
