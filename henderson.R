# Code for calculation of Henderson moving average filter
# See "An Introductory Course on Time Series Analysis" by the Australian
# Bureau of Statistics (PDF is saved as ~simonh/doc/TSAnalysisCourse.pdf) page
# 43ff

# Usage: 
# henderson(N) 
# for an N point Henderson filter 

henderson <-
function(n) {
  r <- (n - 1)/2

  C0 <- (315/8)*(3*r^2 + 12*r - 4)/(64*r^9 + 1152*r^8 + 8592*r^7 + 34272*r^6 +
  78204*r^5 + 99288*r^4 + 57203*r^3 - 4302*r^2 - 19323*r - 5670)
  C1 <- (3465/8)*-1/(64*r^9 + 1152*r^8 + 8592*r^7 + 34272*r^6 +
  78204*r^5 + 99288*r^4 + 57203*r^3 - 4302*r^2 - 19323*r - 5670)

# Hard code for test of n=5
#  C0 <- 1.5617323E-4
#  C1 <- -5.4836106E-5


# Solve for coefficients using matrices as suggested by CW
#  Xa <- 0
#  Xb <- 0
#  Yb <- 0
#  for (i in -r:r){
#    Xa <- Xa + 
#      ((r + 1)^2 - i^2)*((r + 2)^2 - i^2)*((r + 3)^2 - i^2)
#    Xb <- Xb + 
#      ((r + 1)^2 - i^2)*((r + 2)^2 - i^2)*((r + 3)^2 - i^2)*i^2
#    Ya <- Xb
#    Yb <- Yb + 
#      ((r + 1)^2 - i^2)*((r + 2)^2 - i^2)*((r + 3)^2 - i^2)*i^4
#  }
#  Mat1 <- matrix(c(Xa,Xb, Ya,Yb), nrow = 2, ncol=2, byrow=TRUE)
#  Mat2 <- matrix(c(1,0), nrow=2)
#  C <- solve(Mat1,Mat2)

# wk are the weights from -r to r
  wk <- vector(mode="numeric", length=n)
  for (i in -r:r){
    wk[i+(r+1)] <- 
      ((r + 1)^2 - i^2)*((r + 2)^2 - i^2)*((r + 3)^2 - i^2)*(C0 + C1*i^2)
#      ((r + 1)^2 - i^2)*((r + 2)^2 - i^2)*((r + 3)^2 - i^2)*(C[1] + C[2]*i^2)
  }

# Return the weights to the calling function
  wk
}
