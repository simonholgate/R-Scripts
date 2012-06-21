# Code to reproduce the Jet colour palette from Matlab
#
# Copying the numbers from Matlab default
#
#jet.nums <- read.table("jet",col.names=c("h","s","v"),sep="")
#jet.colors <- rgb(jet.nums$h,jet.nums$s,jet.nums$v)
#
# Done as a function
jet.colors <-
function (n)
{
# n must be >= 8 to work
    if ((n <- as.integer(n[1])) >= 8) {
        k <- ceiling(c(1,3,5,7)*n/8)
        r <- array(c(1:n)*0,c(1,n))
        g <- array(c(1:n)*0,c(1,n))
        b <- array(c(1:n)*0,c(1,n))

        r[k[2]:k[3]] <- seq(0,1,length=length(k[2]:k[3]))
        r[k[3]:k[4]] <- 1
        r[k[4]:n] <- seq(1,0.5,length=length(k[4]:n))
        g[k[1]:(k[2]+1)] <- seq(0,1,length=length(k[1]:(k[2]+1)))
        g[(k[2]+1):k[3]] <- 1
        g[k[3]:(k[4]+1)] <- seq(1,0,length=length(k[1]:(k[2]+1)))
        b[1:k[1]] <- seq(0.5,1,length=length(1:k[1]))
        b[k[1]:k[2]] <- 1
        b[k[2]:k[3]] <- seq(1,0,length=length(k[2]:k[3]))

        c(rgb(r,g,b))
    }
    else character(0)
}


