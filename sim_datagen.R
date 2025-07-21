rm(list=ls())

library(fields)
library(sn)
m        <- 25          # Number of spatial locations
s        <- 1:m
mu       <- c(1*s/m+0,  # Temp-obs
              3*s/m+5,  # Precip-obs
              1*s/m+0,  # Temp-model
              2*s/m+1)  # Precip-model    
d        <- rdist(s,s)
Omega    <- exp(-d/5)   # Spatial correlation
rho      <- diag(4)     # Cross correlations
rho[1,2] <- rho[2,1] <- .7  # Cor(temp,prec) for obs
rho[1,3] <- rho[3,1] <- .5
rho[1,4] <- rho[4,1] <- .5
rho[2,3] <- rho[3,2] <- .5
rho[2,4] <- rho[4,2] <- .5
rho[4,3] <- rho[3,4] <- .9 # Cor(temp,prec) for model

skew   <- c(0,100,0,10) # Skewness
df     <- 3             # Degrees of freedom
Smooth <- exp(-d/2)     # Extra smoothing for the model output

Sig    <- kronecker(rho,Omega)
alpha  <- rep(skew,each=m)
dat    <- rmst(n=100, xi=mu, Omega=Sig, alpha=alpha,nu=df)
Temp0  <- dat[,1:m+0*m]
Prec0  <- dat[,1:m+1*m]
Temp1  <- dat[,1:m+2*m]%*%Smooth
Prec1  <- dat[,1:m+3*m]%*%Smooth

par(mfrow=c(2,2))
matplot(s,t(Temp0),type="l",col="gray",main="Temp - obs")
matplot(s,t(Prec0),type="l",col="gray",main="Prec - obs")
matplot(s,t(Temp1),type="l",col="gray",main="Temp - model")
matplot(s,t(Prec1),type="l",col="gray",main="Prec - model")
par(mfrow=c(1,1))
summary(c(Prec1))

Prec0[Prec0<0] = 0
Prec1[Prec1<0] = 0
Prec1 = log(Prec1 + 0.0001)
Prec0 = log(Prec0 + 0.0001)
summary(c(Prec1))
cor(c(Prec0),c(Temp0))
