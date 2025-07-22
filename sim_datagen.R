rm(list=ls())
library(fields)
library(sn)
library(scales)
set.seed(0)
m        <- 25          # Number of spatial locations
s        <- 1:m
mu       <- c(2*s/m+0,  # Temp-obs
              3*s/m+5,  # Precip-obs
              1*s/m+0,  # Temp-model
              2*s/m+1)  # Precip-model  
locs <- expand.grid(1:5,1:5)
d        <- rdist(locs)
Omega    <- exp(-d/5)   # Spatial correlation
rho      <- diag(4)     # Cross correlations
rho[1,2] <- rho[2,1] <- .8  # Cor(temp,prec) for obs
rho[1,3] <- rho[3,1] <- .5
rho[1,4] <- rho[4,1] <- .5
rho[2,3] <- rho[3,2] <- .5
rho[2,4] <- rho[4,2] <- .5
rho[4,3] <- rho[3,4] <- .4 # Cor(temp,prec) for model
ss   <- diag(c(1,-1,1,-1))
rho <- ss%*%rho%*%ss

skew   <- c(0,100,0,10) # Skewness
df     <- 3             # Degrees of freedom
Smooth <- exp(-d/2)     # Extra smoothing for the model output

Sig    <- kronecker(rho,Omega)
alpha  <- rep(skew,each=m)
dat    <- rmst(n=50*30, xi=mu, Omega=Sig, alpha=alpha,nu=df) #1 month for 50 years
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

summary(c(Temp0))
summary(c(Temp1))

#tails go out too far so bringing them back in a bit
t0q <- quantile(c(Temp0),c(0.05,0.999)) 
t1q <- quantile(c(Temp1),c(0.05,0.99))
Temp0 <- rescale(Temp0,to=t0q)
Temp1 <- rescale(Temp1,to=t1q)

Temp0 <- Temp0 + 20
Temp1 <- Temp1 + 20

plot(density(Temp0),xlim=c(15,35))
lines(density(Temp1))

summary(c(Prec1))
summary(c(Prec0))
p0q <- quantile(c(Prec0),0.25) #25% will become 0s
p1q <- quantile(c(Prec1),0.5) #50% will become 0s

# scale shifting prcp
Prec0 <- Prec0 - p0q 
Prec1 <- Prec1 - p1q
summary(c(Prec1))
summary(c(Prec0))

#thresholding to 0
Prec0[Prec0<0] = 0 
Prec1[Prec1<0] = 0

save(locs,Temp0,Temp1,Prec0,Prec1,file = 'data/simdata.RData')

#log transform like we normally do
Prec1 = log(Prec1 + 0.0001)
Prec0 = log(Prec0 + 0.0001)

cor(c(Prec0),c(Temp0)) #cross correlation for obs
cor(c(Prec1),c(Temp1)) #cross correlation for model

crosscors <- matrix(NA,25,2)
for(i in 1:25){
    crosscors[i,1] <- cor(Temp0[,i],Prec0[,i])
    crosscors[i,2] <- cor(Temp1[,i],Prec1[,i])
}

plot(crosscors[,2:1],xlim=c(-0.6,0),ylim=c(-0.6,0))
abline(0,1)

plot(density(Prec0), ylim=c(0,0.4),xlim = c(-15,6)) #obs prcp (more 0s)
lines(density(Prec1)) #model prcp (fewer 0s)
