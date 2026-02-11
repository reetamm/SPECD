rm(list=ls())
library(fields)
library(sn)
library(scales)
sim=0
for(sim in 0:100){
    set.seed(sim)
    m        <- 25          # Number of spatial locations
    s        <- 1:m
    mu       <- c(2*s/m+0,  # Temp-obs
                  3*s/m+5,  # Precip-obs
                  1*s/m+0,  # Temp-model
                  2*s/m+1)  # Precip-model  
    locs <- expand.grid(1:5,1:5)
    d        <- rdist(locs)
    Omega    <- exp(-d/2)   # Spatial correlation
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
    # skep = rep(0,4)
    df     <- 20             # Degrees of freedom
    Smooth <- exp(-d)     # Extra smoothing for the model output
    
    Sig    <- kronecker(rho,Omega)
    alpha  <- rep(skew,each=m)
    dat    <- rmst(n=64*30, xi=mu, Omega=Sig, alpha=alpha,nu=df) #1 month for 64 years
    Temp0  <- dat[,1:m+0*m]
    Prec0  <- dat[,1:m+1*m]
    Temp1  <- dat[,1:m+2*m]%*%Smooth
    Prec1  <- dat[,1:m+3*m]%*%Smooth
    
    crosscors <- matrix(NA,nrow = 25,ncol=2)
    for(i in 1:25){
        crosscors[i,1] <- cor(Temp0[,i],Prec0[,i])
        crosscors[i,2] <- cor(Temp1[,i],Prec1[,i])
    }
    
    plot(crosscors,xlim=c(-0.9,-0.2),ylim=c(-0.9,-0.2))
    abline(0,1)
    
    # par(mfrow=c(2,2))
    # matplot(s,t(Temp0),type="l",col="gray",main="Temp - obs")
    # matplot(s,t(Prec0),type="l",col="gray",main="Prec - obs")
    # matplot(s,t(Temp1),type="l",col="gray",main="Temp - model")
    # matplot(s,t(Prec1),type="l",col="gray",main="Prec - model")
    # par(mfrow=c(1,1))
    
    summary(c(Temp0))
    summary(c(Temp1))
    
    #tails go out too far so bringing them back in a bit
    Temp0 <- rescale(Temp0,to=c(255,285))
    Temp1 <- rescale(Temp1,to=c(250,280))
    
    # Temp0 <- Temp0 + 20
    # Temp1 <- Temp1 + 20
    
    
    summary(c(Prec1))
    summary(c(Prec0))
    
    for(i in 1:25){
        p0q <- quantile(Prec0[,i],0.6) #75% will become 0s
        p1q <- quantile(Prec1[,i],0.3) #50% will become 0s
        
        # scale shifting prcp
        Prec0[,i] <- Prec0[,i] - p0q 
        Prec1[,i] <- Prec1[,i] - p1q 
    }
    
    summary(c(Prec1))
    summary(c(Prec0))
    
    #thresholding to 0
    Prec0[Prec0<0] = 0 
    Prec1[Prec1<0] = 0
    
    pdf(paste0('plots/density_sim.pdf'),width = 8, height = 4)
    par(mfrow=c(1,2),mgp=c(2.25,0.75,0),mar=c(4,4,1,1))
    d0 <-density(Temp0) 
    d1 <-density(Temp1) 
    plotmax.y = max(d0$y,d1$y)
    plotmin.y = min(d0$y,d1$y)
    plotmax.x = max(d0$x,d1$x)
    plotmin.x = min(d0$x,d1$x)
    plot(d1,col=2,ylim=range(c(plotmin.y,plotmax.y)),
         xlim=range(c(plotmin.x,plotmax.x)),ylab="Density",xlab='TMAX',main=paste('TMAX'))
    lines(d0,col=1)
    legend('topleft',c('Mod','Obs'),col=c(2,1),lty = c(1,1),lwd=2)
    
    d0 <-density(log(0.0001+Prec0)) 
    d1 <-density(log(0.0001+Prec1)) 
    plotmax.y = max(d0$y,d1$y)
    plotmin.y = min(d0$y,d1$y)
    plotmax.x = max(d0$x,d1$x)
    plotmin.x = min(d0$x,d1$x)
    plot(d1,col=2,ylim=range(c(plotmin.y,plotmax.y)),
         xlim=range(c(plotmin.x,plotmax.x)),ylab="Density",xlab = 'PRCP',main=paste('PRCP'))
    lines(d0,col=1)
    legend('topright',c('Mod','Obs'),col=c(2,1),lty = c(1,1),lwd=2)
    par(mfrow=c(1,1))
    dev.off()
    # savename <- paste0('data/simdata/',sim,'.RData')
    # save(locs,Temp0,Temp1,Prec0,Prec1,file = savename)
    if(sim==0)
    df <- data.frame(year = rep(1:64, each=30), prcp_obs = c(Prec0), prcp_mod = c(Prec1),
                     tmax_obs = c(Temp0), tmax_mod = c(Temp1),
                     lat = rep(locs[,1],each = 30*64), lon = rep(locs[,2],each = 30*64))
    if(sim>0)
        df <- data.frame(prcp_obs = c(Prec0), prcp_mod = c(Prec1),
                         tmax_obs = c(Temp0), tmax_mod = c(Temp1))
    write.csv(df,file = paste0('data/simdata_csv/',sim,'.csv'))
}
