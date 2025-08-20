rm(list = ls())
library(GpGp)
library(SPQR)
library(lubridate)
load(file = 'data/simdata.RData')

coords = as.matrix(locs)
head(coords)

set.seed(303)
vecchia.order = order_maxmin(coords,lonlat = F)
NNarray <- find_ordered_nn(coords[vecchia.order,],lonlat = F,m=5)

Temp0 <- Temp0[,vecchia.order]
Temp1 <- Temp1[,vecchia.order]
Prec0 <- Prec0[,vecchia.order]
Prec1 <- Prec1[,vecchia.order]
loc = 1

for(loc in 1:6){
    pdfname     <- paste0('plots/sim/fits_l',loc,'.pdf')
    predname1   <- paste0( 'fits/sim/fits_temp_l',loc,'.RDS')
    predname2   <- paste0( 'fits/sim/fits_prcp_l',loc,'.RDS')
    
    # current.loc = vecchia.order[loc]
    
    y1 <- c(Temp0[,loc],Temp1[,loc])
    y2 <- c(Prec0[,loc],Prec1[,loc])
    y2 <- log(0.0001+y2)
    n0 <- n1 <- nrow(Temp0)
    n = n0 + n1
    y0 <- rep(1:0,each=n0)
    
    # plot(y1,y2,col=y0+1)
    
    X1 = matrix(y0, ncol=1)
    
    if(loc>1){
        nns <- NNarray[loc,] # select the correct row
        nns <- nns[complete.cases(nns)] # drop the NAs
        nns <- nns[-1] # drop the response
        for(k in nns){
            # vlocs = vecchia.order[k]
            # x.vec = c(Temp0[,vlocs],Temp1[,vlocs])
            x.vec = c(Temp0[,k],Temp1[,k])
            X1 = cbind(X1,x.vec)
        }    
    }
    
    # head(X1)

    control <- list(iter = 300, batch.size = 100, lr = 0.001, save.name = paste('SPQR.model.temp',loc,'pt',sep='.'))
    fit.y1.mle <- SPQR(X = X1, Y = y1, method = "MLE", control = control, normalize = T, verbose = T,use.GPU=F,
                          n.hidden = c(30,20), activation = 'relu',n.knots = 20,seed = loc)

    # plotGOF(fit.y1.mle)
    cdf.y1.mle = rep(NA,n)
    for(i in 1:n){
        cdf.y1.mle[i] <- predict(fit.y1.mle,   X = X1[i,], Y=y1[i], type = "CDF")    
        if(i%%1000==0)
            print(i)
    }
    
    qout11 <- cdf.y1.mle
    adjust = which(qout11>0.99999)
    qout11[adjust] = 0.99999
    
    
    ###################################
    ###################################
    X2 = cbind(X1,y1)
    nx1 = ncol(X1)+1
    if(loc>1){
        for(k in nns){
            # vlocs = vecchia.order[k]
            # x.vec = c(Prec0[,vlocs],Prec1[,vlocs])
            x.vec = c(Prec0[,k],Prec1[,k])
            x.vec = log(0.0001+x.vec)
            X2 = cbind(X2,x.vec)
        }    
    }
    nx2 = ncol(X2)
    # head(X2)
    # head(y2)
    control <- list(iter = 300, batch.size = 100, lr = 0.001, save.name = paste('SPQR.model.prcp',loc,'pt',sep='.'))
    fit.y2.mle <- SPQR(X = X2, Y = y2, method = "MLE", control = control, normalize = T, verbose = T,use.GPU=F,
                          n.hidden = c(30,20), activation = 'relu',n.knots = 20,seed = loc)
    # plotGOF(fit.y2.mle)
    cdf.y2.mle = rep(NA,n)
    for(i in 1:n){
        cdf.y2.mle[i] <- predict(fit.y2.mle,   X = X2[i,], Y=y2[i], type = "CDF")   
        if(i%%1000==0)
            print(i)
    }
    
    qout21 <- cdf.y2.mle
    adjust = which(qout21>0.99999)
    qout21[adjust] = 0.99999
    
    
############### Predictions temp     
    
    if(loc==1){
        qf.y1.mle = rep(NA,n)
        for(i in 1:n){
            x_pred = 1
            qf.y1.mle[i] <- predict(fit.y1.mle,   X = x_pred, type = "QF",tau=qout11[i])
        }    
    }
    
    if(loc>1){
        X1_pred = X1[,1]
        for(k in nns){
            vecname = paste0('fits/sim/fits_temp_l',k,'.RDS')
            x.vec = readRDS(vecname)
            X1_pred = cbind(X1_pred,x.vec)
        }
        X1_pred[,1] = 1
        qf.y1.mle = rep(NA,n)
        for(i in 1:n){
            if(i%%1000==0)
                print(i)
            qf.y1.mle[i] <- predict(fit.y1.mle,   X = X1_pred[i,], type = "QF",tau=qout11[i])
        }   
    }
    
    saveRDS(qf.y1.mle,file = predname1)
############### Predictions prcp        
    if(loc==1){
        qf.y2.mle = rep(NA,n)
        for(i in 1:n){
            if(i%%1000==0)
                print(i)
            if(i<=n0+1)
                x_pred = c(1,y1[i])
            if(i>n0+1)
                x_pred = c(1,qf.y1.mle[i])
            qf.y2.mle[i] <- predict(fit.y2.mle,   X = x_pred, type = "QF",tau=qout21[i])
        } 
    }
    
    if(loc>1){
        X2_pred = cbind(X1_pred,qf.y1.mle)
        for(k in nns){
            vecname = paste0('fits/sim/fits_prcp_l',k,'.RDS')
            x.vec = readRDS(vecname)
            X2_pred = cbind(X2_pred,x.vec)
        }
        head(X2)
        head(X1)
        qf.y2.mle = rep(NA,n)
  
        for(i in 1:n){
            if(i%%1000==0)
                print(i)
            qf.y2.mle[i] <- predict(fit.y2.mle,   X = X2_pred[i,], type = "QF",tau=qout21[i])
        }   
    }
    saveRDS(qf.y2.mle,file = predname2)
    
    pdf(file = pdfname,width = 6,height = 6)
    par(mfrow=c(2,2))
    # 
    plot(y1,qf.y1.mle,col=y0+1, main = 'MLE-TS',pch=20,cex=0.2)
    abline(0,1)
    
    summary(qf.y1.mle[y0==0])
    summary(qf.y1.mle[y0==1])
    
    summary(y1[y0==1])
    summary(y1[y0==0])
    
    d0 <-density(y1[y0==0]) 
    d1 <-density(y1[y0==1]) 
    d2 <- density(qf.y1.mle[y0==0])
    plotmax.y = max(d0$y,d1$y,d2$y)
    plotmin.y = min(d0$y,d1$y,d2$y)
    plotmax.x = max(d0$x,d1$x,d2$x)
    plotmin.x = min(d0$x,d1$x,d2$x)
    plot(d0,col=1,ylim=range(c(plotmin.y,plotmax.y)),
         xlim=range(c(plotmin.x,plotmax.x)),ylab="Y1",main="Y1")
    lines(d1,col=2)
    lines(d2,col=3,lty=2)

    
    plot(y2,qf.y2.mle,col=y0+1, main = 'MLE-TS',pch=20,cex=0.2)
    abline(0,1)
    
    
    summary(qf.y2.mle[y0==0])
    summary(qf.y2.mle[y0==1])
    
    summary(y2[y0==1])
    summary(y2[y0==0])
    
    d0 <-density(y2[y0==0]) 
    d1 <-density(y2[y0==1]) 
    d2 <- density(qf.y2.mle[y0==0])
    plotmax.y = max(d0$y,d1$y,d2$y)
    plotmin.y = min(d0$y,d1$y,d2$y)
    plotmax.x = max(d0$x,d1$x,d2$x)
    plotmin.x = min(d0$x,d1$x,d2$x)
    plot(d0,col=1,ylim=range(c(plotmin.y,plotmax.y)),
         xlim=range(c(plotmin.x,plotmax.x)),ylab="Y2",main="Y2")
    lines(d1,col=2)
    lines(d2,col=2,lty=2)
    
    dev.off()
    par(mfrow=c(1,1))
}

