rm(list = ls())
library(GpGp)
library(SPQR)
library(lubridate)
region = 'SW'
load(file = 'data/simdata.RData')

coords = locs
head(coords)

set.seed(303)
vecchia.order = order_maxmin(coords,lonlat = T)
loc = 1


# Y.range <- range(Y)
# .Y <- (Y - Y.range[1])/diff(Y.range)

for(loc in 1:25){
    pdfname     <- paste0('plots/sim/fits_temp_l',loc,'.pdf')
    predname1   <- paste0( 'fits/sim/fits_temp_l',loc,'.RDS')
    predname2   <- paste0( 'fits/sim/fits_prcp_l',loc,'.RDS')
    
    y1 <- c(Temp0[,vecchia.order==loc],Temp1[,vecchia.order==loc])
    y2 <- c(Prec0[,vecchia.order==loc],Prec1[,vecchia.order==loc])
    y2 <- log(0.0001+y2)
    n0 <- n1 <- nrow(Temp0)
    n = n0 + n1
    y0 <- rep(1:0,each=n0)
    
    plot(y1,y2,col=y0+1)
    
    X1 = matrix(y0, ncol=1)
    
    if(loc>1){
        k.end = loc-1
        k.start = max(1,loc-5)
        for(k in k.start:k.end){
            x.vec = c(Temp0[,vecchia.order==loc-k],Temp1[vecchia.order==loc-k])
            X1 = cbind(X1,x.vec)
        }    
    }
    
    head(X1)
    # pdf(file = pdfname,width = 6,height = 6)
    # par(mfrow=c(2,2))
    control <- list(iter = 300, batch.size = 100, lr = 0.001, save.name = paste('SPQR.model.temp',loc,'pt',sep='.'))
    fit.y1.mle.ts <- SPQR(X = X1, Y = y1, method = "MLE", control = control, normalize = T, verbose = T,use.GPU=F,
                          n.hidden = c(30,20), activation = 'relu',n.knots = 20,seed = loc)
    # save.SPQR(fit.y1.mle.ts,name = modelname1)
    # plotGOF(fit.y1.mle.ts)
    cdf.y1.mle.ts = rep(NA,n)
    for(i in 1:n){
        cdf.y1.mle.ts[i] <- predict(fit.y1.mle.ts,   X = X1[i,], Y=y1[i], type = "CDF")    
        if(i%%1000==0)
            print(i)
    }
    
    qout11 <- cdf.y1.mle.ts
    adjust = which(qout11>0.99999)
    qout11[adjust] = 0.99999
    
    
    ###################################
    ###################################
    X2 = cbind(X1,y1)
    nx1 = ncol(X1)+1
    if(loc>1){
        k.end = loc-1
        k.start = max(1,loc-5)
        for(k in k.start:k.end){
            x.vec = c(Prec0[vecchia.order==loc-k],Prec1[vecchia.order==loc-k])
            x.vec = log(0.0001+x.vec)
            X2 = cbind(X2,x.vec)
        }    
    }
    nx2 = ncol(X2)
    head(X2)
    head(y2)
    control <- list(iter = 300, batch.size = 100, lr = 0.001, save.name = paste('SPQR.model.prcp',loc,'pt',sep='.'))
    fit.y2.mle.ts <- SPQR(X = X2, Y = y2, method = "MLE", control = control, normalize = T, verbose = T,use.GPU=F,
                          n.hidden = c(30,20), activation = 'relu',n.knots = 20,seed = loc)
    # save.SPQR(fit.y2.mle.ts,name = modelname2)
    # plotGOF(fit.y2.mle.ts)
    cdf.y2.mle.ts = rep(NA,n)
    for(i in 1:n){
        cdf.y2.mle.ts[i] <- predict(fit.y2.mle.ts,   X = X2[i,], Y=y2[i], type = "CDF")   
        if(i%%1000==0)
            print(i)
    }
    
    qout21 <- cdf.y2.mle.ts
    adjust = which(qout21>0.99999)
    qout21[adjust] = 0.99999
    
    
############### Predictions temp     
    
    if(loc==1){
        qf.y1.mle.ts = rep(NA,n)
        for(i in 1:n){
            x_pred = 1
            qf.y1.mle.ts[i] <- predict(fit.y1.mle.ts,   X = x_pred, type = "QF",tau=qout11[i])
        }    
    }
    
    if(loc>1){
        X1_pred = X1[,1]
        k.end = loc-1
        k.start = max(1,loc-5)
        for(k in k.start:k.end){
            vecname = paste0('fits/sim/fits_temp_m',mnth,'_l',loc-k,'.RDS')
            x.vec = readRDS(vecname)
            X1_pred = cbind(X1_pred,x.vec)
        }
        qf.y1.mle.ts = rep(NA,n)
        for(i in 1:n){
            if(i%%1000==0)
                print(i)
            if(i<=n0+1)
                x_pred = c(1,X1[i,2],X1[i,-c(1:2)])
            if(i>n0+1)
                x_pred = c(1,qf.y1.mle.ts[i-1],X1_pred[i,-c(1:2)])
            qf.y1.mle.ts[i] <- predict(fit.y1.mle.ts,   X = x_pred, type = "QF",tau=qout11[i])
        }   
    }
    
    saveRDS(qf.y1.mle.ts,file = predname1)
############### Predictions prcp        
    if(loc==1){
        qf.y2.mle.ts = rep(NA,n)
        for(i in 1:n){
            if(i%%1000==0)
                print(i)
            if(i<=n0+1)
                x_pred = c(1,y1[i])
            if(i>n0+1)
                x_pred = c(1,qf.y1.mle.ts[i])
            qf.y2.mle.ts[i] <- predict(fit.y2.mle.ts,   X = x_pred, type = "QF",tau=qout21[i])
        } 
    }
    
    if(loc>1){
        X2_pred = cbind(X1_pred,qf.y1.mle.ts,X2[,nx1+1])
        k.end = loc-1
        k.start = max(1,loc-5)
        for(k in k.start:k.end){
            vecname = paste0('fits/',region,'/fits_prcp_m',mnth,'_l',loc-k,'.RDS')
            x.vec = readRDS(vecname)
            X2_pred = cbind(X2_pred,x.vec)
        }
        head(X2)
        head(X1)
        qf.y2.mle.ts = rep(NA,n)
  
        for(i in 1:n){
            if(i%%1000==0)
                print(i)
            if(i<=n0+1)
                x_pred = c(1,X1[i,-1],y1[i],X2[i,-c(1:nx1)])
            if(i>n0+1)
                x_pred = c(1,qf.y1.mle.ts[i-1],X1_pred[i,-c(1:2)],qf.y1.mle.ts[i],qf.y2.mle.ts[i-1],X2_pred[i,(nx1+2):nx2])
            qf.y2.mle.ts[i] <- predict(fit.y2.mle.ts,   X = x_pred, type = "QF",tau=qout21[i])
        }   
    }
    saveRDS(qf.y2.mle.ts,file = predname2)
    
    pdf(file = pdfname,width = 6,height = 6)
    par(mfrow=c(2,2))
    
    plot(y1,qf.y1.mle.ts,col=y0+1, main = 'MLE-TS',pch=20,cex=0.2)
    abline(0,1)
    
    summary(qf.y1.mle.ts[y0==0])
    summary(qf.y1.mle.ts[y0==1])
    cor(qf.y1.mle.ts[y0==0][-1],qf.y1.mle.ts[y0==0][-n0])
    cor(qf.y1.mle.ts[y0==1][-1],qf.y1.mle.ts[y0==1][-n0])
    
    summary(y1[y0==1])
    summary(y1[y0==0])
    
    d0 <-density(y1[y0==0]) 
    d1 <-density(y1[y0==1]) 
    d2 <- density(qf.y1.mle.ts[y0==0])
    plotmax.y = max(d0$y,d1$y,d2$y)
    plotmin.y = min(d0$y,d1$y,d2$y)
    plotmax.x = max(d0$x,d1$x,d2$x)
    plotmin.x = min(d0$x,d1$x,d2$x)
    plot(d0,col=1,ylim=range(c(plotmin.y,plotmax.y)),
         xlim=range(c(plotmin.x,plotmax.x)),ylab="Y1",main="Y1")
    lines(d1,col=2)
    lines(d2,col=3,lty=2)
    
    # legend("topright",c("Model","Observations",'MLE-TS'),
    #        col=1:6,lwd=2,bty="n",lty=c(1,1,2,3))
    
    
    
    
    plot(y2,qf.y2.mle.ts,col=y0+1, main = 'MLE-TS',pch=20,cex=0.2)
    abline(0,1)
    
    
    summary(qf.y2.mle.ts[y0==0])
    summary(qf.y2.mle.ts[y0==1])
    cor(qf.y2.mle.ts[y0==0][-1],qf.y2.mle.ts[y0==0][-n0])
    cor(qf.y2.mle.ts[y0==1][-1],qf.y2.mle.ts[y0==1][-n0])
    
    summary(y2[y0==1])
    summary(y2[y0==0])
    
    d0 <-density(y2[y0==0]) 
    d1 <-density(y2[y0==1]) 
    d2 <- density(qf.y2.mle.ts[y0==0])
    plotmax.y = max(d0$y,d1$y,d2$y)
    plotmin.y = min(d0$y,d1$y,d2$y)
    plotmax.x = max(d0$x,d1$x,d2$x)
    plotmin.x = min(d0$x,d1$x,d2$x)
    plot(d0,col=1,ylim=range(c(plotmin.y,plotmax.y)),
         xlim=range(c(plotmin.x,plotmax.x)),ylab="Y2",main="Y2")
    lines(d1,col=2)
    lines(d2,col=2,lty=2)
    
    # legend("topright",c("Model","Observations",'MLE-TS'),
    #        col=1:6,lwd=2,bty="n",lty=c(1,1,2,3))
    
    # plot(qf.y1.mle.ts,qf.y2.mle.ts,col = y0+1)
    
    cor(y1[y0==1],y2[y0==1])
    cor(qf.y1.mle.ts[y0==1],qf.y2.mle.ts[y0==1])
    dev.off()
    par(mfrow=c(1,1))
}
