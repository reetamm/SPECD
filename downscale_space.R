rm(list = ls())
library(GpGp)
library(SPQR)
library(lubridate)
use.gpu = torch::cuda_is_available()
gcm.long = read.csv('selectedgcmdata.csv')
obs.long = read.csv('selectedobsdata.csv')

gcm.months = month(gcm.long[,1])

grid.no = as.factor(gcm.long$lat*gcm.long$lon)
str(grid.no)
coords = aggregate(gcm.long[,2:3],by = list(grid.no), FUN = mean)
coords = coords[,-1]
head(coords)
table(grid.no)
set.seed(303)
vecchia.order = order_maxmin(coords,lonlat = T)
loc = 3
mnth = 9

for(mnth in 1:12)
    for(loc in 1:25){
        pdfname = paste0('plots/fits_temp_m',mnth,'_l',loc,'.pdf')
        modelname1 = paste0('fits/fits_temp_m',mnth,'_l',loc)
        predname1 = paste0('fits/fits_temp_m',mnth,'_l',loc,'.RDS')
        modelname2 = paste0('fits/fits_prcp_m',mnth,'_l',loc)
        predname2 = paste0('fits/fits_prcp_m',mnth,'_l',loc,'.RDS')
        
        y1 <- c(obs.long$tmax[vecchia.order==loc & gcm.months==mnth],gcm.long$tmax[vecchia.order==loc & gcm.months==mnth])
        y2 <- c(obs.long$pr[vecchia.order==loc & gcm.months==mnth],gcm.long$pr[vecchia.order==loc & gcm.months==mnth])
        y2 <- log(1+y2)
        n0 = length(gcm.long$pr[vecchia.order==loc & gcm.months==mnth]); n1 = length(obs.long$pr[vecchia.order==loc & gcm.months==mnth])
        n = n0 + n1
        y0 <- rep(1:0,each=n0)
        
        plot(y1,y2,col=y0+1)
        
        y11 = y1[y0==1]
        y10 = y1[y0==0]
        y1 = c(y11,y10)
        X1 = y0
        if(loc==1)
            X1 = matrix(X1,nrow = n, ncol = 1)
        if(loc>1){
            k.end = loc-1
            k.start = max(1,loc-5)
            for(k in k.start:k.end){
                x.vec = c(obs.long$tmax[vecchia.order==loc-k & gcm.months==mnth],gcm.long$tmax[vecchia.order==loc-k & gcm.months==mnth])
                X1 = cbind(X1,x.vec)
            }    
        }
        
        head(X1)
        # pdf(file = pdfname,width = 6,height = 6)
        # par(mfrow=c(2,2))
        control <- list(iter = 300, batch.size = 100, lr = 0.001)
        fit.y1.mle.ts <- SPQR(X = X1, Y = y1, method = "MLE", control = control, normalize = T, verbose = T,use.GPU=use.gpu,
                              n.hidden = c(30,20), activation = 'relu',n.knots = 15)
        save.SPQR(fit.y1.mle.ts,name = modelname1)
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
        y21 = y2[y0==1]
        y20 = y2[y0==0]
        y2 = c(y21,y20)
        X2 = cbind(X1,y1)

        if(loc>1){
            k.end = loc-1
            k.start = max(1,loc-5)
            for(k in k.start:k.end){
                x.vec = c(obs.long$pr[vecchia.order==loc-k & gcm.months==mnth],gcm.long$pr[vecchia.order==loc-k & gcm.months==mnth])
                x.vec = log(1+x.vec)
                X2 = cbind(X2,x.vec)
            }    
        }
        
        head(X2)
        head(y2)
        control <- list(iter = 300, batch.size = 100, lr = 0.001)
        fit.y2.mle.ts <- SPQR(X = X2, Y = y2, method = "MLE", control = control, normalize = T, verbose = T,use.GPU=use.gpu,
                              n.hidden = c(30,20), activation = 'relu',n.knots = 20)
        save.SPQR(fit.y2.mle.ts,name = modelname2)
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
        
        if(loc==1)
            qf.y1.mle.ts <- predict(fit.y1.mle.ts,   X = 1, type = "QF",tau=qout11)
        
        if(loc>1){
            X1_pred = matrix(X1[,1],nrow = n,ncol = 1)
            k.end = loc-1
            k.start = max(1,loc-5)
            for(k in k.start:k.end){
                vecname = paste0('fits/fits_temp_m',mnth,'_l',loc-k,'.RDS')
                x.vec = readRDS(vecname)
                X1_pred = cbind(X1_pred,x.vec)
            }
            X1_pred[y0==1,-1] = X1[y0==1,-1]
            X1_pred[,1] = 1
            rownames(X1_pred) = NULL
            qf.y1.mle.ts = rep(NA,n)
        
            for(i in 1:n){
                if(i%%1000==0)
                    print(i)
                qf.y1.mle.ts[i] <- predict(fit.y1.mle.ts,   X = X1_pred[i,], type = "QF",tau=qout11[i])
            }   
        }
        
        saveRDS(qf.y1.mle.ts,file = predname1)
 ############### Predictions prcp        
        if(loc==1){
            qf.y2.mle.ts = rep(NA,n)
            X2_pred = X2
            X2_pred[y0==0,2] = qf.y1.mle.ts[y0==0]
            X2_pred[,1] = 1
            for(i in 1:n){
                if(i%%1000==0)
                    print(i)
                qf.y2.mle.ts[i] <- predict(fit.y2.mle.ts,   X = X2_pred[i,], type = "QF",tau=qout21[i])
            } 
        }
        
        if(loc>1){
            X2_pred = cbind(y0,X1_pred[,-1],qf.y1.mle.ts)
            k.end = loc-1
            k.start = max(1,loc-5)
            for(k in k.start:k.end){
                vecname = paste0('fits/fits_prcp_m',mnth,'_l',loc-k,'.RDS')
                x.vec = readRDS(vecname)
                X2_pred = cbind(X2_pred,x.vec)
            }
            head(X2[y0==0,])
            head(X2_pred[y0==0,])
            X2_pred[y0==1,-1] = X2[y0==1,-1]
            X2_pred[,1] = 1
            
            qf.y2.mle.ts = rep(NA,n)
            for(i in 1:n){
                if(i%%1000==0)
                    print(i)
                qf.y2.mle.ts[i] <- predict(fit.y2.mle.ts,   X = X2_pred[i,], type = "QF",tau=qout21[i])
            }   
        }
        saveRDS(qf.y2.mle.ts,file = predname2)
        
        pdf(file = pdfname,width = 6,height = 6)
        par(mfrow=c(2,2))
        
        plot(y1,qf.y1.mle.ts,col=y0+1, main = 'MLE-TS',pch=20,cex=0.2)
        abline(0,1)
        
        summary(qf.y1.mle.ts[y0==0])
        summary(qf.y1.mle.ts[y0==1])
        
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
