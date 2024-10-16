rm(list = ls())
library(GpGp)
library(SPQR)
library(lubridate)
region = 'SE'
gcm.long = read.csv(paste0('data/',region,'_gcm_data.csv'))
obs.long = read.csv(paste0('data/',region,'_obs_data.csv'))

gcm.months = month(gcm.long[,1])
gcm.years = year(gcm.long[,1])
numdays = c(31,28,31,30,31,30,31,31,30,31,30,31)

grid.no = as.factor(gcm.long$lat*gcm.long$lon)
str(grid.no)
coords = aggregate(gcm.long[,2:3],by = list(grid.no), FUN = mean)
coords = coords[,-1]
head(coords)
table(grid.no)
set.seed(303)
vecchia.order = order_maxmin(coords,lonlat = T)
loc = 1
mnth = 1

# Y.range <- range(Y)
# .Y <- (Y - Y.range[1])/diff(Y.range)

for(mnth in 1:12)
    for(loc in 1:25){
        pdfname = paste0('plots/',region,'_validation/fits_temp_m',mnth,'_l',loc,'.pdf')
        modelname1 = paste0('fits/',region,'_validation/fits_temp_m',mnth,'_l',loc)
        predname1 = paste0('fits/',region,'_validation/fits_temp_m',mnth,'_l',loc,'.RDS')
        modelname2 = paste0('fits/',region,'_validation/fits_prcp_m',mnth,'_l',loc)
        predname2 = paste0('fits/',region,'_validation/fits_prcp_m',mnth,'_l',loc,'.RDS')
        
        y2_train <- c(obs.long$tmax[vecchia.order==loc & gcm.months==mnth & gcm.years <= 2000],
                      gcm.long$tmax[vecchia.order==loc & gcm.months==mnth & gcm.years <= 2000])
        y1_train <- c(  obs.long$pr[vecchia.order==loc & gcm.months==mnth & gcm.years <= 2000],
                        gcm.long$pr[vecchia.order==loc & gcm.months==mnth & gcm.years <= 2000])
        y2_test  <- c(obs.long$tmax[vecchia.order==loc & gcm.months==mnth & gcm.years > 2000],
                      gcm.long$tmax[vecchia.order==loc & gcm.months==mnth & gcm.years > 2000])
        y1_test  <- c(  obs.long$pr[vecchia.order==loc & gcm.months==mnth & gcm.years > 2000],
                        gcm.long$pr[vecchia.order==loc & gcm.months==mnth & gcm.years > 2000])
        y1_train <- log(0.0001+y1_train)
        y1_test <- log(0.0001+y1_test)
        
        n0_train <- n1_train <- length(y1_train)/2
        n0_test  <- n1_test  <- length(y1_test)/2
        n_train = n0_train + n1_train
        n_test = n0_test + n1_test
        
        # temperature covariates for training data
        y0_train <- rep(1:0,each=n0_train)
        y11_train = y1_train[y0_train==1]
        x11_train = y11_train[c(n1_train,1:(n1_train-1))]
        cor(cbind(y11_train,x11_train))
        y10_train = y1_train[y0_train==0]
        x10_train = y10_train[c(n0_train,1:(n0_train-1))]
        cor(cbind(y10_train,x10_train))
        y1_train = c(y11_train,y10_train)
        x1_train = c(x11_train,x10_train)
        X1_train = cbind(y0_train,x1_train)
        
        if(loc>1){
            k.end = loc-1
            k.start = max(1,loc-5)
            for(k in k.start:k.end){
                x.vec = c(obs.long$pr[vecchia.order==loc-k & gcm.months==mnth & gcm.years <= 2000],
                          gcm.long$pr[vecchia.order==loc-k & gcm.months==mnth & gcm.years <= 2000])
                x.vec = log(0.0001+x.vec)
                X1_train = cbind(X1_train,x.vec)
            }    
        }
        
        # temperature covariates for testing data
        y0_test  <- rep(1:0,each=n0_test)
        y0_test <- rep(1:0,each=n0_test)
        y11_test = y1_test[y0_test==1]
        x11_test = y11_test[c(n1_test,1:(n1_test-1))]
        cor(cbind(y11_test,x11_test))
        y10_test = y1_test[y0_test==0]
        x10_test = y10_test[c(n0_test,1:(n0_test-1))]
        cor(cbind(y10_test,x10_test))
        y1_test = c(y11_test,y10_test)
        x1_test = c(x11_test,x10_test)
        X1_test = cbind(y0_test,x1_test)
        
        if(loc>1){
            k.end = loc-1
            k.start = max(1,loc-5)
            for(k in k.start:k.end){
                x.vec = c(obs.long$pr[vecchia.order==loc-k & gcm.months==mnth & gcm.years > 2000],
                          gcm.long$pr[vecchia.order==loc-k & gcm.months==mnth & gcm.years > 2000])
                x.vec = log(0.0001+x.vec)
                X1_test = cbind(X1_test,x.vec)
            }    
        }
        
        # normalize variables
        y1_range <- range(y1_train,y1_test)
        y1_train_scaled <- (y1_train - y1_range[1])/diff(y1_range)
        y1_test_scaled <- (y1_test - y1_range[1])/diff(y1_range)
        
        x1_range <- apply(rbind(X1_train,X1_test),2,range)
        X1_train_scaled <- X1_train
        X1_test_scaled <- X1_test
        for(i in 1:ncol(X1_train)){
            X1_train_scaled[,i] <- (X1_train[,i] - x1_range[1,i])/diff(x1_range[,i])
            X1_test_scaled[,i] <- (X1_test[,i] - x1_range[1,i])/diff(x1_range[,i])
        }
        
        # pdf(file = pdfname,width = 6,height = 6)
        # par(mfrow=c(2,2))
        control <- list(iter = 300, batch.size = 100, lr = 0.001, save.name = paste('SPQR.model.prcp',region,mnth,loc,'pt',sep='.'))
        fit.y1.mle.ts <- SPQR(X = X1_train_scaled, Y = y1_train_scaled, method = "MLE", control = control, normalize = F, verbose = T,use.GPU=F,
                              n.hidden = c(30,20), activation = 'relu',n.knots = 20, seed = mnth*loc)
        # save.SPQR(fit.y1.mle.ts,name = modelname1)
        # plotGOF(fit.y1.mle.ts)
        cdf.y1.mle.ts = rep(NA,n_test)
        for(i in 1:n_test){
            cdf.y1.mle.ts[i] <- predict(fit.y1.mle.ts,   X = X1_test_scaled[i,], Y=y1_test_scaled[i], type = "CDF")    
            if(i%%100==0)
                print(i)
        }
        
        qout11 <- cdf.y1.mle.ts
        # adjust = which(qout11>0.99999)
        # qout11[adjust] = 0.99999
        
        
        ###################################
        ###################################
        # Precip covariates for training data
        y21_train = y2_train[y0_train==1]
        x21_train = y21_train[c(n1_train,1:(n1_train-1))]
        cor(cbind(y21_train,x21_train))
        y20_train = y2_train[y0_train==0]
        x20_train = y20_train[c(n0_train,1:(n0_train-1))]
        cor(cbind(y20_train,x20_train))
        y2_train = c(y21_train,y20_train)
        x2_train = c(x21_train,x20_train)
        # X2 = cbind(x2,X1[,2],y1,X1[,-2])
        X2_train = cbind(X1_train,y1_train,x2_train)
        nx1 = ncol(X1_train)+1
        if(loc>1){
            k.end = loc-1
            k.start = max(1,loc-5)
            for(k in k.start:k.end){
                x.vec = c(obs.long$tmax[vecchia.order==loc-k & gcm.months==mnth & gcm.years <= 2000],
                          gcm.long$tmax[vecchia.order==loc-k & gcm.months==mnth & gcm.years <= 2000])
                # x.vec = log(0.0001+x.vec)
                X2_train = cbind(X2_train,x.vec)
            }    
        }
        nx2 = ncol(X2_train)
        
        # Precip covariates for testing data
        y21_test = y2_test[y0_test==1]
        x21_test = y21_test[c(n1_test,1:(n1_test-1))]
        cor(cbind(y21_test,x21_test))
        y20_test = y2_test[y0_test==0]
        x20_test = y20_test[c(n0_test,1:(n0_test-1))]
        cor(cbind(y20_test,x20_test))
        y2_test = c(y21_test,y20_test)
        x2_test = c(x21_test,x20_test)
        # X2 = cbind(x2,X1[,2],y1,X1[,-2])
        X2_test = cbind(X1_test,y1_test,x2_test)
        nx1 = ncol(X1_test)+1
        if(loc>1){
            k.end = loc-1
            k.start = max(1,loc-5)
            for(k in k.start:k.end){
                x.vec = c(obs.long$tmax[vecchia.order==loc-k & gcm.months==mnth & gcm.years > 2000],
                          gcm.long$tmax[vecchia.order==loc-k & gcm.months==mnth & gcm.years > 2000])
                # x.vec = log(0.0001+x.vec)
                X2_test = cbind(X2_test,x.vec)
            }    
        }
        nx2 = ncol(X2_test)
        
        # normalize variables
        y2_range <- range(y2_train,y2_test)
        y2_train_scaled <- (y2_train - y2_range[1])/diff(y2_range)
        y2_test_scaled <- (y2_test - y2_range[1])/diff(y2_range)
        
        x2_range <- apply(rbind(X2_train,X2_test),2,range)
        X2_train_scaled <- X2_train
        X2_test_scaled <- X2_test
        for(i in 1:ncol(X2_train)){
            X2_train_scaled[,i] <- (X2_train[,i] - x2_range[1,i])/diff(x2_range[,i])
            X2_test_scaled[,i] <- (X2_test[,i] - x2_range[1,i])/diff(x2_range[,i])
        }
      
        control <- list(iter = 300, batch.size = 100, lr = 0.001, save.name = paste('SPQR.model.prcp',region,mnth,loc,'pt',sep='.'))
        fit.y2.mle.ts <- SPQR(X = X2_train_scaled, Y = y2_train_scaled, method = "MLE", control = control, normalize = F, verbose = T,use.GPU=F,
                              n.hidden = c(30,20), activation = 'relu',n.knots = 20, seed = mnth*loc)
        # save.SPQR(fit.y2.mle.ts,name = modelname2)
        # plotGOF(fit.y2.mle.ts)
        cdf.y2.mle.ts = rep(NA,n_test)
        for(i in 1:n_test){
            cdf.y2.mle.ts[i] <- predict(fit.y2.mle.ts,   X = X2_test_scaled[i,], Y=y2_test_scaled[i], type = "CDF")   
            if(i%%1000==0)
                print(i)
        }
        
        qout21 <- cdf.y2.mle.ts
        adjust = which(qout21>0.99999)
        qout21[adjust] = 0.99999
        
        
 ############### Predictions temp     
        
        if(loc==1){
            qf.y1.mle.ts = rep(NA,n_test)
            for(i in 1:n_test){
                if(i%%1000==0)
                    print(i)
                if(i<=n0_test+1) #these are the ones which are already observed (+ the first gcm one)
                    x_pred = c(1,X1_test_scaled[i,-1])
                if(i>n0_test+1) #these are the GCM ones which therefore need the time lagged GCM predictions
                    x_pred = c(1,qf.y1.mle.ts[i-1])
                qf.y1.mle.ts[i] <- predict(fit.y1.mle.ts,   X = x_pred, type = "QF",tau=qout11[i])
            }    
        }
        
        if(loc>1){
            X1_pred = X1_test[,1:2]
            k.end = loc-1
            k.start = max(1,loc-5)
            for(k in k.start:k.end){
                vecname = paste0('fits/',region,'_validation/fits_prcp_m',mnth,'_l',loc-k,'.RDS')
                x.vec = readRDS(vecname)
                X1_pred = cbind(X1_pred,x.vec)
            }
            
            # need to scale the x.vec
            X1_pred_scaled <- X1_pred
            for(i in 1:ncol(X1_train)){
                X1_pred_scaled[,i] <- (X1_pred[,i] - x1_range[1,i])/diff(x1_range[,i])
            }
            
            
            qf.y1.mle.ts = rep(NA,n_test)
            for(i in 1:n_test){
                if(i%%100==0)
                    print(i)
                if(i<=n0_test+1)
                    x_pred = c(1,X1_test_scaled[i,2],X1_test_scaled[i,-c(1:2)])
                if(i>n0_test+1)
                    x_pred = c(1,qf.y1.mle.ts[i-1],X1_pred_scaled[i,-c(1:2)])
                qf.y1.mle.ts[i] <- predict(fit.y1.mle.ts,   X = x_pred, type = "QF",tau=qout11[i])
            }   
        }
        
        y1_pred <- qf.y1.mle.ts*diff(y1_range) + y1_range[1]
        saveRDS(y1_pred,file = predname1)
 ############### Predictions prcp        
        if(loc==1){
            qf.y2.mle.ts = rep(NA,n_test)
            for(i in 1:n_test){
                if(i%%1000==0)
                    print(i)
                if(i<=n0_test+1)
                    x_pred = c(1,X1_test_scaled[i,-1],y1_test_scaled[i],X2_test_scaled[i,-c(1:nx1)])
                if(i>n0_test+1)
                    x_pred = c(1,qf.y1.mle.ts[i-1],qf.y1.mle.ts[i],qf.y2.mle.ts[i-1])
                qf.y2.mle.ts[i] <- predict(fit.y2.mle.ts,   X = x_pred, type = "QF",tau=qout21[i])
            } 
        }
        
        if(loc>1){
            X2_pred = cbind(X1_pred_scaled,qf.y1.mle.ts,X2_test_scaled[,nx1+1])
            nx3 <- ncol(X2_pred)+1
            k.end = loc-1
            k.start = max(1,loc-5)
            for(k in k.start:k.end){
                vecname = paste0('fits/',region,'_validation/fits_temp_m',mnth,'_l',loc-k,'.RDS')
                x.vec = readRDS(vecname)
                X2_pred = cbind(X2_pred,x.vec)
            }
            
            # need to scale the x.vec
            X2_pred_scaled <- X2_pred
            for(i in nx3:ncol(X2_train)){
                X2_pred_scaled[,i] <- (X2_pred[,i] - x2_range[1,i])/diff(x2_range[,i])
            }
            
            qf.y2.mle.ts = rep(NA,n_test)
      
            for(i in 1:n_test){
                if(i%%100==0)
                    print(i)
                if(i<=n0_test+1)
                    x_pred = c(1,X1_test_scaled[i,-1],y1_test_scaled[i],X2_test_scaled[i,-c(1:nx1)])
                if(i>n0_test+1)
                    x_pred = c(1,qf.y1.mle.ts[i-1],X1_pred_scaled[i,-c(1:2)],qf.y1.mle.ts[i],qf.y2.mle.ts[i-1],X2_pred_scaled[i,(nx1+2):nx2])
                qf.y2.mle.ts[i] <- predict(fit.y2.mle.ts,   X = x_pred, type = "QF",tau=qout21[i])
            }   
        }
        y2_pred <- qf.y2.mle.ts*diff(y2_range) + y2_range[1]
        saveRDS(y2_pred,file = predname2)
        
        pdf(file = pdfname,width = 6,height = 6)
        par(mfrow=c(2,2))
        
        plot(y1_test,y1_pred,col=y0_test+1, main = 'MLE-TS',pch=20,cex=0.2)
        abline(0,1)
        
        summary(y1_pred[y0_test==0])
        summary(y1_pred[y0_test==1])
        cor(y1_pred[y0_test==0][-1],y1_pred[y0_test==0][-n0_test])
        cor(y1_pred[y0_test==1][-1],y1_pred[y0_test==1][-n0_test])
        
        summary(y1_test[y0_test==1])
        summary(y1_test[y0_test==0])
        
        d0 <-density(y1_test[y0_test==0]) 
        d1 <-density(y1_test[y0_test==1]) 
        d2 <- density(y1_pred[y0_test==0])
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
        
        
        
        # y2_pred <- exp(y2_pred) - 0.0001
        # y2_test <- exp(y2_test) - 0.0001
        plot(y2_test,y2_pred,col=y0_test+1, main = 'MLE-TS',pch=20,cex=0.2)
        abline(0,1)
        
        
        summary(y2_pred[y0_test==0])
        summary(y2_pred[y0_test==1])
        cor(y2_pred[y0_test==0][-1],y2_pred[y0_test==0][-n0_test])
        cor(y2_pred[y0_test==1][-1],y2_pred[y0_test==1][-n0_test])
        
        summary(y2_test[y0_test==1])
        summary(y2_test[y0_test==0])
        
        d0 <-density(y2_test[y0_test==0])
        d1 <-density(y2_test[y0_test==1])
        d2 <- density(y2_pred[y0_test==0])
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
        
        cor(y1_test[y0_test==1],y2_test[y0_test==1])
        cor(y1_pred[y0_test==1],y2_pred[y0_test==1])
        par(mfrow=c(1,1))
        dev.off()
    }
