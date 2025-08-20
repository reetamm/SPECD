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
NNarray <- find_ordered_nn(coords[vecchia.order,],lonlat = T,m=5)
loc = 3
mnth = 1
mnths = 1:3
loc.vector <- 1:25

for(mnth in mnths)
    for(loc in 1:25){
        
        cur.loc <- vecchia.order[loc]
        nns <- NNarray[loc,] # select the correct row
        nns <- nns[complete.cases(nns)] # drop the NAs
        nns <- nns[-1] # drop the response
        nns <- vecchia.order[nns]
        
        pdfname = paste0('plots/',region,'_validation/fits_m',mnth,'_l',loc,'.pdf')
        predname1 = paste0('fits/',region,'_validation/fits_temp_m',mnth,'_l',loc,'.RDS')
        predname2 = paste0('fits/',region,'_validation/fits_prcp_m',mnth,'_l',loc,'.RDS')
        
        y1_train <- c(obs.long$tmax[loc.vector==cur.loc & gcm.months==mnth & gcm.years <= 2000],
                      gcm.long$tmax[loc.vector==cur.loc & gcm.months==mnth & gcm.years <= 2000])
        y2_train <- c(  obs.long$pr[loc.vector==cur.loc & gcm.months==mnth & gcm.years <= 2000],
                        gcm.long$pr[loc.vector==cur.loc & gcm.months==mnth & gcm.years <= 2000])
        y1_test  <- c(obs.long$tmax[loc.vector==cur.loc & gcm.months==mnth & gcm.years > 2000],
                      gcm.long$tmax[loc.vector==cur.loc & gcm.months==mnth & gcm.years > 2000])
        y2_test  <- c(  obs.long$pr[loc.vector==cur.loc & gcm.months==mnth & gcm.years > 2000],
                        gcm.long$pr[loc.vector==cur.loc & gcm.months==mnth & gcm.years > 2000])
        y2_train <- log(0.0001+y2_train)
        y2_test <- log(0.0001+y2_test)
        
        n0_train <- n1_train <- length(y1_train)/2
        n0_test  <- n1_test  <- length(y1_test)/2
        n_train = n0_train + n1_train
        n_test = n0_test + n1_test
        
        # temperature covariates for training data
        y0_train <- rep(1:0,each=n0_train)
        X1_train = matrix(y0_train,ncol=1)
        
        if(loc>1){
            for(k in nns){
                x.vec = c(obs.long$tmax[loc.vector==k & gcm.months==mnth & gcm.years <= 2000],
                          gcm.long$tmax[loc.vector==k & gcm.months==mnth & gcm.years <= 2000])
                X1_train = cbind(X1_train,x.vec)
            }    
        }
        
        # temperature covariates for testing data
        y0_test  <- rep(1:0,each=n0_test)
        X1_test = matrix(y0_test,ncol=1)
        
        if(loc>1){
            for(k in nns){
                x.vec = c(obs.long$tmax[loc.vector==k & gcm.months==mnth & gcm.years > 2000],
                          gcm.long$tmax[loc.vector==k & gcm.months==mnth & gcm.years > 2000])
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
        
        print("Fitting temp")
        control <- list(iter = 300, batch.size = 100, lr = 0.001, save.name = paste('SPQR.model.temp',region,mnth,loc,'v.pt',sep='.'))
        fit.y1.mle <- SPQR(X = X1_train_scaled, Y = y1_train_scaled, method = "MLE", control = control, normalize = F, verbose = T,use.GPU=F,
                              n.hidden = c(30,20), activation = 'relu',n.knots = 20, seed = mnth*loc)

        print("CDF transform of temp")
        cdf.y1.mle = rep(NA,n_test)
        for(i in 1:n_test){
            cdf.y1.mle[i] <- predict(fit.y1.mle,   X = X1_test_scaled[i,], Y=y1_test_scaled[i], type = "CDF")    
            if(i%%100==0)
                print(i)
        }
        
        qout11 <- cdf.y1.mle
        adjust = which(qout11>0.99999)
        qout11[adjust] = 0.99999
        
        
        ###################################
        ###################################
        # Precip covariates for training data
        X2_train = cbind(X1_train,y1_train)
        nx1 = ncol(X1_train)+1
        if(loc>1){
            for(k in nns){
                x.vec = c(obs.long$pr[loc.vector==k & gcm.months==mnth & gcm.years <= 2000],
                          gcm.long$pr[loc.vector==k & gcm.months==mnth & gcm.years <= 2000])
                x.vec = log(0.0001+x.vec)
                X2_train = cbind(X2_train,x.vec)
            }    
        }
        nx2 = ncol(X2_train)
        
        # Precip covariates for testing data
        X2_test = cbind(X1_test,y1_test)
        nx1 = ncol(X1_test)+1
        if(loc>1){
            for(k in nns){
                x.vec = c(obs.long$pr[loc.vector==k & gcm.months==mnth & gcm.years > 2000],
                          gcm.long$pr[loc.vector==k & gcm.months==mnth & gcm.years > 2000])
                x.vec = log(0.0001+x.vec)
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
      
        print("Fitting prcp")
        control <- list(iter = 300, batch.size = 100, lr = 0.001, save.name = paste('SPQR.model.prcp',region,mnth,loc,'v.pt',sep='.'))
        fit.y2.mle <- SPQR(X = X2_train_scaled, Y = y2_train_scaled, method = "MLE", control = control, normalize = F, verbose = T,use.GPU=F,
                              n.hidden = c(30,20), activation = 'relu',n.knots = 20, seed = mnth*loc)

        print("CDF transform of prcp")
        cdf.y2.mle = rep(NA,n_test)
        for(i in 1:n_test){
            cdf.y2.mle[i] <- predict(fit.y2.mle,   X = X2_test_scaled[i,], Y=y2_test_scaled[i], type = "CDF")   
            if(i%%100==0)
                print(i)
        }
        
        qout21 <- cdf.y2.mle
        adjust = which(qout21>0.99999)
        qout21[adjust] = 0.99999
        
        ###### vecchia locs for predictions
        nns <- NNarray[loc,] # select the correct row
        nns <- nns[complete.cases(nns)] # drop the NAs
        nns <- nns[-1] # drop the response
        
 ############### Predictions temp     
        print("Predictions temp")
        if(loc==1){
            qf.y1.mle = rep(NA,n_test)
            for(i in 1:n_test){
                if(i%%100==0)
                    print(i)
               x_pred = 1
                qf.y1.mle[i] <- predict(fit.y1.mle,   X = x_pred, type = "QF",tau=qout11[i])
            }    
        }
        
        if(loc>1){
            X1_pred = X1_test[,1]
            for(k in nns){
                vecname = paste0('fits/',region,'_validation/fits_temp_m',mnth,'_l',k,'.RDS')
                x.vec = readRDS(vecname)
                X1_pred = cbind(X1_pred,x.vec)
            }
            
            # need to scale the x.vec
            X1_pred_scaled <- X1_pred
            for(i in 1:ncol(X1_train)){
                X1_pred_scaled[,i] <- (X1_pred[,i] - x1_range[1,i])/diff(x1_range[,i])
            }
            
            X1_pred_scaled[,1] <- 1
            qf.y1.mle = rep(NA,n_test)
            for(i in 1:n_test){
                if(i%%100==0)
                    print(i)
                qf.y1.mle[i] <- predict(fit.y1.mle,   X = X1_pred_scaled[i,], type = "QF",tau=qout11[i])
            }   
        }
        
        y1_pred <- qf.y1.mle*diff(y1_range) + y1_range[1]
        saveRDS(y1_pred,file = predname1)
 ############### Predictions prcp        
        print("Predictions prcp")
        if(loc==1){
            qf.y2.mle = rep(NA,n_test)
            for(i in 1:n_test){
                if(i%%100==0)
                    print(i)
                if(i<=n0_test)
                    x_pred = c(1,y1_test_scaled[i])
                if(i>n0_test)
                    x_pred = c(1,qf.y1.mle[i])
                qf.y2.mle[i] <- predict(fit.y2.mle,   X = x_pred, type = "QF",tau=qout21[i])
            } 
        }
        
        if(loc>1){
            X2_pred = cbind(X1_pred_scaled,qf.y1.mle)
            nx3 <- ncol(X2_pred)+1
            for(k in nns){
                vecname = paste0('fits/',region,'_validation/fits_prcp_m',mnth,'_l',k,'.RDS')
                x.vec = readRDS(vecname)
                X2_pred = cbind(X2_pred,x.vec)
            }
            
            # need to scale the x.vec
            X2_pred_scaled <- X2_pred
            for(i in nx3:ncol(X2_train)){
                X2_pred_scaled[,i] <- (X2_pred[,i] - x2_range[1,i])/diff(x2_range[,i])
            }
            
            qf.y2.mle = rep(NA,n_test)
            X2_pred_scaled[,1] <- 1

            for(i in 1:n_test){
                if(i%%100==0)
                    print(i)
                qf.y2.mle[i] <- predict(fit.y2.mle,   X = X2_pred_scaled[i,], type = "QF",tau=qout21[i])
            }   
        }
        y2_pred <- qf.y2.mle*diff(y2_range) + y2_range[1]
        saveRDS(y2_pred,file = predname2)
        
        pdf(file = pdfname,width = 6,height = 6)
        par(mfrow=c(2,2))
        
        plot(y1_test,y1_pred,col=y0_test+1, main = 'MLE-TS',pch=20,cex=0.2)
        abline(0,1)
        
        summary(y1_pred[y0_test==0])
        summary(y1_pred[y0_test==1])
        
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
        

        plot(y2_test,y2_pred,col=y0_test+1, main = 'MLE-TS',pch=20,cex=0.2)
        abline(0,1)
        
        
        summary(y2_pred[y0_test==0])
        summary(y2_pred[y0_test==1])
        
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
        
        par(mfrow=c(1,1))
        dev.off()
    }
