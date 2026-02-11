rm(list = ls())
library(ggplot2)
library(scales)
library(lubridate)
library(transport)
library(GpGp)
library(Metrics)
train = T
dataset=8

est_matrix <- array(NA,dim = c(8,3,100))
se_matrix <- array(NA,dim = c(8,3,100))

for(method in 1:3)
for(dataset in 1:100){
    sim = c('sim_marg','sim_cross','sim_SPCDE')[method]
    dir <- file.path('fits/',method,dataset) 
    filename <- paste0('data/simdata/',dataset,'.RData')
    load(file = filename)

coords = as.matrix(locs)
head(coords)

set.seed(303)
vecchia.order = order_maxmin(coords,lonlat = F)
NNarray <- find_ordered_nn(coords[vecchia.order,],lonlat = F,m=5)

Temp0 <- Temp0[,vecchia.order]
Temp1 <- Temp1[,vecchia.order]
Prec0 <- Prec0[,vecchia.order]
Prec1 <- Prec1[,vecchia.order]
train_indices <- c(1:1500,1921:3420)
test_indices <- c(1501:1920,3421:3840)

y1y2.cors.0 = NA
y1y2.cors.1 = NA
y1y2.cors.2 = NA
daysinmonth = c(30)

cal.data = vector('list',1)
for(mnth in 1:1){
    if(train)
    cal.array = array(dim = c(daysinmonth[mnth]*50,6,25))
    if(!train)
        cal.array = array(dim = c(daysinmonth[mnth]*14,6,25))
    for(loc in 1:25){
        if(train){
            y1 <- c(Temp0[,loc],Temp1[,loc])[train_indices]
            y2 <- c(Prec0[,loc],Prec1[,loc])[train_indices]
            envname = paste0( 'fits/',sim,'/',dataset,'/fits_temp_l',loc,'.RDS')
            qf.y1.mle.ts <- readRDS(envname)[train_indices]
            envname = paste0( 'fits/',sim,'/',dataset,'/fits_prcp_l',loc,'.RDS')
            qf.y2.mle.ts <- readRDS(envname)[train_indices]
            qf.y2.mle.ts <- exp(qf.y2.mle.ts) - 0.0001
        }
        
        if(!train){
            y1 <- c(Temp0[,loc],Temp1[,loc])[test_indices]
            y2 <- c(Prec0[,loc],Prec1[,loc])[test_indices]
            envname = paste0( 'fits/',sim,'/',dataset,'/fits_temp_l',loc,'.RDS')
            qf.y1.mle.ts <- readRDS(envname)[test_indices]
            envname = paste0( 'fits/',sim,'/',dataset,'/fits_prcp_l',loc,'.RDS')
            qf.y2.mle.ts <- readRDS(envname)[test_indices]
            qf.y2.mle.ts <- exp(qf.y2.mle.ts) - 0.0001
        }
       
        # y2 <- log(0.0001+y2)
        n0 <- n1 <- length(y1)/2
        n = n0 + n1
        y0 <- rep(1:0,each=n0)


        y1y2.cors.1 = c(y1y2.cors.1,cor(y1[y0==1],y2[y0==1]))
        y1y2.cors.2 = c(y1y2.cors.2,cor(y1[y0==0],y2[y0==0]))
        y1y2.cors.0 = c(y1y2.cors.0,cor(qf.y1.mle.ts[y0==1],qf.y2.mle.ts[y0==1]))

        cal.array[,1,loc] = y1[y0==0]           #model
        cal.array[,2,loc] = qf.y1.mle.ts[y0==0] #calibrated
        cal.array[,3,loc] = y1[y0==1]           #obs
        cal.array[,4,loc] = y2[y0==0]           #model
        cal.array[,5,loc] = qf.y2.mle.ts[y0==0] #calibrated
        cal.array[,6,loc] = y2[y0==1]           #obs
    }
    cal.data[[mnth]] = cal.array
}
metrics_all <- rep(NA,10)
metrics_se <- rep(NA,10)
cal.array = do.call(abind::abind,c(cal.data,along=1))
cal.array2 = apply(cal.array, 2, c)
summary(cal.array2)
# pdf(paste0('plots/density_',region,'_validation.pdf'),width = 8, height = 4)
# par(mfrow=c(1,2),mgp=c(2.25,0.75,0),mar=c(4,4,1,1))
# d0 <-density(cal.array2[,1]) # gcm
# d1 <-density(cal.array2[,3]) # obs 
# d2 <- density(cal.array2[,2]) # pred
# plotmax.y = max(d0$y,d1$y,d2$y)
# plotmin.y = min(d0$y,d1$y,d2$y)
# plotmax.x = max(d0$x,d1$x,d2$x)
# plotmin.x = min(d0$x,d1$x,d2$x)
# plot(d0,col=2,ylim=range(c(plotmin.y,plotmax.y)),
#      xlim=range(c(plotmin.x,plotmax.x)),ylab="Density",xlab='TMAX')
# lines(d1,col=1)
# lines(d2,col=1,lty=2)
# legend('topleft',c('Mod','Obs','Cal'),col=c(2,1,1),lty = c(1,1,2),lwd=2)
# d0 <-density(log(0.0001 + cal.array2[,4])) # gcm
# d1 <-density(log(0.0001 + cal.array2[,6])) # obs 
# d2 <-density(log(0.0001 + cal.array2[,5])) # pred
# plotmax.y = max(d0$y,d1$y,d2$y)
# plotmin.y = min(d0$y,d1$y,d2$y)
# plotmax.x = max(d0$x,d1$x,d2$x)
# plotmin.x = min(d0$x,d1$x,d2$x)
# plot(d0,col=2,ylim=range(c(plotmin.y,plotmax.y)),
#      xlim=range(c(plotmin.x,plotmax.x)),ylab="Density",xlab = 'PRCP')
# lines(d1,col=1)
# lines(d2,col=1,lty=2)
# legend('topright',c('Mod','Obs','Cal'),col=c(2,1,1),lty = c(1,1,2),lwd=2)
# par(mfrow=c(1,1))
# dev.off()

wasdist = array(dim = c(25,2))
for(loc in 1:25){
    wasdist[loc,1] <- wasserstein1d(cal.array[,2,loc],cal.array[,3,loc])
    wasdist[loc,2] <- wasserstein1d(cal.array[,5,loc],cal.array[,6,loc])
}
summary(wasdist)
metrics_all[c(1,5)] = apply(wasdist,2,mean)
metrics_se[c(1,5)] = apply(wasdist,2,sd)

# pdf(paste0('plots/crosscorr_',model.type,'_',region,'_validation.pdf'),width = 5,height = 4)
# par(mfrow=c(1,1),mgp=c(2.25,0.75,0),mar=c(4,4,1,1))
# lim_min = round(min(c(y1y2.cors.0,y1y2.cors.1,y1y2.cors.2),na.rm = T),5)
# lim_max = round(max(c(y1y2.cors.0,y1y2.cors.1,y1y2.cors.2),na.rm = T),5)
# plot(y1y2.cors.0,y1y2.cors.1,pch=20,xlab = 'Model',ylab = 'Observed',cex=0.75,col=1,
#      xlim = c(lim_min,lim_max),ylim = c(lim_min,lim_max))
# points(y1y2.cors.2,y1y2.cors.1,pch=1,col=2,cex=0.75)
# points(y1y2.cors.0,y1y2.cors.1,pch=20,cex=0.75,col=1)
# legend('topleft',c('Uncalibrated','Calibrated'),pch = c(1,20),col=c(2,1))
# abline(0,1)
# dev.off()

## rmse of autocorrelation and cross correlations
# metrics_all[2] = sqrt(mean((y1.cors.0[-1] - y1.cors.1[-1])**2,na.rm = T))
# metrics_all[6] = sqrt(mean((y2.cors.0[-1] - y2.cors.1[-1])**2,na.rm = T))
### rmse of cross correlations
metrics_all[9] = mae(y1y2.cors.0[-1],y1y2.cors.1[-1])
metrics_se[9] = sd(y1y2.cors.0[-1]-y1y2.cors.1[-1])


q1 = 0.95
q2 = 0.95

joint.tail = array(NA,dim = c(25,1,3))
prcp.tail = array(NA,dim = c(25,1,3))
temp.tail = array(NA,dim = c(25,1,3))
for(mnth in 1:1)
    for(loc in 1:25){
        cal.array = cal.data[[mnth]] 
        temp.q = quantile(cal.array[,3,loc],q1)
        prcp.q = quantile(cal.array[,6,loc],q2)
        for(i in 1:3){
            joint.tail[loc,mnth,i] = mean(cal.array[,i,loc]<temp.q & cal.array[,i+3,loc]>prcp.q)    
            prcp.tail[loc,mnth,i] = mean(cal.array[,i+3,loc]>prcp.q)    
            temp.tail[loc,mnth,i] = mean(cal.array[,i,loc]>temp.q)    
        }
        
        
    }

# compare means and upper quantiles
# 1 = obs 2 = pred 
pred_summaries = matrix(NA,25,9)
count = 0
for(mnth in 1:1)
    for(loc in 1:25){
        count = count+1
        cal.array = cal.data[[mnth]] 
        # obs summary
        pred_summaries[count,1] = mean(cal.array[,6,loc])
        pred_summaries[count,4] = quantile(cal.array[,6,loc],q2)
        pred_summaries[count,7] = quantile(cal.array[,3,loc],q1)
        # pred summary
        pred_summaries[count,2] = mean(cal.array[,5,loc])
        pred_summaries[count,5] = quantile(cal.array[,5,loc],q2)
        pred_summaries[count,8] = quantile(cal.array[,2,loc],q1)
        # gcm summary
        pred_summaries[count,3] = mean(cal.array[,4,loc])
        pred_summaries[count,6] = quantile(cal.array[,4,loc],q2)
        pred_summaries[count,9] = quantile(cal.array[,1,loc],q1)
    }
# pdf(paste0('plots/summaries_',region,'_validation.pdf'),width = 8,height = 4)
# par(mfrow=c(1,2),mgp=c(2.25,0.75,0),mar=c(4,4,1,1))
# 
# lim_min = floor(min(pred_summaries[,7:9],na.rm = T))
# lim_max = ceiling(max(pred_summaries[,7:9],na.rm = T))
# plot(pred_summaries[,8:7],xlab='Model',ylab = 'Observed',pch=20,
#      ylim=c(lim_min,lim_max),xlim=c(lim_min,lim_max),cex=0.75)
# points(pred_summaries[,c(7,9)],pch=1,col=2,cex=0.75)
# points(pred_summaries[,8:7],pch=20,col=1,cex=0.75)
# abline(0,1)
# legend('topleft',c('Uncalibrated','Calibrated'),pch = c(1,20),col=c(2,1))
# 
# lim_min = floor(min(pred_summaries[,4:6],na.rm = T))
# lim_max = ceiling(max(pred_summaries[,4:6],na.rm = T))
# plot(pred_summaries[,5:4],xlab='Model',ylab = 'Observed',pch=20,
#      ylim=c(lim_min,lim_max),xlim=c(lim_min,lim_max),cex=0.75)
# points(pred_summaries[,c(6,4)],pch=1,col=2,cex=0.75)
# points(pred_summaries[,5:4],pch=20,cex=0.75,col=1)
# abline(0,1)
# legend('topleft',c('Uncalibrated','Calibrated'),pch = c(1,20),col=c(2,1))
# par(mfrow=c(1,1))
# dev.off()

# rmse of upper quantiles
metrics_all[8] = mae(pred_summaries[,4],pred_summaries[,5]) # prcp
metrics_all[4] = mae(pred_summaries[,7],pred_summaries[,8]) # tmax
metrics_se[8] = sd(pred_summaries[,4]-pred_summaries[,5]) # prcp
metrics_se[4] = sd(pred_summaries[,7]-pred_summaries[,8]) # tmax

correls = matrix(NA,300,8)
count = 0
for(mnth in 1:1){
    cal.array = cal.data[[mnth]]
    # print(mnth)
    for(i in 1:24)
        for(j in (i+1):25){
            count = count+1
            loc1 = coords[vecchia.order==i,-3]
            loc2 = coords[vecchia.order==j,-3]
            correls[count,7] = as.numeric(sqrt((loc1[1]-loc2[1])**2 + (loc1[2]-loc2[2])**2))
            correls[count,8] = mnth
            for(k in 1:6)
                correls[count,k] = cor(cal.array[,k,c(i,j)])[1,2]
        }
}

correls2 = correls
# pdf(paste0('plots/spatcorr_',model.type,'_',region,'_validation.pdf'),width = 8,height = 4)
# par(mfrow=c(1,2),mgp=c(2.25,0.75,0),mar=c(4,4,1,1))
# plot(correls2[,2],correls2[,3],pch=20,col=1,
#      xlab = 'Model',ylab = 'Observed')
# points(correls2[,3],correls2[,1],col=2,pch=1,cex=0.75)
# points(correls2[,2],correls2[,3],pch=20,col=1,cex=0.75)
# abline(0,1)
# legend('topleft',c('Uncalibrated','Calibrated'),pch = c(1,20),col=c(2,1))
# 
# plot(correls2[,5],correls2[,6],col=1,pch=20,
#      xlab = 'Model',ylab = 'Observed')
# points(correls2[,6],correls2[,4],col=2,pch=1,cex=0.75)
# points(correls2[,5],correls2[,6],col=1,pch=20,cex=0.75)
# abline(0,1)
# legend('topleft',c('Uncalibrated','Calibrated'),pch = c(1,20),col = c(2,1))
# dev.off()

# spatial correlations RMSE
metrics_all[3] = mae(correls[,2],correls[,3])
metrics_all[7] = mae(correls[,5],correls[,6])
metrics_se[3] = sd(correls[,2]-correls[,3])
metrics_se[7] = sd(correls[,5]-correls[,6])
metrics_all
round(metrics_all,4)


count = 0
propzero = matrix(NA,25,3)
length(cal.data)
for(mnth in 1:1)
    for(loc in 1:25){
        count = count+1
        tmp <- cal.data[[mnth]][,4:6,loc]
        # tmp <- exp(tmp) - 0.0001
        propzero[count,] <- apply(tmp,2,function(x)mean(round(x,4)==0))        
    }

# pdf(paste0('plots/propzero_',region,'_validation.pdf'),width = 5,height = 4)
# plot(propzero[,c(2,3)],pch=20,col=1,xlab = 'Model',ylab = 'Observed',cex=0.75)
# abline(0,1)
# points(propzero[,c(1,3)],pch=1,col=2,cex=0.75)
# points(propzero[,c(2,3)],pch=20,col=1,cex=0.75)
# legend('bottomright',c('Uncalibrated','Calibrated'),pch = c(1,20),col=c(2,1))
# dev.off()

metrics_all[10] <- mae(propzero[,3],propzero[,2])
metrics_se[10] <- sd(propzero[,3]-propzero[,2])


#wasdist,autocorr,spatcorr,quantile,wasdist,autocorr,spatcorr,quantile,crosscorr,propzero

est_matrix[, method,dataset] <- metrics_all[c(1,4,3,5,8,7,10,9)]
se_matrix[, method,dataset] <- metrics_se[c(1,4,3,5,8,7,10,9)]

print(paste(method,dataset))
}

names(metrics_all) <- c('wasdist','quantile','spatcorr',
                        'wasdist','quantile','spatcorr',
                        'propzero','crosscorr')

dim(est_matrix)
par(mfrow=c(2,4))
for(i in 1:8)
boxplot(t(est_matrix[i,,]),main=names(metrics_all)[i])
par(mfrow=c(1,1))

save(est_matrix,se_matrix,file = 'results/sim_SPCDE_results_train.RData')
