rm(list = ls())
library(ggplot2)
library(scales)
library(lubridate)
library(Metrics)
dataset<-8
est_matrix <- array(NA,dim = c(8,2,100))
se_matrix <- array(NA,dim = c(8,2,100))



train<-T
for(type in 1:2)
for(dataset in 1:100){
    method <- c('QM','CCE')[type]
    sim <- c('sim_marg','sim_cross','sim_SPCDE')[method]
    dir <- file.path('fits/',method,dataset) 
    filename <- paste0('data/simdata/',dataset,'.RData')
    load(file = filename)
    pred.name <- paste0('fits/sim_competing/bias_corrected_file_',sprintf("%03d", dataset),'.csv')
pred.long <- read.csv(pred.name)
coords <- as.matrix(locs)
head(coords)

set.seed(303)
vecchia.order <- order_maxmin(coords,lonlat = F)
NNarray <- find_ordered_nn(coords[vecchia.order,],lonlat = F,m=5)

Temp0 <- Temp0[,vecchia.order]
Temp1 <- Temp1[,vecchia.order]
Prec0 <- Prec0[,vecchia.order]
Prec1 <- Prec1[,vecchia.order]

train_indices   <- c(1:1500,1921:3420)
test_indices    <- c(1501:1920,3421:3840)

y1y2.cors.0 <- NA
y1y2.cors.1 <- NA
y1y2.cors.2 <- NA
daysinmonth <- 30

mnth <- 1
loc  <-1
cal.data <- vector('list',1)
for(mnth in 1:1){
    if(train)
        cal.array   <- array(dim = c(daysinmonth[mnth]*50,6,25))
    if(!train)
        cal.array   <- array(dim = c(daysinmonth[mnth]*14,6,25))
    for(loc in 1:25){
        cur.loc     <- vecchia.order[loc]
       if(train){
           y1   <- c(Temp0[,loc],Temp1[,loc])[train_indices]
           y2   <- c(Prec0[,loc],Prec1[,loc])[train_indices]
       }
        if(!train){
            y1  <- c(Temp0[,loc],Temp1[,loc])[test_indices]
            y2  <- c(Prec0[,loc],Prec1[,loc])[test_indices]
        }
        n0  <- n1 <- length(y1)/2
        n   <- n0 + n1
        y0  <- rep(1:0,each=n0)
        
        y10 <- y1[y0==0]
        y11 <- y1[y0==1]
        x10 <- y10[c(n1,1:(n1-1))]
        x11 <- y11[c(n1,1:(n1-1))]
        y20 <- y2[y0==0]
        y21 <- y2[y0==1]
        x20 <- y20[c(n1,1:(n1-1))]
        x21 <- y21[c(n1,1:(n1-1))]
        
        if(train){
            if(method=='CCA'){
                qf.y1.mle.ts <-(pred.long$tmax_CCA[vecchia.order==loc])[train_indices]
                qf.y2.mle.ts <- (pred.long$prcp_CCA[vecchia.order==loc])[train_indices]
            }
            if(method=='QM'){
                qf.y1.mle.ts <- (pred.long$tmax_QR[vecchia.order==loc])[train_indices]
                qf.y2.mle.ts <- (pred.long$prcp_QR[vecchia.order==loc])[train_indices]
            }
        }
        if(!train){
            if(method=='CCA'){
                qf.y1.mle.ts <-(pred.long$tmax_CCA[vecchia.order==loc])[test_indices]
                qf.y2.mle.ts <- (pred.long$prcp_CCA[vecchia.order==loc])[test_indices]
            }
            if(method=='QM'){
                qf.y1.mle.ts <- (pred.long$tmax_QR[vecchia.order==loc])[test_indices]
                qf.y2.mle.ts <- (pred.long$prcp_QR[vecchia.order==loc])[test_indices]
            }
        }
        
        
        y1y2.cors.1 <- c(y1y2.cors.1,cor(y1[y0==1],y2[y0==1]))
        y1y2.cors.2 <- c(y1y2.cors.2,cor(y1[y0==0],y2[y0==0]))
        y1y2.cors.0 <- c(y1y2.cors.0,cor(qf.y1.mle.ts[y0==1],qf.y2.mle.ts[y0==1]))
        
        cal.array[,1,loc] <- y1[y0==0]
        cal.array[,3,loc] <- y1[y0==1]
        cal.array[,2,loc] <- qf.y1.mle.ts[y0==1]
        cal.array[,4,loc] <- y2[y0==0]
        cal.array[,6,loc] <- y2[y0==1]
        cal.array[,5,loc] <- qf.y2.mle.ts[y0==1]
    }
    cal.data[[mnth]] <- cal.array
}
metrics_all <- rep(NA,10)
metrics_se  <- rep(NA,10)
wasdist     <- array(dim = c(25,2))

cal.array   <- do.call(abind::abind,c(cal.data,along=1))
cal.array2  <- apply(cal.array, 2, c)

for(loc in 1:25){
    wasdist[loc,1] <- wasserstein1d(cal.array[,2,loc],cal.array[,3,loc])
    wasdist[loc,2] <- wasserstein1d(cal.array[,5,loc],cal.array[,6,loc])
}
metrics_all[c(1,5)] <- apply(wasdist,2,mean)
metrics_se[c(1,5)]  <- apply(wasdist,2,sd)

## rmse of cross correlations
metrics_all[9] <- mae(y1y2.cors.0[-1], y1y2.cors.1[-1])
metrics_se[9]  <- sd(y1y2.cors.0[-1]- y1y2.cors.1[-1])

q1 <- 0.95
q2 <- 0.95

joint.tail <- array(NA,dim = c(25,1,3))
prcp.tail  <- array(NA,dim = c(25,1,3))
temp.tail  <- array(NA,dim = c(25,1,3))
for(mnth in 1:1)
    for(loc in 1:25){
        cal.array   <- cal.data[[mnth]] 
        temp.q      <- quantile(cal.array[,3,loc],q1)
        prcp.q      <- quantile(cal.array[,6,loc],q2)
        for(i in 1:3){
            joint.tail[loc,mnth,i]  <- mean(cal.array[,i,loc]<temp.q & cal.array[,i+3,loc]>prcp.q)    
            prcp.tail[loc,mnth,i]   <- mean(cal.array[,i+3,loc]>prcp.q)    
            temp.tail[loc,mnth,i]   <- mean(cal.array[,i,loc]>temp.q)    
        }
    }

# compare means and upper quantiles
# 1 = obs 2 = pred 
pred_summaries <- matrix(NA,25,9)
count <- 0
for(mnth in 1:1)
    for(loc in 1:25){
        count <- count+1
        cal.array <- cal.data[[mnth]] 
        # obs summary
        pred_summaries[count,1] <- mean(cal.array[,6,loc])
        pred_summaries[count,4] <- quantile(cal.array[,6,loc],q2)
        pred_summaries[count,7] <- quantile(cal.array[,3,loc],q1)
        # pred summary
        pred_summaries[count,2] <- mean(cal.array[,5,loc])
        pred_summaries[count,5] <- quantile(cal.array[,5,loc],q2)
        pred_summaries[count,8] <- quantile(cal.array[,2,loc],q1)
        # gcm summary
        pred_summaries[count,3] <- mean(cal.array[,4,loc])
        pred_summaries[count,6] <- quantile(cal.array[,4,loc],q2)
        pred_summaries[count,9] <- quantile(cal.array[,1,loc],q1)
    }

# rmse of upper quantiles
metrics_all[8] <- mae(pred_summaries[,4],pred_summaries[,5]) # prcp
metrics_all[4] <- mae(pred_summaries[,7],pred_summaries[,8]) # tmax
metrics_se[8]  <- sd(pred_summaries[,4]-pred_summaries[,5]) # prcp
metrics_se[4]  <- sd(pred_summaries[,7]-pred_summaries[,8]) # tmax

correls <- matrix(NA,300,8)
count <- 0
for(mnth in 1:1){
    cal.array <- cal.data[[mnth]]
    # print(mnth)
    for(i in 1:24)
        for(j in (i+1):25){
            count <- count+1
            loc1  <- coords[vecchia.order==i,-3]
            loc2  <- coords[vecchia.order==j,-3]
            correls[count,7] <- as.numeric(sqrt((loc1[1]-loc2[1])**2 + (loc1[2]-loc2[2])**2))
            correls[count,8] <- mnth
            for(k in 1:6)
                correls[count,k] <- cor(cal.array[,k,c(i,j)])[1,2]
        }
}
mnth <- rep(1:300,each=12)
correls2 <- correls

# spatial correlations RMSE
metrics_all[3] <- mae(correls[,2], correls[,3])
metrics_all[7] <- mae(correls[,5],correls[,6])
metrics_se[3]  <- sd(correls[,2]- correls[,3])
metrics_se[7]  <- sd(correls[,5]-correls[,6])
metrics_all
round(metrics_all,4)

count <- 0
propzero <- matrix(NA,25,3)
length(cal.data)
for(mnth in 1:1)
    for(loc in 1:25){
        count <- count+1
        tmp <- cal.data[[mnth]][,4:6,loc]
        propzero[count,] <- apply(tmp,2,function(x)mean(round(x,4)==0))        
    }

metrics_all[10] <- mae(propzero[,3],propzero[,2])
metrics_se[10]  <- sd(propzero[,3]-propzero[,2])


est_matrix[,type,dataset] <- metrics_all[c(1,4,3,5,8,7,10,9)]
se_matrix[,type,dataset]  <- metrics_se[c(1,4,3,5,8,7,10,9)]

print(paste(type,dataset))
}
save(est_matrix,se_matrix,file = 'results/sim_competing_results_train.RData')
