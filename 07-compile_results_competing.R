rm(list = ls())
library(ggplot2)
library(scales)
library(lubridate)
library(GpGp)
library(usmap)
library(transport)
library(gridExtra)
region = 'SW'
method = 'CCA'
model.type = 'space'
gcm.long = read.csv(paste0('data/',region,'_gcm_data.csv'))
obs.long = read.csv(paste0('data/',region,'_obs_data.csv'))

if(region=='SW')
    pred.long <- read.csv('data/West_data.csv')
if(region=='SE')
    pred.long <- read.csv('data/East_data.csv')

gcm.months = month(gcm.long[,1])
gcm.years = year(gcm.long[,1])
pred.months = month(pred.long[,1])

grid.no = as.factor(gcm.long$lat*gcm.long$lon)
str(grid.no)
coords = aggregate(gcm.long[,2:3],by = list(grid.no), FUN = mean)
coords = coords[,-1]
head(coords)
coords$lon = coords$lon - 360
GeoLocations <- usmap_transform(coords)

# p1 = plot_usmap(regions ="region", include = c('UT','CO','AZ','NM')) + geom_sf(data = GeoLocations)
# p2 = plot_usmap(regions ="region", include = 
#                     c('Indiana','Ohio','West Virginia','Kentucky',
#                       'Virginia','Tennessee','North Carolina', 'South Carolina',
#                       'Mississippi','Alabama','Georgia')) + geom_sf(data = GeoLocations)

table(grid.no)
set.seed(303)
vecchia.order = order_maxmin(coords,lonlat = T)
y1.cors.0 = NA
y2.cors.0 = NA
y1y2.cors.0 = NA
y1.cors.1 = NA
y2.cors.1 = NA
y1y2.cors.1 = NA
daysinmonth = c(31,28,31,30,31,30,31,31,30,31,30,31)
mnth = 1; loc = 1
cal.data = vector('list',12)
for(mnth in 1:12){
    cal.array = array(dim = c(daysinmonth[mnth]*64,6,25))
    for(loc in 1:25){
        print(paste(mnth,loc))
        y1 <- c(obs.long$tmax[vecchia.order==loc & gcm.months==mnth],
                gcm.long$tmax[vecchia.order==loc & gcm.months==mnth])
        y2 <- c(obs.long$pr[vecchia.order==loc & gcm.months==mnth],
                gcm.long$pr[vecchia.order==loc & gcm.months==mnth])
        # y2 <- log(0.0001+y2)
        n0 = length(gcm.long$pr[vecchia.order==loc & gcm.months==mnth])
        n1 = length(obs.long$pr[vecchia.order==loc & gcm.months==mnth])
        n = n0 + n1
        y0 <- rep(1:0,each=n0)

        y11 = y1[y0==1]
        x11 = y11[c(n1,1:(n1-1))]
        y21 = y2[y0==1]
        x21 = y21[c(n1,1:(n1-1))]


        if(method=='CCA'){
            qf.y1.mle.ts <- c(gcm.long$tmax[vecchia.order==loc & gcm.months==mnth],
                              pred.long$tmax_CCA[vecchia.order==loc & pred.months==mnth])

            qf.y2.mle.ts <- c(gcm.long$pr[vecchia.order==loc & gcm.months==mnth],
                              pred.long$pr_CCA[vecchia.order==loc & pred.months==mnth])
        }
        if(method=='QM'){
            qf.y1.mle.ts <- c(gcm.long$tmax[vecchia.order==loc & gcm.months==mnth],
                              pred.long$tmax_QR[vecchia.order==loc & pred.months==mnth])

            qf.y2.mle.ts <- c(gcm.long$pr[vecchia.order==loc & gcm.months==mnth],
                              pred.long$pr_QR[vecchia.order==loc & pred.months==mnth])
        }

        # qf.y2.mle.ts <- log(0.0001+qf.y2.mle.ts)

        y1.cors.0 = c(y1.cors.0,cor(qf.y1.mle.ts[y0==0][-1],qf.y1.mle.ts[y0==0][-n0]))
        y1.cors.1 = c(y1.cors.1,cor(cbind(y11,x11))[1,2])
        y2.cors.0 = c(y2.cors.0,cor(qf.y2.mle.ts[y0==0][-1],qf.y2.mle.ts[y0==0][-n0]))
        y2.cors.1 = c(y2.cors.1,cor(cbind(y21,x21))[1,2])
        y1y2.cors.1 = c(y1y2.cors.1,cor(y1[y0==1],y2[y0==1]))
        y1y2.cors.0 = c(y1y2.cors.0,cor(qf.y1.mle.ts[y0==1],qf.y2.mle.ts[y0==1]))

        cal.array[,1,loc] = y1[y0==0]
        cal.array[,3,loc] = y1[y0==1]
        cal.array[,2,loc] = qf.y1.mle.ts[y0==0]
        cal.array[,4,loc] = y2[y0==0]
        cal.array[,6,loc] = y2[y0==1]
        cal.array[,5,loc] = qf.y2.mle.ts[y0==0]
    }
    cal.data[[mnth]] = cal.array
}

save(y1.cors.0,y1.cors.1,y2.cors.0,y2.cors.1,y1y2.cors.0,y1y2.cors.1,cal.data,
           file = paste0('summary_',model.type,'_',region,'_',method,'.RData'))
load(paste0('summary_',model.type,'_',region,'_',method,'.RData'))
metrics_all <- rep(NA,10)
eachmonth = rep(NA,12)
# for(i in 1:12){
#     eachmonth[i] = dim(cal.data[[i]])[1]
#     cal.data[[i]][,4:6,] = exp(cal.data[[i]][,4:6,])-1
# }
wasdist = array(dim = c(25,2))
# for(mnth in 1:12){
#     cal.array = cal.data[[mnth]]
#     for(loc in 1:25){
#         wasdist[loc,1,mnth] <- wasserstein1d(cal.array[,2,loc],cal.array[,3,loc])
#         wasdist[loc,2,mnth] <- wasserstein1d(cal.array[,5,loc],cal.array[,6,loc])
#     }
#     
# }

cal.array = do.call(abind::abind,c(cal.data,along=1))
cal.array2 = apply(cal.array, 2, c)

# png(paste0('plots/',region,'_density.png'),width = 800, height = 400)
par(mfrow=c(1,2))
d0 <-density(cal.array2[,1]) # gcm
d1 <-density(cal.array2[,3]) # obs 
d2 <- density(cal.array2[,2]) # pred
plotmax.y = max(d0$y,d1$y,d2$y)
plotmin.y = min(d0$y,d1$y,d2$y)
plotmax.x = max(d0$x,d1$x,d2$x)
plotmin.x = min(d0$x,d1$x,d2$x)
plot(d0,col=1,ylim=range(c(plotmin.y,plotmax.y)),
     xlim=range(c(plotmin.x,plotmax.x)),ylab="Density",main="Temp")
lines(d1,col=2)
lines(d2,col=2,lty=2)
legend('topleft',c('GCM','Obs','Cal'),col=c(1,2,2),lty = c(1,1,2),lwd=2)
d0 <-density(log(0.0001+cal.array2[,4])) # gcm
d1 <-density(log(0.0001+cal.array2[,6])) # obs 
d2 <- density(log(0.0001+cal.array2[,5])) # pred
plotmax.y = max(d0$y,d1$y,d2$y)
plotmin.y = min(d0$y,d1$y,d2$y)
plotmax.x = max(d0$x,d1$x,d2$x)
plotmin.x = min(d0$x,d1$x,d2$x)
plot(d0,col=1,ylim=range(c(plotmin.y,plotmax.y)),
     xlim=range(c(plotmin.x,plotmax.x)),ylab="Density",main="Prcp")
lines(d1,col=2)
lines(d2,col=2,lty=2)
legend('topright',c('GCM','Obs','Cal'),col=c(1,2,2),lty = c(1,1,2),lwd=2)
par(mfrow=c(1,1))
# dev.off()

for(loc in 1:25){
    wasdist[loc,1] <- wasserstein1d(cal.array[,2,loc],cal.array[,3,loc])
    wasdist[loc,2] <- wasserstein1d(cal.array[,5,loc],cal.array[,6,loc])
}
summary(wasdist)
metrics_all[c(1,5)] = apply(wasdist,2,mean)

metrics = data.frame(coords,wasdist,vecchia.order)
ggplot(metrics,aes(x=lon,y=lat,fill=X1)) + geom_raster() + coord_equal() +
    geom_text(aes(label=round(X1,2)),col='white') + ggtitle(paste(region,'temp')) +
    theme(legend.position = 'none')
ggplot(metrics,aes(x=lon,y=lat,fill=X2)) + geom_raster() + coord_equal() +
    geom_text(aes(label=round(X2,2)),col='white') + ggtitle(paste(region,'prcp')) +
    theme(legend.position = 'none')

# mnth=12
# loc=1
# maxy = 1.01*max(cal.array[1:31,4:6,loc])
# miny = 0.99*min(cal.array[1:31,4:6,loc])
# plot(1:31,cal.array[1:31,4,loc],type = 'b',col=1,pch=20,ylim = c(miny,maxy))
# lines(1:31,cal.array[1:31,5,loc],type = 'b',col=2,pch=20)
# lines(1:31,cal.array[1:31,6,loc],type = 'b',col=3,pch=20)

season = rep(1:4,75)
# png(paste0('plots/autocorr_',model.type,'_',region,'_validation.png'),width = 800,height = 400)
par(mfrow=c(1,2))
plot(y1.cors.0[-1],y1.cors.1[-1],pch=20,xlab = 'calibrated',ylab = 'observed',
     col=season,main = 'temp autocorrelations',xlim = c(0,1),ylim = c(0,1))
legend('topleft',c('JFM','AMJ','JAS','OND'),col=1:4,pch=20,pt.cex = 2)
abline(0,1)
plot(y2.cors.0[-1],y2.cors.1[-1],pch=20,xlab = 'calibrated',ylab = 'observed',
     col = season, main = 'prcp autocorrelations',xlim = c(0,1),ylim = c(0,1))
legend('topleft',c('JFM','AMJ','JAS','OND'),col=1:4,pch=20,pt.cex = 2)
abline(0,1)
# plot(y1y2.cors.0,y1y2.cors.1,pch=20,xlab = 'calibrated',ylab = 'observed',main = 'cross correlations')
# abline(0,1)
par(mfrow=c(1,1))
# dev.off()
# png(paste0('plots/crosscorr_',model.type,'_',region,'_validation.png'),width = 400,height = 400)
plot(y1y2.cors.0,y1y2.cors.1,pch=20,xlab = 'calibrated',ylab = 'observed',
     col=season,main = 'cross correlations')
abline(0,1)
# dev.off()
## rmse of autocorrelation and cross correlations
metrics_all[2] = sqrt(mean((y1.cors.0[-1] - y1.cors.1[-1])**2,na.rm = T))
metrics_all[6] = sqrt(mean((y2.cors.0[-1] - y2.cors.1[-1])**2,na.rm = T))
### rmse of cross correlations
metrics_all[9] = sqrt(mean((y1y2.cors.0[-1] - y1y2.cors.1[-1])**2))

q1 = 0.95
q2 = 0.95

joint.tail = array(NA,dim = c(25,12,3))
prcp.tail = array(NA,dim = c(25,12,3))
temp.tail = array(NA,dim = c(25,12,3))
for(mnth in 1:12)
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
pred_summaries = matrix(NA,300,9)
count = 0
for(mnth in 1:12)
    for(loc in 1:25){
        count = count+1
        cal.array = cal.data[[mnth]] 
        # obs summary
        pred_summaries[count,1] = mean(cal.array[,6,loc],na.rm = T)
        pred_summaries[count,4] = quantile(cal.array[,6,loc],q2,na.rm = T)
        pred_summaries[count,7] = quantile(cal.array[,3,loc],q1,na.rm = T)
        # pred summary
        pred_summaries[count,2] = mean(cal.array[,5,loc],na.rm = T)
        pred_summaries[count,5] = quantile(cal.array[,5,loc],q2,na.rm = T)
        pred_summaries[count,8] = quantile(cal.array[,2,loc],q1,na.rm = T)
        # gcm summary
        pred_summaries[count,3] = mean(cal.array[,4,loc],na.rm = T)
        pred_summaries[count,6] = quantile(cal.array[,4,loc],q2,na.rm = T)
        pred_summaries[count,9] = quantile(cal.array[,1,loc],q1,na.rm = T)
    }
# png(paste0('plots/summaries_',region,'_validation.png'),width = 1200,height = 400)
par(mfrow=c(1,3))
plot(pred_summaries[,1:2],col=season,xlab='Observed',ylab = 'Calibrated',pch=20,main = 'Prcp mean')
abline(0,1)
points(pred_summaries[,c(1,3)],pch=3,col=alpha(1,0.4))
legend('topleft',c('Uncalibrated','Calibrated'),pch = c(3,20))

plot(pred_summaries[,4:5],col=season,xlab='Observed',ylab = 'Calibrated',pch=20,main = 'Prcp 0.90 quantile')
abline(0,1)
points(pred_summaries[,c(4,6)],pch=3,col=alpha(1,0.4))
legend('topleft',c('Uncalibrated','Calibrated'),pch = c(3,20))

plot(pred_summaries[,7:8],col=season,xlab='Observed',ylab = 'Calibrated',pch=20,main = 'Temp 0.90 quantile')
abline(0,1)
points(pred_summaries[,c(7,9)],pch=3,col=alpha(1,0.4))
legend('topleft',c('Uncalibrated','Calibrated'),pch = c(3,20))
par(mfrow=c(1,1))
# dev.off()

# rmse of upper quantiles
# metrics_all[8] = sqrt(mean((pred_summaries[,4]-pred_summaries[,5])**2,na.rm = T)) # prcp
# metrics_all[4] = sqrt(mean((pred_summaries[,7]-pred_summaries[,8])**2,na.rm = T)) # tmax
metrics_all[8] = mean(abs((pred_summaries[,4]-pred_summaries[,5]))) # prcp
metrics_all[4] = mean(abs((pred_summaries[,7]-pred_summaries[,8]))) # tmax
# png('tailprob_space.png',width = 1200,height = 400)
# par(mfrow=c(1,3))
# tmp = apply(prcp.tail,3,c)
# dim(tmp)
# xmax = max(tmp)
# xmin = min(tmp)
# plot(density(tmp[,2]),xlab = 'probability',main = 'precip',xlim = c(xmin,xmax))
# lines(density(tmp[,1]),col=2)
# abline(v=1-q2)
# legend('topright',c('Calibrated','GCM'),col=1:2,lwd=2,bty="n")
# 
# 
# tmp = apply(temp.tail,3,c)
# dim(tmp)
# xmax = max(tmp)
# xmin = min(tmp)
# plot(density(tmp[,1]),xlab = 'probability',main = 'temp',xlim = c(xmin,xmax),col=2)
# lines(density(tmp[,2]))
# abline(v=1-q1)
# legend('topright',c('Calibrated','GCM'),col=1:2,lwd=2,bty="n")
# 
# tmp = apply(joint.tail,3,c)
# dim(tmp)
# xmax = max(tmp)
# xmin = min(tmp)
# plot(density(tmp[,2]),xlab = 'probability',main = 'joint',xlim = c(xmin,xmax),col=2)
# lines(density(tmp[,1]))
# lines(density(tmp[,3]))
# abline(v=1-q1)
# legend('topright',c('Calibrated','GCM'),col=1:2,lwd=2,bty="n")
# par(mfrow=c(1,1))
# dev.off()
# png(paste0('plots/tailprob_',model.type,'_',region,'_validation.png'),width = 800,height = 400)
par(mfrow=c(1,2))
tmp = apply(temp.tail,3,c)
plot(tmp[,1],tmp[,2],xlab = 'uncalibrated',ylab = 'calibrated',pch=20,
     main = 'Temp exceedance probability above 0.90 quantile',col=season)
legend('topright',c('JFM','AMJ','JAS','OND'),col=1:4,pch=20,pt.cex = 2)
abline(h=0.1)

tmp = apply(prcp.tail,3,c)
plot(tmp[,1],tmp[,2],xlab = 'uncalibrated',ylab = 'calibrated',pch=20,
     main = 'Prcp exceedance probability above 0.90 quantile',col=season)
legend('topleft',c('JFM','AMJ','JAS','OND'),col=1:4,pch=20,pt.cex = 2)
abline(h=0.1)
# dev.off()
# abline(0,1)
# 
# coords = cbind(coords,vecchia.order)
# names(coords)
# ggplot(coords,aes(x=lon,y=lat,col='white'))+geom_tile()+geom_text(aes(label=vecchia.order))

correls = matrix(NA,300*12,8)
count = 0
for(mnth in 1:12){
cal.array = cal.data[[mnth]]
print(mnth)
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
season = rep(1:4,each=900)
mnth = rep(1:300,each=12)
correls2 = aggregate(correls,by=list(mnth),FUN=mean)
correls2 = correls2[,-1]
# png(paste0('plots/spatcorr_',model.type,'_',region,'_validation.png'),width = 800,height = 400)
par(mfrow=c(1,2))
plot(correls[,2],correls[,1],col=alpha(1,0.4),pch=3,
     xlab = 'Calibrated',ylab = 'Observed',main = 'Temp spatial correlations')
points(correls[,2],correls[,3],col=alpha(season,0.2),pch=20)
abline(0,1)
legend('topleft',c('GCM','Obs'),pch = c(3,20))

plot(correls[,5],correls[,4],col=alpha(1,0.4),pch=3,
     xlab = 'Calibrated',ylab = 'Observed',main = 'Prcp spatial correlations')
points(correls[,5],correls[,6],col=alpha(season,0.2),pch=20)
abline(0,1)
legend('topleft',c('GCM','Obs'),pch = c(3,20))
# dev.off()

# spatial correlations RMSE
# metrics_all[3] = sqrt(mean((correls[,2]-correls[,1])**2,na.rm = T))
# metrics_all[7] = sqrt(mean((correls[,5]-correls[,4])**2,na.rm = T))
metrics_all[3] = mean(abs((correls[,2]-correls[,3])))
metrics_all[7] = mean(abs((correls[,5]-correls[,6])))

count = 0
propzero = matrix(NA,300,3)
length(cal.data)
for(mnth in 1:12)
    for(loc in 1:25){
        count = count+1
        propzero[count,] <- apply(cal.data[[mnth]][,4:6,loc],2,function(x)mean(round(x,4)==0))        
    }


# metrics_all[10] <- sqrt(mean((propzero[,3]-propzero[,2])**2))
metrics_all[10] <- mean(abs((propzero[,3]-propzero[,2])))

names(metrics_all) <- c('wasdist','autocorr','spatcorr','quantile',
                        'wasdist','autocorr','spatcorr','quantile',
                        'crosscorr','propzero')
round(metrics_all,4)
round(metrics_all[c(1,4,2,3,5,8,10,6,7,9)],4)
