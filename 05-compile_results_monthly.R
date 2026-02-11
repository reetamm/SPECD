rm(list = ls())
library(ggplot2)
library(scales)
library(lubridate)
library(GpGp)
library(usmap)
library(transport)
library(gridExtra)
region = 'SW'
method = 'MLE'
model.type = 'space'
gcm.long = read.csv(paste0('data/',region,'_gcm_data.csv'))
obs.long = read.csv(paste0('data/',region,'_obs_data.csv'))

gcm.months = month(gcm.long[,1])
gcm.years = year(gcm.long[,1])

grid.no = as.factor(gcm.long$lat*gcm.long$lon)
str(grid.no)
coords = aggregate(gcm.long[,2:3],by = list(grid.no), FUN = mean)
coords = coords[,-1]
head(coords)
GeoLocations <- usmap_transform(coords)

# p1 = plot_usmap(regions ="region", include = c('UT','CO','AZ','NM')) + geom_sf(data = GeoLocations)
# p2 = plot_usmap(regions ="region", include = 
#                     c('Indiana','Ohio','West Virginia','Kentucky',
#                       'Virginia','Tennessee','North Carolina', 'South Carolina',
#                       'Mississippi','Alabama','Georgia')) + geom_sf(data = GeoLocations)
m=5
table(grid.no)
set.seed(303)
vecchia.order = order_maxmin(coords,lonlat = T)
NNarray <- find_ordered_nn(coords[vecchia.order,],lonlat = T,m=m)
loc = 3
mnth = 1
mnths = 11:12
loc.vector <- 1:25

y1y2.cors.0 = NA
y1y2.cors.1 = NA
y1y2.cors.2 = NA
daysinmonth = c(31,28,31,30,31,30,31,31,30,31,30,31)

cal.data = vector('list',12)
for(mnth in 1:12){
    cal.array = array(dim = c(daysinmonth[mnth]*50,7,25))
    for(loc in 1:25){
        cur.loc <- vecchia.order[loc]
        print(paste(mnth,loc))
        y1 <- c(obs.long$tmax[loc.vector==cur.loc & gcm.months==mnth & gcm.years <= 2000],
                gcm.long$tmax[loc.vector==cur.loc & gcm.months==mnth & gcm.years <= 2000])
        y2 <- c(obs.long$pr[loc.vector==cur.loc & gcm.months==mnth & gcm.years <= 2000],
                gcm.long$pr[loc.vector==cur.loc & gcm.months==mnth & gcm.years <= 2000])
        n0 = length(gcm.long$pr[loc.vector==cur.loc & gcm.months==mnth & gcm.years <= 2000])
        n1 = length(obs.long$pr[loc.vector==cur.loc & gcm.months==mnth & gcm.years <= 2000])
        n = n0 + n1
        y0 <- rep(1:0,each=n0)

        envname = paste0('fits/',region,'_m',m,'/',method,'_temp_m',mnth,'_l',loc,'.RDS')
        qf.y1.mle.ts <- readRDS(envname)
        envname = paste0('fits/',region,'_m',m,'/',method,'_prcp_m',mnth,'_l',loc,'.RDS')
        qf.y2.mle.ts <- readRDS(envname)
        qf.y2.mle.ts <- exp(qf.y2.mle.ts) - 0.0001

        y1y2.cors.1 = c(y1y2.cors.1,cor(y1[y0==1],y2[y0==1]))
        y1y2.cors.2 = c(y1y2.cors.2,cor(y1[y0==0],y2[y0==0]))
        y1y2.cors.0 = c(y1y2.cors.0,cor(qf.y1.mle.ts[y0==1],qf.y2.mle.ts[y0==1]))

        cal.array[,1,loc] = y1[y0==0]
        cal.array[,3,loc] = y1[y0==1]
        cal.array[,2,loc] = qf.y1.mle.ts[y0==0]
        cal.array[,4,loc] = y2[y0==0]
        cal.array[,6,loc] = y2[y0==1]
        cal.array[,5,loc] = qf.y2.mle.ts[y0==0]
        cal.array[,7,loc] = mnth
    }
    cal.data[[mnth]] = cal.array
}

metrics_all <- matrix(NA,10,4)
metrics_se <- matrix(NA,10,4)
eachmonth = rep(NA,12)

cal.array = do.call(abind::abind,c(cal.data,along=1))
cal.array2 = apply(cal.array, 2, c)
mnth.vec = cal.array2[,7]
season.vec <- rep(NA,length(mnth.vec))
season.vec[mnth.vec <= 3] = 1
season.vec[mnth.vec >=4 & mnth.vec <=6] = 2
season.vec[mnth.vec >=7 & mnth.vec <=9] = 3
season.vec[mnth.vec >=10 & mnth.vec <=12] = 4
season.abb <- c('JFM','AMJ','JAS','OND')
summary(mnth.vec)

tmax.mnth.mean <- matrix(NA,12,3)
prcp.mnth.mean <- matrix(NA,12,3)
prcp.mnth.zero <- matrix(NA,12,3)
prcp.mnth.95 <- matrix(NA,12,3)
tmax.mnth.95 <- matrix(NA,12,3)

for(mnth in 1:12){
    tmax.mnth.mean[mnth,] = apply(cal.array2[mnth.vec==mnth,1:3],2,mean)
    tmax.mnth.95[mnth,] = apply(cal.array2[mnth.vec==mnth,1:3],2,quantile,0.95)
    cal.array3 = round(cal.array2[,4:6],4)
    prcp.mnth.mean[mnth,] = apply(cal.array3[mnth.vec==mnth,],2,mean)
    prcp.mnth.95[mnth,] = apply(cal.array3[mnth.vec==mnth,],2,quantile,0.95)
    cal.array3 = ifelse(cal.array3>0,0,1)
    prcp.mnth.zero[mnth,] <- apply(cal.array3[mnth.vec==mnth,],2,mean)
}

pdf(paste0('plots/means_',region,'_monthly.pdf'),width = 9,height = 4)
par(mfrow=c(1,2),mgp=c(2.25,0.75,0),mar=c(4,4,1,1))
plot(1:12,tmax.mnth.mean[,1],'b',ylim = range(tmax.mnth.mean),col='red',
     ylab = 'mean TMAX (K)',xlab = 'Month',xaxt = "n",pch=20,main=paste(region,'TMAX mean'))
axis(1,at=1:12,labels=month.abb)
lines(1:12,tmax.mnth.mean[,3],'b',col='black',pch=20)
lines(1:12,tmax.mnth.mean[,2],'b',col='black',lty=2,pch=4)
legend('topleft',c('Mod','Obs','Cal'),col=c(2,1,1),lty = c(1,1,2),lwd=2,pch=c(20,20,4))

plot(1:12,prcp.mnth.mean[,1],'b',ylim = range(prcp.mnth.mean),col='red',
     ylab = 'mean PRCP (mm)',xlab = 'Month',xaxt = "n",pch=20,main=paste(region,'PRCP mean'))
axis(1,at=1:12,labels=month.abb)
lines(1:12,prcp.mnth.mean[,3],'b',col='black',pch=20)
lines(1:12,prcp.mnth.mean[,2],'b',col='black',lty=2,pch=4)
legend('topleft',c('Mod','Obs','Cal'),col=c(2,1,1),lty = c(1,1,2),lwd=2,pch=c(20,20,4))
dev.off()

pdf(paste0('plots/quantiles_',region,'_monthly.pdf'),width = 13,height = 4)
par(mfrow=c(1,3),mgp=c(2.25,0.75,0),mar=c(4,4,1,1))
plot(1:12,tmax.mnth.95[,1],'b',ylim = range(tmax.mnth.95),col='red',
     ylab = '0.95 quantile TMAX (K)',xlab = 'Month',xaxt = "n",pch=20,main=paste(region,'TMAX quantile'))
axis(1,at=1:12,labels=month.abb)
lines(1:12,tmax.mnth.95[,3],'b',col='black',pch=20)
lines(1:12,tmax.mnth.95[,2],'b',col='black',lty=2,pch=4)
legend('topleft',c('Mod','Obs','Cal'),col=c(2,1,1),lty = c(1,1,2),lwd=2,pch=c(20,20,4))

plot(1:12,prcp.mnth.zero[,1],'b',ylim = range(prcp.mnth.zero),col='red',
     ylab = 'Proportion of zero PRCP',xlab = 'Month',xaxt = "n",pch=20,main=paste(region,'PRCP zeros'))
axis(1,at=1:12,labels=month.abb)
lines(1:12,prcp.mnth.zero[,3],'b',col='black',pch=20)
lines(1:12,prcp.mnth.zero[,2],'b',col='black',lty=2,pch=4)
legend('topleft',c('Mod','Obs','Cal'),col=c(2,1,1),lty = c(1,1,2),lwd=2,pch=c(20,20,4))

plot(1:12,prcp.mnth.95[,1],'b',ylim = range(prcp.mnth.95),col='red',
     ylab = '0.95 quantile PRCP (mm)',xlab = 'Month',xaxt = "n",pch=20,main=paste(region,'PRCP quantile'))
axis(1,at=1:12,labels=month.abb)
lines(1:12,prcp.mnth.95[,3],'b',col='black',pch=20)
lines(1:12,prcp.mnth.95[,2],'b',col='black',lty=2,pch=4)
legend('topleft',c('Mod','Obs','Cal'),col=c(2,1,1),lty = c(1,1,2),lwd=2,pch=c(20,20,4))
dev.off()
# 
# save(y1y2.cors.0,y1y2.cors.1,y1y2.cors.2,cal.data,
#            file = paste0('results/summary_',method,'_',region,'_m',m,'_SPQR.RData'))
# load(paste0('results/summary_',method,'_',region,'_m',m,'_SPQR.RData'))

pdf(paste0('plots/density_',region,'_seasonal.pdf'),width = 10,height = 4)
par(mfrow=c(2,4),mgp=c(2.25,0.75,0),mar=c(4,4,1,1))
for(ssn in 1:4){
    d0 <-density(cal.array2[season.vec==ssn,1]) # gcm
    d1 <-density(cal.array2[season.vec==ssn,3]) # obs 
    d2 <- density(cal.array2[season.vec==ssn,2]) # pred
    plotmax.y = max(d0$y,d1$y,d2$y)
    plotmin.y = min(d0$y,d1$y,d2$y)
    plotmax.x = max(d0$x,d1$x,d2$x)
    plotmin.x = min(d0$x,d1$x,d2$x)
    plot(d0,col=2,ylim=range(c(plotmin.y,plotmax.y)),
         xlim=range(c(plotmin.x,plotmax.x)),ylab="Density",xlab='TMAX',main=paste(season.abb[ssn],region,'TMAX'))
    lines(d1,col=1)
    lines(d2,col=1,lty=2)
    legend('topleft',c('Mod','Obs','Cal'),col=c(2,1,1),lty = c(1,1,2),lwd=2)
    d0 <-density(log(0.0001+cal.array2[season.vec==ssn,4])) # gcm
    d1 <-density(log(0.0001+cal.array2[season.vec==ssn,6])) # obs 
    d2 <- density(log(0.0001+cal.array2[season.vec==ssn,5])) # pred
    plotmax.y = max(d0$y,d1$y,d2$y)
    plotmin.y = min(d0$y,d1$y,d2$y)
    plotmax.x = max(d0$x,d1$x,d2$x)
    plotmin.x = min(d0$x,d1$x,d2$x)
    plot(d0,col=2,ylim=range(c(plotmin.y,plotmax.y)),
         xlim=range(c(plotmin.x,plotmax.x)),ylab="Density",xlab = 'PRCP',main=paste(season.abb[ssn],region,'PRCP'))
    lines(d1,col=1)
    lines(d2,col=1,lty=2)
    legend('topright',c('Mod','Obs','Cal'),col=c(2,1,1),lty = c(1,1,2),lwd=2)
}
par(mfrow=c(1,1)) 
dev.off()