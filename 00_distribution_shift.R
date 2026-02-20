rm(list = ls())
library(GpGp)
library(SPQR)
library(lubridate)
library(tidyr)
region = 'SW'
method = 'MLE'
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
gcm.long$grid.no <- grid.no
obs.long$grid.no <- grid.no
mnth = 1

pdf('plots/dist_shift_SW_TMAX.pdf',width = 9,height = 12)
# temperature correlations GCM
par(mfrow=c(4,3),mgp=c(2.25,0.75,0),mar=c(4,4,1,1))
for(mnth in 1:12){
X10 = gcm.long[gcm.months==mnth & gcm.years <= 2000,c(1,5,6)]
head(X10)
X10_long <- pivot_wider(X10,names_from = grid.no,values_from = tmax)[,-1]
X10_cor = gdata::upperTriangle(cor(X10_long))
length(X10_cor)

X11 = gcm.long[gcm.months==mnth & gcm.years > 2000,c(1,5,6)]
head(X11)
X11_long <- pivot_wider(X11,names_from = grid.no,values_from = tmax)[,-1]
X11_cor = gdata::upperTriangle(cor(X11_long))
length(X11_cor)

plot(X10_cor,X11_cor,xlim=c(0,1),ylim=c(0,1), main = paste(month.abb[mnth],region,'TMAX'),pch=20,cex=0.75,
     xlab = '1951-2000', ylab = '2001-2014',col='red')


X10 = obs.long[gcm.months==mnth & gcm.years <= 2000,c(1,4,6)]
head(X10)
X10_long <- pivot_wider(X10,names_from = grid.no,values_from = tmax)[,-1]
X10_cor = gdata::upperTriangle(cor(X10_long))
length(X10_cor)

X11 = obs.long[gcm.months==mnth & gcm.years > 2000,c(1,4,6)]
head(X11)
X11_long <- pivot_wider(X11,names_from = grid.no,values_from = tmax)[,-1]
X11_cor = gdata::upperTriangle(cor(X11_long))
length(X11_cor)

points(X10_cor,X11_cor,xlim=c(0,1),ylim=c(0,1), main = paste(month.abb[mnth],region,'Obs tmax'),pch=20,cex=0.75,
     xlab = '1951-2000', ylab = '2001-2014')

abline(0,1)
legend('topleft',c('Mod','Obs'),col=c(2,1),lty = c(1,1),lwd=2)
}
par(mfrow=c(1,1))
dev.off()

pdf('plots/dist_shift_SW_PRCP.pdf',width = 9,height = 12)
# precip correlations GCM
par(mfrow=c(4,3),mgp=c(2.25,0.75,0),mar=c(4,4,1,1))
for(mnth in 1:12){
    X10 = gcm.long[gcm.months==mnth & gcm.years <= 2000,c(1,4,6)]
    head(X10)
    X10_long <- pivot_wider(X10,names_from = grid.no,values_from = pr)[,-1]
    X10_cor = gdata::upperTriangle(cor(X10_long))
    length(X10_cor)
    
    X11 = gcm.long[gcm.months==mnth & gcm.years > 2000,c(1,4,6)]
    head(X11)
    X11_long <- pivot_wider(X11,names_from = grid.no,values_from = pr)[,-1]
    X11_cor = gdata::upperTriangle(cor(X11_long))
    length(X11_cor)
    
    plot(X10_cor,X11_cor,xlim=c(0,1),ylim=c(0,1), main = paste(month.abb[mnth],region,'PRCP'),pch=20,cex=0.75,
         xlab = '1951-2000', ylab = '2001-2014',col='red')
   
    X10 = obs.long[gcm.months==mnth & gcm.years <= 2000,c(1,5,6)]
    head(X10)
    X10_long <- pivot_wider(X10,names_from = grid.no,values_from = pr)[,-1]
    X10_cor = gdata::upperTriangle(cor(X10_long))
    length(X10_cor)
    
    X11 = obs.long[gcm.months==mnth & gcm.years > 2000,c(1,5,6)]
    head(X11)
    X11_long <- pivot_wider(X11,names_from = grid.no,values_from = pr)[,-1]
    X11_cor = gdata::upperTriangle(cor(X11_long))
    length(X11_cor)
    
    points(X10_cor,X11_cor,xlim=c(0,1),ylim=c(0,1), main = paste(month.abb[mnth],region,'Obs prcp'),pch=20,cex=0.75,
         xlab = '1951-2000', ylab = '2001-2014')
    abline(0,1)
    legend('topleft',c('Mod','Obs'),col=c(2,1),lty = c(1,1),lwd=2)
}
par(mfrow=c(1,1))
dev.off()

