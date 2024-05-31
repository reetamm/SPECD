rm(list = ls())
library(ggplot2)
library(scales)
library(lubridate)
library(GpGp)
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
y1y2.cors.0 = NA
y1y2.cors.1 = NA
daysinmonth = c(31,28,31,30,31,30,31,31,30,31,30,31)

cal.data = vector('list',12)
for(mnth in 1:12){
    cal.array = array(dim = c(daysinmonth[mnth]*64,6,25))
    for(loc in 1:25){
        print(paste(mnth,loc))
        y1 <- c(obs.long$tmax[vecchia.order==loc & gcm.months==mnth],gcm.long$tmax[vecchia.order==loc & gcm.months==mnth])
        y2 <- c(obs.long$pr[vecchia.order==loc & gcm.months==mnth],gcm.long$pr[vecchia.order==loc & gcm.months==mnth])
        y2 <- log(1+y2)
        n0 = length(gcm.long$pr[vecchia.order==loc & gcm.months==mnth]); n1 = length(obs.long$pr[vecchia.order==loc & gcm.months==mnth])
        n = n0 + n1
        y0 <- rep(1:0,each=n0)
        
        envname = paste0('fits/fits_temp_m',mnth,'_l',loc,'.RDS')
        qf.y1.mle.ts <- readRDS(envname)
        envname = paste0('fits/fits_prcp_m',mnth,'_l',loc,'.RDS')
        qf.y2.mle.ts <- readRDS(envname)
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

save(y1y2.cors.0,y1y2.cors.1,cal.data,
           file = 'summary_space.RData')
load('summary_space.RData')
eachmonth = rep(NA,12)
# for(i in 1:12){
#     eachmonth[i] = dim(cal.data[[i]])[1]
#     cal.data[[i]][,4:6,] = exp(cal.data[[i]][,4:6,])-1
# }
cal.array = cal.data[[6]]
mnth=12
loc=1
maxy = 1.01*max(cal.array[1:31,4:6,loc])
miny = 0.99*min(cal.array[1:31,4:6,loc])
plot(1:31,cal.array[1:31,4,loc],type = 'b',col=1,pch=20,ylim = c(miny,maxy))
lines(1:31,cal.array[1:31,5,loc],type = 'b',col=2,pch=20)
lines(1:31,cal.array[1:31,6,loc],type = 'b',col=3,pch=20)



plot(y1y2.cors.1,y1y2.cors.0)
abline(0,1)
q1 = 0.9
q2 = 0.9

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
tmp = apply(prcp.tail,3,c)
dim(tmp)
xmax = max(tmp)
xmin = min(tmp)
plot(density(tmp[,2]),xlab = 'probability',main = 'precip',xlim = c(xmin,xmax))
lines(density(tmp[,1]),col=2)
abline(v=1-q2)
legend('topright',c('Calibrated','GCM'),col=1:2,lwd=2,bty="n")



tmp = apply(temp.tail,3,c)
dim(tmp)
xmax = max(tmp)
xmin = min(tmp)
plot(density(tmp[,1]),xlab = 'probability',main = 'temp',xlim = c(xmin,xmax),col=2)
lines(density(tmp[,2]))
abline(v=1-q1)
legend('topright',c('Calibrated','GCM'),col=1:2,lwd=2,bty="n")

tmp = apply(joint.tail,3,c)
dim(tmp)
xmax = max(tmp)
xmin = min(tmp)
plot(density(tmp[,2]),xlab = 'probability',main = 'temp',xlim = c(xmin,xmax),col=2)
lines(density(tmp[,1]))
lines(density(tmp[,3]))
abline(v=1-q1)
legend('topright',c('Calibrated','GCM'),col=1:2,lwd=2,bty="n")

plot(tmp[,3],tmp[,2])
plot(tmp[,1],tmp[,2])
abline(0,1)

coords = cbind(coords,vecchia.order)
names(coords)
ggplot(coords,aes(x=lon,y=lat,col='white'))+geom_tile()+geom_text(aes(label=vecchia.order))

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
           correls[count,k] = cor(cal.array[,k,c(1,j)])[1,2]
}
}
par(mfrow=c(2,3))
plot(correls[,7],correls[,1],col=correls[,8],xlim = c(0,8),pch=20,
     xlab = 'Distance',ylab = 'Correlation',main = 'GCM temp')
legend('topright',month.abb,col=1:12,lwd=2,bty="n")

plot(correls[,7],correls[,2],col=correls[,8],xlim = c(0,8),pch=20,
     xlab = 'Distance',ylab = 'Correlation',main = 'Calibrated temp')
legend('topright',month.abb,col=1:12,lwd=2,bty="n")

plot(correls[,7],correls[,3],col=correls[,8],xlim = c(0,8),pch=20,
     xlab = 'Distance',ylab = 'Correlation',main = 'Observed temp')
legend('topright',month.abb,col=1:12,lwd=2,bty="n")


plot(correls[,7],correls[,4],col=correls[,8],xlim = c(0,8),pch=20,
     xlab = 'Distance',ylab = 'Correlation',main = 'GCM prcp')
legend('topright',month.abb,col=1:12,lwd=2,bty="n")

plot(correls[,7],correls[,5],col=correls[,8],xlim = c(0,8),pch=20,
     xlab = 'Distance',ylab = 'Correlation',main = 'Calibrated prcp')
legend('topright',month.abb,col=1:12,lwd=2,bty="n")

plot(correls[,7],correls[,6],col=correls[,8],xlim = c(0,8),pch=20,
     xlab = 'Distance',ylab = 'Correlation',main = 'Observed prcp')
legend('topright',month.abb,col=1:12,lwd=2,bty="n")
par(mfrow=c(1,1))

plot(correls[,2],correls[,3],col=correls[,7]+1,pch=20,
     xlab = 'Calibrated',ylab = 'Observed',main = 'Temp')
abline(0,1)
points(correls[,2],correls[,1],col=alpha(1,0.05),pch=20,
     xlab = 'Calibrated',ylab = 'Observed',main = 'Temp')
legend('topleft','GCM',col=1,lwd=2,bty="n")

plot(correls[,5],correls[,6],col=correls[,7]+1,pch=20,
     xlab = 'Calibrated',ylab = 'Observed',main = 'Prcp')
abline(0,1)
points(correls[,5],correls[,4],col=alpha(1,0.05),pch=20,
       xlab = 'Calibrated',ylab = 'Observed',main = 'Prcp')
legend('topleft','GCM',col=1,lwd=2,bty="n")

# creating correlation matrix
corr_mat1 <- gdata::upperTriangle(cor(cal.array[,2,]))
corr_mat2 <- gdata::upperTriangle(cor(cal.array[,3,]))
plot(corr_mat1,corr_mat2)
abline(0,1)

corr_mat1 <- gdata::upperTriangle(cor(cal.array[,5,]))
corr_mat2 <- gdata::upperTriangle(cor(cal.array[,6,]))
plot(corr_mat1,corr_mat2)
abline(0,1)
