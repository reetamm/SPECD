# Conditional quantile function for precip
rm(list=ls())
library(SPQR)
library(mvtnorm)
q1 <- function(u1,y0){  
    n = length(u)
    
    a <- 3+0.5*y0
    b <- 4-1.0*y0
    y <- qgamma(u1,a,b)
    return(y)}

# Conditional quantile function for temperature
q2 <- function(u2,y0,y1){
    m <- 1.0-0.5*y1 + 0.5*y0*log(y1)
    s <- 0.2*y1*(1-y0)+.2
    y <- qnorm(u2,m,s)
    return(y)}

# data generation
set.seed(919)
n <- 10000                  #define length of simulation
n0 = n*0.4
n1 = n*0.6
Sigma = matrix(c(1,0.5,0.5,1),nrow = 2)
yy = pnorm(rmvnorm(n0,sigma = Sigma))
cor(yy)
aa = matrix(0,nrow = n0,ncol = 2)

aa[1,] = yy[1,]
for(i in 2:n0)
    aa[i,] = c(0.3,0.5)*aa[i-1,] + c(0.7,0.5)*yy[i,]

cor(aa[-1,1],aa[-n0,1])
cor(aa[-1,2],aa[-n0,2])
cor(aa)

Sigma = matrix(c(1,0.6,0.6,1),nrow = 2)
yy = pnorm(rmvnorm(n1,sigma = Sigma))
cor(yy)
bb = matrix(0,nrow = n1,ncol = 2)

bb[1,] = bb[1,]
for(i in 2:n1)
    bb[i,] = c(0.6,0.8)*bb[i-1,] + c(0.4,0.2)*yy[i,]

cor(bb[-1,1],bb[-n1,1])
cor(bb[-1,2],bb[-n1,2])
cor(bb)


y0 <- c(rep(1,n*.6),rep(0,n*.4))
u1 = c(aa[,2],bb[,2])
u2 = c(aa[,1],bb[,1])
summary(u1)
summary(u2)

y1 <- q1(u1,y0)
y2 <- q2(u2,y0,y1)
plot(y1,y2,col=y0+1)

y11 = y1[y0==1]
x11 = y11[c(n1,1:(n1-1))]
cor(cbind(y11,x11))
y10 = y1[y0==0]
x10 = y10[c(n0,1:(n0-1))]
cor(cbind(y10,x10))
y1 = c(y11,y10)
x1 = c(x11,x10)
X1 = cbind(x1,y0)


control <- list(iter = 300)
fit.y1.map.ts <- SPQR(X = X1, Y = y1, method = "MAP", control = control, normalize = T, verbose = T,use.GPU=T,
                      n.hidden = c(30,20), activation = 'relu')
fit.y1.map <- SPQR(X = X1[,2], Y = y1, method = "MAP", control = control, normalize = T, verbose = T,use.GPU=T,
                   n.hidden = c(30,20), activation = 'relu')
fit.y1.mle.ts <- SPQR(X = X1, Y = y1, method = "MLE", control = control, normalize = T, verbose = T,use.GPU=T,
                      n.hidden = c(30,20), activation = 'relu')
fit.y1.mle <- SPQR(X = X1[,2], Y = y1, method = "MLE", control = control, normalize = T, verbose = T,use.GPU=T,
                   n.hidden = c(30,20), activation = 'relu')

cdf.y1.mle.ts = rep(NA,n)
cdf.y1.map.ts = rep(NA,n)
cdf.y1.mle = rep(NA,n)
cdf.y1.map = rep(NA,n)
for(i in 1:n){
    cdf.y1.mle.ts[i] <- predict(fit.y1.mle.ts,   X = X1[i,], Y=y1[i], type = "CDF")    
    cdf.y1.mle[i] <- predict(fit.y1.mle,   X = X1[i,2], Y=y1[i], type = "CDF")    
    cdf.y1.map.ts[i] <- predict(fit.y1.map.ts,   X = X1[i,], Y=y1[i], type = "CDF")    
    cdf.y1.map[i] <- predict(fit.y1.map,   X = X1[i,2], Y=y1[i], type = "CDF")    
    print(i)
}

qout11 <- cdf.y1.mle.ts
qout12 <- cdf.y1.mle
qout13 <- cdf.y1.map.ts
qout14 <- cdf.y1.map
adjust = which(qout11>0.99999 & qout12>0.99999 & qout13>0.99999 & qout14>0.99999)
qout11[adjust] = 0.99999
qout12[adjust] = 0.99999
qout13[adjust] = 0.99999
qout14[adjust] = 0.99999

qf.y1.mle.ts = rep(NA,n)
qf.y1.mle = rep(NA,n)
qf.y1.map.ts = rep(NA,n)
qf.y1.map = rep(NA,n)
x_pred = c(X1[1,1],1)
qf.y1.mle.ts[1] <- predict(fit.y1.mle.ts,   X = x_pred, type = "QF",tau=qout11[1])    
qf.y1.map.ts[1] <- predict(fit.y1.map.ts,   X = x_pred, type = "QF",tau=qout13[1])    


for(i in 2:n){
    print(i)
    if(i==12001)
        x_pred = c(X1[i,1],1)
    if(i!=12001)
        x_pred = c(qf.y1.mle.ts[i-1],1)
        qf.y1.mle.ts[i] <- predict(fit.y1.mle.ts,   X = x_pred, type = "QF",tau=qout11[i])
        qf.y1.map.ts[i] <- predict(fit.y1.map.ts,   X = x_pred, type = "QF",tau=qout13[i])
}
qf.y1.mle <- predict(fit.y1.mle,   X = 1, type = "QF",tau=qout12)    
qf.y1.map <- predict(fit.y1.map,   X = 1, type = "QF",tau=qout14)    
par(mfrow= c(2,2))
plot(y1,qf.y1.mle.ts,col=y0+1, main = 'MLE-TS')
abline(0,1)
plot(y1,qf.y1.mle,col=y0+1, main = 'MLE')
abline(0,1)
plot(y1,qf.y1.map.ts,col=y0+1, main = 'MAP-TS')
abline(0,1)
plot(y1,qf.y1.map,col=y0+1, main = 'MAP')
abline(0,1)
par(mfrow= c(1,1))


summary(qf.y1.mle.ts[y0==0])
summary(qf.y1.mle[y0==0])
summary(qf.y1.mle.ts[y0==1])
summary(qf.y1.mle[y0==1])
cor(qf.y1.mle[y0==0][-1],qf.y1.mle[y0==0][-8000])
cor(qf.y1.mle.ts[y0==0][-1],qf.y1.mle.ts[y0==0][-8000])
cor(qf.y1.mle[y0==1][-1],qf.y1.mle[y0==1][-12000])
cor(qf.y1.mle.ts[y0==1][-1],qf.y1.mle.ts[y0==1][-12000])

summary(y1[y0==1])
summary(y1[y0==0])

d0 <-density(y1[y0==0]) 
d1 <-density(y1[y0==1]) 
d2 <- density(qf.y1.map.ts[y0==0])
d3 <- density(qf.y1.map[y0==0])
d4 <- density(qf.y1.mle.ts[y0==0])
d5 <- density(qf.y1.mle[y0==0])
plot(d0,col=1,ylim=range(c(d0$y,d1$y)),ylab="Y1",main="Y1")
lines(d1,col=2)
lines(d2,col=3,lty=2)
lines(d3,col=4,lty=3)
lines(d4,col=5,lty=2)
lines(d5,col=6,lty=3)

legend("topright",c("Model","Observations",'MAP-TS','MAP','MLE-TS','MLE'),
       col=1:6,lwd=2,bty="n",lty=c(1,1,2,3,2,3))


###################################
###################################
y21 = y2[y0==1]
x21 = y21[c(n1,1:(n1-1))]
cor(cbind(y21,x21))
y20 = y2[y0==0]
x20 = y20[c(n0,1:(n0-1))]
cor(cbind(y20,x20))
y2 = c(y21,y20)
x2 = c(x21,x20)
X2 = cbind(x2,y1,X1)
head(X2)
head(y2)
plot(X2[,3],y2,col=y0+1)
fit.y2.mle.ts <- SPQR(X = X2, Y = y2, method = "MLE", control = control, normalize = T, verbose = T,use.GPU=T,
                      n.hidden = c(30,20), activation = 'relu')
fit.y2.mle <- SPQR(X = X2[,c(2,4)], Y = y2, method = "MLE", control = control, normalize = T, verbose = T,use.GPU=T,
                   n.hidden = c(30,20), activation = 'relu')

cdf.y2.mle.ts = rep(NA,n)
cdf.y2.mle = rep(NA,n)
for(i in 1:n){
    cdf.y2.mle.ts[i] <- predict(fit.y2.mle.ts,   X = X2[i,], Y=y2[i], type = "CDF")    
    cdf.y2.mle[i] <- predict(fit.y2.mle,   X = X2[i,c(2,4)], Y=y2[i], type = "CDF")    
    print(i)
}
qout21 <- cdf.y2.mle.ts
qout22 <- cdf.y2.mle
adjust = which(qout21>0.99999)
qout21[adjust] = 0.99999
adjust = which(qout22>0.99999)
qout22[adjust] = 0.99999

qf.y2.mle.ts = rep(NA,n)
qf.y2.mle = rep(NA,n)
x_pred = c(X2[1,1],qf.y1.mle.ts[1],X2[1,3],1)
qf.y2.mle.ts[1] <- predict(fit.y2.mle.ts,   X = x_pred,             type = "QF",tau=qout21[1])    
qf.y2.mle[1]    <- predict(fit.y2.mle,      X = cbind(X2[1,2],1),   type = "QF",tau=qout22[1])
for(i in 2:n){
    print(i)
    if(i!=12001)
        x_pred = c(qf.y2.mle.ts[i-1],qf.y1.mle.ts[i],qf.y1.mle.ts[i-1],1)
    if(i==12001)
        x_pred = c(X2[i,1],qf.y1.mle.ts[i],qf.y1.mle.ts[i-1],1)
    
    qf.y2.mle.ts[i] <- predict(fit.y2.mle.ts,   X = x_pred,             type = "QF",tau=qout21[i])
    qf.y2.mle[i]    <- predict(fit.y2.mle,      X = cbind(X2[i,2],1),   type = "QF",tau=qout22[i])
    
}
plot(y2,qf.y2.mle.ts,col=y0+1)
abline(0,1)

qout2 <- cdf.y2.mle#[80001:100000]
adjust = which(qout2>0.99999)
qout2[adjust] = 0.99999
max(qout2)
plot(y2,qf.y2.mle,col=y0+1)
abline(0,1)

summary(qf.y2.mle.ts[y0==0])
summary(qf.y2.mle[y0==0])
summary(qf.y2.mle.ts[y0==1])
summary(qf.y2.mle[y0==1])
cor(qf.y2.mle[y0==0][-1],qf.y2.mle[y0==0][-8000])
cor(qf.y2.mle.ts[y0==0][-1],qf.y2.mle.ts[y0==0][-8000])
cor(qf.y2.mle[y0==1][-1],qf.y2.mle[y0==1][-12000])
cor(qf.y2.mle.ts[y0==1][-1],qf.y2.mle.ts[y0==1][-12000])

summary(y2[y0==1])
summary(y2[y0==0])

d0 <-density(y2[y0==0]) 
d1 <-density(y2[y0==1]) 
# d2 <- density(qf.y1.map.ts)
# d3 <- density(qf.y1.map)
d4 <- density(qf.y2.mle.ts[y0==0])
d5 <- density(qf.y2.mle[y0==0])
plot(d0,col=1,ylim=range(c(d0$y,d1$y)),ylab="Y2",main="Y2")
lines(d1,col=2)
# lines(d2,col=3,lty=2)
# lines(d3,col=4,lty=3)
lines(d4,col=5,lty=2)
lines(d5,col=6,lty=3)

legend("topright",c("Model","Observations",'MAP-TS','MAP','MLE-TS','MLE'),
       col=1:6,lwd=2,bty="n",lty=c(1,1,2,3,2,3))

