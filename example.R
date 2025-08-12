# Conditional quantile function for precip
rm(list=ls())
library(SPQR)

q1 <- function(u1,y0){  
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

set.seed(919)
n  <- 20000
n0 = n*0.4
n1 = n*0.6
aa1 = arima.sim(n = n1, list(ar = c(0.6)),
               sd = 1)
aa2 = arima.sim(n = n0, list(ar = c(0.3)),
                sd = 1)
u1 = pnorm(c(aa1,aa2))
aa1 = arima.sim(n = n1, list(ar = c(0.75)),
                sd = 1)
aa2 = arima.sim(n = n0, list(ar = c(0.5)),
                sd = 1)
u2 = pnorm(c(aa1,aa2))


y0 <- c(rep(1,n*.6),rep(0,n*.4))
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
control <- list(iter = 300, batch.size = 100, lr = 0.001)
fit.y1.mle.ts <- SPQR(X = X1, Y = y1, method = "MLE", control = control, normalize = T, verbose = T,use.GPU=T,
                      n.hidden = c(30,20), activation = 'relu',n.knots = 15)
fit.y1.mle <- SPQR(X = X1[,2], Y = y1, method = "MLE", control = control, normalize = T, verbose = T,use.GPU=T,
                   n.hidden = c(30,20), activation = 'relu',n.knots = 15)

cdf.y1.mle.ts = rep(NA,n)
cdf.y1.mle = rep(NA,n)
for(i in 1:n){
    cdf.y1.mle.ts[i] <- predict(fit.y1.mle.ts,   X = X1[i,], Y=y1[i], type = "CDF")    
    cdf.y1.mle[i] <- predict(fit.y1.mle,   X = X1[i,2], Y=y1[i], type = "CDF")    
    print(i)
}

qout11 <- cdf.y1.mle.ts
qout12 <- cdf.y1.mle
adjust = which(qout11>0.99999 & qout12>0.99999)
qout11[adjust] = 0.99999
qout12[adjust] = 0.99999

qf.y1.mle.ts = rep(NA,n)
qf.y1.mle = rep(NA,n)
x_pred = c(X1[1,1],1)
qf.y1.mle.ts[1] <- predict(fit.y1.mle.ts,   X = x_pred, type = "QF",tau=qout11[1])    

for(i in 2:n){
    print(i)
    if(i==12001)
        x_pred = c(X1[i,1],1)
    if(i!=12001)
        x_pred = c(qf.y1.mle.ts[i-1],1)
        qf.y1.mle.ts[i] <- predict(fit.y1.mle.ts,   X = x_pred, type = "QF",tau=qout11[i])
}
qf.y1.mle <- predict(fit.y1.mle,   X = 1, type = "QF",tau=qout12)    
par(mfrow= c(1,2))
plot(y1,qf.y1.mle.ts,col=y0+1, main = 'MLE-TS')
abline(0,1)
plot(y1,qf.y1.mle,col=y0+1, main = 'MLE')
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
d2 <- density(qf.y1.mle.ts[y0==0])
d3 <- density(qf.y1.mle[y0==0])
plot(d0,col=1,ylim=range(c(d0$y,d1$y)),ylab="Y1",main="Y1")
lines(d1,col=2)
lines(d2,col=3,lty=2)
lines(d3,col=4,lty=3)

legend("topright",c("Model","Observations",'MLE-TS','MLE'),
       col=1:6,lwd=2,bty="n",lty=c(1,1,2,3))


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
X2 = cbind(x2,y1,y0)
head(X2)
head(y2)
fit.y2.mle.ts <- SPQR(X = X2, Y = y2, method = "MLE", control = control, normalize = T, verbose = T,use.GPU=T,
                      n.hidden = c(20,10), activation = 'relu',n.knots = 20)
fit.y2.mle <- SPQR(X = X2[,c(2,3)], Y = y2, method = "MLE", control = control, normalize = T, verbose = T,use.GPU=T,
                   n.hidden = c(20,10), activation = 'relu', n.knots = 20)
plotGOF(fit.y2.mle)
plotGOF(fit.y2.mle.ts)
cdf.y2.mle.ts = rep(NA,n)
cdf.y2.mle = rep(NA,n)
for(i in 1:n){
    cdf.y2.mle.ts[i] <- predict(fit.y2.mle.ts,   X = X2[i,], Y=y2[i], type = "CDF")    
    cdf.y2.mle[i] <- predict(fit.y2.mle,   X = X2[i,c(2,3)], Y=y2[i], type = "CDF")    
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
x_pred = c(X2[1,1],qf.y1.mle.ts[1],1)
qf.y2.mle.ts[1] <- predict(fit.y2.mle.ts,   X = x_pred, type = "QF",tau=qout21[1])    
qf.y2.mle[1]    <- predict(fit.y2.mle,      X = c(qf.y1.mle[1],1), type = "QF",tau=qout22[1])
for(i in 2:n){
    print(i)
    if(i!=12001)
        x_pred = c(qf.y2.mle.ts[i-1],qf.y1.mle.ts[i],1)
    if(i==12001)
        x_pred = c(X2[i,1],qf.y1.mle.ts[i],1)
    qf.y2.mle.ts[i] <- predict(fit.y2.mle.ts,   X = x_pred, type = "QF",tau=qout21[i])
    qf.y2.mle[i]    <- predict(fit.y2.mle,      X = c(qf.y1.mle[i],1), type = "QF",tau=qout22[i])
    
}
par(mfrow= c(1,2))
plot(y2,qf.y2.mle.ts,col=y0+1, main = 'MLE-TS')
abline(0,1)
plot(y2,qf.y2.mle,col=y0+1, main = 'MLE')
abline(0,1)
par(mfrow= c(1,1))

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
d2 <- density(qf.y2.mle.ts[y0==0])
d3 <- density(qf.y2.mle[y0==0])
plot(d0,col=1,ylim=range(c(d0$y,d1$y)),ylab="Y2",main="Y2")
lines(d1,col=2)
lines(d2,col=2,lty=2)
lines(d3,col=3,lty=3)

legend("topright",c("Model","Observations",'MLE-TS','MLE'),
       col=1:6,lwd=2,bty="n",lty=c(1,1,2,3))

