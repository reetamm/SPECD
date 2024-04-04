# Conditional quantile function for precip
rm(list = ls())
s1 = c(0.25,0.75)
s2 = c(0.25,0.75)

s = expand.grid(s1,s2)
set.seed(303)
n <- 4  
A <- matrix(runif(n^2,0)*2, ncol=n) 
Sigma <- (t(A) %*% A)
Sigma = cov2cor(Sigma)
Sigma = Sigma/2
diag(Sigma) = 1
Sigma

q1 <- function(s,y0,Sigma){  
    nlocs = nrow(s)
    n = length(y0)
    alldata = matrix(NA,n,4)
    mm <- ss <- alldata
    
    for(i in 1:nlocs){
        mm[,i] = 3*s[i,1]+0.5*y0    
        ss[,i] = 4-1.0*y0
    }
    
    alldata = rmvnorm(n,sigma = Sigma)
    alldata = alldata*ss + mm
    return(alldata)}


y0 = rbinom(100000,1,.8)

y1 = q1(s,y0,Sigma)

set.seed(919)
n  <- 100000
y0 <- rbinom(n,1,.5)
u1 <- runif(n)
y1 <- q1(u1,y0)
u2 <- runif(n)
y2 <- q2(u2,y0,y1)

library(SPQR)
control <- list(iter = 300, warmup = 200, thin = 1)
fit.temp <- SPQR(X = y0, Y = y1, method = "MLE", control = control, normalize = TRUE, verbose = T,use.GPU=T)

qout <- u1[y0==0]

qf.mle1  <- predict(fit.temp, X = 1, type = "QF",tau=qout)
qf.mle2  <- predict(fit.temp, X = 1, type = "QF",tau=runif(10000))
# hist(qf.mle)
# hist(y1[y0==1])

summary(qf.mle1)
summary(y1[y0==1])
summary(y1[y0==0])

d0 <-density(y1[y0==0]) 
d1 <-density(y1[y0==1]) 
d2 <- density(qf.mle1)
plot(d0,col=1,ylim=range(c(d0$y,d1$y)),ylab="Y1",main="Y1")
lines(d1,col=2)
lines(d2,col=3)
legend("topright",c("Model","Observations",'Bias-corrected'),col=1:3,lwd=2,bty="n")


X = cbind(y1,y0)
fit.prcp <- SPQR(X = X, Y = y2, method = "MLE", control = control, normalize = TRUE, verbose = T,use.GPU=T)
X_test <- cbind(as.numeric(qf.mle1),1)
head(X_test)
ntest = nrow(X_test)
u2_test = qout
qf.mle3 = rep(NA,ntest)
for(i in 1:ntest){
    print(i)
    qf.mle3[i] <- predict(fit.prcp,X=X_test[i,],tau = u2_test[i])
}

d0 <-density(y2[y0==0]) 
d1 <-density(y2[y0==1]) 
d2 <- density(qf.mle3)
plot(d0,col=1,ylim=range(c(d0$y,d1$y)),ylab="Y2",main="Y2")
lines(d1,col=2)
lines(d2,col=3)
legend("topright",c("Model","Observations",'Bias-corrected'),col=1:3,lwd=2,bty="n")
