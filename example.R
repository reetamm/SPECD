# Conditional quantile function for precip
q1 <- function(u1,y0){  
    a <- 3+0.5*y0
    b <- 4-1.0*y0
    y <- qgamma(u1,a,b)
    return(y)}

# Conditional inverse quantile function for precip
p1 <- function(y1,y0){
    a <- 3+0.5*y0
    b <- 4-1.0*y0
    u <- pgamma(y1,a,b)
    return(u)}

# Conditional quantile function for temperature
q2 <- function(u2,y0,y1){
    m <- 1.0-0.5*y1 + 0.5*y0*log(y1)
    s <- 0.2*y1*(1-y0)+.2
    y <- qnorm(u2,m,s)
    return(y)}

# Conditional inverse quantile function for temperature
p2 <- function(y2,y0,y1){
    m <- 1.0-0.5*y1 + 0.5*y0*log(y1)
    s <- 0.2*y1*(1-y0)+.2
    u <- pnorm(y2,m,s)
    return(u)}

set.seed(919)
n  <- 100000
y0 <- rbinom(n,1,.8)
u1 <- runif(n)
y1 <- q1(u1,y0)
u2 <- runif(n)
y2 <- q2(u2,y0,y1)
X <- cbind(u1,y0)
library(SPQR)
control <- list(iter = 300, warmup = 200, thin = 1)
fit <- SPQR(X = y0, Y = y1, method = "MLE", control = control, normalize = TRUE, verbose = T)
head(X)

qout <- u1[y0==0]
head(x_test)
x_test[,2] = 1

qf.mle  <- predict(fit, X = 1, type = "QF",tau=qout)
qf.mle  <- predict(fit, X = 1, type = "QF",tau=runif(1000))
# hist(qf.mle)
# hist(y1[y0==1])

summary(qf.mle)
summary(y1[y0==1])
summary(y1[y0==0])

d0 <-density(y1[y0==0]) 
d1 <-density(y1[y0==1]) 
d2 <- density(qf.mle)
plot(d0,col=1,ylim=range(c(d0$y,d1$y)),ylab="Y1",main="Y1")
lines(d1,col=2)
lines(d2,col=3)
legend("topright",c("Model","Observations",'Bias-corrected'),col=1:3,lwd=2,bty="n")


