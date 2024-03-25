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
