library(splines2)

basis <- function(x,L,integral=FALSE){  
  Bknots <- 0:1
  Iknots <- seq(1/(L-3),1-1/(L-3),length=L-4)
  B      <- mSpline(x,knots=Iknots,Boundary.knots=Bknots,
                    intercept=TRUE,integral=integral)
return(B)}

act      <- function(x){pmax(x,0)}
actp     <- function(x){x>0}
Pen      <- function(x){sum(x^2)}
Pen_grad <- function(x){2*x}

SPQR_init <- function(p,L,K,init_mn=0,sigma=.01){
 w       <- list()
 w[[1]]  <- sigma*rnorm(K)
 w[[2]]  <- sigma*matrix(rnorm(p*K),p,K)
 w[[3]]  <- sigma*rnorm(L)
 w[[4]]  <- sigma*matrix(rnorm(L*K),K,L)
return(w)}

dSPQR <- function(y,x,w){
    B   <- basis(y,L=length(w[[3]]))
    p   <- SPQR_probs(w,x)
    out <- rowSums(B*p)
return(out)}

pSPQR <- function(y,x,w){
    B   <- basis(y,L=length(w[[3]]),integral=TRUE)
    p   <- SPQR_probs(w,x)
    out <- rowSums(B*p)
return(out)}

qSPQR <- function(tau,x,w){
   q  <- tau
   fx <- function(y,X,W,P){pSPQR(y,X,W)-P} 
   for(i in 1:length(q)){      
     q[i] <- uniroot(fx, interval=0:1,X=x[i,],W=w,P=tau[i])[[1]]
   }
return(q)}


SPQR_probs <- function(w,x,Z1=NULL){
   if(is.null(Z1)){
     Z1 <- act(sweep(x%*%w[[2]],2,w[[1]],"+")) # n x K
   }
   p <- sweep(Z1%*%w[[4]],2,w[[3]],"+")
   p <- exp(sweep(p,1,p[,1],"-"))
   p <- sweep(p,1,rowSums(p),"/")
return(p)}


SPQR_loss <- function(w,x,B,lam=rep(0,4)){
  -sum(log(rowSums(B*SPQR_probs(w,x)))) + 
   sum(lam*unlist(lapply(w,Pen)))
}

SPQR_grad <- function(w,x,B,lam=rep(0,4)){
  g   <- lapply(w,Pen_grad)
  for(j in 1:4){g[[j]] <- lam[j]*g[[j]]}

  Z1     <- act(sweep(x%*%w[[2]],2,w[[1]],"+"))
  Z1p    <- actp(sweep(x%*%w[[2]],2,w[[1]],"+"))
  p      <- SPQR_probs(w,x,Z1=Z1)
  Bp     <- B*p
  f      <- rowSums(Bp)
  Bp2    <- sweep(Bp,1,f,"/")
  Z1pf   <- sweep(Z1p,1,f,"/")
  dfdg2  <- Bp2-p
  g[[3]] <- g[[3]] - colSums(dfdg2)
  g[[4]] <- g[[4]] - t(Z1)%*%dfdg2

  for(k in 1:length(w[[1]])){
     pw         <- p%*%(w[[4]][k,])
     Bpw        <- Bp%*%w[[4]][k,] - f*pw 
     temp       <- Bpw*Z1pf[,k] 
     g[[1]][k]  <- g[[1]][k] - sum(temp)
     g[[2]][,k] <- g[[2]][,k] - as.vector(t(x)%*%temp)
  }
return(g)}



if(FALSE){

 p   <- 3
 K   <- 4
 L   <- 5
 n   <- 6
 x   <- matrix(rnorm(n*p),n,p)
 y   <- rbeta(n,2,3)
 w   <- SPQR_init(p,L,K,sigma=1)
 lam <- 1:4

 B   <- basis(y,L)
 l0  <- SPQR_loss(w,x,B,lam)
 
 g   <- w

 numgrad2 <- function(w,x,B,lam,k,l0,eps){
   g <- w[[k]]
   if(is.vector(g)){
      for(i in 1:length(g)){
        We         <- w
        We[[k]][i] <- We[[k]][i] + eps
        g[i]       <- SPQR_loss(We,x,B,lam)
      }
   }
   if(is.matrix(g)){
      for(i in 1:nrow(g)){for(j in 1:ncol(g)){
        We           <- w
        We[[k]][i,j] <- We[[k]][i,j] + eps
        g[i,j]       <- SPQR_loss(We,x,B,lam)
      }}
   }
 return((g-l0)/eps)}

 
 eps  <- 10^(-8)
 g1   <- w

 for(k in 1:4){g1[[k]]<-numgrad2(w,x,B,lam,k,l0,eps)}
 g2 <- SPQR_grad(w,x,B,lam)
 for(j in 1:4){
   print(paste0("w",j))
   print(g1[[j]])
   print(g2[[j]])
 }
}

if(FALSE){
  source("C:\\Users\\bjreich\\Desktop\\NormalizingFlows\\Code\\adam.R")
  n   <- 100000
  p   <- 2
  L   <- 20
  K   <- 25
  lam <- rep(.01,4)
  X   <- cbind(runif(n),0)
  y   <- rnorm(n,3*X[,1]^2,1.2-X[,1])
  Y   <- pnorm(y,mean(y),sd(y))
  plot(X[,1],Y)

  initw <- SPQR_init(p,L,K)
  B     <- basis(Y,L)
  fit   <- adam(initw, X, B, loss=SPQR_loss, grad=SPQR_grad,lam=lam,lr=0.001)
  
  ygrid  <- seq(0,1,0.01)
  plot(NA,xlim=0:1,ylim=c(0,5),xlab="y",ylab="PDF",main="PDF estimation")
  for(x0  in c(0.2,0.5,0.8)){
    X0   <- matrix(c(x0,0),length(ygrid),2,byrow=T)
    den  <- density(Y[abs(X[,1]-x0)<0.01])
    pdf  <- dSPQR(ygrid,X0,fit$w)
    lines(den$x,den$y)
    lines(ygrid,pdf,col=2)  
  }

  xgrid  <- seq(0,1,.05)
  Xgrid  <- cbind(xgrid,0)
  XX     <- round(X[,1],1)
  boxplot(Y~XX,xlim=0:1,outline=FALSE,at=sort(unique(XX)),boxwex=.05,
          xlab="X1",ylab="Quantile",main="Quantile estimation")
  for(tau in c(0.05,0.5,0.95)){ 
    q <- qSPQR(rep(tau,nrow(Xgrid)),Xgrid,fit$w)
    lines(xgrid,q,col=2,lwd=2)
  }


}



