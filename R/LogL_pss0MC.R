LogL.pss0MC<- function(parameters, X, Z, data,M,trace)
{
  
  loglik1 <- function(param, X, y,M,bj.m)
  {
    y[is.na(y)]<-(-1)
    y<- as.integer(y)
    n <- as.integer(length(y))
    npar <-as.integer(length(param))
    beta<- as.double(param[1:(npar-1)])
    rho<-as.double(0)
    omega<-as.double(param[npar])
    x<-matrix(as.double(X),nrow=n,ncol=npar-1)
    theta<- work<- as.double(rep(0,n))
    Lik <- as.double(0)
    logL <- as.double(0)
    link <- as.integer(1)
    m <- max(y)
    fact <- rep(1, m + 1)
    if(m > 0)
      for(i in 2:(m + 1))
        fact[i] <- fact[i - 1] * (i - 1)
    fact <- as.double(fact)

    # Monte Carlo method
    for (j in 1:M)
    { bj<-bj.m[j]
    beta[1]<-as.double(param[1] + bj)
    
    results <- .Fortran("psslik0",logL,beta,
                        npar,x,y,theta,work,n,fact,link,PACKAGE="cold")
  
    {if  (results[[1]]=="NaN" ) results[[1]]<-(-Inf)}

    Lik <- Lik + exp( results[[1]] )
    }
    
    L<- (Lik/M)

    return(log(L))
  }
  
  if(trace)	cat(paste("\t",(format(parameters[length(parameters)], digit=4)), collapse=" "), "\t")
  
  param<-parameters
  ti.repl<-data[[1]]
  cumti.repl<-cumsum(ti.repl)
  n.cases<- as.integer(length(ti.repl))
  y<-data[[2]]
  counts<-data[[3]]
  logL<-as.double(0)
  k1<-1
  omega<-as.double(param[length(param)])
  bj.m<-rep(as.double(0),M)
  npar <- as.integer(length(parameters))
  
  for (i in 1:n.cases)
  {
    k2<-cumti.repl[i]
    set.seed(10*i)
    bj.m<-rnorm(M,0,exp(omega/2))
 
    z<- loglik1(param=parameters,X=X[k1:k2,],y=y[k1:k2],M=M, bj.m=bj.m)
    
    logL<-logL+counts[i]*z

    k1<-k2+1
  }
  
  if(trace)	cat(paste("\t",(format( logL,digit=6)), collapse=" "), "\n")
  
  return(-logL)
}
