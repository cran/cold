LogL.pss0MC2<- function(parameters, X, Z, data,M,trace)
{
  
  loglik1 <- function(param, X, y,M,b1j.m,b2j.m,pos.r2)
  {
    npar <-as.integer(length(param)-1)
    beta<- as.double(param[1:(npar-1)])
    y[is.na(y)]<-(-1)
    y<- as.integer(y)
    n <- as.integer(length(y))
    x<-matrix(as.double(X),nrow=n,ncol=npar-1)
    theta<- work<- as.double(rep(0,n))
    logL <- as.double(0)
    Lik <- as.double(0)
    link <- as.integer(1)
    m <- max(y)
    fact <- rep(1, m + 1)
    if(m > 0)
      for(i in 2:(m + 1))
        fact[i] <- fact[i - 1] * (i - 1)
    fact <- as.double(fact)
    
    # Monte Carlo method
    for (j in 1:M)
    {	b1j<-b1j.m[j]
    b2j<-b2j.m[j]
    beta[1]<-as.double(param[1] + b1j)
    beta[pos.r2]<-as.double(param[pos.r2] + b2j)
    
    results <- .Fortran("psslik0",logL,beta,
                        npar,x,y,theta,work,n,fact,link,PACKAGE="cold")
    
    {if  (results[[1]]=="NaN" ) results[[1]]<-(-Inf)}
    
    Lik <- Lik + exp( results[[1]])
    }
    L<- (Lik/M)
    
    {if  (L=="Inf" )  L<-(1e+150)}   

    return(log(L))    }
  
  if(trace)	cat(paste(format(parameters[length(parameters)-1], digit=6), collapse=" "), "\t")
  if(trace)	cat(paste(format(parameters[length(parameters)], digit=6), collapse=" "), "\t")
  
  param<-parameters
  ti.repl<-data[[1]]
  cumti.repl<-cumsum(ti.repl)
  n.cases<- as.integer(length(ti.repl))
  y<-data[[2]]
  counts<-data[[3]]
  logL<-as.double(0)
  k1<-1
  b1j.m<-rep(as.double(0),M)
  b2j.m<-rep(as.double(0),M)
  pos.r2<-as.double(0)
  omega1<-as.double(param[length(param)-1])
  omega2<-as.double(param[length(param)])
  
  names.Z <- dimnames(Z)[[2]]  #new for random
  names.X <- dimnames(X)[[2]]
  
  for (i in 2:ncol(X))
  { if (!is.na(match(names.Z[2],names.X[i]))) pos.r2<-i  }
  
  for (i in 1:n.cases)
  { k2<-cumti.repl[i]
  set.seed(10*i)
  b1j.m<-rnorm(M,0,exp(omega1/2))
  #     set.seed(100*i)
  b2j.m<-rnorm(M,0,exp(omega2/2))   
  
  z<- loglik1(param=parameters,X=X[k1:k2,],y=y[k1:k2],M=M,b1j.m=b1j.m, 
              b2j.m=b2j.m, pos.r2=pos.r2)
  
  logL<-logL+counts[i]*z
  
  k1<-k2+1
  }
  
  if(trace)	cat(paste("\t",(format( logL,digit=6)), collapse=" "), "\n")
  return(-logL)
}
