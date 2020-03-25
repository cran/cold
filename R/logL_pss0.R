logL.pss0 <- function(parameters, X, data, trace)
{
  
  loglik1<- function(param, X, y, trace)
  {#calculate logLi(beta/bi) for each individual
    
    npar <- as.integer(length(param)+1)
    beta <- as.double(param[1:(npar-1)])
    
    y[is.na(y)]<-(-1)
    y <- as.integer(y)
    n <- as.integer(length(y))
    x <-matrix(as.double(X),nrow=n,ncol=npar-1)
    theta <- work <- as.double(rep(0,n))
    logL <- L<-as.double(0)
    link <- as.integer(1)
    m <- max(y)
    fact <- rep(1, m + 1)
    if(m > 0)
    {for(i in 2:(m + 1))
      fact[i] <- fact[i - 1] * (i - 1)}
    fact <- as.double(fact)
    
    results <- .Fortran("psslik0",logL,beta,
                        npar,x,y,theta,work,n,fact,link,PACKAGE="cold")
    
    return(results[[1]])
  }
  
  ti.repl<-data[[1]]
  cumti.repl<-cumsum(ti.repl)
  n.cases<- length(ti.repl)
  y<-data[[2]]
  logL<-0
  k1<-1
  
  for (i in 1:n.cases)
  {
    k2<-cumti.repl[i]
    z<- loglik1(param=parameters,X=X[k1:k2,],y=y[k1:k2],trace=trace)
    
    logL<-logL+ z
    k1<-k2+1
  }
  if(trace)	cat(paste("\t",(format( logL,digit=6)), collapse=" "), "\n")
  
  return(-logL)}