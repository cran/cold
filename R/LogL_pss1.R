LogL.pss1 <- function(parameters, X, data, trace)
{
  
  loglik1<- function(param, X, y, trace)
  {#calculate logLi(beta/bi) for each individual
    
    npar <- as.integer(length(param))
    beta <- as.double(param[1:(npar-1)])
    rho<- as.double(param[npar])
    y[is.na(y)]<-(-1)
    y <- as.integer(y)
    n <- as.integer(length(y))
    x <-matrix(as.double(X),nrow=n,ncol=npar-1)
    theta <- work <-as.double(rep(0,n))
    m.theta <-as.double(rep(0,n))
    vt<-as.double(rep(0,n-1))
    logL <-as.double(0)
    link <- as.integer(1)
    m <- max(y)
    fact <- rep(1, m + 1)
    if(m > 0)
    {for(i in 2:(m + 1)) 
      fact[i] <- fact[i - 1] * (i - 1)}
    fact <- as.double(fact)
    
    
    m.theta<-exp(x%*%beta)
    for (i in 2:n)
    {vt[i-1]<- m.theta[i]-rho*m.theta[i-1]
    if (vt[i-1]<0) return(NaN)}
    
    results <- .Fortran("psslik",logL,beta,rho,
                        npar,x,y,theta,work,n,fact,link,PACKAGE="cold")
    
    return(results[[1]])
  }
  
  if(trace)	cat(paste(format(parameters[length(parameters)], digit=6), collapse=" "), "\t")
  
  ti.repl<-data[[1]]
  cumti.repl<-cumsum(ti.repl)
  n.cases<- length(ti.repl)
  y<-data[[2]]
  counts<-data[[3]]
  logL1<-0
  k1<-1
  
  npar <- as.integer(length(parameters))
  rho<- as.double(parameters[npar])
  
  if (rho <= 0 |  rho >= 1 )
  { logL1<-NaN
  if(trace)	cat(paste("\t",(format( logL1,digit=6)), collapse=" "), "\n")
  return(NaN)}
  
  else if (rho > 0 &  rho < 1 )
  {
    for (i in 1:n.cases)
    {
      k2<-cumti.repl[i]
      z<- loglik1(param=parameters,X=X[k1:k2,],y=y[k1:k2],trace=trace)
      logL1<-logL1+counts[i]*z
      k1<-k2+1
    }
    
    if(trace)	cat(paste("\t",(format( logL1,digit=6)), collapse=" "), "\n")
    return(-logL1)
  }
}