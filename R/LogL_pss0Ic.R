LogL.pss0Ic<- function(parameters, X, data,trace, cublim)
{
  
  loglik<- function(param, X, y)
  {#calculate logLi(beta/bi) for each individual
    
    npar <-as.integer(length(param))
    beta<- as.double(param[1:(npar-1)])
    rho<-as.double(0)
    y[is.na(y)]<-(-1)
    y<- as.integer(y)
    n <- as.integer(length(y))
    x<-matrix(as.double(X),nrow=n,ncol=npar-1)
    theta<- work<- as.double(rep(0,n))
    logL <- as.double(0)
    Lik<-as.double(0)
    link <- as.integer(1)
    m <- max(y)
    fact <- rep(1, m + 1)
    if(m > 0)
    {for(i in 2:(m + 1))
      fact[i] <- fact[i - 1] * (i - 1)}
    
    fact <- as.double(fact)

    results <- .Fortran("psslik0",logL,beta,
                        npar,x,y,theta,work,n,fact,link,PACKAGE="cold")
    
    return(results[[1]]) }
  
 
  int0c<-function(v,parameters,X,y)
  {
    FUN<-get("loglik", inherits=TRUE)
    
    param<- parameters
    nparam<-length(parameters)
    omega1<-as.double(parameters[nparam])
    k<-length(v)
    z<-as.vector(length(v))
    
    #creates a expression to integrate in a vector of length(v)
    for(j in 1:k)
    {param[1]<-as.double(parameters[1]+v[j])
    z[j]<-FUN(param,X,y)  }
    
    a<-exp(z-(v^2)/(2*exp(omega1)))
    
    return(a) }
  
  if(trace)	cat(paste("\t",(format(parameters[length(parameters)], digit=4)), collapse=" "), "\t")
  
  nparam<-length(parameters)
  omega1<-as.double(parameters[nparam])
  ti.repl<-data[[1]]
  cumti.repl<-cumsum(ti.repl)
  n.cases<- as.integer(length(ti.repl))
  y<-data[[2]]
  counts<-data[[3]]
  logL1<-0
  k1<-1
  
  l1i<-as.double(cublim$l1i)
  l1s<-as.double(cublim$l1s)
  
  for (i in 1:n.cases)
  {
    k2<-cumti.repl[i]
    
    z<- hcubature (int0c,lowerLimit=l1i*exp(omega1/2),upperLimit=l1s*exp(omega1/2),
                   parameters=parameters,X=X[k1:k2,], y=y[k1:k2], vectorInterface=TRUE)
    
    {if  (z[[1]]=="Inf" )  z[[1]]<-(1e+150)}
    
    #logL1 gives the log-likelihood
    logL1<-logL1+counts[i]*log(z[[1]]*(1/(sqrt(2*pi)*exp(omega1/2))))
    
    k1<-k2+1
  }
  
  if(trace)	cat(paste("\t",(format( logL1,digit=6)), collapse=" "), "\n")
  return(-logL1)
}
