LogL.pss1Ic2<- function(parameters, X, Z, data, trace,cublim)
{
  
  loglik<- function(param, X, y)
  {#calculate logLi(beta/bi) for each individual
    
    npar <-as.integer(length(param)-2)
    beta<- as.double(param[1:(npar-1)])
    rho<-as.double(param[npar])
    y[is.na(y)]<-(-1)
    y<- as.integer(y)
    n <- as.integer(length(y))
    x<-matrix(as.double(X),nrow=n,ncol=npar-1)
    theta<- work<- as.double(rep(0,n))
    m.theta <-as.double(rep(0,n))
    vt<-as.double(rep(0,n-1))
    logL <- as.double(0)
    Lik <- as.double(0)
    
    link <- as.integer(1)
    m <- max(y)
    fact <- rep(1, m + 1)
    if(m > 0)
    {for(i in 2:(m + 1))
      fact[i] <- fact[i - 1] * (i - 1)}
    fact <- as.double(fact)
    
    m.theta<-exp(x%*%beta)
    for (i in 2:n)
    {vt[i-1]<- m.theta[i]-rho*m.theta[i-1]}

    if ( all(vt > 0) &  all(vt != Inf) &  !anyNA(vt) )
    {results <- .Fortran("psslik",logL,beta,rho,
                        npar,x,y,theta,work,n,fact,link,PACKAGE="cold")
    
    if  (results[[1]]=="NaN" ) results[[1]]<-(-Inf)
    if  (results[[1]]=="Inf" ) results[[1]]<-(-Inf)
   
    Lik<-results[[1]]  }
   
    else Lik<-(-Inf)
    
    return(Lik) }
  
  
  int1c<-function(v,parameters,X,y,pos.r2)
  {
    FUN<-get("loglik", inherits=TRUE)
    
    param<- parameters
    nparam<-length(parameters)
    omega1<-as.double(parameters[nparam-1])
    omega2<-as.double(parameters[nparam])
    k<-length(v)
    k1<-as.double(ncol(v))
    z<-as.vector(length(k1))
    
    #creates a expression to integrate in a vector of length(v)
    for(j in 1:k1)
    {param[1]<-as.double(parameters[1]+v[1,j])
    param[pos.r2]<-as.double(parameters[pos.r2]+v[2,j])
    z[j]<-FUN(param,X,y)
    }

    a<- exp (z- ((v[1,]^2/(2*exp(omega1))) + (v[2,]^2/(2*exp(omega2))) ))
    
    am<- matrix(a,ncol=ncol(v))
    
    return(am)
  }
  
  if(trace)	cat(paste(format(parameters[length(parameters)-2], digit=6), collapse=" "), "\t")
  if(trace)	cat(paste(format(parameters[length(parameters)-1], digit=6), collapse=" "), "\t")
  if(trace)	cat(paste(format(parameters[length(parameters)], digit=6), collapse=" "), "\t")
  
  nparam<-length(parameters)
  omega1<-as.double(parameters[nparam-1])
  omega2<-as.double(parameters[nparam])
  ti.repl<-data[[1]]
  cumti.repl<-cumsum(ti.repl)
  n.cases<- as.integer(length(ti.repl))
  y<-data[[2]]
  logL1<-0
  k1<-1
  pos.r2<-as.double(0)
  npar <- as.integer(length(parameters))
  rho<- as.double(parameters[npar-2])
  
  l1i<-as.double(cublim$l1i)
  l1s<-as.double(cublim$l1s)
  l2i<-as.double(cublim$l2i)
  l2s<-as.double(cublim$l2s)
  
  names.Z <- dimnames(Z)[[2]]  #new for random
  names.X <- dimnames(X)[[2]]
  
  for (i in 2:ncol(X))
  { if (!is.na(match(names.Z[2],names.X[i]))) pos.r2<-i  }
  
  if (omega1 > 10 | omega2 > 10 )
  { logL1<-NaN
  if(trace)	cat(paste("\t",(format( logL1,digit=6)), collapse=" "), "\n")
  return(NaN)}
  
  if (omega1 < -10 | omega2 < -10 )
  { logL1<-NaN
  if(trace)	cat(paste("\t",(format( logL1,digit=6)), collapse=" "), "\n")
  return(NaN)}
  
  if (rho < 0 | rho >= 0.99 )
  { logL1<-NaN
  if(trace)	cat(paste("\t",(format( logL1,digit=6)), collapse=" "), "\n")
  return(NaN)}
  
  else if (rho >= 0 &  rho < 0.99 )
  {
    for (i in 1:n.cases)
    {
      k2<-cumti.repl[i]
      
      z<- hcubature (int1c,lowerLimit=c(l1i*exp(omega1/2),l2i*exp(omega2/2)),
                     upperLimit=c(l1s*exp(omega1/2),l2s*exp(omega2/2)),
                     parameters=parameters,X=X[k1:k2,], y=y[k1:k2],
                     vectorInterface=TRUE, pos.r2=pos.r2)

    {if  (z[[1]]=="0" )  z[[1]]<-(1e-300)}
     if  (z[[1]]=="Inf" )  z[[1]]<-(1e+150)
      if  (z[[1]]=="NaN" )  z[[1]]<-(1e-150)
#      if  (z[[1]]>700)  z[[1]]<-700
      
      if(z[[4]]==0)
        #logL1 gives the log-likelihood
        logL1<-logL1+ log(z[[1]]*(1/(2*pi*exp(omega1/2)*exp(omega2/2))))
      
      k1<-k2+1
    }
    
    if(trace)	cat(paste("\t",(format( logL1,digit=6)), collapse=" "), "\n")
    
    return(-logL1)}
}
