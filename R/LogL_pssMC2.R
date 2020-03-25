LogL.pssMC2<- function(parameters, X, Z, data,M,trace)
{
  
  loglik1 <- function(param, X, y,M,b1j.m,b2j.m,pos.r2)
  {
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
    M1<-1
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
    
    m.theta<-exp(x%*%beta)
    for (i in 2:n)
    {vt[i-1]<- m.theta[i]-rho*m.theta[i-1]}
    
    if ( all(vt > 0) &  all(vt != Inf) &  !anyNA(vt) )
    {
      results <- .Fortran("psslik",logL,beta,rho,
                          npar,x,y,theta,work,n,fact,link,PACKAGE="cold")
      

      if  (results[[1]]=="-Inf" )   M1<-M1-1
      
      else if  (results[[1]]=="NaN" ) 
        {results[[1]]<-(-Inf)
      M1<-M1-1}
      
      else if  (results[[1]] > 700) results[[1]] <- 700 
      
      Lik <- Lik + exp( results[[1]])
      M1<-M1+1 }
    }
    
    if (M1=="1") L<- 1e-150
    else if (M1!="1") L<- (Lik/(M1-1))
    
    {if  (L=="Inf" )  L<-(1e+100)}
    
    #   {if  (L=="Inf" )  L<-(1e+150)}
    #   {if  (L=="0" )  L<-(1e-150)}

    return(log(L))
  }
  
  if(trace)	cat(paste(format(parameters[length(parameters)-2], digit=6), collapse=" "), "\t")
  if(trace)	cat(paste(format(parameters[length(parameters)-1], digit=6), collapse=" "), "\t")
  if(trace)	cat(paste(format(parameters[length(parameters)], digit=6), collapse=" "), "\t")
  
  param<-parameters
  ti.repl<-data[[1]]
  cumti.repl<-cumsum(ti.repl)
  n.cases<- as.integer(length(ti.repl))
  y<-data[[2]]
  logL<-as.double(0)
  k1<-1
  mvcov<-matrix(as.double(0),2,2)
  b1j.m<-rep(as.double(0),M)
  b2j.m<-rep(as.double(0),M)
  pos.r2<-as.double(0)
  omega1<-as.double(param[length(param)-1])
  omega2<-as.double(param[length(param)])
  
  names.Z <- dimnames(Z)[[2]]  #new for random
  names.X <- dimnames(X)[[2]]
  
  for (i in 2:ncol(X))
  { if (!is.na(match(names.Z[2],names.X[i]))) pos.r2<-i  }
  
  npar <- as.integer(length(parameters))
  rho<- as.double(parameters[npar-2])
  
  if (rho <= 0 |  rho >= 1  )
  { logL<-NaN
  if(trace)	cat(paste("\t",(format( logL,digit=6)), collapse=" "), "\n")
  return(NaN)}
  
  else if (rho > 0 &  rho < 1 )
  {
    for (i in 1:n.cases)
    { k2<-cumti.repl[i]
    set.seed(10*i)
    mvcov<- matrix(c(exp(omega1), 0, 0, exp(omega2)),2)
    bi<-mvrnorm(n=M, c(0,0), mvcov)
    b1j.m<-bi[,1]
    b2j.m<-bi[,2] 
    
    z<- loglik1(param=parameters,X=X[k1:k2,],y=y[k1:k2],M=M, b1j.m=b1j.m, 
                b2j.m=b2j.m, pos.r2=pos.r2)
    
    logL<-logL+ z
    
    k1<-k2+1
    }
  
    if(trace)	cat(paste("\t",(format( logL,digit=6)), collapse=" "), "\n")
    return(-logL)}
}
