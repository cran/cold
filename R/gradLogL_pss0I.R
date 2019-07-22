gradLogL.pss0I <- function(parameters, X,data,integrate,trace)
{
  
  gradient1 <-  function(param,X,y,integrate)
  {
    y[is.na(y)]<-(-1)
    y <- as.integer(y)
    n <- as.integer(length(y))
    npar <- as.integer(length(param))
    beta <- as.double(param[1:(npar-1)])
    bt <- as.double(param[1:(npar-1)])
    omega<-as.double(param[npar])
    theta <- work <- as.double(rep(0,n))
    grad<- as.double(rep(0,npar-1))
    gvar<-as.double(0)
    x <- matrix(as.double(X),nrow=n, ncol=npar-1)
    m <- max(y)
    
    li<-as.double(integrate$lig)
    ls<-as.double(integrate$lsg)
    epsabs<-as.double(integrate$epsabs)
    epsrel<-as.double(integrate$epsrel)
    limit<-as.integer(integrate$limit)
    key<-as.integer(integrate$key)
    link <- as.integer(1)
    
    result <- .Fortran("gintp0",grad,gvar,bt,beta,omega,npar,link,
                       m,x,y,theta,work,n,li,ls,epsabs,epsrel,key,limit,PACKAGE="cold")
    
    return(c(result[[1]],result[[2]]))	}
  
  
  loglik1 <- function(param, X, y,integrate)
  {
    npar <-as.integer(length(param))
    beta<- as.double(param[1:(npar-1)])
    bt<- as.double(param[1:(npar-1)])
    omega<-as.double(param[npar])
    y[is.na(y)]<-(-1)
    y<- as.integer(y)
    n <- as.integer(length(y))
    x<-matrix(as.double(X),nrow=n,ncol=npar-1)
    theta<- work<- as.double(rep(0,n))
    logL <- as.double(0)
    m <- max(y)
    
    li<-as.double(integrate$lig)
    ls<-as.double(integrate$lsg)
    epsabs<-as.double(integrate$epsabs)
    epsrel<-as.double(integrate$epsrel)
    limit<-as.integer(integrate$limit)
    key<-as.integer(integrate$key)
    link <- as.integer(1)
    
    results <- .Fortran("intp0",logL,bt,beta,omega,npar,link,m,
                        x,y,theta,work,n,li,ls,epsabs,epsrel,key,limit,PACKAGE="cold")
    
    return(results[[1]])}
  
  nparam <- as.integer(length(parameters))
  omega1<-parameters[nparam]
  ti.repl<-data[[1]]
  cumti.repl<-cumsum(ti.repl)
  n.cases<- length(ti.repl)
  y<-data[[2]]
  counts<-data[[3]]
  dgr<-as.double(rep(0,nparam-1))
  dvar<-0
  k1<-1
  
  for (i in 1:n.cases)
  {
    k2<-cumti.repl[i]
    
    z<-loglik1(param=parameters,X=X[k1:k2,], y=y[k1:k2],integrate=integrate)
    
    num<-gradient1(param=parameters,X=X[k1:k2,], y=y[k1:k2],integrate=integrate)

    {if  (z=="Inf" )  z<-(1e+300)} # New
    
    for (j in 1:(nparam-1))
    { 
      if (is.na(num[j]) | is.na(num[j]-Inf)) num[j]<-0 
      dgr[j]<-dgr[j]+counts[i]*(num[j]/z)
    }
    
    #using the chain rule
    if (is.na(num[nparam]) | is.na(num[nparam]-Inf)) num[nparam]<- 0
    dvar<-dvar+counts[i]*(num[nparam]/z)*exp(omega1)
    k1<-k2+1
  }
  
  gr<-c(dgr,dvar)

  return(-gr)}