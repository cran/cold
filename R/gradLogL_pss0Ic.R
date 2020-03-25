gradLogL.pss0Ic <- function(parameters, X,Z,data,cublim,trace)
{
  
  gradient <- function(param,  X, y)
  {
    npar <- as.integer(length(param))
    beta <- as.double(param[1:(npar-1)])
    rho<-as.double(0)
    y[is.na(y)] <- (-1)
    y <- as.integer(y)
    n <- as.integer(length(y))
    theta <- work <- as.double(rep(0, n))
    grad <- as.double(rep(0, npar))
    x <-matrix(as.double(X),nrow=n,ncol=npar-1)
    m <- max(y)
    fact <- rep(1, m + 1)
    if(m > 0)
    {for(i in 2:m + 1)
      fact[i] <- fact[i - 1] * (i - 1)}
    fact <- as.double(fact)
    link <- as.integer(1)
    
    result <- .Fortran("pssgrd0",grad,beta,
                       npar,x,y,theta,work,n,link,PACKAGE="cold")
    
    for (i in 1:length(grad))
    {if  (result[[1]][i]=="NaN" ) result[[1]][i]<-0}
    
    return(result[[1]])}
  
  loglik<- function(param, X, y)
  {#calculate logLi(beta/bi) for each individual
    
    npar <-as.integer(length(param))
    beta<- as.double(param[1:(npar-1)])
    y[is.na(y)]<-(-1)
    y<- as.integer(y)
    n <- as.integer(length(y))
    x<-matrix(as.double(X),nrow=n,ncol=npar-1)
    theta<- work<- as.double(rep(0,n))
    logL <- as.double(0)
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
    z[j]<-FUN(param,X,y) }
    
    a<-exp(z-(v^2)/(2*exp(omega1)))
    
    return(a) }
  
  
  #int.deriv calculates derivatives of beta 
  int.deriv<-function(v,parameters,X,y,i1)
  {
    FUN<-get("loglik", inherits=TRUE)
    FUN1<-get("gradient", inherits=TRUE)
    param<- parameters
    omega<-as.double(parameters[length(param)])
    k<-length(v)
    grad1<-as.vector(length(param))
    z<-as.vector(length(v))
    d0<-as.vector(length(v))
    
    for(j in 1:k)
    {param[1]<-as.double(parameters[1]+v[j])
    z[j]<-FUN(param,X,y)
    grad1<-FUN1(param,X,y)
    #just derivatives for beta0
    d0[j]<-grad1[i1]     }
    
    a<-exp(z-(v^2)/(2*exp(omega)))*d0
    return(a)   }
  
  #int.deriv.var calculates derivatives of variance
  int.deriv.var<-function(v,parameters,X,y)
  {
    FUN<-get("loglik", inherits=TRUE)
    param<- parameters
    omega<-as.double(parameters[length(param)])
    k<-length(v)
    z<-as.vector(length(v))
    
    for(j in 1:k)
    {param[1]<-as.double(parameters[1]+v[j])
    z[j]<-FUN(param,X,y)    }
    
    a<-exp(z-(v^2)/(2*exp(omega)))*((v^2-exp(omega))/(2*exp(2*omega)))
    
    return(a)   }
  
  nparam <- as.integer(length(parameters)-1)
  omega1<-as.double(parameters[nparam+1])
  ti.repl<-data[[1]]
  cumti.repl<-cumsum(ti.repl)
  n.cases<- length(ti.repl)
  y<-data[[2]]
  dgr<-as.double(rep(0,nparam))
  dvar<-0
  k1<-1
  
  l1i<-as.double(cublim$l1i)
  l1s<-as.double(cublim$l1s)
  
  for (i in 1:n.cases)
  {
    k2<-cumti.repl[i]
    
    z<- hcubature (int0c,lowerLimit=l1i*exp(omega1/2),upperLimit=l1s*exp(omega1/2),
                   parameters=parameters,X=X[k1:k2,], y=y[k1:k2],vectorInterface=TRUE)
    
    {if  (z[[1]]=="Inf" ) z[[1]]<-(1e+150)}
    
    for (i1 in 1:nparam)
    {
      deriv<-hcubature (int.deriv,lowerLimit=l1i*exp(omega1/2),upperLimit=l1s*exp(omega1/2),
                        parameters=parameters,X=X[k1:k2,], y=y[k1:k2],i1=i1,vectorInterface=TRUE)
      
      if (is.na(deriv[[1]]) | is.na(deriv[[1]]-Inf)) deriv[[1]]<-0
      else if (deriv[[1]]=="Inf")  deriv[[1]]<- (1e+150)
      else if (deriv[[1]]=="-Inf") deriv[[1]]<- (-1e+150)
      
      dgr[i1]<-dgr[i1]+(deriv[[1]]/z[[1]])
    }
    
    deriv.var<- hcubature (int.deriv.var,lowerLimit=l1i*exp(omega1/2),upperLimit=l1s*exp(omega1/2),
                           parameters=parameters,X=X[k1:k2,], y=y[k1:k2],vectorInterface=TRUE)
    
    if (is.na(deriv.var[[1]]) | is.na(deriv.var[[1]]-Inf)) deriv.var[[1]]<-0
    if (deriv.var[[1]]=="Inf")  deriv.var[[1]]<- (1e+150)
    if (deriv.var[[1]]=="-Inf") deriv.var[[1]]<- (-1e+150)
    
    dvar<-dvar+(deriv.var[[1]]/z[[1]])*exp(omega1) #using the chain rule
    
    k1<-k2+1
  }
  gr<-c(dgr,dvar)
  
  return(-gr)}
