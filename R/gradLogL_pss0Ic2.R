gradLogL.pss0Ic2 <- function(parameters, X,data,cublim, trace)
{
  
  gradient <- function(param,  X, y)
  {
    npar <- as.integer(length(param)-1)
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
  { #calculate logLi(beta/bi) for each individual
    
    npar <-as.integer(length(param)-1)
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
    
    if  (results[[1]]=="NaN" ) results[[1]]<-(-Inf)
    else if  (results[[1]]>700 )   results[[1]]<-(700)
    else if  (results[[1]]=="Inf" ) results[[1]]<-(-Inf)
    
    return(results[[1]])
  }
  
  
  int0c<-function(v,parameters,X,y)
  {
    FUN<-get("loglik",  inherits=TRUE)
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
    param[2]<-as.double(parameters[2]+v[2,j])
    z[j]<-FUN(param,X,y)  }
    
    a<- exp(z- ((v[1,]^2/(2*exp(omega1))) + (v[2,]^2/(2*exp(omega2))) ))
    
    am<-matrix(a,byrow=TRUE)
    return(am)
  }
  
  #int.deriv calculates derivatives of beta and rho
  int.deriv<-function(v,parameters,X,y,i1)
  {
    FUN<-get("loglik",  inherits=TRUE)
    FUN1<-get("gradient",  inherits=TRUE)
    param<- parameters
    nparam<-length(parameters)
    omega1<-as.double(parameters[nparam-1])
    omega2<-as.double(parameters[nparam])
    k1<-as.double(ncol(v))
    z<-as.vector(length(k1))
    d0<-as.vector(length(k1))
    grad1<-as.vector(length(param)-2)
    
    #creates a expression to integrate in a vector of length(v)
    for(j in 1:k1)
    {param[1]<-as.double(parameters[1]+v[1,j])
    param[2]<-as.double(parameters[2]+v[2,j])
    
    z[j]<-FUN(param,X,y)
    grad1<-FUN1(param,X,y)
    
    #just derivatives for beta
    d0[j]<-grad1[i1]      }
    
    a0<-exp(z- ((v[1,]^2/(2*exp(omega1))) +  (v[2,]^2/(2*exp(omega2))) ))*d0
    
    am<-matrix(a0,byrow=TRUE)
    return(am)    
  }
  
  
  #int.deriv.var calculates derivatives of variance
  int.deriv.var1<-function(v,parameters,X,y)
  {
    FUN<-get("loglik", inherits=TRUE)
    param<- parameters
    nparam<-length(parameters)
    omega1<-as.double(parameters[nparam-1])
    omega2<-as.double(parameters[nparam])
    k1<-as.double(ncol(v))
    z<-as.vector(length(k1))
    d0<-as.vector(length(k1))
    
    for(j in 1:k1)
    {param[1]<-as.double(parameters[1]+v[1,j])
    param[2]<-as.double(parameters[2]+v[2,j])
    z[j]<-FUN(param,X,y)}
    
    a0<- exp (z- ((v[1,]^2/(2*exp(omega1)))   +   (v[2,]^2/(2*exp(omega2))) ))*
      ((v[1,]^2-exp(omega1))/(2*exp(2*omega1)))
    
    am<-matrix(a0,byrow=TRUE)
    return(am)    
  }
  
  #int.deriv.var calculates derivatives of variance
  int.deriv.var2<-function(v,parameters,X,y)
  {
    FUN<-get("loglik",  inherits=TRUE)
    param<- parameters
    nparam<-length(parameters)
    omega1<-as.double(parameters[nparam-1])
    omega2<-as.double(parameters[nparam])
    k1<-as.double(ncol(v))
    z<-as.vector(length(k1))
    d0<-as.vector(length(k1))
    
    for(j in 1:k1)
    {param[1]<-as.double(parameters[1]+v[1,j])
    param[2]<-as.double(parameters[2]+v[2,j])
    z[j]<-FUN(param,X,y)}
    
    a0<- exp (z- ((v[1,]^2/(2*exp(omega1)))   +   (v[2,]^2/(2*exp(omega2))) ))*
      ((v[2,]^2-exp(omega2))/(2*exp(2*omega2)))
    
    am<- matrix(a0,byrow=TRUE)
    return(am)
  }
  
  nparam <- as.integer(length(parameters)-1)
  omega1<-parameters[nparam]
  omega2<-parameters[nparam+1]
  ti.repl<-data[[1]]
  cumti.repl<-cumsum(ti.repl)
  n.cases<- length(ti.repl)
  y<-data[[2]]
  counts<-data[[3]]
  dgr<-as.double(rep(0,nparam-1))
  dvar1<-0
  dvar2<-0
  k1<-1
  
  l1i<-as.double(cublim$l1i)
  l1s<-as.double(cublim$l1s)
  l2i<-as.double(cublim$l2i)
  l2s<-as.double(cublim$l2s)
  
  for (i in 1:n.cases)
  {
    k2<-cumti.repl[i]
    
    z<- hcubature (int0c,lowerLimit=c(l1i*exp(omega1/2),l2i*exp(omega2/2)),
                   upperLimit=c(l1s*exp(omega1/2),l2s*exp(omega2/2)),
                   parameters=parameters,X=X[k1:k2,], y=y[k1:k2], vectorInterface=TRUE)
    
      if  (z[[1]]=="0" ) z[[1]]<-(1e-300)
      else if  (z[[1]]=="Inf" ) z[[1]]<-(1e+150)
      else if  (z[[1]]=="NaN" )  z[[1]]<-(1e-150)
    
    for (i1 in 1:(nparam-1))
    {
      deriv<- hcubature (int.deriv,lowerLimit=c(l1i*exp(omega1/2),l2i*exp(omega2/2)),
                         upperLimit=c(l1s*exp(omega1/2),l2s*exp(omega2/2)),
                         parameters=parameters,X=X[k1:k2,], y=y[k1:k2],i1=i1,vectorInterface=TRUE)
      
      if (is.na(deriv[[1]]) | is.na(deriv[[1]]-Inf)) deriv[[1]]<-0
      if (deriv[[1]]=="Inf")  deriv[[1]]<- (1e+150)
      if (deriv[[1]]=="-Inf") deriv[[1]]<- (-1e+150)
      
      dgr[i1]<-dgr[i1]+counts[i] *(deriv[[1]]/z[[1]])
    }
    
    
    deriv.var1<- hcubature (int.deriv.var1,lowerLimit=c(l1i*exp(omega1/2),l2i*exp(omega2/2)),
                            upperLimit=c(l1s*exp(omega1/2),l2s*exp(omega2/2)),
                            parameters=parameters,X=X[k1:k2,], y=y[k1:k2],vectorInterface=TRUE)
    
    deriv.var2<- hcubature (int.deriv.var2,lowerLimit=c(l1i*exp(omega1/2),l2i*exp(omega2/2)),
                            upperLimit=c(l1s*exp(omega1/2),l2s*exp(omega2/2)),
                            parameters=parameters,X=X[k1:k2,], y=y[k1:k2], vectorInterface=TRUE)
    
    if (is.na(deriv.var1[[1]]) | is.na(deriv.var1[[1]]-Inf)) deriv.var1[[1]]<-0
    if (deriv.var1[[1]]=="Inf")  deriv.var1[[1]]<- (1e+150)
    if (deriv.var1[[1]]=="-Inf") deriv.var1[[1]]<- (-1e+150)
    
    dvar1<-dvar1+counts[i] * (deriv.var1[[1]]/z[[1]])*exp(omega1) #using the chain rule
    
    if (is.na(deriv.var2[[1]]) | is.na(deriv.var2[[1]]-Inf)) deriv.var2[[1]]<-0
    if (deriv.var2[[1]]=="Inf")  deriv.var2[[1]]<- (1e+150)
    if (deriv.var2[[1]]=="-Inf") deriv.var2[[1]]<- (-1e+150)
    
    dvar2<-dvar2+counts[i] * (deriv.var2[[1]]/z[[1]])*exp(omega2) #using the chain rule
    
    k1<-k2+1
  }
  
  gr<-c(dgr,dvar1,dvar2)
  
  return(-gr)}

