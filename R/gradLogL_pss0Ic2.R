gradLogL.pss0Ic2 <- function(parameters, X,Z,data,cublim, trace)
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
  
  

  #int.deriv calculates derivatives 
  int.deriv<-function(v,parameters,X,y,pos.r2)
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
    daux<-matrix(as.double(0),nrow=nparam-2,ncol=k1)
    am<-matrix(as.double(0),nrow=nparam+1,ncol=k1)
    
    #creates a expression to integrate in a vector of length(v)
    for(j in 1:k1)
    {param[1]<-as.double(parameters[1]+v[1,j])
    param[pos.r2]<-as.double(parameters[pos.r2]+v[2,j])
    
    z[j]<-FUN(param,X,y)
    grad1<-FUN1(param,X,y)

    for (i0 in 1:(nparam-2)) 
    daux[i0,j]<-grad1[i0]
    }
    
    #just derivatives for beta
    for (i1 in 1:(nparam-2))
    {
    d0<-daux[i1,]
    
    #just derivatives for beta 
    a<-exp(z- ((v[1,]^2/(2*exp(omega1))) +  (v[2,]^2/(2*exp(omega2))) ))*d0
    
    am[i1,]<- matrix(a,ncol=ncol(v))   
    }
    
     #int.deriv.var1
    a<- exp (z- ((v[1,]^2/(2*exp(omega1)))   +   (v[2,]^2/(2*exp(omega2))) ))*
      ((v[1,]^2-exp(omega1))/(2*exp(2*omega1)))
    
    am[nparam-1,]<- matrix(a,ncol=ncol(v))
    
    #int.deriv.var2
    a<- exp (z- ((v[1,]^2/(2*exp(omega1)))   +   (v[2,]^2/(2*exp(omega2))) ))*
      ((v[2,]^2-exp(omega2))/(2*exp(2*omega2)))
    
    am[nparam,]<- matrix(a,ncol=ncol(v)) 
    
    #loglik
    a<- exp(z- ((v[1,]^2/(2*exp(omega1))) + (v[2,]^2/(2*exp(omega2))) ))
    
    am[nparam+1,]<- matrix(a,ncol=ncol(v))
 
    return(am)    
  }
  
  
  nparam <- as.integer(length(parameters)-1)
  omega1<-parameters[nparam]
  omega2<-parameters[nparam+1]
  ti.repl<-data[[1]]
  cumti.repl<-cumsum(ti.repl)
  n.cases<- length(ti.repl)
  y<-data[[2]]
  dgr<-as.double(rep(0,nparam-1))
  dvar1<-0
  dvar2<-0
  k1<-1
  pos.r2<-as.double(0)
  
  l1i<-as.double(cublim$l1i)
  l1s<-as.double(cublim$l1s)
  l2i<-as.double(cublim$l2i)
  l2s<-as.double(cublim$l2s)
  
  names.Z <- dimnames(Z)[[2]]  #new for random
  names.X <- dimnames(X)[[2]]
  
  for (i in 2:ncol(X))
  { if (!is.na(match(names.Z[2],names.X[i]))) pos.r2<-i  }
  
  for (i in 1:n.cases)
  {
    k2<-cumti.repl[i]
    
       deriv<- hcubature (int.deriv,lowerLimit=c(l1i*exp(omega1/2),l2i*exp(omega2/2)),
                         upperLimit=c(l1s*exp(omega1/2),l2s*exp(omega2/2)), fDim=nparam+2,
                         parameters=parameters,X=X[k1:k2,], y=y[k1:k2], vectorInterface=TRUE,
                         pos.r2=pos.r2)

     for (i0 in 1:(nparam+2))
     {  
      if (is.na(deriv[[1]][i0]) | is.na(deriv[[1]][i0]-Inf)) deriv[[1]][i0]<-0
      else if (deriv[[1]][i0]=="Inf")  deriv[[1]][i0]<- (1e+150)
      else if (deriv[[1]][i0]=="-Inf") deriv[[1]][i0]<- (-1e+150)
     }
       
      #loglik
      z<-deriv[[1]][nparam+2]
      
      for (i1 in 1:(nparam-1))
      {    
      dgr[i1]<-dgr[i1]+(deriv[[1]][i1]/z)
      }
    
    dvar1<-dvar1+ (deriv[[1]][nparam]/z)*exp(omega1) #using the chain rule

    dvar2<-dvar2+ (deriv[[1]][nparam+1]/z)*exp(omega2) #using the chain rule
    
    k1<-k2+1
  }
  
  gr<-c(dgr,dvar1,dvar2)

  return(-gr)}

