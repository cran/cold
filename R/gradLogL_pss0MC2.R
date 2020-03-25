gradLogL.pss0MC2<- function(parameters,X,Z,data,M,trace)
{
  
  gradient <- function(param, X, y,M, b1j.m, b2j.m, pos.r2)
  {
    y[is.na(y)] <- (-1)
    y <- as.integer(y)
    n <-as.integer(length(y))
    npar <- as.integer(length(param)-1)
    beta <- as.double(param[1:(npar-1)])
    rho <- as.double(0)
    theta <- work <- as.double(rep(0, n))
    grad <- as.double(rep(0, npar))
    x <-matrix(as.double(X),nrow=n,ncol=npar-1)
    s.grad <-n.grad<- as.double(rep(0,npar))
    s.prod1<-s.prod2<-n.omega1<-n.omega2<-as.double(0)
    omega1<-as.double(param[length(param)-1])
    omega2<-as.double(param[length(param)])
    m <- max(y)
    fact <- rep(1, m + 1)
    if(m > 0)
      for(i in 2:m + 1)
        fact[i] <- fact[i - 1] * (i - 1)
    fact <- as.double(fact)
    link <- as.integer(1)
    Lik <- as.double(0)
    
    # Monte Carlo method
    for (j in 1:M)
    {	b1j<-b1j.m[j]
    b2j<-b2j.m[j]
    beta[1]<-as.double(param[1] + b1j)
    beta[pos.r2]<-as.double(param[pos.r2] + b2j)
    
    results <- .Fortran("psslik0",logL,beta,
                        npar,x,y,theta,work,n,fact,link,PACKAGE="cold")
    
    results1 <- .Fortran("pssgrd0",grad,beta,
                         npar,x,y,theta,work,n,link,PACKAGE="cold")
    
    if  (results[[1]]=="NaN" )   results[[1]]<-(-Inf)
    else if  (results[[1]]>700 )   results[[1]]<-(700)
    
    for (i in 1:length(grad))
    {if  (results1[[1]][i]=="NaN") results1[[1]][i]<-0}
    
    Lik <- Lik + exp(results[[1]] )
    s.grad<- s.grad + results1[[1]]*exp(results[[1]])
    s.prod1<-s.prod1+exp(results[[1]])*((b1j^2-exp(omega1))/(2*exp(2*omega1)))
    s.prod2<-s.prod2+exp(results[[1]])*((b2j^2-exp(omega2))/(2*exp(2*omega2)))
    }
    
    n.grad<- (s.grad/M)
    n.gradf<-n.grad[1:(npar-1)]
    n.omega1<- exp(omega1)*(s.prod1/M)
    n.omega2<- exp(omega2)*(s.prod2/M)
    L<- (Lik/M)
    
    return(c(n.gradf,n.omega1,n.omega2,L))
  }
  
  param<-parameters
  ti.repl <- data[[1]]
  cumti.repl <- cumsum(ti.repl)
  n.cases <- length(ti.repl)
  y <- data[[2]]
  derivatives<-rep(as.double(0),length(param))
  der.grad<-rep(as.double(0),(length(param)-2))
  der.omega1<-der.omega2<-as.double(0)
  k1 <- 1
  logL<-0
  omega1<-as.double(param[length(param)-1])
  omega2<-as.double(param[length(param)])
  mvcov<-matrix(as.double(0),2,2)
  b1j.m<-rep(as.double(0),M)
  b2j.m<-rep(as.double(0),M)
  pos.r2<-as.double(0)
  
  names.Z <- dimnames(Z)[[2]]  #new for random
  names.X <- dimnames(X)[[2]]
  
  for (i in 2:ncol(X))
  { if (!is.na(match(names.Z[2],names.X[i]))) pos.r2<-i  }
  
  for (i in 1:n.cases)
  {
    k2<-cumti.repl[i]
    set.seed(10*i)
    mvcov<- matrix(c(exp(omega1), 0, 0, exp(omega2)),2)
    bi<-mvrnorm(n=M, c(0,0), mvcov)
    b1j.m<-bi[,1]
    b2j.m<-bi[,2] 
    
    numerator<-gradient(param=parameters, X=X[k1:k2,], y=y[k1:k2], M=M, 
                        b1j.m=b1j.m, b2j.m=b2j.m, pos.r2=pos.r2)
    
    z<-numerator[length(param)+1]
    
    {if  (z=="Inf" ) z<-(1e+150)}
    
    #      logL<-logL+  (z) #log-likelihood for N individuals
    
    for (k in 1:(length(param)-2))
    {
      if (is.na(numerator[k]) | is.na(numerator[k]-Inf)) numerator[k]<-0
      if (numerator[k]=="Inf")  numerator[k]<- (1e+150)
      if (numerator[k]=="-Inf") numerator[k]<- (-1e+150)
      
      der.grad[k]<-der.grad[k] + (numerator[k]/z)
      k<-k+1
    }
    
    if (is.na(numerator[length(param)-1]) | is.na(numerator[length(param)-1]-Inf)) numerator[length(param)-1]<-0
    if (numerator[length(param)-1]=="Inf")  numerator[length(param)-1]<- (1e+150)
    if (numerator[length(param)-1]=="-Inf") numerator[length(param)-1]<- (-1e+150)
    
    
    if (is.na(numerator[length(param)]) | is.na(numerator[length(param)]-Inf)) numerator[length(param)]<-0
    if (numerator[length(param)]=="Inf")  numerator[length(param)]<- (1e+150)
    if (numerator[length(param)]=="-Inf") numerator[length(param)]<- (-1e+150)
    
    der.omega1<-der.omega1 + (numerator[length(param)-1]/z)
    der.omega2<-der.omega2 + (numerator[length(param)]/z)
    
    k1<-k2+1
  }
  
  derivatives<-c(der.grad,der.omega1,der.omega2)

  return(-derivatives)
}