gradLogL.pssMC<- function(parameters,X,Z,data,M,trace)
{
  
  gradient <- function(param, X, y,M, bj.m)
  {
    y[is.na(y)] <- (-1)
    y <- as.integer(y)
    n <- length(y)
    npar <- as.integer(length(param)-1)
    beta <- as.double(param[1:(npar-1)])
    rho <- as.double(param[npar])
    omega<-as.double(param[npar+1])
    theta <- work <- as.double(rep(0, n))
    grad <- as.double(rep(0, npar))
    x <-matrix(as.double(X),nrow=n,ncol=npar-1)
    s.grad <-n.grad<- as.double(rep(0,npar))
    s.prod<-n.omega<-as.double(0)
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
    {	bj<-bj.m[j]
    beta[1]<-as.double(param[1]+bj)
    
    results <- .Fortran("psslik",logL,beta,rho,
                        npar,x,y,theta,work,n,fact,link,PACKAGE="cold")
    
    results1 <- .Fortran("pssgrd",grad,beta,rho,
                         npar,x,y,theta,work,n,fact,link,PACKAGE="cold")
    
    if  (results[[1]]=="NaN" )   results[[1]]<-(-Inf)
    else if  (results[[1]]>700 )   results[[1]]<-(700)
    
    for (i in 1:length(grad))
    {if  (results1[[1]][i]=="NaN") results1[[1]][i]<-0}
    
    Lik <- Lik + exp(results[[1]] )
    s.grad<- s.grad + results1[[1]]*exp(results[[1]])
    s.prod<-s.prod + exp(results[[1]])*((bj^2-exp(omega))/(2*exp(2*omega)))
    }
    
    n.grad<- (s.grad/M)
    n.omega<- exp(omega)*(s.prod/M)
    L<- (Lik/M)
    
    return(c(n.grad,n.omega,L))
  }
  
  param<-parameters
  ti.repl <- data[[1]]
  cumti.repl <- cumsum(ti.repl)
  n.cases <- length(ti.repl)
  y <- data[[2]]
  derivatives<-rep(as.double(0),length(param))
  der.grad<-rep(as.double(0),(length(param)-1))
  der.omega<-as.double(0)
  k1 <- 1
  logL<-0
  omega<-as.double(param[length(param)])
  bj.m<-rep(as.double(0),M)
  
  for (i in 1:n.cases)
  {
    k2<-cumti.repl[i]
    set.seed(10*i)
    bj.m<-rnorm(M,0,exp(omega/2))
    
    numerator<-gradient(param=parameters, X=X[k1:k2,], y=y[k1:k2], M=M, bj.m=bj.m)
    
    z<-numerator[length(param)+1]
    
    {if  (z=="Inf" ) z<-(1e+150)}
    
    #      logL<-logL+  (z) #log-likelihood for N individuals
    
    for (k in 1:(length(param)-1))
    {
      if (is.na(numerator[k]) | is.na(numerator[k]-Inf)) numerator[k]<-0
      if (numerator[k]=="Inf")  numerator[k]<- (1e+150)
      if (numerator[k]=="-Inf") numerator[k]<- (-1e+150)
      
      der.grad[k]<-der.grad[k] + (numerator[k]/z)
      k<-k+1
    }
    
    if (is.na(numerator[length(param)]) | is.na(numerator[length(param)]-Inf)) numerator[length(param)]<-0
    if (numerator[length(param)]=="Inf")  numerator[length(param)]<- (1e+150)
    if (numerator[length(param)]=="-Inf") numerator[length(param)]<- (-1e+150)
    
    der.omega<-der.omega + (numerator[length(param)]/z)
    
    k1<-k2+1
  }
  
  derivatives<-c(der.grad,der.omega)
  
  return(-derivatives)
}