gradLogL.pss1<- function(parameters, X,data, trace)
{
  
  gradient <- function(param,  X, y)
  {
    npar <- as.integer(length(param))
    beta <- as.double(param[1:(npar-1)])
    rho <- as.double(param[npar])
    y[is.na(y)] <- (-1)
    y <- as.integer(y)
    n <- as.integer(length(y))
    theta <- work <- as.double(rep(0, n))
    grad <- as.double(rep(0, npar))
    x <-matrix(as.double(X),nrow=n,ncol=npar-1)
    m <- max(y)
    fact <- rep(1, m + 1)
    if(m > 0)
    {for (i in 2:(m + 1))
      fact[i] <- fact[i - 1] * (i - 1)}
    fact <- as.double(fact)
    link <- as.integer(1)
    
    result <- .Fortran("pssgrd",grad,beta,rho,
                       npar,x,y,theta,work,n,fact,link,PACKAGE="cold")
    
    for (i in 1:length(grad))
    {if  (result[[1]][i]=="NaN" ) result[[1]][i]<-0}
    
    return(result[[1]])
  }
  
  nparam <- length(parameters)
  ti.repl <- data[[1]]
  cumti.repl <- cumsum(ti.repl)
  n.cases <- length(ti.repl)
  y <- data[[2]]
  gr <- rep(as.double(0), nparam)
  k1 <- 1
  
  for (i in 1:n.cases)
  {
    k2<-cumti.repl[i]
    u<-gradient (param=parameters,X=X[k1:k2,], y=y[k1:k2])
    gr<-gr+u
    k1<-k2+1
  }
  
  return(-gr)}