LogL.pss0Ic2.aux<- function(parameters, X, Z, data, trace)
{
  loglik1 <- function(param, X, Z, y, trace)
  {
    
    npar <-as.integer(length(param)-1)
    beta<- as.double(param[1:(npar-1)])
    w1<-param[npar-1]
    w2<-param[npar]
    y<- as.integer(y)
    n <- as.integer(length(y))
    x<-matrix(as.double(X),nrow=n,ncol=npar-1)
    D<-matrix(c(exp(w1),exp(w2)),nrow=2,byrow=TRUE)
    theta<- work<-as.double(rep(0,n))
    logL <- as.double(0)
    eta<-i.fit<- fit<-as.vector(npar-1)
    eta<-x%*%beta

    ui<-y-exp(eta) #residuals
    
    pos.r2<-as.double(0)
    names.Z <- dimnames(Z)[[2]]  #new for random
    names.X <- dimnames(X)[[2]]

    for (i in 2:ncol(X))
    { if (!is.na(match(names.Z[2],names.X[i]))) pos.r2<-i  }
   
   m<-glm(as.numeric(y)~  offset(eta) + Z[,2], family=poisson, maxit=100)
 
    b1i<-coef(m)[1]
    b2i<-coef(m)[2]
    beta[1]<-beta[1]+b1i
    beta[pos.r2]<-beta[pos.r2]+b2i
    y[is.na(y)]<-(-1)
    y<- as.integer(y)
    link <- as.integer(1)
    m <- max(y)
    fact <- rep(1, m + 1)
    bi<-c(b1i,b2i)

    if(m > 0)
    {for(i in 2:(m + 1))
      fact[i] <- fact[i - 1] * (i - 1)}
    fact <- as.double(fact)
    
    results <- .Fortran("psslik0",logL,beta,
                        npar,x,y,theta,work,n,fact,link,PACKAGE="cold")
    
    i.fit<-x%*%beta
    
  return(list(loglik=results[[1]], fit=i.fit, bi.est=bi))}
  
  if(trace)	cat(paste(format(parameters[length(parameters)-1], digit=4), collapse=" "), "\t")
  if(trace)	cat(paste("\t",(format(parameters[length(parameters)], digit=4)), collapse=" "), "\t")
  
  nparam<-length(parameters)
  omega1<-as.double(parameters[nparam]-1)
  omega2<-as.double(parameters[nparam])
  ti.repl<-data[[1]]
  cumti.repl<-cumsum(ti.repl)
  n.cases<- length(ti.repl)
  y<-data[[2]]
#  counts<-data[[3]]
  fitted<-as.double(rep(0,length(y)))
  bi.estimate<-matrix(as.double(0),nrow=n.cases,ncol=2)
  gi.estimate<-as.double(rep(0,length(n.cases)))
  logL1<-as.double(0)
  k1<-1

  for (i in 1:n.cases)
  {
    k2<-cumti.repl[i]
    z1<-loglik1(param=parameters,X=X[k1:k2,], Z=Z[k1:k2,], y=y[k1:k2],trace=trace)
    fitted[k1:k2]<-z1$fit
    bi.estimate[i,]<-z1$bi.est
    k1<-k2+1

  }

return(list(fit=fitted,bi.est=bi.estimate))}
