LogL.pss0I.aux<- function(parameters, X, data, trace)
{
  loglik1 <- function(param, X, y,trace)
  {
    
    npar <-as.integer(length(param))
    beta<- as.double(param[1:(npar-1)])
    bt<- as.double(param[1:(npar-1)])
    omega<-as.double(param[npar])
    y<- as.integer(y)
    n <- as.integer(length(y))
    x<-matrix(as.double(X),nrow=n,ncol=npar-1)
    theta<- work<-as.double(rep(0,n))
    logL <- as.double(0)
    eta<-i.fit<- fit<-as.vector(n)
    eta<-x%*%beta

##### Estimate bi 
    m<-glm(as.numeric(y)~ offset(eta) , family=poisson)

    bi<-coef(m)[1]
    beta[1]<-beta[1]+bi   
    y[is.na(y)]<-(-1)
    y<- as.integer(y)
    link <- as.integer(1)
    m <- max(y)
    fact <- rep(1, m + 1)
    if(m > 0)
    {for(i in 2:(m + 1))
      fact[i] <- fact[i - 1] * (i - 1)}
    
    fact <- as.double(fact)
    
    results <- .Fortran("psslik0",logL,beta,
                        npar,x,y,theta,work,n,fact,link,PACKAGE="cold")
    
   i.fit<-x%*%beta

    return(list(loglik=results[[1]], fit=i.fit, bi.est=bi))}
  
  if(trace)	cat(paste("\t",(format(parameters[length(parameters)], digit=4)), collapse=" "), "\t")
  
  nparam<-length(parameters)
  omega1<-as.double(parameters[nparam])
  ti.repl<-data[[1]]
  cumti.repl<-cumsum(ti.repl)
  n.cases<- length(ti.repl)
  y<-data[[2]]
  fitted<-as.double(rep(0,length(y)))
  res<-as.double(rep(0,length(y)))
  bi.estimate<-as.double(rep(0,length(n.cases)))
  logL1<-as.double(0)
  
  
  k1<-1  
  for (i in 1:n.cases)
  {
    k2<-cumti.repl[i]
    z<-loglik1(param=parameters,X=X[k1:k2,],  y=y[k1:k2],trace=trace)
    fitted[k1:k2]<-z$fit
    bi.estimate[i]<-z$bi.est
    k1<-k2+1
  }

  return(list(fit=fitted,bi.est=bi.estimate))}