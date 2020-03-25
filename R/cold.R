

setClass("summary.cold", representation(coefficients = "matrix", se = "matrix", covariance = "matrix", correlation="matrix",
                                        log.likelihood="numeric", message ="integer",n.cases="numeric", ni.cases="numeric", aic="numeric",call="language"))

setClass("cold", representation(coefficients = "matrix", se = "matrix", covariance = "matrix", correlation="matrix",
                                log.likelihood="numeric", message ="integer",n.cases="numeric", ni.cases="numeric", aic="numeric",
                                Fitted="numeric",  bi.estimate="matrix",Fitted.av="numeric", Time="numeric", model.matrix= "matrix", y.matrix="matrix", 
                                random.matrix="matrix", subset.data="data.frame",final.data="data.frame", y.av="numeric", data.id="numeric",
                                call="language"))

setGeneric("getAIC",def=function(object) standardGeneric("getAIC"))

setGeneric("getLogLik",def=function(object) standardGeneric("getLogLik"))

setGeneric("getcoef",def=function(object) standardGeneric("getcoef"))

setGeneric("getvcov",def=function(object) standardGeneric("getvcov"))

setGeneric("randeff",def=function(object) standardGeneric("randeff"))

setGeneric("fixeff",def=function(object) standardGeneric("fixeff"))

setGeneric("vareff",def=function(object) standardGeneric("vareff"))

setGeneric("model.mat",def=function(object) standardGeneric("model.mat"))

setGeneric("coeftest",def=function(object) standardGeneric("coeftest"))

setGeneric("resid",def=function(object) standardGeneric("resid"))


cold<-function(formula,random=~0,data,id="id",time="time",subSET,
               dependence="ind",start=NULL,method="BFGS", 
               integration="QUADPACK",M=6000,control=coldControl(),
               integrate=coldIntegrate(),cublim=coldcublim(),trace=FALSE)
{
  # *****************DEFINITION OF INTERNAL FUNCTIONS ******************

  na.discrete.replace <- function(frame,  n.times, ti.repl)
  {
    vars <- names(frame)
    names(vars) <- vars
    cumti.repl<-cumsum(ti.repl)
    n.cases<- length(ti.repl)
    for(j in 1:length(vars))
    {k1<-1
    for (i in 1:n.cases)
    {k2<-cumti.repl[i]
    x <- frame[[j]][k1:k2]
    pos <- is.na(x)
    if(any(pos))
      if(j == 1) x[pos] <- -1
    else x[pos] <- x[1]

    frame[[j]][k1:k2]<-x
    k1<-k2+1
    }
    }
    return(data=frame)
  }


# ******************* MAIN PROGRAM *******************************
#
	call <- match.call()
#	vect.time <- F
	if(missing(data) || !is.data.frame(data))
		stop("a data.frame must be supplied")
	if(is.null(names(data)))
		stop("objects in data.frame must have a name")
	expr1 <- terms(formula, data=data)
	expr <- attr(expr1, "variables")
	var.names <- all.vars(expr)
	response <- all.vars(expr)[1]
	
	expr2 <- terms(random, data=data) #new for random
	names.R <- all.vars(expr2)  #new for random

	
if(!missing(time)) {Time<-as.vector(data[[time]])}
if (missing(time)) { if (all (is.na(match(names(data), "time")))) stop ("time must be defined")
			else  Time<-as.vector(data$time)}

if(!missing(id)) {id<-as.vector(data[[id]])}
if (missing(id)){ if (all(is.na(match(names(data), "id")))) stop ("id must be defined")
		else id<-as.vector(data$id)}

if(any(is.na(match(var.names, names(data)))))
	  stop("Variables in formula not contained in the data.frame")
	
# select subset if necessary
	if(!missing(subSET)) {id1 <- eval(substitute(subSET), data)
			      data<-subset(data, id1)}

#returns data of a subset
subset.data<-data

	ti.repl<-as.vector(0)
	i1<-1
	i2<-1
	for (i in 1:(length(data[[response]])-1))
	{
		if (id[i]==id[i+1])
		{ 	i2<-i2+1
			ti.repl[i1]<-i2}
		else {	ti.repl[i1]<-i2
			i1<-i1+1
			i2<-1}
	}

	n.cases <- length(ti.repl)
	n.tot<-cumsum(ti.repl)[n.cases]
	n.time<-length(unique(Time))
	ni.cases <- length(ti.repl)
	pos.ind<-cumsum(ti.repl)

	final.data <- data
	data <- data[var.names]
	n.var <- length(data)
	Y.resp <- as.vector(data[[response]])
	y1 <- Y.resp[!is.na(Y.resp)]
	if((all(y1 >= 0) | all(y1[y1 != 0] == as.integer(y1[y1 != 0]))) == FALSE)
	stop("Unfeasible values of response variable: must be non negative integers or NA"	)


# ********** creation of individual profile according to NA patterns *******************

	data2<-data
	final.data <- na.discrete.replace(frame=data,  n.times=n.time, ti.repl=ti.repl)

	data<-final.data

	
# ********** design matrices creation *******************
	# define a plausible starting point for the optimizer if not given
		data1 <- na.omit(data2)
		data1.resp <- data1[, response]
		data1[, c(response)] <- data1.resp

		if (dependence=="AR1")	 init<-0.5
		else  if (dependence=="AR1R")  init<-c(0.5,0)
		else  if (dependence=="indR")  init<-0
		else  if (dependence=="indR2")  init<-c(0,0)
		else  if (dependence=="AR1R2")  init<-c(0.5,0,0)

		if (dependence=="indR2"  &&  integration=="QUADPACK") 
		  stop ("integration argument must be MC or cubature")
		else if (dependence=="AR1R2"  &&  integration=="QUADPACK") 
      stop ("integration argument must be MC or cubature")

		
		if(is.null(start) && dependence!="ind")
			start <- c(glm(formula, data1,family=poisson)$coefficients, init)
		else if(!is.null(start) && dependence!="ind")
			start <- c(glm(formula, data1,family=poisson, maxit=100)$coefficients, start)
		else if (dependence=="ind") start <- c(glm(formula, data1, family=poisson)$coefficients)

   if (any(is.na(start))) stop("starting values produced by glm contains NA")

	id.not.na<-rep(TRUE,n.tot)
	X <- model.matrix(expr1, data)
	names.output <- dimnames(X)[[2]]
	Z <- model.matrix(expr2, data)  #new for random
	names.Z <- dimnames(Z)[[2]]  #new for random
	sum.ti <- sum(ti.repl)
	data <- list(ti.repl, data[[response]])
	data2<-list(ti.repl, data2[[response]])
	p <- dim(X)[2] + 1
	F.aux<-as.double(rep(0,length(data[[2]])))


	if (dependence=="ind")
	{	if(trace)	cat("\t log.likelihood\n")
	temp<-optim(par= start, fn=logL.pss0, gr=gradlogL.pss0, method=method,
	data = data, X = X, trace=trace,control=control)}
	else  if (dependence=="indR"& integration=="QUADPACK")
	{	if(trace)	cat("\n omega \t log.likelihood\n")
	temp <- optim(par = start, fn =LogL.pss0I,gr = gradLogL.pss0I, method=method,
	data = data, X = X, integrate=integrate, trace=trace,control=control)}
	else  if (dependence=="indR"& integration=="cubature")
	{	if(trace)	cat("\n omega \t log.likelihood\n")
	temp <- optim(par = start, fn =LogL.pss0Ic,gr = gradLogL.pss0Ic, method=method,
	data = data, X = X, Z = Z, trace=trace, cublim=cublim)}
	else  if (dependence=="indR"& integration=="MC")
	{	if(trace)	cat("\n omega \t log.likelihood\n")
	temp <- optim(par = start, fn =LogL.pss0MC,gr = gradLogL.pss0MC,  method=method,
	data = data, X = X, Z = Z, trace=trace,control=control, M=M)}
	else if (dependence=="AR1")
	{	if(trace)	cat("\n rho \t log.likelihood\n")
	temp <-optim(par= start, fn = LogL.pss1,  gr=gradLogL.pss1, method=method,
	data = data, X = X, trace=trace,control=control)}
	else  if (dependence=="AR1R" & integration=="QUADPACK")
	{	if(trace)	cat("\n rho\t omega\t log.likelihood\n")
	temp <- optim(par = start, fn =LogL.pss1I,gr = gradLogL.pss1I,  method=method,
	data = data, X = X, integrate=integrate, trace=trace,control=control)}
	else  if (dependence=="AR1R"& integration=="cubature")
	{	if(trace)	cat("\n rho\t omega\t log.likelihood\n")
	temp <- optim(par = start, fn =LogL.pss1Ic,gr = gradLogL.pss1Ic,  method=method,
	 data = data, X = X, Z = Z, trace=trace, cublim=cublim)}
	else  if (dependence=="AR1R"& integration=="MC")
	{	if(trace)	cat("\n rho\t omega\t log.likelihood\n")
	 temp <- optim(par = start, fn =LogL.pssMC1,gr = gradLogL.pssMC,  method=method,
	 data = data, X = X, Z = Z, trace=trace,control=control, M=M)}
	else  if (dependence=="indR2"& integration=="cubature")
	{	if(trace)	cat("\n  omega1\t  omega2\t log.likelihood\n")
	  temp <- optim(par = start, fn =LogL.pss0Ic2,gr = gradLogL.pss0Ic2, method=method,
	  data = data, X = X, Z = Z, trace=trace, cublim=cublim)}
	else  if (dependence=="indR2"& integration=="MC")
	{	if(trace)	cat("\n  omega1\t  omega2\t log.likelihood\n")
	  temp <- optim(par = start, fn =LogL.pss0MC2, gr = gradLogL.pss0MC2,  method=method,
	  data = data, X = X, Z = Z, trace=trace,control=control, M=M)}
	else  if (dependence=="AR1R2"& integration=="cubature")
	{	if(trace)	cat("\n rho\t  omega1\t  omega2\t  log.likelihood\n")
	  temp <- optim(par = start, fn =LogL.pss1Ic2,gr = gradLogL.pss1Ic2,  method=method,
	  data = data, X = X, Z = Z,trace=trace, cublim=cublim)}
	else  if (dependence=="AR1R2"& integration=="MC")
	{	if(trace)	cat("\n rho\t  omega1\t  omega2\t  log.likelihood\n")
	  temp <- optim(par = start, fn =LogL.pssMC2,gr = gradLogL.pssMC2,  method=method,
	  data = data, X = X, Z = Z, trace=trace,control=control, M=M)}
	
	coefficients <- temp$par
	log.lik <-  - temp$value
	
	if (trace)
	  cat("Convergence reached. Computing the information matrix now\n")

	if (dependence=="ind")
	Info <- num.info(coefficients, "gradlogL.pss0", X, data)
	else  if (dependence=="indR"& integration=="QUADPACK")
	Info <- num.infoI(coefficients, "gradLogL.pss0I", X, data, integrate=integrate)
	else  if (dependence=="indR"& integration=="cubature")
	Info <- num.infoIc(coefficients, "gradLogL.pss0Ic", X, Z, data, cublim=cublim)
	else  if (dependence=="indR"& integration=="MC")
	Info <- num.infoMC(coefficients, "gradLogL.pss0MC", X, Z, data, M=M)
	else if (dependence=="AR1")
	Info <- num.info(coefficients, "gradLogL.pss1", X, data)
	else  if (dependence=="AR1R"& integration=="QUADPACK")
	Info <- num.infoI(coefficients, "gradLogL.pss1I", X, data, integrate=integrate)
	else  if (dependence=="AR1R"& integration=="cubature")
	  Info <- num.infoIc(coefficients, "gradLogL.pss1Ic", X, Z, data, cublim=cublim)
	else  if (dependence=="AR1R"& integration=="MC")
	  Info <- num.infoMC(coefficients, "gradLogL.pssMC", X, Z, data, M=M)
	else  if (dependence=="indR2"& integration=="cubature")
	  Info <- num.infoIc(coefficients, "gradLogL.pss0Ic2", X, Z, data, cublim=cublim)
	else  if (dependence=="indR2"& integration=="MC")
	  Info <- num.infoMC(coefficients, "gradLogL.pss0MC2", X, Z, data, M=M)
	else  if (dependence=="AR1R2"& integration=="cubature")
	  Info <- num.infoIc(coefficients, "gradLogL.pss1Ic2", X, Z, data, cublim=cublim)
	else  if (dependence=="AR1R2"& integration=="MC")
	  Info <- num.infoMC(coefficients, "gradLogL.pssMC2", X, Z, data, M=M)

	se <- matrix(sqrt(diag(solve(Info))), ncol = 1)
	coefficients <- matrix(coefficients, ncol = 1)
	if (dependence=="ind")
	dimnames(coefficients) <- dimnames(se) <- list(names.output, " ")
	else  if (dependence=="indR")
	dimnames(coefficients) <- dimnames(se) <-  list(c(names.output, "omega1"), " ")
	else if (dependence=="AR1")
	dimnames(coefficients) <- dimnames(se) <- list(c(names.output, "rho"), " ")
	else  if (dependence=="AR1R")
	dimnames(coefficients) <- dimnames(se) <-  list(c(names.output, "rho","omega1"), " ")
	else  if (dependence=="indR2")
	dimnames(coefficients) <- dimnames(se) <-  list(c(names.output, "omega1", "omega2"), " ")
	else  if (dependence=="AR1R2")
	dimnames(coefficients) <- dimnames(se) <-  list(c(names.output, "rho","omega1", "omega2"), " ")
	
	covariance <- solve(Info)
	cr<- diag(1/sqrt(diag(covariance)))
	correlation <- cr %*% covariance %*% cr

	if (dependence=="ind")
	{dimnames(covariance) <- list(names.output, names.output)
	dimnames(correlation) <- list(names.output, names.output)}
	else  if (dependence=="indR")
	{dimnames(covariance) <- list(c(names.output, "omega1"), c(names.output, "omega1"))
	dimnames(correlation) <- list(c(names.output, "omega1"), c(names.output, "omega1"))}
	else if (dependence=="AR1")
	{dimnames(covariance) <- list(c(names.output, "rho"), c(names.output, "rho"))
	dimnames(correlation) <- list(c(names.output, "rho"), c(names.output, "rho"))}
	else  if (dependence=="AR1R")
	{dimnames(covariance) <- list(c(names.output, "rho","omega1"), c(names.output, 	"rho","omega1"))
	dimnames(correlation) <- list(c(names.output, "rho","omega1"), c(names.output, 	"rho","omega1"))}
	else  if (dependence=="indR2")
	{dimnames(covariance) <- list(c(names.output, "omega1","omega2"), c(names.output, 	"omega1","omega2"))
	dimnames(correlation) <- list(c(names.output, "omega1","omega2"), c(names.output, 	"omega1","omega2"))}
	else  if (dependence=="AR1R2")
	{dimnames(covariance) <- list(c(names.output, "rho","omega1","omega2"), c(names.output, 	"rho","omega1","omega2"))
	dimnames(correlation) <- list(c(names.output, "rho","omega1","omega2"), c(names.output, 	"rho","omega1","omega2"))}
	

#### To compute fitted values
	Fitted <- rep(NA, n.tot)
 	if (dependence=="ind"|dependence=="AR1")
	{Fitted[id.not.na] <- X %*% coefficients[1:(p - 1)]	 
	bi.estimate<-matrix(NaN, ncol = 1)}
	else  if (dependence=="indR")
	{aux <- LogL.pss0I.aux (parameters=coefficients, X=X, data=data2, trace=trace)
	 Fitted<- aux$fit 
	 bi.estimate<- aux$bi.est
	 bi.estimate<-matrix(bi.estimate, ncol = 1)
	 colnames(bi.estimate)<-names.Z}
	else  if (dependence=="AR1R")
	{aux<- LogL.pss1I.aux (parameters=coefficients, X=X, data=data2, trace=trace)
	 Fitted<- aux$fit 
	 bi.estimate<- aux$bi.est
	 bi.estimate<-matrix(bi.estimate, ncol = 1)
	 colnames(bi.estimate)<-names.Z}
	else  if (dependence=="indR2")
	{aux<- LogL.pss0Ic2.aux (parameters=coefficients, X=X ,Z=Z, data=data2, trace=trace)
	 Fitted<- aux$fit 
	 bi.estimate<- aux$bi.est
   bi.estimate<-matrix(bi.estimate, ncol = 2)
   colnames(bi.estimate)<-names.Z}
	else  if (dependence=="AR1R2")
	{aux<- LogL.pss1Ic2.aux (parameters=coefficients, X=X, Z=Z, data=data2, trace=trace)
	Fitted<- aux$fit 
	bi.estimate<- aux$bi.est
	bi.estimate<-matrix(bi.estimate, ncol = 2)
	colnames(bi.estimate)<-names.Z}

	ncoef<-length(coefficients)
	aic<-(2*temp$value+2*ncoef)
	y<-data2[[2]]

	Fitted <- exp(Fitted)
	Fitted[is.na(y)] <- NA

	y.matrix<-matrix(y,ncol=n.time,byrow=TRUE)
 	y.av<-apply(y.matrix,2,mean,na.rm=TRUE)
	Fitted.matrix<-matrix(Fitted,ncol=n.time,byrow=TRUE)
	Fitted.av<-apply(Fitted.matrix,2,mean,na.rm=TRUE)


cl<- new("cold", coefficients = coefficients, se = se, covariance =covariance, correlation=correlation,
	log.likelihood=- temp$value, message = temp$convergence, n.cases=n.cases, ni.cases=ni.cases, aic=aic,
      	Fitted=Fitted, bi.estimate=bi.estimate,Fitted.av=Fitted.av, Time=Time, 
	      model.matrix=X,y.matrix=y.matrix, random.matrix=Z, subset.data=subset.data, final.data=final.data,
	      y.av=y.av, data.id=id, call=call)

}


