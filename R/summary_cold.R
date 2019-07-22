setMethod("summary",
    signature(object = "cold"),
    function(object,cov=FALSE,cor=FALSE)
	{    
  cat("\nCall:\n")
  print(object@call)
	cat("\nNumber of profiles in the dataset: ", object@ni.cases)
	cat("\nNumber of profiles used in the fit: ", object@n.cases)
	cat("\nDependence structure:",object@call$dependence)
	cat("\nLog likelihood: ", round(object@log.likelihood, 4))
	cat("\nAIC: ", round(object@aic, 4),"\n")
  
	coef <- object@coefficients
	nas <- is.na(coef[, 1])
	cnames <- names(coef[, 1][!nas])
	coef <- matrix(rep(coef[, 1][!nas], 4), ncol = 4)
	coef.aux<-matrix(rep(coef[, 1][!nas], 4), ncol = 4)
#	coef[, 1] <- 1:dim(coef)[[1]]
	coef[, 2] <- object@se[, 1][!nas]
	coef[, 3] <- round(coef[, 1]/coef[, 2], 3)
	coef.aux[,1]<-pnorm(coef[, 1]/coef[, 2])
	coef.aux[,2]<-1-pnorm(coef[, 1]/coef[, 2])
	for (i in 1:dim(coef)[[1]])
	{coef.aux[i,3]<-min(coef.aux[i, 1],coef.aux[i, 2])}
	coef[, 4] <- round(2*coef.aux[,3],6)
  dimnames(coef) <- list(cnames, c("Estimate", "Std. Error", "z value", "p-value"))
	
	Z.aux<-object@random.matrix
	names.Z <- dimnames(Z.aux)[[2]]  #new for random
	
	# for dependence="ind" 
	if(all(is.na(match(cnames, "omega1"))) && all(is.na(match(cnames, "rho"))) )
	 {
	  cat("\nFixed effects:\t\n") 
	  print(coef[ ,  ])
	}
	
	# for dependence="AR1" 
	else if(all(is.na(match(cnames, "omega1"))) && !all(is.na(match(cnames, "rho"))) )
	{
	  cat("\nFixed effects:\t\n") 
	  print(coef[1:dim(coef)[[1]]-1 ,  ])
	  
	  cat("\nEstimated correlation parameter: rho \t\n") 
	  print(coef[dim(coef)[[1]],  ])
	}
	
	# for dependence="indR" 
	else if (!all(is.na(match(cnames, "omega1"))) &&  
	         all(is.na(match(cnames, "omega2"))) && 
	         all(is.na(match(cnames, "rho"))) )
	  {
	  cat("\nFixed effects:\t\n") 
	  print(coef[-dim(coef)[[1]],  ])

	  coef[dim(coef)[[1]], 1] <- exp(coef[dim(coef)[[1]], 1])
	  cat("\nRandom effects:\t\n") 
	  table.r1<-data.frame(coef[dim(coef)[[1]], 1])
	  dimnames(table.r1)<-list("(Intercept)", "Variance")
	  print(table.r1)
	}
	
	# for dependence="AR1R" 
	else if (!all(is.na(match(cnames, "omega1"))) &&  
	         all(is.na(match(cnames, "omega2"))) && 
	         !all(is.na(match(cnames, "rho"))) )
	{
	  cat("\nFixed effects:\t\n") 
	  print(coef[1:(dim(coef)[[1]]-2),  ])
	  
	  cat("\nEstimated correlation parameter: \t\n") 
		table.r0<-matrix(coef[dim(coef)[[1]]-1, ], nrow=1)
  	dimnames(table.r0)<-list("rho", c("Estimate", "Std. Error", "z value", "p-value"))
  	print(table.r0)
		
		
		### for random coef  
	  coef[dim(coef)[[1]], 1] <- exp(coef[dim(coef)[[1]], 1])
	  cat("\nRandom effects:\t\n") 
	  table.r1<-data.frame(coef[dim(coef)[[1]], 1])
	  dimnames(table.r1)<-list("(Intercept)", "Variance")
	  print(table.r1)
	}
	
	# for dependence="indR2" 
	else if (!all(is.na(match(cnames, "omega2"))) && 
	         all(is.na(match(cnames, "rho"))))
	 {
	  cat("\nFixed effects:\t\n") 
	  print(coef[1:(dim(coef)[[1]]-2),  ])
	  
	  coef[dim(coef)[[1]]-1, 1] <- exp(coef[dim(coef)[[1]]-1, 1])
	  coef[dim(coef)[[1]], 1] <- exp(coef[dim(coef)[[1]], 1])
	 
	  cat("\nRandom effects:\t\n") 
	  table.r2<-data.frame(coef[(dim(coef)[[1]]-1):dim(coef)[[1]], 1])
	  dimnames(table.r2)<-list(c("(Intercept)",names.Z[2]), "Variance")
	  print(table.r2)
	  }
	  
	# for dependence="AR1R2" 
	else if (!all(is.na(match(cnames, "omega2"))) && 
	         !all(is.na(match(cnames, "rho"))) )
	{
	  cat("\nFixed effects:\t\n") 
	  print(coef[1:(dim(coef)[[1]]-3),  ])
	  
	  cat("\nEstimated correlation parameter: \t\n") 
	  table.r0<-matrix(coef[dim(coef)[[1]]-2, ],nrow=1)
	  dimnames(table.r0)<-list("rho", c("Estimate", "Std. Error", "z value", "p-value"))
	  print(table.r0)
	  
	  ### for random coef
	  coef[dim(coef)[[1]]-1, 1] <- exp(coef[dim(coef)[[1]]-1, 1])
	  coef[dim(coef)[[1]], 1] <- exp(coef[dim(coef)[[1]], 1])
	  
	  cat("\nRandom effects:\t\n") 
	  table.r2<-data.frame(coef[(dim(coef)[[1]]-1):dim(coef)[[1]], 1])
	  dimnames(table.r2)<-list(c("(Intercept)",names.Z[2]), "Variance")
	  print(table.r2)
	}
	
	  if (cov)
	 {
	  cat("\nCovariance of Coefficients: \n")
		print(object@covariance, digits = 2)}
		if (cor){
		cat("\nCorrelation of Coefficients: \n")
		print(object@correlation, digits = 2)}
	  cat("\nMessage: ", object@message,"\n")
	})

