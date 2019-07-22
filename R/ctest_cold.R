setMethod(f="coeftest", signature(object = "cold"), 
          function (object) 
          {
            coef <- object@coefficients
            nas <- is.na(coef[, 1])
            cnames <- names(coef[, 1][!nas])
            coef <- matrix(rep(coef[, 1][!nas], 4), ncol = 4)
            coef.aux<-matrix(rep(coef[, 1][!nas], 4), ncol = 4)
            coef[, 2] <- object@se[, 1][!nas]
            coef[, 3] <- round(coef[, 1]/coef[, 2], 3)
            coef.aux[,1]<-pnorm(coef[, 1]/coef[, 2])
            coef.aux[,2]<-1-pnorm(coef[, 1]/coef[, 2])
            for (i in 1:dim(coef)[[1]])
            {coef.aux[i,3]<-min(coef.aux[i, 1],coef.aux[i, 2])}
            coef[, 4] <- round(2*coef.aux[,3],6)
            dimnames(coef) <- list(cnames, c("Estimate", "Std. Error", "z value", "p-value"))
            
            # for dependence="ind" 
            if(all(is.na(match(cnames, "omega1"))) && all(is.na(match(cnames, "rho"))) )
            {
            #  cat("\nz test of coefficients:\t\n") 
              return(coef[ ,  ])
            }
            
            # for dependence="AR1" 
            else if(all(is.na(match(cnames, "omega1"))) && !all(is.na(match(cnames, "rho"))) )
            {
             # cat("\nz test of coefficients:\t\n") 
              return(coef[1:dim(coef)[[1]],  ])
            }
            
            # for dependence="indR" 
            else if (!all(is.na(match(cnames, "omega1"))) &&  
                     all(is.na(match(cnames, "omega2"))) && 
                     all(is.na(match(cnames, "rho"))) )
            {
            #  cat("\nz test of coefficients:\t\n") 
              return(coef[-dim(coef)[[1]],  ])
            }
            
            # for dependence="AR1R" 
            else if (!all(is.na(match(cnames, "omega1"))) &&  
                     all(is.na(match(cnames, "omega2"))) && 
                     !all(is.na(match(cnames, "rho"))) )
            {
            #  cat("\nz test of coefficients:\t\n") 
              return(coef[1:(dim(coef)[[1]]-1),  ])
            }
            
            # for dependence="indR2" 
            else if (!all(is.na(match(cnames, "omega2"))) && 
                     all(is.na(match(cnames, "rho"))))
            {
            #  cat("\nz test of coefficients:\t\n") 
              return(coef[1:(dim(coef)[[1]]-2),  ])
            }
            
            # for dependence="AR1R2" 
            else if (!all(is.na(match(cnames, "omega2"))) && 
                     !all(is.na(match(cnames, "rho"))) )
            {
            #  cat("\nz test of coefficients:\t\n") 
              return(coef[1:(dim(coef)[[1]]-2),  ])
            }
            
            
          }
        )