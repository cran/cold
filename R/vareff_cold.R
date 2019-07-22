setMethod("vareff",
          signature(object = "cold"),
          function (object) 
          {
            coef <- object@coefficients
            nas <- is.na(coef[, 1])
            cnames <- names(coef[, 1][!nas])
    
            # for dependence="indR" 
            if (!all(is.na(match(cnames, "omega1"))) &&  
                     all(is.na(match(cnames, "omega2"))) && 
                     all(is.na(match(cnames, "rho"))) )
            {
              coef[dim(coef)[[1]], 1] <- exp(coef[dim(coef)[[1]], 1])
             # cat("\nRandom effects:\t\n") 
              table.r1<-data.frame(coef[dim(coef)[[1]], 1])
              dimnames(table.r1)<-list("(Intercept)", "Variance")
              return(table.r1)
            }
            

            # for dependence="AR1R" 
            else if (!all(is.na(match(cnames, "omega1"))) &&  
                     all(is.na(match(cnames, "omega2"))) && 
                     !all(is.na(match(cnames, "rho"))) )
            {
         
              ### for random coef  
              coef[dim(coef)[[1]], 1] <- exp(coef[dim(coef)[[1]], 1])
             # cat("\nRandom effects:\t\n") 
              table.r1<-data.frame(coef[dim(coef)[[1]], 1])
              dimnames(table.r1)<-list("(Intercept)", "Variance")
              return(table.r1)
            }
            
            # for dependence="indR2" 
            else if (!all(is.na(match(cnames, "omega2"))) && 
                     all(is.na(match(cnames, "rho"))))
            {
              ### for random coef
              coef[dim(coef)[[1]]-1, 1] <- exp(coef[dim(coef)[[1]]-1, 1])
              coef[dim(coef)[[1]], 1] <- exp(coef[dim(coef)[[1]], 1])
              
            #  cat("\nRandom effects:\t\n") 
              table.r2<-data.frame(coef[(dim(coef)[[1]]-1):dim(coef)[[1]], 1])
              dimnames(table.r2)<-list(c("(Intercept)","Time"), "Variance")
              return(table.r2)
            }
            
            # for dependence="AR1R2" 
            else if (!all(is.na(match(cnames, "omega2"))) && 
                     !all(is.na(match(cnames, "rho"))) )
            {
              ### for random coef
              coef[dim(coef)[[1]]-1, 1] <- exp(coef[dim(coef)[[1]]-1, 1])
              coef[dim(coef)[[1]], 1] <- exp(coef[dim(coef)[[1]], 1])
              
            #  cat("\nRandom effects:\t\n") 
              table.r2<-data.frame(coef[(dim(coef)[[1]]-1):dim(coef)[[1]], 1])
              dimnames(table.r2)<-list(c("(Intercept)","Time"), "Variance")
              return(table.r2)
            }
            
            else 
              warning("\nOnly to Random effects model") 
              
          }
        )