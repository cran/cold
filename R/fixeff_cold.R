setMethod("fixeff",
          signature(object = "cold"),
          function (object) 
          {
 
            coef <- object@coefficients
            nas <- is.na(coef[, 1])
            cnames <- names(coef[, 1][!nas])
         
            
            # for dependence="ind" 
            if(all(is.na(match(cnames, "omega1"))) && all(is.na(match(cnames, "rho"))) )
            {  
              return(coef[ ,  ])
            }
            
            # for dependence="AR1" 
            else if(all(is.na(match(cnames, "omega1"))) && !all(is.na(match(cnames, "rho"))) )
            {
              return(coef[1:dim(coef)[[1]]-1 ,  ])
            }
            
            # for dependence="indR" 
            else if (!all(is.na(match(cnames, "omega1"))) &&  
                     all(is.na(match(cnames, "omega2"))) && 
                     all(is.na(match(cnames, "rho"))) )
            {
              return(coef[-dim(coef)[[1]],  ])
            }
            
            # for dependence="AR1R" 
            else if (!all(is.na(match(cnames, "omega1"))) &&  
                     all(is.na(match(cnames, "omega2"))) && 
                     !all(is.na(match(cnames, "rho"))) )
            {
              return(coef[1:(dim(coef)[[1]]-2),  ])
            }
            
            # for dependence="indR2" 
            else if (!all(is.na(match(cnames, "omega2"))) && 
                     all(is.na(match(cnames, "rho"))))
            {
              return(coef[1:(dim(coef)[[1]]-2),  ])
            }
            
            # for dependence="AR1R2" 
            else if (!all(is.na(match(cnames, "omega2"))) && 
                     !all(is.na(match(cnames, "rho"))) )
            {
              return(coef[1:(dim(coef)[[1]]-3),  ])
            } 
            
          }
        )