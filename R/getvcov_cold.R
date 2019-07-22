setMethod("getvcov",
    signature(object = "cold"),
    function (object) 
    {
      cov<-object@covariance
      cnames <- rownames(cov)
      r1<-nrow(cov)
      c1<-ncol(cov)
 
      # for dependence="ind" or dependence="AR1"
       if(all(is.na(match(cnames, "omega1")))  )
        {     cov
          }
      
      
      # for dependence="indR" or dependence="AR1R"
      else if (!all(is.na(match(cnames, "omega1"))) &&  
          all(is.na(match(cnames, "omega2"))) )
        {  cov.aux<-cov[1:(r1-1),1:(c1-1)]
           cov.aux
        }
        
      # for dependence="indR2" or dependence="AR1R2"
      else if (!all(is.na(match(cnames, "omega2"))) )
      { cov.aux<-cov[1:(r1-2),1:(c1-2)]
        cov.aux
      }
      
} )