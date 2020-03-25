
setMethod(f="fitted", signature=c(object="cold"), 
  function (object) { 
    
    fit.i<-object@Fitted
    names(fit.i) <- c(1:length(fit.i))
    round(fit.i ,6)                  
    }
)