setMethod("getcoef",
          signature(object = "cold"),
          function (object) 
          {
            return(object@coefficients)
          }
        )