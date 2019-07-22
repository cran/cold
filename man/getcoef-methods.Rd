\name{getcoef-methods}
\docType{methods}
\alias{getcoef-methods}
\alias{getcoef,cold-method}
\title{Methods for function \code{getcoef}}
\description{Methods for function \code{getcoef} extracting the coefficient estimates of the fitted model object from class \code{\link{cold}}. }

\section{Methods}{
\describe{
    \item{\code{signature(object="cold")}:}{Returns the coefficient estimates of the fitted model.}
}}

\examples{
#####  data = seizure

### AR1
seiz1M<-cold(y~lage+lbase+v4+trt+trt:lbase, data=seizure, start=NULL, 
dependence="AR1")

getcoef(seiz1M)

### indR
seiz0R<-cold(y ~ lage + lbase + trt + trt:lbase + v4, random = ~ 1,
        data = seizure, dependence = "indR")

getcoef(seiz0R)

}

\keyword{methods}
