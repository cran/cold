\name{fixeff-methods}
\docType{methods}
\alias{fixeff-methods}
\alias{fixeff,cold-method}
\title{Methods for function \code{fixeff}}
\description{Methods for function \code{fixeff} extracting fixed effects estimates of a fitted model object from class \code{\link[=cold-class]{cold}}. }

\section{Methods}{
\describe{
\item{\code{signature(object="cold")}:}{ fixed effects estimates of a fitted model object.} 
}}

\examples{
#####  data = seizure

### AR1
seiz1M<-cold(y~lage+lbase+v4+trt+trt:lbase, data=seizure, start=NULL, 
dependence="AR1")

fixeff(seiz1M)

### indR
seiz0R<-cold(y ~ lage + lbase + trt + trt:lbase + v4, random = ~ 1,
        data = seizure, dependence = "indR")

fixeff(seiz0R)
}

\keyword{methods}

