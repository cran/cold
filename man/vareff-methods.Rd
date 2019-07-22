\name{vareff-methods}
\docType{methods}
\alias{vareff-methods}
\alias{vareff,cold-method}
\title{Methods for function \code{vareff}}
\description{ Methods for function \code{vareff} extracting the variance estimates of random effects of a fitted model object. }

\section{Methods}{
\describe{
\item{\code{signature(object="cold")}:}{\code{vareff} for \code{\link{cold}} object.} 
}}

\examples{

#####  data = seizure

### indR
seiz0R<-cold(y ~ lage + lbase + trt + trt:lbase + v4, random = ~ 1,
        data = seizure, dependence = "indR")

vareff(seiz0R)

}

\keyword{methods}

