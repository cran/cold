\name{randeff-methods}
\docType{methods}
\alias{randeff-methods}
\alias{randeff,cold-method}
\title{Methods for function \code{randeff}}
\description{Methods for function \code{randeff} extracting conditional random effects of a fitted model object from class \code{\link[=cold-class]{cold}}. }

\section{Methods}{
\describe{
\item{\code{signature(object="cold")}:}{\code{randeff} for \code{\link{cold}} object.} 
}}

\examples{
#####  data = seizure

### indR
seiz0R<-cold(y ~ lage + lbase + trt + trt:lbase + v4, random = ~ 1,
        data = seizure, dependence = "indR")

randeff(seiz0R)

}


\keyword{methods}

