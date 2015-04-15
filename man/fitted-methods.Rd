\name{fitted-methods}
\docType{methods}
\alias{fitted-methods}
\alias{fitted,cold-method}
\title{Extract "cold" Fitted Values}
\description{Methods for function \code{fitted} extracting fitted values of a fitted model object from class \code{\link[=cold-class]{cold}}. }

\usage{
\S4method{fitted}{cold}(object)
}
\arguments{ 
  \item{object}{an object of class \code{\link[=cold-class]{cold}}.}
    }  

\section{Methods}{
\describe{
\item{\code{signature(object="cold")}:}{fitted for \code{\link{cold}} object.} 
}}

\examples{
#####  data= seizure

### without patient 207
seizure207<-seizure[seizure$id!=207,] 

### indR
seiz0R<-cold(y~lage+lbase+v4+trt+trt:lbase, data=seizure207, 
dependence="indR")

str(fitted(seiz0R))

#####  data= datacold
### AR1
modI<- cold(z~Time*Treatment, data=datacold, time="Time", 
id="Subject", aggregate=Treatment, dependence="AR1")

str(fitted(modI))
}

\keyword{methods}

