\name{resid-methods}
\docType{methods}
\alias{resid-methods}
\alias{resid,cold-method}
\title{Extract "cold" Residuals}
\description{ The residuals for the fitted model represented by \code{object} are extracted.}

\usage{
\S4method{resid}{cold}(object, type = c( "pearson","response","null"),...)
}
\arguments{ 
  \item{object}{an object of class \code{\link[=cold-class]{cold}}.}
  \item{type}{ an optional character string specifying the type of residuals to be used. 
  Two types are allowed: pearson and response. Defaults to "pearson".}
  \item{...}{other arguments.}
   }  

\section{Methods}{
\describe{
\item{\code{signature(object="cold")}:}{residuals for \code{\link{cold}} object.} 
}}

\examples{
#####  data= seizure
### without patient 207
seizure207<-seizure[seizure$id!=207,] 

### indR
seiz0R<-cold(y~lage+lbase+v4+trt+trt:lbase, data=seizure207, 
dependence="indR")

str(resid(seiz0R))

#####  data= datacold
### AR1
modI<- cold(z~Time*Treatment, data=datacold, time="Time", 
id="Subject", aggregate=Treatment, dependence="AR1")

str(resid(modI))
}

\keyword{methods}

