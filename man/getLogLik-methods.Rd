\name{getLogLik-methods}
\docType{methods}
\alias{getLogLik-methods}
\alias{getLogLik,cold-method}
\title{Methods for function \code{getLogLik}}
\description{Extract the Log-Likelihood of a fitted model object from class \code{\link[=cold-class]{cold}}.}
\usage{
\S4method{getLogLik}{cold}(object)
}
\arguments{
\item{object}{an object of class \code{\link[=cold-class]{cold}}.}}
\value{Returns a numeric value corresponding to the log-Likelihood of the fitted model.}
\section{Methods}{
\describe{
    \item{\code{signature(object="cold")}:}{Returns a numeric value corresponding to the log-Likelihood of the fitted model.}
}}

\examples{
#####  data = seizure

### AR1R
seiz1M<-cold(y~lage+lbase+v4+trt+trt:lbase, data=seizure, start=NULL, 
dependence="AR1")

getLogLik(seiz1M)

### indR
seiz0R<-cold(y ~ lage + lbase + trt + trt:lbase + v4, random = ~ 1,
        data = seizure, dependence = "indR")

getLogLik(seiz0R)

}

\keyword{methods}
