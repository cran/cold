\name{getvcov-methods}
\docType{methods}
\alias{getvcov-methods}
\alias{getvcov,cold-method}
\title{Methods for function \code{getvcov}}
\description{Extract the variance-covariance matrix of a fitted model object from class \code{\link[=cold-class]{cold}}.}

\section{Methods}{
\describe{
    \item{\code{signature(object="cold")}:}{Returns a numeric value corresponding  to the variance-covariance matrix of the fixed effect estimates.}
}}

\examples{
#####  data = seizure

### AR1
seiz1M<-cold(y~lage+lbase+v4+trt+trt:lbase, data=seizure, start=NULL, 
dependence="AR1")

getvcov(seiz1M)


}

\keyword{methods}
