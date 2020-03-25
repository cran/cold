\name{resid-methods}
\docType{methods}
\alias{resid-methods}
\alias{resid,cold-method}
\title{Methods for function \code{residd}}
\description{Methods for function \code{resid} extracting residual values of a fitted model object from class \code{\link[=cold-class]{cold}}. }

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

#####  data = seizure

seiz1M <- cold(y ~ lage + lbase + v4 + trt + trt:lbase, data = seizure, 
start = NULL, dependence = "AR1")

resid(seiz1M)[1:16]

}

\keyword{methods}

