\name{fitted-methods}
\docType{methods}
\alias{fitted-methods}
\alias{fitted,cold-method}
\title{Methods for function \code{fitted}}
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

#####  data = seizure

seiz1M<-cold(y~lage+lbase+v4+trt+trt:lbase, data=seizure, start=NULL, 
dependence="AR1")

fitted(seiz1M)[1:16]

}

\keyword{methods}

