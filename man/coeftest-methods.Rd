\name{coeftest-methods}
\docType{methods}
\alias{coeftest-methods}
\alias{coeftest,cold-method}
\title{Methods for function \code{coeftest}}
\description{Extract information from poisson regression model objects of class \code{\link{cold}}.}
\usage{
\S4method{coeftest}{cold}(object)}
\arguments{
\item{object}{an object of class \code{\link[=cold-class]{cold}}.}}
\value{Extract a list of summary statistics from poisson regression model corresponding to the fixed effects coefficients.}

\section{Methods}{
\describe{
\item{\code{signature(object="cold")}:}{list of summary statistics of fixed effects coefficients for \code{\link{cold}} object.} 
}
}

\examples{
#####  data= seizure

seiz1M<-cold(y~lage+lbase+v4+trt+trt:lbase, data=seizure, start=NULL, 
dependence="AR1")

coeftest(seiz1M)
}

\keyword{methods}