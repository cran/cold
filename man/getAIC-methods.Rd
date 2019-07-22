\name{getAIC-methods}
\docType{methods}
\alias{getAIC-methods}
\alias{getAIC,cold-method}
\title{Methods for function \code{getAIC}}
\description{Methods for function \code{getAIC} extracting the Akaike information criterion of the fitted model object from class \code{\link{cold}}. }
\usage{
\S4method{getAIC}{cold}(object)
}
\arguments{
\item{object}{an object of class \code{\link[=cold-class]{cold}}.}}
\value{Returns a numeric value corresponding to the AIC of the fitted model.}
\section{Methods}{
\describe{ 
    \item{getAIC}{\code{signature(object="cold")}: Returns a numeric value corresponding to the AIC of the fitted model.}   
}}

\examples{
#####  data = seizure

### AR1
seiz1M<-cold(y~lage+lbase+v4+trt+trt:lbase, data=seizure, start=NULL, 
dependence="AR1")

getAIC(seiz1M)

### indR
seiz0R<-cold(y ~ lage + lbase + trt + trt:lbase + v4, random = ~ 1,
        data = seizure, dependence = "indR")

getAIC(seiz0R)

}

\keyword{methods}

