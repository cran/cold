\name{model.mat-methods}
\docType{methods}
\alias{model.mat-methods}
\alias{model.mat,cold-method}
\title{Methods for function \code{model.mat}}
\description{Methods for function \code{model.mat} extracting the fixed effects model matrix for a fitted model object from class \code{\link[=cold-class]{cold}}.}

\section{Methods}{
\describe{
    \item{\code{signature(object="cold")}:}{Returns the fixed effects model matrix of the fitted model.}
}}

\examples{
#####  data = seizure

### AR1
seiz1M<-cold(y~lage+lbase+v4+trt+trt:lbase, data=seizure, start=NULL, 
dependence="AR1")

model.mat(seiz1M)[1:20,]


}
\keyword{methods}
