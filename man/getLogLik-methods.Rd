\name{getLogLik-methods}
\docType{methods}
\alias{getLogLik-methods}
\alias{getLogLik,cold-method}
\title{Extract Log-Likelihood}
\description{Methods for function \code{getLogLik} extracting the Log-Likelihood 
for one fitted model object from class \code{\link[=cold-class]{cold}}.}
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
\keyword{methods}