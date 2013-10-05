\name{getAIC-methods}
\docType{methods}
\alias{getAIC-methods}
\alias{getAIC,cold-method}
\title{Extract the Akaike Information Criterion}
\description{Methods for function \code{getAIC} extracting the Akaike information criterion 
for one fitted model object from class \code{\link[=cold-class]{cold}}. }
\usage{
\S4method{getAIC}{cold}(object)
}
\arguments{
\item{object}{an object of class \code{\link[=cold-class]{cold}}.}}
\value{Returns a numeric value corresponding to the AIC of the fitted model.}
\section{Methods}{
\describe{ 
    \item{getAIC}{\code{signature(object="cold")}:Returns a numeric value corresponding to the AIC of the fitted model.}   
}}
\keyword{methods}

