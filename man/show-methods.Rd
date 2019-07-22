\name{show-methods}
\docType{methods}
\alias{show-methods}
\alias{show,cold-method}
\alias{show,summary.cold-method}
\title{Methods for function \code{show} }
\description{Show objects of class \code{cold} and \code{summary.cold}.}
\section{Methods}{
\describe{
\item{\code{signature(object = "cold")}}{Print simple summary of a \code{cold} object, just the call, the number of profiles in the fit, 
the number of coefficients, the value of the log-likelihood and a message giving additional information returned by the optimizer. }

\item{\code{signature(object = "summary.cold")}}{Shows call, the number of profiles in the fit, 
table of coefficients, standard errors and p-values, the log-likelihood, the AIC coefficient, 
and a message giving additional information returned by the optimizer.}
}}
\keyword{methods}

