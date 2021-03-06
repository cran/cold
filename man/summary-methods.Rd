\name{summary-methods}
\docType{methods}
\alias{summary-methods}
\alias{summary,cold-method}
\title{Methods for function \code{summary}}
\description{Summarize objects.}
\usage{
\S4method{summary}{cold}(object, cov=FALSE, cor=FALSE)}
\arguments{ 
  \item{object}{an object of class \code{\link[=cold-class]{cold}}.}
  \item{cov}{if set to TRUE prints the matrix of covariances between parameters estimates. The default is FALSE.}
  \item{cor}{if set to TRUE prints the matrix of correlations between parameters estimates. The default is FALSE.}  
   }  
\details{ Computes and returns a list of summary statistics of the fitted model given a \code{cold} object, using the components
(list elements) "call" and "terms" from its argument, plus 
depending on the structure of the dependence model chosen. In the table for coefficient estimates will appear \code{rho} if the dependence structure of the process corresponds to an \code{AR1}, \code{AR1R} or \code{AR1R2}.
If the structure of the dependence model chosen includes the random intercept (models "\code{indR}" and "\code{AR1R}"), 
the variance estimate for random intercept is given. 
If the structure of the dependence model chosen includes two random effects (models "\code{indR2}"  and "\code{AR1R2}") 
the variance estimates for random intercept and for the slope (time) are given. 
}

\section{Methods}{
\describe{
\item{\code{signature(object = "ANY")}:}{Generic function.}
\item{\code{signature(object = "cold")}:}{Prints a summary as an object of 
class \code{\link[=summary.cold-class]{summary.cold}}, 
containing information about the matched call to \code{cold}, the number of profiles in the data, 
the number of profiles used in the fit, the log-likelihood, the AIC,
a table with estimates, asymptotic SE, t-values and p-values,  
the estimated correlation and variance-covariance matrix for the estimated parameters if the user wishes, 
and a message giving additional information returned by the optimizer.}
}}
\keyword{methods}

