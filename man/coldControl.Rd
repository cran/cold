\name{coldControl}
\alias{coldControl}
\title{Auxiliary for Controlling "cold" Fitting}
\description{Auxiliary function as user interface for \code{cold} fitting }
\usage{coldControl(maxit=100, abstol=-Inf, reltol=sqrt(.Machine$double.eps))}    
\arguments{ 
 \item{maxit}{maximum number of iterations.}      
 \item{abstol}{absolute convergence tolerance.}     
 \item{reltol}{relative convergence tolerance.}
 }      
 \details{          
   See R documentation of \code{\link[stats]{optim}} for details of standard default values 
   for the remaining options not considered in \code{coldControl}.              
  }    
   
\value{A list with the arguments as components.}    

\seealso{\code{\link{cold-class}}, \code{\link[stats]{optim}}}

\keyword{function}
