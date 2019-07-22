\name{coldIntegrate}
\alias{coldIntegrate}
\title{ Auxiliary for controlling "cold" fitting}
\description{Auxiliary function as user interface for \code{cold} fitting.}
\usage{coldIntegrate(li=-4,ls=4, epsabs=.Machine$double.eps^.25, 
epsrel=.Machine$double.eps^.25,limit=100,key=6,lig=-4,lsg=4) }
\arguments{ 
  \item{li}{lower limit of integration for the log-likelihood.}                                                                 
  \item{ls}{upper limit of integration for the log-likelihood.}                                                                 
   \item{epsabs}{absolute accuracy requested.}                                                                                  
   \item{epsrel}{relative accuracy requested.}                                                                                 
    \item{key}{integer from 1 to 6 for choice of local integration rule for number of Gauss-Kronrod quadrature points.          
    A gauss-kronrod pair is used with: \cr                                                                                             
    7 - 15 points if key = 1, \cr                                                                                                      
    10 - 21 points if key = 2,\cr                                                                                                      
    15 - 31 points if key = 3,\cr                                                                                                      
    20 - 41 points if key = 4,\cr                                                                                                      
    25 - 51 points if key = 5 and \cr                                                                                                  
    30 - 61 points if key = 6.}                                                                                                        
    \item{limit}{integer that gives an upperbound on the number of subintervals in the partition of (\code{li},\code{ls}), limit.ge.1.}                                                                                            
    \item{lig}{lower limit of integration for the gradient.}                                                                    
    \item{lsg}{upper limit of integration for the gradient.}                                                                    
   }      
 \details{       
 \code{coldIntegrate} returns a list of constants that are used to compute integrals based on a Fortran-77 subroutine \code{dqage} from a   Fortran-77 subroutine package \code{QUADPACK} for the numerical computation of definite one-dimensional integrals.              
 The subroutine \code{dqage} is a simple globally adaptive integrator in which it is possible to choose between 6 pairs          
 of Gauss-Kronrod quadrature formulae for the rule evaluation component. The source code \code{dqage} was modified and re-named  
 \code{dqager}, the change was the introduction of an extra variable that allow, in our Fortran-77 subroutines when              
 have a call to \code{dqager}, to control for which parameter the integral is computed.
 
For given values of \code{li} and \code{ls}, the above-described
numerical integration is performed over the interval 
(\code{li}*\eqn{\sigma}, \code{ls}*\eqn{\sigma}), where \eqn{\sigma=\exp(\omega)/2}
is  associated to the current parameter value \eqn{\omega} examined by
the \code{optim} function.  In some cases, this integration may
generate an error, and the user must suitably adjust the values of \code{li} 
and \code{ls}. In case different choices of these quantities all
lead to a successful run, it is recommended to retain the one with
largest value of the log-likelihood. Integration of the gradient is
regulated similarly by \code{lig} and \code{lsg}.
                         
For datasets where the individual profiles have a high number of 
observed time points (say, more than 30), 
use \code{coldIntegrate} function to set the integration limits for the 
likelihood and for the gradient to small values 
than the default ones. 
 
When the fitting procedure is complete but the computation of the information matrix 
produces NaNs, changing in \code{coldIntegrate} function the default values for the gradient 
integration limits (\code{lig} and \code{lsg}) might solve this problem.
  }    

\value{A list with the arguments as components.}    

\seealso{\code{\link{cold-class}}}

\examples{ 
#####  data= seizure

Integ <- coldIntegrate(li=-3.5, ls=3.5, lig=-3.5, lsg=3.5)

### AR1R without patient 207
\donttest{ 
seizure207 <- seizure[seizure$id!=207,]

seiz1R1.207 <- cold(y~lage+lbase+v4+trt+trt:lbase, random=~1, 
data=seizure207, dependence="AR1R", integrate=Integ)
summary (seiz1R1.207)
}

}

\keyword{function}
