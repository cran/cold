\name{cold}
\alias{cold}
\title{  Fit of parametric models via likelihood method}
\description{Performs the fit of parametric models via likelihood method. 
Serial dependence and two random effects are allowed according to the stochastic model chosen. 
Missing values are automatically accounted for computing the likelihood function.}
\usage{cold(formula, random, data, id, time, subSET, dependence, 
    start, method, integration, M, control, integrate, cublim,
    aggregate, trace)}
\arguments{ 
  \item{formula}{a description of the model to be fitted of the form response~predictors}.
    \item{random}{the predictos that includes random effects of the form response~predictors}.
  \item{data}{a \code{data} frame containing the variables in the formula. NA values are allowed. 
  If data is missing, an error message is produced. See "Details".}
  \item{id}{a string that matches the name of the \code{id} variable in \code{data}, 
  i.e., the subject identification variable. By default, the program expects a variable named \code{id} 
   to be present in the \code{data.frame}, otherwise the name of the variable playing the role of \code{id} must be declared by assigning \code{id} here. }
  \item{time}{a string that matches the name of the \code{time} variable in data. By default, the program expects a variable named \code{time} 
   to be present in the \code{data.frame}, otherwise the name of the variable playing the role of time must be declared by assigning \code{time} here. }
  \item{subSET}{an optional expression indicating the subset of the rows of \code{data} that should be used in the fit. All observations are included by default.}
  \item{dependence}{expression stating which \code{dependence} structure should be used in the fit. The default is "ind".
  According to the stochastic model chosen serial dependence and random effects are allowed. 
  There are six options: "\code{ind}" (independence), "\code{AR1}" (first order autoregressive),  
  "\code{indR}" (independence with random intercept), "\code{AR1R}" (first order autoregressive with random intercept), "\code{indR2}" (independence with two random effects) or "\code{AR1R2}" (first order autoregressive with two random effects).}
   \item{start}{a vector of initial values for the nuisance parameters of the likelihood. The dimension of the vector is according to the structure of the dependence model.}
  \item{method}{The \code{method} used in the optimization process: "\code{BFGS}", "\code{CG}", "\code{L-BFGS-B}" and "\code{SANN}". The default is "\code{BFGS}". 
    See \code{\link[stats]{optim}} for details.}
  \item{integration}{The \code{integration} code allows the user choose the integration method to solve the integrals: "\code{QUADPACK}" (fortran routines, only for random intercept models), "\code{cubature}" (uses cubature package to compute integrals when the dependence structure includes one or two random effects), "\code{MC}" (uses Monte Carlo methods to compute integrals when the dependence structure includes one or two random effects). The default is "\code{QUADPACK}".}
  \item{M}{Number of random points considered to evaluate the integral when the user choose Monte Carlo methods ("\code{integration}=MC"). The default is set to 6000.}
  \item{control}{a list of algorithmic constants for the optimizer \code{optim}. 
  See R documentation of \code{optim.control} for details and possible control options. By default, \code{cold} sets the maximum number 
  of iterations (\code{maxit}) equal to 100, the absolute convergence tolerance (\code{abstol}) and the relative 
  convergence tolerance (\code{rel.tol}) equal to 1e-6 and uses the \code{\link[stats]{optim}} standard default values for the remaining options.}
  \item{integrate}{a list of algorithmic constants for the computation of a definite integral using a Fortran-77 subroutine. See "Details".}      
 \item{cublim}{a list of algorithmic constants for the computation of a definite integral when the integration argument is set to cubature.}   
   \item{aggregate}{a string that permits the user identify the factor to be used in \code{\link{plot-methods}}.}
   \item{trace}{logical flag: if TRUE, details of the nonlinear optimization are printed. By default the flag is set to FALSE.}
}  

\section{Background}{ 
   Assume that each subject of a given set has been observed at number of 
   successive time points. For each subject and for each time point, a count 
   response variable, and a set of covariates are recorded.  

   Individual random effects, \eqn{b_0}, can be incorporated in the form 
   of a random intercept term of the linear predictor of the logarithmic regression, assuming a normal distribution of mean 0 and variance 
   \eqn{\sigma^2}, parameterized as \eqn{\omega=\log(\sigma^2)}. 
   The combination of serial Markov \code{dependence} with a random intercept corresponds here to the \code{dependence} structures \code{indR}, \code{AR1R}.
  
    Two dimensional randoms effects can also be incorporated the linear predictor of the logarithmic regression. Consider a two-dimensional vector of random effects \eqn{b=(b_0,b_1)} where we assumed to be a random sample from the bivariate normal distribution, \eqn{b\sim N(0,D)} with \eqn{var(b_0)=\sigma^2_{b_0}}, \eqn{var(b_1)=\sigma^2_{b_1}} and \eqn{cov(b_0,b_1)=0}.

The combination of serial Markov \code{dependence} with two random effects corresponds here to the \code{dependence} structures \code{indR2}, \code{AR1R2}.


   Original sources of the above formulation are given by Azzalini (1994), as
   for the AR1, and by \enc{Gonçalves}{Goncalves} (2002) and 
   \enc{Gonçalves}{Goncalves} and Azzalini (2008) for the its extensions.
   }

\details{
\code{data} are contained in a \code{data.frame}. Each element of the \code{data} argument must be identifiable by a name. The simplest situation occurs when all subjects are observed at the same time points. If there are missing values in the response variable \code{NA} values must be inserted.
The response variable represent the individual profiles of each subject, it is expected a variable in the \code{data.frame} that identifies the correspondence of each component of the response variable to the subject that it belongs, 
by default is named \code{id} variable. It is expected a variable named \code{time} to be present in the \code{data.frame}. 
If the \code{time} component has been given a different name, this should be declared. 
The \code{time} variable should identify the time points that each individual profile has been observed.
 

 If a response profile is replicated several times, a variable called 
 \code{counts} must be created accordingly. 
 This vector is used for weighting the response profile indicating for each individual profile the number of times that is replicated. 
 The vector \code{counts} must repeat the number of the observed replications of each individual profile as many times as the number of observed time 
 points for the correspondent profile. The program expect such vector to be named \code{counts}. 
 If each profile has been observed only once, the construction of the vector \code{counts} is not required. 

\code{subSET} is an optional expression indicating the subset of \code{data} that should be 
    used in the fit. This is a logical statement of the type 
    \code{variable 1} == "a" & \code{variable 2} > x  
  which identifies the observations to be selected. All observations are included by default.    

For the models with random intercept \code{indR} and \code{AR1R}, by default 
\code{cold} compute integrals based on a Fortran-77 subroutine package 
\code{QUADPACK}. For some data sets, when the \code{dependence} structure has 
a random intercept term, the user could have the need to do a specification 
of the \code{integrate} argument list changing  
the integration limits in the \code{coldIntegrate} function. 
The \code{\link{coldIntegrate}} is an auxiliary function for controlling \code{cold} fitting. There are more two options to fit models with a random intercept by setting \code{integration="cubature"} or \code{integration="MC"}. 
For the models with two random effects \code{indR2} and \code{AR1R2}, the user has two define the \code{integration} method by setting \code{integration="cubature"} or \code{integration="MC"}. The second random effect is 
considered to be included in the \code{time} argument that plays the role of the time variable in the \code{data.frame}. For the two random effects models we have a random intercept and a random slope.
}

\value{An object of class \code{\link[=cold-class]{cold}}.}

\references{Azzalini, A. (1994). Logistic regression and other discrete data models for serially correlated observations.  
\emph{J. Ital. Stat. Society}, 3 (2), 169-179.

\enc{Gonçalves}{Goncalves}, M. Helena (2002). \emph{ Likelihood methods for discrete longitudinal data}. 
PhD thesis, Faculty of Sciences, University of Lisbon.
  
\enc{Gonçalves}{Goncalves}, M. Helena, Cabral, M. \enc{Salomé}{Salome}, Ruiz de Villa, M. Carme, Escrich, Eduardo and Solanas, Montse. (2007). 
Likelihood approach for count data in longitudinal experiments. \emph{Computational statistics and Data Analysis}, 51, 12, 6511-6520.}
  
\author{M. Helena \enc{Gonçalves}{Goncalves} and M. \enc{Salomé}{Salome} Cabral }

\seealso{\code{\link{cold-class}}, \code{\link{coldControl}}, \code{\link{coldIntegrate}}, \code{\link[stats]{optim}}}

\examples{
#####  data= seizure
str(seizure)

### AR1
seiz1M<-cold(y~lage+lbase+v4+trt+trt:lbase, data=seizure, start=NULL, 
dependence="AR1")

summary(seiz1M)
getAIC(seiz1M)
getLogLik(seiz1M)

### independence
seiz0M<-cold(y~lage+lbase+v4+trt+trt:lbase, data=seizure, start=NULL, 
dependence="ind")

summary(seiz0M)
getAIC(seiz0M)
getLogLik(seiz0M)


#####  data= datacold
str(datacold)

### AR1R with the default integration method
mod1R<- cold(z~Time*Treatment, random=~1, data=datacold, time="Time", 
id="Subject", aggregate=Treatment, dependence="AR1R")

summary (mod1R)
vareff(mod1R)
randeff(mod1R)

### AR1R with integration="cubature"
\donttest{ 
mod1R.c<- cold(z~Time*Treatment, random=~1, data=datacold, time="Time", 
id="Subject", aggregate=Treatment, dependence="AR1R", 
integration="cubature")
summary (mod1R.c)
}

}

\keyword{function}
