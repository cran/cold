\name{cold-package}
\alias{cold-package}
\docType{package}
\title{Count Longitudinal Data}
\description{
Performs Poisson regression analysis for longitudinal count data, allowing for serial dependence among observations from a given individual and two random effects. 
  Estimation is via maximization of the exact likelihood of a suitably defined model. Missing values and unbalanced data are allowed. 
}
\details{This package contains functions to perform the fit of parametric models via likelihood method for count longitudinal data using "S4" classes and methods as implemented in the \code{methods} package.}
\author{M. Helena \enc{Gonçalves}{Goncalves} and M. \enc{Salomé}{Salome} Cabral}

\references{Azzalini, A. (1994). Logistic regression and other discrete data models for serially correlated observations.  \emph{J. Ital. Stat. Society}, 3 (2), 169-179. 
\doi{10.1007/bf02589225}.

\enc{Gonçalves}{Goncalves}, M. Helena (2002). \emph{ Likelihood methods for discrete longitudinal data}. PhD thesis, Faculty of Sciences, University of Lisbon.
  
\enc{Gonçalves}{Goncalves}, M. Helena, Cabral, M. \enc{Salomé}{Salome}, Ruiz de Villa, M. Carme, Escrich, Eduardo and Solanas, Montse. (2007). 
Likelihood approach for count data in longitudinal experiments. \emph{Computational Statistics and Data Analysis}, 51, 12, 6511-6520. \doi{10.1016/j.csda.2007.03.002}.

\enc{Gonçalves}{Goncalves}, M. Helena and Cabral, M. \enc{Salomé}{Salome}. (2021). \code{cold}: An \code{R} Package for the Analysis of Count Longitudinal Data. \emph{Journal of Statistical Software}, 99, 3, 1--24. \doi{10.18637/jss.v099.i03}.

}

\keyword{package}

\seealso{\code{\link{cold-class}}, \code{\link{cold}}, \code{\link[methods]{Methods}}}      

