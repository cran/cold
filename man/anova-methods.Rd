\name{anova-methods}
\docType{methods}
\alias{anova-methods}
\alias{anova,cold-method}
\title{Methods for function anova in package "cold"}
\description{ Compute an analysis deviance table for two fitted model objects. }

\usage{
\S4method{anova}{cold}(object, ..., test = TRUE, correct = FALSE)
}
\arguments{ 
  \item{object}{an object of class \code{\link[=cold-class]{cold}}.}
  \item{...}{an object of class \code{\link[=cold-class]{cold}}.}
  \item{test}{an optional logical value controlling whether likelihood ratio tests 
  should be used to compare the fitted models represented by \code{object} and by \code{y}. The default is TRUE. }  
  \item{correct}{an optional logical value controlling whether the p-value of the likelihood ratio test 
  must be corrected. The default is FALSE.}  
   }  
\details{ 
\code{correct} = TRUE is used to test the presence of a random intercept term and the solution proposed by 
Self and Liang (1987) is adopted only to the p-value. }


\section{Warning}{ 
The comparison between two models by anova will only be valid if they are fitted to the same dataset.}

\section{Methods}{
\describe{
\item{\code{signature(object = "ANY")}:}{Generic function.}
\item{\code{signature(object="cold")}:}{Anova for \code{\link{cold}} object.} 
}}

\references{Self, Steven G. and Liang, Kung-Yee (1987). Asymptotic properties of maximum likelihood estimators 
and likelihood ratio tests under nonstandard conditions. 
\emph{Journal of the American Statistical Association}, 82, 605-610.}

\examples{\donttest{
#####  data= seizure

seiz1<-cold(y~lage+lbase+v4+trt+trt:lbase, data=seizure, start=NULL, 
dependence="AR1")

seiz2<-cold(y~lage+lbase+v4+trt, data=seizure, start=NULL, dependence="AR1")

anova(seiz1,seiz2)

#####  data= datacold

mod0<- cold(z~Time*Treatment, data=datacold, time="Time", id="Subject", 
aggregate=Treatment, dependence="ind")

mod0R<- cold(z~Time*Treatment, data=datacold, time="Time", id="Subject",
aggregate=Treatment, dependence="indR")

anova(mod0, mod0R, correct = TRUE)
}}


\keyword{methods}

