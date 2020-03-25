\name{anova-methods}
\docType{methods}
\alias{anova-methods}
\alias{anova,cold-method}
\title{Methods for function \code{anova} }
\description{Computes an analysis deviance table for two nested fitted model objects of class \code{\link{cold}}.}

\usage{
\S4method{anova}{cold}(object, ...)
}
\arguments{ 
  \item{object}{an object of class \code{\link[=cold-class]{cold}}.}
  \item{...}{an object of class \code{\link[=cold-class]{cold}}.}
   }  

\note{It uses the naive solution of Pinheiro et al. (2000) to calculate the p-value when the difference between the models is the number of random effects.}

\section{Warning}{ 
The comparison between two models by anova will only be valid if they are fitted to the same dataset.}


\section{Methods}{
\describe{
\item{\code{signature(object = "ANY")}:}{Generic function.}
\item{\code{signature(object="cold")}:}{Anova for \code{\link{cold}} object.} 
}
}


\references{Pinheiro, J.C. and Bates, D.M. (2000). Mixed-Effects Models in S and S-PLUS.
\emph{Springer-Verlag.} }

\examples{
#####  data = seizure

seiz1 <- cold(y ~ lage + lbase + v4 + trt + trt:lbase, data = seizure, 
start = NULL, dependence = "AR1")

seiz2 <-cold(y ~ lage + lbase + v4 + trt, data = seizure, start = NULL, 
dependence = "AR1")

anova(seiz1, seiz2)

#####  data = datacold

mod0 <- cold(z ~ Time * Treatment, data = datacold, time = "Time", 
id = "Subject", dependence = "ind")

mod0R <- cold(z  ~ Time * Treatment, random = ~ 1, data = datacold, 
time = "Time", id = "Subject", dependence = "indR")
summary(mod0R)

anova(mod0, mod0R)
}


\keyword{methods}

