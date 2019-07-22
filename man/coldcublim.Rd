\name{coldcublim}
\alias{coldcublim}
\title{ Auxiliary for controlling "cold" fitting}
\description{Auxiliary function as user interface for \code{cold} fitting.}
\usage{coldcublim (l1i=-4,l2i=-4,l1s=4,l2s=4, tol=1e-4, maxEval=100)}
\arguments{
  \item{l1i}{lower limit of integration for the log-likelihood.}
  \item{l1s}{upper limit of integration for the log-likelihood.}
  \item{l2i}{lower limit of integration for the log-likelihood.}
  \item{l2s}{upper limit of integration for the log-likelihood.}
   \item{tol}{The maximum tolerance.}
   \item{maxEval}{The maximum number of function evaluations needed.}                }

\value{A list with the arguments as components.}


\keyword{function}
