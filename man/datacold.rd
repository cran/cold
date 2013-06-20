\name{datacold}
\alias{datacold}
\docType{data}
\title{data}
\description{This example is an artificial data.}
\usage{data(datacold)}
\format{
  A data frame with 650 observations on the following 4 variables.
  \describe{
    \item{\code{Subject}}{identifies de number of the individual profile.This vector contains observations of 50 individual profiles.}
    \item{\code{Treatment}}{a factor with levels \code{0} \code{1}.}
    \item{\code{Time}}{a numeric vector that identifies the number of the time points observed.}
    \item{\code{z}}{a numeric vector representing the response variable.}
  }
}

\examples{\donttest{
data(datacold)

mod0<- cold(z~Time*Treatment, data=datacold, time="Time",id="Subject",
aggregate=Treatment, dependence="ind")
summary (mod0)

modI<- cold(z~Time*Treatment, data=datacold, time="Time",id="Subject", 
aggregate=Treatment, dependence="AR1")
summary (modI)

anova(mod0,modI)

plot(modI,which=1,xlab="Time (weeks)",ylab="Count",main="Model AR1")
}}
\keyword{datasets}
