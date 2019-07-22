\name{datacoldM}
\alias{datacoldM}
\docType{data}
\title{Data with missing values}
\description{ This example is an artificial data with missing values.}
\usage{data("datacoldM")}
\format{
  A data frame with 390 observations on the following 4 variables.
  \describe{
    \item{\code{Subject}}{identifies de number of the individual profile.This vector contains observations of 30 individual profiles.}
    \item{\code{Treatment}}{a factor with levels \code{0} and \code{1}.}
    \item{\code{Time}}{a numeric vector that identifies the number of the time points observed.}
    \item{\code{z}}{a numeric vector representing the response variable.}
  }
}


\examples{
data(datacoldM)
 str(datacoldM)
 
mod0.M<- cold(z~Time*Treatment, data=datacoldM, time="Time", 
id="Subject", aggregate=Treatment, dependence="ind")
summary (mod0.M)

mod1.M<- cold(z~Time*Treatment, data=datacoldM, time="Time", 
id="Subject", aggregate=Treatment, dependence="AR1")
summary (mod1.M)

anova(mod0.M,mod1.M)

}
\keyword{datasets}
