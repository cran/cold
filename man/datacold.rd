\name{datacold}
\alias{datacold}
\docType{data}
\title{Data}
\description{This example is an artificial data.}
\usage{data(datacold)}
\format{
  A data frame with 390 observations on the following 4 variables.
  \describe{
    \item{\code{Subject}}{identifies de number of the individual profile. This vector contains observations of 30 individual profiles.}
    \item{\code{Treatment}}{a factor with levels \code{0} and \code{1}.}
    \item{\code{Time}}{a numeric vector that identifies the number of the time points observed.}
    \item{\code{z}}{a numeric vector representing the response variable.}
  }
}

\examples{
data(datacold)

mod0<- cold(z~Time*Treatment, data=datacold, time="Time", 
id="Subject", dependence="ind")
summary (mod0)

modI<- cold(z~Time*Treatment, data=datacold, time="Time", 
id="Subject",  dependence="AR1")
summary (modI)

anova(mod0,modI)

plot(modI,which=1,factor=Treatment,xlab="Time (weeks)", 
ylab="Count", main="Model AR1")


### independent with random intercept 
mod0R <- cold(z ~ Time * Treatment, random = ~ 1, data = datacold, 
time = "Time", id = "Subject", dependence = "indR")
summary(mod0R)

### independent with random intercept (dependence="indR") 
### using cubature (integration = "cubature")
\donttest{
mod0R.C <- cold(z ~ Time * Treatment, random = ~ 1, data = datacold, 
time = "Time", id = "Subject", dependence = "indR", integration = "cubature")
summary(mod0R.C)

randeff(mod0R.C)

}


### dependence="indR2" 
## It takes a long time to run
\donttest{ 
## Using Monte Carlo method (integration = "MC")
mod0R2MC <- cold(z ~ Time * Treatment, ~ 1 + Time, datacold, time = "Time",
id = "Subject", dependence = "indR2", integration = "MC")

summary (mod0R2MC)

randeff(mod0R2MC)

## Using cubature (integration = "cubature")
mod0R2C<-cold(z ~ Time * Treatment, random = ~ 1 + Time, data = datacold, 
time = "Time", id = "Subject", dependence = "indR2", integration = "cubature")
summary (mod0R2C)
}

}

\keyword{datasets}
