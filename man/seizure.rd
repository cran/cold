\name{seizure}
\alias{seizure}
\docType{data}
\title{Epileptic seizure}
\description{The dataset has the number of epileptic seizures in each of four two-week intervals, 
and in a baseline eight-week interval, for treatment and control groups with a total of 59 individuals.}
\usage{data(seizure)}
\format{
  A data frame with 236 observations on the following 9 variables.
  \describe{
    \item{\code{id}}{identifies de number of the individual profile. This vector contains observations of 59 individual profiles.}
    \item{\code{y}}{a numeric vector with the number of epileptic seizures in the four two-weeks intervals observed.}
    \item{\code{v4}}{a numeric vector indicating the fourth visit.}
    \item{\code{time}}{a numeric vector that identifies the number of the time points observed.}
    \item{\code{trt}}{a numeric vector indicator of treatment.}
    \item{\code{base}}{the number of epileptic seizures in a baseline 8-week interval.}
    \item{\code{age}}{a numeric vector of subject \code{age}.}
    \item{\code{lbase}}{recode the variable \code{base} by log(\code{base}/4).}
    \item{\code{lage}}{recode the variable \code{age} by log(\code{age}).}
  }
}

\source{Thall, P.F., and Vail, S.C. (1990). Some covariance models for longitudinal 
count data with overdispersion. \emph{Biometrics}, 46, 657--671.}

\references{
Diggle, P.J., Liang, K.Y., and Zeger, S.L. (1994). Analysis of Longitudinal Data. Clarendon Press.
}

\examples{
#####  data= seizure
 str(seizure) 

### independence 
seiz0<-cold(y~lage+lbase+v4+trt+trt*lbase, data=seizure, dependence="ind")
summary (seiz0)

### AR1
seiz1<-cold(y~lage+lbase+v4+trt+trt*lbase, data=seizure,dependence="AR1")
summary (seiz1)

### AR1R without patient 207
seizure207<-seizure[seizure$id!=207,]

seiz1R<-cold(y~lage+lbase+v4+trt+trt*lbase, data=seizure207, dependence="AR1R")
summary(seiz1R)

}
\keyword{datasets}
