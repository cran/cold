\name{bolus}
\alias{bolus}
\docType{data}
\title{Bolus data }
\description{The dataset has the number of requests per interval in 12 successive four-hourly intervals following abdominal surgery for 65 patients in a clinical trial to compare two groups (bolus/lock-out combinations).}
\usage{data("bolus")}
\format{
  A data frame with 780 observations on the following 4 variables.
  \describe{
    \item{\code{id}}{identifies de number of the individual profile. This vector contains observations of 65 individual profiles.}
    \item{\code{group}}{a factor with levels \code{1mg} and \code{2mg}.}
    \item{\code{time}}{a numeric vector that identifies the number of the time points observed.}
    \item{\code{y}}{a numeric vector with the number of analgesic doses taken by hospital patients in 12 successive four-hourly intervals.}
  }
}
\details{
The group \code{2mg} has 30 patients and the group \code{1mg} has 35 patients.
}

\source{
Weiss, Robert E. (2005). Modeling Longitudinal Data. Springer

https://robweiss.faculty.biostat.ucla.edu/book-data-sets
}

\references{
Henderson, R. and Shimakura, S. (2003).
A Serially Correlated Gamma Frailty Model for Longitudinal Count Data.
\emph{Biometrika}, vol. 90, No. 2, 355--366

}
\examples{
data(bolus)

## change the reference class
contrasts(bolus$group)
bolus$group<-relevel(factor(bolus$group), ref = "2mg")
contrasts(bolus$group)

## Weiss, Robert E. (2005) pp 353-356, compare with Table 11.2

bol0R <- cold(y ~ time + group, random = ~ 1, data = bolus, 
dependence = "indR")
summary (bol0R)


## reparametrization of time 
bolus$time1 <- bolus$time - 6

\donttest{ 
bol0R1 <- cold(y ~ time1 + group, random = ~ 1,data = bolus, 
dependence = "indR")
summary (bol0R1)

bol1R1 <- cold(y ~ time1 + group, random = ~ 1, data = bolus, 
time = "time1", dependence = "AR1R")
summary (bol1R1)

anova(bol0R1, bol1R1)

plot(bol1R1, which = 1, factor = group, ylab = "Bolus count")

}

}
\keyword{datasets}
