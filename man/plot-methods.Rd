\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,cold,missing-method}
\title{Methods for function \code{plot} }
\description{Three plots (selectable by \code{which}) are currently available:   
a plot for the fitted model (\code{which}=1), a plot for the individual mean profile (\code{which}=2) 
and a plot for the observed data and the corresponded mean profile (\code{which}=3).
}
\usage{
\S4method{plot}{cold,missing}(x,which=c(1:3),xlab=NULL,ylab=NULL,main=NULL,subSET,
ident=FALSE,caption=c("Individual mean profiles"),cex.caption=1)}
\arguments{ 
  \item{x}{an object of class \code{\link{cold}}.}
  \item{which}{if a subset of the plots is required, specify a subset of the numbers 1:3. Default is (\code{which}=1).}
  \item{xlab}{label to plots.}  
  \item{ylab}{label to plots.}  
  \item{main}{title to plots in addition to the caption.} 
  \item{subSET}{logical expression indicating elements to keep in individual mean profile plots: missing values are taken as FALSE 
  for \code{which}=2.}
  \item{ident}{logical expression indicating whether or not to add the number of the subject to individual mean profile plots.
  The \code{ident} argument is for option \code{which}=2.}
  \item{caption}{captions to appear above the plots.}   
  \item{cex.caption}{controls the size of caption.}   
    }
  
  \section{Methods}{
\describe{
 \item{\code{signature(x="ANY", y="ANY")}:}{Generic function.}
    \item{\code{signature(x="cold", y="missing")}:}{Plot diagnostics for \code{\link{cold}} object.}
}}
\examples{
#####  data= datacold
### AR1R
mod1R<- cold(z~Time*Treatment, data=datacold, time="Time", 
id="Subject", aggregate=Treatment, dependence="AR1R")

plot(mod1R, which=1, xlab="Time(weeks)", ylab="Count", 
main="Model AR1R")

plot(mod1R, which=2, xlab="Time(weeks)", ylab="Count", 
main="Model AR1R")

par(mfrow=c(1,2))
plot(mod1R, which=2, ident=TRUE, subSET=Treatment=="0", 
ylab="Count", main="0")
plot(mod1R, which=2, ident=TRUE, subSET=Treatment=="1", 
ylab="Count", main="1")
par(mfrow=c(1,1))

par(mfrow=c(2,2))
plot(mod1R, which=3, subSET=(id==c(2)), xlab="Time (weeks)", 
ylab="count", main="0_Subject2")
plot(mod1R, which=3, subSET=(id==c(10)), xlab="Time (weeks)",
ylab="Count", main="0_Subject10")
plot(mod1R, which=3, subSET=(id==c(26)), xlab="Time (weeks)",
ylab="Count", main="1_Subject26")
plot(mod1R, which=3, subSET=(id==c(30)), xlab="Time (weeks)",
ylab="Count", main="1_Subject32")
par(mfrow=c(1,1))

}

\keyword{methods}
