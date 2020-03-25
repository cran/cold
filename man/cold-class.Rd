\name{cold-class}
\Rdversion{1.0}
\docType{class}
\alias{cold-class}
\title{Class \code{cold} of a maximum likelihood estimation}
\description{This class encapsulates results of a maximum likelihood procedure.}
\section{Objects from the Class}{Objects can be created by calls of the form \code{new("cold", ...)}, 
but most often as the result of a call to \code{\link{cold}}.}
\section{Slots}{
  \describe{
    \item{\code{coefficients}:}{Object of class \code{"matrix"}. Estimated parameters.}
    \item{\code{se}:}{Object of class \code{"matrix"}. Standard errors of estimated parameters.}
    \item{\code{covariance}:}{Object of class \code{"matrix"}. Covariance of estimated parameters.}
    \item{\code{correlation}:}{Object of class \code{"matrix"}. Correlation of estimated parameters.}
    \item{\code{log.likelihood}:}{Object of class \code{"numeric"}. The value of the log likelihood.}
    \item{\code{message}:}{Object of class \code{"integer"}. A character string giving any additional information returned by the optimizer, or NULL. 
    See \code{\link[stats]{optim}} for details.}
    \item{\code{n.cases}:}{Object of class \code{"numeric"}. Number of individual profiles used in the optimization procedure.}
    \item{\code{ni.cases}:}{Object of class \code{"numeric"}. Number of individual profiles in the dataset.}
    \item{\code{aic}:}{Object of class \code{"numeric"}. The Akaike information criterion for a fitted model object.}
    \item{\code{Fitted}:}{Object of class \code{"numeric"}. The fitted values for the estimated parameters.}
    \item{\code{bi.estimate}:}{Object of class \code{"numeric"}. The estimated values for the individual random effects.}
    \item{\code{Fitted.av}:}{Object of class \code{"numeric"}. }
    \item{\code{Time}:}{Object of class \code{"numeric"}. Vector of time points.}
    \item{\code{model.matrix}:}{Object of class \code{"matrix"}. The model matrix.}
    \item{\code{y.matrix}:}{Object of class \code{"matrix"}. The matrix of response values.}
    \item{\code{random.matrix}:}{Object of class \code{"matrix"}. The model  matrix of random effects.}
    \item{\code{subset.data}:}{Object of class \code{"data.frame"}. The data subset if considered.}
    \item{\code{final.data}:}{Object of class \code{"data.frame"}. The data set considered.}
    \item{\code{y.av}:}{Object of class \code{"numeric"}. The average of the response value over an individual profile.}
    \item{\code{data.id}:}{Object of class \code{"numeric"}.  Vector of individual observations.}
    \item{\code{call}:}{Object of class \code{"language"}. The call to \code{"cold"}.} 
  }
}
\section{Methods}{
  \describe{
    \item{anova}{\code{signature(object="cold")}: Anova table.}
     \item{coeftest}{\code{signature(object="cold")}: A list of summary 
     statistics of the fixed effects coefficients.}  
  \item{fitted}{\code{signature(object="cold")}: The fitted values of  a 
    fitted model.}
   \item{fixeff}{\code{signature(object="cold")}: The values corresponding to the fixed effects of a fitted model.}  
    \item{getAIC}{\code{signature(object="cold")}: A numeric value corresponding to the AIC of the fitted model.}
  \item{getcoef}{\code{signature(object="cold")}: The values corresponding to the coefficient estimates of the fitted model.}
    \item{getLogLik}{\code{signature(object="cold")}: A numeric value corresponding to the log-Likelihood of the fitted model.}
  \item{getvcov}{\code{signature(object="cold")}: The variance-covariance matrix of the fitted model.}
    \item{model.mat}{\code{signature(object="cold")}: The fixed effects model matrix of the fitted model.}
    \item{plot}{\code{signature(x="cold", y="missing")}: Plots three type of plots.}
	 \item{randeff}{\code{signature(object="cold")}: A data frame corresponding to the conditional random effects of the fitted model.}
    \item{show}{\code{signature(object="cold")}: Display object briefly.}
    \item{summary}{\code{signature(object="cold")}: Generate object summary.}
    \item{vareff}{\code{signature(object="cold")}: Numeric value(s) corresponding to the estimated random effect(s) variance of the fitted model.}

}
}
\keyword{class}

