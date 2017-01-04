\name{dmbc-extract}
\alias{dmbc-extract}
\alias{dmbc.get.postmean}
\alias{dmbc.get.ml}
\alias{dmbc.get.map}
\title{Estimates Extractors}
\usage{
   dmbc.get.postmean(res)
   dmbc.get.ml(res)
   dmbc.get.map(res)
}
\description{
   These functions extract the parameter estimates from an estimated \code{\link{dmbc-class}} model.
}
\arguments{
   \item{res}{object of class \code{\link{dmbc-class}} containing the results returned by the model fitting function \code{\link{dmbc}}().}
}
\details{
   \code{dmbc.get.postmean} returns the posterior mean of the parameters, \code{dmbc.get.ml} the maximum likelihood estimate and \code{dmbc.get.map} the maximum a posteriori (MAP) estimate.

   All these functions discard the initial burnin iterations provided in the call to the \code{\link{dmbc}}() function for fitting the model.
}
\value{
   \code{dmbc.get.postmean} returns a list of values, that is \code{z.postmean}, \code{alpha.postmean}, \code{eta.postmean}, \code{sigma2.postmean}, \code{lambda.postmean}, \code{prob.postmean} and \code{class.postmean}.

   \code{dmbc.get.ml} returns a list of values, that is \code{z.ml}, \code{alpha.ml}, \code{eta.ml}, \code{sigma2.ml}, \code{lambda.ml} and \code{loglik}.

   \code{dmbc.get.map} returns a list of values, that is \code{z.map}, \code{alpha.map}, \code{eta.map}, \code{sigma2.map}, \code{lambda.map}and \code{logpost}.
}
\seealso{
   \code{\link{dmbc-class}},
   \code{\link{dmbc}}.
}
\keyword{dmbc}
\keyword{extract}
