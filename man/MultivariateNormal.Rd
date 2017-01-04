\name{MultivariateNormal}
\alias{MultivariateNormal}
\alias{dmultinorm}
\alias{rmultinorm}
\title{The Multivariate Normal Distribution}
\usage{
   dmultinorm(x, mean, sigma, log = FALSE)
   rmultinorm(n, mean, sigma)
}
\description{
   Density and random generation for the multivariate normal distribution with mean vector equal to \code{mean} and covariance matrix equal to \code{sigma}.
}
\arguments{
   \item{x}{vector or matrix of quantiles. If \code{x} is a matrix, each row is taken to be a quantile.}
   \item{n}{number of observations.}
   \item{mean}{mean vector.}
   \item{sigma}{covariance matrix.}
   \item{log}{logical; if \code{TRUE}, densities \emph{d} are given as \emph{log(d)}.}
}
\value{
   \code{\link{dmultinorm}} returns an \code{nrow(x)}-vector of densities and \code{\link{rmultinorm}} returns an \code{n}-by-\code{length(mean)}-matrix of random deviates.
}
\seealso{
   \code{\link{dnorm}},
   \code{\link{rnorm}}.
}
\keyword{distribution}
\keyword{multivariate}
\keyword{normal}
