\name{Dirichlet}
\alias{Dirichlet}
\alias{ddirichlet}
\alias{rdirichlet}
\title{The Dirichlet Distribution}
\usage{
   ddirichlet(x, alpha)
   rdirichlet(n, alpha)
}
\description{
   Density and random generation for the Dirichlet distribution with parameter vector equal to \code{alpha}.
}
\arguments{
   \item{x}{vector or matrix of quantiles.}
   \item{n}{number of observations.}
   \item{alpha}{vector of parameters. All elements must be strictly positive.}
}
\value{
   \code{\link{ddirichlet}} returns the \code{nrow(x)}-vector of densities and \code{\link{rdirichlet}} returns an \code{n}-by-\code{length(alpha)}-matrix of random deviates.
}
\keyword{Dirichlet}
\keyword{distribution}
