\name{InverseGamma}
\alias{InverseGamma}
\alias{dinvgamma}
\alias{rinvgamma}
\title{The Inverse Gamma Distribution}
\usage{
   dinvgamma(x, alpha, beta = 1, log = FALSE)
   rinvgamma(n, alpha, beta = 1)
}
\description{
   Density and random generation for the inverse gamma distribution with shape parameter equal to \code{alpha} and rate parameter equal to \code{beta}.
}
\arguments{
   \item{x}{vector of quantiles.}
   \item{n}{number of observations.}
   \item{alpha}{inverse gamma shape parameter. Must be strictly positive.}
   \item{beta}{inverse gamma rate parameter. Must be strictly positive.}
   \item{log}{logical; if \code{TRUE}, densities \emph{d} are given as \emph{log(d)}.}
}
\details{
   If \code{beta} is omitted, it assumes the default value of 1.

   The inverse gamma distribution with parameters \code{alpha} = \eqn{\alpha}{alpha} and \code{beta} = \eqn{\beta}{beta} has density

   \deqn{\frac{\beta^{\alpha}}{\Gamma(\alpha)} x^{-\alpha-1} \exp\left(-\frac{\beta}{x}\right)}{%
      \beta^\alpha/\Gamma(\alpha) x^(-\alpha-1) exp(-\beta/x)}

   for \eqn{x \ge 0}, \eqn{\alpha > 0} and \eqn{\beta > 0}. (Here \eqn{\Gamma(\alpha)}{\Gamma(\alpha)} is the function implemented by R's \code{\link{gamma}()} and defined in its help.)

   The inverse gamma is the distribution for the reciprocal of a gamma distributed random variable with shape parameter \eqn{\alpha} and rate parameter \eqn{\beta}.

   The mean and variance are \eqn{E(X) = \frac{\beta}{\alpha - 1}}{E(X) = \beta/(\alpha - 1)}, for \eqn{\alpha > 1}, and \eqn{Var(X) = \frac{\beta^2}{(\alpha - 1)^2(\alpha - 2)}}{\beta^2/((\alpha - 1)^2(\alpha - 2))}, for \eqn{\alpha > 2}.
}
\value{
   \code{\link{dinvgamma}} returns the \code{length(x)}-vector of densities and \code{\link{rinvgamma}} returns an \code{n}-vector of random deviates.
}
\seealso{
   \code{\link{dgamma}},
   \code{\link{rgamma}}.
}
\keyword{distribution}
\keyword{gamma}
