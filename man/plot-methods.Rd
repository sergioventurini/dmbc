\name{plot-methods}
\docType{methods}
\alias{plot-methods}
\alias{plot,ANY,ANY-method}
\alias{plot,dmbc,missing-method}
\alias{plot,dmbcIC,ANY-method}
\alias{plot,dmbcIC-method}
\alias{plot-dmbc}
\alias{plot-dmbcIC}
\title{Summary Plots of \code{dmbc} and \code{dmbcIC} Objects}
\description{
	\code{plot} methods for classes \code{"dmbc"} and \code{"dmbcIC"}. These functions allow to plot the outputs provided by the \code{\link{dmbc}}(), \code{\link{dmbc.IC}}() and \code{\link{dmbc.IC.more}}() functions.
}
\usage{
\S4method{plot}{dmbc,missing}(x, plot.Z = TRUE, plot.alpha = TRUE, plot.eta = TRUE, plot.sigma2 = TRUE,
   plot.lambda = TRUE, plot.dens = TRUE, \dots)

\S4method{plot}{dmbcIC}(x, \dots)
}
\arguments{
   \item{x}{an object of class \code{"dmbc"} or \code{"dmbcIC"}, as returned by the \code{\link{dmbc}}() function or \code{\link{dmbc.IC}}() and \code{\link{dmbc.IC.more}}() respectively.}
   \item{plot.Z}{logical; if \code{TRUE} (default), a summary plot for the latent positions \eqn{Z} is produced.}
   \item{plot.alpha}{logical; if \code{TRUE} (default), a summary plot for the \eqn{\alpha} parameters is produced.}
   \item{plot.eta}{logical; if \code{TRUE} (default), a summary plot for the \eqn{\eta} parameters is produced.}
   \item{plot.sigma2}{logical; if \code{TRUE} (default), a summary plot for the \eqn{\sigma^2} parameters is produced.}
   \item{plot.lambda}{logical; if \code{TRUE} (default), a summary plot for the \eqn{\lambda} parameters is produced.}
   \item{plot.dens}{logical; if \code{TRUE} (default), a summary plot for the loglikelihood, logprior and logposterior values is produced.}
   \item{\dots}{further arguments passed to or from other methods.}
}
\references{
   Piccarreta, R. and Venturini, S. (2014), "Model-Based Clustering of Several Binary Dissimilarity Matrices". Technical report. %\bold{Volume 2}, Number 2, 756--776.
   %\url{http://projecteuclid.org/euclid.aoas/1215118537}
}
\author{Sergio Venturini \email{sergio.venturini@unibocconi.it} }
\seealso{
   \code{\link{dmbc-class}},
   \code{\link{dmbc.IC}},
   \code{\link{dmbc.IC.more}},
   \code{\link{dmbc}}.
}
\keyword{dissimilarity}
\keyword{methods}
\keyword{mixture}
\examples{
\dontrun{
# ATTENTION: the following code may take some time to run! #

data(simdiss, package = "dmbc")

G <- 3
p <- 2
prm.prop <- list(z = 1.5, alpha = .75)
prm.prior <- list(sigma2 = c(1e-1, 1e-1), lambda = rep(1, G))
burnin <- 20000
nsim <- 10000
seed <- 1406

set.seed(seed)
sim.dmbc <- dmbc(simdiss, p, G, burnin, nsim, prm.prop, prm.prior, random.start = TRUE,
   verbose = TRUE)

# summarize the results
plot(sim.dmbc)

# compute information criteria
pmax <- 3
Gmax <- 4
prm.prop <- list(z = 1.5, alpha = .75)
prm.prior <- c(1e-1, 1e-1)
burnin <- 1000
nsim <- 2000
seed <- 1809

set.seed(seed)
sim.ic <- dmbc.IC(simdiss, pmax = pmax, Gmax = Gmax, burnin = burnin, nsim = nsim,
   prm.prop, prm.prior, random.start = TRUE)

# plot the information criteria results
plot(new.ic)
}
}
