\name{summary-methods}
\docType{methods}
\alias{summary-methods}
\alias{summary,ANY-method}
\alias{summary,dmbc-method}
\title{Summary of \code{dmbc} Objects}
\description{
	\code{summary} method for class \code{"dmbc"}. This function allows to summarize the output of the fitting procedure provided by the \code{\link{dmbc}}() function.
}
\usage{
\S4method{summary}{dmbc}(object, summary.Z = FALSE, summary.alpha = TRUE,
   summary.eta = TRUE, summary.sigma2 = TRUE, summary.lambda = TRUE,
   quantiles = c(0.025, 0.25, 0.5, 0.75, 0.975), \dots)
}
\arguments{
   \item{object}{an object of class \code{"dmbc"}, as returned by the \code{\link{dmbc}}() function.}
   \item{summary.Z}{logical; if \code{TRUE}, a summary for the latent positions \eqn{Z} is returned.}
   \item{summary.alpha}{logical; if \code{TRUE} (default), a summary for the \eqn{\alpha} parameters is returned.}
   \item{summary.eta}{logical; if \code{TRUE} (default), a summary for the \eqn{\eta} parameters is returned.}
   \item{summary.sigma2}{logical; if \code{TRUE} (default), a summary for the \eqn{\sigma^2} parameters is returned.}
   \item{summary.lambda}{logical; if \code{TRUE} (default), a summary for the \eqn{\lambda} parameters is returned.}
   \item{quantiles}{a vector of quantiles to evaluate for each variable.}
   \item{\dots}{further arguments passed to or from other methods.}
}
%%\details{
%%  ~~ If necessary, more details than the description above ~~
%%}
\value{
	The function \code{summary} computes and returns a list of summary statistics of the fitted model in \code{object}, in particular:
   \item{z}{summary of the latent positions \eqn{Z}.}
   \item{alpha}{summary of the \eqn{\alpha} parameters.}
   \item{eta}{summary of the \eqn{\eta} parameters.}
   \item{sigma2}{summary of the \eqn{\sigma^2} parameters.}
   \item{lambda}{summary of the \eqn{\lambda} parameters.}
}
\references{
   Piccarreta, R. and Venturini, S. (2014), "Model-Based Clustering of Several Binary Dissimilarity Matrices". Technical report. %\bold{Volume 2}, Number 2, 756--776.
   %\url{http://projecteuclid.org/euclid.aoas/1215118537}
}
\author{Sergio Venturini \email{sergio.venturini@unibocconi.it} }
\seealso{
   \code{\link{dmbc-class}},
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
summary(sim.dmbc)
}
}
