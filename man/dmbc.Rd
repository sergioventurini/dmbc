\name{dmbc}
\alias{dmbc}
\title{
   Model-Based Clustering of Dissimilarity Matrices
}
\description{
   This function implements the Bayesian procedure discussed in Piccarreta and Venturini (2014) for clustering a set of dissimilarity matrices in homogeneous groups.
}
\usage{
dmbc(D, p = 2, G = 3, burnin = 10000, nsim = 5000, prm.prop, prm.prior,
   random.start = TRUE, verbose = FALSE)
}
\arguments{
  \item{D}{list of observed dissimilarities matrices; the elements of the list must be objects of class \code{\link{dist}}.}
  \item{p}{number of dimensions of the latent space; it must be at least 1.}
  \item{G}{number of clusters; it must be at least 1.}
  \item{burnin}{number of MCMC burnin iterations.}
  \item{nsim}{number of MCMC iterations.}
  \item{prm.prop}{named list of parameters for the proposal distributions of the \eqn{Z} latent positions and \eqn{\alpha} parameters; the elements' names must be \code{"z"} and \code{"alpha"} and both need to be a single number representing the standard deviations of the normal proposal distributions for each \eqn{z_{ih}^{g}}{z_ih^g} and \eqn{\alpha_g} respectively.}
  \item{prm.prior}{named list of hyperparameters for the prior distributions of the \eqn{\sigma^2} and \eqn{\lambda} parameters; the elements' names must be \code{"sigma2"} and \code{"lambda"} with the first one a two-sized vector containing the \eqn{\sigma^2} prior hyperparameters and the second element a vector providing the \eqn{\lambda} hyperparameters.}
  \item{random.start}{logical; if \code{TRUE} (default) allows the algorithm start from a random starting point.}
  \item{verbose}{logical; if \code{TRUE}, messages informing about the progress of the computations will be printed.}
}
% \details{
%    The model implemented by this function is ...
% }
\value{
   \code{dmbc} returns an object of \code{\link{class}} \code{"dmbc-class"}.

   The functions \code{summary} and \code{plot} are used to obtain and print a summary and provide plots of the results. The generic accessor functions \code{\link{dmbc.get.postmean}}, \code{\link{dmbc.get.ml}} and \code{\link{dmbc.get.map}} extract various useful features of the value returned by \code{dmbc}.

   An object of class \code{"dmbc-class"} is an S4 object containing the following components:
   \item{\code{z.chain}}{array; posterior draws from the MCMC algorithm for the (untransformed) latent positions \eqn{Z}.}
   \item{\code{z.chain.p}}{array; posterior draws from the MCMC algorithm for the (Procrustes-transformed) latent positions \eqn{Z}.}
   \item{\code{alpha.chain}}{matrix; posterior draws from the MCMC algorithm for the \eqn{\alpha} parameters.}
   \item{\code{eta.chain}}{matrix; posterior draws from the MCMC algorithm for the \eqn{\eta} parameters.}
   \item{\code{sigma2.chain}}{matrix; posterior draws from the MCMC algorithm for the \eqn{\sigma^2} parameters.}
   \item{\code{lambda.chain}}{matrix; posterior draws from the MCMC algorithm for the \eqn{\lambda} parameters.}
   \item{\code{prob.chain}}{array; posterior draws from the MCMC algorithm for the cluster membership probabilities.}
   \item{\code{x.ind.chain}}{array; posterior draws from the MCMC algorithm for the cluster membership indicators.}
   \item{\code{x.chain}}{matrix; posterior draws from the MCMC algorithm for the cluster membership labels.}
   \item{\code{accept}}{matrix; final acceptance rates for the MCMC algorithm.}
   \item{\code{obsdiss}}{list of the observed binary dissimilarity matrices.}
   \item{\code{dens}}{list of the loglikelihood, logprior and logposterior values at each iteration of the MCMC algorithm.}
   \item{\code{control}}{list of the MCMC algorithm tuning parameters, i.e. the number of MCMC burnin (\code{burnin}) and sample (\code{nsim}) iterations, the proposal parameters (\code{prm.prop}), and prior hyperparameters (\code{prm.prior}, \code{hyper.eta.a}, \code{hyper.eta.b}).}
   \item{\code{dim}}{list of dimensions for the estimated model, i.e. the number of items (\code{n}), the number of latent dimensions (\code{p}), the number of clusters (\code{G}), and the number of subjects (\code{S}).}
}
\references{
   Piccarreta, R. and Venturini, S. (2014), "Model-Based Clustering of Several Binary Dissimilarity Matrices". Technical report. %\bold{Volume 2}, Number 2, 756--776.
   %\url{http://projecteuclid.org/euclid.aoas/1215118537}
}
\author{Sergio Venturini \email{sergio.venturini@unibocconi.it} }
\seealso{
   \code{\link{dmbc-class}},
   \code{\link{dmbcIC-class}},
   \code{\link{dmbc.IC}}.
}
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
sim.dmbc <- dmbc(simdiss, p, G, burnin, nsim, prm.prop, prm.prior,
  random.start = TRUE, verbose = TRUE)

# summarize and plot the results
summary(sim.dmbc)
plot(sim.dmbc)
}
}
\keyword{cluster}
\keyword{dissimilarity}
\keyword{mixture}
