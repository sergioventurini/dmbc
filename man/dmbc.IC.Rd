\name{dmbc.IC}
\alias{dmbc.IC}
\alias{dmbc.IC.more}
\title{
   Model Choice for the \code{dmbc} Procedure
}
\description{
   These functions provide the calculations for some information criteria, namely the DCIC, as discussed in Piccarreta and Venturini (2014).
}
\usage{
dmbc.IC(D, pmax = 3, Gmax = 5, burnin, nsim, prm.prop, prm.prior, random.start,
   est = "post.mean", verbose = TRUE)

dmbc.IC.more(D, pmax = NULL, Gmax = NULL, IC.res, verbose = TRUE)
}
\arguments{
  \item{D}{list of observed dissimilarities matrices; the elements of the list must be objects of class \code{\link{dist}}.}
  \item{pmax}{maximum number of dimensions of the latent space; it must be at least 1.}
  \item{Gmax}{maximum number of clusters; it must be at least 1.}
  \item{burnin}{number of MCMC burnin iterations.}
  \item{nsim}{number of MCMC iterations.}
  \item{prm.prop}{named list of parameters for the proposal distributions of the \eqn{Z} latent positions and \eqn{\alpha} parameters; the elements' names must be \code{"z"} and \code{"alpha"} and both need to be a single number representing the standard deviations of the normal proposal distributions for each \eqn{z_{ih}^{g}}{z_ih^g} and \eqn{\alpha_g} respectively.}
  \item{prm.prior}{named list of hyperparameters for the prior distributions of the \eqn{\sigma^2} and \eqn{\lambda} parameters; the elements' names must be \code{"sigma2"} and \code{"lambda"} with the first one a two-sized vector containing the \eqn{\sigma^2} prior hyperparameters and the second element a vector providing the \eqn{\lambda} hyperparameters.}
  \item{random.start}{logical; if \code{TRUE} (default) allows the algorithm start from a random starting point.}
  \item{est}{character; type of estimates to retain for every fitted model; either one of \code{"postmean"}, \code{"ml"} or \code{"map"}.}
  \item{IC.res}{list of intermediate results for proceeding with additional DCIC and BIC calculations.}
  \item{verbose}{logical; if \code{TRUE}, messages informing about the progress of the computations will be printed.}
}  
% \details{
%    The procedure implemented in this function is ...
% }
\value{
   \code{dmbc.IC} and \code{dmbc.IC.more} return an object of \code{\link{class}} \code{"dmbcIC-class"}.

   The function \code{plot} provides summary plots of the results.

   An object of class \code{"dmbcIC-class"} is an S4 object containing the following components:
   \item{\code{logprior}}{\code{pmax}-by-\code{Gmax} matrix of logprior values (see the references).}
   \item{\code{logmlik}}{\code{pmax}-by-\code{Gmax} matrix of marginal loglikelihood values (see the references).}
   \item{\code{logcorrfact}}{\code{pmax}-by-\code{Gmax} matrix of logarithms of the correction factors (see the references).}
   \item{\code{DCIC}}{\code{pmax}-by-\code{Gmax} matrix of DCIC values.}
   \item{\code{post.mean}}{list of \code{pmax}-times-\code{Gmax} estimates (as indicated by the \code{est} element of the \code{control} slot).}
   \item{\code{control}}{list of tuning parameters, i.e. the number of MCMC burnin (\code{burnin}) and sample (\code{nsim}) iterations, the proposal parameters (\code{prm.prop}), the prior hyperparameters (\code{prm.prior}), a logical (\code{random.start}) indicating whether the MCMC fitting procedures should start from a random starting point, and a character (\code{est}) providing the type of estimates to retain for every fitted model.}
   \item{\code{res_last_p}}{list of intermediate results to use for proceeding with additional DCIC calculations afterwards.}
}
\note{
   Differently from \code{\link{dmbc}}(), in these functions it is assumed that the hyperparameters for the \eqn{\lambda}'s are set to \code{rep(1, G)} for every model. Therefore, the \code{prm.prior} component in the \code{control} argument must include only the vector of hyperparameters for the \eqn{\sigma^2}'s.
}
% \section{Warning}{
%    The procedure implemented in this function are experimental.
% }
\references{
   Piccarreta, R. and Venturini, S. (2014), "Model-Based Clustering of Several Binary Dissimilarity Matrices". Technical report. %\bold{Volume 2}, Number 2, 756--776.
   %\url{http://projecteuclid.org/euclid.aoas/1215118537}
}
\author{Sergio Venturini \email{sergio.venturini@unibocconi.it} }
\seealso{
   \code{\link{dmbc-class}},
   \code{\link{dmbcIC-class}},
   \code{\link{dmbc}}.
}
\examples{
\dontrun{
# ATTENTION: the following code may take some time to run! #

data(simdiss, package = "dmbc")

# model choice
pmax <- 2
Gmax <- 4
prm.prop <- list(z = 1.5, alpha = .75)
prm.prior <- c(1e-1, 1e-1)
burnin <- 1000
nsim <- 2000
seed <- 1809

set.seed(seed)
sim.ic <- dmbc.IC(D = simdiss, pmax = pmax, Gmax = Gmax, burnin = burnin, nsim = nsim,
   prm.prop = prm.prop, prm.prior = prm.prior, random.start = TRUE)

pmax <- pmax + 1
new.ic <- dmbc.IC.more(D = simdiss, pmax = pmax, Gmax = Gmax, sim.ic, verbose = TRUE)

# plot the results
plot(new.ic)
}
}
\keyword{cluster}
\keyword{dissimilarity}
\keyword{mixture}
