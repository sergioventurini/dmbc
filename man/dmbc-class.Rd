\name{dmbc-class}
\Rdversion{1.1}
\docType{class}
\alias{dmbc-class}
\alias{initialize,dmbc-method}

\title{Class \code{"dmbc"}}
\description{
   This class encapsulates results of the estimation procedure implemented in the \code{\link{dmbc}}() function.
}
\section{Objects from the Class}{
   Objects can be created by calls of the form \code{new("dmbc", [all the slots below])}, but most often as the result of a call to \code{\link{dmbc}}().
}
\section{Slots}{
  \describe{
    \item{\code{z.chain}:}{Object of class \code{"array"}; posterior draws from the MCMC algorithm for the (untransformed) latent positions \eqn{Z}.}
    \item{\code{z.chain.p}:}{Object of class \code{"array"}; posterior draws from the MCMC algorithm for the (Procrustes-transformed) latent positions \eqn{Z}.}
    \item{\code{alpha.chain}:}{Object of class \code{"matrix"}; posterior draws from the MCMC algorithm for the \eqn{\alpha} parameters.}
    \item{\code{eta.chain}:}{Object of class \code{"matrix"}; posterior draws from the MCMC algorithm for the \eqn{\eta} parameters.}
    \item{\code{sigma2.chain}:}{Object of class \code{"matrix"}; posterior draws from the MCMC algorithm for the \eqn{\sigma^2} parameters.}
    \item{\code{lambda.chain}:}{Object of class \code{"matrix"}; posterior draws from the MCMC algorithm for the \eqn{\lambda} parameters.}
    \item{\code{prob.chain}:}{Object of class \code{"array"}; posterior draws from the MCMC algorithm for the cluster membership probabilities.}
    \item{\code{x.ind.chain}:}{Object of class \code{"array"}; posterior draws from the MCMC algorithm for the cluster membership indicators.}
    \item{\code{x.chain}:}{Object of class \code{"matrix"}; posterior draws from the MCMC algorithm for the cluster membership labels.}
    \item{\code{accept}:}{Object of class \code{"matrix"}; final acceptance rates for the MCMC algorithm.}
    \item{\code{obsdiss}:}{Object of class \code{"list"}; list of the observed binary dissimilarity matrices.}
    \item{\code{dens}:}{Object of class \code{"list"}; list of the loglikelihood, logprior and logposterior values at each iteration of the MCMC algorithm.}
    \item{\code{control}:}{Object of class \code{"list"}; list of the MCMC algorithm tuning parameters, i.e. the number of MCMC burnin (\code{burnin}) and sample (\code{nsim}) iterations, the proposal parameters (\code{prm.prop}), and prior hyperparameters (\code{prm.prior}, \code{hyper.eta.a}, \code{hyper.eta.b}).}
    \item{\code{dim}:}{Object of class \code{"list"}; list of dimensions for the estimated model, i.e. the number of items (\code{n}), the number of latent dimensions (\code{p}), the number of clusters (\code{G}), and the number of subjects (\code{S}).}
  }
}
\section{Methods}{
  \describe{
    % \item{initialize}{\code{signature(.Object = "dmbc")}: Class initialization.}
    \item{plot}{\code{signature(x = "dmbc", y = "missing")}: Plot the results of a fitted \code{\link{dmbc-class}} model; it exploits the \code{coda}'s \code{\link[coda]{mcmc}} function.}
    \item{summary}{\code{signature(object = "dmbc")}: Generate a summary for the results of a fitted \code{\link{dmbc-class}} model; it exploits the \code{coda}'s \code{\link[coda]{mcmc}} function.}
	 }
}
\references{
   Piccarreta, R. and Venturini, S. (2014), "Model-Based Clustering of Several Binary Dissimilarity Matrices". Technical report. %\bold{Volume 2}, Number 2, 756--776.
   %\url{http://projecteuclid.org/euclid.aoas/1215118537}
}
\author{Sergio Venturini \email{sergio.venturini@unibocconi.it} }
\seealso{
   \code{\link{dmbc}},
   \code{\link{summary-methods}},
   \code{\link{plot-methods}}.
}
% \examples{
% showClass("dmbcIC")
% }
\keyword{classes}
\keyword{dmbc}
\keyword{methods}
