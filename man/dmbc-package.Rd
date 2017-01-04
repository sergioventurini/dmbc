\name{dmbc-package}
\docType{package}
\alias{dmbc-package}
\title{Model-Based Clustering of Dissimilarity Measurements}
\description{
   This package implements a Bayesian approach for clustering a number of dissimilarity matrices into disjoint groups. The unconditional distribution of the dissimilarities are modeled as a mixture of multivariate random vectors. The model parameter estimates are derived using a hybrid Metropolis-Gibbs Markov Chain Monte Carlo algorithm. At the moment the package includes methods for binary dissimilarity data. See Piccarreta and Venturini (2014).
}
\references{
   Piccarreta, R. and Venturini, S. (2014), "Model-Based Clustering of Several Binary Dissimilarity Matrices". Technical report. %\bold{Volume 2}, Number 2, 756--776.
   %\url{http://projecteuclid.org/euclid.aoas/1215118537}
}
\author{Sergio Venturini \email{sergio.venturini@unibocconi.it} }
\seealso{
   \code{\link{dmbc}},
   \code{\link{dmbc.IC}}.
}
\keyword{cluster}
\keyword{dissimilarity}
\keyword{mixture}
