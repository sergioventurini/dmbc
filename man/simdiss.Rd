\name{simdiss}
\alias{simdiss}
\docType{data}
\title{Simulated Data}
\description{
   Simulated binary dissimilarity matrices for \code{S} = 10 subjects about \code{n} = 16 items.
}
\usage{data(simdiss)}
\format{
   A named list of \code{S} = 10 \code{\link{dist}} objects each containing the binary dissimilarities expressed by the \code{S} subjects about \code{n} = 16 items.
}
\references{
   Piccarreta, R. and Venturini, S. (2014), "Model-Based Clustering of Several Binary Dissimilarity Matrices". Technical report. %\bold{Volume 2}, Number 2, 756--776.
   %\url{http://projecteuclid.org/euclid.aoas/1215118537}
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
sim.dmbc <- dmbc(simdiss, p, G, burnin, nsim, prm.prop, prm.prior, random.start = TRUE,
   verbose = TRUE)

# summarize and plot the results
summary(sim.dmbc)
plot(sim.dmbc)
}
}
\keyword{datasets}
\keyword{dissimilarity}
