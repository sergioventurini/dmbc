% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract.R
\name{dmbc_match_groups}
\alias{dmbc_match_groups}
\title{Auxiliary function for realigning the grouping of a fitted DMBC model.}
\usage{
dmbc_match_groups(res, est = "mean", ref = 1)
}
\arguments{
\item{res}{An object of class \code{dmbc_fit_list}.}

\item{est}{A length-one character vector indicating the estimate type to use.}

\item{ref}{A length-one numeric vector indicating the chain number to use as
the reference.}
}
\value{
An object of class \code{dmbc_fit_list}.
}
\description{
\code{dmbc_match_groups()} is an auxiliary function for realigning the
  cluster membership estimates provided by the individual chains of the
  fitted model if they do not agree.
}
\examples{
\dontrun{
data(simdiss, package = "dmbc")

G <- 5
p <- 3
prm.prop <- list(z = 4, alpha = 2)
burnin <- 2000
nsim <- 1000
seed <- 2301

set.seed(seed)

control <- list(burnin = burnin, nsim = nsim, z.prop = prm.prop[["z"]],
  alpha.prop = prm.prop[["alpha"]], random.start = TRUE, verbose = TRUE,
  nchains = 6, store.burnin = TRUE, threads = 2, parallel = "snow")
sim.dmbc <- dmbc(simdiss, p, G, control)

sim.dmbc_new <- dmbc_match_groups(sim.dmbc)
}

}
\references{
Venturini, S., Piccarreta, R. (2021), "A Bayesian Approach for Model-Based
  Clustering of Several Binary Dissimilarity Matrices: the \pkg{dmbc}
  Package in \code{R}", Journal of Statistical Software, 100, 16, 1--35, <10.18637/jss.v100.i16>.
}
\seealso{
\code{\link{dmbc_check_groups}()} for checking the consistency
  of the cluster memberships across chains for a fitted DMBC model.

\code{\link{dmbc_get_configuration}()} for a description of the
  configuration extractor function.

\code{\link{dmbc_fit_list}} for a description of a fitted
  DMBC model.
}
\author{
Sergio Venturini \email{sergio.venturini@unicatt.it}
}
