% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dmbc.R
\name{dmbc}
\alias{dmbc}
\title{Estimation of a DMBC model.}
\usage{
dmbc(
  data,
  p = 2,
  G = 3,
  control = dmbc_control(),
  prior = NULL,
  cl = NULL,
  post_all = FALSE
)
}
\arguments{
\item{data}{An object of class \code{dmbc_data} containing the data
to analyze.}

\item{p}{A length-one numeric vector indicating the number of dimensions of the
latent space.}

\item{G}{A length-one numeric vector indicating the number of cluster to
partition the \emph{S} subjects.}

\item{control}{A list of control parameters that affect the sampling
but do not affect the posterior distribution. See
\code{\link{dmbc_control}()} for more details.}

\item{prior}{A list containing the prior hyperparameters. See
\code{\link{dmbc_prior}()} for more details.}

\item{cl}{An optional \pkg{\link{parallel}} or
\pkg{\href{https://CRAN.R-project.org/package=snow}{snow}}
cluster for use if \code{parallel = "snow"}. If not supplied, a cluster
on the local machine is created for the duration of the \code{dmbc()} call.}

\item{post_all}{A length-one logical vector, which if TRUE applies a further
post-processing to the simulated chains (in case these are more than one).}
}
\value{
A \code{dmbc_fit_list} object.
}
\description{
\code{dmbc()}, the main function of the package, estimates a DMBC model
  for a given set of \emph{S} dissimilarity matrices.
}
\examples{
\dontrun{
data(simdiss, package = "dmbc")

G <- 3
p <- 2
prm.prop <- list(z = 1.5, alpha = .75)
burnin <- 20000
nsim <- 10000
seed <- 2301

set.seed(seed)

control <- list(burnin = burnin, nsim = nsim, z.prop = prm.prop[["z"]],
  alpha.prop = prm.prop[["alpha"]], random.start = TRUE, verbose = TRUE,
  nchains = 2, thin = 10, store.burnin = TRUE, threads = 2,
  parallel = "snow")
sim.dmbc <- dmbc(simdiss, p, G, control)

summary(sim.dmbc, include.burnin = FALSE)

library(bayesplot)
library(ggplot2)
color_scheme_set("teal")
plot(sim.dmbc, what = "trace", regex_pars = "eta")

z <- dmbc_get_configuration(sim.dmbc, chain = 1, est = "mean",
  labels = 1:16)
summary(z)
color_scheme_set("mix-pink-blue")
graph <- plot(z, size = 2, size_lbl = 3)
graph + panel_bg(fill = "gray90", color = NA)
}

}
\references{
Venturini, S., Piccarreta, R. (2021), "A Bayesian Approach for Model-Based
  Clustering of Several Binary Dissimilarity Matrices: the \pkg{dmbc}
  Package in \code{R}", Journal of Statistical Software, 100, 16, 1--35, <10.18637/jss.v100.i16>.
}
\seealso{
\code{\link{bmds}} for Bayesian (metric) multidimensional scaling.

\code{\link{dmbc_data}} for a description of the data format.

\code{\link{dmbc_fit_list}} for a description of the elements
  included in the returned object.
}
\author{
Sergio Venturini \email{sergio.venturini@unicatt.it}
}
