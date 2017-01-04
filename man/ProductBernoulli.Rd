\name{ProductBernoulli}
\alias{ProductBernoulli}
\alias{dprodber}
\title{The Product Bernoulli Distribution}
\usage{
   dprodber(x, probs, log.p = FALSE)
}
\description{
   Probability function and random generation for the product of Bernoulli distribution with probability of successes equal to \code{probs}.
}
\arguments{
   \item{x}{vector of trial outcomes.}
   \item{probs}{vector of probabilities of success on each trial.}
   \item{log.p}{logical; if \code{TRUE}, probabilities \emph{p} are given as \emph{log(p)}.}
}
\details{
   The product of Bernoulli distribution with parameters \code{probs} = \eqn{(\pi_1,\pi_2,\ldots,\pi_m)} is a multivariate discrete distribution that corresponds to a product of \eqn{m} independent Bernoulli distributions with probability of seccess \eqn{\pi_1,\pi_2,\ldots,\pi_m} respectively (see Piccarreta and Venturini (2014)). Its probability function is defined as

   \deqn{\prod_{j = 1}^{m}{\pi_{j}^{x_{j}}(1 - \pi_{j})^{(1 - x_{j})}}}{%
      [\pi_1^x_1 * (1 - \pi_1)^(1 - x_1)]*\dots*[\pi_m^x_m * (1 - \pi_m)^(1 - x_m)]}

   with \eqn{0 < \pi_j < 1} and \eqn{x_j \in \{0, 1\}}{x_j in {0, 1}}, for \eqn{j = 1,2,\ldots,m}.
}
\value{
   \code{\link{dprodber}} returns the probability associated with the given vector of successes \code{x}.
}
\references{
   Piccarreta, R. and Venturini, S. (2014), "Model-Based Clustering of Several Binary Dissimilarity Matrices". Technical report. %\bold{Volume 2}, Number 2, 756--776.
   %\url{http://projecteuclid.org/euclid.aoas/1215118537}
}
\seealso{
   \code{\link{dbinom}}.
}
\author{Sergio Venturini \email{sergio.venturini@unibocconi.it} }
\keyword{binomial}
\keyword{distribution}
