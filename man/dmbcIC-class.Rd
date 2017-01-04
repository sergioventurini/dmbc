\name{dmbcIC-class}
\Rdversion{1.1}
\docType{class}
\alias{dmbcIC-class}
\alias{initialize,dmbcIC-method}

\title{Class \code{"dmbcIC"}}
\description{
   This class encapsulates results of the procedure implemented in the \code{\link{dmbc.IC}} and \code{\link{dmbc.IC.more}} functions.
}
\section{Objects from the Class}{
   Objects can be created by calls of the form \code{new("dmbcIC", [all the slots below])}, but most often as the result of a call to \code{\link{dmbc.IC}} or \code{\link{dmbc.IC.more}}.
}
\section{Slots}{
  \describe{
    \item{\code{logprior}:}{Object of class \code{"matrix"}; \code{pmax}-by-\code{Gmax} matrix of logprior values (see the references).}
    \item{\code{logmlik}:}{Object of class \code{"matrix"}; \code{pmax}-by-\code{Gmax} matrix of marginal loglikelihood values (see the references).}
    \item{\code{logcorrfact}:}{Object of class \code{"matrix"}; \code{pmax}-by-\code{Gmax} matrix of logarithms of the correction factors (see the references).}
    \item{\code{DCIC}:}{Object of class \code{"matrix"}; \code{pmax}-by-\code{Gmax} matrix of DCIC values.}
    \item{\code{post.mean}:}{Object of class \code{"list"}; list of \code{pmax}-times-\code{Gmax} estimates (as indicated by the \code{est} element of the \code{control} slot).}
    \item{\code{control}:}{Object of class \code{"list"}; list of tuning parameters, i.e. the number of MCMC burnin (\code{burnin}) and sample (\code{nsim}) iterations, the proposal parameters (\code{prm.prop}), the prior hyperparameters (\code{prm.prior}), a logical (\code{random.start}) indicating whether the MCMC fitting procedures should start from a random starting point, and a character (\code{est}) providing the type of estimates to retain for every fitted model.}
    \item{\code{res_last_p}:}{Object of class \code{"list"}; list of intermediate results to use for proceeding with additional DCIC calculations afterwards.}
  }
}
\section{Methods}{
  \describe{
    % \item{initialize}{\code{signature(.Object = "dmbcIC")}: ... }
    \item{plot}{\code{signature(x = "dmbcIC", y = "missing")}: Plot the results of a call to \code{\link{dmbc.IC}} or \code{\link{dmbc.IC.more}}.}
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
\keyword{dmbcIC}
\keyword{methods}
