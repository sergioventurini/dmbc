% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dmbc_classes.R
\name{initialize,dmbc_ic-method}
\alias{initialize,dmbc_ic-method}
\alias{dmbc_ic-initialize}
\title{Create an instance of the \code{dmbc_ic} class using new/initialize.}
\usage{
\S4method{initialize}{dmbc_ic}(
  .Object,
  logprior = matrix(),
  logmlik = matrix(),
  logcorrfact = matrix(),
  DCIC = matrix(),
  post.est = list(),
  est = character(),
  res_last_p = list()
)
}
\arguments{
\item{.Object}{Prototype object from the class \code{\link{dmbc_ic}}.}

\item{logprior}{An object of class \code{matrix} providing the
log-prior values corresponding to different values of \emph{p} and
\emph{G}.}

\item{logmlik}{An object of class \code{matrix} providing the
marginal log-likelihood values corresponding to different values of
\emph{p} and \emph{G}.}

\item{logcorrfact}{An object of class \code{matrix} providing the
logarithm of the correction factors corresponding to different values of
\emph{p} and \emph{G}.}

\item{DCIC}{An object of class \code{matrix} providing the values
of the DCIC index corresponding to different values of \emph{p} and
\emph{G}.}

\item{post.est}{An object of class \code{list}; named list with
elements representing the parameter estimates corresponding to different
values of \emph{p} and \emph{G}.}

\item{est}{A length-one character vector representing the estimate
type used in computing the DCIC index. Possible values are \code{mean},
\code{median}, \code{ml} and \code{map}. See \code{\link{dmbc_ic}()} for
more details about these values.}

\item{res_last_p}{An object of class \code{list}; list of
\code{dmbc_fit_list} objects with the results of fitting the DMBC
models corresponding to the last value of \emph{p}. This is needed in case
of an update of the DCIC calculations using additional \emph{p} and/or
\emph{G} values.}
}
\description{
Create an instance of the \code{dmbc_ic} class using new/initialize.
}
\author{
Sergio Venturini \email{sergio.venturini@unicatt.it}
}
