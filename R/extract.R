#' Extractor function for a fitted DMBC model.
#'
#' \code{dmbc_get_postmean()} is an extractor function for extracting the
#'   posterior mean estimates of the parameters for a fitted DMBC model.
#'
#' @param res An object of class \code{dmbc_fit_list}.
#' @param chain A length-one numeric vector indicating the MCMC chain number
#'   to use.
#'
#' @return A named \code{list} with the following elements:
#'   \describe{
#'     \item{\code{z}: }{array of latent coordinates posterior mean estimates}
#'     \item{\code{alpha}: }{numeric vector of alpha posterior mean estimates}
#'     \item{\code{eta}: }{numeric vector of eta posterior mean estimates}
#'     \item{\code{sigma2}: }{numeric vector of sigma2 posterior mean estimates}
#'     \item{\code{lambda}: }{numeric vector of lambda posterior mean estimates}
#'     \item{\code{prob}: }{numeric matrix of probability posterior mean estimates}
#'     \item{\code{cluster}: }{numeric vector of cluster membership posterior
#'       mean estimates}
#'     \item{\code{chain}: }{length-one numeric vector of the MCMC chain number
#'       used}
#'   }
#'
#' @author Sergio Venturini \email{sergio.venturini@@unibocconi.it}
#'
#' @seealso \code{\link{dmbc_data}} for a description of the data format.
#' @seealso \code{\link{dmbc_fit_list}} for a description of the elements
#'   included in the returned object.
#'
#' @references
#'   Venturini, S., Piccarreta, R. (2019), "A Bayesian Approach for Model-Based
#'   Clustering of Several Binary Dissimilarity Matrices: the \pkg{dmbc}
#'   Package in \code{R}", Technical report.
#'
#' @examples
#' data(simdiss, package = "dmbc")
#'
#' G <- 3
#' p <- 2
#' prm.prop <- list(z = 1.5, alpha = .75)
#' burnin <- 2000
#' nsim <- 1000
#' seed <- 2301
#'
#' set.seed(seed)
#'
#' control <- list(burnin = burnin, nsim = nsim, z.prop = prm.prop[["z"]],
#'   alpha.prop = prm.prop[["alpha"]], random.start = TRUE, verbose = TRUE,
#'   nchains = 2, thin = 10, store.burnin = TRUE, threads = 2,
#'   parallel = "multicore")
#' sim.dmbc <- dmbc(simdiss, p, G, control)
#'
#' dmbc_get_postmean(sim.dmbc, chain = 1)
#' @export
dmbc_get_postmean <- function(res, chain = 1) {
  res_chain <- res@results[[chain]]
  control <- res@results[[1]]@control
  burnin <- control[["burnin"]]
  nsim <- control[["nsim"]]
  thin <- control[["thin"]]
  nchains <- control[["nchains"]]
  store.burnin <- control[["store.burnin"]]
  totiter <- burnin + nsim
  
  if (chain > nchains)
    stop("the specified chain is not available.")

  if (store.burnin) {
    todrop <- seq(1, burnin, by = thin)
    tokeep <- seq(1, totiter, by = thin)
    tokeep <- (length(todrop) + 1):length(tokeep)
  } else {
    tokeep <- seq(1, nsim, by = thin)
    tokeep <- 1:length(tokeep)
  }

  z.postmean <- apply(res_chain@z.chain.p[tokeep, , , , drop = FALSE], c(2, 3, 4), mean, na.rm = TRUE)
  alpha.postmean <- colMeans(res_chain@alpha.chain[tokeep, , drop = FALSE], na.rm = TRUE)
  eta.postmean <- colMeans(res_chain@eta.chain[tokeep, , drop = FALSE], na.rm = TRUE)
  sigma2.postmean <- colMeans(res_chain@sigma2.chain[tokeep, , drop = FALSE], na.rm = TRUE)
  lambda.postmean <- colMeans(res_chain@lambda.chain[tokeep, , drop = FALSE], na.rm = TRUE)
  prob.postmean <- apply(res_chain@prob.chain[tokeep, , , drop = FALSE], c(2, 3), mean, na.rm = TRUE)
  class.postmean <- apply(prob.postmean, 1, which.max)
  out <- list(z = z.postmean,
              alpha = alpha.postmean,
              eta = eta.postmean,
              sigma2 = sigma2.postmean,
              lambda = lambda.postmean,
              prob = prob.postmean,
              cluster = class.postmean,
              chain = chain)

  return(out)
}

#' Extractor function for a fitted DMBC model.
#'
#' \code{dmbc_get_postmedian()} is an extractor function for extracting the
#'   posterior median estimates of the parameters for a fitted DMBC model.
#'
#' @param res An object of class \code{dmbc_fit_list}.
#' @param chain A length-one numeric vector indicating the MCMC chain number
#'   to use.
#'
#' @return A named \code{list} with the following elements:
#'   \describe{
#'     \item{\code{z}: }{array of latent coordinates posterior median estimates}
#'     \item{\code{alpha}: }{numeric vector of alpha posterior median estimates}
#'     \item{\code{eta}: }{numeric vector of eta posterior median estimates}
#'     \item{\code{sigma2}: }{numeric vector of sigma2 posterior median estimates}
#'     \item{\code{lambda}: }{numeric vector of lambda posterior median estimates}
#'     \item{\code{prob}: }{numeric matrix of probability posterior median estimates}
#'     \item{\code{cluster}: }{numeric vector of cluster membership posterior
#'       median estimates}
#'     \item{\code{chain}: }{length-one numeric vector of the MCMC chain number
#'       used}
#'   }
#'
#' @author Sergio Venturini \email{sergio.venturini@@unibocconi.it}
#'
#' @seealso \code{\link{dmbc_data}} for a description of the data format.
#' @seealso \code{\link{dmbc_fit_list}} for a description of the elements
#'   included in the returned object.
#'
#' @references
#'   Venturini, S., Piccarreta, R. (2019), "A Bayesian Approach for Model-Based
#'   Clustering of Several Binary Dissimilarity Matrices: the \pkg{dmbc}
#'   Package in \code{R}", Technical report.
#'
#' @examples
#' data(simdiss, package = "dmbc")
#'
#' G <- 3
#' p <- 2
#' prm.prop <- list(z = 1.5, alpha = .75)
#' burnin <- 2000
#' nsim <- 1000
#' seed <- 2301
#'
#' set.seed(seed)
#'
#' control <- list(burnin = burnin, nsim = nsim, z.prop = prm.prop[["z"]],
#'   alpha.prop = prm.prop[["alpha"]], random.start = TRUE, verbose = TRUE,
#'   nchains = 2, thin = 10, store.burnin = TRUE, threads = 2,
#'   parallel = "multicore")
#' sim.dmbc <- dmbc(simdiss, p, G, control)
#'
#' dmbc_get_postmedian(sim.dmbc, chain = 1)
#' @export
dmbc_get_postmedian <- function(res, chain = 1) {
  res_chain <- res@results[[chain]]
  control <- res@results[[1]]@control
  burnin <- control[["burnin"]]
	nsim <- control[["nsim"]]
  thin <- control[["thin"]]
  nchains <- control[["nchains"]]
  store.burnin <- control[["store.burnin"]]
	totiter <- burnin + nsim
	
  if (chain > nchains)
    stop("the specified chain is not available.")

  if (store.burnin) {
    todrop <- seq(1, burnin, by = thin)
    tokeep <- seq(1, totiter, by = thin)
    tokeep <- (length(todrop) + 1):length(tokeep)
  } else {
    tokeep <- seq(1, nsim, by = thin)
    tokeep <- 1:length(tokeep)
  }

	z.postmedian <- apply(res_chain@z.chain.p[tokeep, , , , drop = FALSE], c(2, 3, 4), median, na.rm = TRUE)
	alpha.postmedian <- colMedians(res_chain@alpha.chain[tokeep, , drop = FALSE], na.rm = TRUE)
	eta.postmedian <- colMedians(res_chain@eta.chain[tokeep, , drop = FALSE], na.rm = TRUE)
	sigma2.postmedian <- colMedians(res_chain@sigma2.chain[tokeep, , drop = FALSE], na.rm = TRUE)
	lambda.postmedian <- colMedians(res_chain@lambda.chain[tokeep, , drop = FALSE], na.rm = TRUE)
	prob.postmedian <- apply(res_chain@prob.chain[tokeep, , , drop = FALSE], c(2, 3), median, na.rm = TRUE)
	class.postmedian <- apply(prob.postmedian, 1, which.max)
	out <- list(z = z.postmedian,
              alpha = alpha.postmedian,
              eta = eta.postmedian,
              sigma2 = sigma2.postmedian,
		          lambda = lambda.postmedian,
              prob = prob.postmedian,
              cluster = class.postmedian,
              chain = chain)

	return(out)
}

#' Extractor function for a fitted DMBC model.
#'
#' \code{dmbc_get_ml()} is an extractor function for extracting the
#'   maximum likelihood estimates of the parameters for a fitted DMBC model.
#'
#' @param res An object of class \code{dmbc_fit_list}.
#' @param chain A length-one numeric vector indicating the MCMC chain number
#'   to use.
#'
#' @return A named \code{list} with the following elements:
#'   \describe{
#'     \item{\code{z}: }{array of latent coordinates posterior mean estimates}
#'     \item{\code{alpha}: }{numeric vector of alpha posterior mean estimates}
#'     \item{\code{eta}: }{numeric vector of eta posterior mean estimates}
#'     \item{\code{sigma2}: }{numeric vector of sigma2 posterior mean estimates}
#'     \item{\code{lambda}: }{numeric vector of lambda posterior mean estimates}
#'     \item{\code{prob}: }{numeric matrix of probability posterior mean estimates}
#'     \item{\code{cluster}: }{numeric vector of cluster membership posterior
#'       mean estimates}
#'     \item{\code{loglik}: }{length-one numeric vector of the maximum
#'       log-likelihood value}
#'     \item{\code{chain}: }{length-one numeric vector of the MCMC chain number
#'       used}
#'   }
#'
#' @author Sergio Venturini \email{sergio.venturini@@unibocconi.it}
#'
#' @seealso \code{\link{dmbc_data}} for a description of the data format.
#' @seealso \code{\link{dmbc_fit_list}} for a description of the elements
#'   included in the returned object.
#'
#' @references
#'   Venturini, S., Piccarreta, R. (2019), "A Bayesian Approach for Model-Based
#'   Clustering of Several Binary Dissimilarity Matrices: the \pkg{dmbc}
#'   Package in \code{R}", Technical report.
#'
#' @examples
#' data(simdiss, package = "dmbc")
#'
#' G <- 3
#' p <- 2
#' prm.prop <- list(z = 1.5, alpha = .75)
#' burnin <- 2000
#' nsim <- 1000
#' seed <- 2301
#'
#' set.seed(seed)
#'
#' control <- list(burnin = burnin, nsim = nsim, z.prop = prm.prop[["z"]],
#'   alpha.prop = prm.prop[["alpha"]], random.start = TRUE, verbose = TRUE,
#'   nchains = 2, thin = 10, store.burnin = TRUE, threads = 2,
#'   parallel = "multicore")
#' sim.dmbc <- dmbc(simdiss, p, G, control)
#'
#' dmbc_get_ml(sim.dmbc, chain = 1)
#' @export
dmbc_get_ml <- function(res, chain = 1) {
  res_chain <- res@results[[chain]]
  control <- res@results[[1]]@control
  burnin <- control[["burnin"]]
  nsim <- control[["nsim"]]
  thin <- control[["thin"]]
  nchains <- control[["nchains"]]
  store.burnin <- control[["store.burnin"]]
  totiter <- burnin + nsim
	ll <- res_chain@dens$loglik
  
  if (chain > nchains)
    stop("the specified chain is not available.")

  if (store.burnin) {
    todrop <- seq(1, burnin, by = thin)
    tokeep <- seq(1, totiter, by = thin)
    tokeep <- (length(todrop) + 1):length(tokeep)
  } else {
    tokeep <- seq(1, nsim, by = thin)
    tokeep <- 1:length(tokeep)
  }

  i.ml <- which.max(ll[tokeep])
  if (store.burnin) {
    i.ml <- length(todrop) + i.ml
  }

  prob.ml <- res_chain@prob.chain[i.ml, , ]
	out <- list(z = res_chain@z.chain.p[i.ml, , , ],
              alpha = res_chain@alpha.chain[i.ml, ],
          		eta = res_chain@eta.chain[i.ml, ],
              sigma2 = res_chain@sigma2.chain[i.ml, ],
              lambda = res_chain@lambda.chain[i.ml, ],
              prob = prob.ml,
          		cluster = apply(prob.ml, 1, which.max),
              loglik = ll[i.ml],
              chain = chain)

	return(out)
}

#' Extractor function for a fitted DMBC model.
#'
#' \code{dmbc_get_map()} is an extractor function for extracting the
#'   maximum-a-posterior estimates of the parameters for a fitted DMBC model.
#'
#' @param res An object of class \code{dmbc_fit_list}.
#' @param chain A length-one numeric vector indicating the MCMC chain number
#'   to use.
#'
#' @return A named \code{list} with the following elements:
#'   \describe{
#'     \item{\code{z}: }{array of latent coordinates posterior mean estimates}
#'     \item{\code{alpha}: }{numeric vector of alpha posterior mean estimates}
#'     \item{\code{eta}: }{numeric vector of eta posterior mean estimates}
#'     \item{\code{sigma2}: }{numeric vector of sigma2 posterior mean estimates}
#'     \item{\code{lambda}: }{numeric vector of lambda posterior mean estimates}
#'     \item{\code{prob}: }{numeric matrix of probability posterior mean estimates}
#'     \item{\code{cluster}: }{numeric vector of cluster membership posterior
#'       mean estimates}
#'     \item{\code{logpost}: }{length-one numeric vector of the maximum
#'       log-posterior value}
#'     \item{\code{chain}: }{length-one numeric vector of the MCMC chain number
#'       used}
#'   }
#'
#' @author Sergio Venturini \email{sergio.venturini@@unibocconi.it}
#'
#' @seealso \code{\link{dmbc_data}} for a description of the data format.
#' @seealso \code{\link{dmbc_fit_list}} for a description of the elements
#'   included in the returned object.
#'
#' @references
#'   Venturini, S., Piccarreta, R. (2019), "A Bayesian Approach for Model-Based
#'   Clustering of Several Binary Dissimilarity Matrices: the \pkg{dmbc}
#'   Package in \code{R}", Technical report.
#'
#' @examples
#' data(simdiss, package = "dmbc")
#'
#' G <- 3
#' p <- 2
#' prm.prop <- list(z = 1.5, alpha = .75)
#' burnin <- 2000
#' nsim <- 1000
#' seed <- 2301
#'
#' set.seed(seed)
#'
#' control <- list(burnin = burnin, nsim = nsim, z.prop = prm.prop[["z"]],
#'   alpha.prop = prm.prop[["alpha"]], random.start = TRUE, verbose = TRUE,
#'   nchains = 2, thin = 10, store.burnin = TRUE, threads = 2,
#'   parallel = "multicore")
#' sim.dmbc <- dmbc(simdiss, p, G, control)
#'
#' dmbc_get_map(sim.dmbc, chain = 1)
#' @export
dmbc_get_map <- function(res, chain = 1) {
  res_chain <- res@results[[chain]]
  control <- res@results[[1]]@control
  burnin <- control[["burnin"]]
  nsim <- control[["nsim"]]
  thin <- control[["thin"]]
  store.burnin <- control[["store.burnin"]]
  nchains <- control[["nchains"]]
  totiter <- burnin + nsim
	lpost <- res_chain@dens$logpost
  
  if (chain > nchains)
    stop("the specified chain is not available.")

  if (store.burnin) {
    todrop <- seq(1, burnin, by = thin)
    tokeep <- seq(1, totiter, by = thin)
    tokeep <- (length(todrop) + 1):length(tokeep)
  } else {
    tokeep <- seq(1, nsim, by = thin)
    tokeep <- 1:length(tokeep)
  }

  i.map <- which.max(lpost[tokeep])
  if (store.burnin) {
    i.map <- length(todrop) + i.map
  }

	prob.map <- res_chain@prob.chain[i.map, , ]
  out <- list(z = res_chain@z.chain.p[i.map, , , ],
              alpha = res_chain@alpha.chain[i.map, ],
          		eta = res_chain@eta.chain[i.map, ],
              sigma2 = res_chain@sigma2.chain[i.map, ],
          		lambda = res_chain@lambda.chain[i.map, ],
              prob = prob.map,
              cluster = apply(prob.map, 1, which.max),
              logpost = lpost[i.map],
              chain = chain)

	return(out)
}

#' Extractor function for a fitted DMBC model.
#'
#' \code{dmbc_get_configuration()} is an extractor function for extracting the
#'   latent configuration estimates of a fitted DMBC model.
#'
#' @param res An object of class \code{dmbc_fit_list}.
#' @param chain A length-one numeric vector indicating the MCMC chain number
#'   to use.
#' @param est A length-one character vector indicating the estimate type to use.
#' @param labels An optional character vector with the object labels.
#'
#' @return A \code{\link{dmbc_config}} object.
#'
#' @author Sergio Venturini \email{sergio.venturini@@unibocconi.it}
#'
#' @seealso \code{\link{dmbc_data}} for a description of the data format.
#' @seealso \code{\link{dmbc_fit_list}} for a description of the elements
#'   included in the returned object.
#'
#' @references
#'   Venturini, S., Piccarreta, R. (2019), "A Bayesian Approach for Model-Based
#'   Clustering of Several Binary Dissimilarity Matrices: the \pkg{dmbc}
#'   Package in \code{R}", Technical report.
#'
#' @examples
#' data(simdiss, package = "dmbc")
#'
#' G <- 3
#' p <- 2
#' prm.prop <- list(z = 1.5, alpha = .75)
#' burnin <- 2000
#' nsim <- 1000
#' seed <- 2301
#'
#' set.seed(seed)
#'
#' control <- list(burnin = burnin, nsim = nsim, z.prop = prm.prop[["z"]],
#'   alpha.prop = prm.prop[["alpha"]], random.start = TRUE, verbose = TRUE,
#'   nchains = 2, thin = 10, store.burnin = TRUE, threads = 2,
#'   parallel = "multicore")
#' sim.dmbc <- dmbc(simdiss, p, G, control)
#'
#' z <- dmbc_get_configuration(sim.dmbc, chain = 1, est = "mean")
#' summary(z)
#'
#' library(bayesplot)
#' library(ggplot2)
#' color_scheme_set("mix-pink-blue")
#' graph <- plot(z, size = 2, size_lbl = 3, label_object = TRUE, adjust = .1)
#' graph <- graph + panel_bg(fill = "gray90", color = NA)
#' graph + geom_text(aes(label = lbl), nudge_x = .75, nudge_y = 0, size = 3)
#' @export
dmbc_get_configuration <- function(res, chain = 1, est = "mean", labels = character(0)) {
  res_chain <- res@results[[chain]]
  control <- res@results[[1]]@control
  burnin <- control[["burnin"]]
  nsim <- control[["nsim"]]
  thin <- control[["thin"]]
  nchains <- control[["nchains"]]
  store.burnin <- control[["store.burnin"]]
  totiter <- burnin + nsim
  
  if (chain > nchains)
    stop("the specified chain is not available.")
  if (!(est %in% c("mean", "median", "ml", "map")))
    stop("the estimate type specified is not available.")
  labels <- as.character(labels)
  if (length(labels) && (length(labels) != res_chain@dim[["n"]]))
    stop("the number of labels provided must be equal to the number of objects in the data.")

  if (store.burnin) {
    todrop <- seq(1, burnin, by = thin)
    tokeep <- seq(1, totiter, by = thin)
    tokeep <- (length(todrop) + 1):length(tokeep)
  } else {
    tokeep <- seq(1, nsim, by = thin)
    tokeep <- 1:length(tokeep)
  }

  res.est <- switch(est,
    mean = dmbc_get_postmean(res, chain = chain),
    median = dmbc_get_postmedian(res, chain = chain),
    ml = dmbc_get_ml(res, chain = chain),
    map = dmbc_get_map(res, chain = chain))
  Z.est <- res.est$z
  Z.sd <- apply(res_chain@z.chain.p[tokeep, , , , drop = FALSE], c(2, 3, 4), sd, na.rm = TRUE)
  cl <- res.est$cluster

  out <- new("dmbc_config",
    Z.est = Z.est,
    Z.sd = Z.sd,
    cluster = cl,
    est = est,
    n = res_chain@dim[["n"]],
    p = res_chain@dim[["p"]],
    S = res_chain@dim[["S"]],
    G = res_chain@dim[["G"]],
    family = res_chain@model@family,
    chain = chain,
    labels = labels)

  return(out)
}
