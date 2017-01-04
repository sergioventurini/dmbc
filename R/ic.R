dmbc.IC <- function(D, pmax = 3, Gmax = 5, burnin, nsim, prm.prop, prm.prior, random.start, est = "post.mean", verbose = TRUE) {
	logprior <- logmlik <- logcorrfact <- dcic <- matrix(NA, nrow = pmax, ncol = Gmax)
	res_list <- list()
	res_save <- list()
	res_all <- list()
	res.i <- 1
	for (G.i in 1:Gmax) {
		for (p.i in 1:pmax) {
			if (verbose)
				cat("*** p = ", p.i, " -- G = ", G.i, " ***\n", sep = "")
			
			res <- dmbc(D = D, p = p.i, G = G.i, burnin = burnin, nsim = nsim, prm.prop = prm.prop, prm.prior = list(sigma2 = prm.prior,
							lambda = rep(1, G.i)), random.start = random.start, verbose = verbose)
			res_list[[p.i]] <- res
			if (p.i == pmax)
				res_save[[G.i]] <- res
			
			if (est == "post.mean") {
				est.tmp <- dmbc.get.postmean(res)
				z.m <- est.tmp$z.postmean
				alpha.m <- est.tmp$alpha.postmean
				eta.m <- est.tmp$eta.postmean
				sigma2.m <- est.tmp$sigma2.postmean
				lambda.m <- est.tmp$lambda.postmean
			} else if (est == "ml") {
				est.tmp <- dmbc.get.ml(res)
				z.m <- est.tmp$z.ml
				alpha.m <- est.tmp$alpha.ml
				eta.m <- est.tmp$eta.ml
				sigma2.m <- est.tmp$sigma2.ml
				lambda.m <- est.tmp$lambda.ml
			} else if (est == "map") {
				est.tmp <- dmbc.get.map(res)
				z.m <- est.tmp$z.map
				alpha.m <- est.tmp$alpha.map
				eta.m <- est.tmp$eta.map
				sigma2.m <- est.tmp$sigma2.map
				lambda.m <- est.tmp$lambda.map
			}
			class.m <- dmbc.get.postmean(res)$class.postmean
			res_all[[res.i]] <- list(z.m = z.m, alpha.m = alpha.m, eta.m = eta.m, sigma2.m = sigma2.m, lambda.m = lambda.m, class.m = class.m)
			names(res_all)[res.i] <- paste("p = ", p.i, " -- G = ", G.i, sep = "")

			res.i <- res.i + 1
			
			logprior[p.i, G.i] <- log_marg_prior(res, z.m)
			logmlik[p.i, G.i] <- log_marg_lik(res, z.m)
			if (p.i > 1)
				logcorrfact[p.i, G.i] <- log_corr_fact(res_list[[p.i - 1]], z.m)
			dcic[p.i, G.i] <-  -2*(logprior[p.i, G.i] + logmlik[p.i, G.i] + ifelse(p.i > 1, sum(logcorrfact[2:p.i, G.i], na.rm = TRUE), 0))
			
			if (verbose) {
				cat("\n")
			}
		}
		res_list <- list()
	}
	
	out <- new("dmbcIC",
		logprior = logprior,
		logmlik = logmlik,
		logcorrfact = logcorrfact,
		DCIC = dcic,
		post.mean = res_all,
		control = list(burnin = burnin, nsim = nsim, prm.prop = prm.prop, prm.prior = prm.prior, random.start = random.start, est = est),
		res_last_p = res_save
	)
	return(out)
}

dmbc.IC.more <- function(D, pmax = NULL, Gmax = NULL, IC.res, verbose = TRUE) {
	pmax.old <- nrow(IC.res@DCIC)
	Gmax.old <- ncol(IC.res@DCIC)
	if (is.null(pmax) & is.null(Gmax))
		stop("pmax and Gmax cannot be both null.")
	if (is.null(pmax)) {
		warning("you did not provide a pmax value, so it is set to its previous value.", call. = FALSE, immediate. = TRUE)
		pmax <- pmax.old
	}
	if (is.null(Gmax)) {
		warning("you did not provide a Gmax value, so it is set to its previous value.", call. = FALSE, immediate. = TRUE)
		Gmax <- Gmax.old
	}
	if ((pmax <= pmax.old) & (Gmax <= Gmax.old))
		stop("at least one of the new values for pmax and Gmax must be larger than the previous ones.")
	if (pmax < pmax.old) {
		warning("the pmax value (", pmax, ") is smaller than the previous one (", pmax.old, ") and hence it is set at that value.",
					call. = FALSE, immediate. = TRUE)
		pmax <- pmax.old
	}
	if (Gmax < Gmax.old) {
		warning("the Gmax value (", Gmax, ") is smaller than the previous one (", Gmax.old, ") and hence it is set at that value.",
					call. = FALSE, immediate. = TRUE)
		Gmax <- Gmax.old
	}
	burnin <- IC.res@control$burnin
	nsim <- IC.res@control$nsim
	prm.prop <- IC.res@control$prm.prop
	prm.prior <- IC.res@control$prm.prior
	random.start <- IC.res@control$random.start
	est <- IC.res@control$est
	logprior <- logmlik <- logcorrfact <- dcic <- matrix(NA, nrow = pmax, ncol = Gmax)
	logprior[1:pmax.old, 1:Gmax.old] <- IC.res@logprior
	logmlik[1:pmax.old, 1:Gmax.old] <- IC.res@logmlik
	logcorrfact[1:pmax.old, 1:Gmax.old] <- IC.res@logcorrfact
	dcic[1:pmax.old, 1:Gmax.old] <- IC.res@DCIC
	res_list <- list()
	res_save <- list()
	res_all <- IC.res@post.mean
	res_last_p <- IC.res@res_last_p
	res.i <- pmax.old*Gmax.old + 1
	p.seq <- if (pmax == pmax.old) pmax.old else seq(from = (pmax.old + 1), to = pmax, by = 1)
	G.seq <- if (Gmax == Gmax.old) Gmax.old else seq(from = (Gmax.old + 1), to = Gmax, by = 1)
	
	if (Gmax > Gmax.old) {
		for (G.i in G.seq) {
			for (p.i in 1:pmax.old) {
				if (verbose)
					cat("*** p = ", p.i, " -- G = ", G.i, " ***\n", sep = "")
				
				res <- dmbc(D = D, p = p.i, G = G.i, burnin = burnin, nsim = nsim, prm.prop = prm.prop, prm.prior = list(sigma2 = prm.prior,
								lambda = rep(1, G.i)), random.start = random.start, verbose = verbose)
				res_list[[p.i]] <- res
				if (p.i == pmax.old)
					res_last_p[[G.i]] <- res
				
				if (est == "post.mean") {
					est.tmp <- dmbc.get.postmean(res)
					z.m <- est.tmp$z.postmean
					alpha.m <- est.tmp$alpha.postmean
					eta.m <- est.tmp$eta.postmean
					sigma2.m <- est.tmp$sigma2.postmean
					lambda.m <- est.tmp$lambda.postmean
				} else if (est == "ml") {
					est.tmp <- dmbc.get.ml(res)
					z.m <- est.tmp$z.ml
					alpha.m <- est.tmp$alpha.ml
					eta.m <- est.tmp$eta.ml
					sigma2.m <- est.tmp$sigma2.ml
					lambda.m <- est.tmp$lambda.ml
				} else if (est == "map") {
					est.tmp <- dmbc.get.map(res)
					z.m <- est.tmp$z.map
					alpha.m <- est.tmp$alpha.map
					eta.m <- est.tmp$eta.map
					sigma2.m <- est.tmp$sigma2.map
					lambda.m <- est.tmp$lambda.map
				}
				class.m <- dmbc.get.postmean(res)$class.postmean
				res_all[[res.i]] <- list(z.m = z.m, alpha.m = alpha.m, eta.m = eta.m, sigma2.m = sigma2.m, lambda.m = lambda.m, class.m = class.m)
				names(res_all)[res.i] <- paste("p = ", p.i, " -- G = ", G.i, sep = "")

				res.i <- res.i + 1
				
				logprior[p.i, G.i] <- log_marg_prior(res, z.m)
				logmlik[p.i, G.i] <- log_marg_lik(res, z.m)
				if (p.i > 1)
					logcorrfact[p.i, G.i] <- log_corr_fact(res_list[[p.i - 1]], z.m)
				dcic[p.i, G.i] <-  -2*(logprior[p.i, G.i] + logmlik[p.i, G.i] + ifelse(p.i > 1, sum(logcorrfact[2:p.i, G.i], na.rm = TRUE), 0))
				
				if (verbose) {
					cat("\n")
				}
			}
			res_list <- list()
		}
	}

	if (pmax > pmax.old) {
		for (G.i in 1:Gmax) {
			for (p.i in p.seq) {
				if (verbose)
					cat("*** p = ", p.i, " -- G = ", G.i, " ***\n", sep = "")
				
				res <- dmbc(D = D, p = p.i, G = G.i, burnin = burnin, nsim = nsim, prm.prop = prm.prop, prm.prior = list(sigma2 = prm.prior,
								lambda = rep(1, G.i)), random.start = random.start, verbose = verbose)
				res_list[[p.i]] <- res
				if (p.i == p.seq[1])
					res_list[[p.i - 1]] <- res_last_p[[G.i]]
				if (p.i == pmax)
					res_save[[G.i]] <- res
				
				if (est == "post.mean") {
					est.tmp <- dmbc.get.postmean(res)
					z.m <- est.tmp$z.postmean
					alpha.m <- est.tmp$alpha.postmean
					eta.m <- est.tmp$eta.postmean
					sigma2.m <- est.tmp$sigma2.postmean
					lambda.m <- est.tmp$lambda.postmean
				} else if (est == "ml") {
					est.tmp <- dmbc.get.ml(res)
					z.m <- est.tmp$z.ml
					alpha.m <- est.tmp$alpha.ml
					eta.m <- est.tmp$eta.ml
					sigma2.m <- est.tmp$sigma2.ml
					lambda.m <- est.tmp$lambda.ml
				} else if (est == "map") {
					est.tmp <- dmbc.get.map(res)
					z.m <- est.tmp$z.map
					alpha.m <- est.tmp$alpha.map
					eta.m <- est.tmp$eta.map
					sigma2.m <- est.tmp$sigma2.map
					lambda.m <- est.tmp$lambda.map
				}
				class.m <- dmbc.get.postmean(res)$class.postmean
				res_all[[res.i]] <- list(z.m = z.m, alpha.m = alpha.m, eta.m = eta.m, sigma2.m = sigma2.m, lambda.m = lambda.m, class.m = class.m)
				names(res_all)[res.i] <- paste("p = ", p.i, " -- G = ", G.i, sep = "")

				res.i <- res.i + 1
				
				logprior[p.i, G.i] <- log_marg_prior(res, z.m)
				logmlik[p.i, G.i] <- log_marg_lik(res, z.m)
				logcorrfact[p.i, G.i] <- log_corr_fact(res_list[[p.i - 1]], z.m)
				dcic[p.i, G.i] <-  -2*(logprior[p.i, G.i] + logmlik[p.i, G.i] + ifelse(p.i > 1, sum(logcorrfact[2:p.i, G.i], na.rm = TRUE), 0))
				
				if (verbose) {
					cat("\n")
				}
			}
			res_list <- list()
		}
	}
	
	out <- new("dmbcIC",
		logprior = logprior,
		logmlik = logmlik,
		logcorrfact = logcorrfact,
		DCIC = dcic,
		post.mean = res_all,
		control = list(burnin = burnin, nsim = nsim, prm.prop = prm.prop, prm.prior = prm.prior, random.start = random.start, est = est),
		res_last_p = res_save
	)
	return(out)
}

log_marg_lik <- function(res, Z) {
	D <- res@obsdiss
	burnin <- res@control[[1]]
	nsim <- res@control[[2]]
	totiter <- burnin + nsim
	x <- dmbc.get.postmean(res)$class.postmean
	n <- dim(Z)[1]
	p <- dim(Z)[2]
	G <- dim(Z)[3]
	if (G > 1) {
		q <- 3*G
		theta <- cbind(res@alpha.chain[(burnin + 1):totiter, ], res@sigma2.chain[(burnin + 1):totiter, ], res@lambda.chain[(burnin + 1):totiter, ])
		theta.star <- l1median(theta)
		theta.star[(2*G + 1):(3*G)] <- colMeans(theta)[(2*G + 1):(3*G)]   # needed cause the 'ddirichlet' function returns a 0 when the sum of lambdas != 1
		H.star <- covMcd(theta[, -q])$cov   # we removed the last dimension otherwise the hessian is always singular
	} else {
		q <- 2
		theta <- cbind(res@alpha.chain[(burnin + 1):totiter, 1], res@sigma2.chain[(burnin + 1):totiter, 1])
		theta.star <- l1median(theta)
		H.star <- covMcd(theta)$cov
	}
	if (G > 1) {
		loglik <- dmbc.logLik(D, Z, theta.star[1:G], theta.star[(2*G + 1):(3*G)], x)
	} else {
		loglik <- dmbc.logLik(D, Z, theta.star[1], theta.star[2], x)
	}
	alpha <- res@control[[4]]$sigma2[1]
	beta <- res@control[[4]]$sigma2[2]
	if (G > 1) {
		lambda.hyp <- res@control[[4]]$lambda
	}
	logprior <- 0
	for (g in 1:G) {
		logprior <- logprior + dnorm(theta.star[g], sd = sqrt(theta.star[G + g]), log = TRUE)
		logprior <- logprior + dinvgamma(theta.star[G + g], alpha = alpha, beta = beta, log = TRUE)
	}
	if (G > 1) {
		logprior <- logprior + log(ddirichlet(theta.star[(2*G + 1):(3*G)], lambda.hyp))
	}
	logmlik <- q*log(2*pi)/2 + log(det(H.star))/2 + loglik + logprior

	return(logmlik)
}

log_marg_prior <- function(res, Z) {
	n <- dim(Z)[1]
	p <- dim(Z)[2]
	G <- dim(Z)[3]
	a_g <- res@control[[5]]
	b_g <- res@control[[6]]
	logprior <- 0
	for (g in 1:G) {
		logprior <- logprior + lgamma(a_g[g] + n*p/2) - lgamma(a_g[g]) + a_g[g]*log(b_g[g]) - (a_g[g] + n*p/2)*log(b_g[g] + sum(Z[, , g]^2)/2)
	}
	logprior <- logprior - n*p*G*log(2*pi)/2

	return(logprior)
}

log_corr_fact <- function(res, Z) {
	n <- dim(Z)[1]
	p <- dim(Z)[2]
	a <- res@control[[5]][1]
	b <- res@control[[6]][1]
	logcorrfact <- n*log(2*pi)/2 + lgamma(a + n*p/2) - lgamma(a + n*(p + 1)/2) + n*log(b + sum(Z^2)/2)/2

	return(logcorrfact)
}
